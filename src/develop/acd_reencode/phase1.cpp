// phase1.cpp -- Phase 1 of ../plan/acd_cut_reencoding.md: compressibility of a pipeline cut.
//
//   acd_phase1 [options] <file.aig> <cut.txt>
//     cut.txt: one cut signal per line: "<node> [arrival downstream]" (node = aig_mincut node id)
//   --K k        window leaf limit (default 16)
//   --global     also run the global sampling metric (3.1)
//   --samples n  log2 of the number of random patterns for 3.1 (default 20)
//   --slack T    slack-aware grouping (not in the plan): only group signals of the same class
//                A (upstream and downstream slack >= 1), B (only downstream slack >= 1);
//                class C (downstream slack 0, cannot be decoded) stays as singletons
//
// Output (CSV-ish lines on stdout):
//   global,<N>,<ceil_log2_N>,<saturated 0/1>,<samples>,<exhaustive 0/1>
//   group,<id>,<|G|>,<|L|>,<|Img|>,<r>,<saving>,<fd_size>,<max_arrival>,<max_downstream>,<r_eq_L>,<members>,<fd subset>,<leaves>
#include "../aig_mincut/aig.h"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

using std::vector;
typedef unsigned long long u64;

struct CutSig { int node, arrival = 0, downstream = 0; };

// ---------------------------------------------------------------------------------------------
// Reconvergence-driven window: start from the roots as leaves and repeatedly replace the leaf
// whose expansion adds the fewest new leaves (ties: the higher one). Expansions that do not grow
// the leaf set are always taken; growing ones only while |L| <= K. Expanding never enlarges the
// image, so we expand as far as the limit allows.

struct Window { vector<int> leaves, cone; bool ok = false; };

static Window findWindow(const Aig& g, const vector<int>& roots, int K) {
    std::set<int> L(roots.begin(), roots.end());
    int steps = 0;
    while (steps++ < 200000) {
        int best = -1, bestCost = 1 << 30;
        for (int v : L) {
            if (!g.isAnd(v)) continue;
            int f0 = g.fanin0[v] >> 1, f1 = g.fanin1[v] >> 1, cost = -1;
            if (f0 && !L.count(f0)) cost++;
            if (f1 && f1 != f0 && !L.count(f1)) cost++;
            if (cost < bestCost || (cost == bestCost && v > best)) { bestCost = cost; best = v; }
        }
        if (best < 0) break;
        if (bestCost > 0 && (int)L.size() + bestCost > K) break;
        L.erase(best);
        for (int lit : { g.fanin0[best], g.fanin1[best] }) if (lit >> 1) L.insert(lit >> 1);
    }
    Window w;
    w.leaves.assign(L.begin(), L.end());
    w.ok = (int)L.size() <= K;
    // cone: nodes between leaves and roots (inclusive of roots that are not leaves)
    std::set<int> seen(L.begin(), L.end());
    vector<int> st(roots.begin(), roots.end());
    while (!st.empty()) {
        int v = st.back(); st.pop_back();
        if (seen.count(v)) continue;
        seen.insert(v);
        w.cone.push_back(v);
        for (int lit : { g.fanin0[v], g.fanin1[v] }) if ((lit >> 1) && !seen.count(lit >> 1)) st.push_back(lit >> 1);
    }
    std::sort(w.cone.begin(), w.cone.end());   // node ids are topological
    return w;
}

// ---------------------------------------------------------------------------------------------
// Distinct-vector counting with a 128-bit linear hash: h(x) = XOR of key_j over set bits x_j.
// Two distinct vectors collide with probability 2^-128 per pair.

struct Keys { vector<u64> k1, k2; };
static Keys makeKeys(size_t n, u64 seed) {
    std::mt19937_64 rng(seed);
    Keys k;
    for (size_t i = 0; i < n; i++) { k.k1.push_back(rng()); k.k2.push_back(rng()); }
    return k;
}

static size_t countDistinct(vector<std::pair<u64, u64>> h) {
    std::sort(h.begin(), h.end());
    return std::unique(h.begin(), h.end()) - h.begin();
}

static int ceilLog2(u64 n) { int r = 0; while ((1ULL << r) < n) r++; return r; }

// Truth-table pattern of variable i over 64 * w ... patterns (word index w).
static u64 varWord(int i, size_t w) {
    static const u64 m[6] = { 0xAAAAAAAAAAAAAAAAULL, 0xCCCCCCCCCCCCCCCCULL, 0xF0F0F0F0F0F0F0F0ULL,
                              0xFF00FF00FF00FF00ULL, 0xFFFF0000FFFF0000ULL, 0xFFFFFFFF00000000ULL };
    if (i < 6) return m[i];
    return ((w >> (i - 6)) & 1) ? ~0ULL : 0ULL;
}

struct GroupResult {
    vector<int> members, fdSet;
    int nLeaves = 0, r = 0, fd = 0, maxArr = 0, maxDown = 0;
    u64 img = 0;
};

// Image and functional-dependence analysis of G over its window (all 2^|L| leaf values).
static GroupResult analyzeGroup(const Aig& g, const vector<int>& G, const Window& w) {
    int nL = (int)w.leaves.size();
    size_t P = 1ULL << nL, nW = (P + 63) / 64;
    std::unordered_map<int, int> idx;
    vector<vector<u64>> val;
    for (int i = 0; i < nL; i++) {
        idx[w.leaves[i]] = (int)val.size();
        vector<u64> t(nW);
        for (size_t k = 0; k < nW; k++) t[k] = varWord(i, k);
        val.push_back(std::move(t));
    }
    for (int v : w.cone) {
        if (idx.count(v)) continue;
        vector<u64> t(nW);
        int a = g.fanin0[v], b = g.fanin1[v];
        const vector<u64>* va = (a >> 1) ? &val[idx.at(a >> 1)] : nullptr;
        const vector<u64>* vb = (b >> 1) ? &val[idx.at(b >> 1)] : nullptr;
        for (size_t k = 0; k < nW; k++) {
            u64 x = va ? (*va)[k] : 0, y = vb ? (*vb)[k] : 0;
            if (a & 1) x = ~x;
            if (b & 1) y = ~y;
            t[k] = x & y;
        }
        idx[v] = (int)val.size();
        val.push_back(std::move(t));
    }
    Keys keys = makeKeys(G.size(), 0x5eed + G.size());
    vector<std::pair<u64, u64>> h(P, { 0, 0 });
    auto bitOf = [&](int j, size_t p) { return (val[idx.at(G[j])][p >> 6] >> (p & 63)) & 1; };
    for (size_t j = 0; j < G.size(); j++)
        for (size_t p = 0; p < P; p++)
            if (bitOf((int)j, p)) { h[p].first ^= keys.k1[j]; h[p].second ^= keys.k2[j]; }
    GroupResult res;
    res.members = G;
    res.nLeaves = nL;
    res.img = countDistinct(h);
    res.r = ceilLog2(res.img);
    // greedy deletion: drop a signal if it is a function of the remaining ones
    // (removing it does not merge any image points). Deepest signals are tried first.
    vector<int> order(G.size());
    for (size_t j = 0; j < G.size(); j++) order[j] = (int)j;
    std::sort(order.begin(), order.end(), [&](int x, int y) { return G[x] > G[y]; });
    vector<char> dropped(G.size(), 0);
    for (int j : order) {
        vector<std::pair<u64, u64>> h2 = h;
        for (size_t p = 0; p < P; p++)
            if (bitOf(j, p)) { h2[p].first ^= keys.k1[j]; h2[p].second ^= keys.k2[j]; }
        if (countDistinct(h2) == res.img) { h.swap(h2); dropped[j] = 1; }
    }
    for (size_t j = 0; j < G.size(); j++) if (!dropped[j]) res.fdSet.push_back(G[j]);
    res.fd = (int)res.fdSet.size();
    return res;
}

// ---------------------------------------------------------------------------------------------
// 3.1 global sampling over the transitive fanin of the cut.

static void globalSampling(const Aig& g, const vector<CutSig>& cut, int logSamples) {
    int n = g.numNodes();
    vector<char> inCone(n, 0);
    for (auto& c : cut) inCone[c.node] = 1;
    for (int v = n - 1; v > g.nPi; v--)
        if (inCone[v]) inCone[g.fanin0[v] >> 1] = inCone[g.fanin1[v] >> 1] = 1;
    vector<int> ands;
    for (int v = g.nPi + 1; v < n; v++) if (inCone[v]) ands.push_back(v);
    bool exhaustive = g.nPi <= logSamples;
    u64 samples = exhaustive ? (1ULL << g.nPi) : (1ULL << logSamples);
    size_t nW = (samples + 63) / 64;
    Keys keys = makeKeys(cut.size(), 0xC0FFEE);
    std::mt19937_64 rng(42);
    vector<u64> val(n, 0);
    vector<std::pair<u64, u64>> h;
    h.reserve(samples);
    for (size_t w = 0; w < nW; w++) {
        for (int i = 0; i < g.nPi; i++) val[i + 1] = exhaustive ? varWord(i, w) : rng();
        for (int v : ands) {
            u64 x = val[g.fanin0[v] >> 1], y = val[g.fanin1[v] >> 1];
            val[v] = ((g.fanin0[v] & 1) ? ~x : x) & ((g.fanin1[v] & 1) ? ~y : y);
        }
        std::pair<u64, u64> hw[64] = {};
        for (size_t j = 0; j < cut.size(); j++)
            for (u64 bits = val[cut[j].node]; bits; bits &= bits - 1) {
                int p = __builtin_ctzll(bits);
                hw[p].first ^= keys.k1[j]; hw[p].second ^= keys.k2[j];
            }
        size_t valid = std::min<u64>(64, samples - 64 * w);
        for (size_t p = 0; p < valid; p++) h.push_back(hw[p]);
    }
    u64 N = countDistinct(std::move(h));
    printf("global,%llu,%d,%d,%llu,%d\n", N, ceilLog2(N), (!exhaustive && N == samples) ? 1 : 0, samples, exhaustive ? 1 : 0);
}

// ---------------------------------------------------------------------------------------------

int main(int argc, char** argv) {
    int K = 16, logSamples = 20, slackT = 0;
    bool global = false;
    vector<std::string> pos;
    for (int i = 1; i < argc; i++) {
        std::string a = argv[i];
        if (a == "--K" && i + 1 < argc) K = std::stoi(argv[++i]);
        else if (a == "--samples" && i + 1 < argc) logSamples = std::stoi(argv[++i]);
        else if (a == "--global") global = true;
        else if (a == "--slack" && i + 1 < argc) slackT = std::stoi(argv[++i]);
        else pos.push_back(a);
    }
    if (pos.size() != 2 || K < 1 || K > 24) {
        fprintf(stderr, "usage: %s [--K k] [--global] [--samples n] [--slack T] <file.aig> <cut.txt>\n", argv[0]);
        return 2;
    }
    try {
        Aig g = readAiger(pos[0]);
        vector<CutSig> cut;
        std::ifstream cf(pos[1]);
        if (!cf) throw std::runtime_error("cannot open " + pos[1]);
        std::string line;
        while (std::getline(cf, line)) {
            std::istringstream ss(line);
            CutSig c;
            if (!(ss >> c.node)) continue;
            ss >> c.arrival >> c.downstream;
            if (c.node <= 0 || c.node >= g.numNodes()) throw std::runtime_error("bad cut node " + std::to_string(c.node));
            cut.push_back(c);
        }
        if (global) globalSampling(g, cut, logSamples);

        // per-signal windows, used for the overlap score
        int nc = (int)cut.size();
        vector<vector<int>> single(nc);
        for (int j = 0; j < nc; j++) single[j] = findWindow(g, { cut[j].node }, K).leaves;

        // 3.2 step 1-2: greedy grouping by leaf overlap
        vector<char> used(nc, 0);
        int gid = 0;
        vector<int> cls(nc, 0);   // 0 = A, 1 = B, 2 = C (only with --slack)
        if (slackT > 0)
            for (int j = 0; j < nc; j++)
                cls[j] = cut[j].downstream >= slackT ? 2 : (cut[j].arrival >= slackT ? 1 : 0);
        for (int seed = 0; seed < nc; seed++) {
            if (used[seed]) continue;
            used[seed] = 1;
            vector<int> G{ cut[seed].node }, Gi{ seed };
            Window w = findWindow(g, G, K);
            while (cls[seed] != 2) {
                std::set<int> Ls(w.leaves.begin(), w.leaves.end());
                int best = -1, bestOv = 0;
                for (int j = 0; j < nc; j++) {
                    if (used[j] || cls[j] != cls[seed]) continue;
                    int ov = 0;
                    for (int l : single[j]) ov += Ls.count(l);
                    if (ov > bestOv) { bestOv = ov; best = j; }
                }
                if (best < 0) break;                       // nothing shares a leaf
                vector<int> G2 = G; G2.push_back(cut[best].node);
                Window w2 = findWindow(g, G2, K);
                if (!w2.ok) break;                         // leaves would exceed K
                used[best] = 1; G.swap(G2); Gi.push_back(best); w = w2;
            }
            if (!w.ok) throw std::logic_error("single-signal window exceeds K");
            GroupResult r = analyzeGroup(g, G, w);
            for (int j : Gi) { r.maxArr = std::max(r.maxArr, cut[j].arrival); r.maxDown = std::max(r.maxDown, cut[j].downstream); }
            printf("group,%d,%zu,%d,%llu,%d,%d,%d,%d,%d,%d,", gid++, G.size(), r.nLeaves, r.img, r.r,
                   (int)G.size() - r.r, r.fd, r.maxArr, r.maxDown, r.r == r.nLeaves ? 1 : 0);
            for (size_t k = 0; k < G.size(); k++) printf("%s%d", k ? " " : "", G[k]);
            printf(",");
            for (size_t k = 0; k < r.fdSet.size(); k++) printf("%s%d", k ? " " : "", r.fdSet[k]);
            printf(",");
            for (size_t k = 0; k < w.leaves.size(); k++) printf("%s%d", k ? " " : "", w.leaves[k]);
            printf("\n");
        }
        return 0;
    } catch (const std::exception& e) {
        fprintf(stderr, "error: %s\n", e.what());
        return 2;
    }
}
