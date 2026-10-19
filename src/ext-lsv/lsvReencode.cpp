// lsvReencode.cpp -- lsv_reencode: window-based random re-encoding of output equivalence classes
// (src/develop/plan/window_reencode_pareto.md).
//
// Window (root, K): reconvergence-driven leaves L (|L| <= K) above the root (findWindow), then the
//   cone C = every AND whose fanins are all in L u C (forward closure of L, at most -C nodes, so C
//   may contain nodes that do not feed the root), and the outputs G = nodes of C used outside C or
//   by a PO. Replacing the window by any circuit L -> G removes all of C.
// Classes: the 2^|L| leaf patterns are grouped by the value vector of G (analyzeGroup); N classes,
//   r = ceil(log2 N).
// Re-encoding: an injection phi: class -> r-bit codeword. Encoder E: L -> r bits (truth table bit p
//   = phi(class(p))), decoder D: r bits -> G with unused codewords as don't cares; both synthesized
//   from truth tables (synthesize: ISOP / DC fills x {raw, dc2, syn2}). Two phi that differ only by
//   permuting or complementing code bits are isomorphic; every phi is reduced to a canonical form
//   (class 0 -> 0, code-bit columns sorted) and costed once.
// Samplers: base (assignCode, deterministic), U (uniform phi), A (affine A*base, A invertible and not
//   a permutation), W1/W4 (1 or 4 random swaps/moves from base), X (all canonical phi, N <= 8), and
//   P0 (direct resynthesis of L -> G, no code bottleneck), and N (kicks only: null kick, C put back as
//   it was; the AIG is the same up to node order).
// gain = |C| - AND(E) - AND(D); dlev = max over G of (new arrival - old arrival); ddep = increase of
//   the circuit depth if only this window changed.
// Phase 1 (-a): analysis only, one CSV row per (window, sampler); the current GIA is not changed.
// Phase 2 (-k): kicks = window replacement(s) + script R, compared with R(x); -m / -P apply several
//   windows per kick (each window is picked on the AIG left by the previous one).
// Phase 3a (-i): iterated local search, kick + R from the current AIG, keep improvements; -L runs R only on
//   a region around the replacement.
// Both -k and -i report how many code bits survive R (simulation signatures); -b keeps them (R cannot
//   optimize across the code bits).
#include "lsvAcdCore.h"
#include "lsvPipe.h"

#include <climits>
#include <cstdlib>
#include <ctime>
#include <map>
#include <numeric>
#include <queue>
#include <random>
#include <set>

#include "base/abc/abc.h"
#include "base/main/main.h"
#include "base/main/mainInt.h"
#include "misc/extra/extra.h"

namespace lsvpipe {
namespace {

using std::string;
using std::vector;
typedef vector<vector<word>> Tts;

// ---------------------------------------------------------------------------------------------
// Levels of the whole AIG

struct Levels {
    vector<int> lev, down;   // down: longest path (in ANDs) from the node to a PO; -1 if none
    vector<char> poDrv;
    int depth = 0;
};

Levels computeLevels(const PipeAig& g) {
    int n = g.numNodes();
    Levels L;
    L.lev.assign(n, 0); L.down.assign(n, -1); L.poDrv.assign(n, 0);
    for (int v = g.nPi + 1; v < n; v++) L.lev[v] = 1 + std::max(L.lev[g.fanin0[v] >> 1], L.lev[g.fanin1[v] >> 1]);
    for (int lit : g.poLit) { L.poDrv[lit >> 1] = 1; L.depth = std::max(L.depth, L.lev[lit >> 1]); }
    for (int v = n - 1; v >= 1; v--) {
        if (L.poDrv[v]) L.down[v] = 0;
        for (int k = g.foStart[v]; k < g.foStart[v + 1]; k++) {
            int x = g.foList[k];
            if (L.down[x] >= 0) L.down[v] = std::max(L.down[v], L.down[x] + 1);
        }
    }
    return L;
}

// ---------------------------------------------------------------------------------------------
// Windows: leaves, forward closure, output boundary

struct RWin {
    int root = 0;
    vector<int> leaves, cone, outs;   // cone sorted (topological)
    bool trunc = false;               // closure stopped at the cone limit
    bool shrunk = false;              // nodes feeding a leaf were removed from the closure
};

struct Marks {   // stamp-based node sets, reset in O(1)
    vector<int> leaf, cone, queued, vis;
    int cur = 0;
    explicit Marks(int n) : leaf(n, 0), cone(n, 0), queued(n, 0), vis(n, 0) {}
};

enum WinStatus { WIN_OK, WIN_NONE, WIN_CYCLIC };

// WIN_NONE: no window within K leaves, or the cone limit was hit before the root.
// WIN_CYCLIC: some leaf is in the transitive fanout of C; replacing C would then create a loop
// (E reads that leaf, the leaf reads an output of C).
WinStatus buildWindow(const PipeAig& g, const Levels& LV, int root, int K, int maxCone, Marks& M, RWin& w) {
    Window f = findWindow(g, { root }, K);
    if (!f.ok) return WIN_NONE;
    int cur = ++M.cur;
    w = RWin();
    w.root = root;
    w.leaves = f.leaves;
    for (int l : w.leaves) M.leaf[l] = cur;
    auto inW = [&](int v) { return v == 0 || M.leaf[v] == cur || M.cone[v] == cur; };
    std::priority_queue<int, vector<int>, std::greater<int>> pq;   // smallest id first = topological
    auto pushFanouts = [&](int v) {
        for (int k = g.foStart[v]; k < g.foStart[v + 1]; k++) {
            int x = g.foList[k];
            if (M.queued[x] != cur && M.leaf[x] != cur) { M.queued[x] = cur; pq.push(x); }
        }
    };
    for (int l : w.leaves) pushFanouts(l);
    while (!pq.empty()) {
        if ((int)w.cone.size() >= maxCone) { w.trunc = true; break; }
        int v = pq.top(); pq.pop();
        if (!inW(g.fanin0[v] >> 1) || !inW(g.fanin1[v] >> 1)) continue;
        M.cone[v] = cur;
        w.cone.push_back(v);
        pushFanouts(v);
    }
    if (M.cone[root] != cur) return WIN_NONE;   // truncated before the root
    // A leaf must not depend on C: E reads the leaves and D drives G, so a leaf that reads a node of C
    // would close a loop. Drop every node of C in the transitive fanin of a leaf (search only above the
    // smallest cone node; nothing below it can depend on C), then every node of C that reads a dropped
    // node, so that C stays closed. If the root is dropped the window is cyclic.
    int minC = w.cone.front(), nDrop = 0;
    vector<int> st;
    for (int l : w.leaves) if (l > minC && g.isAnd(l)) { M.vis[l] = cur; st.push_back(l); }
    while (!st.empty()) {
        int v = st.back(); st.pop_back();
        for (int lit : { g.fanin0[v], g.fanin1[v] }) {
            int x = lit >> 1;
            if (x <= minC - 1 || !g.isAnd(x) || M.vis[x] == cur) continue;
            M.vis[x] = cur;
            if (M.cone[x] == cur) { M.cone[x] = 0; nDrop++; }
            st.push_back(x);
        }
    }
    if (nDrop) {
        vector<int> kept;
        for (int v : w.cone) {   // topological order: fanins are decided first
            if (M.cone[v] != cur) continue;
            if (!inW(g.fanin0[v] >> 1) || !inW(g.fanin1[v] >> 1)) { M.cone[v] = 0; continue; }
            kept.push_back(v);
        }
        w.cone.swap(kept);
        w.shrunk = true;
        if (M.cone[root] != cur) return WIN_CYCLIC;
    }
    for (int v : w.cone) {
        bool out = LV.poDrv[v];
        for (int k = g.foStart[v]; k < g.foStart[v + 1] && !out; k++) out = M.cone[g.foList[k]] != cur;
        if (out) w.outs.push_back(v);
    }
    return WIN_OK;
}

// Truth tables (over the leaf patterns) of the window outputs, by direct simulation of the cone.
Tts outputTruths(const PipeAig& g, const RWin& w) {
    int nL = (int)w.leaves.size(), nw = nWordsOf(nL);
    std::map<int, vector<word>> val;
    val[0] = vector<word>(nw, 0);
    for (int i = 0; i < nL; i++) {
        vector<word> t(nw);
        for (int k = 0; k < nw; k++) t[k] = varWord(i, k);
        val[w.leaves[i]] = t;
    }
    for (int v : w.cone) {
        const vector<word>& a = val.at(g.fanin0[v] >> 1);
        const vector<word>& b = val.at(g.fanin1[v] >> 1);
        vector<word> t(nw);
        for (int k = 0; k < nw; k++) t[k] = ((g.fanin0[v] & 1) ? ~a[k] : a[k]) & ((g.fanin1[v] & 1) ? ~b[k] : b[k]);
        val[v] = t;
    }
    Tts out;
    for (int v : w.outs) out.push_back(val.at(v));
    return out;
}

// Equivalence classes of the window outputs over the leaf patterns.
Group groupOf(const PipeAig& g, const RWin& w) {
    Group G;
    G.members = w.outs;
    Window cw; cw.leaves = w.leaves; cw.cone = w.cone; cw.ok = true;
    analyzeGroup(g, G, cw, false);
    return G;
}

// Window roots: the ANDs that reach a PO.
vector<int> liveRoots(const PipeAig& g, const Levels& LV) {
    vector<int> roots;
    for (int v = g.nPi + 1; v < g.numNodes(); v++) if (LV.down[v] >= 0) roots.push_back(v);
    return roots;
}

// ---------------------------------------------------------------------------------------------
// Small GIA helpers

// Arrival time of every CO when CI i arrives at ciArr[i].
vector<int> outArrivals(Gia_Man_t* p, const vector<int>& ciArr) {
    vector<int> a(Gia_ManObjNum(p), 0);
    Gia_Obj_t* pObj;
    int i;
    Gia_ManForEachCi(p, pObj, i) a[Gia_ObjId(p, pObj)] = i < (int)ciArr.size() ? ciArr[i] : 0;
    Gia_ManForEachAnd(p, pObj, i) a[i] = 1 + std::max(a[Gia_ObjFaninId0(pObj, i)], a[Gia_ObjFaninId1(pObj, i)]);
    vector<int> out;
    Gia_ManForEachCo(p, pObj, i) out.push_back(a[Gia_ObjFaninId0p(p, pObj)]);
    return out;
}

// ---------------------------------------------------------------------------------------------
// Codes: phi as code[q] = codeword of class q

typedef vector<int> Code;

// Canonical representative under permutation / complementation of the r code bits: XOR with
// code[0] (class 0 -> 0), then sort the r columns (bit b of every class) lexicographically.
// For N > 2^(r-1) no two columns are equal, so the representative is unique.
Code canonical(const Code& c, int r) {
    int N = (int)c.size();
    vector<vector<char>> cols(r, vector<char>(N));
    for (int b = 0; b < r; b++) for (int q = 0; q < N; q++) cols[b][q] = ((c[q] ^ c[0]) >> b) & 1;
    std::sort(cols.begin(), cols.end());
    Code out(N, 0);
    for (int b = 0; b < r; b++) for (int q = 0; q < N; q++) out[q] |= cols[b][q] << b;
    return out;
}

Code sampleUniform(int N, int r, std::mt19937_64& rng) {
    vector<int> perm(1 << r);
    std::iota(perm.begin(), perm.end(), 0);
    std::shuffle(perm.begin(), perm.end(), rng);
    return Code(perm.begin(), perm.begin() + N);
}

bool invertible(vector<int> rows, int r) {
    for (int col = 0; col < r; col++) {
        int piv = -1;
        for (int i = col; i < r; i++) if ((rows[i] >> col) & 1) { piv = i; break; }
        if (piv < 0) return false;
        std::swap(rows[col], rows[piv]);
        for (int i = 0; i < r; i++) if (i != col && ((rows[i] >> col) & 1)) rows[i] ^= rows[col];
    }
    return true;
}

// A * base over GF(2), A random invertible and not a permutation matrix (rows = output bits).
Code sampleAffine(const Code& base, int r, std::mt19937_64& rng) {
    vector<int> rows(r);
    while (true) {
        for (auto& x : rows) x = (int)(rng() & ((1u << r) - 1));
        bool perm = true;
        for (int x : rows) perm &= __builtin_popcount(x) == 1;
        if (!perm && invertible(rows, r)) break;
    }
    Code c(base.size(), 0);
    for (size_t q = 0; q < base.size(); q++)
        for (int b = 0; b < r; b++) c[q] |= (__builtin_popcount(rows[b] & base[q]) & 1) << b;
    return c;
}

// s random steps from base: swap the codewords of two classes, or move a class to an unused codeword.
Code sampleWalk(const Code& base, int r, int s, std::mt19937_64& rng) {
    int N = (int)base.size();
    Code canonBase = canonical(base, r), c;
    for (int attempt = 0; attempt < 20; attempt++) {
        c = base;
        vector<char> used(1 << r, 0);
        for (int x : c) used[x] = 1;
        for (int step = 0; step < s; step++) {
            int q = (int)(rng() % N);
            if (N < (1 << r) && (rng() & 1)) {
                int w;
                do w = (int)(rng() % (1 << r)); while (used[w]);
                used[c[q]] = 0; used[w] = 1; c[q] = w;
            } else {
                int q2;
                do q2 = (int)(rng() % N); while (q2 == q);
                std::swap(c[q], c[q2]);
            }
        }
        if (canonical(c, r) != canonBase) break;
    }
    return c;
}

// All canonical codes (class 0 -> 0, the remaining classes on distinct nonzero codewords).
void enumerateCanonical(int N, int r, vector<Code>& out) {
    Code c(N, 0);
    vector<char> used(1 << r, 0);
    used[0] = 1;
    std::function<void(int)> rec = [&](int q) {
        if (q == N) { if (canonical(c, r) == c) out.push_back(c); return; }
        for (int w = 1; w < (1 << r); w++) {
            if (used[w]) continue;
            used[w] = 1; c[q] = w;
            rec(q + 1);
            used[w] = 0;
        }
    };
    rec(1);
}

u64 numCanonical(int N, int r) {   // (2^r)! / ((2^r - N)! r! 2^r)
    u64 x = 1;
    for (int i = 0; i < N; i++) x *= (u64)((1 << r) - i);
    for (int i = 2; i <= r; i++) x /= i;
    return x >> r;
}

enum Sampler { S_U, S_A, S_W1, S_W4, S_BASE, S_P0, S_NULL };
const char* kMethods[] = { "U", "A", "W1", "W4", "base", "P0", "N" };   // indexed by Sampler

// The deterministic code of assignCode.
Code baseCode(const Group& G) {
    vector<int> reused, idx;
    int nnew;
    Code base;
    assignCode(G, reused, nnew, idx, base);
    return base;
}

// A code of sampler m (S_U .. S_BASE).
Code sampleCode(int m, const Group& G, const Code& base, std::mt19937_64& rng) {
    switch (m) {
    case S_U: return sampleUniform(G.nImg, G.r, rng);
    case S_A: return sampleAffine(base, G.r, rng);
    case S_W1: return sampleWalk(base, G.r, 1, rng);
    case S_W4: return sampleWalk(base, G.r, 4, rng);
    default: return base;
    }
}

// ---------------------------------------------------------------------------------------------
// Costing one code / the direct resynthesis

struct Pt { int gain = 0, dlev = 0, ddep = 0; };

struct Cost {
    bool ok = false;
    Pt A, D;   // min-AND encoder + min-AND decoder / min-depth encoder + min-depth decoder
};

struct WinCtx {
    const PipeAig& g;
    const Levels& LV;
    const RWin& w;
    const Group& G;
    const Tts& outTt;     // G over the leaf patterns
    int budget;
    int nErrors = 0;
};

void delayOf(const WinCtx& X, const vector<int>& arr, int& dlev, int& ddep) {
    dlev = INT_MIN; int newDepth = 0;
    for (size_t j = 0; j < X.w.outs.size(); j++) {
        int v = X.w.outs[j];
        dlev = std::max(dlev, arr[j] - X.LV.lev[v]);
        if (X.LV.down[v] >= 0) newDepth = std::max(newDepth, arr[j] + X.LV.down[v]);
    }
    ddep = std::max(0, newDepth - X.LV.depth);
}

vector<int> leafArrivals(const WinCtx& X) {
    vector<int> a;
    for (int l : X.w.leaves) a.push_back(X.LV.lev[l]);
    return a;
}

// Encoder (L -> r code bits) and decoder (r bits -> G) of a code; both empty if over budget.
// Each holds the min-AND candidate first and the min-depth one last. The caller frees them.
bool synthCode(const WinCtx& X, const Code& code, vector<Synth>& enc, vector<Synth>& dec) {
    int nL = (int)X.w.leaves.size(), nw = nWordsOf(nL), r = X.G.r;
    size_t P = (size_t)1 << nL;
    Tts on(r, vector<word>(nw, 0)), dc(r, vector<word>(nw, 0));
    for (size_t p = 0; p < P; p++) {
        int cw = code[X.G.inv[p]];
        for (int b = 0; b < r; b++) if ((cw >> b) & 1) on[b][p >> 6] |= (word)1 << (p & 63);
    }
    enc = synthesize(on, dc, nL, X.budget);
    if (enc.empty()) return false;
    vector<int> outs(X.w.outs.size());
    std::iota(outs.begin(), outs.end(), 0);
    Tts don, ddc;
    decoderTruths(X.G, code, r, outs, don, ddc);
    dec = synthesize(don, ddc, r, X.budget);
    if (dec.empty()) { freeSynth(enc); return false; }
    return true;
}

// end-to-end check: D(E(L)) (or P(L) when d is null) equals the window outputs on every leaf pattern
bool composesTo(const WinCtx& X, Gia_Man_t* e, Gia_Man_t* d) {
    int nL = (int)X.w.leaves.size(), nw = nWordsOf(nL);
    Tts y = simulate(e, leafPatterns(nL), nw);
    if (d) y = simulate(d, y, nw);
    return sameTruths(y, X.outTt, nL);
}

Cost costCode(WinCtx& X, const Code& code) {
    Cost c;
    vector<Synth> enc, dec;
    if (!synthCode(X, code, enc, dec)) return c;
    int nC = (int)X.w.cone.size();
    for (int pick = 0; pick < 2; pick++) {
        Synth& e = pick ? enc.back() : enc[0];
        Synth& d = pick ? dec.back() : dec[0];
        if (!composesTo(X, e.p, d.p)) {
            X.nErrors++;
            freeSynth(enc); freeSynth(dec);
            return c;
        }
        Pt& t = pick ? c.D : c.A;
        t.gain = nC - e.ands - d.ands;
        delayOf(X, outArrivals(d.p, outArrivals(e.p, leafArrivals(X))), t.dlev, t.ddep);
    }
    freeSynth(enc); freeSynth(dec);
    c.ok = true;
    return c;
}

Cost costDirect(WinCtx& X) {
    Cost c;
    int nL = (int)X.w.leaves.size(), nw = nWordsOf(nL);
    Tts dc(X.outTt.size(), vector<word>(nw, 0));
    vector<Synth> s = synthesize(X.outTt, dc, nL, X.budget);
    if (s.empty()) return c;
    int nC = (int)X.w.cone.size();
    for (int pick = 0; pick < 2; pick++) {
        Synth& e = pick ? s.back() : s[0];
        Pt& t = pick ? c.D : c.A;
        t.gain = nC - e.ands;
        delayOf(X, outArrivals(e.p, leafArrivals(X)), t.dlev, t.ddep);
    }
    freeSynth(s);
    c.ok = true;
    return c;
}

// ---------------------------------------------------------------------------------------------
// Per-sampler statistics

struct MStat {
    int samples = 0, fail = 0, nOk = 0, nPos = 0;
    std::set<Code> seen;
    int best = INT_MIN, bestDlev = 0, bestDdep = 0;
    double sumGain = 0;
    int minDlev = INT_MAX, gainAtMinDlev = 0;
    void add(const Code* canon, const Cost& c) {
        samples++;
        if (canon) seen.insert(*canon);
        if (!c.ok) { fail++; return; }
        nOk++;
        sumGain += c.A.gain;
        nPos += c.A.gain > 0;
        if (c.A.gain > best || (c.A.gain == best && c.A.dlev < bestDlev)) { best = c.A.gain; bestDlev = c.A.dlev; bestDdep = c.A.ddep; }
        if (c.D.dlev < minDlev || (c.D.dlev == minDlev && c.D.gain > gainAtMinDlev)) { minDlev = c.D.dlev; gainAtMinDlev = c.D.gain; }
    }
};

struct Options {
    int K = 8, nRoots = 1000, nTry = 16, xCap = 30, maxCone = 500, minN = 5, nWin = 1, coverPct = 0, budgetMul = 3, timeLimit = 0, localHops = 0;
    // AND budget of each synthesized circuit (E, D, or the direct L -> G): budgetMul * (|C| + 10); 0 = none
    int budget(int nC) const { return budgetMul ? budgetMul * (nC + 10) : INT_MAX; }
    unsigned seed = 1;
    bool analysis = false, kick = false, ils = false, sideways = false, keepCode = false, cec = true, fault = false, selfTest = false,
         verbose = false;
    string out, method = "U", script = "&put; dc2; &get -n", bestOut, kickDir;
};

int selfTest() {
    int fail = 0;
    for (int N = 3; N <= 8; N++) {
        int r = ceilLog2(N);
        vector<Code> all;
        enumerateCanonical(N, r, all);
        std::mt19937_64 rng(N);
        std::set<Code> hit;
        for (int i = 0; i < 20000; i++) hit.insert(canonical(sampleUniform(N, r, rng), r));
        bool ok = (u64)all.size() == numCanonical(N, r) && hit == std::set<Code>(all.begin(), all.end());
        printf("  N=%d r=%d: %zu canonical codes enumerated (formula %llu), uniform sampling reaches %zu: %s\n", N, r,
               all.size(), numCanonical(N, r), hit.size(), ok ? "ok" : "FAIL");
        fail |= !ok;
    }
    return fail;
}

int runAnalysis(Gia_Man_t* pGia, const Options& opt) {
    PipeAig g;
    string msg;
    if (!pipeAigFromGia(pGia, g, msg)) { Abc_Print(-1, "%s\n", msg.c_str()); return 1; }
    Levels LV = computeLevels(g);
    const char* name = pGia->pName ? pGia->pName : "circuit";
    FILE* f = opt.out.empty() ? stdout : fopen(opt.out.c_str(), "w");
    if (!f) { Abc_Print(-1, "Cannot open \"%s\" for writing.\n", opt.out.c_str()); return 1; }
    fprintf(f, "circuit,K,win,root,nL,nC,nG,N,r,trunc,shrunk,out_lev,method,samples,distinct,fail,best_gain,best_dlev,best_ddep,"
               "mean_gain,frac_pos,min_dlev,gain_at_min_dlev,secs\n");
    vector<int> roots = liveRoots(g, LV);
    std::mt19937_64 rng(opt.seed);
    std::shuffle(roots.begin(), roots.end(), rng);
    if ((int)roots.size() > opt.nRoots) roots.resize(opt.nRoots);
    Marks M(g.numNodes());
    std::set<vector<int>> seenLeaves;
    int nWin = 0, nDup = 0, nBad = 0, nCyclic = 0, nX = 0, nErrors = 0;
    clock_t t0 = clock();
    for (int root : roots) {
        RWin w;
        WinStatus ws = buildWindow(g, LV, root, opt.K, opt.maxCone, M, w);
        if (ws == WIN_NONE) { nBad++; continue; }
        if (ws == WIN_CYCLIC) { nCyclic++; continue; }
        if (!seenLeaves.insert(w.leaves).second) { nDup++; continue; }
        clock_t tw = clock();
        Group G = groupOf(g, w);
        Tts outTt = outputTruths(g, w);
        int outLev = 0;
        for (int v : w.outs) outLev = std::max(outLev, LV.lev[v]);
        char head[256];
        sprintf(head, "%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d", name, opt.K, nWin, root, (int)w.leaves.size(), (int)w.cone.size(),
                (int)w.outs.size(), G.nImg, G.r, (int)w.trunc, (int)w.shrunk, outLev);
        nWin++;
        if (G.nImg < 3) { fprintf(f, "%s,none,0,0,0,,,,,,,,0\n", head); continue; }
        WinCtx X{ g, LV, w, G, outTt, opt.budget((int)w.cone.size()) };
        std::map<Code, Cost> cache;
        auto cost = [&](const Code& code, Code& canon) {
            canon = canonical(code, G.r);
            auto it = cache.find(canon);
            if (it != cache.end()) return it->second;
            Cost c = costCode(X, canon);
            cache[canon] = c;
            return c;
        };
        std::mt19937_64 wr(opt.seed * 1000003ULL + root);
        Code base = baseCode(G);
        vector<std::pair<string, MStat>> stats;
        Code canon;
        { MStat s; s.add(nullptr, costDirect(X)); stats.push_back({ "P0", s }); }
        { MStat s; Cost c = cost(base, canon); s.add(&canon, c); stats.push_back({ "base", s }); }
        for (int m = S_U; m <= S_W4; m++) {
            MStat s;
            for (int t = 0; t < opt.nTry; t++) {
                Cost k = cost(sampleCode(m, G, base, wr), canon);
                s.add(&canon, k);
            }
            stats.push_back({ kMethods[m], s });
        }
        if (G.nImg <= 8 && nX < opt.xCap) {
            nX++;
            vector<Code> all;
            enumerateCanonical(G.nImg, G.r, all);
            MStat s;
            for (const Code& c : all) { Cost k = cost(c, canon); s.add(&canon, k); }
            stats.push_back({ "X", s });
        }
        nErrors += X.nErrors;
        double secs = (double)(clock() - tw) / CLOCKS_PER_SEC;
        for (auto& ms : stats) {
            const MStat& s = ms.second;
            if (s.nOk)
                fprintf(f, "%s,%s,%d,%zu,%d,%d,%d,%d,%.2f,%.3f,%d,%d,%.3f\n", head, ms.first.c_str(), s.samples, s.seen.size(),
                        s.fail, s.best, s.bestDlev, s.bestDdep, s.sumGain / s.nOk, (double)s.nPos / s.nOk, s.minDlev,
                        s.gainAtMinDlev, secs);
            else
                fprintf(f, "%s,%s,%d,%zu,%d,,,,,,,,%.3f\n", head, ms.first.c_str(), s.samples, s.seen.size(), s.fail, secs);
        }
        if (opt.verbose)
            printf("  win %4d root %6d: |L|=%2zu |C|=%3zu |G|=%3zu N=%4d r=%2d  P0 %+d  base %+d  U best %+d  (%.2fs)\n", nWin - 1,
                   root, w.leaves.size(), w.cone.size(), w.outs.size(), G.nImg, G.r, stats[0].second.best, stats[1].second.best,
                   stats[2].second.best, secs);
    }
    if (f != stdout) fclose(f);
    printf("lsv_reencode -a: %s K=%d: %d windows (%d duplicate, %d without a window, %d cyclic), %d with X, %d internal errors, %.1fs\n",
           name, opt.K, nWin, nDup, nBad, nCyclic, nX, nErrors, (double)(clock() - t0) / CLOCKS_PER_SEC);
    return nErrors ? 1 : 0;
}

// ---------------------------------------------------------------------------------------------
// Kick mode (phase 2): replace one window (or, with -m / -P, several windows in sequence) by a
// re-encoded implementation, then run the script R. Every kick starts from the same AIG x and is
// compared with R(x).

// Copy of a small GIA into p with its CIs driven by ins; returns its CO literals.
vector<int> instantiate(Gia_Man_t* p, Gia_Man_t* s, const vector<int>& ins) {
    vector<int> m(Gia_ManObjNum(s), 0);
    Gia_Obj_t* pObj;
    int i;
    Gia_ManForEachCi(s, pObj, i) m[Gia_ObjId(s, pObj)] = i < (int)ins.size() ? ins[i] : 0;
    Gia_ManForEachAnd(s, pObj, i)
        m[i] = Gia_ManHashAnd(p, Abc_LitNotCond(m[Gia_ObjFaninId0(pObj, i)], Gia_ObjFaninC0(pObj)),
                              Abc_LitNotCond(m[Gia_ObjFaninId1(pObj, i)], Gia_ObjFaninC1(pObj)));
    vector<int> out;
    Gia_ManForEachCo(s, pObj, i) out.push_back(Abc_LitNotCond(m[Gia_ObjFaninId0p(s, pObj)], Gia_ObjFaninC0(pObj)));
    return out;
}

Gia_Man_t* dupNamed(Gia_Man_t* p, Gia_Man_t* pNames) {
    Gia_Man_t* q = Gia_ManDup(p);
    copyNames(q, pNames);
    return q;
}

// Writes p as an AIGER file, with the names of pNames if p has none of its own.
void writeAig(Gia_Man_t* p, const string& fn, Gia_Man_t* pNames = nullptr) {
    Gia_Man_t* b = pNames ? dupNamed(p, pNames) : p;
    Gia_AigerWrite(b, (char*)fn.c_str(), 0, 0, 0);
    if (b != p) Gia_ManStop(b);
}

// Origin of every AND of a kicked AIG (written next to it by -W, for residue analysis): X = outside
// the transitive fanout of C (copied from x), E = encoder, D = decoder (or the direct L -> G circuit of
// P0), T = transitive fanout of C rebuilt on top of the replacement. leaves / code / gout: literals of
// L, of the code bits Y and of the new G in the kicked AIG (-1 if the node was swept away).
struct KickTrace {
    string lab;
    vector<int> leaves, code, gout;
};

// The AIG with C replaced: first maps L to Y and second maps Y to G (or first maps L to G when second
// is null). Nodes outside the transitive fanout of C are built first (every leaf is among them, see
// buildWindow), then the replacement, then the rest. Structural hashing is on, so the new logic may
// share nodes with the old. Returns nullptr if a removed node is still referenced (internal error).
Gia_Man_t* replaceWindow(const PipeAig& g, const RWin& w, Gia_Man_t* first, Gia_Man_t* second, Gia_Man_t* pNames,
                         bool fault = false, KickTrace* tr = nullptr) {
    int n = g.numNodes();
    vector<char> inC(n, 0), tfo(n, 0);
    for (int v : w.cone) inC[v] = 1;
    for (int v = g.nPi + 1; v < n; v++) tfo[v] = inC[v] || tfo[g.fanin0[v] >> 1] || tfo[g.fanin1[v] >> 1];
    Gia_Man_t* p = Gia_ManStart(n + (first ? 2 * Gia_ManObjNum(first) : 0) + (second ? 2 * Gia_ManObjNum(second) : 0) + 16);
    Gia_ManHashAlloc(p);
    vector<int> lit(n, -1);
    lit[0] = 0;
    for (int i = 0; i < g.nPi; i++) lit[i + 1] = Gia_ManAppendCi(p);
    bool bad = false;
    auto litOf = [&](int l) { if (lit[l >> 1] < 0) { bad = true; return 0; } return Abc_LitNotCond(lit[l >> 1], l & 1); };
    auto build = [&](int v) { lit[v] = Gia_ManHashAnd(p, litOf(g.fanin0[v]), litOf(g.fanin1[v])); };
    for (int v = g.nPi + 1; v < n; v++) if (!tfo[v]) build(v);
    vector<int> in;
    for (int l : w.leaves) in.push_back(litOf(2 * l));
    int n0 = Gia_ManObjNum(p);
    vector<int> mid, out;
    if (first) {
        mid = instantiate(p, first, in);
    } else {   // null kick: C itself, i.e. the same AIG with the nodes renumbered as for a real kick
        for (int v : w.cone) build(v);
        for (int v : w.outs) mid.push_back(lit[v]);
    }
    int n1 = Gia_ManObjNum(p);
    out = second ? instantiate(p, second, mid) : mid;
    int n2 = Gia_ManObjNum(p);
    if (fault) out[0] = Abc_LitNot(out[0]);   // test hook (-F): a wrong replacement that cec must catch
    for (size_t j = 0; j < w.outs.size(); j++) lit[w.outs[j]] = out[j];
    for (int v = g.nPi + 1; v < n; v++) if (tfo[v] && !inC[v]) build(v);
    for (int l : g.poLit) Gia_ManAppendCo(p, litOf(l));
    Gia_ManHashStop(p);
    Gia_Man_t* q = Gia_ManCleanup(p);   // sets Value of every kept object of p to its literal in q
    if (tr && !bad) {
        auto mapLit = [&](int l) { unsigned v = Gia_ManObj(p, Abc_Lit2Var(l))->Value; return v == ~0u ? -1 : Abc_LitNotCond((int)v, Abc_LitIsCompl(l)); };
        tr->lab.assign(Gia_ManObjNum(q), '-');
        Gia_Obj_t* pObj;
        int i;
        Gia_ManForEachAnd(p, pObj, i) {
            if (pObj->Value == ~0u) continue;
            tr->lab[Abc_Lit2Var(pObj->Value)] = i < n0 ? 'X' : i < n1 ? (second ? 'E' : 'D') : i < n2 ? 'D' : 'T';
        }
        for (int l : in) tr->leaves.push_back(mapLit(l));
        if (second) for (int l : mid) tr->code.push_back(mapLit(l));
        for (int l : out) tr->gout.push_back(mapLit(l));
    }
    Gia_ManStop(p);
    if (bad) { Gia_ManStop(q); return nullptr; }
    copyNames(q, pNames);
    return q;
}

// Runs script with a copy of p as the current GIA; returns a copy of the resulting GIA, or nullptr.
Gia_Man_t* runScript(Abc_Frame_t* pAbc, Gia_Man_t* p, const string& script) {
    Abc_FrameUpdateGia(pAbc, Gia_ManDup(p));
    vector<char> cmd(script.begin(), script.end());
    cmd.push_back(0);
    if (Cmd_CommandExecute(pAbc, cmd.data())) return nullptr;
    Gia_Man_t* q = Abc_FrameReadGia(pAbc);
    return q ? Gia_ManDup(q) : nullptr;
}

// R on a region of y only (-L hops, ILS): the region is every AND within hops fanin / fanout steps of the
// replacement (the nodes labelled E or D by the kick trace), made convex (a node that depends on the region
// through a node outside it is dropped, so the optimized region can be put back as one block). The
// region becomes a separate AIG (inputs: its fanins from outside, outputs: its nodes used outside or by a
// PO), R runs on it, and the result (or the region as it was, if R did not make it smaller) replaces it;
// the rest of y is copied. Returns nullptr if R failed.
Gia_Man_t* localR(Abc_Frame_t* pAbc, Gia_Man_t* y, const string& lab, int hops, const string& script, Gia_Man_t* pNames,
                  int& regionSize) {
    PipeAig g;
    string msg;
    if (!pipeAigFromGia(y, g, msg)) return nullptr;
    int n = g.numNodes();
    vector<char> inR(n, 0);
    vector<int> front;
    for (int v = g.nPi + 1; v < n && v < (int)lab.size(); v++)
        if (lab[v] == 'E' || lab[v] == 'D') { inR[v] = 1; front.push_back(v); }
    for (int h = 0; h < hops && !front.empty(); h++) {
        vector<int> next;
        auto add = [&](int u) { if (g.isAnd(u) && !inR[u]) { inR[u] = 1; next.push_back(u); } };
        for (int v : front) {
            add(g.fanin0[v] >> 1);
            add(g.fanin1[v] >> 1);
            for (int k = g.foStart[v]; k < g.foStart[v + 1]; k++) add(g.foList[k]);
        }
        front.swap(next);
    }
    // convexity: a region node with a fanin outside the region that depends on the region would close a
    // loop when the region is put back as one block, so it leaves the region (in topological order, so
    // the nodes depending on it leave too). The region only shrinks: it stays within the hops neighbourhood.
    vector<char> down(n, 0);   // outside the region, in its transitive fanout
    for (int v = g.nPi + 1; v < n; v++) {
        int a = g.fanin0[v] >> 1, b = g.fanin1[v] >> 1;
        if (inR[v] && ((!inR[a] && down[a]) || (!inR[b] && down[b]))) inR[v] = 0;
        down[v] = !inR[v] && (inR[a] || down[a] || inR[b] || down[b]);
    }
    // the region as an AIG
    vector<char> isOut(n, 0);
    for (int l : g.poLit) if (inR[l >> 1]) isOut[l >> 1] = 1;
    for (int v = g.nPi + 1; v < n; v++) {
        if (inR[v]) continue;
        if (inR[g.fanin0[v] >> 1]) isOut[g.fanin0[v] >> 1] = 1;
        if (inR[g.fanin1[v] >> 1]) isOut[g.fanin1[v] >> 1] = 1;
    }
    RWin rw;   // the region as a window: leaves = its fanins from outside, cone = its nodes, outs = isOut
    vector<int> lit(n, -1);
    Gia_Man_t* s = Gia_ManStart(2 * n + 16);
    Gia_ManHashAlloc(s);
    for (int v = g.nPi + 1; v < n; v++) {
        if (!inR[v]) continue;
        rw.cone.push_back(v);
        for (int f : { g.fanin0[v] >> 1, g.fanin1[v] >> 1 })
            if (!inR[f] && lit[f] < 0) { lit[f] = f == 0 ? 0 : Gia_ManAppendCi(s); if (f) rw.leaves.push_back(f); }
        lit[v] = Gia_ManHashAnd(s, Abc_LitNotCond(lit[g.fanin0[v] >> 1], g.fanin0[v] & 1), Abc_LitNotCond(lit[g.fanin1[v] >> 1], g.fanin1[v] & 1));
        if (isOut[v]) rw.outs.push_back(v);
    }
    regionSize = (int)rw.cone.size();
    for (int v : rw.outs) Gia_ManAppendCo(s, lit[v]);
    Gia_ManHashStop(s);
    Gia_Man_t* s0 = Gia_ManCleanup(s);
    Gia_ManStop(s);
    Gia_Man_t* z = runScript(pAbc, s0, script);
    if (!z) { Gia_ManStop(s0); return nullptr; }
    // y with the region replaced by blk (the region is convex, so its leaves are built before it)
    Gia_Man_t* q = replaceWindow(g, rw, Gia_ManAndNum(z) < Gia_ManAndNum(s0) ? z : s0, nullptr, pNames);
    Gia_ManStop(z);
    Gia_ManStop(s0);
    return q;
}

// Simulation signatures of every object of p: nW random words per CI (fixed seed, by CI index, so AIGs
// with the same CIs get the same patterns), canonical up to complement (bit 0 = 0).
vector<vector<word>> simSigs(Gia_Man_t* p, int nW = 32) {
    vector<vector<word>> s(Gia_ManObjNum(p), vector<word>(nW, 0));
    std::mt19937_64 r(12345);
    Gia_Obj_t* o;
    int i;
    Gia_ManForEachCi(p, o, i) for (auto& w : s[Gia_ObjId(p, o)]) w = r();
    Gia_ManForEachAnd(p, o, i) {
        const vector<word>& a = s[Gia_ObjFaninId0(o, i)];
        const vector<word>& b = s[Gia_ObjFaninId1(o, i)];
        word ca = Gia_ObjFaninC0(o) ? ~(word)0 : 0, cb = Gia_ObjFaninC1(o) ? ~(word)0 : 0;
        for (int k = 0; k < nW; k++) s[i][k] = (a[k] ^ ca) & (b[k] ^ cb);
    }
    for (auto& v : s) if (v[0] & 1) for (auto& w : v) w = ~w;
    return s;
}

std::set<vector<word>> sigSet(Gia_Man_t* p) {
    vector<vector<word>> s = simSigs(p);
    return std::set<vector<word>>(s.begin(), s.end());
}

// Code bits of a traced kick y (lits = K.tr.code) after R (z): bits = distinct non-constant code
// functions, fresh = those that are not a function of x (Fx = sigSet(x)), surv / freshSurv = those still
// computed by some node of z.
struct CodeSurv { int bits = 0, fresh = 0, surv = 0, freshSurv = 0; };

CodeSurv codeSurvival(Gia_Man_t* y, const vector<int>& lits, const std::set<vector<word>>& Fx, Gia_Man_t* z) {
    CodeSurv c;
    if (lits.empty()) return c;
    vector<vector<word>> sy = simSigs(y);
    std::set<vector<word>> Fz = sigSet(z), seen;
    vector<word> zero(sy[0].size(), 0);
    for (int l : lits) {
        if (l < 0) continue;
        const vector<word>& f = sy[Abc_Lit2Var(l)];
        if (f == zero || !seen.insert(f).second) continue;
        bool fresh = !Fx.count(f), surv = Fz.count(f) > 0;
        c.bits++; c.fresh += fresh; c.surv += surv; c.freshSurv += fresh && surv;
    }
    return c;
}

// R with the code bits kept (-b): every code node of y (an AND; code bits that are CIs or constants are
// left alone) is cut: its uses read a new CI and a new CO exports it, so R cannot optimize across the code
// bits; afterwards each new CI is reconnected to the optimized code bit. nCut = number of cut nodes.
// Returns nullptr if R failed or the reconnection is cyclic (an optimized code bit depending on itself).
Gia_Man_t* runScriptCut(Abc_Frame_t* pAbc, Gia_Man_t* y, const vector<int>& code, const string& script, Gia_Man_t* pNames,
                        int& nCut) {
    int n = Gia_ManObjNum(y), nPo = Gia_ManCoNum(y), nPi = Gia_ManCiNum(y);
    vector<int> cutIdx(n, -1), cuts;
    for (int l : code) {
        if (l < 0) continue;
        int v = Abc_Lit2Var(l);
        if (!Gia_ObjIsAnd(Gia_ManObj(y, v)) || cutIdx[v] >= 0) continue;
        cutIdx[v] = (int)cuts.size();
        cuts.push_back(v);
    }
    nCut = (int)cuts.size();
    Gia_Man_t* p = Gia_ManStart(n + 2 * nCut + 16);
    vector<int> lit(n, -1), pc(nCut);
    lit[0] = 0;
    Gia_Obj_t* o;
    int i;
    Gia_ManForEachCi(y, o, i) lit[Gia_ObjId(y, o)] = Gia_ManAppendCi(p);
    for (int j = 0; j < nCut; j++) pc[j] = Gia_ManAppendCi(p);
    auto use = [&](int v, int c) { return Abc_LitNotCond(cutIdx[v] >= 0 ? pc[cutIdx[v]] : lit[v], c); };
    Gia_ManHashAlloc(p);
    Gia_ManForEachAnd(y, o, i) lit[i] = Gia_ManHashAnd(p, use(Gia_ObjFaninId0(o, i), Gia_ObjFaninC0(o)), use(Gia_ObjFaninId1(o, i), Gia_ObjFaninC1(o)));
    Gia_ManForEachCo(y, o, i) Gia_ManAppendCo(p, use(Gia_ObjFaninId0p(y, o), Gia_ObjFaninC0(o)));
    for (int v : cuts) Gia_ManAppendCo(p, lit[v]);
    Gia_ManHashStop(p);
    Gia_Man_t* p0 = Gia_ManCleanup(p);
    Gia_ManStop(p);
    Gia_Man_t* zc = runScript(pAbc, p0, script);
    Gia_ManStop(p0);
    if (!zc) return nullptr;
    // reconnect: resolve the cut CIs in passes (a code bit may be computed from another one)
    int m = Gia_ManObjNum(zc);
    vector<int> zl(m, -1);
    zl[0] = 0;
    Gia_Man_t* q = Gia_ManStart(m + 16);
    Gia_ManHashAlloc(q);
    for (int j = 0; j < nPi; j++) zl[Gia_ObjId(zc, Gia_ManCi(zc, j))] = Gia_ManAppendCi(q);
    auto zlit = [&](int v, int c) { return zl[v] < 0 ? -1 : Abc_LitNotCond(zl[v], c); };
    for (int left = nCut, pass = 0; ; pass++) {
        Gia_ManForEachAnd(zc, o, i) {
            if (zl[i] >= 0) continue;
            int a = zlit(Gia_ObjFaninId0(o, i), Gia_ObjFaninC0(o)), b = zlit(Gia_ObjFaninId1(o, i), Gia_ObjFaninC1(o));
            if (a >= 0 && b >= 0) zl[i] = Gia_ManHashAnd(q, a, b);
        }
        if (!left) break;
        int before = left;
        for (int j = 0; j < nCut; j++) {
            int ci = Gia_ObjId(zc, Gia_ManCi(zc, nPi + j));
            if (zl[ci] >= 0) continue;
            Gia_Obj_t* co = Gia_ManCo(zc, nPo + j);
            int d = zlit(Gia_ObjFaninId0p(zc, co), Gia_ObjFaninC0(co));
            if (d >= 0) { zl[ci] = d; left--; }
        }
        if (left == before) { Gia_ManHashStop(q); Gia_ManStop(q); Gia_ManStop(zc); return nullptr; }
    }
    for (int j = 0; j < nPo; j++) {
        Gia_Obj_t* co = Gia_ManCo(zc, j);
        Gia_ManAppendCo(q, zlit(Gia_ObjFaninId0p(zc, co), Gia_ObjFaninC0(co)));
    }
    Gia_ManHashStop(q);
    Gia_ManStop(zc);
    Gia_Man_t* r = Gia_ManCleanup(q);
    Gia_ManStop(q);
    copyNames(r, pNames);
    return r;
}

// One random window of g with at least max(3, minN) classes, re-encoded with sampler m (s1 = E or the
// direct L -> G circuit for P0, s2 = D). Returns 1 on success, 0 if no window was found in 200 tries,
// -1 on an internal error.
int pickWindow(const PipeAig& g, const Levels& LV, Marks& M, int m, const Options& opt, std::mt19937_64& rng,
               RWin& w, Group& G, vector<Synth>& s1, vector<Synth>& s2, int& tries) {
    vector<int> roots = liveRoots(g, LV);
    if (roots.empty()) return 0;
    for (tries = 1; tries <= 200; tries++) {
        freeSynth(s1); freeSynth(s2);
        int root = roots[rng() % roots.size()];
        if (buildWindow(g, LV, root, opt.K, opt.maxCone, M, w) != WIN_OK) continue;
        G = groupOf(g, w);
        if (G.nImg < std::max(3, opt.minN)) continue;
        Tts outTt = outputTruths(g, w);
        WinCtx X{ g, LV, w, G, outTt, opt.budget((int)w.cone.size()) };
        if (m == S_P0) {
            Tts dc(outTt.size(), vector<word>(nWordsOf((int)w.leaves.size()), 0));
            s1 = synthesize(outTt, dc, (int)w.leaves.size(), X.budget);
            if (s1.empty()) continue;
            if (!composesTo(X, s1[0].p, nullptr)) return -1;
        } else if (m != S_NULL) {   // the null kick needs no synthesis
            Code code = sampleCode(m, G, baseCode(G), rng);
            if (!synthCode(X, canonical(code, G.r), s1, s2)) continue;
            if (!composesTo(X, s1[0].p, s2[0].p)) return -1;
        }
        return 1;
    }
    freeSynth(s1); freeSynth(s2);
    return 0;
}

int methodIndex(const Options& opt) {
    for (int i = S_U; i <= S_NULL; i++) if (opt.method == kMethods[i]) return i;
    Abc_Print(-1, "Unknown sampler \"%s\" (use U, A, W1, W4, base, P0 or N).\n", opt.method.c_str());
    return -1;
}

// One kick of x: windows re-encoded one after another, each picked on the AIG left by the previous one,
// until -m windows were replaced or (-P) the replaced cones add up to P% of AND(x). With one window
// per kick, that window is traced. Returns 1 with the kicked AIG in y, 0 if no
// window was found, -1 on an internal error (y = nullptr unless 1).
struct Kick {
    RWin w0;
    int nImg0 = 0, r0 = 0;   // first window: classes, code bits
    int tries0 = 0, nWin = 0, cov = 0;
    int rootLev = 0, rootDown = 0, depth = 0;   // first window: level of the root, its distance to a PO, AIG depth
    KickTrace tr;
};

int makeKick(Gia_Man_t* x, int m, const Options& opt, std::mt19937_64& rng, Gia_Man_t*& y, Kick& K) {
    int aX = Gia_ManAndNum(x);
    y = dupNamed(x, x);
    while (opt.coverPct ? 100LL * K.cov < (long long)opt.coverPct * aX : K.nWin < opt.nWin) {
        PipeAig g;
        string msg;
        if (!pipeAigFromGia(y, g, msg)) { Abc_Print(-1, "%s\n", msg.c_str()); Gia_ManStop(y); y = nullptr; return -1; }
        Levels LV = computeLevels(g);
        Marks M(g.numNodes());
        RWin w;
        Group G;
        vector<Synth> s1, s2;
        int tries = 0, st = pickWindow(g, LV, M, m, opt, rng, w, G, s1, s2, tries);
        if (st < 0) { Gia_ManStop(y); y = nullptr; return -1; }
        if (st == 0) break;
        bool tw = K.nWin == 0 && opt.nWin == 1 && !opt.coverPct;
        Gia_Man_t* y2 = replaceWindow(g, w, s1.empty() ? nullptr : s1[0].p, s2.empty() ? nullptr : s2[0].p, x, opt.fault && K.nWin == 0,
                                      tw ? &K.tr : nullptr);
        freeSynth(s1); freeSynth(s2);
        if (!y2) { Gia_ManStop(y); y = nullptr; return -1; }
        Gia_ManStop(y);
        y = y2;
        if (K.nWin == 0) { K.w0 = w; K.nImg0 = G.nImg; K.r0 = G.r; K.tries0 = tries; K.rootLev = LV.lev[w.root]; K.rootDown = LV.down[w.root]; K.depth = LV.depth; }
        K.nWin++;
        K.cov += (int)w.cone.size();
    }
    if (K.nWin == 0) { Gia_ManStop(y); y = nullptr; return 0; }
    return 1;
}

int runKicks(Abc_Frame_t* pAbc, Gia_Man_t* pCur, const Options& opt) {
    int m = methodIndex(opt);
    if (m < 0) return 1;
    Gia_Man_t* x = dupNamed(pCur, pCur);   // own copy: the current GIA is replaced while R runs
    string name = x->pName ? x->pName : "circuit";
    int aX = Gia_ManAndNum(x), lX = Gia_ManLevelNum(x);
    Gia_Man_t* Rx = runScript(pAbc, x, opt.script);
    if (!Rx) { Abc_Print(-1, "The script \"%s\" failed.\n", opt.script.c_str()); Abc_FrameUpdateGia(pAbc, x); return 1; }
    int aRx = Gia_ManAndNum(Rx), lRx = Gia_ManLevelNum(Rx);
    if (!opt.kickDir.empty()) {
        writeAig(Rx, opt.kickDir + "/Rx.aig", x);
        writeAig(x, opt.kickDir + "/x.aig");   // x as numbered by the window ids below (CIs, then ANDs)
    }
    FILE* f = opt.out.empty() ? stdout : fopen(opt.out.c_str(), "w");
    if (!f) { Abc_Print(-1, "Cannot open \"%s\" for writing.\n", opt.out.c_str()); Abc_FrameUpdateGia(pAbc, x); Gia_ManStop(Rx); return 1; }
    fprintf(f, "circuit,K,method,kick,root,nL,nC,nG,N,r,shrunk,tries,and_x,lev_x,and_kick,lev_kick,and_R,lev_R,"
               "and_Rx,lev_Rx,d_R,d_x,undo,cec,secs,nwin,cov,code_bits,code_fresh,code_surv,code_fresh_surv,ncut\n");
    std::mt19937_64 rng(opt.seed * 7919ULL + m);
    std::set<vector<word>> Fx = sigSet(x);   // functions of x, for the code-bit survival columns
    Gia_Man_t* best = nullptr;
    int bestA = aRx, nKicks = 0, nUndo = 0, nBetter = 0, nErr = 0, nNoWin = 0, nCecFail = 0, nUndec = 0;
    clock_t t0 = clock();
    for (int k = 0; k < opt.nRoots && !nCecFail && !nErr; k++) {
        clock_t tk = clock();
        Kick K;
        Gia_Man_t* y = nullptr;
        int st = makeKick(x, m, opt, rng, y, K);
        if (st < 0) { nErr++; break; }
        if (st == 0) { nNoWin++; continue; }
        const KickTrace& tr = K.tr;
        const RWin& w0 = K.w0;
        int cec = opt.cec ? cecGiaSafe(y, x) : -2;
        nUndec += cec == -1;
        if (cec == 0) {
            nCecFail++;
            printf("lsv_reencode -k: kick %d (root %d) is NOT equivalent to x%s\n", k, w0.root, opt.fault ? " (fault injected by -F)" : "");
            Gia_ManStop(y);
            break;
        }
        int nCut = 0;
        Gia_Man_t* z = opt.keepCode && !tr.code.empty() ? runScriptCut(pAbc, y, tr.code, opt.script, x, nCut) : runScript(pAbc, y, opt.script);
        if (!z) { Gia_ManStop(y); Abc_Print(-1, "The script \"%s\" failed.\n", opt.script.c_str()); nErr++; break; }
        if (opt.keepCode && cecGiaSafe(z, x) == 0) {   // the reconnected AIG is new code: check it too
            nCecFail++;
            printf("lsv_reencode -k: kick %d (root %d): R with the code bits kept is NOT equivalent to x\n", k, w0.root);
            Gia_ManStop(y); Gia_ManStop(z);
            break;
        }
        CodeSurv cs = codeSurvival(y, tr.code, Fx, z);
        int aY = Gia_ManAndNum(y), lY = Gia_ManLevelNum(y), aZ = Gia_ManAndNum(z), lZ = Gia_ManLevelNum(z);
        bool undo = aZ == aRx && lZ == lRx;
        nKicks++; nUndo += undo; nBetter += aZ < aRx;
        if (!opt.kickDir.empty() && !tr.lab.empty()) {   // the kicked AIG y and the origin of its nodes
            writeAig(y, opt.kickDir + "/y" + std::to_string(k) + ".aig");
            string fn = opt.kickDir + "/y" + std::to_string(k) + ".lab";
            FILE* fl = fopen(fn.c_str(), "w");
            if (fl) {
                auto list = [&](const char* key, const vector<int>& v) {
                    fprintf(fl, "%s", key);
                    for (int x : v) fprintf(fl, " %d", x);
                    fprintf(fl, "\n");
                };
                list("xleaves", w0.leaves);
                list("xcone", w0.cone);
                list("xouts", w0.outs);
                list("leaves", tr.leaves);
                list("code", tr.code);
                list("gout", tr.gout);
                fprintf(fl, "lab %s\n", tr.lab.c_str());
                fclose(fl);
            }
        }
        if (!opt.kickDir.empty()) writeAig(z, opt.kickDir + "/k" + std::to_string(k) + ".aig", x);
        if (aZ < bestA) { if (best) Gia_ManStop(best); best = z; bestA = aZ; }
        else Gia_ManStop(z);
        Gia_ManStop(y);
        fprintf(f, "%s,%d,%s,%d,%d,%zu,%zu,%zu,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%.3f,%d,%d,%d,%d,%d,%d,%d\n", name.c_str(),
                opt.K, kMethods[m], k, w0.root, w0.leaves.size(), w0.cone.size(), w0.outs.size(), K.nImg0, K.r0,
                (int)w0.shrunk, K.tries0, aX, lX, aY, lY, aZ, lZ, aRx, lRx, aZ - aRx, aZ - aX, (int)undo, cec,
                (double)(clock() - tk) / CLOCKS_PER_SEC, K.nWin, K.cov, cs.bits, cs.fresh, cs.surv, cs.freshSurv, nCut);
        fflush(f);
    }
    if (f != stdout) fclose(f);
    int bestCec = -2;
    if (best) {
        bestCec = cecGiaSafe(best, x);
        if (!opt.bestOut.empty() && bestCec == 1) writeAig(best, opt.bestOut, x);
    }
    printf("lsv_reencode -k: %s K=%d M=%s: x %d/%d, R(x) %d/%d (and/lev); %d kicks, undo %d, better %d, best %+d (cec %s), "
           "no window %d, cec undecided %d, cec failures %d, internal errors %d, %.1fs\n",
           name.c_str(), opt.K, kMethods[m], aX, lX, aRx, lRx, nKicks, nUndo, nBetter, bestA - aRx,
           bestCec == 1 ? "ok" : bestCec == -2 ? "-" : bestCec == -1 ? "undecided" : "FAIL", nNoWin, nUndec, nCecFail, nErr, (double)(clock() - t0) / CLOCKS_PER_SEC);
    Abc_FrameUpdateGia(pAbc, x);   // the frame owns x again
    if (best) Gia_ManStop(best);
    Gia_ManStop(Rx);
    return (nCecFail || nErr || bestCec == 0) ? 1 : 0;
}

// Iterated local search (-i): cur = R(x) (or x, if R makes it worse); every iteration kicks cur (as -k does), runs R and accepts the
// result if it has fewer ANDs, or as many and fewer levels (-e: also as many ANDs and no more levels).
// Stops after -n iterations or -T seconds (wall clock). One CSV row per iteration; -W writes every
// accepted AIG as dir/a<iter>.aig; -w writes the final AIG if it passes cec against x.
int runIls(Abc_Frame_t* pAbc, Gia_Man_t* pCur, const Options& opt) {
    int m = methodIndex(opt);
    if (m < 0) return 1;
    abctime t0 = Abc_Clock();
    auto elapsed = [&]() { return (double)(Abc_Clock() - t0) / CLOCKS_PER_SEC; };
    Gia_Man_t* x = dupNamed(pCur, pCur);   // own copy: the current GIA is replaced while R runs
    string name = x->pName ? x->pName : "circuit";
    int aX = Gia_ManAndNum(x), lX = Gia_ManLevelNum(x);
    Gia_Man_t* cur = runScript(pAbc, x, opt.script);
    if (!cur) { Abc_Print(-1, "The script \"%s\" failed.\n", opt.script.c_str()); Abc_FrameUpdateGia(pAbc, x); return 1; }
    int aRx = Gia_ManAndNum(cur), lRx = Gia_ManLevelNum(cur);
    if (aX < aRx || (aX == aRx && lX < lRx)) { Gia_ManStop(cur); cur = dupNamed(x, x); }   // R made x worse: start from x
    int aCur = Gia_ManAndNum(cur), lCur = Gia_ManLevelNum(cur);
    FILE* f = opt.out.empty() ? stdout : fopen(opt.out.c_str(), "w");
    if (!f) { Abc_Print(-1, "Cannot open \"%s\" for writing.\n", opt.out.c_str()); Gia_ManStop(cur); Abc_FrameUpdateGia(pAbc, x); return 1; }
    fprintf(f, "circuit,K,method,iter,root,nL,nC,nG,N,r,nwin,cov,and_cur,lev_cur,and_kick,and_R,lev_R,d,acc,cec,secs,elapsed,"
               "lev_root,down_root,depth,region,cec_acc,code_bits,code_fresh,code_surv,code_fresh_surv,ncut\n");
    fprintf(f, "%s,%d,%s,-1,-1,0,0,0,0,0,0,0,%d,%d,%d,%d,%d,0,1,-2,0,%.2f,0,0,%d,0,-2,0,0,0,0,0\n", name.c_str(), opt.K, kMethods[m], aX, lX, aX, aCur, lCur,
            elapsed(), lCur);
    std::mt19937_64 rng(opt.seed * 7919ULL + m);
    std::set<vector<word>> Fcur;   // functions of cur, for the code-bit survival columns (built when needed)
    bool FcurOk = false;
    int nIter = 0, nAcc = 0, nErr = 0, nNoWin = 0, nCecFail = 0, nUndec = 0;
    for (int k = 0; k < opt.nRoots && !nCecFail && !nErr && (!opt.timeLimit || elapsed() < opt.timeLimit); k++) {
        abctime tk = Abc_Clock();
        Kick K;
        Gia_Man_t* y = nullptr;
        int st = makeKick(cur, m, opt, rng, y, K);
        if (st < 0) { nErr++; break; }
        if (st == 0) { nNoWin++; continue; }
        int cec = opt.cec ? cecGiaSafe(y, cur) : -2;
        nUndec += cec == -1;
        if (cec == 0) {
            nCecFail++;
            printf("lsv_reencode -i: iteration %d (root %d) is NOT equivalent to the current AIG%s\n", k, K.w0.root, opt.fault ? " (fault injected by -F)" : "");
            Gia_ManStop(y);
            break;
        }
        int region = 0;   // -L: R on the region around the replacement only (whole AIG if the kick was not traced)
        int nCut = 0;   // -b: R with the code bits kept
        Gia_Man_t* z = opt.localHops > 0 && !K.tr.lab.empty() ? localR(pAbc, y, K.tr.lab, opt.localHops, opt.script, x, region)
                     : opt.keepCode && !K.tr.code.empty() ? runScriptCut(pAbc, y, K.tr.code, opt.script, x, nCut)
                                                          : runScript(pAbc, y, opt.script);
        if (!z) { Gia_ManStop(y); Abc_Print(-1, "The script \"%s\" failed.\n", opt.script.c_str()); nErr++; break; }
        if (!K.tr.code.empty() && !FcurOk) { Fcur = sigSet(cur); FcurOk = true; }
        CodeSurv cs = codeSurvival(y, K.tr.code, Fcur, z);
        int aY = Gia_ManAndNum(y), aZ = Gia_ManAndNum(z), lZ = Gia_ManLevelNum(z);
        Gia_ManStop(y);
        bool acc = aZ < aCur || (aZ == aCur && (lZ < lCur || (opt.sideways && lZ <= lCur)));
        int cecAcc = -2;   // without per-iteration cec (-c) or with -L, every accepted result is checked against cur
        if (acc && (!opt.cec || opt.localHops > 0 || opt.keepCode)) {
            cecAcc = cecGiaSafe(z, cur);
            nUndec += cecAcc == -1;
            if (cecAcc == 0) {
                nCecFail++;
                printf("lsv_reencode -i: iteration %d (root %d): the accepted AIG is NOT equivalent to the current AIG\n", k, K.w0.root);
                Gia_ManStop(z);
                break;
            }
        }
        nIter++;
        fprintf(f, "%s,%d,%s,%d,%d,%zu,%zu,%zu,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%.3f,%.2f,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n", name.c_str(), opt.K, kMethods[m], k,
                K.w0.root, K.w0.leaves.size(), K.w0.cone.size(), K.w0.outs.size(), K.nImg0, K.r0, K.nWin, K.cov,
                aCur, lCur, aY, aZ, lZ, aZ - aCur, (int)acc, cec, (double)(Abc_Clock() - tk) / CLOCKS_PER_SEC, elapsed(),
                K.rootLev, K.rootDown, K.depth, region, cecAcc, cs.bits, cs.fresh, cs.surv, cs.freshSurv, nCut);
        fflush(f);
        if (!acc) { Gia_ManStop(z); continue; }
        nAcc++;
        Gia_ManStop(cur);
        cur = dupNamed(z, x);
        Gia_ManStop(z);
        aCur = aZ; lCur = lZ;
        FcurOk = false;
        if (!opt.kickDir.empty()) writeAig(cur, opt.kickDir + "/a" + std::to_string(k) + ".aig");
    }
    if (f != stdout) fclose(f);
    int fin = cecGiaSafe(cur, x);
    if (!opt.bestOut.empty() && fin == 1) writeAig(cur, opt.bestOut);
    printf("lsv_reencode -i: %s K=%d M=%s: x %d/%d, R(x) %d/%d, final %d/%d (and/lev, %+d vs R(x), cec %s); %d iterations, "
           "accepted %d, no window %d, cec undecided %d, cec failures %d, internal errors %d, %.1fs\n",
           name.c_str(), opt.K, kMethods[m], aX, lX, aRx, lRx, aCur, lCur, aCur - aRx,
           fin == 1 ? "ok" : fin == -1 ? "undecided" : "FAIL", nIter, nAcc, nNoWin, nUndec, nCecFail, nErr, elapsed());
    Abc_FrameUpdateGia(pAbc, x);   // the frame owns x again
    Gia_ManStop(cur);
    return (nCecFail || nErr || fin == 0) ? 1 : 0;
}

// ---------------------------------------------------------------------------------------------
// Command

int Lsv_CommandReencode(Abc_Frame_t* pAbc, int argc, char** argv) {
    Options opt;
    int c, seed = 1;
    Extra_UtilGetoptReset();
    while ((c = Extra_UtilGetopt(argc, argv, "KntXCSNmPBTLoMRwWakiebcFsvh")) != EOF) {
        switch (c) {
        case 'K': if (!parseInt(argc, argv, opt.K, 2, "K")) goto usage; if (opt.K > 14) { Abc_Print(-1, "-K must be <= 14.\n"); goto usage; } break;
        case 'n': if (!parseInt(argc, argv, opt.nRoots, 1, "n")) goto usage; break;
        case 't': if (!parseInt(argc, argv, opt.nTry, 1, "t")) goto usage; break;
        case 'X': if (!parseInt(argc, argv, opt.xCap, 0, "X")) goto usage; break;
        case 'C': if (!parseInt(argc, argv, opt.maxCone, 1, "C")) goto usage; break;
        case 'S': if (!parseInt(argc, argv, seed, 0, "S")) goto usage; break;
        case 'N': if (!parseInt(argc, argv, opt.minN, 3, "N")) goto usage; break;
        case 'm': if (!parseInt(argc, argv, opt.nWin, 1, "m")) goto usage; break;
        case 'P': if (!parseInt(argc, argv, opt.coverPct, 0, "P")) goto usage; break;
        case 'B': if (!parseInt(argc, argv, opt.budgetMul, 0, "B")) goto usage; break;
        case 'T': if (!parseInt(argc, argv, opt.timeLimit, 0, "T")) goto usage; break;
        case 'L': if (!parseInt(argc, argv, opt.localHops, 0, "L")) goto usage; break;
        case 'o':
            if (globalUtilOptind >= argc) { Abc_Print(-1, "Command line switch \"-o\" should be followed by a file name.\n"); goto usage; }
            opt.out = argv[globalUtilOptind++];
            break;
        case 'M': case 'R': case 'w': case 'W':
            if (globalUtilOptind >= argc) { Abc_Print(-1, "Command line switch \"-%c\" should be followed by a string.\n", c); goto usage; }
            (c == 'M' ? opt.method : c == 'R' ? opt.script : c == 'w' ? opt.bestOut : opt.kickDir) = argv[globalUtilOptind++];
            break;
        case 'a': opt.analysis = !opt.analysis; break;
        case 'k': opt.kick = !opt.kick; break;
        case 'i': opt.ils = !opt.ils; break;
        case 'e': opt.sideways = !opt.sideways; break;
        case 'b': opt.keepCode = !opt.keepCode; break;
        case 'c': opt.cec = !opt.cec; break;
        case 'F': opt.fault = !opt.fault; break;
        case 's': opt.selfTest = !opt.selfTest; break;
        case 'v': opt.verbose = !opt.verbose; break;
        default: goto usage;
        }
    }
    opt.seed = (unsigned)seed;
    if (opt.selfTest) return selfTest();
    if (opt.analysis + opt.kick + opt.ils != 1) { Abc_Print(-1, "Choose exactly one of -a (analysis), -k (kicks) and -i (ILS).\n"); return 1; }
    if (opt.keepCode && (opt.localHops || opt.nWin != 1 || opt.coverPct)) { Abc_Print(-1, "-b needs one window per kick and no -L.\n"); return 1; }
    {
        Gia_Man_t* pGia = currentGia(pAbc, false);
        if (!pGia) { Abc_Print(-1, "There is no AIG (read a circuit, or use &r).\n"); return 1; }
        return opt.kick ? runKicks(pAbc, pGia, opt) : opt.ils ? runIls(pAbc, pGia, opt) : runAnalysis(pGia, opt);
    }
usage:
    Abc_Print(-2, "usage: lsv_reencode [-KntXCSNmPBTL num] [-o file] [-M sampler] [-R script] [-w file] [-W dir] [-akiebcsvh]\n");
    Abc_Print(-2, "\t         window-based re-encoding of output equivalence classes (see lsvReencode.cpp)\n");
    Abc_Print(-2, "\t-K num  : window leaf limit, <= 14 [default = %d]\n", opt.K);
    Abc_Print(-2, "\t-n num  : number of random window roots (-a), of kicks (-k) or of iterations (-i) [default = %d]\n", opt.nRoots);
    Abc_Print(-2, "\t-t num  : random codes per sampler (U, A, W1, W4) and window [default = %d]\n", opt.nTry);
    Abc_Print(-2, "\t-X num  : max windows (N <= 8) with all canonical codes enumerated [default = %d]\n", opt.xCap);
    Abc_Print(-2, "\t-C num  : max ANDs in a window [default = %d]\n", opt.maxCone);
    Abc_Print(-2, "\t-B num  : AND budget of each synthesized circuit = num * (|C| + 10), 0 = no limit [default = %d]\n", opt.budgetMul);
    Abc_Print(-2, "\t-S num  : random seed [default = 1]\n");
    Abc_Print(-2, "\t-o file : CSV output, one row per (window, sampler) [default = stdout]\n");
    Abc_Print(-2, "\t-M str  : sampler used by -k: U, A, W1, W4, base, P0 or N (null kick: C itself, only the node\n");
    Abc_Print(-2, "\t          order changes) [default = %s]\n", opt.method.c_str());
    Abc_Print(-2, "\t-R str  : re-optimization script run after each kick [default = \"%s\"]\n", opt.script.c_str());
    Abc_Print(-2, "\t-w file : write the best R(kick) AIG, if it beats R(x) and passes cec\n");
    Abc_Print(-2, "\t-W dir  : write R(x) as dir/Rx.aig and every R(kick) as dir/k<i>.aig; with one window per kick\n");
    Abc_Print(-2, "\t          also x (dir/x.aig), each kicked AIG (dir/y<i>.aig) and the origin of its nodes (dir/y<i>.lab)\n");
    Abc_Print(-2, "\t-a      : analysis only (phase 1): cost every sampled code, do not change the AIG\n");
    Abc_Print(-2, "\t-k      : kicks (phase 2): -n single kicks from the current AIG x, each followed by R, vs. R(x)\n");
    Abc_Print(-2, "\t-i      : iterated local search: kick the current AIG, run R, keep the result if it has fewer ANDs\n");
    Abc_Print(-2, "\t          (or as many and fewer levels); -W dir writes every accepted AIG as dir/a<iter>.aig\n");
    Abc_Print(-2, "\t-e      : -i also accepts results with as many ANDs and no more levels [default = %s]\n", opt.sideways ? "yes" : "no");
    Abc_Print(-2, "\t-T num  : -i time limit in seconds, 0 = none [default = %d]\n", opt.timeLimit);
    Abc_Print(-2, "\t-L num  : -i runs R only on the ANDs within num fanin/fanout steps of the replacement (minus the\n");
    Abc_Print(-2, "\t          nodes that depend on it through outside nodes),\n");
    Abc_Print(-2, "\t          0 = R on the whole AIG; accepted results are cec'd against the current AIG [default = %d]\n", opt.localHops);
    Abc_Print(-2, "\t-b      : keep the code bits: R runs with every code bit cut (exported as a CO, read back as a CI), so it\n");
    Abc_Print(-2, "\t          cannot optimize across them; reconnected afterwards [default = %s]\n", opt.keepCode ? "yes" : "no");
    Abc_Print(-2, "\t-N num  : kick only windows with at least num classes [default = %d]\n", opt.minN);
    Abc_Print(-2, "\t-m num  : windows re-encoded (one after another) per kick [default = %d]\n", opt.nWin);
    Abc_Print(-2, "\t-P num  : instead of -m, re-encode windows until their cones add up to num%% of AND(x) [default = off]\n");
    Abc_Print(-2, "\t-c      : toggle cec of every kicked AIG against x [default = yes]\n");
    Abc_Print(-2, "\t-F      : test hook: complement one replaced output, so cec must fail\n");
    Abc_Print(-2, "\t-s      : self-test of the code enumeration / canonical form\n");
    Abc_Print(-2, "\t-v      : toggle verbose output\n");
    Abc_Print(-2, "\t-h      : print the command usage\n");
    return 1;
}

void init(Abc_Frame_t* pAbc) { Cmd_CommandAdd(pAbc, "LSV", "lsv_reencode", Lsv_CommandReencode, 0); }
void destroy(Abc_Frame_t* pAbc) {}
Abc_FrameInitializer_t frame_initializer = { init, destroy };
struct PackageRegistrationManager {
    PackageRegistrationManager() { Abc_FrameAddInitializer(&frame_initializer); }
} lsvReencodeRegistrationManager;

}  // namespace
}  // namespace lsvpipe
