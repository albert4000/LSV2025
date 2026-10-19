// lsvPipeAcd.cpp -- ACD-style re-encoding of a 2-stage pipeline cut (src/develop/plan/acd_cut_reencoding.md).
//
// Phase 1 (analysis): group the cut signals, take a reconvergence-driven window (<= K leaves) per
//   group, enumerate the image Img(G) over the leaves, r = ceil(log2 |Img|), and a minimal
//   functionally-dependent subset S (greedy deletion).
// Phase 2 (synthesis): per group, (1) register only S and decode G \ S after the register, else
//   (2) re-encode G into r bits (existing G signals reused as code bits where possible, the rest
//   computed from the leaves). Encoders/decoders are synthesized inside ABC: truth table ->
//   (ISOP with don't-cares | DC=0 | DC=1) -> Kit_TruthToGia -> {none, dc2, syn2}, smallest kept.
//   A candidate is accepted iff AND + W*reg strictly decreases and both stage depths stay <= T.
// The result is checked with cec (registers as wires) against the original.
// Ported from src/develop/acd_reencode (phase1.cpp / phase2.py).
#include "lsvPipe.h"
#include "lsvAcdCore.h"

#include <random>
#include <set>
#include <unordered_map>

#include "aig/gia/giaAig.h"
#include "misc/util/utilTruth.h"
#include "proof/fra/fra.h"


namespace lsvpipe {

using std::string;
using std::vector;

// ---------------------------------------------------------------------------------------------
// Pipeline netlist (stage 0 / registers / stage 1)

struct Net {
    enum { DEAD = 0, PI = 1, AND = 2, REG = 3 };
    vector<char> kind, st;
    vector<int> a, b;
    vector<int> pos, regOf;   // regOf: original node -> its register id (-1 if none)
    int nPi = 0;

    int add(char k, int x, int y, char s) {
        kind.push_back(k); a.push_back(x); b.push_back(y); st.push_back(s);
        return (int)kind.size() - 1;
    }
    vector<int> topo() const {
        vector<int> order;
        vector<char> state(kind.size(), 0);
        vector<std::pair<int, int>> stack;
        for (int lit : pos) {
            int root = lit >> 1;
            if (!root || state[root]) continue;
            stack.push_back({ root, 0 });
            while (!stack.empty()) {
                int v = stack.back().first, i = stack.back().second;
                if (i == 0) { if (state[v]) { stack.pop_back(); continue; } state[v] = 1; }
                int nf = kind[v] == AND ? 2 : (kind[v] == REG ? 1 : 0);
                if (i < nf) {
                    stack.back().second++;
                    int u = (i == 0 ? a[v] : b[v]) >> 1;
                    if (u && !state[u]) stack.push_back({ u, 0 });
                } else { order.push_back(v); stack.pop_back(); }
            }
        }
        return order;
    }
    void sweep() {
        vector<char> live(kind.size(), 0);
        for (int v : topo()) live[v] = 1;
        for (size_t v = 1; v < kind.size(); v++) if (!live[v] && kind[v] != PI) kind[v] = DEAD;
        for (auto& r : regOf) if (r >= 0 && kind[r] == DEAD) r = -1;
    }
    void stats(int& ands, int& regs) const {
        ands = regs = 0;
        for (char k : kind) { ands += k == AND; regs += k == REG; }
    }
    // stage depths; returns false on an illegal structure
    bool depths(int& d0, int& d1) const {
        vector<int> lev(kind.size(), 0);
        d0 = d1 = 0;
        for (int v : topo()) {
            if (kind[v] == REG) {
                int x = a[v] >> 1;
                if (x && st[x] != 0) return false;
                continue;
            }
            if (kind[v] != AND) continue;
            for (int x : { a[v] >> 1, b[v] >> 1 }) {
                if (!x) continue;
                if (st[v] == 0 && (kind[x] == REG || st[x] != 0)) return false;
                if (st[v] == 1 && kind[x] != REG && st[x] != 1) return false;
            }
            lev[v] = 1 + std::max(lev[a[v] >> 1], lev[b[v] >> 1]);
            (st[v] ? d1 : d0) = std::max(st[v] ? d1 : d0, lev[v]);
        }
        for (int l : pos) {
            int x = l >> 1;
            if (x && kind[x] != REG && st[x] != 1) return false;
        }
        return true;
    }
    vector<int> instantiate(Gia_Man_t* p, const vector<int>& inputs, char stage) {
        vector<int> m(Gia_ManObjNum(p), 0);
        Gia_Obj_t* pObj;
        int i;
        Gia_ManForEachCi(p, pObj, i) m[Gia_ObjId(p, pObj)] = i < (int)inputs.size() ? inputs[i] : 0;
        auto lit = [&](int id, int c) { return m[id] ^ c; };
        Gia_ManForEachAnd(p, pObj, i)
            m[i] = 2 * add(AND, lit(Gia_ObjFaninId0(pObj, i), Gia_ObjFaninC0(pObj)),
                           lit(Gia_ObjFaninId1(pObj, i), Gia_ObjFaninC1(pObj)), stage);
        vector<int> outs;
        Gia_ManForEachCo(p, pObj, i) outs.push_back(lit(Gia_ObjFaninId0p(p, pObj), Gia_ObjFaninC0(pObj)));
        return outs;
    }
    void substituteRegs(const vector<std::pair<int, int>>& subst) {
        vector<int> to(kind.size(), -1);
        for (auto& s : subst) to[s.first] = s.second;
        auto f = [&](int l) { int x = l >> 1; return (x < (int)to.size() && to[x] >= 0) ? (to[x] ^ (l & 1)) : l; };
        for (size_t v = 0; v < kind.size(); v++)
            if (kind[v] == AND && st[v] == 1) { a[v] = f(a[v]); b[v] = f(b[v]); }
        for (auto& l : pos) l = f(l);
        for (auto& s : subst) kind[s.first] = DEAD;
        sweep();
    }
};

static Net buildNet(const Problem& P, const Solution& sol, const vector<int>& cut) {
    const PipeAig& g = P.g;
    int n = g.numNodes();
    Net N;
    N.nPi = g.nPi;
    N.kind.assign(n, Net::DEAD); N.st.assign(n, 0); N.a.assign(n, 0); N.b.assign(n, 0);
    N.regOf.assign(n, -1);
    for (int v = 1; v <= g.nPi; v++) N.kind[v] = Net::PI;
    for (int u : cut) N.regOf[u] = N.add(Net::REG, 2 * u, 0, 1);
    auto litIn = [&](int l, int s) {
        int u = l >> 1;
        return (u && s == 1 && sol.s[u] == 0) ? 2 * N.regOf[u] + (l & 1) : l;
    };
    for (int v = g.nPi + 1; v < n; v++) {
        if (!P.live[v]) continue;
        N.kind[v] = Net::AND; N.st[v] = (char)sol.s[v];
        N.a[v] = litIn(g.fanin0[v], sol.s[v]); N.b[v] = litIn(g.fanin1[v], sol.s[v]);
    }
    for (int l : g.poLit) N.pos.push_back(litIn(l, 1));
    N.sweep();
    return N;
}

// Sequential GIA (registers = flops, reset to their all-zero-input value) or combinational GIA
// (registers = wires).
static Gia_Man_t* netToGia(const Net& N, bool comb, Gia_Man_t* pOrig) {
    vector<int> order = N.topo(), regs;
    for (int v : order) if (N.kind[v] == Net::REG) regs.push_back(v);
    vector<char> z(N.kind.size(), 0);   // all-zero-input values
    for (int v : order) {
        if (N.kind[v] == Net::AND) z[v] = (z[N.a[v] >> 1] ^ (N.a[v] & 1)) & (z[N.b[v] >> 1] ^ (N.b[v] & 1));
        else if (N.kind[v] == Net::REG) z[v] = z[N.a[v] >> 1] ^ (N.a[v] & 1);
    }
    Gia_Man_t* p = Gia_ManStart((int)N.kind.size() + 1);
    Gia_ManHashAlloc(p);
    vector<int> lit(N.kind.size(), 0);
    for (int i = 1; i <= N.nPi; i++) lit[i] = Gia_ManAppendCi(p);
    vector<int> ro;
    if (!comb) for (size_t r = 0; r < regs.size(); r++) ro.push_back(Gia_ManAppendCi(p));
    std::unordered_map<int, int> roOf;
    for (size_t r = 0; r < regs.size(); r++) roOf[regs[r]] = (int)r;
    auto L = [&](int l) { return Abc_LitNotCond(lit[l >> 1], l & 1); };
    for (int v : order) {
        if (N.kind[v] == Net::AND) lit[v] = Gia_ManHashAnd(p, L(N.a[v]), L(N.b[v]));
        else if (N.kind[v] == Net::REG) lit[v] = comb ? L(N.a[v]) : Abc_LitNotCond(ro[roOf[v]], z[v]);
    }
    for (int l : N.pos) Gia_ManAppendCo(p, L(l));
    if (!comb) {
        for (int v : regs) Gia_ManAppendCo(p, Abc_LitNotCond(L(N.a[v]), z[v]));
        Gia_ManSetRegNum(p, (int)regs.size());
    }
    Gia_ManHashStop(p);
    copyNames(p, pOrig);
    return p;
}

// ---------------------------------------------------------------------------------------------
// Phase 2: one group

enum Status { ACCEPTED, REJ_AREA, REJ_DELAY, SKIPPED };

struct Ctx {
    const Problem& P;
    const AcdOptions& opt;
    int T;
};

static ll areaOf(const Net& N, ll W) { int a, r; N.stats(a, r); return a + W * r; }

// Evaluate a candidate built on a copy of the net. Returns 0 ok, 1 area, 2 delay, 3 illegal.
static int tryCandidate(const Ctx& C, const Net& N, const std::function<void(Net&)>& build, Net& out, string& note) {
    Net M = N;
    build(M);
    int d0, d1;
    ll a0 = areaOf(N, C.opt.regWeight), a1 = areaOf(M, C.opt.regWeight);
    if (!M.depths(d0, d1)) { note += "illegal "; return 3; }
    char buf[96];
    sprintf(buf, "%+lld:%d/%d ", a1 - a0, d0, d1);
    note += buf;
    if (a1 >= a0) return 1;
    if (std::max(d0, d1) > C.T && !C.opt.ignoreDelay) return 2;
    out = M;
    return 0;
}

static Status processGroup(const Ctx& C, Net& N, Group& G, string& note, string& method) {
    int m = (int)G.members.size();
    vector<int> regs(m);
    for (int j = 0; j < m; j++) regs[j] = N.regOf[G.members[j]];
    bool tried = false, sawDelay = false;
    // A candidate can only pay off if enc + dec < W*saving + removed ANDs; the window cone (+100)
    // is a loose bound on the removed logic. Raw ISOP/factored results beyond 8x that are dropped.
    auto budget = [&](int saving) { return 8 * (int)(C.opt.regWeight * saving + G.coneAnds + 100); };
    auto record = [&](int st) { tried = true; if (st == 2) sawDelay = true; };

    // 1. functional dependence
    int nS = (int)G.fdSet.size();
    if (nS < m && nS <= 16) {
        vector<int> code(G.nImg, 0);
        for (int q = 0; q < G.nImg; q++)
            for (int k = 0; k < nS; k++) code[q] |= G.U[G.fdSet[k]][q] << k;
        vector<int> rest;
        for (int j = 0; j < m; j++) if (std::find(G.fdSet.begin(), G.fdSet.end(), j) == G.fdSet.end()) rest.push_back(j);
        vector<vector<word>> on, dc;
        decoderTruths(G, code, nS, rest, on, dc);
        vector<Synth> dec = synthesize(on, dc, nS, budget(m - nS));
        if (dec.empty()) { note += "fd[dec over budget] "; record(1); }
        for (auto& d : dec) {
            Net out;
            note += "fd[dec=" + std::to_string(d.ands) + "/d" + std::to_string(d.depth) + "]";
            int st = tryCandidate(C, N, [&](Net& M) {
                vector<int> ins;
                for (int k : G.fdSet) ins.push_back(2 * regs[k]);
                vector<int> o = M.instantiate(d.p, ins, 1);
                vector<std::pair<int, int>> sub;
                for (size_t i = 0; i < rest.size(); i++) sub.push_back({ regs[rest[i]], o[i] });
                M.substituteRegs(sub);
            }, out, note);
            record(st);
            if (st == 0) { N = out; method = "fd"; freeSynth(dec); return ACCEPTED; }
        }
        freeSynth(dec);
    }
    // 2. general re-encoding
    if (G.r < m) {
        vector<int> reused, idx, code;
        int nnew;
        assignCode(G, reused, nnew, idx, code);
        vector<int> rest;
        for (int j = 0; j < m; j++) if (std::find(reused.begin(), reused.end(), j) == reused.end()) rest.push_back(j);
        vector<vector<word>> on, dc;
        decoderTruths(G, code, G.r, rest, on, dc);
        vector<Synth> dec = synthesize(on, dc, G.r, budget(m - G.r)), enc;
        int nL = (int)G.leaves.size();
        if (nnew) {
            int nw = nWordsOf(nL);
            vector<vector<word>> et(nnew, vector<word>(nw, 0)), none(nnew, vector<word>(nw, 0));
            for (size_t p = 0; p < G.inv.size(); p++)
                for (int t = 0; t < nnew; t++)
                    if ((idx[G.inv[p]] >> t) & 1) et[t][p >> 6] |= 1ULL << (p & 63);
            enc = synthesize(et, none, nL, budget(m - G.r));
        }
        if (dec.empty() || (nnew && enc.empty())) { note += "enc[over budget] "; record(1); }
        vector<std::pair<int, int>> combos;
        if (!dec.empty() && (!nnew || !enc.empty())) combos.push_back({ 0, 0 });
        if (!combos.empty() && (enc.size() > 1 || dec.size() > 1)) combos.push_back({ (int)enc.size() - 1, (int)dec.size() - 1 });
        for (auto& cb : combos) {
            Synth* e = nnew ? &enc[std::max(cb.first, 0)] : nullptr;
            Synth& d = dec[cb.second];
            note += "enc[r=" + std::to_string(G.r) + ",reuse=" + std::to_string(reused.size()) +
                    ",enc=" + std::to_string(e ? e->ands : 0) + ",dec=" + std::to_string(d.ands) + "/d" + std::to_string(d.depth) + "]";
            Net out;
            int st = tryCandidate(C, N, [&](Net& M) {
                vector<int> ins;
                for (int j : reused) ins.push_back(2 * regs[j]);
                if (e) {
                    vector<int> lv;
                    for (int l : G.leaves) lv.push_back(2 * l);
                    for (int o : M.instantiate(e->p, lv, 0)) ins.push_back(2 * M.add(Net::REG, o, 0, 1));
                }
                vector<int> o = M.instantiate(d.p, ins, 1);
                vector<std::pair<int, int>> sub;
                for (size_t i = 0; i < rest.size(); i++) sub.push_back({ regs[rest[i]], o[i] });
                M.substituteRegs(sub);
            }, out, note);
            record(st);
            if (st == 0) { N = out; method = "encode"; freeSynth(dec); freeSynth(enc); return ACCEPTED; }
        }
        freeSynth(dec); freeSynth(enc);
    }
    if (!tried) return SKIPPED;
    return sawDelay ? REJ_DELAY : REJ_AREA;
}

// ---------------------------------------------------------------------------------------------

static void globalSampling(const PipeAig& g, const vector<int>& cut, int logSamples) {
    int n = g.numNodes();
    vector<char> inCone(n, 0);
    for (int u : cut) inCone[u] = 1;
    for (int v = n - 1; v > g.nPi; v--)
        if (inCone[v]) inCone[g.fanin0[v] >> 1] = inCone[g.fanin1[v] >> 1] = 1;
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
        for (int v = g.nPi + 1; v < n; v++) {
            if (!inCone[v]) continue;
            u64 x = val[g.fanin0[v] >> 1], y = val[g.fanin1[v] >> 1];
            val[v] = ((g.fanin0[v] & 1) ? ~x : x) & ((g.fanin1[v] & 1) ? ~y : y);
        }
        std::pair<u64, u64> hw[64] = {};
        for (size_t j = 0; j < cut.size(); j++)
            for (u64 bits = val[cut[j]]; bits; bits &= bits - 1) {
                int p = __builtin_ctzll(bits);
                hw[p].first ^= keys.k1[j]; hw[p].second ^= keys.k2[j];
            }
        size_t valid = std::min<u64>(64, samples - 64 * w);
        for (size_t p = 0; p < valid; p++) h.push_back(hw[p]);
    }
    u64 N = countDistinct(h);
    printf("  global sampling : %llu distinct cut vectors over %llu %s patterns, ceil(log2 N) = %d%s\n",
           N, samples, exhaustive ? "exhaustive" : "random", ceilLog2(N), (!exhaustive && N == samples) ? " (saturated)" : "");
}


Gia_Man_t* runAcd(Problem& P, const Solution& sol, const AcdOptions& opt, Gia_Man_t* pOrig, string& msg) {
    const PipeAig& g = P.g;
    int n = g.numNodes(), T = P.T;
    if (P.S != 2) { msg = "ACD re-encoding supports 2 stages only"; return nullptr; }
    vector<int> cut = cutSets(P, sol)[1];
    // per-signal arrival (depth inside stage 0) and downstream depth (through stage 1)
    vector<int> arrival = stageDepth(P, sol.s), down = stageDown(P, sol.s), arr(n, 0), dn(n, 0);
    for (int u : cut) {
        arr[u] = g.isAnd(u) ? arrival[u] : 0;
        for (int k = g.foStart[u]; k < g.foStart[u + 1]; k++) {
            int w = g.foList[k];
            if (P.live[w] && sol.s[w] == 1) dn[u] = std::max(dn[u], down[w]);
        }
    }
    int Tref = T > 0 ? T : std::max(*std::max_element(arrival.begin(), arrival.end()), 1);
    if (opt.global) globalSampling(g, cut, opt.logSamples);

    // Phase 1: grouping (plan: leaf overlap; slack: also within slack classes)
    int nc = (int)cut.size();
    vector<int> cls(nc, 0);
    if (opt.slackGrouping)
        for (int j = 0; j < nc; j++) cls[j] = dn[cut[j]] >= Tref ? 2 : (arr[cut[j]] >= Tref ? 1 : 0);
    vector<vector<int>> single(nc);
    for (int j = 0; j < nc; j++) single[j] = findWindow(g, { cut[j] }, opt.K).leaves;
    vector<char> used(nc, 0);
    vector<Group> groups;
    for (int seed = 0; seed < nc; seed++) {
        if (used[seed]) continue;
        used[seed] = 1;
        Group G;
        G.members = { cut[seed] };
        Window w = findWindow(g, G.members, opt.K);
        while (cls[seed] != 2) {
            std::set<int> Ls(w.leaves.begin(), w.leaves.end());
            int best = -1, bestOv = 0;
            for (int j = 0; j < nc; j++) {
                if (used[j] || cls[j] != cls[seed]) continue;
                int ov = 0;
                for (int l : single[j]) ov += (int)Ls.count(l);
                if (ov > bestOv) { bestOv = ov; best = j; }
            }
            if (best < 0) break;
            vector<int> M2 = G.members; M2.push_back(cut[best]);
            Window w2 = findWindow(g, M2, opt.K);
            if (!w2.ok) break;
            used[best] = 1; G.members = M2; w = w2;
        }
        if (!w.ok) { msg = "single-signal window exceeds K"; return nullptr; }
        analyzeGroup(g, G, w);
        G.coneAnds = (int)w.cone.size();
        for (int u : G.members) { G.maxArr = std::max(G.maxArr, arr[u]); G.maxDown = std::max(G.maxDown, dn[u]); }
        groups.push_back(std::move(G));
    }
    int sumSav = 0, sumFd = 0, sumFeas = 0, ge2 = 0, reqL = 0;
    for (auto& G : groups) {
        int m = (int)G.members.size(), up = Tref - G.maxArr, dw = Tref - G.maxDown;
        int sav = m - G.r, fd = m - (int)G.fdSet.size();
        int feas = (up >= 1 && dw >= 1) ? sav : (dw >= 1 ? fd : 0);
        sumSav += sav; sumFd += fd; sumFeas += feas; ge2 += (up >= 2 || dw >= 2); reqL += G.r == (int)G.leaves.size();
        if (opt.verbose && sav > 0)
            printf("    group |G|=%-3d |L|=%-2d |Img|=%-6d r=%-2d saving=%-3d fd=%-3d slack up/down=%d/%d%s\n", m,
                   (int)G.leaves.size(), G.nImg, G.r, sav, (int)G.fdSet.size(), up, dw, G.r == (int)G.leaves.size() ? " r==|L|" : "");
    }
    bool go = sumFeas >= std::max(0.05 * nc, 4.0) && ge2 > 0;
    printf("  phase 1 (K=%d, %s grouping): cut %d, %d groups, saving %d (fd-only %d), slack-feasible %d, r==|L| groups %d, go=%s\n",
           opt.K, opt.slackGrouping ? "slack" : "plan", nc, (int)groups.size(), sumSav, sumFd, sumFeas, reqL, go ? "yes" : "no");
    if (opt.analysisOnly) return nullptr;

    // Phase 2
    Net N = buildNet(P, sol, cut);
    int and0, reg0, d0, d1;
    N.stats(and0, reg0);
    if (!N.depths(d0, d1)) { msg = "internal error: initial pipeline is not legal"; return nullptr; }
    Ctx C{ P, opt, T > 0 ? T : std::max(d0, d1) };
    vector<int> order;
    for (size_t i = 0; i < groups.size(); i++) if ((int)groups[i].members.size() > groups[i].r) order.push_back((int)i);
    std::stable_sort(order.begin(), order.end(), [&](int x, int y) {
        return (int)groups[x].members.size() - groups[x].r > (int)groups[y].members.size() - groups[y].r; });
    int cnt[4] = { 0, 0, 0, 0 };
    for (int i : order) {
        string note, method;
        ll before = areaOf(N, opt.regWeight);
        Status st = processGroup(C, N, groups[i], note, method);
        cnt[st]++;
        if (opt.verbose) {
            static const char* names[4] = { "accepted", "rejected (area)", "rejected (delay)", "skipped" };
            printf("    group |G|=%-3d r=%-2d fd=%-3d -> %-16s %-6s area %+lld  [%s]\n", (int)groups[i].members.size(), groups[i].r,
                   (int)groups[i].fdSet.size(), names[st], method.c_str(), areaOf(N, opt.regWeight) - before, note.c_str());
        }
    }
    int and1, reg1;
    N.stats(and1, reg1);
    N.depths(d0, d1);
    Gia_Man_t* pComb = netToGia(N, true, pOrig);
    int eq = cecGiaSafe(pOrig, pComb);
    Gia_ManStop(pComb);
    if (eq != 1) { msg = eq == 0 ? "CEC FAILED: re-encoded circuit is not equivalent; result discarded"
                                  : "CEC undecided; result discarded"; return nullptr; }
    Gia_Man_t* pSeq = netToGia(N, false, pOrig);
    and1 = Gia_ManAndNum(pSeq);      // after structural hashing of the final netlist
    reg1 = Gia_ManRegNum(pSeq);
    ll a0 = and0 + opt.regWeight * reg0, a1 = and1 + opt.regWeight * reg1;
    printf("  phase 2: groups tried %d, accepted %d, rejected area %d, rejected delay %d\n",
           (int)order.size(), cnt[ACCEPTED], cnt[REJ_AREA], cnt[REJ_DELAY]);
    printf("  registers %d -> %d, ANDs %d -> %d, area (AND + %lld*reg) %lld -> %lld (%+.3f%%), stage depths %d/%d%s\n",
           reg0, reg1, and0, and1, opt.regWeight, a0, a1, a0 ? 100.0 * (a1 - a0) / a0 : 0.0, d0, d1,
           (T > 0 && std::max(d0, d1) > T) ? " (exceeds T)" : "");
    printf("  cec (registers as wires) vs. original: equivalent\n");
    return pSeq;
}

}  // namespace lsvpipe
