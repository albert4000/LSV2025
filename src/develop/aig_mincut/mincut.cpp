// mincut.cpp -- fixed-stage pipeline register minimization on a combinational AIG.
//
//   Stage A: S = 2, node-split max-flow / min-cut (exact).
//   Stage B: general S, difference-constraint LP solved through its min-cost-flow dual (exact).
//   Stage C: per-stage depth limit T, level windows + lazy path constraints.
//   Stage D: stage_assign.json, cuts.json, report, verification, BLIF dump for ABC dsec.
//
// See ../plan/aig_mincut_spec.md.
#include "aig.h"
#include "flow.h"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

using std::string;
using std::vector;

struct Options {
    string input;
    int S = 2;
    ll W = 1;
    string side = "source";
    int T = 0;                    // 0 = no depth limit
    string method = "auto";       // auto | flow | lp
    string outDir = ".";
    bool writeJson = true;
    string blif, refBlif;
    int brute = 0;                // max #ANDs for brute force, 0 = off
    bool specGraph = false;       // stage A with the original spec 4.1 v_in->u_out arc (over-constrained)
    int simCycles = 32;
    bool quiet = false;
};

struct Problem {
    const Aig& g;
    int S;
    ll W;
    int T;
    vector<char> live, poDrv;
    vector<int> L, R;                     // AND level from PI side / to PO side (excluding the node)
    vector<int> lo, hi;                   // stage window
    vector<std::pair<int, int>> lazy;     // (u, v): s(v) >= s(u) + 1
    Problem(const Aig& g_, int S_, ll W_, int T_) : g(g_), S(S_), W(W_), T(T_) {}
};

struct Solution {
    vector<int> s, m;   // -1 for nodes not in the problem
    ll regs = 0;        // number of registers, sum (m - s)
};

static double nowSec() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

// ---------------------------------------------------------------------------------------------
// Problem setup

static void buildProblem(Problem& P) {
    const Aig& g = P.g;
    int n = g.numNodes();
    P.live.assign(n, 0); P.poDrv.assign(n, 0);
    for (int lit : g.poLit) if (lit >> 1) P.poDrv[lit >> 1] = P.live[lit >> 1] = 1;
    for (int v = n - 1; v > g.nPi; v--)        // dangling removal: transitive fanin of POs
        if (P.live[v]) { P.live[g.fanin0[v] >> 1] = P.live[g.fanin1[v] >> 1] = 1; }
    P.live[0] = 0;

    P.L.assign(n, 0); P.R.assign(n, 0);
    for (int v = g.nPi + 1; v < n; v++) P.L[v] = 1 + std::max(P.L[g.fanin0[v] >> 1], P.L[g.fanin1[v] >> 1]);
    for (int v = n - 1; v > 0; v--) {
        if (!P.live[v]) continue;
        for (int k = g.foStart[v]; k < g.foStart[v + 1]; k++) {
            int w = g.foList[k];
            if (P.live[w]) P.R[v] = std::max(P.R[v], P.R[w] + 1);
        }
    }
    P.lo.assign(n, 0); P.hi.assign(n, P.S - 1);
    for (int v = 1; v <= g.nPi; v++) P.hi[v] = 0;
    if (P.T > 0) {
        for (int v = g.nPi + 1; v < n; v++) {
            P.lo[v] = (P.L[v] + P.T - 1) / P.T - 1;                       // ASAP
            P.hi[v] = P.S - (P.R[v] + 1 + P.T - 1) / P.T;                  // ALAP
        }
    }
}

static void computeM(const Problem& P, Solution& sol) {
    const Aig& g = P.g;
    int n = g.numNodes();
    sol.m.assign(n, -1);
    sol.regs = 0;
    for (int u = 1; u < n; u++) {
        if (!P.live[u]) continue;
        int m = sol.s[u];
        if (P.poDrv[u]) m = P.S - 1;
        for (int k = g.foStart[u]; k < g.foStart[u + 1]; k++)
            if (P.live[g.foList[k]]) m = std::max(m, sol.s[g.foList[k]]);
        sol.m[u] = m;
        sol.regs += m - sol.s[u];
    }
}

// Combinational depth inside each stage (PIs and registered fanins count 0).
static vector<int> stageDepth(const Problem& P, const vector<int>& s, vector<int>* pred = nullptr) {
    const Aig& g = P.g;
    int n = g.numNodes();
    vector<int> d(n, 0);
    if (pred) pred->assign(n, -1);
    for (int v = g.nPi + 1; v < n; v++) {
        if (!P.live[v]) continue;
        d[v] = 1;
        for (int lit : { g.fanin0[v], g.fanin1[v] }) {
            int u = lit >> 1;
            if (g.isAnd(u) && s[u] == s[v] && d[u] + 1 > d[v]) { d[v] = d[u] + 1; if (pred) (*pred)[v] = u; }
        }
    }
    return d;
}

// ---------------------------------------------------------------------------------------------
// Stage A: S = 2 via node-split min-cut
//
// Source side X = stage 0. v_in in X <=> s(v) = 0; v is registered iff v_in in X, v_out not in X.
// The closure arc goes v_in -> u_in (fanin u must be in stage 0 if v is). The original spec 4.1 wrote
// v_in -> u_out, which additionally forbids a registered u from also feeding a stage-0 fanout
// directly, so it over-counts whenever such fanout sharing is optimal (specArc = true keeps it).

static bool solveFlow(const Problem& P, const string& side, bool specArc, Solution& sol, ll& cutValue) {
    const Aig& g = P.g;
    int n = g.numNodes();
    auto vin = [](int v) { return 2 * v; };
    auto vout = [](int v) { return 2 * v + 1; };
    int src = 2 * n, snk = 2 * n + 1;
    FlowGraph G(2 * n + 2);
    for (int v = 1; v < n; v++) {
        if (!P.live[v]) continue;
        G.addArc(vin(v), vout(v), P.W);
        if (g.isPi(v)) G.addArc(src, vin(v), FLOW_INF);
        if (P.poDrv[v]) G.addArc(vout(v), snk, FLOW_INF);
        if (g.isAnd(v)) {
            for (int lit : { g.fanin0[v], g.fanin1[v] }) {
                int u = lit >> 1;
                if (u == 0) continue;
                G.addArc(vout(u), vin(v), FLOW_INF);
                G.addArc(vin(v), specArc ? vout(u) : vin(u), FLOW_INF);   // closed cut
            }
        }
    }
    cutValue = G.maxFlow(src, snk);
    if (cutValue < 0) return false;
    vector<char> X;
    if (side == "sink") { X = G.reachTo(snk); for (auto& x : X) x = !x; }
    else X = G.reachFrom(src);
    sol.s.assign(n, -1);
    for (int v = 1; v < n; v++) if (P.live[v]) sol.s[v] = X[vin(v)] ? 0 : 1;
    computeM(P, sol);
    return true;
}

// ---------------------------------------------------------------------------------------------
// Stage B: LP over s(v), m(u) solved as the dual of a min-cost flow.
//
// Every LP constraint x_b - x_a <= c becomes an uncapacitated arc a->b of cost c. The objective
// sum W*(m(u) - s(u)) puts supply W on m(u) and demand W on s(u). Optimal potentials of the
// final residual graph are optimal LP solutions; shortest distances from the reference node z
// give the pointwise-maximal one (cut nearest the PIs), distances to z the minimal one.

struct LpInfo { int phases = 0; ll mcfCost = 0; };

static int solveLp(const Problem& P, const string& side, Solution& sol, LpInfo& info) {
    const Aig& g = P.g;
    int n = g.numNodes();
    vector<int> idx(n, -1), nodes;
    for (int v = 1; v < n; v++) if (P.live[v]) { idx[v] = (int)nodes.size(); nodes.push_back(v); }
    int K = (int)nodes.size();
    auto sI = [&](int v) { return 1 + 2 * idx[v]; };
    auto mI = [&](int v) { return 2 + 2 * idx[v]; };
    const int z = 0, SS = 1 + 2 * K, ST = 2 + 2 * K;
    const int S1 = P.S - 1;

    // Initial feasible solution = ASAP schedule under all forward constraints; it also
    // decides feasibility since every lower-bound constraint points forward in topo order.
    vector<vector<int>> lazyPred(n);
    for (auto& c : P.lazy) lazyPred[c.second].push_back(c.first);
    vector<ll> pot(2 * K + 3, 0);
    vector<int> s0(n, 0);
    for (int v : nodes) {
        int s = P.lo[v];
        if (g.isAnd(v)) {
            for (int lit : { g.fanin0[v], g.fanin1[v] }) if (lit >> 1) s = std::max(s, s0[lit >> 1]);
            for (int u : lazyPred[v]) s = std::max(s, s0[u] + 1);
        }
        if (s > P.hi[v]) return 1;
        s0[v] = s;
        pot[sI(v)] = s; pot[mI(v)] = S1;
    }
    pot[SS] = S1; pot[ST] = 0;

    FlowGraph G(2 * K + 3);
    auto le = [&](int a, int b, ll c) { G.addArc(a, b, FLOW_INF, c); };   // x_b - x_a <= c
    for (int v : nodes) {
        le(z, sI(v), P.hi[v]);
        le(sI(v), z, -P.lo[v]);
        le(mI(v), sI(v), 0);                 // s(v) <= m(v)
        le(z, mI(v), S1);                    // m(v) <= S-1 (never binding at the optimum)
        if (P.poDrv[v]) le(mI(v), z, -S1);   // m(v) >= S-1
        if (g.isAnd(v)) {
            for (int lit : { g.fanin0[v], g.fanin1[v] }) {
                int u = lit >> 1;
                if (u == 0) continue;
                le(sI(v), sI(u), 0);         // s(u) <= s(v)
                le(mI(u), sI(v), 0);         // s(v) <= m(u)
            }
        }
        G.addArc(SS, mI(v), P.W);
        G.addArc(sI(v), ST, P.W);
    }
    for (auto& c : P.lazy) le(sI(c.second), sI(c.first), -1);   // s(u) - s(v) <= -1
    for (int e = 0; e < G.numArcs(); e++)
        if (G.cap[e] > 0 && G.cost[e] + pot[G.tail(e)] - pot[G.to[e]] < 0)
            throw std::logic_error("initial potential is not feasible");

    FlowGraph::McfResult r = G.minCostFlow(SS, ST, pot);
    if (r.unbounded || r.flow != (ll)K * P.W) throw std::logic_error("min-cost flow did not route all supply");
    info.phases = r.phases;
    info.mcfCost = r.cost;

    vector<char> skip(G.n, 0);
    skip[SS] = skip[ST] = 1;
    bool rev = (side == "sink");
    vector<ll> d = G.dijkstra(z, pot, rev, &skip);
    auto x = [&](int i) -> ll {
        if (d[i] >= FLOW_INF) throw std::logic_error("LP node unreachable from z");
        return rev ? -(d[i] - pot[i] + pot[z]) : d[i] + pot[i] - pot[z];
    };
    sol.s.assign(n, -1);
    vector<int> mLp(n, -1);
    for (int v : nodes) { sol.s[v] = (int)x(sI(v)); mLp[v] = (int)x(mI(v)); }
    computeM(P, sol);
    for (int v : nodes)
        if (mLp[v] != sol.m[v]) throw std::logic_error("LP m(u) is not tight at node " + std::to_string(v));
    if (sol.regs * P.W != -r.cost) throw std::logic_error("primal/dual objective mismatch");
    return 0;
}

// Stage C wrapper: re-solve with lazily added path constraints until every stage has depth <= T.
// Each added constraint is necessary (a same-stage path of T+1 ANDs is illegal), so the final
// solution is optimal for the depth-constrained problem and infeasibility is exact.
static int solveLpTimed(Problem& P, const string& side, Solution& sol, LpInfo& info, int& rounds) {
    rounds = 0;
    while (true) {
        rounds++;
        if (solveLp(P, side, sol, info)) return 1;
        if (P.T <= 0) return 0;
        vector<int> pred;
        vector<int> d = stageDepth(P, sol.s, &pred);
        size_t before = P.lazy.size();
        for (int v = P.g.nPi + 1; v < P.g.numNodes(); v++) {
            if (!P.live[v] || d[v] <= P.T) continue;
            int u = v;
            for (int k = 0; k < P.T; k++) u = pred[u];
            P.lazy.push_back({ u, v });
        }
        if (P.lazy.size() == before) return 0;
    }
}

// ---------------------------------------------------------------------------------------------
// Brute force (small circuits): enumerate all monotone stage assignments.

static ll bruteForce(const Problem& P, Solution& best) {
    const Aig& g = P.g;
    int n = g.numNodes();
    vector<int> ands;
    for (int v = g.nPi + 1; v < n; v++) if (P.live[v]) ands.push_back(v);
    Solution cur;
    cur.s.assign(n, -1);
    for (int v = 1; v <= g.nPi; v++) if (P.live[v]) cur.s[v] = 0;
    vector<int> d(n, 0);
    ll bestRegs = -1;
    std::function<void(size_t)> rec = [&](size_t i) {
        if (i == ands.size()) {
            computeM(P, cur);
            if (bestRegs < 0 || cur.regs < bestRegs) { bestRegs = cur.regs; best = cur; }
            return;
        }
        int v = ands[i], lo = P.lo[v];
        for (int lit : { g.fanin0[v], g.fanin1[v] }) if (lit >> 1) lo = std::max(lo, cur.s[lit >> 1]);
        for (int s = lo; s <= P.hi[v]; s++) {
            cur.s[v] = s;
            d[v] = 1;
            for (int lit : { g.fanin0[v], g.fanin1[v] }) {
                int u = lit >> 1;
                if (g.isAnd(u) && cur.s[u] == s) d[v] = std::max(d[v], d[u] + 1);
            }
            if (P.T > 0 && d[v] > P.T) continue;
            rec(i + 1);
        }
        cur.s[v] = -1;
    };
    rec(0);
    return bestRegs;
}

// ---------------------------------------------------------------------------------------------
// Verification

// Structural checks: range, monotonicity, chain lengths, balanced register count on every path.
static string verifyStructure(const Problem& P, const Solution& sol) {
    const Aig& g = P.g;
    int n = g.numNodes();
    vector<int> regsTo(n, 0);   // registers between the PIs and v's combinational output
    for (int v = 1; v < n; v++) {
        if (!P.live[v]) continue;
        int s = sol.s[v];
        if (s < 0 || s > P.S - 1) return "stage out of range at node " + std::to_string(v);
        if (g.isPi(v) && s != 0) return "PI not in stage 0: node " + std::to_string(v);
        if (sol.m[v] < s) return "negative chain length at node " + std::to_string(v);
        if (!g.isAnd(v)) continue;
        int r = -1;
        for (int lit : { g.fanin0[v], g.fanin1[v] }) {
            int u = lit >> 1;
            if (u == 0) continue;
            if (s < sol.s[u]) return "edge goes backwards in stage: " + std::to_string(u) + "->" + std::to_string(v);
            if (s > sol.m[u]) return "chain too short for edge " + std::to_string(u) + "->" + std::to_string(v);
            int ru = regsTo[u] + (s - sol.s[u]);
            if (r >= 0 && ru != r) return "unbalanced paths into node " + std::to_string(v);
            r = ru;
        }
        regsTo[v] = r < 0 ? s : r;
    }
    for (int lit : g.poLit) {
        int u = lit >> 1;
        if (!u) continue;
        if (sol.m[u] < P.S - 1) return "PO driver chain too short at node " + std::to_string(u);
        if (regsTo[u] + (P.S - 1 - sol.s[u]) != P.S - 1) return "PO path does not cross S-1 registers";
    }
    if (P.T > 0) {
        vector<int> d = stageDepth(P, sol.s);
        for (int v = 1; v < n; v++) if (P.live[v] && d[v] > P.T) return "stage depth exceeds T at node " + std::to_string(v);
    }
    return "";
}

// Values of all nodes with all PIs at 0 (the register reset values that make the
// pipeline equivalent to the delayed original from cycle 0).
static vector<char> zeroValues(const Aig& g) {
    vector<char> z(g.numNodes(), 0);
    for (int v = g.nPi + 1; v < g.numNodes(); v++)
        z[v] = (z[g.fanin0[v] >> 1] ^ (g.fanin0[v] & 1)) & (z[g.fanin1[v] >> 1] ^ (g.fanin1[v] & 1));
    return z;
}

// Cycle-accurate 64-way simulation of the pipeline against the original AIG fed with inputs
// delayed by S-1 cycles. Returns "" on success.
static string simulate(const Problem& P, const Solution& sol, int cycles, unsigned seed) {
    const Aig& g = P.g;
    int n = g.numNodes(), S = P.S;
    typedef unsigned long long u64;
    std::mt19937_64 rng(seed);
    vector<char> z0 = zeroValues(g);
    vector<vector<u64>> hist(S, vector<u64>(n));   // hist[t % S][v] = pipeline value of v at cycle t
    vector<vector<u64>> inHist(S, vector<u64>(g.nPi + 1, 0));
    vector<u64> ref(n);
    for (int v = 0; v < n; v++) for (int t = 0; t < S; t++) hist[t][v] = z0[v] ? ~0ULL : 0;
    auto lit = [](u64 x, int l) { return (l & 1) ? ~x : x; };
    for (int t = 0; t < cycles; t++) {
        vector<u64>& cur = hist[t % S];
        cur[0] = 0;
        for (int v = 1; v <= g.nPi; v++) cur[v] = inHist[t % S][v] = rng();
        for (int v = g.nPi + 1; v < n; v++) {
            if (!P.live[v]) continue;
            u64 a[2];
            int ls[2] = { g.fanin0[v], g.fanin1[v] };
            for (int k = 0; k < 2; k++) {
                int u = ls[k] >> 1, back = u ? sol.s[v] - sol.s[u] : 0;
                u64 x = (t - back >= 0) ? hist[(t - back) % S][u] : (z0[u] ? ~0ULL : 0);
                a[k] = lit(x, ls[k]);
            }
            cur[v] = a[0] & a[1];
        }
        // reference: original logic on the inputs of cycle t-(S-1)
        int tr = t - (S - 1);
        ref[0] = 0;
        for (int v = 1; v <= g.nPi; v++) ref[v] = tr >= 0 ? inHist[tr % S][v] : 0;
        for (int v = g.nPi + 1; v < n; v++) ref[v] = lit(ref[g.fanin0[v] >> 1], g.fanin0[v]) & lit(ref[g.fanin1[v] >> 1], g.fanin1[v]);
        for (size_t o = 0; o < g.poLit.size(); o++) {
            int l = g.poLit[o], u = l >> 1;
            int back = u ? (S - 1 - sol.s[u]) : 0;
            u64 x = (t - back >= 0) ? hist[(t - back) % S][u] : (z0[u] ? ~0ULL : 0);
            if (lit(x, l) != lit(ref[u], l))
                return "simulation mismatch at cycle " + std::to_string(t) + " on PO " + g.poName[o];
        }
    }
    return "";
}

// ---------------------------------------------------------------------------------------------
// Output

static string sanitize(const string& s) {
    string r;
    for (char c : s) r += (c == ' ' || c == '\t' || c == '#' || c == '=') ? '_' : c;
    return r;
}

static string jsonStr(const string& s) {
    string r = "\"";
    for (char c : s) {
        if (c == '"' || c == '\\') r += '\\';
        if ((unsigned char)c < 0x20) { r += ' '; continue; }
        r += c;
    }
    return r + "\"";
}

// Pipeline as BLIF: registers reset to the node's all-zero-input value.
static void writePipelineBlif(const Problem& P, const Solution& sol, const string& path) {
    const Aig& g = P.g;
    std::ofstream f(path);
    if (!f) throw std::runtime_error("cannot write " + path);
    vector<char> z0 = zeroValues(g);
    auto base = [&](int v) { return g.isPi(v) ? sanitize(g.piName[v - 1]) : (v ? "n" + std::to_string(v) : string("const0")); };
    auto sig = [&](int u, int k) { return k ? base(u) + "_r" + std::to_string(k) : base(u); };
    f << ".model pipeline\n.inputs";
    for (int i = 0; i < g.nPi; i++) f << " " << sanitize(g.piName[i]);
    f << "\n.outputs";
    for (auto& o : g.poName) f << " " << sanitize(o);
    f << "\n.names const0\n";
    for (int u = 1; u < g.numNodes(); u++) {
        if (!P.live[u]) continue;
        for (int k = 1; k <= sol.m[u] - sol.s[u]; k++)
            f << ".latch " << sig(u, k - 1) << " " << sig(u, k) << " " << (int)z0[u] << "\n";
        if (!g.isAnd(u)) continue;
        int a = g.fanin0[u], b = g.fanin1[u];
        f << ".names " << sig(a >> 1, (a >> 1) ? sol.s[u] - sol.s[a >> 1] : 0) << " "
          << sig(b >> 1, (b >> 1) ? sol.s[u] - sol.s[b >> 1] : 0) << " " << base(u) << "\n"
          << ((a & 1) ? '0' : '1') << ((b & 1) ? '0' : '1') << " 1\n";
    }
    for (size_t o = 0; o < g.poLit.size(); o++) {
        int l = g.poLit[o], u = l >> 1;
        f << ".names " << sig(u, u ? P.S - 1 - sol.s[u] : 0) << " " << sanitize(g.poName[o]) << "\n"
          << ((l & 1) ? "0 1\n" : "1 1\n");
    }
    f << ".end\n";
}

// Reference: original logic with every PI delayed by S-1 registers (reset 0).
static void writeReferenceBlif(const Problem& P, const string& path) {
    const Aig& g = P.g;
    std::ofstream f(path);
    if (!f) throw std::runtime_error("cannot write " + path);
    int D = P.S - 1;
    auto sig = [&](int v) {
        if (g.isPi(v)) return sanitize(g.piName[v - 1]) + (D ? "_d" + std::to_string(D) : string());
        return v ? "n" + std::to_string(v) : string("const0");
    };
    f << ".model reference\n.inputs";
    for (int i = 0; i < g.nPi; i++) f << " " << sanitize(g.piName[i]);
    f << "\n.outputs";
    for (auto& o : g.poName) f << " " << sanitize(o);
    f << "\n.names const0\n";
    for (int i = 0; i < g.nPi; i++) {
        string p = sanitize(g.piName[i]);
        for (int k = 1; k <= D; k++)
            f << ".latch " << (k == 1 ? p : p + "_d" + std::to_string(k - 1)) << " " << p << "_d" << k << " 0\n";
    }
    for (int v = g.nPi + 1; v < g.numNodes(); v++) {
        if (!P.live[v]) continue;
        int a = g.fanin0[v], b = g.fanin1[v];
        f << ".names " << sig(a >> 1) << " " << sig(b >> 1) << " " << sig(v) << "\n"
          << ((a & 1) ? '0' : '1') << ((b & 1) ? '0' : '1') << " 1\n";
    }
    for (size_t o = 0; o < g.poLit.size(); o++) {
        int l = g.poLit[o];
        f << ".names " << sig(l >> 1) << " " << sanitize(g.poName[o]) << "\n" << ((l & 1) ? "0 1\n" : "1 1\n");
    }
    f << ".end\n";
}

static void writeStageAssign(const Problem& P, const Solution& sol, const string& path) {
    const Aig& g = P.g;
    std::ofstream f(path);
    if (!f) throw std::runtime_error("cannot write " + path);
    f << "{\n  \"stages\": " << P.S << ",\n  \"nodes\": [";
    bool first = true;
    for (int v = 1; v < g.numNodes(); v++) {
        if (!P.live[v]) continue;
        f << (first ? "\n" : ",\n") << "    {\"id\": " << v << ", \"aiger_var\": " << g.origVar[v]
          << ", \"type\": \"" << (g.isPi(v) ? "pi" : "and") << "\""
          << (g.isPi(v) ? ", \"name\": " + jsonStr(g.piName[v - 1]) : string())
          << ", \"s\": " << sol.s[v] << ", \"m\": " << sol.m[v] << ", \"po_driver\": " << (P.poDrv[v] ? "true" : "false") << "}";
        first = false;
    }
    f << "\n  ]\n}\n";
}

// Signals on boundary i (between stage i-1 and stage i): every u with s(u) < i <= m(u).
static vector<vector<int>> cutSets(const Problem& P, const Solution& sol) {
    vector<vector<int>> cuts(P.S);
    for (int u = 1; u < P.g.numNodes(); u++)
        if (P.live[u])
            for (int i = sol.s[u] + 1; i <= sol.m[u]; i++) cuts[i].push_back(u);
    return cuts;
}

static void writeCuts(const Problem& P, const Solution& sol, const string& path) {
    const Aig& g = P.g;
    std::ofstream f(path);
    if (!f) throw std::runtime_error("cannot write " + path);
    vector<vector<int>> cuts = cutSets(P, sol);
    int n = g.numNodes();
    // arrival = depth of u inside its own stage; down[w] = ANDs on the longest path from w
    // (inclusive) to the end of w's stage (a PO or a register).
    vector<int> arrival = stageDepth(P, sol.s), down(n, 0);
    for (int w = n - 1; w > g.nPi; w--) {
        if (!P.live[w]) continue;
        down[w] = 1;
        for (int k = g.foStart[w]; k < g.foStart[w + 1]; k++) {
            int x = g.foList[k];
            if (P.live[x] && sol.s[x] == sol.s[w]) down[w] = std::max(down[w], down[x] + 1);
        }
    }
    f << "{\n  \"stages\": " << P.S << ",\n  \"reg_weight\": " << P.W << ",\n  \"cuts\": [";
    for (int i = 1; i < P.S; i++) {
        f << (i > 1 ? ",\n" : "\n") << "    {\"boundary\": " << i << ", \"width\": " << cuts[i].size() << ", \"signals\": [";
        bool first = true;
        for (int u : cuts[i]) {
            // downstream: longest AND path after the register at boundary i (through stage i) to its end;
            // sameStageUses: distinct fanouts that read u combinationally in u's own stage
            int downstream = 0, sameStageUses = 0, prev = -1;
            for (int k = g.foStart[u]; k < g.foStart[u + 1]; k++) {
                int w = g.foList[k];
                if (!P.live[w] || w == prev) continue;
                prev = w;
                if (sol.s[w] == i) downstream = std::max(downstream, down[w]);
                if (sol.s[w] == sol.s[u]) sameStageUses++;
            }
            f << (first ? "\n" : ",\n") << "      {\"node\": " << u << ", \"aiger_var\": " << g.origVar[u]
              << ", \"lit\": " << 2 * g.origVar[u] << ", \"src_stage\": " << sol.s[u]
              << ", \"chain_pos\": " << i - sol.s[u] << ", \"chain_len\": " << sol.m[u] - sol.s[u]
              << ", \"balancing\": " << (g.isPi(u) ? "true" : "false")
              << ", \"po_driver\": " << (P.poDrv[u] ? "true" : "false")
              << ", \"level\": " << P.L[u] << ", \"arrival\": " << (g.isAnd(u) ? arrival[u] : 0)
              << ", \"downstream\": " << downstream << ", \"level_to_po\": " << P.R[u]
              << ", \"same_stage_uses\": " << sameStageUses << "}";
            first = false;
        }
        f << (cuts[i].empty() ? "]}" : "\n    ]}");
    }
    f << "\n  ]\n}\n";
}

static void printReport(const Problem& P, const Solution& sol) {
    const Aig& g = P.g;
    vector<vector<int>> cuts = cutSets(P, sol);
    ll balancing = 0;
    for (int v = 1; v <= g.nPi; v++) if (P.live[v]) balancing += sol.m[v] - sol.s[v];
    vector<int> d = stageDepth(P, sol.s);
    vector<int> depth(P.S, 0), size(P.S, 0);
    for (int v = g.nPi + 1; v < g.numNodes(); v++)
        if (P.live[v]) { depth[sol.s[v]] = std::max(depth[sol.s[v]], d[v]); size[sol.s[v]]++; }
    printf("  cut widths      :");
    for (int i = 1; i < P.S; i++) printf(" R%d=%zu", i, cuts[i].size());
    printf("\n  registers       : %lld (balancing/PI-delay: %lld)\n", sol.regs, balancing);
    printf("  register area   : %lld (weight %lld)\n", sol.regs * P.W, P.W);
    printf("  stage ANDs      :");
    for (int s = 0; s < P.S; s++) printf(" %d", size[s]);
    printf("\n  stage depth     :");
    for (int s = 0; s < P.S; s++) printf(" %d", depth[s]);
    printf("  (max %d)\n", *std::max_element(depth.begin(), depth.end()));
}

// ---------------------------------------------------------------------------------------------

static void usage(const char* prog) {
    fprintf(stderr,
        "usage: %s [options] <file.aig|file.aag>\n"
        "  --stages S            number of pipeline stages (default 2)\n"
        "  --mode count|area     W = 1 (count) or W = --reg-weight (area, default 10)\n"
        "  --reg-weight W        register area weight in NAND2 units (implies --mode area)\n"
        "  --cut-side source|sink  min-cut nearest the PIs (source) or the POs (sink)\n"
        "  --period T            max AND depth per stage (stage C), 0 = none\n"
        "  --method auto|flow|lp flow: stage A only (S=2, no T); lp: stage B/C;\n"
        "                        auto: lp, plus flow cross-check when S=2 and no T\n"
        "  --out-dir D           where stage_assign.json / cuts.json go (default .)\n"
        "  --no-json             do not write JSON files\n"
        "  --write-blif F        write the pipelined circuit (BLIF with latches)\n"
        "  --write-ref-blif F    write the original logic with PIs delayed by S-1 latches\n"
        "  --brute N             brute-force check when the circuit has <= N live ANDs\n"
        "  --spec-graph          stage A with the original spec 4.1 v_in->u_out arc (for comparison)\n"
        "  --sim N               simulation cycles for latency check (default 32, 0 = off)\n"
        "  -q                    quiet: one summary line\n", prog);
    exit(2);
}

static Options parseArgs(int argc, char** argv) {
    Options o;
    bool area = false, weightSet = false;
    for (int i = 1; i < argc; i++) {
        string a = argv[i];
        auto next = [&]() -> string { if (i + 1 >= argc) usage(argv[0]); return argv[++i]; };
        if (a == "--stages") o.S = std::stoi(next());
        else if (a == "--mode") { string m = next(); if (m == "area") area = true; else if (m == "count") area = false; else usage(argv[0]); }
        else if (a == "--reg-weight") { o.W = std::stoll(next()); weightSet = area = true; }
        else if (a == "--cut-side") o.side = next();
        else if (a == "--period") o.T = std::stoi(next());
        else if (a == "--method") o.method = next();
        else if (a == "--out-dir") o.outDir = next();
        else if (a == "--no-json") o.writeJson = false;
        else if (a == "--write-blif") o.blif = next();
        else if (a == "--write-ref-blif") o.refBlif = next();
        else if (a == "--brute") o.brute = std::stoi(next());
        else if (a == "--spec-graph") o.specGraph = true;
        else if (a == "--sim") o.simCycles = std::stoi(next());
        else if (a == "-q") o.quiet = true;
        else if (a == "-h" || a == "--help") usage(argv[0]);
        else if (!a.empty() && a[0] == '-') { fprintf(stderr, "unknown option %s\n", a.c_str()); usage(argv[0]); }
        else o.input = a;
    }
    if (area && !weightSet) o.W = 10;
    if (!area) o.W = 1;
    if (o.input.empty() || o.S < 1 || o.W < 1 || o.T < 0) usage(argv[0]);
    if (o.side == "balanced") { fprintf(stderr, "--cut-side balanced is not implemented yet\n"); exit(2); }
    if (o.side != "source" && o.side != "sink") usage(argv[0]);
    if (o.method != "auto" && o.method != "flow" && o.method != "lp") usage(argv[0]);
    if (o.method == "flow" && (o.S != 2 || o.T)) { fprintf(stderr, "--method flow needs --stages 2 and no --period\n"); exit(2); }
    return o;
}

int main(int argc, char** argv) {
    Options opt = parseArgs(argc, argv);
    try {
        double t0 = nowSec();
        Aig g = readAiger(opt.input);
        Problem P(g, opt.S, opt.W, opt.T);
        buildProblem(P);
        int nLive = 0, nLiveAnd = 0, maxLevel = 0;
        for (int v = 1; v < g.numNodes(); v++)
            if (P.live[v]) { nLive++; if (g.isAnd(v)) { nLiveAnd++; maxLevel = std::max(maxLevel, P.L[v]); } }
        double tParse = nowSec() - t0;
        if (!opt.quiet)
            printf("%s: %d PI, %d PO, %d AND (%d live), level %d; S=%d W=%lld T=%d side=%s\n",
                   opt.input.c_str(), g.nPi, (int)g.poLit.size(), g.numAnds(), nLiveAnd, maxLevel,
                   opt.S, opt.W, opt.T, opt.side.c_str());

        Solution sol, flowSol;
        ll flowValue = -1;
        bool runFlow = opt.method == "flow" || (opt.method == "auto" && opt.S == 2 && opt.T == 0);
        bool runLp = opt.method != "flow";
        double tFlow = 0, tLp = 0;
        LpInfo info;
        int rounds = 0;
        if (opt.S == 1) {   // single stage: no registers at all
            sol.s.assign(g.numNodes(), -1);
            for (int v = 1; v < g.numNodes(); v++) if (P.live[v]) sol.s[v] = 0;
            computeM(P, sol);
            runFlow = runLp = false;
            if (opt.T > 0 && maxLevel > opt.T) { printf("INFEASIBLE: depth %d > T=%d with 1 stage\n", maxLevel, opt.T); return 3; }
        }
        if (runFlow) {
            double t = nowSec();
            if (!solveFlow(P, opt.side, opt.specGraph, flowSol, flowValue)) throw std::logic_error("unbounded max-flow");
            tFlow = nowSec() - t;
            if (!opt.quiet) printf("[A] max-flow %s: min-cut = %lld (assignment cost %lld), %.3fs\n",
                                   opt.specGraph ? "(spec graph)" : "", flowValue, flowSol.regs * opt.W, tFlow);
            sol = flowSol;
        }
        if (runLp) {
            double t = nowSec();
            if (solveLpTimed(P, opt.side, sol, info, rounds)) {
                printf("INFEASIBLE: no assignment with S=%d and stage depth <= %d (circuit depth %d)\n", opt.S, opt.T, maxLevel);
                if (opt.brute > 0 && nLiveAnd <= opt.brute) {
                    Solution best;
                    if (bruteForce(P, best) >= 0) { printf("  check brute force FAIL: found a feasible assignment\n"); return 1; }
                }
                return 3;
            }
            tLp = nowSec() - t;
            if (!opt.quiet) printf("[B] LP/min-cost-flow: %lld registers, %d mcf phases, %d round(s), %zu lazy constraints, %.3fs\n",
                                   sol.regs, info.phases, rounds, P.lazy.size(), tLp);
        }

        int status = 0;
        auto check = [&](const char* what, bool ok, const string& msg) {
            if (!opt.quiet || !ok) printf("  check %-26s %s%s%s\n", what, ok ? "PASS" : "FAIL", msg.empty() ? "" : ": ", msg.c_str());
            if (!ok) status = 1;
        };
        if (runFlow) {
            string e = verifyStructure(P, flowSol);
            check("A: legal assignment", e.empty(), e);
            check("A: cut value == cost", flowValue == flowSol.regs * opt.W,
                  std::to_string(flowValue) + " vs " + std::to_string(flowSol.regs * opt.W));
            if (runLp) {
                check("A == B (S=2 cross-check)", flowValue == sol.regs * opt.W,
                      std::to_string(flowValue) + " vs " + std::to_string(sol.regs * opt.W));
                check("A, B same extreme cut", flowSol.s == sol.s, "");
            }
        }
        string e = verifyStructure(P, sol);
        check("balanced / legal", e.empty(), e);
        if (opt.simCycles > 0 && e.empty()) {
            string se = simulate(P, sol, opt.simCycles + opt.S, 12345);
            check("latency simulation", se.empty(), se);
        }
        if (opt.brute > 0) {
            if (nLiveAnd <= opt.brute) {
                Solution best;
                ll b = bruteForce(P, best);
                if (b < 0) check("brute force", !runLp, "brute force finds no feasible assignment");
                else check("brute-force optimum", b == sol.regs, std::to_string(b) + " vs " + std::to_string(sol.regs));
            } else if (!opt.quiet) printf("  check brute force skipped (%d live ANDs > %d)\n", nLiveAnd, opt.brute);
        }

        if (!opt.quiet) printReport(P, sol);
        if (opt.writeJson) {
            writeStageAssign(P, sol, opt.outDir + "/stage_assign.json");
            writeCuts(P, sol, opt.outDir + "/cuts.json");
        }
        if (!opt.blif.empty()) writePipelineBlif(P, sol, opt.blif);
        if (!opt.refBlif.empty()) writeReferenceBlif(P, opt.refBlif);
        if (opt.quiet)
            printf("%s S=%d T=%d regs=%lld area=%lld flow=%lld time=%.3fs %s\n", opt.input.c_str(), opt.S, opt.T,
                   sol.regs, sol.regs * opt.W, flowValue, tParse + tFlow + tLp, status ? "FAIL" : "ok");
        else
            printf("  time            : parse %.3fs, flow %.3fs, lp %.3fs\n", tParse, tFlow, tLp);
        return status;
    } catch (const std::exception& ex) {
        fprintf(stderr, "error: %s\n", ex.what());
        return 2;
    }
}
