// lsvPipeMinCut.cpp -- fixed-stage pipeline register minimization on a combinational GIA.
//
//   Stage A: S = 2, node-split max-flow / min-cut (exact).
//   Stage B: general S, difference-constraint LP solved through its min-cost-flow dual (exact).
//   Stage C: per-stage depth limit T, level windows + lazy path constraints (exact).
// Ported from src/develop/aig_mincut (standalone tool); see src/develop/plan/aig_mincut_spec.md.
#include "lsvPipe.h"

#include <chrono>
#include <random>

namespace lsvpipe {

using std::string;
using std::vector;

static double nowSec() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

// ---------------------------------------------------------------------------------------------
// PipeAig

void PipeAig::buildFanout() {
    int n = numNodes();
    foStart.assign(n + 1, 0);
    for (int v = nPi + 1; v < n; v++) { foStart[(fanin0[v] >> 1) + 1]++; foStart[(fanin1[v] >> 1) + 1]++; }
    for (int i = 0; i < n; i++) foStart[i + 1] += foStart[i];
    foList.assign(foStart[n], 0);
    vector<int> pos(foStart.begin(), foStart.end() - 1);
    for (int v = nPi + 1; v < n; v++) { foList[pos[fanin0[v] >> 1]++] = v; foList[pos[fanin1[v] >> 1]++] = v; }
}

bool pipeAigFromGia(Gia_Man_t* p, PipeAig& g, string& msg) {
    if (Gia_ManRegNum(p) > 0) { msg = "the GIA is sequential; a combinational AIG is required"; return false; }
    Gia_Obj_t* pObj;
    int i;
    vector<int> map(Gia_ManObjNum(p), -1);
    map[0] = 0;
    g = PipeAig();
    g.nPi = Gia_ManCiNum(p);
    g.fanin0.assign(g.nPi + 1, 0); g.fanin1.assign(g.nPi + 1, 0);
    Gia_ManForEachCi(p, pObj, i) {
        map[Gia_ObjId(p, pObj)] = i + 1;
        const char* nm = Gia_ObjCiName(p, i);
        g.piName.push_back(nm ? nm : "pi" + std::to_string(i));
    }
    Gia_ManForEachAnd(p, pObj, i) {
        if (Gia_ObjIsXor(pObj) || Gia_ObjIsMux(p, pObj)) { msg = "the GIA has XOR/MUX nodes; an AIG of ANDs is required"; return false; }
        map[i] = g.numNodes();
        g.fanin0.push_back(2 * map[Gia_ObjFaninId0(pObj, i)] + Gia_ObjFaninC0(pObj));
        g.fanin1.push_back(2 * map[Gia_ObjFaninId1(pObj, i)] + Gia_ObjFaninC1(pObj));
    }
    Gia_ManForEachCo(p, pObj, i) {
        g.poLit.push_back(2 * map[Gia_ObjFaninId0p(p, pObj)] + Gia_ObjFaninC0(pObj));
        const char* nm = Gia_ObjCoName(p, i);
        g.poName.push_back(nm ? nm : "po" + std::to_string(i));
    }
    g.buildFanout();
    return true;
}

void copyNames(Gia_Man_t* pNew, Gia_Man_t* pOrig) {
    ABC_FREE(pNew->pName);   // Gia_ManDup / Gia_ManCleanup copy the name of their input
    if (pOrig->pName) pNew->pName = Abc_UtilStrsav(pOrig->pName);
    if (!pOrig->vNamesIn || !pOrig->vNamesOut) return;
    int nRegs = Gia_ManRegNum(pNew), i;
    pNew->vNamesIn = Vec_PtrAlloc(Gia_ManCiNum(pNew));
    pNew->vNamesOut = Vec_PtrAlloc(Gia_ManCoNum(pNew));
    for (i = 0; i < Gia_ManPiNum(pNew); i++) Vec_PtrPush(pNew->vNamesIn, Abc_UtilStrsav((char*)Vec_PtrEntry(pOrig->vNamesIn, i)));
    for (i = 0; i < nRegs; i++) { char b[32]; sprintf(b, "pipe_r%d", i); Vec_PtrPush(pNew->vNamesIn, Abc_UtilStrsav(b)); }
    for (i = 0; i < Gia_ManPoNum(pNew); i++) Vec_PtrPush(pNew->vNamesOut, Abc_UtilStrsav((char*)Vec_PtrEntry(pOrig->vNamesOut, i)));
    for (i = 0; i < nRegs; i++) { char b[32]; sprintf(b, "pipe_r%d_in", i); Vec_PtrPush(pNew->vNamesOut, Abc_UtilStrsav(b)); }
}

// ---------------------------------------------------------------------------------------------
// FlowGraph

ll FlowGraph::maxFlow(int s, int t, const vector<ll>* pot) {
    ll total = 0;
    vector<int> level(n), it(n), path, q(n);
    while (true) {
        std::fill(level.begin(), level.end(), -1);
        int qh = 0, qt = 0;
        level[s] = 0; q[qt++] = s;
        while (qh < qt) {
            int u = q[qh++];
            for (int e = head[u]; e != -1; e = nxt[e])
                if (level[to[e]] < 0 && admissible(e, u, pot)) { level[to[e]] = level[u] + 1; q[qt++] = to[e]; }
        }
        if (level[t] < 0) break;
        it = head;
        path.clear();
        int u = s;
        while (true) {
            if (u == t) {
                ll f = FLOW_INF;
                for (int e : path) f = std::min(f, cap[e]);
                if (f >= FLOW_INF) return -1;
                for (int e : path) { cap[e] -= f; cap[e ^ 1] += f; }
                total += f;
                size_t k = 0;
                while (cap[path[k]] > 0) k++;
                path.resize(k);
                u = k ? to[path[k - 1]] : s;
                continue;
            }
            bool adv = false;
            for (int& e = it[u]; e != -1; e = nxt[e]) {
                int v = to[e];
                if (level[v] == level[u] + 1 && admissible(e, u, pot)) { path.push_back(e); u = v; adv = true; break; }
            }
            if (adv) continue;
            level[u] = -1;
            if (u == s) break;
            u = tail(path.back()); path.pop_back();
        }
    }
    return total;
}

vector<char> FlowGraph::reachFrom(int s) const {
    vector<char> r(n, 0);
    vector<int> st{ s }; r[s] = 1;
    while (!st.empty()) {
        int u = st.back(); st.pop_back();
        for (int e = head[u]; e != -1; e = nxt[e])
            if (cap[e] > 0 && !r[to[e]]) { r[to[e]] = 1; st.push_back(to[e]); }
    }
    return r;
}

vector<char> FlowGraph::reachTo(int t) const {
    vector<char> r(n, 0);
    vector<int> st{ t }; r[t] = 1;
    while (!st.empty()) {
        int v = st.back(); st.pop_back();
        for (int e = head[v]; e != -1; e = nxt[e]) {
            int u = to[e];
            if (cap[e ^ 1] > 0 && !r[u]) { r[u] = 1; st.push_back(u); }
        }
    }
    return r;
}

vector<ll> FlowGraph::dijkstra(int src, const vector<ll>& pot, bool reverse, const vector<char>* skip) const {
    vector<ll> d(n, FLOW_INF);
    typedef std::pair<ll, int> P;
    std::priority_queue<P, vector<P>, std::greater<P>> pq;
    d[src] = 0; pq.push({ 0, src });
    while (!pq.empty()) {
        P p = pq.top(); pq.pop();
        int u = p.second;
        if (p.first != d[u]) continue;
        for (int e = head[u]; e != -1; e = nxt[e]) {
            int v = to[e];
            if (skip && (*skip)[v]) continue;
            int arc = reverse ? (e ^ 1) : e;
            if (cap[arc] <= 0) continue;
            ll rc = reverse ? cost[arc] + pot[v] - pot[u] : cost[arc] + pot[u] - pot[v];
            if (d[u] + rc < d[v]) { d[v] = d[u] + rc; pq.push({ d[v], v }); }
        }
    }
    return d;
}

FlowGraph::McfResult FlowGraph::minCostFlow(int s, int t, vector<ll>& pot) {
    McfResult r;
    while (true) {
        vector<ll> d = dijkstra(s, pot);
        if (d[t] >= FLOW_INF) break;
        for (int v = 0; v < n; v++) pot[v] += std::min(d[v], d[t]);
        ll f = maxFlow(s, t, &pot);
        if (f < 0) { r.unbounded = true; return r; }
        if (f == 0) break;
        r.flow += f;
        r.phases++;
    }
    for (int e = 0; e < numArcs(); e += 2) r.cost += cost[e] * flowOn(e);
    return r;
}

// ---------------------------------------------------------------------------------------------
// Problem setup

void buildProblem(Problem& P) {
    const PipeAig& g = P.g;
    int n = g.numNodes();
    P.live.assign(n, 0); P.poDrv.assign(n, 0);
    for (int lit : g.poLit) if (lit >> 1) P.poDrv[lit >> 1] = P.live[lit >> 1] = 1;
    for (int v = n - 1; v > g.nPi; v--)
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
    if (P.T > 0)
        for (int v = g.nPi + 1; v < n; v++) {
            P.lo[v] = (P.L[v] + P.T - 1) / P.T - 1;
            P.hi[v] = P.S - (P.R[v] + 1 + P.T - 1) / P.T;
        }
}

void computeM(const Problem& P, Solution& sol) {
    const PipeAig& g = P.g;
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

vector<int> stageDepth(const Problem& P, const vector<int>& s, vector<int>* pred) {
    const PipeAig& g = P.g;
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

// Longest same-stage path (in ANDs) from each live AND down through its fanouts, counting itself.
vector<int> stageDown(const Problem& P, const vector<int>& s) {
    const PipeAig& g = P.g;
    int n = g.numNodes();
    vector<int> down(n, 0);
    for (int w = n - 1; w > g.nPi; w--) {
        if (!P.live[w]) continue;
        down[w] = 1;
        for (int k = g.foStart[w]; k < g.foStart[w + 1]; k++) {
            int x = g.foList[k];
            if (P.live[x] && s[x] == s[w]) down[w] = std::max(down[w], down[x] + 1);
        }
    }
    return down;
}

// ---------------------------------------------------------------------------------------------
// Stage A: S = 2 via node-split min-cut. Source side X = stage 0; v is registered iff
// v_in in X and v_out not in X. Closure arc v_in -> u_in (not u_out: a registered node may also
// feed stage-0 logic directly).

static bool solveFlow(const Problem& P, bool sinkSide, Solution& sol, ll& cutValue) {
    const PipeAig& g = P.g;
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
        if (g.isAnd(v))
            for (int lit : { g.fanin0[v], g.fanin1[v] }) {
                int u = lit >> 1;
                if (u == 0) continue;
                G.addArc(vout(u), vin(v), FLOW_INF);
                G.addArc(vin(v), vin(u), FLOW_INF);
            }
    }
    cutValue = G.maxFlow(src, snk);
    if (cutValue < 0) return false;
    vector<char> X;
    if (sinkSide) { X = G.reachTo(snk); for (auto& x : X) x = !x; }
    else X = G.reachFrom(src);
    sol.s.assign(n, -1);
    for (int v = 1; v < n; v++) if (P.live[v]) sol.s[v] = X[vin(v)] ? 0 : 1;
    computeM(P, sol);
    return true;
}

// ---------------------------------------------------------------------------------------------
// Stage B: LP over s(v), m(u) as the dual of a min-cost flow. Constraint x_b - x_a <= c is an
// uncapacitated arc a->b of cost c; supply W at m(u), demand W at s(u). Shortest distances from
// the reference node z in the final residual graph give the pointwise-maximal optimum (cut
// nearest the PIs), distances to z the minimal one.
// Returns 0 = ok, 1 = infeasible, 2 = internal error (err set).

static int solveLp(const Problem& P, bool sinkSide, Solution& sol, int& phases, string& err) {
    const PipeAig& g = P.g;
    int n = g.numNodes();
    vector<int> idx(n, -1), nodes;
    for (int v = 1; v < n; v++) if (P.live[v]) { idx[v] = (int)nodes.size(); nodes.push_back(v); }
    int K = (int)nodes.size();
    auto sI = [&](int v) { return 1 + 2 * idx[v]; };
    auto mI = [&](int v) { return 2 + 2 * idx[v]; };
    const int z = 0, SS = 1 + 2 * K, ST = 2 + 2 * K, S1 = P.S - 1;

    vector<vector<int>> lazyPred(n);
    for (auto& c : P.lazy) lazyPred[c.second].push_back(c.first);
    vector<ll> pot(2 * K + 3, 0);
    vector<int> s0(n, 0);
    for (int v : nodes) {           // ASAP schedule: initial feasible potential + exact feasibility
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
    auto le = [&](int a, int b, ll c) { G.addArc(a, b, FLOW_INF, c); };
    for (int v : nodes) {
        le(z, sI(v), P.hi[v]);
        le(sI(v), z, -P.lo[v]);
        le(mI(v), sI(v), 0);
        le(z, mI(v), S1);
        if (P.poDrv[v]) le(mI(v), z, -S1);
        if (g.isAnd(v))
            for (int lit : { g.fanin0[v], g.fanin1[v] }) {
                int u = lit >> 1;
                if (u == 0) continue;
                le(sI(v), sI(u), 0);
                le(mI(u), sI(v), 0);
            }
        G.addArc(SS, mI(v), P.W);
        G.addArc(sI(v), ST, P.W);
    }
    for (auto& c : P.lazy) le(sI(c.second), sI(c.first), -1);
    for (int e = 0; e < G.numArcs(); e++)
        if (G.cap[e] > 0 && G.cost[e] + pot[G.tail(e)] - pot[G.to[e]] < 0) { err = "initial potential is not feasible"; return 2; }

    FlowGraph::McfResult r = G.minCostFlow(SS, ST, pot);
    if (r.unbounded || r.flow != (ll)K * P.W) { err = "min-cost flow did not route all supply"; return 2; }
    phases = r.phases;
    vector<char> skip(G.n, 0);
    skip[SS] = skip[ST] = 1;
    vector<ll> d = G.dijkstra(z, pot, sinkSide, &skip);
    sol.s.assign(n, -1);
    vector<int> mLp(n, -1);
    for (int v : nodes) {
        for (int i : { sI(v), mI(v) }) if (d[i] >= FLOW_INF) { err = "LP node unreachable"; return 2; }
        auto x = [&](int i) { return (int)(sinkSide ? -(d[i] - pot[i] + pot[z]) : d[i] + pot[i] - pot[z]); };
        sol.s[v] = x(sI(v)); mLp[v] = x(mI(v));
    }
    computeM(P, sol);
    for (int v : nodes) if (mLp[v] != sol.m[v]) { err = "LP m(u) is not tight"; return 2; }
    if (sol.regs * P.W != -r.cost) { err = "primal/dual objective mismatch"; return 2; }
    return 0;
}

// Stage C: lazily add s(v) >= s(u) + 1 for every same-stage path of T+1 ANDs. Each such
// constraint is necessary, so the final solution is optimal and infeasibility is exact.
static int solveLpTimed(Problem& P, bool sinkSide, Solution& sol, int& phases, int& rounds, string& err) {
    rounds = 0;
    while (true) {
        rounds++;
        int r = solveLp(P, sinkSide, sol, phases, err);
        if (r || P.T <= 0) return r;
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
// Brute force and verification

static ll bruteForce(const Problem& P) {
    const PipeAig& g = P.g;
    int n = g.numNodes();
    vector<int> ands;
    for (int v = g.nPi + 1; v < n; v++) if (P.live[v]) ands.push_back(v);
    Solution cur;
    cur.s.assign(n, -1);
    for (int v = 1; v <= g.nPi; v++) if (P.live[v]) cur.s[v] = 0;
    vector<int> d(n, 0);
    ll best = -1;
    std::function<void(size_t)> rec = [&](size_t i) {
        if (i == ands.size()) {
            computeM(P, cur);
            if (best < 0 || cur.regs < best) best = cur.regs;
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
    return best;
}

static string verifyStructure(const Problem& P, const Solution& sol) {
    const PipeAig& g = P.g;
    int n = g.numNodes();
    vector<int> regsTo(n, 0);
    for (int v = 1; v < n; v++) {
        if (!P.live[v]) continue;
        int s = sol.s[v];
        if (s < 0 || s > P.S - 1) return "stage out of range at node " + std::to_string(v);
        if (g.isPi(v) && s != 0) return "PI not in stage 0";
        if (sol.m[v] < s) return "negative chain length at node " + std::to_string(v);
        if (!g.isAnd(v)) continue;
        int r = -1;
        for (int lit : { g.fanin0[v], g.fanin1[v] }) {
            int u = lit >> 1;
            if (u == 0) continue;
            if (s < sol.s[u]) return "edge goes backwards in stage at node " + std::to_string(v);
            if (s > sol.m[u]) return "chain too short at node " + std::to_string(v);
            int ru = regsTo[u] + (s - sol.s[u]);
            if (r >= 0 && ru != r) return "unbalanced paths into node " + std::to_string(v);
            r = ru;
        }
        regsTo[v] = r < 0 ? s : r;
    }
    for (int lit : g.poLit) {
        int u = lit >> 1;
        if (!u) continue;
        if (sol.m[u] < P.S - 1) return "PO driver chain too short";
        if (regsTo[u] + (P.S - 1 - sol.s[u]) != P.S - 1) return "PO path does not cross S-1 registers";
    }
    if (P.T > 0) {
        vector<int> d = stageDepth(P, sol.s);
        for (int v = 1; v < n; v++) if (P.live[v] && d[v] > P.T) return "stage depth exceeds T";
    }
    return "";
}

static vector<char> zeroValues(const PipeAig& g) {
    vector<char> z(g.numNodes(), 0);
    for (int v = g.nPi + 1; v < g.numNodes(); v++)
        z[v] = (z[g.fanin0[v] >> 1] ^ (g.fanin0[v] & 1)) & (z[g.fanin1[v] >> 1] ^ (g.fanin1[v] & 1));
    return z;
}

// Cycle-accurate 64-way simulation of the pipeline vs. the original fed with inputs delayed S-1.
static string simulate(const Problem& P, const Solution& sol, int cycles) {
    const PipeAig& g = P.g;
    int n = g.numNodes(), S = P.S;
    typedef unsigned long long u64;
    std::mt19937_64 rng(12345);
    vector<char> z0 = zeroValues(g);
    vector<vector<u64>> hist(S, vector<u64>(n)), inHist(S, vector<u64>(g.nPi + 1, 0));
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
        int tr = t - (S - 1);
        ref[0] = 0;
        for (int v = 1; v <= g.nPi; v++) ref[v] = tr >= 0 ? inHist[tr % S][v] : 0;
        for (int v = g.nPi + 1; v < n; v++) ref[v] = lit(ref[g.fanin0[v] >> 1], g.fanin0[v]) & lit(ref[g.fanin1[v] >> 1], g.fanin1[v]);
        for (size_t o = 0; o < g.poLit.size(); o++) {
            int l = g.poLit[o], u = l >> 1, back = u ? (S - 1 - sol.s[u]) : 0;
            u64 x = (t - back >= 0) ? hist[(t - back) % S][u] : (z0[u] ? ~0ULL : 0);
            if (lit(x, l) != lit(ref[u], l)) return "mismatch at cycle " + std::to_string(t) + " on PO " + g.poName[o];
        }
    }
    return "";
}

// ---------------------------------------------------------------------------------------------

MinCutResult solveMinCut(Problem& P, const MinCutOptions& opt, Solution& sol) {
    MinCutResult res;
    const PipeAig& g = P.g;
    int nLiveAnd = 0;
    for (int v = g.nPi + 1; v < g.numNodes(); v++) nLiveAnd += P.live[v];
    bool runFlow = opt.method == 1 || (opt.method == 0 && P.S == 2 && P.T == 0);
    bool runLp = opt.method != 1;
    Solution flowSol;
    if (P.S == 1) {
        sol.s.assign(g.numNodes(), -1);
        for (int v = 1; v < g.numNodes(); v++) if (P.live[v]) sol.s[v] = 0;
        computeM(P, sol);
        if (P.T > 0 && !verifyStructure(P, sol).empty()) res.infeasible = true;
        return res;
    }
    if (runFlow) {
        double t = nowSec();
        if (!solveFlow(P, opt.sinkSide, flowSol, res.flowValue)) { res.error = "unbounded max-flow"; res.ok = false; return res; }
        res.tFlow = nowSec() - t;
        sol = flowSol;
    }
    if (runLp) {
        double t = nowSec();
        int r = solveLpTimed(P, opt.sinkSide, sol, res.lpPhases, res.lpRounds, res.error);
        res.tLp = nowSec() - t;
        res.lazy = (int)P.lazy.size();
        if (r == 1) { res.infeasible = true; return res; }
        if (r == 2) { res.ok = false; return res; }
    }
    auto check = [&](const char* what, bool good, const string& msg) {
        if (opt.verbose || !good) printf("  check %-26s %s%s%s\n", what, good ? "PASS" : "FAIL", msg.empty() ? "" : ": ", msg.c_str());
        if (!good) res.ok = false;
    };
    if (runFlow) {
        string e = verifyStructure(P, flowSol);
        check("A: legal assignment", e.empty(), e);
        check("A: cut value == cost", res.flowValue == flowSol.regs * P.W, "");
        if (runLp) {
            check("A == B (S=2 cross-check)", res.flowValue == sol.regs * P.W,
                  std::to_string(res.flowValue) + " vs " + std::to_string(sol.regs * P.W));
            check("A, B same extreme cut", flowSol.s == sol.s, "");
        }
    }
    string e = verifyStructure(P, sol);
    check("balanced / legal", e.empty(), e);
    if (opt.simCycles > 0 && e.empty()) {
        string se = simulate(P, sol, opt.simCycles + P.S);
        check("latency simulation", se.empty(), se);
    }
    if (opt.brute > 0) {
        if (nLiveAnd <= opt.brute) {
            ll b = bruteForce(P);
            check("brute-force optimum", b == sol.regs, std::to_string(b) + " vs " + std::to_string(sol.regs));
        } else if (opt.verbose)
            printf("  check brute force skipped (%d live ANDs > %d)\n", nLiveAnd, opt.brute);
    }
    return res;
}

vector<vector<int>> cutSets(const Problem& P, const Solution& sol) {
    vector<vector<int>> cuts(P.S);
    for (int u = 1; u < P.g.numNodes(); u++)
        if (P.live[u])
            for (int i = sol.s[u] + 1; i <= sol.m[u]; i++) cuts[i].push_back(u);
    return cuts;
}

void printReport(const Problem& P, const Solution& sol) {
    const PipeAig& g = P.g;
    vector<vector<int>> cuts = cutSets(P, sol);
    ll balancing = 0;
    for (int v = 1; v <= g.nPi; v++) if (P.live[v]) balancing += sol.m[v] - sol.s[v];
    vector<int> d = stageDepth(P, sol.s);
    vector<int> depth(P.S, 0), size(P.S, 0);
    for (int v = g.nPi + 1; v < g.numNodes(); v++)
        if (P.live[v]) { depth[sol.s[v]] = std::max(depth[sol.s[v]], d[v]); size[sol.s[v]]++; }
    printf("  cut widths      :");
    for (int i = 1; i < P.S; i++) printf(" R%d=%d", i, (int)cuts[i].size());
    printf("\n  registers       : %lld (balancing/PI-delay: %lld)\n", sol.regs, balancing);
    printf("  register area   : %lld (weight %lld)\n", sol.regs * P.W, P.W);
    printf("  stage ANDs      :");
    for (int s = 0; s < P.S; s++) printf(" %d", size[s]);
    printf("\n  stage depth     :");
    for (int s = 0; s < P.S; s++) printf(" %d", depth[s]);
    printf("  (max %d)\n", *std::max_element(depth.begin(), depth.end()));
}

bool writeJson(const Problem& P, const Solution& sol, const string& dir) {
    const PipeAig& g = P.g;
    int n = g.numNodes();
    FILE* f = fopen((dir + "/stage_assign.json").c_str(), "w");
    if (!f) return false;
    fprintf(f, "{\n  \"stages\": %d,\n  \"nodes\": [", P.S);
    bool first = true;
    for (int v = 1; v < n; v++) {
        if (!P.live[v]) continue;
        fprintf(f, "%s    {\"id\": %d, \"type\": \"%s\", \"s\": %d, \"m\": %d, \"po_driver\": %s}", first ? "\n" : ",\n",
                v, g.isPi(v) ? "pi" : "and", sol.s[v], sol.m[v], P.poDrv[v] ? "true" : "false");
        first = false;
    }
    fprintf(f, "\n  ]\n}\n");
    fclose(f);

    f = fopen((dir + "/cuts.json").c_str(), "w");
    if (!f) return false;
    vector<vector<int>> cuts = cutSets(P, sol);
    vector<int> arrival = stageDepth(P, sol.s), down = stageDown(P, sol.s);
    fprintf(f, "{\n  \"stages\": %d,\n  \"reg_weight\": %lld,\n  \"cuts\": [", P.S, P.W);
    for (int i = 1; i < P.S; i++) {
        fprintf(f, "%s    {\"boundary\": %d, \"width\": %d, \"signals\": [", i > 1 ? ",\n" : "\n", i, (int)cuts[i].size());
        first = true;
        for (int u : cuts[i]) {
            int downstream = 0, same = 0, prev = -1;
            for (int k = g.foStart[u]; k < g.foStart[u + 1]; k++) {
                int w = g.foList[k];
                if (!P.live[w] || w == prev) continue;
                prev = w;
                if (sol.s[w] == i) downstream = std::max(downstream, down[w]);
                if (sol.s[w] == sol.s[u]) same++;
            }
            fprintf(f, "%s      {\"node\": %d, \"src_stage\": %d, \"chain_pos\": %d, \"chain_len\": %d, "
                       "\"balancing\": %s, \"po_driver\": %s, \"level\": %d, \"arrival\": %d, \"downstream\": %d, "
                       "\"level_to_po\": %d, \"same_stage_uses\": %d}",
                    first ? "\n" : ",\n", u, sol.s[u], i - sol.s[u], sol.m[u] - sol.s[u], g.isPi(u) ? "true" : "false",
                    P.poDrv[u] ? "true" : "false", P.L[u], g.isAnd(u) ? arrival[u] : 0, downstream, P.R[u], same);
            first = false;
        }
        fprintf(f, "%s", cuts[i].empty() ? "]}" : "\n    ]}");
    }
    fprintf(f, "\n  ]\n}\n");
    fclose(f);
    return true;
}

Gia_Man_t* buildPipelineGia(const Problem& P, const Solution& sol, Gia_Man_t* pOrig) {
    const PipeAig& g = P.g;
    int n = g.numNodes();
    vector<char> z0 = zeroValues(g);
    vector<int> regBase(n, -1);   // index of the first register of u's chain
    int nRegs = 0;
    for (int u = 1; u < n; u++)
        if (P.live[u] && sol.m[u] > sol.s[u]) { regBase[u] = nRegs; nRegs += sol.m[u] - sol.s[u]; }
    Gia_Man_t* pNew = Gia_ManStart(Gia_ManObjNum(pOrig) + 2 * nRegs + 1);
    Gia_ManHashAlloc(pNew);
    vector<int> lit(n, 0), ro(nRegs);
    for (int i = 0; i < g.nPi; i++) lit[i + 1] = Gia_ManAppendCi(pNew);
    for (int r = 0; r < nRegs; r++) ro[r] = Gia_ManAppendCi(pNew);
    // value of u delayed by k registers (registers hold the complement when the reset value is 1)
    auto sig = [&](int u, int k) {
        if (u == 0 || k == 0) return lit[u];
        return Abc_LitNotCond(ro[regBase[u] + k - 1], z0[u]);
    };
    for (int v = g.nPi + 1; v < n; v++) {
        if (!P.live[v]) continue;
        int a = g.fanin0[v], b = g.fanin1[v];
        int la = Abc_LitNotCond(sig(a >> 1, (a >> 1) ? sol.s[v] - sol.s[a >> 1] : 0), a & 1);
        int lb = Abc_LitNotCond(sig(b >> 1, (b >> 1) ? sol.s[v] - sol.s[b >> 1] : 0), b & 1);
        lit[v] = Gia_ManHashAnd(pNew, la, lb);
    }
    for (int l : g.poLit) {
        int u = l >> 1;
        Gia_ManAppendCo(pNew, Abc_LitNotCond(sig(u, u ? P.S - 1 - sol.s[u] : 0), l & 1));
    }
    for (int u = 1; u < n; u++)
        for (int k = 1; regBase[u] >= 0 && k <= sol.m[u] - sol.s[u]; k++)
            Gia_ManAppendCo(pNew, Abc_LitNotCond(sig(u, k - 1), z0[u]));
    Gia_ManHashStop(pNew);
    Gia_ManSetRegNum(pNew, nRegs);
    copyNames(pNew, pOrig);
    return pNew;
}

}  // namespace lsvpipe
