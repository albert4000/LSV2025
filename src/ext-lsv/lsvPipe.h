// lsvPipe.h -- pipeline register minimization (min-cut) and ACD cut re-encoding on ABC's GIA.
//
// Commands (lsvPipeCmd.cpp):
//   lsv_pipe      fixed-stage min-register pipelining (src/develop/plan/aig_mincut_spec.md)
//   lsv_pipe_acd  ACD re-encoding of the 2-stage cut     (src/develop/plan/acd_cut_reencoding.md)
//
// Everything works on a private, flat copy of a combinational GIA (PipeAig): node 0 is const 0,
// nodes 1..nPi are the PIs, the rest are ANDs in topological order; literal = 2*node + compl.
// ABC is compiled with -fno-exceptions, so errors are reported through return values.
#ifndef LSV_PIPE_H
#define LSV_PIPE_H

#include <algorithm>
#include <cstdio>
#include <functional>
#include <queue>
#include <string>
#include <utility>
#include <vector>

#include "aig/gia/gia.h"
#include "base/main/main.h"

namespace lsvpipe {

typedef long long ll;
static const ll FLOW_INF = (ll)1e15;

// ---------------------------------------------------------------------------------------------
// Flat AIG

struct PipeAig {
    int nPi = 0;
    std::vector<int> fanin0, fanin1;      // literals, valid for AND nodes
    std::vector<int> poLit;
    std::vector<std::string> piName, poName;
    std::vector<int> foStart, foList;     // fanout (AND fanouts only)

    int numNodes() const { return (int)fanin0.size(); }
    bool isPi(int v) const { return v >= 1 && v <= nPi; }
    bool isAnd(int v) const { return v > nPi; }
    int numAnds() const { return numNodes() - nPi - 1; }
    void buildFanout();
};

// Returns false (with msg) if the GIA is sequential.
bool pipeAigFromGia(Gia_Man_t* p, PipeAig& g, std::string& msg);

// ---------------------------------------------------------------------------------------------
// Residual graph: Dinic max-flow and primal-dual min-cost flow

struct FlowGraph {
    int n;
    std::vector<int> head, nxt, to;
    std::vector<ll> cap, cost;

    explicit FlowGraph(int n_) : n(n_), head(n_, -1) {}
    int addArc(int a, int b, ll c, ll w = 0) {
        int e = (int)to.size();
        to.push_back(b); cap.push_back(c); cost.push_back(w);  nxt.push_back(head[a]); head[a] = e;
        to.push_back(a); cap.push_back(0); cost.push_back(-w); nxt.push_back(head[b]); head[b] = e + 1;
        return e;
    }
    int tail(int e) const { return to[e ^ 1]; }
    ll flowOn(int e) const { return cap[e ^ 1]; }
    int numArcs() const { return (int)to.size(); }
    bool admissible(int e, int u, const std::vector<ll>* pot) const {
        if (cap[e] <= 0) return false;
        return !pot || cost[e] + (*pot)[u] - (*pot)[to[e]] == 0;
    }
    ll maxFlow(int s, int t, const std::vector<ll>* pot = nullptr);   // -1 if unbounded
    std::vector<char> reachFrom(int s) const;
    std::vector<char> reachTo(int t) const;
    std::vector<ll> dijkstra(int src, const std::vector<ll>& pot, bool reverse = false,
                             const std::vector<char>* skip = nullptr) const;
    struct McfResult { ll flow = 0, cost = 0; int phases = 0; bool unbounded = false; };
    McfResult minCostFlow(int s, int t, std::vector<ll>& pot);
};

// ---------------------------------------------------------------------------------------------
// Min-cut pipelining (stages A/B/C)

struct Problem {
    const PipeAig& g;
    int S;
    ll W;
    int T;
    std::vector<char> live, poDrv;
    std::vector<int> L, R;                     // level from the PIs / ANDs on the longest path to a PO
    std::vector<int> lo, hi;                   // stage window
    std::vector<std::pair<int, int>> lazy;     // (u, v): s(v) >= s(u) + 1
    Problem(const PipeAig& g_, int S_, ll W_, int T_) : g(g_), S(S_), W(W_), T(T_) {}
};

struct Solution {
    std::vector<int> s, m;   // stage and last stage that reads the node; -1 for unused nodes
    ll regs = 0;             // sum (m - s)
};

struct MinCutOptions {
    int S = 2, T = 0;
    ll W = 1;
    bool sinkSide = false;    // cut nearest the POs instead of the PIs
    int method = 0;           // 0 auto (LP, + max-flow cross-check if S=2 and no T), 1 flow only, 2 LP only
    int brute = 0;            // brute-force check if <= brute live ANDs
    int simCycles = 32;
    bool verbose = false;
};

struct MinCutResult {
    bool infeasible = false;
    bool ok = true;            // all checks passed
    std::string error;         // internal error (solver)
    ll flowValue = -1;
    int lpPhases = 0, lpRounds = 0, lazy = 0;
    double tFlow = 0, tLp = 0;
};

void buildProblem(Problem& P);
void computeM(const Problem& P, Solution& sol);
std::vector<int> stageDepth(const Problem& P, const std::vector<int>& s, std::vector<int>* pred = nullptr);
std::vector<int> stageDown(const Problem& P, const std::vector<int>& s);
// Runs stage A and/or B/C, then all checks. Prints check lines when verbose or on failure.
MinCutResult solveMinCut(Problem& P, const MinCutOptions& opt, Solution& sol);
std::vector<std::vector<int>> cutSets(const Problem& P, const Solution& sol);
void printReport(const Problem& P, const Solution& sol);
bool writeJson(const Problem& P, const Solution& sol, const std::string& dir);
// Sequential GIA of the pipeline (flops reset to the node value under all-zero inputs, realized
// with init-0 flops by complementing where needed). Equivalent to the original with every PI
// delayed by S-1 cycles.
Gia_Man_t* buildPipelineGia(const Problem& P, const Solution& sol, Gia_Man_t* pOrig);

// ---------------------------------------------------------------------------------------------
// ACD re-encoding of a 2-stage cut

struct AcdOptions {
    int K = 16;
    bool slackGrouping = false;   // group only within slack classes A/B
    bool ignoreDelay = false;     // accept candidates on area alone
    bool analysisOnly = false;    // phase 1 only
    bool global = false;          // global sampling metric
    int logSamples = 20;
    ll regWeight = 10;
    bool verbose = false;
};

// Returns the pipelined sequential GIA (caller owns), or nullptr if nothing to return
// (analysis only or error, msg set).
Gia_Man_t* runAcd(Problem& P, const Solution& sol, const AcdOptions& opt, Gia_Man_t* pOrig, std::string& msg);

// helpers shared by both parts
void copyNames(Gia_Man_t* pNew, Gia_Man_t* pOrig);

// command helpers (lsvPipeCmd.cpp), also used by lsv_reencode
Gia_Man_t* currentGia(Abc_Frame_t* pAbc, bool fromNtk);
bool parseInt(int argc, char** argv, int& value, int minValue, const char* flag);

}  // namespace lsvpipe

#endif
