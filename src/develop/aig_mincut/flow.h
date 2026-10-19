// flow.h -- residual graph with Dinic max-flow and primal-dual min-cost flow.
#pragma once
#include <algorithm>
#include <cstdint>
#include <functional>
#include <queue>
#include <vector>

typedef long long ll;
static const ll FLOW_INF = (ll)1e15;

struct FlowGraph {
    int n;
    std::vector<int> head, nxt, to;
    std::vector<ll> cap, cost;

    explicit FlowGraph(int n_) : n(n_), head(n_, -1) {}

    // Arc a->b; its reverse arc is (index ^ 1).
    int addArc(int a, int b, ll c, ll w = 0) {
        int e = (int)to.size();
        to.push_back(b); cap.push_back(c); cost.push_back(w);  nxt.push_back(head[a]); head[a] = e;
        to.push_back(a); cap.push_back(0); cost.push_back(-w); nxt.push_back(head[b]); head[b] = e + 1;
        return e;
    }
    int tail(int e) const { return to[e ^ 1]; }
    ll  flowOn(int e) const { return cap[e ^ 1]; }   // flow on forward arc e
    int numArcs() const { return (int)to.size(); }

    // Arc usable for augmentation. With pot != nullptr, only zero-reduced-cost arcs.
    bool admissible(int e, int u, const std::vector<ll>* pot) const {
        if (cap[e] <= 0) return false;
        if (!pot) return true;
        return cost[e] + (*pot)[u] - (*pot)[to[e]] == 0;
    }

    // Dinic (iterative DFS, safe for very deep graphs). Returns flow value, or -1 if unbounded.
    ll maxFlow(int s, int t, const std::vector<ll>* pot = nullptr) {
        ll total = 0;
        std::vector<int> level(n), it(n), path;
        std::vector<int> q(n);
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
                level[u] = -1;                       // dead end
                if (u == s) break;
                u = tail(path.back()); path.pop_back();
            }
        }
        return total;
    }

    // Nodes reachable from s in the residual graph.
    std::vector<char> reachFrom(int s) const {
        std::vector<char> r(n, 0);
        std::vector<int> st{ s }; r[s] = 1;
        while (!st.empty()) {
            int u = st.back(); st.pop_back();
            for (int e = head[u]; e != -1; e = nxt[e])
                if (cap[e] > 0 && !r[to[e]]) { r[to[e]] = 1; st.push_back(to[e]); }
        }
        return r;
    }
    // Nodes that can reach t in the residual graph.
    std::vector<char> reachTo(int t) const {
        std::vector<char> r(n, 0);
        std::vector<int> st{ t }; r[t] = 1;
        while (!st.empty()) {
            int v = st.back(); st.pop_back();
            for (int e = head[v]; e != -1; e = nxt[e]) {   // e: v->u, residual arc u->v is e^1
                int u = to[e];
                if (cap[e ^ 1] > 0 && !r[u]) { r[u] = 1; st.push_back(u); }
            }
        }
        return r;
    }

    // Dijkstra on residual arcs with reduced costs (cost + pot[u] - pot[v] >= 0).
    // reverse=true walks arcs backwards (distance *to* src). skip[] nodes are ignored.
    std::vector<ll> dijkstra(int src, const std::vector<ll>& pot, bool reverse = false,
                             const std::vector<char>* skip = nullptr) const {
        std::vector<ll> d(n, FLOW_INF);
        typedef std::pair<ll, int> P;
        std::priority_queue<P, std::vector<P>, std::greater<P>> pq;
        d[src] = 0; pq.push({ 0, src });
        while (!pq.empty()) {
            P p = pq.top(); pq.pop();
            int u = p.second;
            if (p.first != d[u]) continue;
            for (int e = head[u]; e != -1; e = nxt[e]) {
                int v = to[e];
                if (skip && (*skip)[v]) continue;
                int arc = reverse ? (e ^ 1) : e;        // reverse: residual arc v->u
                if (cap[arc] <= 0) continue;
                ll rc = reverse ? cost[arc] + pot[v] - pot[u] : cost[arc] + pot[u] - pot[v];
                if (d[u] + rc < d[v]) { d[v] = d[u] + rc; pq.push({ d[v], v }); }
            }
        }
        return d;
    }

    // Primal-dual min-cost flow from s to t. pot must be a feasible potential
    // (all residual arcs have non-negative reduced cost); it is updated in place and
    // stays feasible. Each phase: Dijkstra, then Dinic on the zero-reduced-cost subgraph.
    struct McfResult { ll flow = 0, cost = 0; int phases = 0; bool unbounded = false; };
    McfResult minCostFlow(int s, int t, std::vector<ll>& pot) {
        McfResult r;
        while (true) {
            std::vector<ll> d = dijkstra(s, pot);
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
};
