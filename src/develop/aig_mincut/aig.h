// aig.h -- combinational AIG (AIGER .aag / .aig reader) with CSR fanout.
#pragma once
#include <string>
#include <vector>

// Node 0 is constant false, nodes 1..nPi are PIs, nodes nPi+1..n-1 are ANDs in
// topological order. A literal is 2*node + complement.
struct Aig {
    int nPi = 0;
    std::vector<int> fanin0, fanin1;       // literals, valid for AND nodes
    std::vector<int> poLit;                // literal driving each PO
    std::vector<int> origVar;              // node -> AIGER variable index
    std::vector<std::string> piName, poName;
    // CSR fanout (node -> AND fanout nodes; a node appears twice if used by both fanins)
    std::vector<int> foStart, foList;

    int numNodes() const { return (int)fanin0.size(); }
    bool isConst(int v) const { return v == 0; }
    bool isPi(int v) const { return v >= 1 && v <= nPi; }
    bool isAnd(int v) const { return v > nPi; }
    int numAnds() const { return numNodes() - nPi - 1; }

    void buildFanout();
};

// Throws std::runtime_error on malformed input or sequential AIGs.
Aig readAiger(const std::string& path);
