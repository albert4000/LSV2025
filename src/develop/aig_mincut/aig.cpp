// aig.cpp -- AIGER reader (ASCII and binary, combinational only).
#include "aig.h"

#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>

void Aig::buildFanout() {
    int n = numNodes();
    foStart.assign(n + 1, 0);
    for (int v = nPi + 1; v < n; v++) { foStart[(fanin0[v] >> 1) + 1]++; foStart[(fanin1[v] >> 1) + 1]++; }
    for (int i = 0; i < n; i++) foStart[i + 1] += foStart[i];
    foList.assign(foStart[n], 0);
    std::vector<int> pos(foStart.begin(), foStart.end() - 1);
    for (int v = nPi + 1; v < n; v++) { foList[pos[fanin0[v] >> 1]++] = v; foList[pos[fanin1[v] >> 1]++] = v; }
}

namespace {

struct Reader {
    const std::string& buf;
    size_t p = 0;
    explicit Reader(const std::string& b) : buf(b) {}
    void fail(const std::string& msg) { throw std::runtime_error("AIGER parse error: " + msg); }
    unsigned readUInt() {
        while (p < buf.size() && (buf[p] == ' ' || buf[p] == '\t')) p++;
        if (p >= buf.size() || buf[p] < '0' || buf[p] > '9') fail("expected number");
        unsigned long long x = 0;
        while (p < buf.size() && buf[p] >= '0' && buf[p] <= '9') x = x * 10 + (buf[p++] - '0');
        return (unsigned)x;
    }
    bool atEol() {
        while (p < buf.size() && (buf[p] == ' ' || buf[p] == '\t' || buf[p] == '\r')) p++;
        return p >= buf.size() || buf[p] == '\n';
    }
    void skipEol() { if (!atEol()) fail("expected end of line"); if (p < buf.size()) p++; }
    unsigned readVarint() {
        unsigned x = 0, i = 0;
        while (true) {
            if (p >= buf.size()) fail("truncated binary AND section");
            unsigned char ch = (unsigned char)buf[p++];
            x |= (unsigned)(ch & 0x7f) << (7 * i++);
            if (!(ch & 0x80)) return x;
        }
    }
};

}  // namespace

Aig readAiger(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("cannot open " + path);
    std::stringstream ss; ss << in.rdbuf();
    std::string buf = ss.str();
    Reader r(buf);

    bool binary;
    if (buf.compare(0, 4, "aag ") == 0) binary = false;
    else if (buf.compare(0, 4, "aig ") == 0) binary = true;
    else throw std::runtime_error("not an AIGER file: " + path);
    r.p = 4;
    unsigned M = r.readUInt(), I = r.readUInt(), L = r.readUInt(), O = r.readUInt(), A = r.readUInt();
    unsigned extra[4] = { 0, 0, 0, 0 };
    for (int k = 0; k < 4 && !r.atEol(); k++) extra[k] = r.readUInt();
    r.skipEol();
    if (L) throw std::runtime_error("sequential AIG (L > 0) is not supported");
    if (extra[0] || extra[1] || extra[2] || extra[3]) throw std::runtime_error("AIGER 1.9 B/C/J/F sections are not supported");
    if (M < I + A) throw std::runtime_error("bad AIGER header (M < I + A)");

    // varDef[var]: -1 undefined, 0 const, 1 PI, 2 AND
    std::vector<int> varDef(M + 1, -1);
    std::vector<unsigned> andL(M + 1, 0), andR(M + 1, 0), piVars;
    varDef[0] = 0;
    for (unsigned i = 0; i < I; i++) {
        unsigned lit = binary ? 2 * (i + 1) : r.readUInt();
        if (!binary) r.skipEol();
        if ((lit & 1) || lit / 2 > M || varDef[lit / 2] != -1) r.fail("bad input literal");
        varDef[lit / 2] = 1; piVars.push_back(lit / 2);
    }
    std::vector<unsigned> outs(O);
    for (unsigned i = 0; i < O; i++) { outs[i] = r.readUInt(); r.skipEol(); if (outs[i] / 2 > M) r.fail("bad output literal"); }
    for (unsigned i = 0; i < A; i++) {
        unsigned lhs, r0, r1;
        if (binary) {
            lhs = 2 * (I + i + 1);
            unsigned d0 = r.readVarint(), d1 = r.readVarint();
            if (d0 > lhs || d1 > lhs - d0) r.fail("bad delta");
            r0 = lhs - d0; r1 = r0 - d1;
        } else {
            lhs = r.readUInt(); r0 = r.readUInt(); r1 = r.readUInt(); r.skipEol();
        }
        if ((lhs & 1) || lhs / 2 > M || varDef[lhs / 2] != -1 || r0 / 2 > M || r1 / 2 > M) r.fail("bad AND definition");
        varDef[lhs / 2] = 2; andL[lhs / 2] = r0; andR[lhs / 2] = r1;
    }

    // Topological renumbering (ASCII files may define ANDs in any order).
    Aig g;
    g.nPi = (int)I;
    std::vector<int> node(M + 1, -1);
    g.origVar.push_back(0); node[0] = 0;
    for (unsigned v : piVars) { node[v] = (int)g.origVar.size(); g.origVar.push_back(v); }
    g.fanin0.assign(I + 1, 0); g.fanin1.assign(I + 1, 0);
    auto mapLit = [&](unsigned lit) { return 2 * node[lit / 2] + (int)(lit & 1); };
    std::vector<char> state(M + 1, 0);   // 1 = on stack, 2 = done
    for (unsigned root = 1; root <= M; root++) {
        if (varDef[root] != 2 || node[root] >= 0) continue;
        std::vector<unsigned> st{ root };
        while (!st.empty()) {
            unsigned v = st.back();
            if (node[v] >= 0) { st.pop_back(); continue; }
            unsigned f[2] = { andL[v] / 2, andR[v] / 2 };
            bool ready = true;
            for (unsigned c : f) {
                if (varDef[c] == -1) r.fail("undefined literal used as AND fanin");
                if (varDef[c] == 2 && node[c] < 0) {
                    if (state[c] == 1) r.fail("combinational cycle");
                    ready = false; st.push_back(c);
                }
            }
            if (ready) {
                node[v] = (int)g.origVar.size(); g.origVar.push_back(v);
                g.fanin0.push_back(mapLit(andL[v])); g.fanin1.push_back(mapLit(andR[v]));
                st.pop_back(); state[v] = 2;
            } else state[v] = 1;
        }
    }
    for (unsigned lit : outs) {
        if (varDef[lit / 2] == -1) r.fail("undefined output literal");
        g.poLit.push_back(mapLit(lit));
    }

    // Symbol table.
    g.piName.resize(I); g.poName.resize(O);
    for (unsigned i = 0; i < I; i++) g.piName[i] = "pi" + std::to_string(i);
    for (unsigned i = 0; i < O; i++) g.poName[i] = "po" + std::to_string(i);
    while (r.p < buf.size()) {
        size_t eol = buf.find('\n', r.p);
        if (eol == std::string::npos) eol = buf.size();
        std::string line = buf.substr(r.p, eol - r.p);
        r.p = eol + 1;
        if (line.empty()) continue;
        if (line[0] == 'c') break;
        size_t sp = line.find(' ');
        if ((line[0] != 'i' && line[0] != 'o') || sp == std::string::npos) continue;
        unsigned idx = (unsigned)std::stoul(line.substr(1, sp - 1));
        std::string name = line.substr(sp + 1);
        if (line[0] == 'i' && idx < I) g.piName[idx] = name;
        if (line[0] == 'o' && idx < O) g.poName[idx] = name;
    }
    g.buildFanout();
    return g;
}
