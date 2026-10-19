// lsvAcdCore.cpp -- windows, image (equivalence classes), code assignment and truth-table synthesis
// shared by lsv_pipe_acd (lsvPipeAcd.cpp) and lsv_reencode (lsvReencode.cpp). See lsvAcdCore.h.
#include "lsvAcdCore.h"

#include <random>
#include <set>
#include <unordered_map>

#include <sys/wait.h>
#include <unistd.h>

#include "aig/gia/giaAig.h"
#include "misc/util/utilTruth.h"
#include "proof/cec/cec.h"
#include "proof/fra/fra.h"

extern "C" {
int Kit_TruthToGia(Gia_Man_t* pMan, unsigned* pTruth, int nVars, Vec_Int_t* vMemory, Vec_Int_t* vLeaves, int fHash);
word Abc_IsopCheck(word* pOn, word* pOnDc, word* pRes, int nVars, word CostLim, int* pCover);
}

namespace lsvpipe {

using std::string;
using std::vector;
u64 varWord(int i, size_t w) {
    static const u64 m[6] = { 0xAAAAAAAAAAAAAAAAULL, 0xCCCCCCCCCCCCCCCCULL, 0xF0F0F0F0F0F0F0F0ULL,
                              0xFF00FF00FF00FF00ULL, 0xFFFF0000FFFF0000ULL, 0xFFFFFFFF00000000ULL };
    return i < 6 ? m[i] : (((w >> (i - 6)) & 1) ? ~0ULL : 0ULL);
}

Keys makeKeys(size_t n, u64 seed) {
    std::mt19937_64 rng(seed);
    Keys k;
    for (size_t i = 0; i < n; i++) { k.k1.push_back(rng()); k.k2.push_back(rng()); }
    return k;
}
size_t countDistinct(vector<std::pair<u64, u64>> h) {
    std::sort(h.begin(), h.end());
    return std::unique(h.begin(), h.end()) - h.begin();
}

// ---------------------------------------------------------------------------------------------
// Windows


Window findWindow(const PipeAig& g, const vector<int>& roots, int K) {
    std::set<int> L(roots.begin(), roots.end());
    for (int steps = 0; steps < 200000; steps++) {
        int best = -1, bestCost = 1 << 30;
        for (int v : L) {
            if (!g.isAnd(v)) continue;
            int f0 = g.fanin0[v] >> 1, f1 = g.fanin1[v] >> 1, cost = -1;
            if (f0 && !L.count(f0)) cost++;
            if (f1 && f1 != f0 && !L.count(f1)) cost++;
            if (cost < bestCost || (cost == bestCost && v > best)) { bestCost = cost; best = v; }
        }
        if (best < 0 || (bestCost > 0 && (int)L.size() + bestCost > K)) break;
        L.erase(best);
        for (int lit : { g.fanin0[best], g.fanin1[best] }) if (lit >> 1) L.insert(lit >> 1);
    }
    Window w;
    w.leaves.assign(L.begin(), L.end());
    w.ok = (int)L.size() <= K;
    std::set<int> seen(L.begin(), L.end());
    vector<int> st(roots.begin(), roots.end());
    while (!st.empty()) {
        int v = st.back(); st.pop_back();
        if (seen.count(v)) continue;
        seen.insert(v);
        w.cone.push_back(v);
        for (int lit : { g.fanin0[v], g.fanin1[v] }) if ((lit >> 1) && !seen.count(lit >> 1)) st.push_back(lit >> 1);
    }
    std::sort(w.cone.begin(), w.cone.end());
    return w;
}

// ---------------------------------------------------------------------------------------------
// Group analysis: image points, pattern -> point map, functional-dependence subset


void analyzeGroup(const PipeAig& g, Group& G, const Window& w, bool computeFd) {
    int nL = (int)w.leaves.size(), m = (int)G.members.size();
    size_t P = 1ULL << nL, nW = (P + 63) / 64;
    G.leaves = w.leaves;
    std::unordered_map<int, int> idx;
    vector<vector<u64>> val;
    for (int i = 0; i < nL; i++) {
        idx[w.leaves[i]] = (int)val.size();
        vector<u64> t(nW);
        for (size_t k = 0; k < nW; k++) t[k] = varWord(i, k);
        val.push_back(t);
    }
    for (int v : w.cone) {
        if (idx.count(v)) continue;
        vector<u64> t(nW);
        int a = g.fanin0[v], b = g.fanin1[v];
        const vector<u64>* va = (a >> 1) ? &val[idx.at(a >> 1)] : nullptr;
        const vector<u64>* vb = (b >> 1) ? &val[idx.at(b >> 1)] : nullptr;
        for (size_t k = 0; k < nW; k++) {
            u64 x = va ? (*va)[k] : 0, y = vb ? (*vb)[k] : 0;
            t[k] = ((a & 1) ? ~x : x) & ((b & 1) ? ~y : y);
        }
        idx[v] = (int)val.size();
        val.push_back(t);
    }
    vector<const u64*> mv(m);   // member truth tables, looked up once (bit() runs 2^|L| * m times)
    for (int j = 0; j < m; j++) mv[j] = val[idx.at(G.members[j])].data();
    auto bit = [&](int j, size_t p) { return (char)((mv[j][p >> 6] >> (p & 63)) & 1); };
    Keys keys = makeKeys(m, 0x5eed + m);
    struct PairHash { size_t operator()(const std::pair<u64, u64>& h) const { return h.first ^ (h.second * 0x9e3779b97f4a7c15ULL); } };
    std::unordered_map<std::pair<u64, u64>, int, PairHash> seen;   // 128-bit hash of the value vector -> point
    vector<std::pair<u64, u64>> ph;
    vector<size_t> rep;
    G.inv.assign(P, 0);
    for (size_t p = 0; p < P; p++) {
        u64 h1 = 0, h2 = 0;
        for (int j = 0; j < m; j++) if (bit(j, p)) { h1 ^= keys.k1[j]; h2 ^= keys.k2[j]; }
        auto ins = seen.insert({ { h1, h2 }, (int)rep.size() });
        if (ins.second) { rep.push_back(p); ph.push_back({ h1, h2 }); }
        G.inv[p] = ins.first->second;
    }
    G.nImg = (int)rep.size();
    G.r = ceilLog2(G.nImg);
    G.U.assign(m, vector<char>(G.nImg));
    for (int j = 0; j < m; j++) for (int q = 0; q < G.nImg; q++) G.U[j][q] = bit(j, rep[q]);
    G.fdSet.clear();
    if (!computeFd) return;
    // greedy deletion, deepest members first
    vector<int> order(m);
    for (int j = 0; j < m; j++) order[j] = j;
    std::sort(order.begin(), order.end(), [&](int x, int y) { return G.members[x] > G.members[y]; });
    vector<char> dropped(m, 0);
    for (int j : order) {
        vector<std::pair<u64, u64>> h = ph;
        for (int k = 0; k < m; k++)
            if (dropped[k] || k == j)
                for (int q = 0; q < G.nImg; q++) if (G.U[k][q]) { h[q].first ^= keys.k1[k]; h[q].second ^= keys.k2[k]; }
        if ((int)countDistinct(h) == G.nImg) dropped[j] = 1;
    }
    for (int j = 0; j < m; j++) if (!dropped[j]) G.fdSet.push_back(j);
}

// ---------------------------------------------------------------------------------------------
// Synthesis of small multi-output functions from truth tables


// SOP -> AIG from an irredundant cover (cube bit 2v: literal !x_v, bit 2v+1: literal x_v).
// Returns -1 if the ISOP has more than maxCubes cubes.
static int sopToGia(Gia_Man_t* p, word* t, int nVars, const vector<int>& leaves, int maxCubes) {
    static vector<int> cover(1 << 16);
    vector<word> res(nWordsOf(nVars));
    word limit = (word)(maxCubes + 1) << 32;
    word cost = Abc_IsopCheck(t, t, res.data(), nVars, limit, cover.data());
    int nCubes = (int)(cost >> 32);
    if (cost >= limit) { std::fill(cover.begin(), cover.end(), 0); return -1; }
    auto tree = [&](vector<int> v, bool isAnd) {
        if (v.empty()) return isAnd ? 1 : 0;
        while (v.size() > 1) {
            vector<int> w;
            for (size_t i = 0; i + 1 < v.size(); i += 2) w.push_back(isAnd ? Gia_ManHashAnd(p, v[i], v[i + 1]) : Gia_ManHashOr(p, v[i], v[i + 1]));
            if (v.size() & 1) w.push_back(v.back());
            v.swap(w);
        }
        return v[0];
    };
    vector<int> cubes;
    for (int c = 0; c < nCubes; c++) {
        vector<int> lits;
        for (int v = 0; v < nVars; v++) {
            if (cover[c] & (1 << (2 * v))) lits.push_back(Abc_LitNot(leaves[v]));
            if (cover[c] & (1 << (2 * v + 1))) lits.push_back(leaves[v]);
        }
        cubes.push_back(tree(lits, true));
        cover[c] = 0;
    }
    return tree(cubes, false);
}

vector<vector<word>> leafPatterns(int nVars) {
    int nw = nWordsOf(nVars);
    vector<vector<word>> in(nVars, vector<word>(nw));
    for (int i = 0; i < nVars; i++) for (int k = 0; k < nw; k++) in[i][k] = varWord(i, k);
    return in;
}

vector<vector<word>> simulate(Gia_Man_t* p, const vector<vector<word>>& in, int nw) {
    vector<vector<word>> val(Gia_ManObjNum(p), vector<word>(nw, 0));
    Gia_Obj_t* pObj;
    int i;
    Gia_ManForEachCi(p, pObj, i) if (i < (int)in.size()) val[Gia_ObjId(p, pObj)] = in[i];
    Gia_ManForEachAnd(p, pObj, i) {
        const vector<word>& x = val[Gia_ObjFaninId0(pObj, i)];
        const vector<word>& y = val[Gia_ObjFaninId1(pObj, i)];
        for (int k = 0; k < nw; k++) val[i][k] = (Gia_ObjFaninC0(pObj) ? ~x[k] : x[k]) & (Gia_ObjFaninC1(pObj) ? ~y[k] : y[k]);
    }
    vector<vector<word>> out;
    Gia_ManForEachCo(p, pObj, i) {
        vector<word> x = val[Gia_ObjFaninId0p(p, pObj)];
        if (Gia_ObjFaninC0(pObj)) for (auto& t : x) t = ~t;
        out.push_back(x);
    }
    return out;
}

bool sameTruths(const vector<vector<word>>& a, const vector<vector<word>>& b, int nVars) {
    if (a.size() != b.size()) return false;
    word mask = nVars >= 6 ? ~(word)0 : Abc_Tt6Mask(1 << nVars);
    for (size_t j = 0; j < a.size(); j++)
        for (size_t k = 0; k < a[j].size(); k++) if ((a[j][k] ^ b[j][k]) & mask) return false;
    return true;
}

bool checkTruths(Gia_Man_t* p, const vector<vector<word>>& tts, int nVars) {
    return sameTruths(simulate(p, leafPatterns(nVars), nWordsOf(nVars)), tts, nVars);
}

// Returns nullptr if the result would exceed the AND budget (or on an internal error).
Gia_Man_t* giaFromTruths(const vector<vector<word>>& tts, int nVars, int budget) {
    Gia_Man_t* p = Gia_ManStart(1000);
    Gia_ManHashAlloc(p);
    Vec_Int_t* vLeaves = Vec_IntAlloc(nVars);
    vector<int> leaves;
    for (int i = 0; i < std::max(nVars, 1); i++) { int l = Gia_ManAppendCi(p); Vec_IntPush(vLeaves, l); leaves.push_back(l); }
    for (auto t : tts) {
        int nw = nWordsOf(nVars), lit;
        if (nVars < 6) t[0] = Abc_Tt6Stretch(t[0], nVars);
        if (Abc_TtIsConst0(t.data(), nw)) lit = 0;
        else if (Abc_TtIsConst1(t.data(), nw)) lit = 1;
        else if (nVars < 16) lit = Kit_TruthToGia(p, (unsigned*)t.data(), nVars, NULL, vLeaves, 1);
        else lit = sopToGia(p, t.data(), nVars, leaves, budget);   // Kit's factoring is limited to < 16 vars
        if (lit < 0 || Gia_ManAndNum(p) > budget) { Vec_IntFree(vLeaves); Gia_ManStop(p); return nullptr; }
        Gia_ManAppendCo(p, lit);
    }
    Vec_IntFree(vLeaves);
    Gia_ManHashStop(p);
    Gia_Man_t* q = Gia_ManCleanup(p);
    Gia_ManStop(p);
    if (!checkTruths(q, tts, nVars)) {
        printf("lsv acd: internal error: synthesized function does not match its truth table (%d vars)\n", nVars);
        Gia_ManStop(q);
        return nullptr;
    }
    return q;
}

// Candidates: for each don't-care fill, the raw ISOP-factored AIG, dc2 and syn2. Returns the
// min-AND candidate and (if different) the min-depth one; all others are freed.
// budget: AND count beyond which a raw (unoptimized) result is not worth optimizing.
vector<Synth> synthesize(const vector<vector<word>>& on, const vector<vector<word>>& dc, int nVars, int budget) {
    int nw = nWordsOf(nVars);
    vector<vector<vector<word>>> fills{ on };
    bool hasDc = false;
    for (auto& d : dc) hasDc |= !Abc_TtIsConst0((word*)d.data(), nw);
    if (hasDc) {
        vector<vector<word>> one = on, isop = on;
        for (size_t j = 0; j < on.size(); j++) {
            vector<word> ondc(nw), res(nw);
            for (int k = 0; k < nw; k++) ondc[k] = on[j][k] | dc[j][k];
            one[j] = ondc;
            vector<word> o = on[j], od = ondc;
            if (nVars < 6) { o[0] = Abc_Tt6Stretch(o[0], nVars); od[0] = Abc_Tt6Stretch(od[0], nVars); }
            Abc_IsopCheck(o.data(), od.data(), res.data(), nVars, ~(word)0 >> 1, NULL);
            bool okRes = true;
            for (int k = 0; k < nw; k++) if ((o[k] & ~res[k]) || (res[k] & ~od[k])) okRes = false;
            if (okRes) isop[j] = res;
        }
        fills.push_back(one);
        fills.push_back(isop);
    }
    vector<Synth> all;
    for (size_t i = 0; i < fills.size(); i++) {
        const vector<vector<word>>& f = fills[i];
        // a fill equal to an earlier one gives the same candidates, which never win a tie against the earlier ones
        bool dup = false;
        for (size_t j = 0; j < i && !dup; j++) dup = sameTruths(f, fills[j], nVars);
        if (dup) continue;
        Gia_Man_t* base = giaFromTruths(f, nVars, budget);
        if (!base) continue;
        Gia_Man_t* vars[3] = { base, Gia_ManCompress2(base, 1, 0), Gia_ManAigSyn2(base, 0, 1, 0, 20, 0, 0, 0) };
        for (Gia_Man_t* q : vars) {
            Synth s; s.p = q; s.ands = Gia_ManAndNum(q); s.depth = Gia_ManLevelNum(q);
            all.push_back(s);
        }
    }
    if (all.empty()) return {};
    size_t ia = 0, id = 0;
    for (size_t i = 1; i < all.size(); i++) {
        if (std::make_pair(all[i].ands, all[i].depth) < std::make_pair(all[ia].ands, all[ia].depth)) ia = i;
        if (std::make_pair(all[i].depth, all[i].ands) < std::make_pair(all[id].depth, all[id].ands)) id = i;
    }
    vector<Synth> res{ all[ia] };
    if (all[id].ands != all[ia].ands || all[id].depth != all[ia].depth) res.push_back(all[id]);
    for (size_t i = 0; i < all.size(); i++) if (i != ia && (res.size() < 2 || i != id)) Gia_ManStop(all[i].p);
    return res;
}

void freeSynth(vector<Synth>& v) { for (auto& s : v) if (s.p) Gia_ManStop(s.p); v.clear(); }

// Truth tables (bit p = minterm p) of decoder outputs: code[q] -> U[j][q]; unused codes are DC.
void decoderTruths(const Group& G, const vector<int>& code, int nVars, const vector<int>& outs,
                          vector<vector<word>>& on, vector<vector<word>>& dc) {
    int nw = nWordsOf(nVars);
    vector<word> used(nw, 0);
    for (int c : code) used[c >> 6] |= 1ULL << (c & 63);
    vector<word> d(nw);
    for (int k = 0; k < nw; k++) d[k] = ~used[k];
    if (nVars < 6) d[0] &= Abc_Tt6Mask(1 << nVars);
    on.clear(); dc.clear();
    for (int j : outs) {
        vector<word> t(nw, 0);
        for (int q = 0; q < G.nImg; q++) if (G.U[j][q]) t[code[q] >> 6] |= 1ULL << (code[q] & 63);
        on.push_back(t); dc.push_back(d);
    }
}

// Greedy code assignment: reuse G columns as code bits while the total stays r.
void assignCode(const Group& G, vector<int>& reused, int& nnew, vector<int>& idx, vector<int>& code) {
    int n = G.nImg, m = (int)G.members.size(), r = G.r;
    vector<int> cls(n, 0);
    int nCls = 1;
    reused.clear();
    while (true) {
        int bestJ = -1, bestN = -1, bestMax = 0;
        for (int j = 0; j < m; j++) {
            if (std::find(reused.begin(), reused.end(), j) != reused.end()) continue;
            vector<int> cnt(2 * nCls, 0);
            for (int q = 0; q < n; q++) cnt[2 * cls[q] + G.U[j][q]]++;
            int nc = 0, mx = 0;
            for (int c : cnt) { nc += c > 0; mx = std::max(mx, c); }
            if ((int)reused.size() + 1 + ceilLog2(mx) > r) continue;
            if (nc > bestN || (nc == bestN && mx < bestMax)) { bestJ = j; bestN = nc; bestMax = mx; }
        }
        if (bestJ < 0) break;
        reused.push_back(bestJ);
        vector<int> remap(2 * nCls, -1);
        int k = 0;
        for (int q = 0; q < n; q++) {
            int c = 2 * cls[q] + G.U[bestJ][q];
            if (remap[c] < 0) remap[c] = k++;
            cls[q] = remap[c];
        }
        nCls = k;
    }
    vector<int> cnt(nCls, 0);
    idx.assign(n, 0);
    int mx = 0;
    for (int q = 0; q < n; q++) { idx[q] = cnt[cls[q]]++; mx = std::max(mx, cnt[cls[q]]); }
    nnew = ceilLog2(mx);
    code.assign(n, 0);
    for (int q = 0; q < n; q++) {
        int c = 0;
        for (size_t k = 0; k < reused.size(); k++) c |= G.U[reused[k]][q] << k;
        code[q] = c | (idx[q] << reused.size());
    }
}

// Combinational equivalence of two GIAs (XOR miter, fraiging + SAT, as in dcec): 1 / 0 / -1.
// The GIA-based CEC engine is much slower on circuits like div, so the Fraig engine is used.
int cecGia(Gia_Man_t* p0, Gia_Man_t* p1) {
    Gia_Man_t* pMiter = Gia_ManMiter(p0, p1, 0, 0, 0, 0, 0);
    if (!pMiter) return -1;
    Aig_Man_t* pAig = Gia_ManToAig(pMiter, 0);
    Gia_ManStop(pMiter);
    int r = Fra_FraigCec(&pAig, 100000, 0);
    Aig_ManStop(pAig);
    return r;
}

// Same with SAT sweeping only (the first engine of &cec), never falling back to the Fraig engine,
// which hits an internal assertion (Dar_ManRwsat -> Aig_ManDupDfs) on some sqrt miters.
// 1: every miter output reduced to constant 0; 0: counter-example; -1: undecided.
int cecGiaCec(Gia_Man_t* p0, Gia_Man_t* p1) {
    Gia_Man_t* pMiter = Gia_ManMiter(p0, p1, 0, 1, 0, 0, 0);
    if (!pMiter) return -1;
    Cec_ParFra_t pars;
    Cec_ManFraSetDefaultParams(&pars);
    pars.nItersMax = 1000;
    pars.fCheckMiter = 1;
    pars.fDualOut = 1;
    Gia_Man_t* pRed = Cec_ManSatSweeping(pMiter, &pars, 1);
    int r = -1;
    if (!pRed) r = pMiter->pCexComb ? 0 : -1;
    else {
        r = 1;
        Gia_Obj_t* pObj;
        int i;
        Gia_ManForEachCo(pRed, pObj, i) if (Gia_ObjFaninLit0p(pRed, pObj) != 0) r = -1;
        Gia_ManStop(pRed);
    }
    if (pMiter->pCexComb) { ABC_FREE(pMiter->pCexComb); r = 0; }
    Gia_ManStop(pMiter);
    return r;
}

// cecGia (fast Fraig engine) in a forked child, so that an internal assertion of ABC cannot take down
// the caller; if the child crashes or is undecided, SAT sweeping (cecGiaCec) decides. 1 / 0 / -1.
int cecGiaSafe(Gia_Man_t* p0, Gia_Man_t* p1) {
    fflush(stdout);
    fflush(stderr);
    pid_t pid = fork();
    if (pid == 0) {
        int r = cecGia(p0, p1);
        _exit(r == 1 ? 0 : r == 0 ? 1 : 2);
    }
    if (pid > 0) {
        int status = 0;
        if (waitpid(pid, &status, 0) == pid && WIFEXITED(status)) {
            if (WEXITSTATUS(status) == 0) return 1;
            if (WEXITSTATUS(status) == 1) return 0;
        }
    }
    return cecGiaCec(p0, p1);
}

}  // namespace lsvpipe
