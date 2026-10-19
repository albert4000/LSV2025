// lsvAcdCore.h -- building blocks shared by lsv_pipe_acd and lsv_reencode:
//   windows on a PipeAig, the image of a window's signals over all leaf patterns (= the output
//   equivalence classes), r-bit code assignment, and synthesis of small multi-output functions
//   from truth tables inside ABC.
#ifndef LSV_ACD_CORE_H
#define LSV_ACD_CORE_H

#include <utility>
#include <vector>

#include "lsvPipe.h"

namespace lsvpipe {

typedef unsigned long long u64;

inline int ceilLog2(u64 n) { int r = 0; while ((1ULL << r) < n) r++; return r; }
inline int nWordsOf(int nVars) { return nVars <= 6 ? 1 : 1 << (nVars - 6); }
// Word w of the truth table of variable i (bit p = minterm p).
u64 varWord(int i, size_t w);

// Random 128-bit keys for hashing bit vectors (XOR of the keys of the set bits).
struct Keys { std::vector<u64> k1, k2; };
Keys makeKeys(size_t n, u64 seed);
size_t countDistinct(std::vector<std::pair<u64, u64>> h);

// ---------------------------------------------------------------------------------------------
// Windows

struct Window { std::vector<int> leaves, cone; bool ok = false; };

// Reconvergence-driven leaves (<= K) above roots; cone = the nodes between the leaves and the roots
// (sorted, i.e. topological). ok = false if the roots alone already need more than K leaves.
Window findWindow(const PipeAig& g, const std::vector<int>& roots, int K);

// ---------------------------------------------------------------------------------------------
// Image of a group of signals over all 2^|leaves| leaf patterns

struct Group {
    std::vector<int> members, leaves, fdSet;   // fdSet: member indices
    int nImg = 0, r = 0, coneAnds = 0;
    int maxArr = 0, maxDown = 0;
    std::vector<std::vector<char>> U;          // U[j][point]: value of member j at image point
    std::vector<int> inv;                      // leaf pattern -> image point (= equivalence class)
};

// w.cone must contain every AND between w.leaves and G.members, in topological order.
// computeFd: also compute a minimal functionally-dependent subset (greedy deletion; O(|G|^2 N)).
void analyzeGroup(const PipeAig& g, Group& G, const Window& w, bool computeFd = true);

// Greedy code assignment: reuse G columns as code bits while the total stays r.
// code[q] = r-bit codeword of image point q.
void assignCode(const Group& G, std::vector<int>& reused, int& nnew, std::vector<int>& idx, std::vector<int>& code);

// ---------------------------------------------------------------------------------------------
// Synthesis of small multi-output functions from truth tables

struct Synth { Gia_Man_t* p = nullptr; int ands = 0, depth = 0; };

// Truth tables of the nVars variables (CI i driven by variable i).
std::vector<std::vector<word>> leafPatterns(int nVars);
// CO functions when CI i is driven by in[i] (nw words each; CIs beyond in are 0).
std::vector<std::vector<word>> simulate(Gia_Man_t* p, const std::vector<std::vector<word>>& in, int nw);
// Equal on the 2^nVars minterms (bits above them in a single word are ignored).
bool sameTruths(const std::vector<std::vector<word>>& a, const std::vector<std::vector<word>>& b, int nVars);
// Simulate all outputs over the 2^nVars minterms and compare with the requested truth tables.
bool checkTruths(Gia_Man_t* p, const std::vector<std::vector<word>>& tts, int nVars);
// Returns nullptr if the result would exceed the AND budget (or on an internal error).
Gia_Man_t* giaFromTruths(const std::vector<std::vector<word>>& tts, int nVars, int budget);
// Candidates: for each don't-care fill (DC=0 | DC=1 | ISOP), the raw factored AIG, dc2 and syn2.
// Returns the min-AND candidate and (if different) the min-depth one; empty if all exceed budget.
std::vector<Synth> synthesize(const std::vector<std::vector<word>>& on, const std::vector<std::vector<word>>& dc,
                              int nVars, int budget);
void freeSynth(std::vector<Synth>& v);
// Decoder truth tables over nVars = r code bits: output j at code[q] is U[outs[j]][q]; unused codes are DC.
void decoderTruths(const Group& G, const std::vector<int>& code, int nVars, const std::vector<int>& outs,
                   std::vector<std::vector<word>>& on, std::vector<std::vector<word>>& dc);

// Combinational equivalence of two GIAs (XOR miter, fraiging + SAT): 1 / 0 / -1 (undecided).
int cecGia(Gia_Man_t* p0, Gia_Man_t* p1);
// Same with SAT sweeping only (never falls back to the Fraig engine, which can hit an internal assertion).
int cecGiaCec(Gia_Man_t* p0, Gia_Man_t* p1);
// cecGia in a forked child (a crash there is survivable); falls back to cecGiaCec. 1 / 0 / -1.
int cecGiaSafe(Gia_Man_t* p0, Gia_Man_t* p1);

}  // namespace lsvpipe

#endif
