// lsvPipeCmd.cpp -- ABC commands lsv_pipe and lsv_pipe_acd (see lsvPipe.h).
#include "lsvPipe.h"

#include <climits>
#include <cstdlib>
#include <string>

#include "base/abc/abc.h"
#include "base/main/main.h"
#include "base/main/mainInt.h"
#include "misc/extra/extra.h"

namespace lsvpipe {

// The current GIA; derived from the current network (&get -n) if there is none or if forced.
Gia_Man_t* currentGia(Abc_Frame_t* pAbc, bool fromNtk) {
    if ((fromNtk || !Abc_FrameReadGia(pAbc)) && Abc_FrameReadNtk(pAbc)) {
        char cmd[] = "&get -n";
        if (Cmd_CommandExecute(pAbc, cmd)) return nullptr;
    }
    return Abc_FrameReadGia(pAbc);
}

bool parseInt(int argc, char** argv, int& value, int minValue, const char* flag) {
    if (globalUtilOptind >= argc) { Abc_Print(-1, "Command line switch \"-%s\" should be followed by an integer.\n", flag); return false; }
    char* end;
    long v = strtol(argv[globalUtilOptind], &end, 10);
    if (end == argv[globalUtilOptind] || *end || v < INT_MIN || v > INT_MAX) {
        Abc_Print(-1, "Command line switch \"-%s\" should be followed by an integer, not \"%s\".\n", flag, argv[globalUtilOptind]);
        return false;
    }
    globalUtilOptind++;
    value = (int)v;
    if (value < minValue) { Abc_Print(-1, "The value of \"-%s\" should be >= %d.\n", flag, minValue); return false; }
    return true;
}

}  // namespace lsvpipe

namespace {

using namespace lsvpipe;

int Lsv_CommandPipe(Abc_Frame_t* pAbc, int argc, char** argv);
int Lsv_CommandPipeAcd(Abc_Frame_t* pAbc, int argc, char** argv);

void init(Abc_Frame_t* pAbc) {
    Cmd_CommandAdd(pAbc, "LSV", "lsv_pipe", Lsv_CommandPipe, 0);
    Cmd_CommandAdd(pAbc, "LSV", "lsv_pipe_acd", Lsv_CommandPipeAcd, 0);
}
void destroy(Abc_Frame_t* pAbc) {}
Abc_FrameInitializer_t frame_initializer = { init, destroy };
struct PackageRegistrationManager {
    PackageRegistrationManager() { Abc_FrameAddInitializer(&frame_initializer); }
} lsvPipeRegistrationManager;

int Lsv_CommandPipe(Abc_Frame_t* pAbc, int argc, char** argv) {
    MinCutOptions opt;
    int c, W = 1, fReplace = 0, fFromNtk = 0;
    std::string jsonDir;
    Extra_UtilGetoptReset();
    while ((c = Extra_UtilGetopt(argc, argv, "STWbsjflkpNvh")) != EOF) {
        switch (c) {
        case 'S': if (!parseInt(argc, argv, opt.S, 1, "S")) goto usage; break;
        case 'T': if (!parseInt(argc, argv, opt.T, 0, "T")) goto usage; break;
        case 'W': if (!parseInt(argc, argv, W, 1, "W")) goto usage; break;
        case 'b': if (!parseInt(argc, argv, opt.brute, 0, "b")) goto usage; break;
        case 's': if (!parseInt(argc, argv, opt.simCycles, 0, "s")) goto usage; break;
        case 'j':
            if (globalUtilOptind >= argc) { Abc_Print(-1, "Command line switch \"-j\" should be followed by a directory.\n"); goto usage; }
            jsonDir = argv[globalUtilOptind++];
            break;
        case 'f': opt.method = 1; break;
        case 'l': opt.method = 2; break;
        case 'k': opt.sinkSide = !opt.sinkSide; break;
        case 'p': fReplace ^= 1; break;
        case 'N': fFromNtk ^= 1; break;
        case 'v': opt.verbose = !opt.verbose; break;
        default: goto usage;
        }
    }
    if (opt.method == 1 && (opt.S != 2 || opt.T)) { Abc_Print(-1, "-f (max-flow only) needs -S 2 and no -T.\n"); return 1; }
    opt.W = W;
    {
        Gia_Man_t* pGia = currentGia(pAbc, fFromNtk);
        if (!pGia) { Abc_Print(-1, "There is no AIG (read a circuit, or use &r).\n"); return 1; }
        PipeAig g;
        std::string msg;
        if (!pipeAigFromGia(pGia, g, msg)) { Abc_Print(-1, "%s\n", msg.c_str()); return 1; }
        Problem P(g, opt.S, opt.W, opt.T);
        buildProblem(P);
        int nLive = 0, depth = 0;
        for (int v = g.nPi + 1; v < g.numNodes(); v++) if (P.live[v]) { nLive++; depth = std::max(depth, P.L[v]); }
        printf("lsv_pipe: %d PI, %d PO, %d AND (%d live), depth %d; S=%d T=%d W=%lld cut side=%s\n", g.nPi,
               (int)g.poLit.size(), g.numAnds(), nLive, depth, opt.S, opt.T, opt.W, opt.sinkSide ? "sink" : "source");
        Solution sol;
        MinCutResult r = solveMinCut(P, opt, sol);
        if (!r.error.empty()) { Abc_Print(-1, "Internal error: %s\n", r.error.c_str()); return 1; }
        if (r.infeasible) { printf("INFEASIBLE: no assignment with S=%d and stage depth <= %d (circuit depth %d)\n", opt.S, opt.T, depth); return 0; }
        if (r.flowValue >= 0) printf("  [A] max-flow min-cut = %lld (%.3fs)\n", r.flowValue, r.tFlow);
        if (opt.method != 1 && opt.S > 1)
            printf("  [B] LP/min-cost-flow: %lld registers, %d mcf phases, %d round(s), %d lazy constraints (%.3fs)\n",
                   sol.regs, r.lpPhases, r.lpRounds, r.lazy, r.tLp);
        printReport(P, sol);
        printf("  checks          : %s\n", r.ok ? "all passed" : "FAILED (see above)");
        if (!jsonDir.empty()) {
            if (writeJson(P, sol, jsonDir)) printf("  wrote %s/stage_assign.json and %s/cuts.json\n", jsonDir.c_str(), jsonDir.c_str());
            else Abc_Print(-1, "Cannot write JSON files to \"%s\".\n", jsonDir.c_str());
        }
        if (fReplace && !r.ok) printf("  current GIA not replaced: the checks failed\n");
        if (fReplace && r.ok) {
            Gia_Man_t* pNew = buildPipelineGia(P, sol, pGia);
            printf("  current GIA replaced by the pipeline: %d flops, %d ANDs (equivalent to the original with PIs delayed by %d)\n",
                   Gia_ManRegNum(pNew), Gia_ManAndNum(pNew), opt.S - 1);
            Abc_FrameUpdateGia(pAbc, pNew);
        }
        return r.ok ? 0 : 1;
    }
usage:
    Abc_Print(-2, "usage: lsv_pipe [-STWbsj num] [-flkpNvh]\n");
    Abc_Print(-2, "\t         balanced S-stage pipelining of the current (combinational) AIG with the\n");
    Abc_Print(-2, "\t         minimum number of registers (exact: max-flow for S=2, LP/min-cost flow otherwise)\n");
    Abc_Print(-2, "\t-S num : number of pipeline stages [default = 2]\n");
    Abc_Print(-2, "\t-T num : max AND depth per stage, 0 = unconstrained [default = 0]\n");
    Abc_Print(-2, "\t-W num : register weight (area per register) [default = 1]\n");
    Abc_Print(-2, "\t-b num : brute-force optimality check if <= num live ANDs, 0 = off [default = 0]\n");
    Abc_Print(-2, "\t-s num : cycles of latency simulation check, 0 = off [default = 32]\n");
    Abc_Print(-2, "\t-j dir : write stage_assign.json and cuts.json to dir\n");
    Abc_Print(-2, "\t-f     : max-flow only (S=2, no -T) [default = LP, + max-flow cross-check when possible]\n");
    Abc_Print(-2, "\t-l     : LP only, no max-flow cross-check\n");
    Abc_Print(-2, "\t-k     : toggle cut nearest the POs instead of the PIs [default = PIs]\n");
    Abc_Print(-2, "\t-p     : replace the current GIA by the pipelined sequential AIG [default = no]\n");
    Abc_Print(-2, "\t-N     : derive the AIG from the current network (&get -n) even if a GIA exists\n");
    Abc_Print(-2, "\t-v     : toggle verbose output (all check results)\n");
    Abc_Print(-2, "\t-h     : print the command usage\n");
    return 1;
}

int Lsv_CommandPipeAcd(Abc_Frame_t* pAbc, int argc, char** argv) {
    AcdOptions acd;
    MinCutOptions mc;
    int c, W = 10, fReplace = 0, fFromNtk = 0;
    Extra_UtilGetoptReset();
    while ((c = Extra_UtilGetopt(argc, argv, "TKWkgxaGpNvh")) != EOF) {
        switch (c) {
        case 'T': if (!parseInt(argc, argv, mc.T, 0, "T")) goto usage; break;
        case 'K': if (!parseInt(argc, argv, acd.K, 1, "K")) goto usage; if (acd.K > 16) { Abc_Print(-1, "-K must be <= 16.\n"); goto usage; } break;
        case 'W': if (!parseInt(argc, argv, W, 1, "W")) goto usage; break;
        case 'k': mc.sinkSide = !mc.sinkSide; break;
        case 'g': acd.slackGrouping = !acd.slackGrouping; break;
        case 'x': acd.ignoreDelay = !acd.ignoreDelay; break;
        case 'a': acd.analysisOnly = !acd.analysisOnly; break;
        case 'G': acd.global = !acd.global; break;
        case 'p': fReplace ^= 1; break;
        case 'N': fFromNtk ^= 1; break;
        case 'v': acd.verbose = !acd.verbose; break;
        default: goto usage;
        }
    }
    acd.regWeight = W;
    {
        Gia_Man_t* pGia = currentGia(pAbc, fFromNtk);
        if (!pGia) { Abc_Print(-1, "There is no AIG (read a circuit, or use &r).\n"); return 1; }
        PipeAig g;
        std::string msg;
        if (!pipeAigFromGia(pGia, g, msg)) { Abc_Print(-1, "%s\n", msg.c_str()); return 1; }
        Problem P(g, 2, 1, mc.T);
        buildProblem(P);
        int depth = 0;
        for (int v = g.nPi + 1; v < g.numNodes(); v++) if (P.live[v]) depth = std::max(depth, P.L[v]);
        printf("lsv_pipe_acd: %d PI, %d PO, %d AND, depth %d; T=%d K=%d W=%d cut side=%s\n", g.nPi, (int)g.poLit.size(),
               g.numAnds(), depth, mc.T, acd.K, W, mc.sinkSide ? "sink" : "source");
        Solution sol;
        mc.simCycles = 16;
        MinCutResult r = solveMinCut(P, mc, sol);
        if (!r.error.empty()) { Abc_Print(-1, "Internal error: %s\n", r.error.c_str()); return 1; }
        if (r.infeasible) { printf("INFEASIBLE: no 2-stage assignment with stage depth <= %d\n", mc.T); return 0; }
        if (!r.ok) { Abc_Print(-1, "Min-cut checks failed.\n"); return 1; }
        std::vector<int> d = stageDepth(P, sol.s);
        int d0 = 0, d1 = 0;
        for (int v = g.nPi + 1; v < g.numNodes(); v++) if (P.live[v]) (sol.s[v] ? d1 : d0) = std::max(sol.s[v] ? d1 : d0, d[v]);
        printf("  min-cut: %lld registers, stage depths %d/%d\n", sol.regs, d0, d1);
        Gia_Man_t* pNew = runAcd(P, sol, acd, pGia, msg);
        if (!pNew) {
            if (!msg.empty()) { Abc_Print(-1, "%s\n", msg.c_str()); return 1; }
            return 0;
        }
        if (fReplace) {
            printf("  current GIA replaced by the re-encoded pipeline: %d flops, %d ANDs\n", Gia_ManRegNum(pNew), Gia_ManAndNum(pNew));
            Abc_FrameUpdateGia(pAbc, pNew);
        } else Gia_ManStop(pNew);
        return 0;
    }
usage:
    Abc_Print(-2, "usage: lsv_pipe_acd [-TKW num] [-kgxaGpNvh]\n");
    Abc_Print(-2, "\t         2-stage min-register pipelining followed by ACD re-encoding of the cut\n");
    Abc_Print(-2, "\t         (functional dependence / r-bit encoders and decoders synthesized in ABC)\n");
    Abc_Print(-2, "\t-T num : max AND depth per stage, 0 = unconstrained [default = 0]\n");
    Abc_Print(-2, "\t-K num : window leaf limit, <= 16 [default = 16]\n");
    Abc_Print(-2, "\t-W num : register area in AND (NAND2) units [default = 10]\n");
    Abc_Print(-2, "\t-k     : toggle cut nearest the POs instead of the PIs [default = PIs]\n");
    Abc_Print(-2, "\t-g     : toggle grouping within slack classes [default = leaf overlap only]\n");
    Abc_Print(-2, "\t-x     : toggle accepting candidates on area alone (stage depth may exceed T)\n");
    Abc_Print(-2, "\t-a     : toggle analysis only (phase 1, no circuit change)\n");
    Abc_Print(-2, "\t-G     : toggle global sampling of the cut (2^20 patterns or exhaustive)\n");
    Abc_Print(-2, "\t-p     : replace the current GIA by the re-encoded sequential AIG [default = no]\n");
    Abc_Print(-2, "\t-N     : derive the AIG from the current network (&get -n) even if a GIA exists\n");
    Abc_Print(-2, "\t-v     : toggle verbose output (per-group details)\n");
    Abc_Print(-2, "\t-h     : print the command usage\n");
    return 1;
}

}  // namespace
