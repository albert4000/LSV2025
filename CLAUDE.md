# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repository is

This is `LSV-PA`, a course fork of [berkeley-abc/abc](https://github.com/berkeley-abc/abc) (the ABC logic
synthesis/verification system) used for programming assignments in NTU's "Logic Synthesis and Verification"
course. It tracks upstream ABC closely (merged via `Merge berkeley-abc/master into course fork` commits) and
layers course-specific work on top in two places:

- `src/ext-lsv/` — student-written ABC commands (`lsv_*`), compiled as an ABC module.
- `lsv/` — assignment materials: PDFs, benchmarks, test scripts, admin scripts, and an `archive/` of past
  semesters' assignments.

Assignments are also sometimes implemented as modifications to core ABC source (e.g. adding a new flag to an
existing command like `&nf` in `src/aig/gia/giaNf.c`), not only as new `lsv_*` commands, so don't assume all
course work lives under `src/ext-lsv/` or `lsv/`.

## Build

```sh
make                          # builds the ./abc executable (default target)
make libabc.a                 # static library build
make ABC_USE_NAMESPACE=xxx    # compile as C++ in a namespace instead of C
make ABC_USE_NO_READLINE=1    # if libreadline is unavailable
make ABC_USE_NO_PTHREADS=1    # if pthreads are unavailable
make clean
```

CMake is also supported (used for the `test/` CTest suite):

```sh
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

If a build fails to start due to a cyclic-dependency check, run `find . -type f -exec touch "{}" \;` and rebuild.

Currently compiled into `abc` (see `src/ext-lsv/module.make`): `lsvUnateBdd.cpp`, `lsvUnateSat.cpp`, the
pipelining commands `lsvPipeMinCut.cpp`, `lsvPipeAcd.cpp`, `lsvPipeCmd.cpp` (header `lsvPipe.h`), the shared
window / code / truth-table synthesis core `lsvAcdCore.cpp` (header `lsvAcdCore.h`, used by `lsv_pipe_acd` and
`lsv_reencode`), and `lsvReencode.cpp` (`lsv_reencode`: window re-encoding analysis, kicks and ILS). Other files in
that directory (`kFeasibleCut.cpp`, `lsvCmdKLMOCut.cpp`, `lsvCmdHello.cpp`, `lsvCmd.cpp`) are commented out of
`SRC` there; add a new `.cpp` to that list to have it built. The module list that makes up `abc` overall is
`MODULES` in the top-level `Makefile`. ABC compiles C++ with `-fno-exceptions`: report errors by return value.

`lsv_pipe` (min-register S-stage pipelining) and `lsv_pipe_acd` (ACD re-encoding of the 2-stage cut) work on
the current GIA (`&r`, or derived via `&get -n`); `-p` replaces it with the pipelined sequential AIG. Their
regression test is `src/ext-lsv/tests/lsv_pipe_test.sh`. `src/develop/` holds the original standalone
prototypes (`aig_mincut/`, `acd_reencode/`) and the plans they implement (`plan/`). ABC in this tree cannot
read ASCII `.aag` files; convert to binary `.aig` first (e.g. with `src/develop/acd_reencode/aiger.py`).

## Running / testing

There is no single top-level test runner for assignments; each PA has its own ad hoc shell scripts that pipe
a batch of ABC commands into `./abc` and diff/inspect the output, e.g. `lsv/pa1/test/pa1.sh`,
`lsv/pa1/test/pa1_4.sh`, `lsv/pa1/test/compare_logs.sh`. Typical pattern:

```sh
./abc <<'EOF' > out.log 2>&1
read <benchmark.blif>
strash
lsv_print_nodes        # or whatever custom command is being exercised
quit
EOF
```

Interactive use: `./abc` drops into a REPL; `./abc -c "cmd1; cmd2; ..."` runs a batch of `;`-separated
commands non-interactively.

The CMake build additionally provides a small CTest regression suite under `test/` (e.g. `test/gia/`,
`test/write_cnf_dangling.sh`) unrelated to the course assignments — run via `ctest --test-dir build`.

EPFL benchmark circuits used for grading PA1/PA2 are expected as a symlink `./benchmarks` pointing at a
checkout of https://github.com/lsils/benchmarks (not included in this repo).

## Architecture: how a new command gets added

ABC is organized as `Abc_Frame_t`-scoped modules that register CLI commands at startup via
`Cmd_CommandAdd`. `src/ext-lsv/lsvCmd.cpp` is the minimal template for this pattern:

1. Declare the command handler: `int Lsv_CommandXxx(Abc_Frame_t* pAbc, int argc, char** argv)`.
2. In an `init(Abc_Frame_t* pAbc)` function, call
   `Cmd_CommandAdd(pAbc, "<group>", "<command_name>", Lsv_CommandXxx, <fChanges>)`.
3. Register `init`/`destroy` via a static `Abc_FrameInitializer_t` and a file-scope
   `Abc_FrameAddInitializer(&frame_initializer)` call (runs at static-init time, so just instantiating a
   global struct with this side effect is enough — see `PackageRegistrationManager` in `lsvCmd.cpp`).
4. Command handlers parse `argv` with `Extra_UtilGetopt`/`Extra_UtilGetoptReset` (getopt-style, ABC's own
   wrapper), fetch the current network with `Abc_FrameReadNtk(pAbc)`, print errors via `Abc_Print(-1, ...)`,
   and print usage via `Abc_Print(-2, ...)` on a `usage:` goto-label, matching the existing commands' style.
5. Add the new `.cpp` to `SRC` in `src/ext-lsv/module.make` so it's actually compiled in.

Network/graph data structures commonly touched by assignment code:
- `Abc_Ntk_t` / `Abc_Obj_t` (`src/base/abc/`) — the general logic network representation, iterated with
  macros like `Abc_NtkForEachNode`, `Abc_ObjForEachFanin`.
- BDD-based work goes through CUDD via `bdd/extrab/extraBdd.h` after running ABC's `collapse` command.
- SAT-based work goes through ABC's own SAT solvers under `src/sat/`.
- `Abc_Ntk_t` (AIG-oriented) vs. `Gia_Man_t` (`src/aig/gia/`) are two different network representations in
  this codebase; mapping/technology-mapping work (e.g. `&nf`) operates on `Gia_Man_t`, not `Abc_Ntk_t`.

## Submission workflow (only relevant if asked to prepare/submit an assignment)

Students fork this repo publicly, develop on a branch, and submit via pull request to their own
student-ID-named branch (not `master`). See `README.md` and `private-fork.md` for the full workflow,
including the private-repo trick to avoid exposing solutions before the deadline.
