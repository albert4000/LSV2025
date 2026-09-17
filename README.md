# LSV-PA, Fall 2025
This is a repository to host programming assignments of the course **Logic Synthesis and Verification** at National Taiwan University.
It is forked from the repository [ABC](https://github.com/berkeley-abc/abc) of UC Berkeley.

## Submission Workflow
We will follow a common [git feature branch workflow](https://www.atlassian.com/git/tutorials/comparing-workflows/feature-branch-workflow) for submission.
All enrolled students will have their own branches named by their students' ID numbers.
You are asked to **fork** this repository and develop your program on the forked repository.
After finishing your programming assignments, the assignments **must** be submitted via a pull request to a student's **own branch** (please do not modify the master branch).
If you cannot find your own branch, please contact the TA.
If you don't know how to create a pull request, please read through this [document](https://guides.github.com/activities/forking/).

### To avoid plagiarism ...
Please note that your fork of this public repository will also be public,
which means that if you push your code to the fork, it is visible to everyone.
In case you want to prevent other students from copying your solution,
an easy way is to push and create a pull request **at the last moment before the deadline**.

Another complicated way is to create a private repository to develop your solutions,
pull your code to the public fork after an assignment is finished,
and create a pull request via the public fork.
The benefit of this method is that you can push your code during the development and keep it private.
The drawback is again you need to create a pull request close to the deadline, as PRs are visible to everyone.
The detailed steps are documented [here](./private-fork.md).

## Participants
Please register your student IDs and GitHub accounts in this [table](./lsv/admin/participants-id.csv).

## Assignments
### [PA1](./lsv/pa1/pa1.pdf): Multi-output Cut Enumeration
Submission deadline:
- Exercises 1-3: 2023/09/21 23:59
- Exercises 4: 2023/10/05 23:59

## Contact
TA: Kuo-Wei Ho (f11943109@ntu.edu.tw)

For questions, you are encouraged to open an [issue](https://github.com/NTU-ALComLab/LSV-PA/issues).
As other students might have the same questions, discussing in an issue will benefit everyone.
Note that you can set labels, e.g., `PA0`, `PA1`, etc, to classify your questions.

---

# Upstream ABC Documentation

The sections below are from the original [berkeley-abc/abc](https://github.com/berkeley-abc/abc) README and describe how to build and use ABC itself.

## Compiling

To compile the ABC executable, run `make`. To compile ABC as a static library,
run `make libabc.a`.

When ABC is used as a static library, `Abc_Start()` and `Abc_Stop()` are
provided for starting and quitting a single ABC session in the calling
application. A simple demo program (file src/demo.c) shows how to
create a stand-alone program performing DAG-aware AIG rewriting, by calling 
APIs of ABC compiled as a static library.

Independent ABC sessions can run concurrently in different threads on Linux
and macOS. Existing single-threaded applications using `Abc_Start()` and
`Abc_Stop()` continue to work unchanged. Concurrent applications should create
one `Abc_Frame_t` per session with `Abc_FrameCreate()`, execute commands using
that frame, and release it with `Abc_FrameDestroy()`; the same frame must not be
used concurrently by multiple threads. Scoped calls to
`Abc_FrameEnter()`/`Abc_FrameLeave()` establish the current frame when invoking
APIs outside command execution that depend on it. See `src/runabc.c` for
single-threaded and concurrent embedding examples. Some obsolete commands still
use shared process state and are not safe to run in concurrent sessions.

Build and run the single-threaded demo from the repository root with:

```sh
gcc -Wall -O2 -c src/demo.c -o demo.o
g++ -o demo demo.o libabc.a -lm -ldl -lreadline -lpthread
./demo i10.aig
```

Build and run the MiniAIG/MiniLUT embedding example with:

```sh
gcc -Wall -O2 -Isrc -c src/runabc.c -o runabc.o
g++ -o runabc runabc.o libabc.a -lm -ldl -lreadline -lpthread
./runabc i10.aig
./runabc -t 4 i10.aig
```

The last command compares four concurrent sessions with a single-threaded
baseline. If ABC is built with `ABC_USE_NO_READLINE=1`, omit `-lreadline`.
Explicitly single-threaded configurations may disable thread-local frame
selection with `ABC_USE_NO_THREAD_LOCAL=1`.

To run the demo program, give it a file with the logic network in AIGER or BLIF. For example:

    [...] ~/abc> demo i10.aig
    i10          : i/o =  257/  224  lat =    0  and =   2396  lev = 37
    i10          : i/o =  257/  224  lat =    0  and =   1851  lev = 35
    Networks are equivalent.
    Reading =   0.00 sec   Rewriting =   0.18 sec   Verification =   0.41 sec

The same can be produced by running the binary in the command-line mode:

    [...] ~/abc> ./abc
    UC Berkeley, ABC 1.01 (compiled Oct  6 2012 19:05:18)
    abc 01> r i10.aig; b; ps; b; rw -l; rw -lz; b; rw -lz; b; ps; cec
    i10          : i/o =  257/  224  lat =    0  and =   2396  lev = 37
    i10          : i/o =  257/  224  lat =    0  and =   1851  lev = 35
    Networks are equivalent.

or in the batch mode:

    [...] ~/abc> ./abc -c "r i10.aig; b; ps; b; rw -l; rw -lz; b; rw -lz; b; ps; cec"
    ABC command line: "r i10.aig; b; ps; b; rw -l; rw -lz; b; rw -lz; b; ps; cec".
    i10          : i/o =  257/  224  lat =    0  and =   2396  lev = 37
    i10          : i/o =  257/  224  lat =    0  and =   1851  lev = 35
    Networks are equivalent.

### Compiling as C or C++

ABC can be compiled with a C or C++ compiler.

- To compile as C code, use the default Make configuration.
- To compile as C++ in a namespace, run `make ABC_USE_NAMESPACE=xxx`.

### Building a shared library

Build the shared library as position-independent code with:

```sh
make ABC_USE_PIC=1 libabc.so
```

### Building with CMake

ABC also supports CMake builds. The default CMake configuration builds the
executable and regression tests:

```sh
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

## Bug Reporting

Please reproduce bugs and unexpected behavior using the latest version of ABC
available from https://github.com/berkeley-abc/abc.

If the problem persists, provide the following information:

1. ABC version or Git commit.
2. Operating system, distribution, architecture, and compiler version.
3. The exact command line and complete error message.
4. The output of `ldd abc` on Linux, or the corresponding dependency report
   on another platform, when relevant.
5. Versions of relevant tools and libraries.

## Troubleshooting

1. If compilation does not start because of the cyclic dependency check, run
   `find . -type f -exec touch "{}" \;` and rebuild.
2. If readline is unavailable, install its development package or build with
   `make ABC_USE_NO_READLINE=1`.
3. If pthreads are unavailable, build without pthread support using
   `make ABC_USE_NO_PTHREADS=1`. Concurrent embedding requires pthread support.
4. On systems where readline depends on curses, include the appropriate curses
   library in `ABC_READLINE_LIBRARIES`.

The following tutorial is kindly offered by Ana Petkovska from EPFL:
https://www.dropbox.com/s/qrl9svlf0ylxy8p/ABC_GettingStarted.pdf

## Final Remarks

ABC includes CMake-based regression tests, which run in continuous integration
on supported configurations. The suite is not exhaustive, so bug reports
should include a small reproducer whenever possible.

This system is maintained by Alan Mishchenko <alanmi@berkeley.edu>. Consider also 
using ZZ framework developed by Niklas Een: https://bitbucket.org/niklaseen/abc-zz (or https://github.com/berkeley-abc/abc-zz)
