CXX      ?= g++
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra

acd_phase1: phase1.cpp ../aig_mincut/aig.cpp ../aig_mincut/aig.h
	$(CXX) $(CXXFLAGS) -o $@ phase1.cpp ../aig_mincut/aig.cpp

clean:
	rm -f acd_phase1

.PHONY: clean
