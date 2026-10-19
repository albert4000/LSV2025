CXX      ?= g++
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra

aig_mincut: mincut.cpp aig.cpp aig.h flow.h
	$(CXX) $(CXXFLAGS) -o $@ mincut.cpp aig.cpp

clean:
	rm -f aig_mincut

.PHONY: clean
