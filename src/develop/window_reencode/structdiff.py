#!/usr/bin/env python3
"""Structural comparison of two combinational AIGs with the same PIs (binary AIGER).

Every node is simulated on the same 2048 random PI patterns; a node's function key is its signature
up to complement. Comparing B against a reference A:
  new_frac   fraction of B's ANDs whose function does not occur anywhere in A (PIs, const and ANDs)
  kept_frac  fraction of A's ANDs whose function still occurs in B
Plus per-AIG structure statistics: ANDs, depth, structural XOR/XNOR and MUX counts (3-AND patterns),
mean fanout, fraction of ANDs with fanout >= 2, and the mean normalized level (level / depth) of B's new
nodes, i.e. where in the circuit the new functions sit.
usage: structdiff.py A.aig B.aig [B2.aig ...]
"""
import random, sys

NBITS = 2048
MASK = (1 << NBITS) - 1


class Aig:
    def __init__(self, path):
        data = open(path, "rb").read()
        nl = data.index(b"\n")
        h = data[:nl].split()
        if h[0] != b"aig":
            raise ValueError(f"{path}: not a binary AIGER file")
        M, I, L, O, A = map(int, h[1:6])
        if L:
            raise ValueError(f"{path}: latches are not supported")
        pos = nl + 1
        self.outs = []
        for _ in range(O):
            e = data.index(b"\n", pos)
            self.outs.append(int(data[pos:e]))
            pos = e + 1
        self.nPi, self.nAnd = I, A
        self.f0, self.f1 = [0] * (I + A + 1), [0] * (I + A + 1)

        def varint():
            nonlocal pos
            x = sh = 0
            while True:
                b = data[pos]
                pos += 1
                x |= (b & 0x7F) << sh
                if not b & 0x80:
                    return x
                sh += 7

        for i in range(A):
            v = I + 1 + i
            r0 = 2 * v - varint()
            r1 = r0 - varint()
            self.f0[v], self.f1[v] = r0, r1
        n = I + A + 1
        self.lev = [0] * n
        self.fo = [0] * n
        for v in range(I + 1, n):
            self.lev[v] = 1 + max(self.lev[self.f0[v] >> 1], self.lev[self.f1[v] >> 1])
            self.fo[self.f0[v] >> 1] += 1
            self.fo[self.f1[v] >> 1] += 1
        for o in self.outs:
            self.fo[o >> 1] += 1
        self.depth = max((self.lev[o >> 1] for o in self.outs), default=0)

    def simulate(self, seed=12345):
        rng = random.Random(seed)
        n = self.nPi + self.nAnd + 1
        val = [0] * n
        for i in range(1, self.nPi + 1):
            val[i] = rng.getrandbits(NBITS)
        for v in range(self.nPi + 1, n):
            a, b = self.f0[v], self.f1[v]
            x = val[a >> 1] ^ (MASK if a & 1 else 0)
            y = val[b >> 1] ^ (MASK if b & 1 else 0)
            val[v] = x & y
        self.keys = [w ^ MASK if w & 1 else w for w in val]   # canonical up to complement
        return self.keys

    def patterns(self):
        """Structural XOR/XNOR and MUX roots: v = AND(!p, !q) with p, q ANDs sharing a variable in
        opposite polarity (MUX); XOR if both shared variables appear in opposite polarities."""
        nx = nm = 0
        for v in range(self.nPi + 1, self.nPi + self.nAnd + 1):
            a, b = self.f0[v], self.f1[v]
            if not (a & 1 and b & 1):
                continue
            p, q = a >> 1, b >> 1
            if p <= self.nPi or q <= self.nPi:
                continue
            P, Q = {self.f0[p], self.f1[p]}, {self.f0[q], self.f1[q]}
            opp = [l for l in P if l ^ 1 in Q]
            if len(opp) == 2:
                nx += 1
            elif len(opp) == 1:
                nm += 1
        return nx, nm

    def stats(self):
        ands = range(self.nPi + 1, self.nPi + self.nAnd + 1)
        nx, nm = self.patterns()
        return dict(ands=self.nAnd, depth=self.depth, xor=nx, mux=nm,
                    mean_fo=round(sum(self.fo[v] for v in ands) / max(1, self.nAnd), 3),
                    frac_shared=round(sum(self.fo[v] >= 2 for v in ands) / max(1, self.nAnd), 3))


def compare(A, B):
    """A, B: simulated Aig objects with the same PIs."""
    ka = set(A.keys)
    kb = set(B.keys)
    andsB = range(B.nPi + 1, B.nPi + B.nAnd + 1)
    newB = [v for v in andsB if B.keys[v] not in ka]
    keptA = sum(A.keys[v] in kb for v in range(A.nPi + 1, A.nPi + A.nAnd + 1))
    return dict(new_frac=round(len(newB) / max(1, B.nAnd), 4),
                kept_frac=round(keptA / max(1, A.nAnd), 4),
                n_new=len(newB),
                new_lev=round(sum(B.lev[v] for v in newB) / max(1, len(newB)) / max(1, B.depth), 3) if newB else "")


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 1
    A = Aig(sys.argv[1])
    A.simulate()
    print("A", sys.argv[1], A.stats())
    for p in sys.argv[2:]:
        B = Aig(p)
        if B.nPi != A.nPi:
            print(p, "different PI count")
            continue
        B.simulate()
        print("B", p, B.stats(), compare(A, B))


if __name__ == "__main__":
    sys.exit(main())
