#!/usr/bin/env python3
"""Small hand-made / random AIGs (ASCII AIGER) for testing aig_mincut.

usage: gen_aig.py adder N | parity N | compare N | random PIS ANDS POS SEED   > out.aag
"""
import random
import sys


class Builder:
    def __init__(self, n_pi):
        self.n_pi = n_pi
        self.ands = []          # (lhs, rhs0, rhs1)
        self.outs = []

    def pi(self, i):
        return 2 * (i + 1)

    def AND(self, a, b):
        lhs = 2 * (self.n_pi + len(self.ands) + 1)
        self.ands.append((lhs, a, b))
        return lhs

    def OR(self, a, b):
        return self.AND(a ^ 1, b ^ 1) ^ 1

    def XOR(self, a, b):
        return self.OR(self.AND(a, b ^ 1), self.AND(a ^ 1, b))

    def dump(self):
        M = self.n_pi + len(self.ands)
        lines = [f"aag {M} {self.n_pi} 0 {len(self.outs)} {len(self.ands)}"]
        lines += [str(self.pi(i)) for i in range(self.n_pi)]
        lines += [str(o) for o in self.outs]
        lines += [f"{l} {a} {b}" for l, a, b in self.ands]
        return "\n".join(lines) + "\n"


def adder(n):
    b = Builder(2 * n)
    c = 0
    for i in range(n):
        x, y = b.pi(i), b.pi(n + i)
        s = b.XOR(b.XOR(x, y), c) if c else b.XOR(x, y)
        c = b.OR(b.AND(x, y), b.AND(c, b.XOR(x, y))) if c else b.AND(x, y)
        b.outs.append(s)
    b.outs.append(c)
    return b


def parity(n):
    b = Builder(n)
    layer = [b.pi(i) for i in range(n)]
    while len(layer) > 1:
        nxt = [b.XOR(layer[i], layer[i + 1]) for i in range(0, len(layer) - 1, 2)]
        if len(layer) % 2:
            nxt.append(layer[-1])
        layer = nxt
    b.outs.append(layer[0])
    return b


def compare(n):
    """a > b, MSB first, plus a == b."""
    b = Builder(2 * n)
    gt, eq = 0, 1
    for i in reversed(range(n)):
        x, y = b.pi(i), b.pi(n + i)
        gt = b.OR(gt, b.AND(eq, b.AND(x, y ^ 1))) if gt else b.AND(x, y ^ 1)
        e = b.XOR(x, y) ^ 1
        eq = b.AND(eq, e) if eq != 1 else e
    b.outs += [gt, eq]
    return b


def rand_aig(n_pi, n_and, n_po, seed):
    rng = random.Random(seed)
    b = Builder(n_pi)
    lits = [b.pi(i) for i in range(n_pi)]
    for _ in range(n_and):
        x, y = rng.sample(lits, 2) if len(lits) > 1 else (lits[0], lits[0])
        lits.append(b.AND(x ^ rng.randint(0, 1), y ^ rng.randint(0, 1)))
    # last AND is always an output so nothing is trivially dangling at the top
    b.outs.append(lits[-1] ^ rng.randint(0, 1))
    for _ in range(n_po - 1):
        b.outs.append(rng.choice(lits) ^ rng.randint(0, 1))
    return b


if __name__ == "__main__":
    kind, args = sys.argv[1], [int(a) for a in sys.argv[2:]]
    gen = {"adder": adder, "parity": parity, "compare": compare, "random": rand_aig}[kind]
    sys.stdout.write(gen(*args).dump())
