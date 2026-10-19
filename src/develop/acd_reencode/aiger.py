"""Minimal AIGER reader/writer (binary and ASCII; latches supported, init 0 only)."""


class Aiger:
    """Variables: 0 = const, 1..I inputs, I+1..I+L latches, then ANDs (topological)."""

    def __init__(self):
        self.I = 0
        self.latches = []     # next-state literals, one per latch
        self.ands = []        # (rhs0, rhs1) for variable I+L+1+k
        self.outputs = []     # literals
        self.in_names, self.out_names = [], []

    @property
    def L(self):
        return len(self.latches)

    def and_var(self, k):
        return self.I + self.L + 1 + k


def _uint(data, pos):
    x, i = 0, 0
    while True:
        b = data[pos]; pos += 1
        x |= (b & 0x7F) << (7 * i); i += 1
        if not b & 0x80:
            return x, pos


def read(path):
    data = open(path, "rb").read()
    nl = data.index(b"\n")
    hdr = data[:nl].split()
    fmt = hdr[0].decode()
    M, I, L, O, A = (int(x) for x in hdr[1:6])
    if len(hdr) > 6 and any(int(x) for x in hdr[6:]):
        raise ValueError("B/C/J/F sections not supported")
    g = Aiger(); g.I = I
    pos = nl + 1

    def line():
        nonlocal pos
        e = data.index(b"\n", pos); s = data[pos:e].decode(); pos = e + 1
        return s.split()

    if fmt == "aag":
        lits = [int(line()[0]) for _ in range(I)]
        if lits != [2 * (i + 1) for i in range(I)]:
            raise ValueError("ASCII inputs must be 2,4,...")
    for k in range(L):
        f = line()
        if fmt == "aag":
            f = f[1:]
        if len(f) > 1 and int(f[1]) != 0:
            raise ValueError("only zero-initialized latches are supported")
        g.latches.append(int(f[0]))
    g.outputs = [int(line()[0]) for _ in range(O)]
    if fmt == "aig":
        for k in range(A):
            lhs = 2 * (I + L + 1 + k)
            d0, pos = _uint(data, pos)
            d1, pos = _uint(data, pos)
            g.ands.append((lhs - d0, lhs - d0 - d1))
    else:
        defs = {}
        for _ in range(A):
            l, a, b = (int(x) for x in line())
            defs[l] = (a, b)
        for k in range(A):
            g.ands.append(defs[2 * (I + L + 1 + k)])
    g.in_names = [f"pi{i}" for i in range(I)]
    g.out_names = [f"po{i}" for i in range(O)]
    while pos < len(data):
        e = data.find(b"\n", pos)
        e = len(data) if e < 0 else e
        s = data[pos:e].decode(errors="replace"); pos = e + 1
        if not s or s[0] == "c":
            break
        if s[0] in "io" and " " in s:
            idx, name = s[1:].split(" ", 1)
            (g.in_names if s[0] == "i" else g.out_names)[int(idx)] = name
    return g


def _enc(x):
    out = bytearray()
    while x & ~0x7F:
        out.append((x & 0x7F) | 0x80); x >>= 7
    out.append(x)
    return bytes(out)


def write(g, path, comment=""):
    """Binary AIGER. ANDs must be topological with rhs0 >= rhs1 not required (normalized here)."""
    M = g.I + g.L + len(g.ands)
    out = [f"aig {M} {g.I} {g.L} {len(g.outputs)} {len(g.ands)}\n".encode()]
    out += [f"{n}\n".encode() for n in g.latches]
    out += [f"{o}\n".encode() for o in g.outputs]
    body = bytearray()
    for k, (a, b) in enumerate(g.ands):
        lhs = 2 * g.and_var(k)
        a, b = max(a, b), min(a, b)
        assert lhs > a, "ANDs are not in topological order"
        body += _enc(lhs - a) + _enc(a - b)
    out.append(bytes(body))
    out += [f"i{i} {n}\n".encode() for i, n in enumerate(g.in_names)]
    out += [f"o{i} {n}\n".encode() for i, n in enumerate(g.out_names)]
    if comment:
        out.append(b"c\n" + comment.encode() + b"\n")
    open(path, "wb").write(b"".join(out))
