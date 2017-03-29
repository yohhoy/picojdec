#!/usr/bin/env python3
"""
pico baseline JPEG decoder

Copyright(c) 2016 yohhoy
"""
import math
import sys


# Marker symbols (X'FFxx')
MSYM = {'SOF0': 0xC0, 'SOF1': 0xC1, 'SOF2': 0xC2, 'SOF3': 0xC3, 'SOF5': 0xC5,
        'SOF6': 0xC6, 'SOF7': 0xC7, 'JPG': 0xC8, 'SOF9': 0xC9, 'SOF10': 0xCA,
        'SOF11': 0xCB, 'SOF13': 0xCD, 'SOF14': 0xCE, 'SOF15': 0xCF,
        'DHT': 0xC4, 'DAC': 0xCC, 'SOI': 0xD8, 'EOI': 0xD9, 'SOS': 0xDA,
        'DQT': 0xDB, 'DNL': 0xDC, 'DRI': 0xDD, 'DHP': 0xDE, 'EXP': 0xDF, 'COM': 0xFE}
MSYM.update(dict([(f'RST{m}', 0xD0 + m) for m in range(8)]))   # RST0..RST7
MSYM.update(dict([(f'APP{n}', 0xE0 + n) for n in range(16)]))  # APP0..APP15
MSYM.update(dict([(f'JPG{n}', 0xF0 + n) for n in range(14)]))  # APP0..APP13

# Zig-zag sequence
ZZ = [ 0,  1,  5,  6, 14, 15, 27, 28,
       2,  4,  7, 13, 16, 26, 29, 42,
       3,  8, 12, 17, 25, 30, 41, 43,
       9, 11, 18, 24, 31, 40, 44, 53,
      10, 19, 23, 32, 39, 45, 52, 54,
      20, 22, 33, 38, 46, 51, 55, 60,
      21, 34, 37, 47, 50, 56, 59, 61,
      35, 36, 48, 49, 57, 58, 62, 63]

# 8x8 IDCT matrix
IDCT = [0.0] * 64
for i in range(64):
    y, x = i / 8, i % 8
    PI = math.pi / 16.0
    for j in range(64):
        u, v = j / 8, j % 8
        cu = 1.0 / math.sqrt(2) if u == 0 else 1.0
        cv = 1.0 / math.sqrt(2) if v == 0 else 1.0
        IDCT[i] += cu * cv * math.cos((2 * x + 1) * u * PI) * math.cos((2 * y + 1) * v * PI)
    IDCT[i] /= 4.0


class NoMoreData(Exception):
    pass

class Reader():
    def __init__(self, filename):
        self.filename = self
        self.fs = open(filename, 'rb')
    def byte_raw(self, n = 1):
        b = self.fs.read(n)
        if len(b) < n:
            raise NoMoreData()
        return b
    def byte(self, n = 1):
        return int.from_bytes(self.byte_raw(n), 'big')
    def __enter__(self):
        return self
    def __exit__(self, type, value, traceback):
        self.fs.close()


def decode_hufftable(l, v):
    huffval = []
    huffsize = []
    for i, vi in enumerate(v):
        huffval.extend(vi)
        huffsize += [i + 1] * len(vi)
    huffsize += [0]

    huffcode = []
    si = huffsize[0]
    code = 0
    k = 0
    while True:
        while True:
            huffcode.append(code)
            k += 1
            code += 1
            if huffsize[k] != si:
                break
        if huffsize[k] == 0:
            break
        while True:
            code <<= 1
            si += 1
            if huffsize[k] == si:
                break
    return (huffval, huffsize, huffcode)


def parse_SOFn(r, n):
    lf = r.byte(2)
    p = r.byte(1)
    y = r.byte(2)
    x = r.byte(2)
    nf = r.byte(1)
    print(f'SOF{n}: Lf={lf} P={p} Y={y} X={x} Nf={nf}')
    for i in range(nf):
        c = r.byte(1)
        b = r.byte(1)
        h, v = b >> 4, b & 0b1111
        tq = r.byte(1)
        print(f'  Component[{i}]: C={c} H={h} V={v} Tq={tq}')
    lf -= 8 + 3 * nf
    assert lf == 0, "invalid SOF payload"


def parse_SOS(r):
    ls = r.byte(2)
    ns = r.byte(1)
    print(f'SOS: Ls={ls} Ns={ns}')
    for j in range(ns):
        cs = r.byte(1)
        b = r.byte(1)
        td, ta = b >> 4, b & 0b1111
        print(f'  Component[{j}]: Cs={cs} Td={td} Ta={ta}')
    ss = r.byte(1)
    se = r.byte(1)
    b = r.byte(1)
    ah, al = b >> 4, b & 0b1111
    ls -= 6 + 2 * ns
    print(f'  Ss={ss} Se={se} Ah={ah} Al={al}')
    assert ls == 0, "invalid SOS payload"


def parse_DQT(r):
    lq = r.byte(2)
    print(f'DQT: Lq={lq}')
    lq -= 2
    t = 0
    while 0 < lq:
        b = r.byte(1)
        pq, tq = b >> 4, b & 0b1111
        n = 1 if pq == 0 else 2
        q = [0] * 64
        for k in range(64):
            q[k] = r.byte(n)
        print(f'  DQT[{t}]: Pq={pq} Tq={tq} Q_k={q}')
        lq -= 1 + 64 * n
        t += 1
    assert lq == 0, "invalid DQT payload"


def parse_DHT(r):
    lh = r.byte(2)
    print(f'DHT: Lh={lh}')
    lh -= 2
    t = 0
    while 0 < lh:
        b = r.byte(1)
        tc, th = b >> 4, b & 0b1111
        l = [0] * 16
        for i in range(16):
            l[i] = r.byte(1)
        v = [[]] * 16
        for i in range(16):
            if l[i] == 0:
                continue
            v[i] = [0] * l[i]
            for j in range(l[i]):
                v[i][j] = r.byte(1)
        print(f'  DHT[{t}]: Tc={tc} Th={th} L_i={l} V_ij={v}')
        lh -= 1 + 16 + sum(l)
        # decode Huffman table
        hv, hs, hc = decode_hufftable(l, v)
        if tc == 0:
            hv_s = ', '.join([str(v) for v in hv])
        else:
            def achv2str(v):
                rrrr, ssss = v >> 4, v & 0b1111
                return ('ZRL' if rrrr == 0 else f'EOB{rrrr}') if ssss == 0 else f'{rrrr}/{ssss}'
            hv_s = ', '.join([achv2str(v) for v in hv])
        hc_s = ', '.join(['{n:0{s}b}'.format(n=c, s=hs[i]) for i, c in enumerate(hc)])
        print(f'  DHT[{t}]: HUFFVAL=[{hv_s}]')
        print(f'  DHT[{t}]: HUFFCODE=[{hc_s}]')
        t += 1
    assert lh == 0, "invalid DHT payload"


def parse_APPn(r, n):
    la = r.byte(2)
    ap = r.byte_raw(la - 2)
    print(f'APP{n}: La={la} Ap={ap}')


def parse_stream(r):
    M = {v: k for k, v in MSYM.items()}  # X'FFxx' -> marker symbol
    try:
        while True:
            # search marker
            b = r.byte()
            if b != 0xFF:
                continue
            b = r.byte()
            if b == 0x00:
                continue
            m = M.get(b, None)
            # parse marker
            if m == 'SOI' or m == 'EOI':
                # 'Start Of Image', 'End Image' markers
                print(f'{m}')
            elif m[:3] == 'SOF':
                # 'Start Of Frame' markers
                n = b - MSYM['SOF0']
                parse_SOFn(r, n)
                #assert n == 0, "support only SOF0/Baseline DCT"
            elif m == 'DQT':
                # 'Define quantization tables' markers
                parse_DQT(r)
            elif m == 'DHT':
                # 'Define Huffman tables' markers
                parse_DHT(r)
            elif m == 'DAC':
                # 'Define arithmetic coding conditionings' markers
                assert False, "Arithmetic coding is not supported"
            elif m == 'SOS':
                # 'Start of scan' markers
                parse_SOS(r)
            elif m[:3] == 'APP':
                # 'Reserved for application segments' markers
                n = b - MSYM['APP0']
                parse_APPn(r, n)
            else:
                print(f'(ignore {m})')
    except NoMoreData:
        pass


def main(filename):
    with Reader(filename) as r:
        parse_stream(r)


if __name__ == '__main__':
    args = sys.argv[1:];
    if len(args) != 1:
        print('usage: picojdec.py <input.jpg>')
        exit(1)
    main(args[0])
