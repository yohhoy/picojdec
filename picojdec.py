#!/usr/bin/env python3
"""
pico baseline JPEG decoder

Copyright (c) 2017 yohhoy
"""
import math
import sys


# Marker symbols (X'FFxx')
MSYM = {'SOF0': 0xC0, 'SOF1': 0xC1, 'SOF2': 0xC2, 'SOF3': 0xC3, 'SOF5': 0xC5,
        'SOF6': 0xC6, 'SOF7': 0xC7, 'JPG': 0xC8, 'SOF9': 0xC9, 'SOF10': 0xCA,
        'SOF11': 0xCB, 'SOF13': 0xCD, 'SOF14': 0xCE, 'SOF15': 0xCF,
        'DHT': 0xC4, 'DAC': 0xCC, 'SOI': 0xD8, 'EOI': 0xD9, 'SOS': 0xDA,
        'DQT': 0xDB, 'DNL': 0xDC, 'DRI': 0xDD, 'DHP': 0xDE, 'EXP': 0xDF, 'COM': 0xFE}
MSYM.update({f'RST{m}': 0xD0 + m for m in range(8)})   # RST0..RST7
MSYM.update({f'APP{n}': 0xE0 + n for n in range(16)})  # APP0..APP15
MSYM.update({f'JPG{n}': 0xF0 + n for n in range(14)})  # APP0..APP13

# Zig-zag sequence
ZZ = [0,  1,  5,  6,  14, 15, 27, 28,
      2,  4,  7,  13, 16, 26, 29, 42,
      3,  8,  12, 17, 25, 30, 41, 43,
      9,  11, 18, 24, 31, 40, 44, 53,
      10, 19, 23, 32, 39, 45, 52, 54,
      20, 22, 33, 38, 46, 51, 55, 60,
      21, 34, 37, 47, 50, 56, 59, 61,
      35, 36, 48, 49, 57, 58, 62, 63]

# 8x8 Inverse DCT matrix
IDCT = [0.0] * 64
for i in range(64):
    y, x = i // 8, i % 8
    PI = math.pi / 16.0
    for j in range(64):
        v, u = j // 8, j % 8
        cu = 1.0 / math.sqrt(2) if u == 0 else 1.0
        cv = 1.0 / math.sqrt(2) if v == 0 else 1.0
        IDCT[i] += cu * cv * math.cos((2 * x + 1) * u * PI) * math.cos((2 * y + 1) * v * PI)
    IDCT[i] /= 4.0


class NoMoreData(Exception):
    pass

# byte/bit stream reader
class Reader():
    def __init__(self, filename):
        self.filename = self
        self.fs = open(filename, 'rb')
        self.blen = 0
        self.bbuf = 0
    # read n-bytes (as byte type)
    def byte_raw(self, n):
        b = self.fs.read(n)
        if len(b) < n:
            raise NoMoreData()
        return b
    # read n-bytes
    def byte(self, n):
        return int.from_bytes(self.byte_raw(n), 'big')
    # read n-bits
    def bits(self, n):
        ret = 0
        while 0 < n:
            if self.blen == 0:
                b = self.fs.read(1)
                if b == b'\xFF':
                    # X'FF00' -> 0xff(256)
                    stuff = self.fs.read(1)
                    if stuff != b'\x00':
                        raise NoMoreData()
                self.bbuf = int.from_bytes(b, 'big')
                self.blen = 8
            m = min(n, self.blen)
            lb = (self.bbuf >> (self.blen - m)) & ((1 << m) - 1)
            ret = (ret << m) | lb
            self.blen -= m
            n -= m
        return ret
    def __enter__(self):
        return self
    def __exit__(self, type, value, traceback):
        self.fs.close()


# Huffman decoder
class HuffmanDecoder():
    def __init__(self, hufftbl):
        self.huffval, self.huffsize, self.huffcode = hufftbl
    # decode single value
    def decode(self, r):
        code, sz = 0, 0
        for i, n in enumerate(self.huffsize):
            if sz < n:
                m, sz = n - sz, n
                code = (code << m) | r.bits(m)
            if self.huffcode[i] == code:
                print(f'Huffman: {code:0{sz}b} -> {self.huffval[i]}')
                return self.huffval[i]
        assert False, "broken Huffman code"


# decode HUFFVAL/HUFFSIZE/HUFFCODE
def decode_hufftable(v):
    huffval = []
    huffsize = []
    for i, w in enumerate(v):
        huffval.extend(w)
        huffsize += [i + 1] * len(w)
    huffsize += [0]

    huffcode = [0] * len(huffval)
    si = huffsize[0]
    code = 0
    k = 0
    while True:
        while True:
            huffcode[k] = code
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


# DC code table
def decode_dccode(r, hdec):
    ssss = hdec.decode(r)
    if ssss == 0:
        print('DC val=0')
        return 0
    b = r.bits(ssss)
    if b < (1 << (ssss - 1)):
        val = b - ((1 << ssss) - 1)
    else:
        val = b
    print(f'DC s={ssss} b={b:0{ssss}b} val={val:+d}')
    return val


# AC code table
def decode_accode(r, hdec):
    rs = hdec.decode(r)
    rrrr, ssss = rs >> 4, rs & 0b1111
    if ssss == 0:
        print('AC ' + 'EOB' if rrrr == 0 else 'ZRL')
        return (rrrr, 0)  # EOB/ZRL
    b = r.bits(ssss)
    if b < (1 << (ssss - 1)):
        val = b - ((1 << ssss) - 1)
    else:
        val = b
    print(f'AC r={rrrr} s={ssss} b={b:0{ssss}b} val={val:+d}')
    return (rrrr, val)


# inverse DCT
def idct(coeff):
    block = [0] * 64
    cos = lambda x: math.cos(x * math.pi / 16)
    for i in range(64):
        y, x = i // 8, i % 8
        s = 0
        for j in range(64):
            v, u = j // 8, j % 8
            cu = 1 / math.sqrt(2) if u == 0 else 1
            cv = 1 / math.sqrt(2) if v == 0 else 1
            s += cu * cv * coeff[j] * cos((2 * x + 1) * u) * cos((2 * y + 1) * v)
        block[i] = int(s / 4)
    return block


# Entropy-coded data segments
def decode_scan(r, image, scan):
    assert scan['SS'] == (0, 63), "SpectralSelection is not supported"
    assert scan['SA'] == (0, 0), "SuccessiveApproximation is not supported"
    # component list in MCU
    mcu_cs = []
    for sc in scan['C']:
        mcu_cs += [c for c in image['C'] if c['C'] == sc['Cs']]
    # component index of blocks in MCU
    mcu_ci = []
    for c in mcu_cs:
        mcu_ci += [c['i']] * (c['H'] * c['V'])
    # set QuantizationTables
    qtbl = []
    for c in mcu_cs:
        qtbl.append(image['QT'][c['Tq']])
    # initalize HuffmanDecoders
    hdec = []
    for c in scan['C']:
        hdec.append(HuffmanDecoder(image['HT'][c['Td'] * 2    ]))  # DC
        hdec.append(HuffmanDecoder(image['HT'][c['Ta'] * 2 + 1]))  # AC
    # DC predictors
    pred = [0] * len(scan)

    shift = 1 << (image['F']['bit'] - 1)
    maxval = (1 << image['F']['bit']) - 1
    rec = image['I']

    # decode MCU
    nmcu = image['F']['nmcu']
    for mcu_idx in range(nmcu[0] * nmcu[1]):
        mcu_y, mcu_x = mcu_idx // nmcu[0], mcu_idx % nmcu[0]
        ci = 0  # TODO: support multiple component
        sq = [0] * 64
        # decode DC coefficient
        sq[0] = decode_dccode(r, hdec[ci * 2]) + pred[ci]
        pred[ci] = sq[0]
        # decode AC coefficients (in Zig-Zag order)
        k = 1
        while k < 64:
            run, val = decode_accode(r, hdec[ci * 2 + 1])
            if (run, val) == (0, 0):  # EOB
                break
            k += run
            assert k < 64, f'out of range: k={k} run/val={run}/{val}'
            sq[k] = val
            k += 1
        print(f'{mcu_idx:4d}: Sq={sq[:k]}...')
        # dequantize
        for i in range(64):
            sq[i] *= qtbl[ci][i]
        print(f'{mcu_idx:4d}: Rz={sq[:k]}...')
        # inverse Zig-Zag scan
        coeff = [0] * 64
        for i, z in enumerate(ZZ):
            coeff[i] = sq[z]
        print(f'{mcu_idx:4d}: Ri=' + str([coeff[i:i + 8] for i in range(0, 64, 8)]))
        # inverse DCT
        block = idct(coeff)
        # level shift
        for i in range(64):
            block[i] = max(0, min(block[i] + shift, maxval))
        print(f'{mcu_idx:4d}:  I=' + str([block[i:i + 8] for i in range(0, 64, 8)]))
        # dump PGM
        for i in range(64):
            rec_x = mcu_x * 8 + i % 8
            rec_y = mcu_y * 8 + i // 8
            rec[rec_y * image['F']['size'][0] + rec_x] = block[i]
        if mcu_x == nmcu[0] - 1:
            with open(f'mcu{mcu_idx:04d}.pgm', 'wb') as f:
                f.write('P5\n{0[0]} {0[1]}\n255\n'.format(image['F']['size']).encode('ascii'))
                for px in rec:
                    f.write(px.to_bytes(1, 'big'))


# Frame header
def parse_SOFn(r, n, image):
    lf = r.byte(2)
    p = r.byte(1)
    y = r.byte(2)
    x = r.byte(2)
    nf = r.byte(1)
    print(f'SOF{n}: Lf={lf} P={p} Y={y} X={x} Nf={nf}')
    component = [{}] * nf
    for i in range(nf):
        c = r.byte(1)
        b = r.byte(1)
        h, v = b >> 4, b & 0b1111
        tq = r.byte(1)
        component[i] = {'i': i, 'C': c, 'H': h, 'V': v, 'Tq': tq}
        print(f'  Component[{i}]: C={c} H={h} V={v} Tq={tq}')
    lf -= 8 + 3 * nf
    assert lf == 0, "invalid SOF payload"
    # check QuantizationTable
    for c in component:
        assert image['QT'][c['Tq']], "QuantizationTable not found"
    # register frame/component
    mcu_h = max(c['H'] * 8 for c in component)
    mcu_v = max(c['V'] * 8 for c in component)
    nmcu = (math.ceil(x / mcu_h), math.ceil(y / mcu_v))
    image['F'] = {'bit': p, 'size': (x, y), 'nmcu': nmcu}
    image['C'] = component
    image['I'] = [0] * (mcu_h * nmcu[0] * mcu_v * nmcu[1])


# Scan header
def parse_SOS(r, image):
    ls = r.byte(2)
    ns = r.byte(1)
    print(f'SOS: Ls={ls} Ns={ns}')
    component = [{}] * ns
    for j in range(ns):
        cs = r.byte(1)
        b = r.byte(1)
        td, ta = b >> 4, b & 0b1111
        component[j] = {'Cs': cs, 'Td': td, 'Ta': ta}
        print(f'  Component[{j}]: Cs={cs} Td={td} Ta={ta}')
    ss = r.byte(1)
    se = r.byte(1)
    b = r.byte(1)
    ah, al = b >> 4, b & 0b1111
    ls -= 6 + 2 * ns
    print(f'  Ss={ss} Se={se} Ah={ah} Al={al}')
    assert ls == 0, "invalid SOS payload"
    # check HuffmanTables
    for sc in component:
        assert image['HT'][sc['Td'] * 2    ], "HuffmanTable/DC does not found"
        assert image['HT'][sc['Ta'] * 2 + 1], "HuffmanTable/AC does not found"
    # SS=spectral selection, SA=successive approximation
    return {'SS': (ss, se), 'SA': (ah, al), 'C': component}


# Quantization table-specification
def parse_DQT(r, image):
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
    # register QuantizationTable
    image['QT'][tq] = q


# Huffman table-specification
def parse_DHT(r, image):
    lh = r.byte(2)
    print(f'DHT: Lh={lh}')
    lh -= 2
    t = 0
    while 0 < lh:
        b = r.byte(1)
        tc, th = b >> 4, b & 0b1111
        l = [0] * 16  # BITS
        for i in range(16):
            l[i] = r.byte(1)
        v = [[]] * 16  # HUFFVAL
        for i in range(16):
            if l[i] == 0:
                continue
            v[i] = [0] * l[i]
            for j in range(l[i]):
                v[i][j] = r.byte(1)
        print(f'  DHT[{t}]: Tc={tc} Th={th} L_i={l} V_ij={v}')
        lh -= 1 + 16 + sum(l)
        # decode HuffmanTable
        hv, hs, hc = decode_hufftable(v)
        if tc == 0:
            hv_s = ', '.join([str(v) for v in hv])
        else:
            def achv2str(rs):
                r, s = rs >> 4, rs & 0b1111
                return ('ZRL' if r == 0 else f'EOB{r}') if s == 0 else f'{r}/{s}'
            hv_s = ', '.join([achv2str(v) for v in hv])
        hc_s = ', '.join(['{n:0{s}b}'.format(n=c, s=hs[i]) for i, c in enumerate(hc)])
        print(f'  DHT[{t}]: HUFFVAL=[{hv_s}]')
        print(f'  DHT[{t}]: HUFFCODE=[{hc_s}]')
        # register HuffmanTable
        image['HT'][th * 2 + tc] = (hv, hs, hc)
        t += 1
    assert lh == 0, "invalid DHT payload"


# Application data
def parse_APPn(r, n):
    la = r.byte(2)
    if n == 0:
        identifier = r.byte_raw(5)
        if identifier == b'JFIF\x00':
            version = r.byte(2)
            units = r.byte(1)
            h_dens = r.byte(2)
            v_dens = r.byte(2)
            h_thumb = r.byte(1)
            v_thumb = r.byte(1)
            print(f'APP0: La={la} {identifier} version={version:04X}'
                  f' units={units} density={h_dens},{v_dens} thumbnail={h_thumb},{v_thumb}')
        else:
            print(f'APP0: La={la} {identifier} ...')
    else:
        ap = r.byte_raw(la - 2)
        print(f'APP{n}: La={la} Ap={ap}')


def parse_stream(r):
    M = {v: k for k, v in MSYM.items()}  # X'FFxx' -> marker symbol
    try:
        image = None
        while True:
            # search marker
            b = r.byte(1)
            if b != 0xFF:
                continue
            b = r.byte(1)
            if b == 0x00:
                continue
            m = M.get(b, None)
            # parse marker
            if m == 'SOI':
                # 'Start of image' marker
                print(f'{m}')
                image = {'QT': [None] * 4,  # QuantizationTable (0..3)
                         'HT': [None] * 8}  # HuffmanTable (DC/AC, 0..3)
            elif m == 'EOI':
                # 'End of image' marker
                print(f'{m}')
                image = None
            elif m[:3] == 'SOF':
                # 'Start of frame' markers
                n = b - MSYM['SOF0']
                frame = parse_SOFn(r, n, image)
                #assert n == 0, "support only SOF0/Baseline DCT"
            elif m == 'DQT':
                # 'Define quantization tables' marker
                parse_DQT(r, image)
            elif m == 'DHT':
                # 'Define Huffman tables' marker
                parse_DHT(r, image)
            elif m == 'DAC':
                # 'Define arithmetic coding conditionings' marker
                assert False, "Arithmetic coding is not supported"
            elif m == 'SOS':
                # 'Start of scan' marker
                scan = parse_SOS(r, image)
                decode_scan(r, image, scan)
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
    args = sys.argv[1:]
    if len(args) != 1:
        print('usage: picojdec.py <input.jpg>')
        exit(1)
    main(args[0])
