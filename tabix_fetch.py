#!/usr/bin/env python3
"""Fetch summary statistics at specific positions using the tabix linear index.

Uses small parallel HTTP range requests around each locus, decompressing BGZF
blocks directly. Far cheaper than streaming a multi-GB file.
"""
import sys, struct, zlib, subprocess, urllib.request
from concurrent.futures import ThreadPoolExecutor

UA = {"User-Agent": "Mozilla/5.0"}
WINDOW = 131072          # bytes fetched per locus (~0.8 Mb of genome)
THREADS = 16


def get_range(url, start, end, tries=4):
    for t in range(tries):
        try:
            req = urllib.request.Request(url, headers={**UA, "Range": f"bytes={start}-{end}"})
            with urllib.request.urlopen(req, timeout=90) as r:
                return r.read()
        except Exception:
            if t == tries - 1:
                return b""
    return b""


def parse_tbi(data):
    if data[:2] == b"\x1f\x8b":
        data = subprocess.run(["bgzip", "-dc"], input=data, capture_output=True).stdout
    assert data[:4] == b"TBI\x01"
    off = 4
    n_ref, = struct.unpack_from("<i", data, off); off += 4
    fmt, col_seq, col_beg, col_end, meta, skip, l_nm = struct.unpack_from("<7i", data, off); off += 28
    names = data[off:off + l_nm].split(b"\x00")[:-1]
    off += l_nm
    refs = []
    for _ in range(n_ref):
        n_bin, = struct.unpack_from("<i", data, off); off += 4
        for _ in range(n_bin):
            bin_id, n_chunk = struct.unpack_from("<Ii", data, off); off += 8
            off += 16 * n_chunk
        n_intv, = struct.unpack_from("<i", data, off); off += 4
        intv = list(struct.unpack_from("<%dQ" % n_intv, data, off)); off += 8 * n_intv
        refs.append(intv)
    return [n.decode() for n in names], refs


def decompress_blocks(data):
    out = []
    pos = 0
    while pos + 18 <= len(data):
        try:
            bsize, = struct.unpack_from("<H", data, pos + 16)
        except struct.error:
            break
        total = bsize + 1
        if total <= 18 or pos + total > len(data):
            break
        try:
            d = zlib.decompressobj(-15)
            out.append(d.decompress(data[pos + 18:pos + total]).decode("utf-8", "replace"))
        except Exception:
            pass
        pos += total
    return "".join(out)


def detect_cols(url):
    head = decompress_blocks(get_range(url, 0, 300000))
    hdr = head.split("\n")[0].split("\t")
    idx = {}
    for name in ["af_EUR", "beta_EUR", "se_EUR", "neglog10_pval_EUR"]:
        idx[name] = hdr.index(name)
    return idx


def main():
    url, keyfile, outpath, trait = sys.argv[1:5]
    global AF, BETA, SE, PCOL
    ci = detect_cols(url)
    AF, BETA, SE, PCOL = ci["af_EUR"], ci["beta_EUR"], ci["se_EUR"], ci["neglog10_pval_EUR"]
    keys = {}
    for line in open(keyfile):
        line = line.strip()
        if line:
            c, p = line.split(":")
            keys.setdefault(c, set()).add(int(p))
    names, refs = parse_tbi(get_range(url + ".tbi", 0, 10 ** 9))
    name2intv = dict(zip(names, refs))

    tasks = {}   # (chrom, coffset) -> (start_byte, end_byte)
    for chrom, positions in keys.items():
        intv = name2intv.get(chrom)
        if not intv:
            continue
        for p in sorted(positions):
            i = p >> 14
            if i >= len(intv):
                continue
            voff = intv[i]
            while voff == 0 and i > 0:
                i -= 1
                voff = intv[i]
            coff = (voff >> 16) if voff else 0
            tasks.setdefault((chrom, coff), [coff, coff + WINDOW])

    def work(item):
        (chrom, base), (s, e) = item
        data = get_range(url, s, e)
        return chrom, decompress_blocks(data)

    found = {}
    with ThreadPoolExecutor(THREADS) as ex:
        for chrom, text in ex.map(work, tasks.items()):
            positions = keys[chrom]
            for line in text.split("\n"):
                parts = line.split("\t")
                if len(parts) < 8:
                    continue
                try:
                    p = int(parts[1])
                except ValueError:
                    continue
                if parts[0] == chrom and p in positions:
                    found[(chrom, p)] = parts

    with open(outpath, "w") as w:
        w.write("trait\tchr\tpos\tea\tnea\tbeta\tse\teaf\tp\n")
        for p in found.values():
            pv = p[PCOL] if len(p) > PCOL else "NA"
            try:
                pv = 10 ** (-float(pv))
            except ValueError:
                pv = "NA"
            w.write("\t".join([trait, p[0], p[1], p[3], p[2], p[BETA], p[SE], p[AF], str(pv)]) + "\n")
    print(f"{trait}: {len(found)}/{sum(len(v) for v in keys.values())} found "
          f"({len(tasks)} fetches)", flush=True)


if __name__ == "__main__":
    layout = sys.argv[5] if len(sys.argv) > 5 else "cont"
    if layout == "bio":
        AF, BETA, SE, PCOL = 13, 19, 25, 31      # 0-based indices for biomarkers files
    else:
        AF, BETA, SE, PCOL = 18, 24, 30, 36
    main()
