#!/usr/bin/env python3
"""Extract requested variant positions from a GWAS summary-statistics table.

Usage: extract_positions.py <sumstats.tsv.gz> <keys.txt> <trait> <out.tsv>

The table may be plain or gzipped and must have a header row. Column names are
matched against common GWAS Catalog spellings, so both the harmonised
(GCST...) layout and the original-study layout are accepted. Output columns
are trait, chr, pos, ea, nea, beta, se, eaf, p -- the layout used by the R
scripts.
"""
import gzip
import sys

CHR = ("chromosome", "chr", "chrom")
POS = ("base_pair_location", "pos", "position", "bp")
EA = ("effect_allele", "ea", "alt")
NEA = ("other_allele", "nea", "oa", "reference_allele", "ref")
BETA = ("beta", "effect", "b")
SE = ("standard_error", "se", "stderr")
EAF = ("effect_allele_frequency", "eaf", "af", "alt_af", "freq")
P = ("p_value", "pval", "p", "p_value_meta")


def pick(header, names):
    for n in names:
        for i, h in enumerate(header):
            if h.strip().lower() == n:
                return i
    return None


def main():
    path, keyfile, trait, outpath = sys.argv[1:5]
    keys = {}
    with open(keyfile) as k:
        for line in k:
            line = line.strip()
            if line:
                c, p = line.split(":")
                keys.setdefault(c, set()).add(int(p))

    opener = gzip.open if path.endswith(".gz") else open
    n = 0
    with opener(path, "rt") as f, open(outpath, "w") as w:
        header = f.readline().rstrip("\n").split("\t")
        ic, ip = pick(header, CHR), pick(header, POS)
        ie, inea, ib = pick(header, EA), pick(header, NEA), pick(header, BETA)
        ise, iaf, ipv = pick(header, SE), pick(header, EAF), pick(header, P)
        if None in (ic, ip, ie, inea, ib, ise, ipv):
            sys.exit("could not map the required columns in " + path)
        w.write("trait\tchr\tpos\tea\tnea\tbeta\tse\teaf\tp\n")
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(ic, ip, ie, inea, ib, ise, ipv):
                continue
            try:
                p = int(float(parts[ip]))
            except ValueError:
                continue
            c = parts[ic].replace("chr", "")
            if p not in keys.get(c, ()):
                continue
            eaf = parts[iaf] if iaf is not None and iaf < len(parts) else "NA"
            w.write("\t".join([trait, c, str(p), parts[ie], parts[inea], parts[ib],
                               parts[ise], eaf, parts[ipv]]) + "\n")
            n += 1
    print(f"{trait}: {n} variants written to {outpath}")


if __name__ == "__main__":
    main()
