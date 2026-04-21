"""Build amplicons for data/variants.json using:
  1) Primer-based in-silico PCR to extract each amplicon from the DPYD reference FASTA
  2) A short 5' "context anchor" (unique flanking sequence) to locate the variant inside
     the extracted amplicon.

The FASTA is in transcript/coding orientation (5'→3' of mRNA). Uppercase/lowercase reflect
soft-masked repeats (Ensembl convention), NOT exon/intron structure — so we don't use case
to infer CDS positions. Anchors were derived by visual inspection of each amplicon and the
known biology of each variant (splice donor motif for *2A, codon context for the rest).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

from Bio.Seq import Seq

PROJECT = Path(__file__).resolve().parent.parent
FASTA_PATH = PROJECT / "data" / "reference" / "dpyd-001-enst00000370192-ccds30777.fasta"
VARIANTS_PATH = PROJECT / "data" / "variants.json"
AMPLICONS_DIR = PROJECT / "data" / "amplicons"


CONTEXT_ANCHORS = {
    "DPYD_2A":       {"anchor_5": "CAGACAAC",  "anchor_3": "TAAGT"},
    "DPYD_13":       {"anchor_5": "CACTCCTA",  "anchor_3": "TGATC"},
    "DPYD_c2846AT":  {"anchor_5": "ATGATTG",   "anchor_3": "TGAAG"},
    "DPYD_HapB3":    {"anchor_5": "GGACAGA",   "anchor_3": "CAAGA"},
}


def load_fasta(path: Path) -> str:
    return "".join(line.strip() for line in path.read_text().splitlines() if not line.startswith(">"))


def rc(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def extract_amplicon(seq_upper: str, primer_fwd: str, primer_rev: str):
    pf = primer_fwd.upper()
    pr = primer_rev.upper()
    pr_rc = rc(pr)
    fwd_pos = seq_upper.find(pf)
    rev_pos = seq_upper.find(pr_rc)
    if fwd_pos < 0 or rev_pos < 0:
        return None
    start, end = fwd_pos, rev_pos + len(pr)
    if end <= start:
        return None
    return seq_upper[start:end], start, end


def locate_variant(amplicon: str, anchor_5: str, ref_allele: str, anchor_3: str):
    template = (anchor_5 + ref_allele + anchor_3).upper()
    matches = []
    start = 0
    amp = amplicon.upper()
    while True:
        i = amp.find(template, start)
        if i < 0:
            break
        matches.append(i)
        start = i + 1
    if len(matches) != 1:
        return None, len(matches)
    offset = matches[0] + len(anchor_5)
    return offset, 1


def main() -> int:
    seq = load_fasta(FASTA_PATH).upper()
    variants = json.loads(VARIANTS_PATH.read_text())
    AMPLICONS_DIR.mkdir(parents=True, exist_ok=True)

    print(f"FASTA: {len(seq):,} bp (uppercased)")
    print()
    failures = 0

    for vid, v in variants.items():
        anchors = CONTEXT_ANCHORS.get(vid)
        if anchors is None:
            print(f"[SKIP]  {vid}: sin anchors definidos")
            failures += 1
            continue
        if not v.get("primer_fwd") or not v.get("primer_rev"):
            print(f"[SKIP]  {vid}: primers vacíos")
            failures += 1
            continue

        ext = extract_amplicon(seq, v["primer_fwd"], v["primer_rev"])
        if ext is None:
            print(f"[FAIL]  {vid}: primers no encontrados en el FASTA")
            failures += 1
            continue
        amplicon, amp_start, amp_end = ext

        offset, n_matches = locate_variant(
            amplicon,
            anchors["anchor_5"],
            v["ref_allele_sense"],
            anchors["anchor_3"],
        )
        if offset is None:
            print(
                f"[FAIL]  {vid}: anchor_5+ref+anchor_3 "
                f"('{anchors['anchor_5']}{v['ref_allele_sense']}{anchors['anchor_3']}') "
                f"matches={n_matches} (debe ser 1)"
            )
            failures += 1
            continue

        base_at = amplicon[offset].upper()
        if base_at != v["ref_allele_sense"].upper():
            print(f"[FAIL]  {vid}: base en offset {offset} = '{base_at}', esperada '{v['ref_allele_sense']}'")
            failures += 1
            continue

        v["amplicon_sense"] = amplicon
        v["variant_offset_in_amplicon"] = offset
        v["context_anchor_5"] = anchors["anchor_5"]
        v["context_anchor_3"] = anchors["anchor_3"]

        fa_path = AMPLICONS_DIR / f"{vid}.fasta"
        header = (
            f">{vid} {v['name']} | {v['hgvs_coding']} | amp={len(amplicon)}bp | "
            f"variant at offset {offset} (ref={base_at}, alt={v['alt_allele_sense']})"
        )
        fa_path.write_text(f"{header}\n{amplicon}\n")

        context_start = max(0, offset - 10)
        context_end = min(len(amplicon), offset + 11)
        ctx = amplicon[context_start:offset] + f"[{base_at}]" + amplicon[offset + 1:context_end]
        print(
            f"[OK]    {vid}: amp={len(amplicon)}bp, offset={offset}, "
            f"ref={base_at} alt={v['alt_allele_sense']}  context: ...{ctx}..."
        )

    VARIANTS_PATH.write_text(json.dumps(variants, indent=2) + "\n")
    print()
    print(f"{len(variants) - failures}/{len(variants)} variantes OK.")
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
