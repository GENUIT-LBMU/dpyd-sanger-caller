"""Diagnostic: extract amplicons with primers and print their content + context around
positions that match ref_allele_sense. Use the output to manually locate variants.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

from Bio.Seq import Seq

PROJECT = Path(__file__).resolve().parent.parent
FASTA_PATH = PROJECT / "data" / "reference" / "dpyd-001-enst00000370192-ccds30777.fasta"
VARIANTS_PATH = PROJECT / "data" / "variants.json"


def load_fasta(path: Path) -> str:
    return "".join(line.strip() for line in path.read_text().splitlines() if not line.startswith(">"))


def find_all(seq: str, sub: str) -> list[int]:
    positions, start = [], 0
    while True:
        i = seq.find(sub, start)
        if i < 0:
            return positions
        positions.append(i)
        start = i + 1


def main() -> int:
    seq = load_fasta(FASTA_PATH).upper()
    variants = json.loads(VARIANTS_PATH.read_text())

    for vid, v in variants.items():
        print("=" * 80)
        print(f"{vid} — {v['name']} — {v['hgvs_coding']} — ref={v['ref_allele_sense']} alt={v['alt_allele_sense']}")
        pf = v["primer_fwd"].upper()
        pr = v["primer_rev"].upper()
        pr_rc = str(Seq(pr).reverse_complement())

        fwd_positions = find_all(seq, pf)
        rev_positions = find_all(seq, pr_rc)
        print(f"primer_fwd matches in FASTA: {len(fwd_positions)} @ {fwd_positions[:5]}")
        print(f"primer_rev(RC) matches in FASTA: {len(rev_positions)} @ {rev_positions[:5]}")

        if not fwd_positions or not rev_positions:
            print("NO se puede extraer amplicón — revisar primers.")
            continue

        # Pick the first pair that gives a positive amplicon length
        best = None
        for f in fwd_positions:
            for r in rev_positions:
                length = (r + len(pr)) - f
                if 50 <= length <= 2000:
                    best = (f, r, length)
                    break
            if best:
                break
        if not best:
            print("Ningún par fwd/rev da un amplicón de tamaño razonable.")
            continue

        f, r, length = best
        amplicon = seq[f : r + len(pr)]
        print(f"\nAmplicón: {length} bp, FASTA [{f}..{r + len(pr)}]\n")
        # Print with 60 char wrap
        for i in range(0, len(amplicon), 60):
            print(f"  {i:4d} {amplicon[i:i+60]}")

        # Candidate positions of ref_allele
        ref = v["ref_allele_sense"].upper()
        ref_positions = [i for i, b in enumerate(amplicon) if b == ref]
        print(f"\nPosiciones de '{ref}' en el amplicón: {len(ref_positions)} (primeras 20: {ref_positions[:20]})")
        print()

    return 0


if __name__ == "__main__":
    sys.exit(main())
