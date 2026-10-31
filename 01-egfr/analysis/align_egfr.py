
# --- path resolution: files live in <repo>/target regardless of where this is run from
import os as _os
_T = _os.path.join(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))), 'target')
def _t(n): return _os.path.join(_T, n)
# ---

#!/usr/bin/env python3
"""
Align human (P00533) and mouse (Q01279) EGFR and report conservation
across the extracellular region, with emphasis on domain III.

Needleman-Wunsch global alignment, BLOSUM62-like simple scoring.
"""

def read_fasta(path):
    seq = []
    for line in open(path):
        line = line.strip()
        if not line or line.startswith('>'):
            continue
        seq.append(line)
    return ''.join(seq)


def nw_align(a, b, match=2, mismatch=-1, gap=-6):
    n, m = len(a), len(b)
    # score matrix rows, keep traceback compactly
    prev = [j * gap for j in range(m + 1)]
    tb = [[0] * (m + 1) for _ in range(n + 1)]
    for j in range(1, m + 1):
        tb[0][j] = 1  # left
    for i in range(1, n + 1):
        cur = [i * gap] + [0] * m
        tb[i][0] = 2  # up
        ai = a[i - 1]
        for j in range(1, m + 1):
            s = match if ai == b[j - 1] else mismatch
            diag = prev[j - 1] + s
            left = cur[j - 1] + gap
            up = prev[j] + gap
            best = diag
            ptr = 0
            if left > best:
                best, ptr = left, 1
            if up > best:
                best, ptr = up, 2
            cur[j] = best
            tb[i][j] = ptr
        prev = cur
    # traceback
    i, j = n, m
    A, B = [], []
    while i > 0 or j > 0:
        p = tb[i][j]
        if i > 0 and j > 0 and p == 0:
            A.append(a[i - 1]); B.append(b[j - 1]); i -= 1; j -= 1
        elif j > 0 and p == 1:
            A.append('-'); B.append(b[j - 1]); j -= 1
        else:
            A.append(a[i - 1]); B.append('-'); i -= 1
    return ''.join(reversed(A)), ''.join(reversed(B))


human = read_fasta(_t('human_EGFR_P00533.fasta'))
mouse = read_fasta(_t('mouse_EGFR_Q01279.fasta'))

# Restrict to the mature extracellular region to keep the alignment clean.
# UniProt numbering: signal peptide 1-24; competition target is residues 25-645.
H_START, H_END = 25, 645
h_ecd = human[H_START - 1:H_END]
m_ecd = mouse[H_START - 1:H_END]

ah, am = nw_align(h_ecd, m_ecd)

# Map alignment columns back to human UniProt numbering
rows = []
hi = H_START - 1
mi = H_START - 1
for ca, cb in zip(ah, am):
    if ca != '-':
        hi += 1
    if cb != '-':
        mi += 1
    rows.append((hi if ca != '-' else None, ca, cb,
                 mi if cb != '-' else None))

ident = sum(1 for _, x, y, _ in rows if x == y and x != '-')
aligned = sum(1 for _, x, y, _ in rows if x != '-' and y != '-')
print(f"Extracellular region {H_START}-{H_END}")
print(f"  aligned columns: {aligned}")
print(f"  identical:       {ident}  ({100*ident/aligned:.1f}% identity)")
print()

# Domain boundaries. The EGFR literature numbers the MATURE protein
# (signal peptide 1-24 removed), so mature position X = UniProt precursor X+24.
# Mature: I 1-165, II 166-310, III 311-480, IV 481-620.
OFFSET = 24
DOMAINS = {
    'I  (L1)':   (1 + OFFSET, 165 + OFFSET),
    'II (CR1)':  (166 + OFFSET, 310 + OFFSET),
    'III (L2)':  (311 + OFFSET, 480 + OFFSET),
    'IV (CR2)':  (481 + OFFSET, 620 + OFFSET),
}
print("Per-domain human/mouse identity:")
for name, (s, e) in DOMAINS.items():
    sub = [(h, x, y) for h, x, y, _ in rows if h is not None and s <= h <= e]
    idn = sum(1 for _, x, y in sub if x == y)
    print(f"  Domain {name}  {s}-{e}: {idn}/{len(sub)} = {100*idn/len(sub):.1f}%")
print()

# Differences within domain III
print("Human/mouse differences in domain III (mature 311-480 = precursor 335-504):")
diffs = [(h, x, y) for h, x, y, _ in rows
         if h is not None and 335 <= h <= 504 and x != y]
line = []
for h, x, y in diffs:
    line.append(f"{x}{h}{y}")
for k in range(0, len(line), 8):
    print("   " + "  ".join(line[k:k + 8]))
print(f"   total: {len(diffs)} differences")
print()

# Native histidines in the extracellular region (relevant to pH behaviour)
his_h = [h for h, x, y, _ in rows if h is not None and x == 'H']
print("Histidines in human EGFR extracellular region:", his_h)
his3 = [h for h in his_h if 310 <= h <= 480]
print("  of which in domain III:", his3)

with open(_t('human_mouse_ecd_alignment.txt'), 'w') as fh:
    for k in range(0, len(ah), 60):
        fh.write(f"human {ah[k:k+60]}\n")
        fh.write(f"mouse {am[k:k+60]}\n")
        fh.write("      " + ''.join('*' if a == b else ' '
                 for a, b in zip(ah[k:k+60], am[k:k+60])) + "\n\n")
print("\nalignment written to human_mouse_ecd_alignment.txt")
