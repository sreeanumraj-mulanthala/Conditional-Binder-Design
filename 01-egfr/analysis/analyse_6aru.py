#!/usr/bin/env python3
"""
Analyse PDB 6ARU (EGFR extracellular region + cetuximab Fab mutant).

Goals:
  1. Identify chains and what they are.
  2. Derive the Fab epitope on EGFR empirically from the coordinates
     (rather than trusting the literature list).
  3. Cross-reference that epitope against human/mouse conservation.
  4. Measure burial of H433 (precursor numbering), the proposed pH handle.
  5. Check the proposed hotspot set forms a contiguous, exposed patch.
"""
import math
from collections import defaultdict
import numpy as np

# --- path resolution: files live in <repo>/target regardless of where this is run from
import os as _os
_T = _os.path.join(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))), 'target')
def _t(n): return _os.path.join(_T, n)
# ---


CIF = _t('6aru.cif')

# ---------------------------------------------------------------- parse
def parse_atom_site(path):
    """Minimal mmCIF atom_site loop parser."""
    cols, rows = [], []
    in_loop = False
    reading = False
    for line in open(path):
        s = line.strip()
        if s.startswith('loop_'):
            in_loop, cols, reading = True, [], False
            continue
        if in_loop and s.startswith('_atom_site.'):
            cols.append(s.split('.', 1)[1])
            reading = True
            continue
        if reading:
            if s.startswith('#') or s.startswith('loop_') or s.startswith('_'):
                break
            if not s:
                continue
            parts = s.split()
            if len(parts) == len(cols):
                rows.append(parts)
        elif in_loop and cols and not s.startswith('_atom_site.'):
            in_loop = False
    idx = {c: i for i, c in enumerate(cols)}
    return idx, rows


idx, rows = parse_atom_site(CIF)
need = ['group_PDB', 'label_atom_id', 'label_comp_id', 'auth_asym_id',
        'auth_seq_id', 'Cartn_x', 'Cartn_y', 'Cartn_z', 'type_symbol']
missing = [c for c in need if c not in idx]
if missing:
    raise SystemExit(f'missing columns: {missing}')

atoms = []
for r in rows:
    if r[idx['group_PDB']] != 'ATOM':
        continue
    try:
        xyz = (float(r[idx['Cartn_x']]), float(r[idx['Cartn_y']]),
               float(r[idx['Cartn_z']]))
        seq = int(r[idx['auth_seq_id']])
    except ValueError:
        continue
    atoms.append(dict(
        name=r[idx['label_atom_id']], res=r[idx['label_comp_id']],
        ch=r[idx['auth_asym_id']], seq=seq, xyz=xyz,
        el=r[idx['type_symbol']]))

print(f'parsed {len(atoms)} protein atoms')

# ---------------------------------------------------------------- chains
chains = defaultdict(set)
for a in atoms:
    chains[a['ch']].add((a['seq'], a['res']))
print('\nchains:')
for ch in sorted(chains):
    res = sorted(chains[ch])
    print(f'  {ch}: {len(res):4d} residues, auth_seq {res[0][0]}–{res[-1][0]}')

THREE2ONE = {
 'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E',
 'GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F',
 'PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}

# EGFR chain = the long one
egfr_ch = max(chains, key=lambda c: len(chains[c]))
other = [c for c in chains if c != egfr_ch]
print(f'\nEGFR chain (longest) = {egfr_ch}; other chains = {other}')

# ------------------------------------------------- empirical epitope
coords = np.array([a['xyz'] for a in atoms])
chs = np.array([a['ch'] for a in atoms])
seqs = np.array([a['seq'] for a in atoms])

egfr_mask = chs == egfr_ch
CUT = 4.5

for partner in other:
    pm = chs == partner
    if pm.sum() == 0:
        continue
    A = coords[egfr_mask]
    B = coords[pm]
    # chunked distance
    contact_res = set()
    sa = seqs[egfr_mask]
    for k in range(0, len(B), 2000):
        d = np.linalg.norm(A[:, None, :] - B[None, k:k+2000, :], axis=-1)
        hit = (d < CUT).any(axis=1)
        contact_res.update(sa[hit].tolist())
    if contact_res:
        print(f'\nEGFR residues within {CUT} A of chain {partner} '
              f'({len(contact_res)}):')
        print('  ', sorted(contact_res))
