#!/usr/bin/env python3
"""
Epitope analysis on 6ARU chain A (EGFR).

NOTE ON NUMBERING: 6ARU auth_seq_id uses MATURE EGFR numbering
(verified: 607/609 residues match human P00533 at offset +24, i.e.
precursor = mature + 24). Literature epitope numbers can therefore be
used directly against this file.

Outputs:
  - empirical Fab epitope (contacts <4.5 A)
  - human/mouse conservation for each epitope residue
  - per-residue SASA on the unbound receptor (Shrake-Rupley)
  - burial and electrostatic environment of H409
  - contiguity of the proposed conserved hotspot patch
"""
import numpy as np
from collections import defaultdict

# --- path resolution: files live in <repo>/target regardless of where this is run from
import os as _os
_T = _os.path.join(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))), 'target')
def _t(n): return _os.path.join(_T, n)
# ---


T = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E',
     'GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F',
     'PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}

VDW = {'C':1.70,'N':1.55,'O':1.52,'S':1.80,'SE':1.90,'P':1.80}
PROBE = 1.4
OFFSET = 24  # precursor = mature + OFFSET

# Gly-X-Gly reference SASA (A^2), Tien et al. 2013 theoretical
REF = {'A':129,'R':274,'N':195,'D':193,'C':167,'Q':225,'E':223,'G':104,
       'H':224,'I':197,'L':201,'K':236,'M':224,'F':240,'P':159,'S':155,
       'T':172,'W':285,'Y':263,'V':174}


def parse_atom_site(path):
    cols, rows, in_loop, reading = [], [], False, False
    for line in open(path):
        s = line.strip()
        if s.startswith('loop_'):
            in_loop, cols, reading = True, [], False
            continue
        if in_loop and s.startswith('_atom_site.'):
            cols.append(s.split('.', 1)[1]); reading = True; continue
        if reading:
            if s.startswith('#') or s.startswith('loop_') or s.startswith('_'):
                break
            if s:
                p = s.split()
                if len(p) == len(cols):
                    rows.append(p)
    return {c: i for i, c in enumerate(cols)}, rows


idx, rows = parse_atom_site(_t('6aru.cif'))
atoms = []
for r in rows:
    if r[idx['group_PDB']] != 'ATOM':
        continue
    try:
        seq = int(r[idx['auth_seq_id']])
        xyz = [float(r[idx['Cartn_x']]), float(r[idx['Cartn_y']]),
               float(r[idx['Cartn_z']])]
    except ValueError:
        continue
    atoms.append(dict(name=r[idx['label_atom_id']], res=r[idx['label_comp_id']],
                      ch=r[idx['auth_asym_id']], seq=seq, xyz=xyz,
                      el=r[idx['type_symbol']].upper()))

A = [a for a in atoms if a['ch'] == 'A']
FAB = [a for a in atoms if a['ch'] in ('B', 'C')]
resname = {}
for a in A:
    resname[a['seq']] = T.get(a['res'], 'X')

# ------------------------------------------------ empirical epitope
ac = np.array([a['xyz'] for a in A]); aseq = np.array([a['seq'] for a in A])
fc = np.array([a['xyz'] for a in FAB])
epi = set()
for k in range(0, len(fc), 1500):
    d = np.linalg.norm(ac[:, None, :] - fc[None, k:k+1500, :], axis=-1)
    epi.update(aseq[(d < 4.5).any(axis=1)].tolist())
epitope = sorted(epi)

# ------------------------------------------------ conservation
def rd(p):
    return ''.join(l.strip() for l in open(p) if l.strip() and not l.startswith('>'))
hum, mou = rd(_t('human_EGFR_P00533.fasta')), rd(_t('mouse_EGFR_Q01279.fasta'))

# ------------------------------------------------ SASA (Shrake-Rupley)
def sphere_points(n=192):
    pts, inc = [], np.pi * (3 - np.sqrt(5))
    off = 2.0 / n
    for i in range(n):
        y = i * off - 1 + off / 2
        r = np.sqrt(max(0.0, 1 - y * y))
        phi = i * inc
        pts.append([np.cos(phi) * r, y, np.sin(phi) * r])
    return np.array(pts)

SP = sphere_points(192)

def sasa(atom_list):
    xyz = np.array([a['xyz'] for a in atom_list])
    rad = np.array([VDW.get(a['el'], 1.70) for a in atom_list]) + PROBE
    out = np.zeros(len(atom_list))
    # neighbour search via grid-free pairwise in chunks
    for i in range(len(atom_list)):
        d = np.linalg.norm(xyz - xyz[i], axis=1)
        nb = np.where((d < rad + rad[i]) & (d > 0))[0]
        pts = xyz[i] + SP * rad[i]
        if len(nb):
            dd = np.linalg.norm(pts[:, None, :] - xyz[nb][None, :, :], axis=-1)
            acc = (dd >= rad[nb][None, :]).all(axis=1)
        else:
            acc = np.ones(len(pts), bool)
        out[i] = 4 * np.pi * rad[i] ** 2 * acc.mean()
    return out

print('computing SASA on unbound chain A ...')
sa = sasa(A)
res_sasa = defaultdict(float)
for a, v in zip(A, sa):
    res_sasa[a['seq']] += v

print('\n=== Empirical Fab epitope on EGFR (6ARU, mature numbering) ===')
print(f'{len(epitope)} residues: {epitope}\n')
print(f"{'mature':>6} {'prec':>5} {'aa':>3} {'mouse':>5} {'cons':>5} "
      f"{'SASA':>7} {'rel%':>6}")
cons_patch, div_patch = [], []
for m in epitope:
    p = m + OFFSET
    if p > len(hum):
        continue
    h, mo = hum[p-1], mou[p-1]
    aa = resname.get(m, '?')
    rel = 100 * res_sasa[m] / REF.get(aa, 200)
    same = (h == mo)
    (cons_patch if same else div_patch).append(m)
    print(f'{m:>6} {p:>5} {aa:>3} {mo:>5} {"yes" if same else "NO":>5} '
          f'{res_sasa[m]:7.1f} {rel:6.1f}')

print(f'\nconserved epitope positions: {cons_patch}')
print(f'divergent epitope positions: {div_patch}')

# ------------------------------------------------ H409 environment
print('\n=== H409 environment (the proposed pH handle) ===')
h409 = [a for a in A if a['seq'] == 409]
if h409:
    print(f'  residue: {h409[0]["res"]}409   SASA {res_sasa[409]:.1f} A^2 '
          f'({100*res_sasa[409]/REF["H"]:.1f} % of Gly-H-Gly reference)')
    side = [a for a in h409 if a['name'] in ('ND1', 'NE2', 'CE1', 'CD2', 'CG')]
    sc = np.array([a['xyz'] for a in side])
    # charged groups nearby
    CH = {'ASP': ('OD1', 'OD2'), 'GLU': ('OE1', 'OE2'),
          'LYS': ('NZ',), 'ARG': ('NH1', 'NH2', 'NE')}
    near = []
    for a in A:
        if a['res'] in CH and a['name'] in CH[a['res']] and a['seq'] != 409:
            d = np.linalg.norm(sc - np.array(a['xyz']), axis=1).min()
            if d < 10:
                near.append((round(float(d), 2), f"{T[a['res']]}{a['seq']}", a['name']))
    near.sort()
    print('  charged groups within 10 A of the imidazole:')
    seen = set()
    for d, lab, nm in near:
        if lab in seen:
            continue
        seen.add(lab)
        sign = '-' if lab[0] in 'DE' else '+'
        print(f'    {lab:>7} ({sign})  {d:5.2f} A')

# ------------------------------------------------ patch contiguity
print('\n=== Contiguity of the conserved epitope patch ===')
def cb(seq):
    cands = [a for a in A if a['seq'] == seq and a['name'] in ('CB', 'CA')]
    if not cands:
        return None
    pref = [a for a in cands if a['name'] == 'CB'] or cands
    return np.array(pref[0]['xyz'])

pts = {m: cb(m) for m in cons_patch}
pts = {k: v for k, v in pts.items() if v is not None}
ks = sorted(pts)
print('    pairwise CB distances (A):')
hdr = '        ' + ''.join(f'{k:>7}' for k in ks)
print(hdr)
for i in ks:
    row = f'    {i:>4}'
    for j in ks:
        row += f'{np.linalg.norm(pts[i]-pts[j]):7.1f}'
    print(row)
