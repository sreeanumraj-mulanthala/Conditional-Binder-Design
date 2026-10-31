
# --- path resolution: files live in <repo>/target regardless of where this is run from
import os as _os
_T = _os.path.join(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))), 'target')
def _t(n): return _os.path.join(_T, n)
# ---

#!/usr/bin/env python3
"""
Write BindCraft-ready target PDBs from 6ARU.

Outputs (mature EGFR numbering, as 6ARU uses):
  EGFR_chainA_full.pdb    - chain A, whole modelled ectodomain
  EGFR_domainIII.pdb      - chain A residues 310-480 only (domain III)

The Fab chains B and C are discarded: the binder must be designed against
the free receptor, not against a surface already occupied by an antibody.
"""

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

def pdb_line(serial, name, altloc, resn, chain, resi, icode, x, y, z,
             occ, bfac, elem):
    # PDB atom-name column rule: 1-3 char names start in column 14
    nm = name if len(name) == 4 else ' ' + name
    return (f"ATOM  {serial:>5} {nm:<4}{altloc:1}{resn:>3} {chain:1}"
            f"{resi:>4}{icode:1}   {x:8.3f}{y:8.3f}{z:8.3f}"
            f"{occ:6.2f}{bfac:6.2f}          {elem:>2}\n")


def write_pdb(out, lo=None, hi=None):
    serial = 0
    seen = set()
    with open(out, 'w') as fh:
        for r in rows:
            if r[idx['group_PDB']] != 'ATOM':
                continue
            if r[idx['auth_asym_id']] != 'A':
                continue
            try:
                resi = int(r[idx['auth_seq_id']])
            except ValueError:
                continue
            if lo is not None and not (lo <= resi <= hi):
                continue
            alt = r[idx['label_alt_id']] if 'label_alt_id' in idx else '.'
            if alt not in ('.', '?', 'A'):
                continue           # keep only the primary conformer
            serial += 1
            seen.add(resi)
            fh.write(pdb_line(
                serial,
                r[idx['label_atom_id']].strip('"'),
                ' ',
                r[idx['label_comp_id']],
                'A',
                resi,
                ' ',
                float(r[idx['Cartn_x']]),
                float(r[idx['Cartn_y']]),
                float(r[idx['Cartn_z']]),
                1.00,
                float(r[idx['B_iso_or_equiv']]) if 'B_iso_or_equiv' in idx else 0.0,
                r[idx['type_symbol']]))
        fh.write('TER\nEND\n')
    return serial, sorted(seen)


n, res = write_pdb(_t('EGFR_chainA_full.pdb'))
print(f'EGFR_chainA_full.pdb  : {n} atoms, {len(res)} residues, '
      f'{res[0]}-{res[-1]}')

n, res = write_pdb(_t('EGFR_domainIII.pdb'), 310, 480)
print(f'EGFR_domainIII.pdb    : {n} atoms, {len(res)} residues, '
      f'{res[0]}-{res[-1]}')

HOT = [409, 411, 412, 408]
print(f'\nhotspots present in domain III file: '
      f'{[h for h in HOT if h in res]}')
print('hotspot string for BindCraft: A' + ',A'.join(str(h) for h in HOT))
