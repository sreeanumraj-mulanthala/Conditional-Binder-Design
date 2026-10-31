#!/usr/bin/env python3
"""
Rank BindCraft designs for the Anthropic x Adaptyv EGFR challenge.

The competition ranks pH-selectivity first, mouse cross-reactivity second,
affinity third. BindCraft optimises only affinity/confidence, so this script
adds the two objectives BindCraft knows nothing about:

  1. pH switch   - does the binder place a carboxylate (Asp/Glu) close enough
                   to the target's H409 imidazole to form a salt bridge that
                   exists at pH 6.5 and is lost at pH 7.4?
  2. mouse x-react - does the binder footprint avoid the 7 epitope positions
                   that differ between human and mouse EGFR?

plus liability checks (free cysteines, length, composition) before submission.

Usage:
    python3 rank_designs.py <dir-with-pdbs> [final_design_stats.csv]

Target numbering is MATURE EGFR (as in 6ARU and in the trimmed target file).
"""
import sys
import os
import glob
import math
import csv
from collections import defaultdict

import numpy as np

# ---------------------------------------------------------------- constants
HIS_TARGET = 409              # the pH handle
HOTSPOTS = [408, 409, 411, 412]
CONSERVED_PATCH = [382, 408, 409, 411, 412, 415, 438]
DIVERGENT = [353, 418, 443, 467, 468, 471, 473]   # differ in mouse

SALT_BRIDGE_CUT = 4.0    # N...O, a real charged H-bond
NEARBY_CUT = 6.0         # close enough to matter with modest rearrangement
CONTACT_CUT = 4.5        # general interface contact

THREE2ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
    'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
    'TYR': 'Y', 'VAL': 'V',
}

CARBOXYL = {'ASP': ('OD1', 'OD2'), 'GLU': ('OE1', 'OE2')}
HIS_N = ('ND1', 'NE2')

KD = {  # Kyte-Doolittle, for a crude hydrophobicity sanity check
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5,
    'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9,
    'M': 1.9, 'F': 2.8, 'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9,
    'Y': -1.3, 'V': 4.2,
}


# ---------------------------------------------------------------- pdb io
def read_pdb(path):
    atoms = []
    for line in open(path):
        if not line.startswith(('ATOM', 'HETATM')):
            continue
        try:
            atoms.append(dict(
                name=line[12:16].strip(),
                res=line[17:20].strip(),
                ch=line[21],
                seq=int(line[22:26]),
                xyz=(float(line[30:38]), float(line[38:46]),
                     float(line[46:54])),
            ))
        except ValueError:
            continue
    return atoms


# Domain III of the design target, mature numbering 310-480, taken from the
# trimmed target file. Used to recover the numbering offset in design outputs.
TARGET_REF_START = 310


def _load_target_ref():
    """Sequence of the trimmed target, keyed by mature residue number."""
    ref = {}
    for cand in ('EGFR_domainIII.pdb', '../EGFR_domainIII.pdb',
                 '../target/EGFR_domainIII.pdb',
                 os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              'EGFR_domainIII.pdb')):
        if os.path.exists(cand):
            for line in open(cand):
                if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                    ref[int(line[22:26])] = THREE2ONE.get(
                        line[17:20].strip(), 'X')
            break
    return ref


TARGET_REF = _load_target_ref()


def detect_offset(target_atoms):
    """BindCraft renumbers the target chain from 1, discarding the original
    numbering. Recover the shift by matching sequence against the reference.

    Returns `off` such that  mature_number = pdb_number + off.
    """
    if not TARGET_REF:
        return 0
    obs = {}
    for a in target_atoms:
        if a['name'] == 'CA':
            obs[a['seq']] = THREE2ONE.get(a['res'], 'X')
    if not obs:
        return 0
    best, best_score = 0, -1
    lo = min(TARGET_REF) - max(obs)
    hi = max(TARGET_REF) - min(obs)
    for off in range(lo, hi + 1):
        score = sum(1 for k, v in obs.items()
                    if TARGET_REF.get(k + off) == v)
        if score > best_score:
            best, best_score = off, score
    # require a convincing match, else assume numbering is already mature
    return best if best_score >= 0.8 * len(obs) else 0


def split_chains(atoms):
    """Target = the chain matching the reference domain III sequence.
    Binder = the other chain.
    Returns (target_atoms, binder_atoms, (target_id, binder_ids), offset)."""
    bych = defaultdict(list)
    for a in atoms:
        bych[a['ch']].append(a)

    # the target is the chain whose sequence matches the reference
    target_ch, best = None, -1
    for ch, al in bych.items():
        off = detect_offset(al)
        obs = {a['seq']: THREE2ONE.get(a['res'], 'X')
               for a in al if a['name'] == 'CA'}
        score = sum(1 for k, v in obs.items() if TARGET_REF.get(k + off) == v)
        if score > best:
            target_ch, best = ch, score
    if target_ch is None or best <= 0:
        target_ch = max(bych, key=lambda c: len({a['seq'] for a in bych[c]}))

    offset = detect_offset(bych[target_ch])
    binder_chs = [c for c in bych if c != target_ch]
    binder = [a for c in binder_chs for a in bych[c]]
    return bych[target_ch], binder, (target_ch, binder_chs), offset


def sequence_of(atoms):
    seen, seq = {}, []
    for a in atoms:
        if a['name'] == 'CA':
            seen[a['seq']] = THREE2ONE.get(a['res'], 'X')
    for k in sorted(seen):
        seq.append(seen[k])
    return ''.join(seq)


# ---------------------------------------------------------------- metrics
def min_dist(a_xyz, b_xyz):
    if not len(a_xyz) or not len(b_xyz):
        return float('inf')
    A = np.asarray(a_xyz)
    B = np.asarray(b_xyz)
    return float(np.linalg.norm(A[:, None, :] - B[None, :, :], axis=-1).min())


def analyse(path):
    atoms = read_pdb(path)
    if not atoms:
        return None
    target, binder, ids, off = split_chains(atoms)
    if not binder:
        return None

    seq = sequence_of(binder)
    n = len(seq)

    # mature numbering = pdb numbering + off, so a mature target number M
    # is found in this file at (M - off)
    def at(mature):
        return mature - off

    # --- H409 imidazole nitrogens
    his_n = [a['xyz'] for a in target
             if a['seq'] == at(HIS_TARGET) and a['name'] in HIS_N]
    if not his_n:
        print(f"  WARNING: H409 not located in {os.path.basename(path)} "
              f"(offset {off})")

    # --- binder carboxylates
    carbox = defaultdict(list)
    for a in binder:
        if a['res'] in CARBOXYL and a['name'] in CARBOXYL[a['res']]:
            carbox[(a['ch'], a['seq'], a['res'])].append(a['xyz'])

    best_d, best_res = float('inf'), None
    for key, xyzs in carbox.items():
        d = min_dist(xyzs, his_n)
        if d < best_d:
            best_d, best_res = d, key

    # how many distinct acidic residues are within reach
    n_near = sum(1 for key, xyzs in carbox.items()
                 if min_dist(xyzs, his_n) <= NEARBY_CUT)

    # --- interface footprint on the target
    B = np.asarray([a['xyz'] for a in binder])
    T = np.asarray([a['xyz'] for a in target])
    tseq = np.asarray([a['seq'] for a in target])
    contacts = set()
    for k in range(0, len(B), 1500):          # chunk to bound memory
        d = np.linalg.norm(T[:, None, :] - B[None, k:k + 1500, :], axis=-1)
        contacts.update((tseq[(d < CONTACT_CUT).any(axis=1)] + off).tolist())

    hit_hot = sorted(set(contacts) & set(HOTSPOTS))
    hit_cons = sorted(set(contacts) & set(CONSERVED_PATCH))
    hit_div = sorted(set(contacts) & set(DIVERGENT))

    # --- liabilities
    comp = {aa: seq.count(aa) / n for aa in set(seq)} if n else {}
    ncys = seq.count('C')
    hydroph = sum(KD.get(c, 0) for c in seq) / n if n else 0
    net_charge = (seq.count('K') + seq.count('R')
                  - seq.count('D') - seq.count('E'))

    return dict(
        name=os.path.splitext(os.path.basename(path))[0],
        path=path,
        seq=seq,
        length=n,
        target_chain=ids[0],
        binder_chains=','.join(ids[1]),
        numbering_offset=off,
        h409_min_dist=round(best_d, 2) if math.isfinite(best_d) else None,
        h409_partner=(f"{THREE2ONE.get(best_res[2], '?')}{best_res[1]}"
                      if best_res else None),
        n_acidic_near_h409=n_near,
        salt_bridge=bool(math.isfinite(best_d) and best_d <= SALT_BRIDGE_CUT),
        n_contacts=len(contacts),
        hotspots_hit=hit_hot,
        conserved_hit=hit_cons,
        divergent_hit=hit_div,
        n_cys=ncys,
        net_charge=net_charge,
        mean_hydropathy=round(hydroph, 2),
        max_aa_frac=round(max(comp.values()), 2) if comp else None,
    )


def ph_score(r):
    """Higher = better pH switch prospect."""
    d = r['h409_min_dist']
    if d is None or not math.isfinite(d):
        return 0.0
    if d <= SALT_BRIDGE_CUT:
        s = 1.0
    elif d <= NEARBY_CUT:
        s = (NEARBY_CUT - d) / (NEARBY_CUT - SALT_BRIDGE_CUT) * 0.6
    else:
        s = 0.0
    if r['n_acidic_near_h409'] > 1:
        s = min(1.0, s + 0.15 * (r['n_acidic_near_h409'] - 1))
    return round(s, 3)


# Human -> mouse substitutions at the domain III epitope, with a severity
# weight. Counting divergent contacts alone is too crude: a conservative
# substitution that preserves charge and polarity is unlikely to abolish
# binding, whereas one that changes the chemistry of the surface may.
#
#   0.0 = effectively harmless, 1.0 = likely to break the contact
DIVERGENT_PENALTY = {
    353: 0.10,   # R->K  basic to basic, BLOSUM62 +2
    418: 0.45,   # S->G  polar to polar but the side chain is lost entirely
    443: 0.10,   # K->R  basic to basic, +2
    467: 0.15,   # I->M  hydrophobic to hydrophobic, +1
    468: 0.25,   # S->N  polar to polar, +1, but larger
    471: 0.50,   # G->A  polar to hydrophobic
    473: 0.60,   # N->K  polar to basic, a charge appears
}


def mouse_score(r):
    """Higher = safer for mouse cross-reactivity.

    Penalises contacts on divergent positions in proportion to how chemically
    different the mouse residue is, and rewards coverage of the conserved
    patch.
    """
    nc = len(r['conserved_hit'])
    if nc == 0:
        return 0.0
    penalty = sum(DIVERGENT_PENALTY.get(p, 0.5) for p in r['divergent_hit'])
    coverage = min(1.0, nc / 5.0)
    return round(max(0.0, 1.0 - 0.4 * penalty) * coverage, 3)


def flags(r):
    out = []
    if not (10 <= r['length'] <= 250):
        out.append('LENGTH-INVALID')
    if set(r['seq']) - set(THREE2ONE.values()):
        out.append('NON-STANDARD-AA')
    if r['n_cys'] == 1:
        out.append('single-Cys')
    elif r['n_cys'] > 2:
        out.append(f"{r['n_cys']}-Cys")
    if r['max_aa_frac'] and r['max_aa_frac'] > 0.30:
        out.append('low-complexity')
    if r['mean_hydropathy'] > 0.5:
        out.append('hydrophobic')
    if abs(r['net_charge']) > 12:
        out.append(f"charge{r['net_charge']:+d}")
    if not r['hotspots_hit']:
        out.append('MISSES-HOTSPOTS')
    return out


# ---------------------------------------------------------------- main
def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    d = sys.argv[1]
    stats_csv = sys.argv[2] if len(sys.argv) > 2 else None

    pdbs = sorted(glob.glob(os.path.join(d, '**', '*.pdb'), recursive=True))
    if not pdbs:
        print(f'no .pdb files under {d}')
        sys.exit(1)

    # optional BindCraft metrics, keyed loosely by design name
    bc = {}
    if stats_csv and os.path.exists(stats_csv):
        with open(stats_csv) as fh:
            for row in csv.DictReader(fh):
                key = (row.get('Design') or row.get('design')
                       or row.get('Name') or row.get('name') or '')
                if key:
                    bc[key.strip()] = row

    rows = []
    for p in pdbs:
        r = analyse(p)
        if r is None:
            print(f'  skipped (unparsable or single chain): {p}')
            continue
        r['ph_score'] = ph_score(r)
        r['mouse_score'] = mouse_score(r)
        r['flags'] = flags(r)

        m = bc.get(r['name'], {})
        for src, dst in (('i_pTM', 'i_ptm'), ('i_ptm', 'i_ptm'),
                         ('pLDDT', 'plddt'), ('plddt', 'plddt'),
                         ('i_pAE', 'i_pae'), ('i_pae', 'i_pae'),
                         ('ShapeComplementarity', 'sc'),
                         ('Binder_Energy_Score', 'energy'),
                         ('dSASA', 'dsasa'), ('n_InterfaceHbonds', 'hbonds')):
            if src in m and m[src] not in ('', None):
                try:
                    r[dst] = float(m[src])
                except ValueError:
                    pass

        # composite: pH-selectivity first, then mouse, then confidence
        conf = r.get('i_ptm', 0.0) or 0.0
        r['total'] = round(2.0 * r['ph_score'] + 1.2 * r['mouse_score']
                           + 0.8 * conf, 3)
        rows.append(r)

    rows.sort(key=lambda x: -x['total'])

    print(f'\n{len(rows)} designs analysed\n')
    hdr = (f"{'rank':>4} {'name':<28} {'len':>4} {'H409':>6} {'partner':>8} "
           f"{'pH':>5} {'mouse':>6} {'i_ptm':>6} {'div':>4} {'total':>6}  flags")
    print(hdr)
    print('-' * len(hdr))
    for i, r in enumerate(rows, 1):
        print(f"{i:>4} {r['name'][:28]:<28} {r['length']:>4} "
              f"{(r['h409_min_dist'] if r['h409_min_dist'] is not None else -1):>6.2f} "
              f"{str(r['h409_partner'] or '-'):>8} "
              f"{r['ph_score']:>5.2f} {r['mouse_score']:>6.2f} "
              f"{r.get('i_ptm', float('nan')):>6.2f} "
              f"{len(r['divergent_hit']):>4} {r['total']:>6.2f}  "
              f"{','.join(r['flags'])}")

    out = os.path.join(d, 'ranked_designs.csv')
    cols = ['name', 'length', 'h409_min_dist', 'h409_partner',
            'n_acidic_near_h409', 'salt_bridge', 'ph_score', 'mouse_score',
            'i_ptm', 'plddt', 'i_pae', 'n_contacts', 'hotspots_hit',
            'conserved_hit', 'divergent_hit', 'n_cys', 'net_charge',
            'mean_hydropathy', 'total', 'flags', 'seq']
    with open(out, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for r in rows:
            w.writerow([r.get(c, '') for c in cols])
    print(f'\nwritten: {out}')


if __name__ == '__main__':
    main()
