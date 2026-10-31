#!/usr/bin/env python3
"""
Build and validate the Proteinbase submission CSV.

Competition rules checked here:
  - columns: name, sequence, molecule_class
  - up to 20 designs
  - each chain 10-250 aa
  - standard amino acids only
  - no duplicate sequences, no duplicate names
  - single_chain molecule class (de novo mini-binders)
"""
import csv
import os
import sys
import glob
import importlib.util

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    'rd', os.path.join(HERE, 'rank_designs.py'))
rd = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rd)

STANDARD = set('ACDEFGHIKLMNPQRSTVWY')
MIN_LEN, MAX_LEN, MAX_DESIGNS = 10, 250, 20
PREFIX = 'AMR_EGFR'          # initials + target, so names are identifiable


def collect(pdb_dir):
    rows = []
    for p in sorted(glob.glob(os.path.join(pdb_dir, '*.pdb'))):
        r = rd.analyse(p)
        if not r:
            continue
        r['ph_score'] = rd.ph_score(r)
        r['mouse_score'] = rd.mouse_score(r)
        r['flags'] = rd.flags(r)
        rows.append(r)
    # one entry per unique sequence, keeping the best H409 geometry
    uniq = {}
    for r in rows:
        k = r['seq']
        if k not in uniq or r['h409_min_dist'] < uniq[k]['h409_min_dist']:
            uniq[k] = r
    return sorted(uniq.values(),
                  key=lambda x: (-(2 * x['ph_score'] + 1.2 * x['mouse_score']),
                                 x['h409_min_dist']))


def validate(entries):
    """Hard checks. Anything here would be rejected at upload."""
    problems = []
    seen_seq, seen_name = set(), set()
    for name, seq in entries:
        if not (MIN_LEN <= len(seq) <= MAX_LEN):
            problems.append(f'{name}: length {len(seq)} outside {MIN_LEN}-{MAX_LEN}')
        bad = set(seq) - STANDARD
        if bad:
            problems.append(f'{name}: non-standard residues {sorted(bad)}')
        if seq in seen_seq:
            problems.append(f'{name}: duplicate sequence')
        if name in seen_name:
            problems.append(f'{name}: duplicate name')
        if not name.replace('_', '').replace('-', '').isalnum():
            problems.append(f'{name}: unusual characters in name')
        seen_seq.add(seq)
        seen_name.add(name)
    if len(entries) > MAX_DESIGNS:
        problems.append(f'{len(entries)} designs, maximum is {MAX_DESIGNS}')
    return problems


def main():
    pdb_dir = sys.argv[1] if len(sys.argv) > 1 else 'run1/all'
    out = sys.argv[2] if len(sys.argv) > 2 else 'submission_01.csv'

    sel = collect(pdb_dir)[:MAX_DESIGNS]
    entries = [(f'{PREFIX}_{i:02d}', r['seq']) for i, r in enumerate(sel, 1)]

    problems = validate(entries)
    print(f'{len(entries)} designs prepared\n')
    print(f"{'name':<14}{'len':>4}{'H409 A':>8}{'via':>6}{'net':>5}"
          f"{'Cys':>4}  sequence")
    print('-' * 100)
    for (name, seq), r in zip(entries, sel):
        print(f'{name:<14}{len(seq):>4}{r["h409_min_dist"]:>8.2f}'
              f'{str(r["h409_partner"]):>6}{r["net_charge"]:>5}{r["n_cys"]:>4}  '
              f'{seq[:46]}{"..." if len(seq) > 46 else ""}')

    print()
    if problems:
        print('VALIDATION FAILED:')
        for p in problems:
            print('   ', p)
        sys.exit(1)
    print('validation passed:')
    print(f'   {len(entries)} designs, all {MIN_LEN}-{MAX_LEN} aa')
    print('   standard amino acids only')
    print('   no duplicate sequences or names')

    with open(out, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['name', 'sequence', 'molecule_class'])
        for name, seq in entries:
            w.writerow([name, seq, 'single_chain'])
    print(f'\nwritten: {out}')

    # a companion table recording why each design was chosen
    meta = os.path.splitext(out)[0] + '_rationale.csv'
    cols = ['name', 'source', 'length', 'h409_min_dist', 'h409_partner',
            'n_acidic_near_h409', 'ph_score', 'mouse_score', 'hotspots_hit',
            'conserved_hit', 'divergent_hit', 'n_contacts', 'net_charge',
            'n_cys', 'mean_hydropathy', 'sequence']
    with open(meta, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for (name, seq), r in zip(entries, sel):
            w.writerow([name, os.path.basename(r['path']).replace('.pdb', ''),
                        r['length'], r['h409_min_dist'], r['h409_partner'],
                        r['n_acidic_near_h409'], r['ph_score'],
                        r['mouse_score'], ' '.join(map(str, r['hotspots_hit'])),
                        ' '.join(map(str, r['conserved_hit'])),
                        ' '.join(map(str, r['divergent_hit'])),
                        r['n_contacts'], r['net_charge'], r['n_cys'],
                        r['mean_hydropathy'], seq])
    print(f'written: {meta}')


if __name__ == '__main__':
    main()
