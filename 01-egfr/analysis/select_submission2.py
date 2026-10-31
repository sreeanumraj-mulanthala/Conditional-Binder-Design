#!/usr/bin/env python3
"""
Build submission 2: maximise distinct binding modes carrying the H409 salt
bridge, with a small number of mechanism-negative controls.

Selection rules
  - a design qualifies if a binder carboxylate sits within 4.0 A of the H409
    imidazole nitrogens
  - it must also have been fully scored by BindCraft (so we know it is not
    structurally poor), with i_pTM >= 0.65 and pLDDT >= 0.80
  - at most 2 per binding mode, and only where the two differ enough in
    sequence to count as distinct (<= 85% identity)
  - modes are ordered by the merit of their best design; stronger modes get
    the second slot
  - the last slots go to controls: well-scored designs with NO carboxylate
    near H409, to give an internal comparison
"""
import csv
import os
from collections import defaultdict

SALT = 4.0
CONTROL_MIN = 7.0        # a control must be clearly non-interacting
MIN_IPTM, MIN_PLDDT = 0.65, 0.80
MAX_TOTAL = 20
N_CONTROLS = 2
MAX_IDENTITY = 0.85


def load():
    stats = {r['Design']: r for r in
             csv.DictReader(open('run3/mpnn_design_stats.csv')) if r.get('Design')}
    geom = list(csv.DictReader(open('run3/flat/ranked_designs.csv')))
    out = {}
    for g in geom:
        if g['h409_min_dist'] in ('', 'None'):
            continue
        d = g['name'].rsplit('_model', 1)[0]
        s = stats.get(d)
        if not s:
            continue
        dist = float(g['h409_min_dist'])
        # keep the model with the best geometry for each design
        if d in out and out[d]['dist'] <= dist:
            continue
        out[d] = dict(
            design=d, traj='_'.join(d.split('_')[:3]), dist=dist,
            partner=g['h409_partner'], seq=g['seq'], length=int(g['length']),
            mouse=float(g['mouse_score']),
            divergent=g['divergent_hit'], conserved=g['conserved_hit'],
            hotspots=g['hotspots_hit'], n_contacts=int(g['n_contacts']),
            net=int(g['net_charge']), cys=int(g['n_cys']),
            iptm=float(s['Average_i_pTM']), plddt=float(s['Average_pLDDT']),
            phob=float(s['Average_Surface_Hydrophobicity']),
            sc=float(s['Average_ShapeComplementarity']),
            dG=float(s['Average_dG']),
        )
    return list(out.values())


def identity(a, b):
    n = min(len(a), len(b))
    return sum(1 for i in range(n) if a[i] == b[i]) / n


def merit(d):
    """Rank within a mode: close salt bridge, conserved footprint, confident."""
    prox = max(0.0, (SALT - d['dist']) / SALT)     # 0..1, closer is better
    return 2.0 * prox + 1.2 * d['mouse'] + 1.0 * d['iptm'] + 0.5 * d['plddt']


def main():
    alld = load()
    qual = [d for d in alld if d['iptm'] >= MIN_IPTM and d['plddt'] >= MIN_PLDDT]

    sb = [d for d in qual if d['dist'] <= SALT]
    ctl = [d for d in qual if d['dist'] >= CONTROL_MIN]

    modes = defaultdict(list)
    for d in sb:
        modes[d['traj']].append(d)
    for v in modes.values():
        v.sort(key=merit, reverse=True)

    order = sorted(modes, key=lambda t: merit(modes[t][0]), reverse=True)
    print(f"{len(sb)} qualifying salt-bridge designs across {len(order)} modes")
    print(f"{len(ctl)} qualifying control designs (>= {CONTROL_MIN} A)\n")

    picked = []
    # one per mode, best first
    for t in order:
        picked.append(modes[t][0])
    # second from the strongest modes, if sequence-distinct
    slots = MAX_TOTAL - N_CONTROLS - len(picked)
    for t in order:
        if slots <= 0:
            break
        for cand in modes[t][1:]:
            if identity(cand['seq'], modes[t][0]['seq']) <= MAX_IDENTITY:
                picked.append(cand)
                slots -= 1
                break
    picked = picked[:MAX_TOTAL - N_CONTROLS]

    # controls: best-scoring, each from a different trajectory, and from modes
    # not already represented, so they are independent of the test set
    used = {d['traj'] for d in picked}
    ctl.sort(key=lambda d: -(d['iptm'] + 0.5 * d['plddt']))
    controls, ctl_trajs = [], set()
    for c in ctl:
        if len(controls) >= N_CONTROLS:
            break
        if c['traj'] in used or c['traj'] in ctl_trajs:
            continue
        controls.append(c)
        ctl_trajs.add(c['traj'])

    final = picked + controls
    print(f"{'#':>3} {'design':<26}{'len':>4}{'H409':>7}{'via':>6}"
          f"{'i_pTM':>7}{'pLDDT':>7}{'mouse':>7}  role")
    print('-' * 80)
    for i, d in enumerate(final, 1):
        role = 'salt bridge' if d['dist'] <= SALT else 'CONTROL (no switch)'
        print(f"{i:>3} {d['design'].replace('EGFR_',''):<26}{d['length']:>4}"
              f"{d['dist']:>7.2f}{d['partner']:>6}{d['iptm']:>7.2f}"
              f"{d['plddt']:>7.2f}{d['mouse']:>7.2f}  {role}")

    # pairwise identity check across the whole set
    print("\nmaximum pairwise identity within the submission:")
    worst = (0, None, None)
    for i in range(len(final)):
        for j in range(i + 1, len(final)):
            idt = identity(final[i]['seq'], final[j]['seq'])
            if idt > worst[0]:
                worst = (idt, final[i]['design'], final[j]['design'])
    print(f"  {100*worst[0]:.1f}%  {worst[1]} vs {worst[2]}")

    import json
    json.dump(final, open('run3/selection.json', 'w'), indent=1)
    print(f"\n{len(final)} designs selected -> run3/selection.json")


if __name__ == '__main__':
    main()
