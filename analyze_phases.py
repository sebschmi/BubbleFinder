#!/usr/bin/env python3
"""Per-phase analysis of SP-Compress instrumentation CSV.

Usage: python3 analyze_phases.py <csv_path>"""
import csv, sys

def fmt(us):
    if us < 1000: return f"{us:.0f}us"
    if us < 1e6: return f"{us/1000:.1f}ms"
    return f"{us/1e6:.2f}s"

def main():
    if len(sys.argv) < 2:
        print(__doc__); sys.exit(1)
    with open(sys.argv[1]) as f:
        rows = list(csv.DictReader(f))
    n = len(rows)
    if n == 0: print("Empty CSV."); return

    legacy_cols = ['t_compress_us', 't_build_spqr_core_us', 't_reconstruct_us', 't_normalize_us', 't_canonicalize_us']
    prod_reconstruct = [
        't_reconstruct_build_builder_us',
        't_reconstruct_normalize_in_place_us',
        't_reconstruct_finalize_us',
        't_reconstruct_defensive_normalize_us',
    ]
    prod_cols = ['t_compress_us', 't_build_spqr_core_us'] + prod_reconstruct + ['t_canonicalize_us']
    sub_canon = ['t_canon_root_us', 't_canon_node_order_us', 't_canon_edge_orient_us', 't_canon_move_root_us']

    if legacy_cols[0] not in rows[0]:
        print("ERROR: this CSV doesn't have per-phase timings.")
        sys.exit(2)

    has_prod_reconstruct = all(c in rows[0] for c in prod_reconstruct)
    cols = prod_cols if has_prod_reconstruct else legacy_cols
    has_sub_canon = sub_canon[0] in rows[0]
    labels = {
        't_compress_us': 'compress',
        't_build_spqr_core_us': 'build spqr core',
        't_reconstruct_us': 'reconstruct',
        't_normalize_us': 'normalize',
        't_reconstruct_build_builder_us': 'reconstruct build builder',
        't_reconstruct_normalize_in_place_us': 'normalize in place',
        't_reconstruct_finalize_us': 'finalize CSR',
        't_reconstruct_defensive_normalize_us': 'defensive normalize',
        't_canonicalize_us': 'canonicalize',
    }

    print(f"=== Global ({n} blocks) ===")
    total = {c: sum(int(r[c]) for r in rows) for c in cols}
    if has_sub_canon:
        for c in sub_canon:
            total[c] = sum(int(r[c]) for r in rows)
    total['t_baseline'] = sum(int(r['t_baseline_us']) for r in rows)
    total['t_spcompress'] = sum(int(r['t_spcompress_us']) for r in rows)
    print(f"  {'baseline (build_spqr direct)':<35} {fmt(total['t_baseline']):>10}")
    print(f"  {'sp-compress total':<35} {fmt(total['t_spcompress']):>10}")
    if has_prod_reconstruct:
        print(f"  --- production-path breakdown of sp-compress ---")
    else:
        print(f"  --- legacy breakdown of sp-compress ---")
    for c in cols:
        share = 100 * total[c] / total['t_spcompress'] if total['t_spcompress'] > 0 else 0
        label = labels.get(c, c.replace('_us', '').replace('t_', '').replace('_', ' '))
        print(f"  {label:<35} {fmt(total[c]):>10}   {share:>5.1f}%")
    if has_sub_canon:
        print(f"  --- canonicalize sub-breakdown ---")
        for c in sub_canon:
            share_total = 100 * total[c] / total['t_spcompress'] if total['t_spcompress'] > 0 else 0
            share_canon = 100 * total[c] / total['t_canonicalize_us'] if total['t_canonicalize_us'] > 0 else 0
            label = c.replace('_us', '').replace('t_canon_', '').replace('_', ' ')
            print(f"  └─ {label:<32} {fmt(total[c]):>10}   {share_canon:>5.1f}% of canon, {share_total:>5.1f}% of total")
    sum_phases = sum(total[c] for c in cols)
    gap = total['t_spcompress'] - sum_phases
    print(f"  (sum of phases: {fmt(sum_phases)}; unaccounted wrapper/timer gap: {fmt(gap)})")
    print()

    # Big block focus
    big = max(rows, key=lambda r: int(r['n_edges']))
    if int(big['n_edges']) > 10000:
        print(f"=== The biggest block (block_idx={big['block_idx']}, n_edges={int(big['n_edges']):,}) ===")
        print(f"  baseline (build_spqr Gblk)     {fmt(int(big['t_baseline_us'])):>10}")
        print(f"  sp-compress total              {fmt(int(big['t_spcompress_us'])):>10}")
        print(f"  --- {'production path' if has_prod_reconstruct else 'main phases'} ---")
        for c in cols:
            v = int(big[c])
            share = 100 * v / int(big['t_spcompress_us'])
            label = labels.get(c, c.replace('_us', '').replace('t_', '').replace('_', ' '))
            print(f"  {label:<28} {fmt(v):>10}  {share:>5.1f}%")
        if has_sub_canon:
            print(f"  --- canonicalize sub-breakdown ---")
            for c in sub_canon:
                v = int(big[c])
                share_total = 100 * v / int(big['t_spcompress_us'])
                share_canon = 100 * v / int(big['t_canonicalize_us']) if int(big['t_canonicalize_us']) > 0 else 0
                label = c.replace('_us', '').replace('t_canon_', '').replace('_', ' ')
                print(f"  └─ {label:<25} {fmt(v):>10}  {share_canon:>5.1f}% canon, {share_total:>5.1f}% total")
        print()

        # Honest verdict
        b_spqr = int(big['t_build_spqr_core_us'])
        b_base = int(big['t_baseline_us'])
        ratio = b_spqr / b_base if b_base > 0 else 0
        print(f"  KEY INSIGHTS:")
        print(f"   build_spqr(coreG) / build_spqr(Gblk) = {ratio:.2f}")
        if ratio > 0.85: print(f"   → R-component dominates HT")
        elif ratio < 0.65: print(f"   → build_spqr is roughly linear in m")
        else: print(f"   → mixed")

        norm = int(big['t_reconstruct_normalize_in_place_us'] if has_prod_reconstruct else big['t_normalize_us'])
        canon = int(big['t_canonicalize_us'])
        if norm + canon > int(big['t_spcompress_us']) * 0.5:
            print(f"   → normalize+canonicalize = {fmt(norm+canon)} dominates {100*(norm+canon)/int(big['t_spcompress_us']):.0f}% of sp-compress total")

        if has_prod_reconstruct:
            sr = [(labels[c], int(big[c])) for c in prod_reconstruct]
            sr.sort(key=lambda x: -x[1])
            print(f"   → top reconstruct-path culprit: {sr[0][0]} ({fmt(sr[0][1])})")

        if has_sub_canon and int(big['t_canonicalize_us']) > 1000:
            sc = [(c, int(big[c])) for c in sub_canon]
            sc.sort(key=lambda x: -x[1])
            top = sc[0]
            top_label = top[0].replace('_us', '').replace('t_canon_', '').replace('_', ' ')
            top_share = 100 * top[1] / int(big['t_canonicalize_us'])
            print(f"   → top canon culprit: {top_label} ({fmt(top[1])} = {top_share:.0f}% of canon)")

if __name__ == '__main__':
    main()
