import argparse
import random
import subprocess
import tempfile
from pathlib import Path
import sys


def generate_random_gfa(gfa, nn, ne, seed=None):
    rng = random.Random(seed)

    nodes = list(range(1, nn + 1))
    links = set()

    with open(gfa, "w") as of:
        for n in nodes:
            of.write(f"S {n} x\n")

        while len(links) < ne:
            a = rng.choice(nodes)
            b = rng.choice(nodes)
            if a == b:
                continue
            sa = rng.choice(["+", "-"])
            sb = rng.choice(["+", "-"])
            link = (a, sa, b, sb)
            if link in links:
                continue
            links.add(link)

        for (a, sa, b, sb) in sorted(links):
            of.write(f"L {a} {sa} {b} {sb} x\n")


def make_canonical_snarl(epts):
    def flip(s): return "+" if s == "-" else "-"
    norm = tuple(sorted(epts))
    inv = tuple(sorted((n, flip(s)) for n, s in epts))
    return min(norm, inv)


def parse_snarls_bruteforce_output(txt):
    sns = set()
    for line in txt.splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) != 2:
            continue
        e1_raw, e2_raw = parts
        if len(e1_raw) < 2 or len(e2_raw) < 2:
            continue
        try:
            n1, s1 = int(e1_raw[:-1]), e1_raw[-1]
            n2, s2 = int(e2_raw[:-1]), e2_raw[-1]
        except ValueError:
            continue
        sn = make_canonical_snarl([(n1, s1), (n2, s2)])
        sns.add(sn)
    return sns


def parse_snarls_bubblefinder_file(path: Path):
    sns = set()

    with open(path) as f:
        lines = [l.strip() for l in f if l.strip()]

    if not lines:
        return sns

    start = 0
    first = lines[0].split()
    if len(first) == 1 and first[0].isdigit():
        start = 1

    for line in lines[start:]:
        parts = line.split()
        ep = []

        for token in parts:
            token = token.strip()
            if not token:
                continue
            if len(token) < 2:
                continue
            try:
                node = int(token[:-1])
                sign = token[-1]
            except ValueError:
                continue
            if sign not in {"+", "-"}:
                continue
            ep.append((node, sign))

        k = len(ep)
        if k < 2:
            continue
        for i in range(k):
            for j in range(i + 1, k):
                sns.add(make_canonical_snarl((ep[i], ep[j])))

    return sns


def parse_superbubbles_bruteforce_output(txt):
    bbs = set()
    lines = [l.strip() for l in txt.splitlines() if l.strip()]
    if not lines:
        return bbs

    start = 0
    if lines and len(lines[0].split()) == 1 and lines[0].split()[0].isdigit():
        start = 1

    for line in lines[start:]:
        parts = line.split()
        if len(parts) != 2:
            continue
        a, b = parts
        try:
            val_a = a[:-1] if a[-1] in "+-" else a
            val_b = b[:-1] if b[-1] in "+-" else b
            u = int(val_a)
            v = int(val_b)
        except ValueError:
            continue
        if u == v:
            continue
        if u > v:
            u, v = v, u
        bbs.add((u, v))

    return bbs


def parse_superbubbles_bubblefinder_file(path: Path):
    bbs = set()

    with open(path) as f:
        lines = [l.strip() for l in f if l.strip()]

    if not lines:
        return bbs

    start = 0
    if lines and len(lines[0].split()) == 1 and lines[0].split()[0].isdigit():
        start = 1

    for line in lines[start:]:
        parts = line.split()
        if len(parts) != 2:
            continue
        a, b = parts
        try:
            u = int(a)
            v = int(b)
        except ValueError:
            continue
        if u == v:
            continue
        if u > v:
            u, v = v, u
        bbs.add((u, v))

    return bbs


def _decode_bytes(blob: bytes) -> str:
    if not blob:
        return ""
    return blob.decode("utf-8", errors="replace")


def run_snarls_bf(bin_path, gfa_path):
    proc = subprocess.run(
        [bin_path, str(gfa_path)],
        capture_output=True,
        check=True,
    )
    return parse_snarls_bruteforce_output(_decode_bytes(proc.stdout))


def run_superbubbles_bf(bin_path, gfa_path):
    proc = subprocess.run(
        [bin_path, str(gfa_path)],
        capture_output=True,
        check=True,
    )
    return parse_superbubbles_bruteforce_output(_decode_bytes(proc.stdout))


def run_bubblefinder(bin_path, gfa_path, out_path, threads, mode="snarls"):
    cmd = [
        bin_path,
        mode,
        "-g",
        str(gfa_path),
        "-o",
        str(out_path),
        "--gfa",
        "-j",
        str(threads),
    ]
    proc = subprocess.run(cmd, capture_output=True)

    sout = _decode_bytes(proc.stdout)
    serr = _decode_bytes(proc.stderr)

    if proc.returncode != 0:
        raise subprocess.CalledProcessError(
            proc.returncode,
            proc.args,
            output=sout,
            stderr=serr,
        )

    if mode == "superbubbles":
        return parse_superbubbles_bubblefinder_file(out_path)
    return parse_snarls_bubblefinder_file(out_path)


def compare_snarls(a, b):
    return a - b, b - a


def _has_divergence_reason(reasons):
    return any(r.startswith("divergence_") for r in reasons)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--bruteforce-bin",
        default="./snarls_bf",
        help="Path to brute-force binary (auto-detects mode if 'ultrabubbles' in name)",
    )
    parser.add_argument("--bubblefinder-bin", default="./BubbleFinder",
                        help="Path to BubbleFinder binary")
    parser.add_argument("--n-graphs", type=int, default=10000,
                        help="Number of random graphs to test")
    parser.add_argument("--min-nodes", type=int, default=3)
    parser.add_argument("--max-nodes", type=int, default=10)
    parser.add_argument("--min-edges", type=int, default=2)
    parser.add_argument("--max-edges", type=int, default=20)
    parser.add_argument("--seed", type=int, default=440)
    parser.add_argument(
        "--threads",
        type=int,
        nargs="+",
        default=[1],
        help="List of thread counts to test for BubbleFinder (e.g. --threads 1 7)"
    )
    parser.add_argument("--keep-failing", action="store_true",
                        help="Save failing .gfa files (and BF outputs) in /tmp")
    parser.add_argument(
        "--max-report-graphs",
        type=int,
        default=5,
        help="Max number of failing graphs whose full GFA is printed at the end"
    )
    parser.add_argument(
        "--superbubbles",
        action="store_true",
        help="Test superbubbles instead of snarls/ultrabubbles",
    )

    args = parser.parse_args()

    if args.superbubbles:
        mode = "superbubbles"
        feat = "superbubbles"
    elif "ultrabubbles" in Path(args.bruteforce_bin).name.lower():
        mode = "ultrabubbles"
        feat = "ultrabubbles"
    else:
        mode = "snarls"
        feat = "snarls"

    Feat = feat.capitalize()
    print(f"Tests on {Feat} (Mode detected: {mode}).")

    brute = args.bruteforce_bin
    bf_bin = args.bubblefinder_bin
    rng = random.Random(args.seed)
    tlist = sorted(set(args.threads))

    fails = 0
    notes = []

    for idx in range(args.n_graphs):
        nn = rng.randint(args.min_nodes, args.max_nodes)
        ne = rng.randint(args.min_edges, args.max_edges)

        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            gfa_path = tmp / "input.gfa"
            generate_random_gfa(gfa_path, nn, ne, seed=rng.randint(0, 10**9))

            bad = False
            reasons = set()

            try:
                if mode == "superbubbles":
                    bf_snarls = run_superbubbles_bf(brute, gfa_path)
                else:
                    bf_snarls = run_snarls_bf(brute, gfa_path)
            except subprocess.CalledProcessError as e:
                bad = True
                reasons.add("bruteforce_error")
                print(f"[Graph {idx}] ERROR in brute-force (returncode={e.returncode})")
                if e.output:
                    print("  [bruteforce stdout]:")
                    print("  <<<")
                    print(e.output.rstrip())
                    print("  >>>")
                if e.stderr:
                    print("  [bruteforce stderr]:")
                    print("  <<<")
                    print(e.stderr.rstrip())
                    print("  >>>")

                sample = None
                if len(notes) < args.max_report_graphs:
                    sample = gfa_path.read_text()
                if args.keep_failing:
                    dst = Path("/tmp") / f"fail_graph_{idx}.gfa"
                    dst.write_text(gfa_path.read_text())
                    dstout = Path("/tmp") / f"fail_graph_{idx}_bruteforce.stdout.txt"
                    dsterr = Path("/tmp") / f"fail_graph_{idx}_bruteforce.stderr.txt"
                    dstout.write_text(e.output or "")
                    dsterr.write_text(e.stderr or "")

                notes.append({
                    "index": idx,
                    "n_nodes": nn,
                    "n_edges": ne,
                    "reasons": sorted(reasons),
                    "threads": tlist,
                    "gfa": sample,
                    "has_divergence": _has_divergence_reason(reasons),
                })
                fails += 1
                continue

            per_thread = {}
            bf_err = False
            for th in tlist:
                out = tmp / f"out_t{th}.sbfind"
                try:
                    per_thread[th] = run_bubblefinder(
                        bf_bin, gfa_path, out, threads=th, mode=mode
                    )
                except subprocess.CalledProcessError as e:
                    bad = True
                    bf_err = True
                    reasons.add("bubblefinder_error")
                    print(f"[Graph {idx}] ERROR in BubbleFinder (threads={th}, returncode={e.returncode})")
                    if e.output:
                        print("  [BubbleFinder stdout]:")
                        print("  <<<")
                        print(e.output.rstrip())
                        print("  >>>")
                    if e.stderr:
                        print("  [BubbleFinder stderr]:")
                        print("  <<<")
                        print(e.stderr.rstrip())
                        print("  >>>")

                    sample = None
                    if len(notes) < args.max_report_graphs:
                        sample = gfa_path.read_text()
                    if args.keep_failing:
                        dstgfa = Path("/tmp") / f"fail_graph_{idx}.gfa"
                        dstgfa.write_text(gfa_path.read_text())
                        dstout = Path("/tmp") / f"fail_graph_{idx}_t{th}.bubblefinder.stdout.txt"
                        dsterr = Path("/tmp") / f"fail_graph_{idx}_t{th}.bubblefinder.stderr.txt"
                        dstout.write_text(e.output or "")
                        dsterr.write_text(e.stderr or "")

                    notes.append({
                        "index": idx,
                        "n_nodes": nn,
                        "n_edges": ne,
                        "reasons": sorted(reasons),
                        "threads": tlist,
                        "gfa": sample,
                        "has_divergence": _has_divergence_reason(reasons),
                    })
                    fails += 1
                    break

            if bf_err:
                continue

            divergent = False
            for th in tlist:
                missing, extra = compare_snarls(bf_snarls, per_thread[th])
                if missing or extra:
                    divergent = True
                    reasons.add("divergence_bruteforce_vs_bubblefinder")
                    print(f"[Graph {idx}] DIVERGENCE brute-force vs BubbleFinder (threads={th})")
                    print(f"  nodes={nn}, edges={ne}")
                    if missing:
                        print(f"  {Feat} missing in BubbleFinder:")
                        for sn in sorted(missing):
                            print("    ", sn)
                    if extra:
                        print(f"  Extra {feat} in BubbleFinder:")
                        for sn in sorted(extra):
                            print("    ", sn)

            if len(tlist) > 1:
                ref = tlist[0]
                ref_set = per_thread[ref]
                for th in tlist[1:]:
                    missing, extra = compare_snarls(ref_set, per_thread[th])
                    if missing or extra:
                        divergent = True
                        reasons.add("divergence_between_threads")
                        print(f"[Graph {idx}] DIVERGENCE between threads={ref} and threads={th}")
                        if missing:
                            print(f"  {Feat} missing in threads={th} vs {ref}:")
                            for sn in sorted(missing):
                                print("    ", sn)
                        if extra:
                            print(f"  Extra {feat} in threads={th} vs {ref}:")
                            for sn in sorted(extra):
                                print("    ", sn)

            if divergent:
                bad = True

            if bad:
                sample = None
                if len(notes) < args.max_report_graphs:
                    sample = gfa_path.read_text()
                if args.keep_failing:
                    dst = Path("/tmp") / f"fail_graph_{idx}.gfa"
                    dst.write_text(gfa_path.read_text())
                    print(f"  GFA saved to {dst}")

                notes.append({
                    "index": idx,
                    "n_nodes": nn,
                    "n_edges": ne,
                    "reasons": sorted(reasons),
                    "threads": tlist,
                    "gfa": sample,
                    "has_divergence": _has_divergence_reason(reasons),
                })
                fails += 1
            else:
                print(f"[Graph {idx}] OK (nodes={nn}, edges={ne}, threads={tlist})")

    print(f"Total number of graphs tested: {args.n_graphs}")
    print(f"Number of graphs with divergence or errors: {fails}")

    if fails == 0:
        print(f"\nAll {args.n_graphs} graphs passed: no logical divergences and no runtime errors detected.")
        return 0

    rc = {}
    for rec in notes:
        for reason in rec["reasons"]:
            rc[reason] = rc.get(reason, 0) + 1

    if rc:
        print("\nFailure causes summary:")
        for r, c in sorted(rc.items(), key=lambda x: -x[1]):
            print(f"  - {r}: {c} graphs")

    fail_idxs = [rec["index"] for rec in notes]
    print("\nIndices of failing graphs:")
    print("  ", ", ".join(str(i) for i in fail_idxs))

    sample = notes[:args.max_report_graphs]
    print(f"\nDetails for {len(sample)} failing graph(s):")
    for rec in sample:
        print("\n------------------------------")
        print(f"Graph {rec['index']} (nodes={rec['n_nodes']}, edges={rec['n_edges']})")
        print(f"Reasons: {', '.join(rec['reasons'])}")
        print(f"Thread counts tested: {rec['threads']}")
        if rec["gfa"] is not None:
            print("GFA:")
            print("<<<")
            print(rec["gfa"].rstrip())
            print(">>>")
        else:
            print("GFA not shown (beyond max-report-graphs).")

    divs = sum(1 for rec in notes if rec.get("has_divergence"))
    runtime = len(notes) - divs

    print("\n\n\nSummary of issues:")
    print(f"  Graph with logical divergences: {divs}")
    print(f"  - Graphs with runtime errors only: {runtime}")

    if runtime > 0:
        print("")
        print("  Some graphs failed due to segfaults or other behavior")

    if divs > 0:
        print("Divergences were detected.")
        return 1
    else:
        print("\nNo divergences detected (only runtime errors)")
        return 0


if __name__ == "__main__":
    sys.exit(main())