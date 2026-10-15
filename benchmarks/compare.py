#!/usr/bin/env python3
"""Compare baseline vs current benchmark results."""

import csv
import pathlib

BASE = pathlib.Path(__file__).parent / "base_results.csv"
CUR = pathlib.Path(__file__).parent / "current_results.csv"


def main():
    base = load(BASE)
    cur = load(CUR)

    tiers = sorted(set(base.keys()) | set(cur.keys()))

    header = f"{'Tier':<4} {'Name':<38} {'Baseline (ms)':>13} {'Current (ms)':>12} {'Speedup':>9}"
    print(header)
    print("-" * len(header))

    for t in tiers:
        b = base.get(t, {})
        c = cur.get(t, {})

        name = c.get("name", b.get("name", "?"))
        base_ms = float(b["median_ms"]) if "median_ms" in b else None
        cur_ms = float(c["median_ms"]) if "median_ms" in c else None

        speedup = (
            f"{base_ms / cur_ms:.2f}x"
            if (base_ms is not None and cur_ms is not None and cur_ms > 0)
            else "-"
        )
        base_str = f"{base_ms:>11.1f}" if base_ms is not None else "    -"
        cur_str = f"{cur_ms:>10.1f}" if cur_ms is not None else "    -"

        print(f"{t:<4} {name:<38} {base_str:>13} {cur_str:>12} {speedup:>9}")


def load(path):
    if not path.exists():
        raise SystemExit(f"missing results file: {path} (run benchmarks/benchmark.R <label>)")
    with open(path) as f:
        return {row["tier"]: row for row in csv.DictReader(f)}
        

if __name__ == "__main__":
    main()
