"""Step 4: solve every mesh with FEniCSx (P2), and the beams with P1 too.

Runs inside the DOLFINx container (CI); reads the .msh files and JSON sidecars
made in step 2. Writes out/results/fenicsx_*.json.
"""
import argparse

import _common  # noqa: F401
from cadfea import pipeline
from cadfea.params import CASES, LEVELS, levels


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", default=",".join(CASES))
    ap.add_argument("--levels", default="", help="default: every level of each case")
    ap.add_argument("--degree", type=int, default=2, choices=[1, 2])
    ap.add_argument("--with-p1", action="store_true", help="also run the beams with P1")
    ap.add_argument("--solver", default="lu", choices=["lu", "cg"])
    args = ap.parse_args()
    runs = [(c, lv, args.degree) for c in args.cases.split(",")
            for lv in (args.levels.split(",") if args.levels else levels(c))]
    if args.with_p1:
        runs += [(c, lv, 1) for c in ("beam", "stocky") for lv in LEVELS]
    print(f"{'case':8s} {'level':7s} {'space':5s} {'dofs':>8s} {'seconds':>8s}  tip displacement per load case (mm)")
    for case, level, degree in runs:
        r = pipeline.ensure_fenicsx(case, level, degree, args.solver)
        tips = "  ".join(f"{lc}: {v['tip_u'][0]:+.6e} {v['tip_u'][1]:+.6e} {v['tip_u'][2]:+.6e}"
                         for lc, v in r["load_cases"].items())
        print(f"{case:8s} {level:7s} P{degree:<4d} {r['n_dofs']:8d} {r['seconds']:8.2f}  {tips}", flush=True)


if __name__ == "__main__":
    main()
