"""Step 3: solve every mesh with CalculiX (C3D10), and the beams with C3D4 too.

Needs the ccx executable on PATH (CI: apt-get install calculix-ccx).
Writes out/ccx/* (decks, .dat, .frd) and out/results/ccx_*.json.
"""
import argparse

import _common  # noqa: F401
from cadfea import pipeline
from cadfea.params import CASES, LEVELS, NU0, STEEL, levels


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", default=",".join(CASES))
    ap.add_argument("--levels", default="", help="default: every level of each case")
    ap.add_argument("--element", default="C3D10", choices=["C3D10", "C3D4"])
    ap.add_argument("--with-c3d4", action="store_true", help="also run the beams with C3D4")
    ap.add_argument("--with-nu0", action="store_true", help="also run the beams with Poisson's ratio 0")
    args = ap.parse_args()
    runs = [(c, lv, args.element, STEEL) for c in args.cases.split(",")
            for lv in (args.levels.split(",") if args.levels else levels(c))]
    if args.with_c3d4:
        runs += [(c, lv, "C3D4", STEEL) for c in ("beam", "stocky") for lv in LEVELS]
    if args.with_nu0:
        runs += [(c, lv, "C3D10", NU0) for c in ("beam", "stocky") for lv in LEVELS]
    print(f"{'case':8s} {'level':7s} {'element':7s} {'nu':>4s} {'nodes':>7s} {'seconds':>8s}  tip displacement per load case (mm)")
    for case, level, element, mat in runs:
        r = pipeline.ensure_ccx(case, level, element, mat)
        tips = "  ".join(f"{lc}: {v['tip_u'][0]:+.6e} {v['tip_u'][1]:+.6e} {v['tip_u'][2]:+.6e}"
                         for lc, v in r["load_cases"].items())
        print(f"{case:8s} {level:7s} {element:7s} {mat.nu:4.1f} {r['n_nodes']:7d} {r['seconds']:8.2f}  {tips}", flush=True)


if __name__ == "__main__":
    main()
