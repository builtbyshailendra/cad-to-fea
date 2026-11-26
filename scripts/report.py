"""Step 5: put the solver results next to the known answers.

Reads out/results/*.json (from steps 1-4) and writes results/summary.md,
results/cantilever.txt, results/bracket.txt and two plots.
"""
import json
import os

import _common  # noqa: F401
from cadfea import paths, pipeline, theory
from cadfea.params import BEAM, BRACKET, LEVELS, NU0, STOCKY, THIN

SOLVERS = [("calculix", "C3D10"), ("fenicsx", "P2"), ("calculix", "C3D4"), ("fenicsx", "P1")]
LABEL = {"C3D10": "CalculiX C3D10", "P2": "FEniCSx P2", "C3D4": "CalculiX C3D4", "P1": "FEniCSx P1"}


def result(solver, case, level, element):
    return pipeline.load_result(solver, case, level, element)


def cantilever_rows():
    rows = []
    for p in (BEAM, STOCKY):
        eb, ts = theory.euler_bernoulli_tip(p), theory.timoshenko_tip(p)
        for solver, el in SOLVERS:
            for lv in LEVELS:
                r = result(solver, p.name, lv, el)
                if r is None:
                    continue
                d = -r["load_cases"]["tip"]["tip_u"][2]
                rows.append({"case": p.name, "solver": LABEL[el], "level": lv, "dofs": r.get("n_dofs", 3 * r["n_nodes"]),
                             "tip": d, "err_timo": d / ts - 1, "err_eb": d / eb - 1, "seconds": r["seconds"]})
    return rows


def bracket_rows():
    rows = []
    kt_ref, s_net = theory.heywood_kt(BRACKET), theory.net_stress(BRACKET)
    for lv in LEVELS:
        a, b = result("calculix", "bracket", lv, "C3D10"), result("fenicsx", "bracket", lv, "P2")
        for name, r in (("CalculiX C3D10", a), ("FEniCSx P2", b)):
            if r is None:
                continue
            sh, pu = r["load_cases"]["shelf"], r["load_cases"]["pull"]
            rows.append({
                "solver": name, "level": lv, "dofs": r.get("n_dofs", 3 * r["n_nodes"]),
                "feature": r.get("feature_size"),
                "shelf_tip_uz": sh["tip_u"][2], "pull_tip_ux": pu["tip_u"][0],
                "kt": pu["hole_sxx_max"] / s_net, "kt_face": pu["hole_sxx_face_max"] / s_net,
                "kt_err": pu["hole_sxx_max"] / s_net / kt_ref - 1,
                "shelf_fillet_s1": sh["fillet_s1_max"], "shelf_hole_s1": sh["hole_s1_max"],
                "profile": pu["hole_sxx_profile"], "seconds": r["seconds"],
            })
    return rows


def poisson_rows():
    """The beams again with nu = 0: is the gap to Timoshenko the clamped root's Poisson effect?"""
    rows = []
    for p in (BEAM, STOCKY):
        for lv in LEVELS:
            r = pipeline.load_result("calculix", p.name, lv, "C3D10", NU0)
            if r is None:
                continue
            d = -r["load_cases"]["tip"]["tip_u"][2]
            ts, eb = theory.timoshenko_tip(p, NU0), theory.euler_bernoulli_tip(p, NU0)
            rows.append({"case": p.name, "level": lv, "tip": d, "timo": ts, "err_timo": d / ts - 1, "err_eb": d / eb - 1})
    return rows


def thin_rows():
    """Kt against plate thickness: the 2 mm bracket next to the 8 mm one (fine meshes)."""
    rows = []
    for p in (THIN, BRACKET):
        for solver, el in (("calculix", "C3D10"), ("fenicsx", "P2")):
            r = result(solver, p.name, "fine", el)
            if r is None:
                continue
            pu, s_net = r["load_cases"]["pull"], theory.net_stress(p)
            prof = [s for _, s in pu["hole_sxx_profile"]]
            rows.append({"case": p.name, "t_over_d": p.t / p.d, "solver": LABEL[el],
                         "kt": pu["hole_sxx_max"] / s_net, "kt_face": pu["hole_sxx_face_max"] / s_net,
                         "kt_mean": sum(prof) / len(prof) / s_net, "heywood": theory.heywood_kt(p)})
    return rows


def same_mesh_gaps():
    gaps = []
    for case in ("beam", "stocky", "bracket", "thin"):
        for lv in (LEVELS if case != "thin" else ["fine"]):
            a, b = result("calculix", case, lv, "C3D10"), result("fenicsx", case, lv, "P2")
            if a is None or b is None:
                continue
            for lc in a["load_cases"]:
                ua, ub = a["load_cases"][lc]["tip_u"], b["load_cases"][lc]["tip_u"]
                big = max(range(3), key=lambda i: abs(ub[i]))
                gaps.append({"case": case, "level": lv, "load": lc, "ccx": ua[big], "fenicsx": ub[big],
                             "rel": abs(ua[big] - ub[big]) / abs(ub[big])})
    return gaps


def plots(cant, brk, nu0, thin) -> list[str]:
    try:
        import matplotlib
    except ImportError:
        return []
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullFormatter

    made = []
    colors = {"CalculiX C3D10": "#2a6fdb", "FEniCSx P2": "#e8702a", "CalculiX C3D4": "#7fa7ea", "FEniCSx P1": "#f2b48a"}
    marks = {"CalculiX C3D10": "o", "FEniCSx P2": "x", "CalculiX C3D4": "s", "FEniCSx P1": "+"}
    sizes = {"CalculiX C3D10": 7, "FEniCSx P2": 9, "CalculiX C3D4": 7, "FEniCSx P1": 11}

    # Cantilevers: everything (top), and the quadratic elements close up (bottom).
    fig, axes = plt.subplots(2, 2, figsize=(10, 7.5), dpi=90)
    for col, p in enumerate((BEAM, STOCKY)):
        ts, eb = theory.timoshenko_tip(p), theory.euler_bernoulli_tip(p)
        for row, labs in enumerate((list(colors), ["CalculiX C3D10", "FEniCSx P2"])):
            ax = axes[row, col]
            for lab in labs:
                pts = [(r["dofs"], r["tip"] / ts) for r in cant if r["case"] == p.name and r["solver"] == lab]
                if pts:
                    x, y = zip(*pts)
                    ax.plot(x, y, marker=marks[lab], color=colors[lab], label=f"{lab}, nu = 0.3", lw=1.2,
                            ms=sizes[lab], mfc="none" if lab.startswith("CalculiX") else None)
            if row == 1:
                dofs = {r["level"]: r["dofs"] for r in cant if r["case"] == p.name and r["solver"] == "CalculiX C3D10"}
                pts = [(dofs.get(r["level"]), r["tip"] / r["timo"]) for r in nu0 if r["case"] == p.name and r["level"] in dofs]
                if pts:
                    x, y = zip(*pts)
                    ax.plot(x, y, marker="D", color="#2a9d6f", label="CalculiX C3D10, nu = 0\n(vs Timoshenko at nu = 0)", lw=1.2)
            ax.axhline(1.0, color="#333", lw=1, label="Timoshenko (Cowper k)")
            ax.axhline(eb / ts, color="#333", lw=1, ls="--", label="Euler-Bernoulli")
            ax.set_xscale("log")
            ax.xaxis.set_minor_formatter(NullFormatter())
            ax.set_ylabel("tip deflection / Timoshenko")
            ax.grid(alpha=0.3)
            if row == 0:
                ax.set_title(f"{p.name}: L/h = {p.L / p.h:.0f}, all elements")
            else:
                ax.set_title(f"{p.name}: quadratic elements, close up")
                ax.set_xlabel("degrees of freedom")
        axes[0, 0].legend(fontsize=7, loc="lower right")
        axes[1, 0].legend(fontsize=7, loc="center right")
    fig.tight_layout()
    f = paths.REPORT_DIR / "cantilever_convergence.png"
    fig.savefig(f)
    plt.close(fig)
    made.append(f.name)

    if brk:
        fig, axes = plt.subplots(1, 2, figsize=(10, 4), dpi=90)
        kt_ref = theory.heywood_kt(BRACKET)
        for lab in ("CalculiX C3D10", "FEniCSx P2"):
            pts = [(r["feature"], r["kt"]) for r in brk if r["solver"] == lab]
            if pts:
                x, y = zip(*pts)
                axes[0].plot(x, y, marker=marks[lab], color=colors[lab], label=lab, ms=sizes[lab],
                             mfc="none" if lab.startswith("CalculiX") else None)
        axes[0].axhline(kt_ref, color="#333", lw=1, ls="--", label=f"Heywood (2D): {kt_ref:.3f}")
        axes[0].invert_xaxis()
        axes[0].set_xlabel("element size on the hole (mm)")
        axes[0].set_ylabel("Kt, peak (net section)")
        axes[0].set_title("8 mm bracket, pull load: Kt vs mesh")
        axes[0].grid(alpha=0.3)
        axes[0].legend(fontsize=8)
        for case, p, style in (("bracket", BRACKET, "-"), ("thin", THIN, ":")):
            s_net = theory.net_stress(p)
            for solver, el in (("calculix", "C3D10"), ("fenicsx", "P2")):
                r = result(solver, case, "fine", el)
                if r is None:
                    continue
                lab = LABEL[el]
                z, s = zip(*r["load_cases"]["pull"]["hole_sxx_profile"])
                axes[1].plot([v / s_net for v in s], [v / p.t for v in z], ls=style, marker=marks[lab],
                             color=colors[lab], ms=sizes[lab], mfc="none" if lab.startswith("CalculiX") else None,
                             label=f"{lab}, t/d = {p.t / p.d:.1f}")
        axes[1].axvline(kt_ref, color="#333", lw=1, ls="--", label="Heywood (2D)")
        axes[1].set_xlabel("peak hoop stress / net-section stress")
        axes[1].set_ylabel("z / t (0 and 1: the plate faces)")
        axes[1].set_title("fine meshes: through the thickness")
        axes[1].grid(alpha=0.3)
        axes[1].legend(fontsize=7)
        fig.tight_layout()
        f = paths.REPORT_DIR / "bracket_hole_kt.png"
        fig.savefig(f)
        plt.close(fig)
        made.append(f.name)
    return made


def main() -> None:
    paths.ensure(paths.REPORT_DIR)
    cant, brk, gaps = cantilever_rows(), bracket_rows(), same_mesh_gaps()
    where = os.environ.get("CADFEA_WHERE", "local")
    lines = [f"# Results ({where})", ""]
    th = theory.summary()
    lines += ["## Known answers", "",
              "| | L/h | Euler-Bernoulli (mm) | Timoshenko, Cowper k (mm) | Timoshenko, k = 5/6 (mm) | shear share |",
              "|---|---|---|---|---|---|"]
    for name in ("beam", "stocky"):
        t = th[name]
        lines.append(f"| {name} | {t['L_over_h']:.0f} | {t['euler_bernoulli']:.6f} | {t['timoshenko_cowper']:.6f} | "
                     f"{t['timoshenko_5_6']:.6f} | {100 * t['shear_share']:.2f}% |")
    b = th["bracket"]
    lines += ["", f"Bracket, pull load: net-section stress {b['net_stress']:.3f} MPa, Heywood Kt {b['heywood_kt']:.4f}.", ""]

    lines += ["## Cantilevers: tip deflection (centre of the tip face)", "",
              "| case | solver | mesh | dofs | tip (mm) | vs Timoshenko | vs Euler-Bernoulli | solve (s) |",
              "|---|---|---|---|---|---|---|---|"]
    for r in cant:
        lines.append(f"| {r['case']} | {r['solver']} | {r['level']} | {r['dofs']:,} | {r['tip']:.6f} | "
                     f"{100 * r['err_timo']:+.3f}% | {100 * r['err_eb']:+.3f}% | {r['seconds']:.2f} |")
    lines += ["", "## Bracket", "",
              "| solver | mesh | dofs | hole size (mm) | shelf: tip uz (mm) | pull: tip ux (mm) | Kt | Kt vs Heywood | Kt at the faces | shelf: fillet s1 (MPa) | shelf: hole s1 (MPa) | solve (s) |",
              "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in brk:
        lines.append(f"| {r['solver']} | {r['level']} | {r['dofs']:,} | {r['feature']} | {r['shelf_tip_uz']:.6f} | "
                     f"{r['pull_tip_ux']:.7f} | {r['kt']:.4f} | {100 * r['kt_err']:+.2f}% | {r['kt_face']:.4f} | "
                     f"{r['shelf_fillet_s1']:.2f} | {r['shelf_hole_s1']:.2f} | {r['seconds']:.2f} |")
    nu0 = poisson_rows()
    if nu0:
        lines += ["", "## Experiment: the beams with Poisson's ratio 0 (CalculiX C3D10)", "",
                  "| case | mesh | tip (mm) | Timoshenko, nu = 0 (mm) | vs Timoshenko | vs Euler-Bernoulli |", "|---|---|---|---|---|---|"]
        for r in nu0:
            lines.append(f"| {r['case']} | {r['level']} | {r['tip']:.6f} | {r['timo']:.6f} | {100 * r['err_timo']:+.3f}% | {100 * r['err_eb']:+.3f}% |")
    thin = thin_rows()
    if thin:
        lines += ["", "## Experiment: hole Kt against plate thickness (pull load, fine meshes)", "",
                  "| case | t/d | solver | Kt, peak | Kt, mean of the 8 through-thickness bins | Kt at the faces | Heywood |",
                  "|---|---|---|---|---|---|---|"]
        for r in thin:
            lines.append(f"| {r['case']} | {r['t_over_d']:.1f} | {r['solver']} | {r['kt']:.4f} ({100 * (r['kt'] / r['heywood'] - 1):+.1f}%) | "
                         f"{r['kt_mean']:.4f} ({100 * (r['kt_mean'] / r['heywood'] - 1):+.1f}%) | {r['kt_face']:.4f} | {r['heywood']:.4f} |")
    lines += ["", "## Same mesh, two solvers: tip displacement (largest component)", "",
              "| case | mesh | load | CalculiX (mm) | FEniCSx (mm) | relative difference |", "|---|---|---|---|---|---|"]
    for g in gaps:
        lines.append(f"| {g['case']} | {g['level']} | {g['load']} | {g['ccx']:.9f} | {g['fenicsx']:.9f} | {g['rel']:.1e} |")
    made = plots(cant, brk, poisson_rows(), thin_rows())
    if made:
        lines += ["", "Plots: " + ", ".join(made)]
    text = "\n".join(lines) + "\n"
    (paths.REPORT_DIR / "summary.md").write_text(text)
    (paths.RESULTS_DIR / "report.json").write_text(json.dumps({"cantilever": cant, "bracket": brk, "same_mesh": gaps, "nu0": nu0, "thickness": thin}, indent=1))
    print(text)


if __name__ == "__main__":
    main()
