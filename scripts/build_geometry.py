"""Step 1: build every model in CadQuery and in build123d, export STEP, compare.

Writes out/results/geometry.json and results/geometry.txt.
"""
import json
import time

import _common  # noqa: F401
from cadfea import geom_b123d, geom_cq, paths, props
from cadfea.params import CASES


def _z(x: float) -> str:
    """Six decimals, without printing -0.000000 for a value like -5e-17."""
    return f"{round(x, 9) + 0.0:.6f}"


def main() -> None:
    paths.ensure(paths.STEP_DIR, paths.RESULTS_DIR, paths.REPORT_DIR)
    report, lines = {}, []
    for name, p in CASES.items():
        t0 = time.perf_counter()
        cq_model = geom_cq.build(p)
        t_cq = time.perf_counter() - t0
        t0 = time.perf_counter()
        bd_model = geom_b123d.build(p)
        t_bd = time.perf_counter() - t0
        a, b = geom_cq.properties(cq_model), geom_b123d.properties(bd_model)
        sa = geom_cq.export_step(cq_model, paths.step_path(name, "cadquery"))
        sb = geom_b123d.export_step(bd_model, paths.step_path(name, "build123d"))
        ra, rb = props.step_properties(sa), props.step_properties(sb)
        exact = props.analytic(p)
        r = {
            "cadquery": a, "build123d": b, "step_cadquery": ra, "step_build123d": rb,
            "closed_form": exact,
            "diff_cq_vs_b123d": props.rel_diff(a, b),
            "diff_step_vs_step": props.rel_diff(ra, rb),
            "diff_cq_vs_closed_form": props.rel_diff(a, exact),
            "diff_b123d_vs_closed_form": props.rel_diff(b, exact),
            "diff_cq_vs_own_step": props.rel_diff(a, ra),
            "seconds_cadquery": round(t_cq, 3), "seconds_build123d": round(t_bd, 3),
            "step_bytes": {"cadquery": sa.stat().st_size, "build123d": sb.stat().st_size},
        }
        report[name] = r
        worst = lambda d: max(d.values())  # noqa: E731
        lines += [
            f"{name}: volume {a['volume']:.6f} mm3, area {a['area']:.6f} mm2, "
            f"centre of mass ({_z(a["com"][0])}, {_z(a["com"][1])}, {_z(a["com"][2])}) mm, "
            f"{a['n_faces']} faces, {a['n_edges']} edges",
            f"  closed form: volume {exact['volume']:.6f}, area {exact['area']:.6f}",
            f"  largest relative difference, CadQuery vs build123d (own APIs): {worst(r['diff_cq_vs_b123d']):.1e}",
            f"  largest relative difference, the two STEP files read back:     {worst(r['diff_step_vs_step']):.1e}",
            f"  largest relative difference, CadQuery vs closed form:           {worst(r['diff_cq_vs_closed_form']):.1e}",
            f"  largest relative difference, build123d vs closed form:          {worst(r['diff_b123d_vs_closed_form']):.1e}",
            f"  build time: CadQuery {t_cq*1000:.0f} ms, build123d {t_bd*1000:.0f} ms; "
            f"STEP {sa.stat().st_size} and {sb.stat().st_size} bytes",
        ]
    (paths.RESULTS_DIR / "geometry.json").write_text(json.dumps(report, indent=1))
    text = "\n".join(lines) + "\n"
    (paths.REPORT_DIR / "geometry.txt").write_text(text)
    print(text)


if __name__ == "__main__":
    main()
