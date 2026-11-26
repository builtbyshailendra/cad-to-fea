"""Step 2: mesh every case at three levels with Gmsh (from the CadQuery STEP files).

Writes out/mesh/* and results/meshes.txt.
"""
import json

import _common  # noqa: F401
from cadfea import paths, pipeline, props
from cadfea.params import CASES, LEVELS


def main() -> None:
    paths.ensure(paths.MESH_DIR, paths.RESULTS_DIR, paths.REPORT_DIR)
    rows, out = [], {}
    for case in CASES:
        exact = props.analytic(CASES[case])["volume"]
        for level in LEVELS:
            info = pipeline.ensure_mesh(case, level)
            keep = {k: v for k, v in info.items() if k not in ("nsets", "loads", "load_tris", "load_xyz", "feature_xyz")}
            keep["nset_sizes"] = {k: len(v) for k, v in info["nsets"].items()}
            out[f"{case}_{level}"] = keep
            rows.append(
                f"{case:8s} {level:7s} {info['size']:5.1f} {info['feature_size']:5.1f} "
                f"{info['n_nodes']:8d} {info['n_tets']:8d} {info['n_dofs']:9d} "
                f"{abs(info['volume_straight'] / exact - 1):9.1e} {abs(info['volume_curved'] / exact - 1):9.1e} "
                f"{info['min_quality_sicn']:6.3f} {info['seconds']:6.2f}"
            )
    (paths.RESULTS_DIR / "meshes.json").write_text(json.dumps(out, indent=1))
    head = ("case     level    size  feat    nodes     tets      dofs  vol err   vol err  min    mesh\n"
            "                  (mm)  (mm)                                straight  curved   SICN   (s)")
    text = head + "\n" + "\n".join(rows) + "\n"
    (paths.REPORT_DIR / "meshes.txt").write_text(text)
    print(text)


if __name__ == "__main__":
    main()
