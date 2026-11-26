"""Where generated files go. Big outputs (STEP, meshes, solver files) stay out of git."""
from __future__ import annotations

import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = Path(os.environ.get("CADFEA_OUT", ROOT / "out"))
STEP_DIR = OUT / "step"
MESH_DIR = OUT / "mesh"
CCX_DIR = OUT / "ccx"
FENICS_DIR = OUT / "fenicsx"
RESULTS_DIR = OUT / "results"  # small JSON summaries per run
REPORT_DIR = ROOT / "results"  # committed: text summaries and plots


def ensure(*dirs: Path) -> None:
    for d in dirs:
        d.mkdir(parents=True, exist_ok=True)


def step_path(case: str, tool: str) -> Path:
    return STEP_DIR / f"{case}_{tool}.step"


def mesh_stem(case: str, level: str) -> str:
    return f"{case}_{level}"
