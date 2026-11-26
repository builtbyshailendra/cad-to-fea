"""FEniCSx (DOLFINx 0.10): linear elasticity on the same Gmsh mesh.

CalculiX is told "these nodes, these forces". FEniCSx is told the equations:
find u such that  integral(sigma(u) : eps(v)) dx = integral(T . v) ds(LOAD)
for every test function v, with u = 0 on FIXED. The form compiler (FFCx) turns
that into element kernels; PETSc solves the linear system.

Only importable where DOLFINx is installed (the CI container); everything in
here imports it lazily.
"""
from __future__ import annotations

import time
from pathlib import Path

import numpy as np

from .mesh import GROUP_TAGS, load_vectors
from .params import CASES, STEEL, Material


def _read_mesh(path: Path):
    from mpi4py import MPI

    try:  # DOLFINx >= 0.10
        from dolfinx.io import gmsh as gmshio
    except ImportError:  # 0.9 and older
        from dolfinx.io import gmshio
    data = gmshio.read_from_msh(str(path), MPI.COMM_WORLD, 0, gdim=3)
    if hasattr(data, "mesh"):
        return data.mesh, data.facet_tags
    return data[0], data[2]


def solve(msh_path: Path, info: dict, degree: int = 2, mat: Material = STEEL,
          solver: str = "lu") -> dict:
    import ufl
    from dolfinx import default_scalar_type, fem, geometry
    from dolfinx.fem.petsc import LinearProblem, assemble_vector

    from . import post

    p = CASES[info["case"]]
    t0 = time.perf_counter()
    domain, ft = _read_mesh(msh_path)
    t_read = time.perf_counter() - t0
    V = fem.functionspace(domain, ("Lagrange", degree, (3,)))

    E, nu = mat.E, mat.nu
    mu = fem.Constant(domain, default_scalar_type(E / (2 * (1 + nu))))
    lmbda = fem.Constant(domain, default_scalar_type(E * nu / ((1 + nu) * (1 - 2 * nu))))

    def eps(w):
        return ufl.sym(ufl.grad(w))

    def sigma(w):
        return 2 * mu * eps(w) + lmbda * ufl.tr(eps(w)) * ufl.Identity(3)

    fdim = domain.topology.dim - 1
    fixed = ft.find(GROUP_TAGS["FIXED"])
    fixed_dofs = fem.locate_dofs_topological(V, fdim, fixed)  # block (node) indices
    bc = fem.dirichletbc(np.zeros(3, dtype=default_scalar_type), fixed_dofs, V)

    ds = ufl.Measure("ds", domain=domain, subdomain_data=ft)
    load_area = fem.assemble_scalar(fem.form(1.0 * ds(GROUP_TAGS["LOAD"])))
    T = fem.Constant(domain, np.zeros(3, dtype=default_scalar_type))
    u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
    a = ufl.inner(sigma(u), eps(v)) * ufl.dx
    L = ufl.dot(T, v) * ds(GROUP_TAGS["LOAD"])

    if solver == "lu":
        opts = {"ksp_type": "preonly", "pc_type": "lu", "pc_factor_mat_solver_type": "mumps"}
    else:  # conjugate gradients with algebraic multigrid
        opts = {"ksp_type": "cg", "pc_type": "gamg", "ksp_rtol": 1e-12, "ksp_max_it": 2000}
    opts["ksp_error_if_not_converged"] = True
    problem = LinearProblem(a, L, bcs=[bc], petsc_options=opts,
                            petsc_options_prefix=f"elastic_{info['case']}_")
    if solver != "lu":
        _set_near_nullspace(problem, V)

    # Tip point: a mesh node (embedded by Gmsh), so evaluation is exact.
    tip = np.array([p.tip_point])
    tree = geometry.bb_tree(domain, domain.topology.dim)
    cells = geometry.compute_colliding_cells(domain, geometry.compute_collisions_points(tree, tip), tip)
    tip_cell = cells.links(0)[0]

    out = {"solver": "fenicsx", "element": f"P{degree}", "case": info["case"],
           "level": info["level"], "n_nodes": info["n_nodes"], "n_tets": info["n_tets"],
           "n_dofs": int(V.dofmap.index_map.size_global * V.dofmap.index_map_bs),
           "load_area": float(load_area), "seconds_read": round(t_read, 2), "load_cases": {}}
    t_solve = 0.0
    for lc, total in load_vectors(p).items():
        T.value[:] = np.array(total) / load_area
        t1 = time.perf_counter()
        uh = problem.solve()
        t_solve += time.perf_counter() - t1
        tip_u = uh.eval(tip, np.array([tip_cell], dtype=np.int32))[0]
        # Consistent reaction: the residual a(u_h, v) - L(v) on the fixed dofs.
        r = assemble_vector(fem.form(ufl.inner(sigma(uh), eps(v)) * ufl.dx - L))
        r.assemble()
        reaction = r.array.reshape(-1, 3)[fixed_dofs].sum(axis=0)
        nodal = nodal_stress(domain, uh, sigma)
        out["load_cases"][lc] = {
            "tip_u": tip_u.tolist(),
            "reaction": reaction.tolist(),
            "max_abs_u": float(np.abs(uh.x.array).max()),
            "max_von_mises": float(post.von_mises(nodal["stress"]).max()),
            **post.feature_stresses(info, nodal),
        }
    out["seconds"] = round(t_solve, 2)
    return out


def _set_near_nullspace(problem, V) -> None:
    """Rigid body modes for algebraic multigrid (3 translations, 3 rotations)."""
    from dolfinx import la
    from dolfinx.fem import Function
    from petsc4py import PETSc

    x = V.tabulate_dof_coordinates()
    basis = [la.vector(V.dofmap.index_map, bs=3) for _ in range(6)]
    for i in range(3):
        basis[i].array.reshape(-1, 3)[:, i] = 1.0
    b = [bb.array.reshape(-1, 3) for bb in basis[3:]]
    b[0][:, 0], b[0][:, 1] = -x[:, 1], x[:, 0]
    b[1][:, 1], b[1][:, 2] = -x[:, 2], x[:, 1]
    b[2][:, 0], b[2][:, 2] = x[:, 2], -x[:, 0]
    la.orthonormalize(basis)
    vecs = [Function(V) for _ in range(6)]
    for f, bb in zip(vecs, basis):
        f.x.array[:] = bb.array
    ns = PETSc.NullSpace().create(vectors=[f.x.petsc_vec for f in vecs])
    problem.A.setNearNullSpace(ns)


def nodal_stress(domain, uh, sigma) -> dict:
    """Stress at every mesh node, averaged over the cells around it.

    sigma(u_h) is interpolated into a discontinuous quadratic tensor space, whose
    10 points per cell are the cell's 10 geometry nodes (corners and mid-edges).
    For a straight 10-node tet the stress of a P2 displacement is linear, so these
    values are exact; CalculiX extrapolates its 4 integration points to the
    corners the same way. Averaging the cells' values at each node gives the
    smoothed nodal stress both programs report.
    """
    from dolfinx import fem

    W = fem.functionspace(domain, ("Discontinuous Lagrange", 2, (3, 3)))
    pts = W.element.interpolation_points
    pts = pts() if callable(pts) else pts
    s = fem.Function(W)
    s.interpolate(fem.Expression(sigma(uh), pts))
    vals = s.x.array.reshape(-1, 9)
    geo = domain.geometry.dofmap  # (cells, 10) geometry node per cell node
    ncell = domain.topology.index_map(domain.topology.dim).size_local
    cell_dofs = np.array([W.dofmap.cell_dofs(c) for c in range(ncell)])
    # Check the claim above: DG2 point i of a cell is that cell's geometry node i.
    gap = np.abs(W.tabulate_dof_coordinates()[cell_dofs] - domain.geometry.x[geo[:ncell]]).max()
    if gap > 1e-8:
        raise RuntimeError(f"DG2 points don't match geometry nodes (gap {gap})")
    nnode = domain.geometry.x.shape[0]
    acc = np.zeros((nnode, 9))
    cnt = np.zeros(nnode)
    np.add.at(acc, geo[:ncell].ravel(), vals[cell_dofs.ravel()])
    np.add.at(cnt, geo[:ncell].ravel(), 1.0)
    avg = acc / np.maximum(cnt, 1)[:, None]
    # Gmsh node tag = original input index + 1 (dolfinx sorts nodes by tag).
    tags = np.asarray(domain.geometry.input_global_indices, dtype=np.int64) + 1
    six = avg[:, [0, 4, 8, 1, 5, 2]]  # xx, yy, zz, xy, yz, zx
    return {"nodes": tags, "stress": six, "xyz": domain.geometry.x.copy()}
