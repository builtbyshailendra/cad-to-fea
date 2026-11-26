# Learning notes: cad-to-fea

One idea runs through the whole repo: push one engineering result through several tools, and check it against an answer you already know. Two tools agreeing with each other, and with theory, is the evidence. When they don't agree, the gap is where the learning is.

```text
params.py ─► CadQuery (sketch + extrude) ─► STEP ─┐
          └► build123d (box + box)       ─► STEP ─┴─► same volume, area, centre of mass, inertia, bbox (1e-16)
CadQuery STEP ─► Gmsh: named faces, size field, 10-node tets ─► .msh ─► FEniCSx (P2) ─┐
                                                              └► .inp ─► CalculiX (C3D10) ─┴─► same tip deflection
                                                                                            ─► vs Timoshenko, vs Heywood Kt
```

## 1. CadQuery

**What it is.** A Python library for parametric CAD. You write a script; it calls OpenCASCADE (OCCT), the open-source geometry kernel, through the OCP bindings. The result is a B-rep solid you can export as STEP.

**The concepts that matter.**
- **Workplanes.** Every sketch lives on a plane. The named planes have fixed orientations: `"XZ"` maps local (x, y) to global (X, Z) and its normal is **-Y**, so extrusions go towards -Y unless you use `both=True`. Getting this wrong mirrors your part.
- **The stack and selectors.** Each call returns a new `Workplane` holding a stack of objects. Selectors pick from it: `edges("|Y")` (edges parallel to Y), `faces(">Z")` (the face furthest along +Z), or selector objects such as `NearestToPointSelector((t, 0, t))`.
- **Operations on solids.** `extrude`, `fillet`, `cut` and the other booleans work on the B-rep: faces, edges and vertices with exact surfaces (planes, cylinders) underneath.
- **B-rep, not mesh.** A STEP file stores exact surfaces. The fillet is a true cylinder; the mesh comes later.
- **Units.** STEP stores its length unit. Here it is millimetres (`SI_UNIT(.MILLI.,.METRE.)`), and Gmsh is told to read millimetres.

**Where it is used.** `cadfea/geom_cq.py`: the beam (`YZ` workplane, rectangle, extrude) and the bracket (L profile on `XZ`, symmetric extrude, fillet on the inner corner edge, cut a cylinder). `props.py` reads the STEP back with plain OCP.

**Pitfalls met.**
- A cutting cylinder that ends exactly on the part's faces gives coplanar faces in the boolean. Overshoot by 1 mm instead.
- `Shape.BoundingBox()` in CadQuery 2.8 takes a tolerance, not `optimal=True` (build123d does). Version differences like this are why the dependency versions are pinned.

## 2. build123d

**What it is.** Another Python CAD library on the same OCCT kernel, with a more Pythonic API: builder mode (context managers that collect what you create) and algebra mode (`part = box + box - cylinder`).

**The concepts that matter.**
- **Builders.** `with BuildPart() as bp:` collects solids; `mode=Mode.SUBTRACT` turns the next object into a cut.
- **Alignment.** `Box(L, W, t, align=(Align.MIN, Align.CENTER, Align.MIN))` says which corner or centre sits at the origin. The default is centred on all three axes, which is rarely what an engineering drawing means.
- **ShapeList filters.** `bp.edges().filter_by(Axis.Y).sort_by_distance((t, 0, t))[0]` is the inner corner edge: select by geometry, not by index.
- **Locations.** `with Locations((xh, 0, t/2)):` places the next object.
- **Same kernel, different route.** Here the bracket is two fused boxes, not an extruded profile, on purpose: if two constructions give the same solid, the parameters mean what I think.

**Where it is used.** `cadfea/geom_b123d.py`.

**Pitfalls met.**
- CadQuery depends on `cadquery-ocp` and build123d on `cadquery-ocp-novtk`; both install the same `OCP` module. In this venv they coexist (pip also installed `cadquery-ocp-proxy`), but it is the first thing to suspect if imports break.
- The two STEP files describe the same solid but list the faces in a different order with different orientations. Gmsh meshes them differently (6,308 vs 6,330 nodes for the coarse bracket, locally). That is why the pipeline meshes one file, the CadQuery one.

## 3. Gmsh

**What it is.** An open-source mesh generator with a CAD engine (OCC) and a Python API. It reads STEP, meshes in 1D/2D/3D, and writes many formats.

**The concepts that matter.**
- **Entities vs physical groups.** Entities are the CAD's points, curves, surfaces and volumes. Physical groups are named collections of entities: they are what solvers see (`FIXED`, `LOAD`, `HOLE`, `FILLET`, `SOLID`). If any physical group exists, Gmsh only saves elements that belong to one.
- **Find faces by geometry.** Entity tags depend on how the CAD tool ordered the B-rep. `classify_faces` uses each face's type (`Plane`, `Cylinder`) and bounding box instead.
- **Size fields.** A `Distance` field (distance to the hole and fillet surfaces) feeds a `Threshold` field (small size near, large size far), set as the background mesh. Turn off the other size sources (`MeshSizeFromPoints`, `MeshSizeFromCurvature`, `MeshSizeExtendFromBoundary`) or they fight the field.
- **Element order.** `setOrder(2)` turns 4-node tets into 10-node tets and moves the new mid-edge nodes onto the CAD surface, so the hole is curved in the mesh. Check the Jacobian determinant stays positive.
- **Embedding.** `mesh.embed(0, [point], 2, face)` forces a node at an exact point: both solvers report the tip deflection at the same node.
- **Formats.** `.msh` 4.1 keeps physical groups (FEniCSx reads it). `.inp` is Abaqus format (CalculiX reads it), and Gmsh reorders the 10-node tet's mid-edge nodes into Abaqus order on the way out.

**Where it is used.** `cadfea/mesh.py`; `scripts/build_meshes.py`.

**Pitfalls met.**
- Surface physical groups would be written to the `.inp` as 2D elements with no section, which CalculiX rejects. The `.inp` is written with only the volume group; the node sets go into the deck from the JSON sidecar.
- The same Gmsh version, same STEP and same options gave slightly different meshes on Windows and on Linux (the fine beam: 39,702 nodes locally, 39,626 in CI). Repeatable on one machine, not across machines.

## 4. CalculiX (ccx)

**What it is.** An open-source finite element solver whose input deck follows Abaqus conventions. You describe the model as cards in a text file, run `ccx -i job`, and read `job.dat` (printed tables) and `job.frd` (results for post-processing).

**The concepts that matter.**
- **The deck.** `*NODE`, `*ELEMENT` (from Gmsh), `*NSET`, `*MATERIAL` + `*ELASTIC`, `*SOLID SECTION` (ties elements to a material), `*BOUNDARY` (fixed DOFs), `*STEP` + `*STATIC`, `*CLOAD` (nodal forces), and output requests.
- **Output requests.** `*NODE PRINT` writes tables to `.dat` (easy to parse: tip displacement, reaction totals). `*NODE FILE` / `*EL FILE` write `.frd` (all nodes; stresses extrapolated from integration points to nodes and averaged).
- **Element choice.** C3D10 (quadratic tet) is accurate in bending. C3D4 (linear tet) has constant strain per element and is far too stiff in bending ("locking").
- **Consistent loads.** A uniform traction on a 6-node triangle puts **zero** force on the corner nodes and a third on each mid-edge node. Spreading the force evenly over all six nodes is a different (wrong) load.
- **Steps.** Loads carry over from one `*STEP` to the next unless you write `*CLOAD, OP=NEW`.
- **Equilibrium check.** The reaction total on the fixed face must equal minus the applied load. It did in every run; the other components stayed below 4e-7 N.

**Where it is used.** `cadfea/ccx.py`; `scripts/run_ccx.py`. CI only (`apt-get install calculix-ccx`, version 2.21 on Ubuntu 24.04).

**Pitfalls met.**
- `.frd` numbers have no separators: `1.00000E-03-2.50000E-01` is two numbers. Parse fixed columns (12 characters), never `split()`.
- For C3D4, the mid-edge nodes are still in the node list but belong to no element; they must leave the node sets, or CalculiX constrains and loads nodes that are not connected to anything.
- The default solver (SPOOLES) took 95-124 s for the 279,000-DOF bracket, against 36-48 s for FEniCSx with MUMPS: fine here, slow for bigger models.

## 5. FEniCSx (DOLFINx)

**What it is.** A finite element framework where you write the weak form of the equations in Python (UFL), and the form compiler (FFCx) generates the element code. DOLFINx is the library, Basix the element tabulator, PETSc the linear algebra.

**The concepts that matter.**
- **Weak form.** Find u with u = 0 on FIXED such that ∫ σ(u) : ε(v) dx = ∫ T · v ds(LOAD) for every v. `a` and `L` in `fenics.py` are exactly these two integrals.
- **Function spaces.** `functionspace(mesh, ("Lagrange", 2, (3,)))` is a vector P2 field: the same unknowns as C3D10 on the same mesh.
- **Boundary conditions from tags.** Gmsh's physical groups become facet tags; `locate_dofs_topological(V, 2, ft.find(11))` finds every DOF on the fixed face. The traction acts on `ds(12)`.
- **Solving.** `LinearProblem` with PETSc options: here direct LU (MUMPS). CG with algebraic multigrid also works for elasticity if you give GAMG the six rigid-body modes (the near-nullspace).
- **Post-processing.** Stress is not a solution field; interpolate the UFL expression σ(u_h) into a discontinuous space, then average at the nodes to compare with CalculiX.

**Where it is used.** `cadfea/fenics.py`; `scripts/run_fenicsx.py`. CI only, in the `ghcr.io/fenics/dolfinx/dolfinx:v0.10.0-r1` image (DOLFINx 0.10.0, Gmsh 4.14).

**Pitfalls met.**
- The API moves between versions: in 0.10 `dolfinx.io.gmshio` became `dolfinx.io.gmsh`, `read_from_msh` returns a `MeshData` tuple, `LinearProblem` needs a `petsc_options_prefix`, and `Function.eval` returned a flat array where the code expected one row per point.
- P1 (linear) elements lock in bending exactly like C3D4.
- DOLFINx renumbers nodes; `mesh.geometry.input_global_indices` maps them back to Gmsh node tags.

## 6. The checks against known answers

- **Cantilever, Euler-Bernoulli:** δ = PL³/(3EI). Bending only.
- **Timoshenko:** δ = PL³/(3EI) + PL/(kGA). Adds shear deformation; k is the shear correction factor, here Cowper's k = 10(1 + ν)/(12 + 11ν) = 0.850 for a rectangle at ν = 0.3. The shear term grows like (h/L)²: 0.76% of the Timoshenko deflection at L/h = 10, 4.6% at L/h = 4.
- **What "clamped" means.** In beam theory: no deflection, no rotation at the root. In the 3D model: every node of the root face is held, which also stops the steel contracting sideways (Poisson's ratio). That extra restraint made the 3D beams 0.71% and 1.71% stiffer than Timoshenko; with ν = 0 the gap fell to 0.004% and 0.06%.
- **Heywood's Kt** for a hole in a strip under tension: Kt = 2 + (1 - d/W)³ on the net-section stress P/((W - d)t). It is a plane-stress (2D) result. In a plate as thick as the hole is wide, the hoop stress is higher at mid-thickness and lower at the faces (8 mm bracket: +5.4% and -8%); a thin plate comes close to the 2D number.
- **Same mesh, two codes.** On straight-sided tets C3D10 and P2 are the same element, so with the same loads and constraints the answers must match to the digits printed. That is a stronger test of the toolchain (mesh export, node order, loads, BCs, parsing) than either code against theory.

## 7. What broke while building

1. The 3D cantilevers were stiffer than Timoshenko, and the slender one matched Euler-Bernoulli better. Not the mesh, not the solver: the clamp's Poisson restraint. ν = 0 proved it.
2. The hole Kt failed a 5% check at +5.4%. Plotting through the thickness and thinning the plate showed the 3D effect; the tolerance was replaced with checks of that pattern.
3. The first thin-plate mesh left 3 mm curved tets in the 2 mm fillet, and CalculiX and FEniCSx disagreed by 0.34%. Refining the fillet closed it.
4. The C3D4 runs first reported the C3D10 mesh's DOF count.
5. DOLFINx 0.10 API changes broke the first FEniCSx run (`Function.eval` returned a flat array).
6. Gmsh's `.inp` would have included the surface groups as 2D elements; only the volume group is kept in that file.
7. Meshes from the same inputs differ between Windows and Linux.

## 8. Explain it without notes

1. Why build the bracket two different ways, and why does their agreement not test the geometry kernel?
2. What is a B-rep, and what does a STEP file store?
3. What is a physical group in Gmsh, and why find faces by position rather than by tag?
4. How does a Distance + Threshold size field work, and what other size sources must be turned off?
5. What changes when you call `setOrder(2)`, and how would you check the curved elements are valid?
6. What are the consistent nodal forces for a uniform traction on a 6-node triangle, and why not one sixth per node?
7. Walk through a CalculiX deck card by card. What does `OP=NEW` do?
8. Where does CalculiX write the reaction forces, and where the stresses? Why does the `.frd` parser read fixed columns?
9. Write the weak form of linear elasticity. Which term carries the traction, and where does the Dirichlet condition go?
10. Why do C3D4 and P1 give the same (wrong) answer, and why do C3D10 and P2 give the same (right) one?
11. Why were the 3D cantilevers stiffer than Timoshenko, and what experiment showed it?
12. Why is the 8 mm bracket's Kt above Heywood's, and what happens as the plate gets thinner?
13. Why did the two solvers disagree on the first thin-plate mesh?
14. How does the CI make sure no test is silently skipped in both jobs?

## 9. Your turn

One exercise per tool, 20-40 minutes each. Each one changes something real, says what to expect (measured where I could measure it), and is what makes "I have worked with this tool" true. Write your prediction down before you run anything.

CadQuery, build123d and Gmsh run on the laptop:
```bash
cd cad-to-fea
D:/cache/venv-sim/Scripts/activate        # or: pip install -r requirements-ci.txt
export CADFEA_OUT=D:/cache/cadfea-out     # keep generated files off drive C
```

CalculiX and FEniCSx run only in CI. To run an experiment there, push a branch and read the log:
```bash
git switch -c try-something
# ...edit...
git commit -am "try: something"
git push -u origin try-something
gh run watch                                          # pick the run; about 9 minutes
gh run view --log | grep -A 25 "step 3, CalculiX"     # or "step 4, FEniCSx", or "FAILED"
gh run download -n results -D D:/cache/try-results    # summary.md and the plots
# afterwards: git switch main; git branch -D try-something; git push origin --delete try-something
```

### 1. CadQuery: change the fillet, then break the agreement on purpose (local, 25 min)
1. In `cadfea/params.py`, change the bracket's `r: float = 6.0` to `3.0`. Predict the new volume by hand: the fillet adds a spandrel of area r²(1 - π/4) over the 40 mm width, so the volume drops by 40 × 0.2146 × (36 - 9) = 231.770 mm³, to **67,288.938 mm³**.
2. Run `python scripts/build_geometry.py`. `results/geometry.txt` should show volume 67,288.938130 mm³ and the CadQuery vs build123d difference still around 1e-16. `python -m pytest tests/test_geometry.py` passes: both tools followed the parameter.
3. Put `r` back to 6.0. Now make CadQuery alone wrong: in `cadfea/geom_cq.py` change `.fillet(p.r)` to `.fillet(p.r * 1.001)`. Run the tests again: `test_cadquery_and_build123d_agree[bracket]` and `test_matches_closed_form[bracket]` fail, with a volume difference of about **9.2e-6**. Undo the change.
4. Explain to yourself: why does the build123d model not need to change in step 1, and which selector found the corner edge?

### 2. build123d: algebra mode, and a chamfer in both tools (local, 40 min)
1. Rewrite `bracket()` in `cadfea/geom_b123d.py` in algebra mode instead of builder mode:
   ```python
   from build123d import Pos  # add to the imports
   arm = Box(p.L, p.W, p.t, align=(MIN, CEN, MIN))
   upright = Box(p.t, p.W, p.H, align=(MIN, CEN, MIN))
   part = arm + upright
   corner = part.edges().filter_by(Axis.Y).sort_by_distance((p.t, 0, p.t))[0]
   part = fillet(corner, radius=p.r)
   part -= Pos(p.xh, 0, p.t / 2) * Cylinder(p.d / 2, p.t + 2)
   return part
   ```
   Expect `pytest tests/test_geometry.py` to pass unchanged (I measured 3.5e-16 against CadQuery).
2. Add a 1 mm × 45° chamfer to the top edge of the hole, in **both** tools. CadQuery: select the edge with `NearestToPointSelector((p.xh, 0, p.t))` (a circle's "centre" is its centre point) and `.chamfer(1.0)`. build123d: `chamfer(edge, length=1.0)`.
3. Trap: in build123d, `edges().filter_by(GeomType.CIRCLE).sort_by(Axis.Z)[-1]` does **not** pick the hole's edge; the fillet's end arcs sit higher. That chamfers the wrong edge and the volume drops by 101 mm³. Also filter by radius (`abs(e.radius - p.d / 2) < 1e-9`).
4. Expect both models to lose **16.755 mm³** (a ring of triangular section 0.5 mm² at mean radius 5.333 mm: 2π × 5.333 × 0.5) and to agree to about 4.5e-16, with 11 faces each. `test_matches_closed_form` will fail until you subtract the chamfer in `props.analytic` too: that is the point of an independent formula.

### 3. Gmsh: switch the refinement off, then double it (local; the Kt part in CI, 30 min)
1. Look at a mesh first:
   ```bash
   python -c "import gmsh; gmsh.initialize(); gmsh.open('D:/cache/cadfea-out/mesh/bracket_medium.msh'); gmsh.fltk.run()"
   ```
   Tools > Visibility lists the physical groups; hide `SOLID` to see the small triangles on the hole and fillet.
2. In `cadfea/params.py`, set the bracket's medium level to `MeshLevel("medium", 5.0, 5.0)` (no refinement). Run `python scripts/build_meshes.py`. On my laptop the medium bracket went from **27,251 nodes to 6,165**, the hole from 1,280 nodes to 74, and the curved-volume error from 3.2e-8 to 1.7e-5. `tests/test_mesh.py::test_size_field_refines_the_hole` now fails: that test is doing its job.
3. Now `MeshLevel("medium", 5.0, 0.5, fillet=1.0)`: expect about **56,885 nodes** and 4,866 on the hole (Linux counts differ slightly).
4. Push each version to a branch and read Kt for `bracket medium` in `results/summary.md`. With 0.5 mm on the hole, expect it within about 0.5% of the fine mesh's 2.552 (medium and fine already agree to 0.01%). With no refinement, predict first: the 2 mm hole mesh read 2.63-2.66, so 5 mm elements around a 10 mm hole will be further off, and you will see by how much.

### 4. CalculiX: forget OP=NEW, then go linear (CI, 35 min)
1. In `cadfea/ccx.py`, change `"*CLOAD" + (", OP=NEW" if i else "")` to `"*CLOAD"`. Before pushing, predict the bracket's second step: without `OP=NEW` the pull step carries the shelf load too, and the model is linear, so its tip displacement is the sum of the two load cases. On the fine mesh: **ux = 0.0239958 mm, uz = -0.642666 mm** (instead of 0.0244522 and 0.0228137).
2. Push. In the "step 3, CalculiX" log, the `bracket fine ... pull:` numbers should match your prediction to about 6 digits, and `test_reactions_balance_the_load[bracket]` fails: the pull step's reaction now has +200 N in z. The reaction check is what catches a load leaking from one step into the next. Revert.
3. Linear elements on the bracket: add a CI line `python scripts/run_ccx.py --cases bracket --element C3D4`. On the beams, C3D4 was 51%, 17% and 5.5% too stiff (coarse, medium, fine). Predict the bracket's shelf deflection against the C3D10 value of -0.6655 mm before you look; the arm is a plate in bending, so expect the same kind of error.
4. Write a deck locally (no ccx needed) and read it card by card:
   `python -c "from cadfea import pipeline, ccx, paths; i = pipeline.ensure_mesh('bracket', 'coarse'); print(ccx.write_deck(i, paths.MESH_DIR, paths.CCX_DIR)[0])"`

### 5. FEniCSx: P1, then ν = 0 in your own code (CI, 40 min)
1. `scripts/run_fenicsx.py --degree 1` solves with linear elements. Already in CI for the beams (`--with-p1`): -51.3%, -16.8%, -5.5% on the slender beam, identical to CalculiX C3D4. Explain why identical.
2. Make FEniCSx run the ν = 0 experiment, which so far only CalculiX runs. Give `pipeline.ensure_fenicsx` a `mat: Material = STEEL` argument, pass it to `fenics.solve(..., mat=mat)` and into the result path (as `ensure_ccx` does), add `--nu0` to `scripts/run_fenicsx.py`, and a CI line `python3 scripts/run_fenicsx.py --cases beam,stocky --levels fine --nu0`.
3. Expect the fine beam at **0.958054 mm** and the stocky one at **0.063198 mm**: the CalculiX ν = 0 values, to about 7 digits, because it is the same mesh and the same element. Against Timoshenko at ν = 0: -0.004% and -0.06%.
4. Optional: `--solver cg` (CG with algebraic multigrid). The test `test_cg_with_multigrid_matches_lu` already shows it matches LU to 1e-8 on the medium beam; compare the iteration count and time on the fine bracket. Then delete the `_set_near_nullspace` call and predict what happens to the iteration count (it may not converge within the 2,000-iteration limit at all).
