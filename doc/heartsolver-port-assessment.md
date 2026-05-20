# HeartSolver Port Assessment: dolfin-hpc 0.8 → 0.9

**Date:** 2026-05-07  
**Source assessed:** `fenicsheartsolver-heartsolver-f3923f758c3d` (porous branch snapshot)  
**Reference headers:** `~/repos/dolfin-hpc/include/` (dolfin-hpc 0.9, dev branch)

---

## Architecture of the current codebase

The HeartSolver snapshot is largely self-contained. It embeds the unicorn-hpc 0.1.2-hpc FSI
infrastructure directly in its own `src/` and `include/` tree — there is no separate `.so`
dependency to install. The full dependency surface is:

```
HeartSolver (this repo)
  ├── dolfin-hpc 0.8          (external, pkg-config)
  ├── UFC 1.0                 (external, pkg-config as ufc-1)
  └── [embedded unicorn 0.1.2-hpc]
        TimeDependentPDE, ALESolver, ElasticSmoother,
        LaplacianSmoother, ErrorEstimate, NodeNormal,
        unicorn_init, unicorn_solve, Checkpoint (0.8 API)
```

The `.h` files in `include/` that appear in `.py`/`.form`/`.h` triplets are FFC 0.5.1 + UFC 1.0
generated form code. Every one of them must be regenerated for 0.9.

---

## Blocking issues (must resolve before compilation)

### 1. UFC 1.0 → ufc2-hpc (form recompilation)

All 14 generated form headers (`NavierStokes3D`, `NavierStokesPressure3D`, `BF`,
`Elasticity2D/3D`, `ElasticityJac2D/3D`, `ElasticityStress2D/3D`, `Laplacian2D/3D`,
`Projection2D/3D`, `MeanPressure`) conform to `ufc::finite_element` / `ufc::form` from UFC 1.0.
ufc2-hpc has a different interface.

The `.py` source files use `ffc 0.5.1` syntax: `from ffc import *`,
`compile([...], ..., {'language': 'dolfin', ...})`. This must be ported to ffc-hpc 0.9's
UFL-based pipeline. The form language changed significantly between the two eras — old-style FFC
used tensor-product notation, while the 0.9 pipeline uses UFL expressions.

The `Makefile` links against `pkg-config --cflags ufc-1`. This needs to change to the ufc2-hpc
cmake install.

### 2. Checkpoint API — complete redesign

The 0.8 Checkpoint is used with a procedural, state-mutating interface:

```cpp
// 0.8 pattern (used in init.cpp, ALESolver.cpp, main.cpp):
chkp.restart(path);     // set restart path
chkp.load(mesh);        // load mesh from checkpoint
chkp.restart()          // returns bool: are we restarting?
chkp.restart_time()     // returns real: restart time
```

The 0.9 Checkpoint has a completely different, explicit map-based interface:

```cpp
// 0.9 pattern:
Checkpoint::MeshMap meshes;
Checkpoint::FunctionMap funcs;
chkp.load_header(filename);
chkp.load(filename, meshes);  // explicit filename each time
chkp.time();                  // replaces restart_time()
```

There is no `restart()` bool method. The entire checkpoint/restart logic in `ALESolver.cpp` and
`main.cpp` (the conditional branch at `main.cpp:142`) must be redesigned.

### 3. BoundaryMesh API

The 0.8 patterns used pervasively throughout `main.cpp` and ALESolver:

```cpp
// 0.8:
BoundaryMesh boundary_mesh;
boundary_mesh.init(mesh);                          // default ctor + init
MeshFunction<uint>* vmap = bmesh.data().meshFunction("vertex map");
MeshFunction<uint>* cmap = bmesh.data().meshFunction("cell map");
vmap->get(vertex->index())                         // indirect lookup
```

The 0.9 BoundaryMesh requires a type at construction and exposes direct index lookup:

```cpp
// 0.9:
BoundaryMesh bmesh(mesh, BoundaryMesh::exterior);  // type required
bmesh.vertex_index(boundary_vertex_index)           // direct
bmesh.facet_index(boundary_cell_index)              // direct
// MeshData/meshFunction("vertex map") does not exist
```

This pattern appears in `main.cpp:206–276` and throughout ALESolver.

---

## High-priority API changes (significant but mechanical)

### 4. Function initialization

The 0.8 two-step `init(mesh, vec, form, index)` is used for every working field:

```cpp
// 0.8:
Function boolbc;
Vector mf;
boolbc.init(mesh, mf, *fform, 0);
```

The 0.9 `Function` is initialized via a `FiniteElementSpace` or `Form`:

```cpp
// 0.9:
Function boolbc;
boolbc.init(*fform, 0);  // mesh is implicit in form; vector managed internally
```

Every `Function` in `main.cpp` and the entire ALESolver uses the 0.8 pattern — roughly 30+ call
sites.

### 5. DirichletBC with MeshFunction markers

The 0.8 code uses integer-marker MeshFunctions:

```cpp
// 0.8:
DirichletBC bc_mom_sys_wall(bcf_mom_wall, sub_domains_facets, WALL_BOUNDARY);
```

The 0.9 `DirichletBC` takes a `SubDomain const&`, not a `MeshFunction<uint>` + integer. The
entire subdomain-marking approach must be replaced by `SubDomain` subclasses that implement
`inside()` for each physical boundary.

This is a design-level change — the `SurfaceDistribution` class that populates `sub_domains_facets`
via integer tags must either be kept as a mesh-marking utility (and the BCs reconstructed from it
using helper SubDomains) or replaced entirely.

**Additional: time-dependent BCs — `update()` does not exist in 0.9.**

The 0.8 HeartSolver required a custom `update()` method added manually to `DirichletBC.h`. That
method was never upstreamed and is absent from 0.9. The correct 0.9 replacement is the
`bc(t)` pattern:

1. Wrap each time-varying boundary value in a `TimeDependent` subclass
   (`include/dolfin/evolution/TimeDependent.h`) that overrides `sync(Time const& t)` to update
   its value from `clock()`.
2. Pass an instance to `DirichletBC` as the `Coefficient` argument.
3. In the time loop, call `bc(t)` before `bc.apply(A, b)`. This propagates time through
   `BoundaryCondition::operator()(Time const& t)` → `DirichletBC::sync(t)` → each stored
   `Coefficient::operator()(t)` (see `include/dolfin/fem/BoundaryCondition.h`).

There is no demo in the shipped codebase that demonstrates this — the `demo/pde/bcs/` demo uses
only static coefficients. Every time-dependent BC call site in `main.cpp` and `ALESolver.cpp`
must be converted; this work is included in the 1.0-week estimate for this item.

### 6. MPI API rename

| 0.8 call | 0.9 equivalent |
|----------|----------------|
| `dolfin::MPI::processNumber()` | `MPI::rank()` |
| `dolfin::MPI::numProcesses()` | `MPI::size()` |
| `dolfin::MPI::DOLFIN_COMM` | `MPI::DOLFIN_COMM` (still a static member) |

`MPI_Barrier(dolfin::MPI::DOLFIN_COMM)` and `MPI_Allreduce(...)` calls in `main.cpp` are fine —
they use raw MPI, just need to extract the communicator correctly.

### 7. MeshDistributedData API

The 0.8 pattern is used heavily in `main.cpp` vertex loops:

```cpp
// 0.8:
mesh.distdata().get_global(vertex->index(), 0)
mesh.distdata().is_ghost(vertex.index(), 0)
mesh.distdata().global_numVertices()
mesh.distdata().global_numCells()
```

The 0.9 `MeshDistributedData` is accessed as `mesh.distdata()[dim]` returning a `DistributedData&`
— a completely different structure. The `get_global`, `is_ghost` methods need to be located in the
0.9 `DistributedData` interface.

---

## Medium-priority changes (isolated, local)

### 8. NodeNormal

- 0.8: `NodeNormal nn(mesh)` — extends nothing
- 0.9: `NodeNormal(Mesh& mesh, Type, real alpha)` — extends `BoundaryNormal`

Since HeartSolver has its own `NodeNormal.cpp`, the implementation itself must be ported to the 0.9
internal APIs (mesh iteration, dof access). The interface change is modest; the implementation
churn may be larger.

### 9. dolfin_set / dolfin_get

`dolfin_set`/`dolfin_get` exist in 0.9 as templates, so the call syntax works. The parameter keys
must be audited — keys like `"output destination"`, `"Mesh read in serial"`,
`"PDE reassemble matrix"`, `"beta"`, `"adapt_tol"`, `"output_format"`, `"solution file name"` must
exist in the 0.9 parameter system or be added.

### 10. `dolfin/fem/UFC.h`

`ALESolver.cpp` and `ALESolver.h` include `<dolfin/fem/UFC.h>`. This header does not appear in the
0.9 include tree — the UFC assembly helper is replaced by the ufc2-hpc mechanism. Every usage of
the `UFC` class in ALESolver must be found and replaced.

### 11. `Array<BoundaryCondition*>` container

0.8 uses a custom `Array<T>` throughout. If 0.9 has dropped this in favour of `std::vector`, all
container types need updating. `Array<BoundaryCondition*>` is used as the BC storage type in
`TimeDependentPDE`, `ALESolver`, `Heart_Solver`.

---

## What maps directly (minimal change)

| Component | Status |
|-----------|--------|
| `Assembler` | API stable, minor changes |
| `KrylovSolver` / `LUSolver` | API stable |
| `Matrix` / `Vector` | API stable |
| `File` I/O for `.bin` | Compatible in 0.9 |
| `MeshFunction<T>` basic usage | Stable |
| `CellIterator`, `VertexIterator`, `FacetIterator` | Stable |
| `SubDomain` | Stable |
| `dolfin_set` / `dolfin_get` | Present in 0.9 (template form) |
| `message()` | Likely stable |
| `dolfin_init` / `dolfin_finalize` | Check parameter count |

The cardiac physics — the NS formulations, the ALE kinematics, the leaflet motion data
interpolation — are not affected. The `heart.h` parameter file, `leaflet_motion_*.h` data tables,
and `MitralValve`/`Leaflet` geometry classes require no dolfin-hpc API changes.

---

## What is absent from dolfin-hpc 0.9

| Missing item | Notes |
|---|---|
| `dolfin/fem/UFC.h` | No direct equivalent; removed with UFC 1.0 |
| `Checkpoint` 0.8 API (`restart()`, `load(mesh)`) | Replaced by explicit map-based API |
| `BoundaryMesh::data().meshFunction("vertex map")` | Replaced by `vertex_index()`/`facet_index()` |
| `MPI::processNumber()` / `numProcesses()` | Renamed to `rank()` / `size()` |
| UFC 1.0 `ufc::` base classes | Replaced by ufc2-hpc |

---

## Build system

The current `Makefile` uses `pkg-config --cflags ufc-1`. The 0.9 build requires:

- ufc2-hpc (cmake install) — must be built first
- dolfin-hpc 0.9 configured with `--with-petsc=...` as documented in CLAUDE.md
- The Makefile or a replacement CMakeLists.txt must find ufc2-hpc headers and link against
  dolfin-hpc 0.9

---

## Effort estimate (focused MSc student)

| Task | Est. weeks |
|---|---|
| Environment setup, build system migration, codebase orientation | 0.5 |
| Form recompilation: port 14 `.py` FFC 0.5.1 forms to ffc-hpc 0.9 / UFL | 2.0 |
| Checkpoint API redesign throughout ALESolver + main | 1.0 |
| BoundaryMesh API + MeshDistributedData update | 1.0 |
| Function init API update (~30 call sites) | 0.5 |
| DirichletBC marker → SubDomain refactor | 1.0 |
| MPI rename + UFC.h removal + Array→vector | 0.5 |
| NodeNormal port to 0.9 internals | 0.5 |
| ElasticSmoother / LaplacianSmoother port | 0.5 |
| TimeDependentPDE / ALESolver port | 1.5 |
| Mesh format conversion (xml → bin) + parameter key audit | 0.5 |
| First build + link errors | 0.5 |
| Functional testing (compare against 0.8 on reference case) | 2.0 |
| **Total** | **~12 weeks** |

**Confidence interval: 10–16 weeks.** The lower bound assumes the form recompilation is
straightforward (the FFC 0.5.1 `.py` sources are clean and the UFL equivalent is direct). The
upper bound accounts for the DirichletBC subdomain redesign revealing unforeseen assumptions about
how `sub_domains_facets` is built, and for test-case infrastructure to be set up from scratch.

---

## Recommended entry point

The logical sequencing:

1. Get `ufc2-hpc` installed and `ffc-hpc 0.9` working on Dardel — recompile the `.py` forms first;
   compilation failures here expose any form-level incompatibilities early.
2. Port `TimeDependentPDE.cpp` — it is the base of the whole class hierarchy and the simplest
   piece.
3. Work up: `ALESolver` → `Heart_Solver` → `main.cpp`.
4. Checkpoint and DirichletBC redesign last, since they require the most design judgment.

Note: the `KNOWN_ISSUES.md` entry about `MeshQualityFunction` linking errors on Dardel is directly
relevant — `MeshQuality qual(mesh)` appears in `ALESolver.h` and that same module may surface the
same linker issue during the 0.9 port.
