# Architecture Overview — dolfin-hpc 0.9.x

**Target audience:** New MSc students with FEM background but no FEniCS experience.  
**Goal:** Understand how a variational form written in UFL becomes an assembled PETSc matrix.

---

## 1. The four repos and their roles

**ufl-hpc** is the form language. You write your PDE as a variational problem in a
Python-based domain-specific language — elements, trial/test functions, integrals — and save it
in a `.ufl` file. ufl-hpc is a fork of UFL (Unified Form Language) tuned for the HPC stack.
It defines the symbolic objects (`FiniteElement`, `TrialFunction`, `inner`, `grad`, `dx`, etc.)
and produces a symbolic representation of the form. Nothing numerical happens here; ufl-hpc is a
pure math layer. Everything downstream depends on it through ffc-hpc.

**ffc-hpc** is the form compiler. It takes a `.ufl` file, applies symbolic differentiation and
optimisation (quadrature point selection, tensor contraction optimisation), and emits a
self-contained C++ header. You invoke it as `ffc -l dolfin MyForm.ufl`. The output header
contains all the element-level mathematics, pre-compiled into static C++ functions. ffc-hpc
depends on ufl-hpc to parse the form language and on ufc2-hpc to know which interface the
generated code must implement. The generated header is the only ffc-hpc artefact that ends up in
your solver build.

**ufc2-hpc** is the interface contract. It defines a small set of abstract C++ base classes —
`ufc::finite_element`, `ufc::dof_map`, `ufc::form`, `ufc::cell_integral`,
`ufc::exterior_facet_integral` — in a header-only cmake install. ffc-hpc generates classes
that *implement* these interfaces; dolfin-hpc's assembler calls through them at runtime without
knowing which form is being assembled. ufc2-hpc is the glue that decouples the compiler output
from the runtime library. Both ffc-hpc (as a code-generation target) and dolfin-hpc (as a
consumer) depend on it.

**dolfin-hpc** is the runtime library. It reads meshes, manages degree-of-freedom maps, drives
the assembler loop over mesh cells, applies boundary conditions, and wraps PETSc for distributed
linear algebra. A user solver links against dolfin-hpc and includes a generated form header;
dolfin-hpc provides the `Mesh`, `Function`, `Assembler`, `KrylovSolver`, and all supporting
infrastructure. It depends on ufc2-hpc (to call through the `ufc::` interfaces), PETSc (linear
algebra), and ParMETIS (mesh partitioning).

---

## 2. Data flow: from .ufl to PETSc Mat

The path from a variational form to an assembled sparse matrix has five stages.

```
  You write                ffc-hpc compiles         dolfin-hpc runs
  ──────────               ───────────────          ───────────────

  Poisson.ufl              Poisson.h                main.cpp
  ─────────                ─────────                ───────
  element = FiniteElement  class poisson_           Mesh mesh("...bin");
    ("Lagrange","tri",1)     cell_integral_0        Poisson::BilinearForm a(mesh);
  v = TestFunction(elem)     : ufc::cell_integral   Poisson::LinearForm   L(mesh,f,g);
  u = TrialFunction(elem)  {                        Matrix A;  Vector b;
  f = Coefficient(elem)      tabulate_tensor(       a.assemble(A, true);
  a = inner(grad(v),           A, w,                L.assemble(b, true);
        grad(u))*dx            coord_dofs,          bc.apply(A, b, a);
  L = v*f*dx + v*g*ds          orientation)         KrylovSolver s(...);
                             { /* Gaussian          s.solve(A, u.vector(), b);
                               quadrature,
                               Jacobian,
                               basis fns */
                             }
                           };
                           class BilinearForm
                             : dolfin::BilinearForm
                           { /* wraps above */ };

  STEP 1          STEP 2            STEP 3           STEP 4         STEP 5
  .ufl source  →  ffc-hpc       →   Form object  →   Assembler  →   PETSc Mat
  (symbolic)      generates         (holds ufc::     loops cells,    (distributed,
                  Poisson.h         form ref +        calls           ready for
                                    DofMap)           tabulate_       solve)
                                                      tensor()
```

**Step 1 — Write the form.** `Poisson.ufl` (12 lines) expresses the bilinear form
`a(v,u) = ∫ ∇v·∇u dx` and linear form `L(v) = ∫ vf dx + ∫ vg ds` symbolically.

**Step 2 — Compile with ffc.** Running `ffc -l dolfin Poisson.ufl` produces `Poisson.h`
(~3500 lines). All quadrature points, basis function tabulations, and Jacobian transformations
are baked into static C++ at this step. The generated file includes `<ufc.h>` and the
DOLFIN-HPC wrappers at its base.

**Step 3 — Construct Form objects.** In `main.cpp`, `Poisson::BilinearForm a(mesh)` creates an
object that holds a reference to the mesh, instantiates the generated `ufc::finite_element` and
`ufc::dof_map` objects, and builds `FiniteElementSpace` and `DofMap` tables for the trial and
test spaces. `Form` (`include/dolfin/fem/Form.h`) is the abstract base; the generated subclass
provides `form()`, `coefficients()`, and `cell_integrals()`.

**Step 4 — Assemble.** `a.assemble(A, true)` delegates to `Assembler::assemble(A, form, reset)`
(`include/dolfin/fem/Assembler.h`). The assembler loops over all local cells, fetches geometric
data and coefficient values, calls `ufc::cell_integral::tabulate_tensor(A, w, coordinate_dofs, cell_orientation)` for
each cell to compute the local element stiffness matrix, then scatters `A_local` into the global
`GenericTensor` using the DOF map.

**Step 5 — Distributed linear algebra.** `Matrix A` is backed by `PETScMatrix`
(`include/dolfin/la/petsc/PETScMatrix.h`), which wraps PETSc's `Mat`. Each MPI process
contributes its local rows; PETSc handles the parallel sparse format. After assembly the matrix
is ready for `KrylovSolver` or `LUSolver`.

---

## 3. Key abstractions in dolfin-hpc

**`Mesh`** (`include/dolfin/mesh/Mesh.h`) — stores the computational domain as a combination of
`MeshTopology` (cell-vertex-edge connectivity tables) and `MeshGeometry` (vertex coordinates);
loaded from `.bin` binary format via `Mesh mesh("file.bin")`.

**`MeshFunction<T>`** (`include/dolfin/mesh/MeshFunction.h`) — a field of values of type `T`
indexed by mesh entities (vertices, edges, facets, or cells); used to tag boundary regions and
mark subdomains for selective assembly or boundary condition application.

**`FiniteElementSpace`** (`include/dolfin/fem/FiniteElementSpace.h`) — pairs a `Mesh` with a
`FiniteElement` and a `DofMap`, providing the map from local cell DOFs to the global DOF index
vector; created automatically when a `Form` is constructed.

**`Function`** (`include/dolfin/function/Function.h`) — a finite element function: a
`FiniteElementSpace` plus a `Vector` of DOF coefficients; in `main.cpp`,
`Function u(a.trial_space())` allocates the solution vector in the trial space of form `a`.

**`Form`** (`include/dolfin/fem/Form.h`) — abstract base wrapping a `ufc::form const&`; holds
`cell_integrals_`, `exterior_facet_integrals_`, and vectors of `FiniteElementSpace` and `DofMap`
objects; subclassed by the generated `BilinearForm` and `LinearForm` in each `MyForm.h`.

**`Assembler`** (`include/dolfin/fem/Assembler.h`) — a namespace (not a class) exposing
`assemble(GenericTensor&, Form&, bool)` and subdomain-restricted variants; the core loop that
drives `tabulate_tensor` over all cells and scatters into the global tensor.

**`BoundaryCondition` / `DirichletBC`** (`include/dolfin/fem/DirichletBC.h`) — applies essential
boundary conditions; `bc.apply(A, b, a)` modifies assembled rows of `A` and entries of `b` to
enforce the constraint, using a `SubDomain::inside()` predicate to identify DOFs on the boundary.

---

## 4. Where parallelism enters

**Mesh partitioning.** When a mesh is loaded in parallel (N MPI processes), `MeshPartition`
(using ParMETIS via `LoadBalancer`) assigns each cell to one process. Each process then holds its
local partition plus a one-cell-thick **ghost layer** — copies of neighbouring cells owned by
adjacent processes. Ghost cells are needed so that assembly near partition boundaries can access
the full element geometry and coefficient values.

**Distributed assembly.** `Assembler::assemble` only iterates over locally owned cells. Because
ghost cells are present, all integrals contributing to locally owned DOFs are computed correctly
without communication during the assembly loop itself. The `DofMap` maps local DOF indices to
global indices in the distributed numbering.

**Distributed linear algebra.** `PETScMatrix` and `PETScVector` wrap PETSc's parallel `Mat` and
`Vec`. Each process owns a contiguous block of global rows (determined by `DofMap`). After
assembly, PETSc's standard solvers (`KrylovSolver` with `bjacobi`, `PETScLUSolver` with MUMPS)
operate on the distributed system without further changes to user code.

**Key implication for solvers.** The in-built mesh generators (`UnitSquare`, `UnitCube`) do not
work in parallel — meshes must be generated in serial and loaded in parallel from `.bin` files.
See `KNOWN_ISSUES.md` §4.

---

## 5. What ffc generates — inside Poisson.h

The generated header `Poisson.h` (produced by `ffc -l dolfin Poisson.ufl`, FFC 1.2.0, UFC 2.4.0)
contains three layers.

**Element classes.** `poisson_finite_element_0 : public ufc::finite_element` implements the
Lagrange P1 basis: `space_dimension()` returns 3 (triangle), `evaluate_basis()` evaluates φ_i(x),
and `evaluate_basis_derivatives()` evaluates ∂φ_i/∂x_j. These methods are compiled from symbolic
expressions derived from the UFL element definition.

**DOF map class.** `poisson_dof_map_0 : public ufc::dof_map` implements `tabulate_dofs()`, which
maps local DOF indices (0, 1, 2 for a triangle) to global DOF indices given the mesh topology.
The assembler calls this for each cell to know where to scatter the local tensor.

**Cell integral class.** `poisson_cell_integral_0 : public ufc::cell_integral` implements
`tabulate_tensor(double* A, const double* const* w, const double* coordinate_dofs, int cell_orientation)` (UFC 2.4.0 signature). This single function
is the computational core: it evaluates the bilinear form `a(v,u)` on one cell by Gaussian
quadrature, computing the Jacobian of the geometric mapping, pulling back the basis function
gradients, and accumulating the 3×3 local stiffness matrix into `A`. All constants and
quadrature weights are hard-coded; no symbolic evaluation happens at runtime.

**DOLFIN-HPC wrappers.** At the end of the file, `Poisson::BilinearForm : public
dolfin::BilinearForm` and `Poisson::LinearForm : public dolfin::LinearForm` register the above
classes with the dolfin-hpc runtime, expose `trial_space()` and `test_space()`, and wire up
named coefficients (`f`, `g`) to their `Coefficient` objects.

**Why this design?** Separating symbolic compilation (ffc-hpc, at build time) from numerical
execution (dolfin-hpc, at run time) means that assembly is as fast as hand-written C++. The
assembler is a generic loop; all PDE-specific mathematics lives in `tabulate_tensor`. Adding a
new PDE requires writing a `.ufl` file and running `ffc` — dolfin-hpc does not change.
