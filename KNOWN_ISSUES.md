# Known Issues

This file documents known bugs and limitations in dolfin-hpc 0.9.x.
Where a fix is not yet available, a workaround is described.
Issues are tracked here until a proper fix is committed.

---

## 1. Interpolation operator

**Status:** Known bug, no fix available  
**Affected versions:** 0.8.x, 0.9.x  
**Description:** The interpolation operator `<<` does not work correctly for all meshes and polynomial degrees. The precise conditions under which it fails are not fully characterised.  
**Workaround:** Use projection instead of interpolation. Replace expressions of the form:

```cpp
u << expr;
```

with an explicit projection:

```cpp
project(expr, V, u);
```

**Resolution:** Under investigation. Contributions welcome.

---

## 2. Periodic boundary conditions

**Status:** Known bug, no fix available  
**Affected versions:** 0.9.x  
**Description:** Periodic boundary conditions do not work correctly in all cases. The failure mode is not fully characterised.  
**Workaround:** Avoid periodic BCs where possible. Use alternative boundary condition formulations (e.g. Dirichlet or Neumann) that achieve an equivalent problem setup.  
**Resolution:** Under investigation. Contributions welcome.

---

## 3. In-built mesh generators in parallel

**Status:** Known limitation  
**Affected versions:** 0.9.x  
**Description:** In-built mesh generators do not work correctly in parallel computing environments.  
**Workaround:** Generate meshes in serial (e.g. using an external mesh generator), then load the mesh in parallel. Use `dolfin-convert` to convert from `.xml` to `.bin` format as required by 0.9.x.  
**Resolution:** No fix planned in the near term.

---

## 4. MeshQualityFunction linking error

**Status:** Known build issue on Dardel  
**Affected versions:** 0.9.x  
**Description:** A linking error involving `MeshQualityFunction` has been observed when building on Dardel with certain module environments. The class is defined as a header-only implementation in `include/dolfin/mesh/utilities/MeshQualityFunction.h`.  
**Workaround:** If the linking error occurs, remove the include of `MeshQualityFunction.h` from `include/dolfin/mesh/dolfin_mesh.h` and the corresponding entry from `include/dolfin/Makefile.am`.  
**Resolution:** Under investigation. Will be revisited when Dardel build environment is fully documented.

---

## 5. dolfin-convert requires Python 2

**Status:** Known limitation  
**Affected versions:** 0.9.x  
**Description:** `dolfin-convert` (in `misc/utils/convert/attic/`) is used to convert Gmsh `.msh`
files to the dolfin-hpc `.bin` format required by 0.9.x. The script requires Python 2, which is
unavailable on Dardel (Python 3.11 only) and on most modern local machines. This means new mesh
generation from Gmsh is currently blocked on these systems.  
**Workaround:** Use a system with Python 2 installed to run the conversion, then transfer the
resulting `.bin` file. Alternatively, use the pre-converted `.bin` meshes provided with the
installed demos (e.g. `UnitSquareMesh_32x32.bin` in the Poisson demo).  
**Resolution:** Port `dolfin-convert` to Python 3, or replace with
[meshio](https://github.com/nschloe/meshio) (a modern Python 3 mesh conversion library that
supports both Gmsh `.msh` and DOLFIN XML output).

---

## 6. --enable-optimize-p1 blocks Taylor-Hood demos (Stokes)

**Status:** Known limitation of the installed Dardel module  
**Affected versions:** dolfin-hpc/0.9.5 module on Dardel (software environment 24.11)  
**Description:** The installed module was compiled with `--enable-optimize-p1`. This flag
enables a P1-specific optimisation path that is incompatible with higher-order elements.
Any demo or solver using Taylor-Hood P2/P1 elements (e.g. the shipped Stokes demo) will
abort immediately with:

```
This demo cannot be run with p1 optimizations enabled!
```

The Poisson demo (pure P1) is unaffected and runs correctly.

**Note:** `--enable-optimize-p1` is the correct build choice for HeartSolver, which uses
P1/P1 elements with stabilization throughout. No module rebuild is needed for HeartSolver
use. The limitation applies only to demos and solvers that use higher-order elements
(P2/P1 Taylor-Hood or above).  
**Workaround:** None without rebuilding. Do not attempt the Stokes demo or any P2/P1
solver with the current module.  
**Resolution:** Rebuild the module without `--enable-optimize-p1`, or request an updated
module from PDC. When requesting, specify that P2/P1 (Taylor-Hood) support is required.
This is only needed if running solvers that use higher-order elements; HeartSolver is unaffected.

---

## Reporting new issues

Please open a GitHub issue at https://github.com/dolfin-forge/dolfin-hpc/issues with:
- A minimal reproducer (mesh, form file, driver code)
- The Dardel module environment (`module list` output)
- The error message or unexpected behaviour observed
