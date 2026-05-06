# Known Issues

This file documents known bugs and limitations in dolfin-hpc 0.9.x.
Where a fix is not yet available, a workaround is described.
Issues are tracked here until a proper fix is committed.

---

## 1. Memory leaks / pointer allocation

**Status:** Fix in progress (patch expected from Joel Kronborg)  
**Affected versions:** 0.9.x  
**Description:** Memory leaks have been identified in pointer allocation in the core library. The exact locations are being identified via Valgrind4hpc on Dardel.  
**Workaround:** None required for correctness; affects long-running simulations and memory-constrained environments.  
**Resolution:** Patch pending — will be committed to `dev` once received and reviewed.

---

## 2. Interpolation operator

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

## 3. Periodic boundary conditions

**Status:** Known bug, no fix available  
**Affected versions:** 0.9.x  
**Description:** Periodic boundary conditions do not work correctly in all cases. The failure mode is not fully characterised.  
**Workaround:** Avoid periodic BCs where possible. Use alternative boundary condition formulations (e.g. Dirichlet or Neumann) that achieve an equivalent problem setup.  
**Resolution:** Under investigation. Contributions welcome.

---

## 4. In-built mesh generators in parallel

**Status:** Known limitation  
**Affected versions:** 0.9.x  
**Description:** In-built mesh generators do not work correctly in parallel computing environments.  
**Workaround:** Generate meshes in serial (e.g. using an external mesh generator), then load the mesh in parallel. Use `dolfin-convert` to convert from `.xml` to `.bin` format as required by 0.9.x.  
**Resolution:** No fix planned in the near term.

---

## 5. MeshQualityFunction linking error

**Status:** Known build issue on Dardel  
**Affected versions:** 0.9.x  
**Description:** A linking error involving `MeshQualityFunction` has been observed when building on Dardel with certain module environments. The class is defined as a header-only implementation in `include/dolfin/mesh/utilities/MeshQualityFunction.h`.  
**Workaround:** If the linking error occurs, remove the include of `MeshQualityFunction.h` from `include/dolfin/mesh/dolfin_mesh.h` and the corresponding entry from `include/dolfin/Makefile.am`.  
**Resolution:** Under investigation. Will be revisited when Dardel build environment is fully documented.

---

## Reporting new issues

Please open a GitHub issue at https://github.com/dolfin-forge/dolfin-hpc/issues with:
- A minimal reproducer (mesh, form file, driver code)
- The Dardel module environment (`module list` output)
- The error message or unexpected behaviour observed
