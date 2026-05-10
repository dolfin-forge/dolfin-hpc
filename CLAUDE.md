# Project Context for Claude Code

This file provides context for Claude Code sessions working on the dolfin-forge project.
Read this before doing anything else.

---

## Who I am

Johan Hoffman, Professor of Numerical Analysis at KTH Royal Institute of Technology, Stockholm.
Co-founder of FEniCS and one of the two original DOLFIN developers.
Research focus: Direct FEM Simulation (DFS), adaptive FEM, turbulent flows, Navier-Stokes.

I have deep domain knowledge — no need to explain FEM basics. Be direct and technical.
Flag anything that may break existing HPC build configurations.
Prefer conservative, well-justified changes over clever modernization.

---

## The codebase

All repos live under the GitHub organization: https://github.com/dolfin-forge

### Core library stack (dolfin-hpc 0.9.x — active development target)
- `dolfin-hpc` — core C++ FEM library, HPC-optimized fork of DOLFIN
- `ufl-hpc` — Unified Form Language
- `ffc-hpc` — FEniCS Form Compiler  
- `ufc2-hpc` — Unified Form-assembly Code interface

### HeartSolver stack (currently runs on dolfin-hpc 0.8 only)
- `heartsolver` — main heart solver (FSI, cardiac mechanics)
- `unicorn-hpc-heartsolver` — modified unicorn-hpc v0.1.2 FSI infrastructure
- `heart-valves` — FSI solver for valve geometries (not yet compatible with heartsolver)
- `heartsolver-scripts` — visualization/patient scripts (private, IP status unclear)

### Local clones on this machine
- `~/repos/dolfin-hpc` — working clone of dolfin-hpc (this is the active dev repo)
- `~/repos/heartsolver_repos/` — flat snapshots of HeartSolver repos downloaded from Bitbucket

---

## Current state of dolfin-hpc

### Branches
- `master` — protected, stable, mirrors original Bitbucket codebase + migration fixes
- `dev` — integration branch, two commits ahead of master:
  1. `fix: Dardel build patches for modern PETSc and C++ compatibility`
  2. `docs: add KNOWN_ISSUES.md`
- Feature branches follow `fix/`, `feature/`, `docs/` prefix convention
- All work goes to `dev` via PR, never directly to `master`

### Known patches already applied to dev
- `m4/ax_petsc.m4`: `variables` → `petscvariables` (4 occurrences, PETSc 3.19+ compatibility)
- `include/dolfin/fem/UFCCell.h`: added `#include <utility>`
- `src/la/petsc/PETScMatrix.cpp`: `PETSC_NULL` → `NULL`
- `src/la/petsc/PETScLUSolver.cpp`: `PETSC_NULL` → `NULL`
- `src/main/SubSystemsManager.cpp`: `PETSC_NULL` → `NULL`
- `include/dolfin/la/petsc/PETScVector.h`: `PETSC_NULL` → `NULL`

### Known issues (see KNOWN_ISSUES.md in dolfin-hpc)
- Memory leaks in pointer allocation (patch expected from Joel Kronborg)
- Interpolation operator `<<` broken for some meshes/degrees — workaround: use projection
- Periodic BCs broken in some cases — workaround: use alternative BC formulations
- In-built mesh generators don't work in parallel
- MeshQualityFunction linking error on Dardel in some module environments

### Build environment (Dardel, HPE Cray EX)
- Currently installed as module: `ml PDC; ml dolfin-hpc/0.9.5`
- PETSc 3.21.1 at `/pdc/software/23.12/other/petsc/3.21.1/petsc-3.21.1/`
- Configure command:
  ```
  ./configure CC=cc CXX=CC MPICXX=CC CXXFLAGS="-std=c++14 -g" \
    --with-parmetis --with-petsc=/pdc/software/23.12/other/petsc/3.21.1/petsc-3.21.1/ \
    --with-gts --host=x86_64-unknown-linux-gnu \
    --enable-mpi --enable-mpi-io --enable-function-cache \
    --enable-optimize-p1 --disable-progress-bar --enable-static \
    --prefix=<your-install-path>
  ```
- ufc2-hpc built with cmake before dolfin-hpc

---

## Agreed priorities

1. ✅ Consolidate codebase — largely complete (pending Joel's memory patch)
2. ⬜ Dardel build guide — waiting for Dardel to come back up
3. ⬜ Architecture overview — how the four repos relate, data flow UFL → assembled matrix
4. ⬜ Minimal working example — Stokes or simple NS demo runnable on Dardel in week 1
5. ⬜ Simple CI / regression check

### Parallel tracks
- FEniCSx on MacBook M5 Pro (local development, teaching)
- HeartSolver port from 0.8 to 0.9 (6-month target)
- Novel framework exploration (JAX-FEM, Firedrake, Julia/Gridap)

---

## HeartSolver port — current task

**Goal:** Assess the scope of porting HeartSolver from dolfin-hpc 0.8 to 0.9.

**Why:** dolfin-hpc 0.8 is maintenance-only. HeartSolver is the only reason 0.8 still exists.
Completing the port allows us to retire 0.8 entirely.

**Key differences between 0.8 and 0.9:**
- UFC version: 0.8 uses UFC 1.1; 0.9 uses ufc2-hpc (cmake-based)
- Mesh format: 0.8 reads `.xml`; 0.9 requires `.bin` (use dolfin-convert)
- Form compilation: ffc-hpc 0.9 pipeline differs from 0.8
- API changes: significant C++ API restructuring between versions
- DirichletBC: 0.8 required a custom `update()` method added for HeartSolver
- DofMap.h: 0.8 required a manual fix to line 114 (not needed in 0.9)

**HeartSolver source location:**
- Flat snapshot: `~/repos/heartsolver_repos/fenicsheartsolver-heartsolver-f3923f758c3d/`
- Key files: `main.cpp` (21KB), `src/`, `include/`, `build/Makefile`, `data/`
- The `porous` branch was the last active Dardel branch (Joel Kronborg's work)

**Assessment task for this Claude Code session:**
1. Read `main.cpp` in full at `~/repos/heartsolver_repos/fenicsheartsolver-heartsolver-f3923f758c3d/main.cpp`
2. Read all files in `~/repos/heartsolver_repos/fenicsheartsolver-heartsolver-f3923f758c3d/src/` and `include/`
3. Read `~/repos/heartsolver_repos/fenicsheartsolver-heartsolver-f3923f758c3d/build/Makefile`
4. Identify every dolfin-hpc API call, header include, and UFC/UFL dependency
5. Cross-reference against the dolfin-hpc 0.9 headers in `~/repos/dolfin-hpc/include/`
6. Produce a structured report: what maps directly, what needs adaptation, what is missing
7. Estimate effort in person-weeks for a focused MSc student

---

## Development philosophy

- Conservative modernization — HPC environments are traditional
- Do not break working functionality
- GPU and Apple Silicon support are a later phase
- Kokkos is on the radar for future HPC portability but not a current priority
- All commits via PR to dev, never direct push to master
- Commit messages follow conventional commits format: `fix:`, `feat:`, `docs:`, `refactor:`

---

## Hardware targets

1. **Dardel** (PDC/KTH, HPE Cray EX) — primary HPC target, AMD CPUs
2. **LUMI** (CSC, AMD GPUs) — secondary HPC target
3. **MacBook Pro M5 Pro** (Apple Silicon) — local development via FEniCSx/conda
4. **Google Colab** — student access via Jupyter notebooks (legacy FEniCS + FEniCSx migration planned)

---

## Key contacts

- Joel Kronborg — former student, Dardel build expert, memory patch pending
- Ashish Bhole — contributed demo codes at https://github.com/ashishbhole/dolfin-hpc-codes

## FEniCSx local environment (MacBook M5)

- Installed via conda: `conda activate fenicsx`
- Version: FEniCSx 0.10.0, PETSc 3.25.1, MPICH 5.0
- Working notebooks in `~/repos/fenicsx-notebooks/`:
  - `poisson.ipynb` — P1, exact solution verification
  - `stokes.ipynb` — Taylor-Hood P2/P1, lid-driven cavity
  - `navier-stokes.ipynb` — backward Euler, Re=100, regularized lid

### MUMPS note for Apple Silicon
On M5, set `OMP_NUM_THREADS=1` before `import petsc4py` and set
`MUMPS CNTL(1)=0.1` (pivot threshold) to avoid non-deterministic
factorization errors. See `fenicsx-notebooks` commit `1686ac2`.
This is Mac-specific; on Dardel (x86, MPI-distributed) MUMPS behaves
differently and these workarounds are not needed.

## FEniCSx local environment (MacBook M5)

- Installed via conda: `conda activate fenicsx`
- Version: FEniCSx 0.10.0, PETSc 3.25.1, MPICH 5.0
- Working notebooks in ~/repos/fenicsx-notebooks/:
  - poisson.ipynb — P1, exact solution verification
  - stokes.ipynb — Taylor-Hood P2/P1, lid-driven cavity
  - navier-stokes.ipynb — backward Euler, Re=100, regularized lid

### MUMPS note for Apple Silicon
On M5, set OMP_NUM_THREADS=1 before import petsc4py and set
MUMPS CNTL(1)=0.1 (pivot threshold) to avoid non-deterministic
factorization errors. See fenicsx-notebooks commit 1686ac2.
This is Mac-specific; on Dardel (x86, MPI-distributed) MUMPS behaves
differently and these workarounds are not needed.
