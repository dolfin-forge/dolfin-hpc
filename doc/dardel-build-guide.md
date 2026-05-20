# Dardel Build Guide — dolfin-hpc 0.9.5

**Target audience:** New MSc students with basic Linux skills.  
**Goal:** Build and run the Poisson demo on Dardel within one week of onboarding.  
**System:** Dardel (PDC/KTH, HPE Cray EX), software environment 24.11.

---

## 1. Prerequisites

You need:
- A PDC account (apply via SUPR: https://supr.naiss.se)
- Membership in the project allocation (ask Johan)
- SSH access to Dardel: `ssh <username>@dardel.pdc.kth.se`

---

## 2. Environment setup

Dardel uses the Lmod module system. The dolfin-hpc installation is provided
as a user module. Load it with exactly two commands:

```bash
ml PDC
ml dolfin-hpc/0.9.5
```

For linking against ParMETIS (required for mesh partitioning), also load:

```bash
ml parmetis/4.0.3-cpeCray-24.11
```

Verify the environment:

```bash
pkg-config --modversion dolfin
```

Expected output: `0.9.5-hpc`

**Note:** The pkg-config file is named `dolfin.pc`, not `dolfin-hpc.pc`.
Always use `pkg-config --cflags dolfin` and `pkg-config --libs dolfin`.

---

## 3. Build the Poisson demo

The installed module includes pre-compiled demo sources with pre-generated
form headers (`Poisson.h`). Copy the Poisson demo to your home directory:

```bash
cp -r /pdc/software/24.11/other/dolfin_hpc/dolfin-hpc/demo/pde/poisson ~/poisson_test
cd ~/poisson_test
ls
```

You should see:
```
Makefile  Poisson.h  Poisson.ufl  README  UnitSquareMesh_32x32.bin  main.cpp
```

Build with make:

```bash
make
```

Expected output: two deprecation warnings about `PETSC_NULL` (harmless,
these are fixed in the dolfin-hpc source but not in the installed headers),
followed by successful linking. The binary `demo` is produced.

Verify:
```bash
ls -la demo
```

---

## 4. Run the demo on a compute node

Dardel requires all MPI jobs to run on compute nodes, not the login node.
Submit an interactive allocation (replace `<project-id>` with your NAISS
project, e.g. `naiss2026-3-301`):

```bash
salloc -n 4 -A <project-id> -t 00:30:00 -p shared
```

> **Note:** After `salloc` succeeds, the shell moves to the compute node and all loaded
> modules are reset. You must reload them before running anything.

```bash
ml PDC
ml dolfin-hpc/0.9.5
ml parmetis/4.0.3-cpeCray-24.11
```

Run the demo:

```bash
cd ~/poisson_test
srun -n 4 ./demo
```

Expected output: DOLFIN solver log showing mesh partitioning across 4 MPI
ranks, assembly, and linear solve.

---

## 5. Module summary

| Module | Purpose |
|---|---|
| `PDC` | Enables PDC software stack |
| `dolfin-hpc/0.9.5` | Core library, headers, pkg-config |
| `parmetis/4.0.3-cpeCray-24.11` | Mesh partitioning (linking) |

The dolfin-hpc module sets:
- `PKG_CONFIG_PATH` to include `dolfin.pc`
- Compiler `CC` (Cray C++ wrapper, wraps clang++)
- Flags: `-std=c++14 -g -O2`
- PETSc 3.22.2 include and library paths

---

## 6. Mesh generation

New meshes cannot currently be generated from Gmsh directly on Dardel.
`dolfin-convert` (the tool that converts Gmsh `.msh` files to the `.bin`
format required by dolfin-hpc 0.9.x) requires Python 2, which is not
available on Dardel (Python 3.11 only).

For demos and initial work, use the pre-converted `.bin` meshes provided
with the installed module. The Poisson demo includes `UnitSquareMesh_32x32.bin`
and can be run immediately without any mesh conversion step.

If you need a custom mesh, generate and convert it on a machine with Python 2
available, then transfer the `.bin` file to Dardel via `scp`.

See `KNOWN_ISSUES.md` §5 for details and the long-term fix (Python 3 port or
replacement with `meshio`).

---

## 7. Writing your own solver

For a new solver `mysolver.cpp` that uses an existing pre-compiled form
header `MyForm.h`:

```bash
CC $(pkg-config --cflags dolfin) -c mysolver.cpp
CC -o mysolver mysolver.o $(pkg-config --cflags dolfin) $(pkg-config --libs dolfin)
```

Or add a `Makefile` following the pattern in the Poisson demo.

---

## 8. Compiling new UFL forms

To compile a new `.ufl` form file into a C++ header, `ffc` (FEniCS Form
Compiler) is needed. The ffc-hpc installation is at:

```
/pdc/software/24.11/other/dolfin_hpc/ffc-hpc/bin/ffc
```

However, ffc requires FIAT (Finite element Automatic Tabulator) which is
**not currently installed** on Dardel. This means new UFL forms cannot be
compiled without first installing FIAT.

**Workaround:** Use the pre-compiled form headers provided with the demos,
or compile forms locally on your local machine and transfer the generated `.h`
file to Dardel via `scp`.

To compile a form locally, `ffc` is needed. Two options:

- **FEniCSx conda environment** (if you have it set up):
  ```bash
  conda activate fenicsx
  ffc -l dolfin MyForm.ufl
  ```
- **Legacy FEniCS** (if installed separately):
  ```bash
  ffc -l dolfin MyForm.ufl
  ```

Both produce `MyForm.h`. Transfer it to Dardel with:
```bash
scp MyForm.h <username>@dardel.pdc.kth.se:~/your-project/
```

The file can then be used directly in compilation on Dardel.

**TODO:** Install FIAT on Dardel to enable on-system form compilation.
Track in KNOWN_ISSUES.md.

---

## 9. Batch job script

For production runs, use a batch script. Save as `run.sh`:

```bash
#!/bin/bash
#SBATCH -A <project-id>
#SBATCH -J poisson_demo
#SBATCH -p main
#SBATCH -n 128
#SBATCH -t 01:00:00

ml PDC
ml dolfin-hpc/0.9.5
ml parmetis/4.0.3-cpeCray-24.11

srun ./demo
```

> **Note:** Replace `<project-id>` with your NAISS allocation ID (e.g. `naiss2025-5-152`).
> Ask your PI for the correct ID, or find it yourself by running `projinfo` or `groups`
> on the Dardel login node.

Submit with:
```bash
sbatch run.sh
```

Monitor with:
```bash
squeue -u $USER
```

---

## 10. Troubleshooting

**`pkg-config: dolfin not found`**  
You forgot to load the modules. Run `ml PDC && ml dolfin-hpc/0.9.5`.

**`unable to find library -lparmetis`**  
Load `ml parmetis/4.0.3-cpeCray-24.11` before linking.

**`ffc: No module named FIAT`**  
FIAT is not installed. Compile forms locally and transfer the `.h` file.
See section 8.

**`AssocGrpSubmitJobsLimit`**  
Your allocation is exhausted or not yet activated. Contact your PI or
check SUPR for allocation status.

**PETSC_NULL deprecation warnings**  
Harmless. These are fixed in the dolfin-hpc 0.9.5 source but not in the
installed headers. They do not affect correctness.

---

## 11. Installed paths reference

| Component | Path |
|---|---|
| Library | `/pdc/software/24.11/other/dolfin_hpc/install/lib/` |
| Headers | `/pdc/software/24.11/other/dolfin_hpc/install/include/` |
| pkg-config | `/pdc/software/24.11/other/dolfin_hpc/install/lib/pkgconfig/dolfin.pc` |
| Demos | `/pdc/software/24.11/other/dolfin_hpc/dolfin-hpc/demo/` |
| ffc | `/pdc/software/24.11/other/dolfin_hpc/ffc-hpc/bin/ffc` |
| ufl-hpc egg | `/pdc/software/24.11/other/dolfin_hpc/ufl-hpc/lib/python3.11/site-packages/` |
| PETSc | `/cfs/klemming/pdc/software/dardel/24.11/other/petsc/3.22.2/` |
| ParMETIS | `/pdc/software/24.11/eb/software/parmetis/4.0.3-cpeCray-24.11/` |

