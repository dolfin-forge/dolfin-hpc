// Copyright (C) 2005 Johan Jansson.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_PETSC

#include <dolfin/la/petsc/PETScPreconditioner.h>

#include <dolfin/la/petsc/PETScVector.h>

#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 6
#include <petsc/private/pcimpl.h>
#elif PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 3
#include <petsc-private/pcimpl.h>
#else
#include <private/pcimpl.h>
#endif


using namespace dolfin;

//-----------------------------------------------------------------------------
void PETScPreconditioner::setup(const KSP ksp, PETScPreconditioner &pc)
{
  PC petscpc;
  KSPGetPC(ksp, &petscpc);

  PETScPreconditioner::PCCreate(petscpc);

  petscpc->data = &pc;
  petscpc->ops->apply = PETScPreconditioner::PCApply;
  petscpc->ops->applytranspose = PETScPreconditioner::PCApply;
  petscpc->ops->applysymmetricleft = PETScPreconditioner::PCApply;
  petscpc->ops->applysymmetricright = PETScPreconditioner::PCApply;
}
//-----------------------------------------------------------------------------
auto PETScPreconditioner::PCApply(PC pc, Vec x, Vec y) -> int
{
  // Convert vectors to DOLFIN wrapper format and pass to DOLFIN preconditioner

  PETScPreconditioner* newpc = (PETScPreconditioner*)pc->data;

  PETScVector dolfinx(x), dolfiny(y);

  newpc->solve(dolfiny, dolfinx);

  return 0;
}
//-----------------------------------------------------------------------------
auto PETScPreconditioner::PCCreate(PC pc) -> int
{
  // Initialize function pointers to 0

  pc->ops->setup               = nullptr;
  pc->ops->apply               = nullptr;
  pc->ops->applyrichardson     = nullptr;
  pc->ops->applytranspose      = nullptr;
  pc->ops->applysymmetricleft  = nullptr;
  pc->ops->applysymmetricright = nullptr;
  pc->ops->setfromoptions      = nullptr;
  pc->ops->view                = nullptr;
  pc->ops->destroy             = nullptr;

  // Set PETSc name of preconditioner
  PetscObjectChangeTypeName((PetscObject)pc, "DOLFIN");

  return 0;
}
//-----------------------------------------------------------------------------
auto PETScPreconditioner::getType(PreconditionerType pc) -> PCType
{
  switch (pc)
  {
  case default_pc:
    return "default";
  case amg:
    return PCHYPRE;
  case icc:
    return PCICC;
  case ilu:
    return PCILU;
  case jacobi:
    return PCJACOBI;
  case bjacobi:
    return PCBJACOBI;
  case sor:
    return PCSOR;
  case none:
    return PCNONE;
  default:
    warning("Requested preconditioner unkown. Using incomplete LU.");
    return PCILU;
  }
}
//-----------------------------------------------------------------------------

#endif
