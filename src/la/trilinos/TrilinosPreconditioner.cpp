// Copyright (C) 2021 Julian Hornich
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_TRILINOS

#include <dolfin/la/trilinos/TrilinosPreconditioner.h>

#include <dolfin/la/PreconditionerType.h>
#include <dolfin/la/trilinos/TrilinosVector.h>


namespace dolfin
{

namespace trilinos
{

//-----------------------------------------------------------------------------

// void Preconditioner::setup(const KSP ksp, Preconditioner &pc)
// {
// }

//-----------------------------------------------------------------------------

// auto Preconditioner::PCApply(PC pc, Vec x, Vec y) -> int
// {
// }

//-----------------------------------------------------------------------------

// auto Preconditioner::PCCreate(PC pc) -> int
// {
// }

//-----------------------------------------------------------------------------

// auto Preconditioner::getType(PreconditionerType pc) -> PCType
// {
// }

//-----------------------------------------------------------------------------

} // namespace trilinos

} // namespace dolfin

#endif // HAVE_TRILINOS
