// Copyright (C) 2020 Julian Hornich
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_TRILINOS_MATRIX_H
#define __DOLFIN_TRILINOS_MATRIX_H

#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_TRILINOS

#include <dolfin/common/Variable.h>
#include <dolfin/la/GenericMatrix.h>
#include <dolfin/la/trilinos/TrilinosObject.h>

namespace dolfin
{

namespace trilinos
{

class Matrix : public GenericMatrix, public Object, public Variable
{
};

} // end namespace trilinos

} // end namespace dolfin

#endif // HAVE_TRILINOS

#endif // __DOLFIN_TRILINOS_MATRIX_H
