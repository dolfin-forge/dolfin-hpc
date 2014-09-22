// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-02-27
// Last changed: 2014-02-27

#include <dolfin/fem/Functional.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
Functional::Functional(Mesh& mesh) :
    Form(mesh)
{
}

//-----------------------------------------------------------------------------
Functional::~Functional()
{
}

}
