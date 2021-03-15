// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/fem/Functional.h>

namespace dolfin
{

//-----------------------------------------------------------------------------

Functional::Functional( Mesh & mesh )
  : Form( mesh )
{
}

//-----------------------------------------------------------------------------

} // namespace dolfin
