// Copyright (C) 2007 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/fem/LinearForm.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
LinearForm::LinearForm( Mesh & mesh )
  : Form( mesh )
{
}

//-----------------------------------------------------------------------------
LinearForm::~LinearForm()
{
}

} // namespace dolfin
