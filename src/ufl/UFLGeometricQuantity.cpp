// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#include <dolfin/ufl/UFLGeometricQuantity.h>

namespace ufl
{

//-----------------------------------------------------------------------------
GeometricQuantity::GeometricQuantity(Cell const& cell):
    cell_(cell)
{
}

//-----------------------------------------------------------------------------
GeometricQuantity::~GeometricQuantity()
{
}

//-----------------------------------------------------------------------------
Cell const& GeometricQuantity::cell()
{
  return cell_;
}

}
