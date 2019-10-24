// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/ufl/UFLGeometricQuantity.h>

namespace ufl
{

//-----------------------------------------------------------------------------
GeometricQuantity::GeometricQuantity(std::string const& name, Cell const& cell) :
    Class(name),
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

//-----------------------------------------------------------------------------
void GeometricQuantity::display() const
{
  Class::display();
  std::cout << std::endl;
}

}
