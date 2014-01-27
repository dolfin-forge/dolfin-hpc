// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#include <dolfin/ufl/UFLSpace.h>

#include <iomanip>
#include <sstream>

namespace ufl
{

//-----------------------------------------------------------------------------
Space::Space(uint const& dim) :
    Class("Space"),
    dimension_(dim)
{
  std::stringstream ssrepr;
  ssrepr << "Space(" << dimension_ << ")";
  repr_ = ssrepr.str();

  std::stringstream ssstr;
  ssstr << "R" << dimension_;
  str_ = ssstr.str();
}

//-----------------------------------------------------------------------------
Space::~Space()
{
}

//-----------------------------------------------------------------------------
uint Space::dimension() const
{
  return dimension_;
}

//-----------------------------------------------------------------------------
Object::repr_t const Space::repr() const
{
  return repr_;
}

//-----------------------------------------------------------------------------
std::string const Space::str() const
{
  return str_;
}

//-----------------------------------------------------------------------------
void Space::display() const
{
  Class::display();
  std::cout << std::setw(16) << "dimension" << " = " << this->dimension() << std::endl;
  std::cout << std::endl;
}

}
