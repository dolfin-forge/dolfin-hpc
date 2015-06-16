// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#include <dolfin/ufl/UFLSpace.h>

#include <dolfin/log/log.h>

#include <iomanip>
#include <sstream>

namespace ufl
{

//-----------------------------------------------------------------------------
Space::Space(dolfin::uint const& dim) :
    Class("Space"),
    dimension_(dim),
    repr_(*this, dimension_)
{
  std::stringstream ssstr;
  ssstr << "R" << dimension_;
  str_ = ssstr.str();

  if(dimension_ > 3)
  {
    dolfin::error("Number of space dimension is greater than 3");
  }
}

//-----------------------------------------------------------------------------
Space::Space(repr_t const& repr) :
    Class("Space", repr),
    dimension_(arg(0)),
    repr_(repr)
{
  std::stringstream ssstr;
  ssstr << "R" << dimension_;
  str_ = ssstr.str();

  if(dimension_ > 3)
  {
    dolfin::error("Number of space dimension is greater than 3");
  }
}

//-----------------------------------------------------------------------------
Space::~Space()
{
}

//-----------------------------------------------------------------------------
dolfin::uint Space::dimension() const
{
  return dimension_;
}

//-----------------------------------------------------------------------------
Object::repr_t const& Space::repr() const
{
  return repr_;
}

//-----------------------------------------------------------------------------
std::string const& Space::str() const
{
  return str_;
}

//-----------------------------------------------------------------------------
void Space::display() const
{
  Class::display();
  std::cout << std::setw(24) << "dimension" << " = " << this->dimension()
      << std::endl;
  std::cout << std::endl;
}

}
