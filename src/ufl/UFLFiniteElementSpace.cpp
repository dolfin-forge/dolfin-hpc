// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/ufl/UFLFiniteElementSpace.h>
#include <dolfin/ufl/UFLFiniteElement.h>
#include <dolfin/ufl/UFLVectorElement.h>
#include <dolfin/ufl/UFLMixedElement.h>
#include <dolfin/ufl/UFLEnrichedElement.h>
#include <dolfin/ufl/UFLRestrictedElement.h>
#include <dolfin/ufl/UFLTensorElement.h>

namespace ufl
{

using dolfin::error;

//-----------------------------------------------------------------------------
FiniteElementSpace * FiniteElementSpace::create(Object::repr_t const repr)
{
  std::string name = Class::make_name(repr);
  if (name == "FiniteElement")
  {
    return new FiniteElement(repr);
  }
  else if (name == "VectorElement")
  {
    return new VectorElement(repr);
  }
  else if (name == "TensorElement")
  {
    return new TensorElement(repr);
  }
  else if (name == "MixedElement")
  {
    return new MixedElement(repr);
  }
  else if (name == "EnrichedElement")
  {
    error("UFL implementation for EnrichedElement is not available yet.");
//FIXME:    return new EnrichedElement(repr);
  }
  else if (name == "RestrictedElement")
  {
    error("UFL implementation for RestrictedElement is not available yet.");
//FIXME:    return new RestrictedElement(repr);
  }
  else
  {
    error("Unknown type of ufl::FiniteElementSpace: '" + name + "'");
  }
  
  return NULL;
}

//-----------------------------------------------------------------------------
FiniteElementSpace::FiniteElementSpace(std::string const& name,
                                     QuadratureScheme quad_scheme) :
    Class(name),
    quad_scheme_(quad_scheme)
{
}

//-----------------------------------------------------------------------------
FiniteElementSpace::FiniteElementSpace(std::string const& name, repr_t repr) :
    Class(name, repr),
    quad_scheme_()
{
}

//-----------------------------------------------------------------------------
FiniteElementSpace::~FiniteElementSpace()
{
}

//-----------------------------------------------------------------------------
QuadratureScheme const& FiniteElementSpace::quadrature_scheme() const
{
  return quad_scheme_;
}

//-----------------------------------------------------------------------------
dolfin::uint FiniteElementSpace::value_size() const
{
  return std::max(this->value_shape().prod(), (dolfin::uint) 1);
}

//-----------------------------------------------------------------------------
bool FiniteElementSpace::component_is_valid(ValueArray const& i) const
{
  dolfin::uint r = this->value_shape().size();
  bool range_ok = true;
  for (size_t idx = 0; idx < this->value_shape().size(); ++idx)
  {
    range_ok = range_ok && (i[idx] < this->value_shape()[idx]);
  }
  return (i.size() == r && range_ok);
}

//-----------------------------------------------------------------------------
void FiniteElementSpace::check_component(ValueArray const& i) const
{
  if (!component_is_valid(i))
  {
    error("Requested component is invalid");
  }
}

//-----------------------------------------------------------------------------
Cell FiniteElementSpace::get_cell(List const& elements) const
{
  if (elements.size() < 1)
  {
    error("Invalid element list given as argument to get_cell");
  }
  List::const_iterator it = elements.begin();
  Cell ret = (*it)->cell();
  for (++it; it != elements.end(); ++it)
  {
    if (ret != (*it)->cell())
    {
      error("All subelements of mixed element should have the same cell.");
    }
  }
  return ret;
}

//-----------------------------------------------------------------------------
dolfin::uint FiniteElementSpace::get_degree_max(List const& elements) const
{
  dolfin::uint ret = 0;
  for (List::const_iterator it = elements.begin(); it != elements.end(); ++it)
  {
    ret = std::max((dolfin::uint) (*it)->degree(), ret);
  }
  return ret;
}

//-----------------------------------------------------------------------------
dolfin::uint FiniteElementSpace::get_value_size(List const& elements) const
{
  // Compute value size
  dolfin::uint ret = 0;
  for (List::const_iterator it = elements.begin(); it != elements.end(); ++it)
  {
    ret += (*it)->value_size();
  }
  return ret;
}

//-----------------------------------------------------------------------------
void FiniteElementSpace::display() const
{
  Class::display();
  std::cout << std::setw(24) << "family" << " = " << this->family().str()
            << std::endl;
  std::cout << std::setw(24) << "cell" << " = " << this->cell().str()
            << std::endl;
  std::cout << std::setw(24) << "degree" << " = " << this->degree()
            << std::endl;
  std::cout << std::setw(24) << "value_shape" << " = "
            << this->value_shape().str() << std::endl;
  std::cout << std::setw(24) << "quadrature_scheme" << " = "
            << this->quadrature_scheme().str() << std::endl;
  std::cout << std::endl;
}

}

