// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#include <dolfin/ufl/UFLFiniteElementBase.h>
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
FiniteElementBase * FiniteElementBase::create(Object::repr_t const repr)
{
  FiniteElementBase * ret = NULL;

  std::string name = Class::make_name(repr);
  if (name == "FiniteElement")
  {
    ret = new FiniteElement(repr);
  }
  else if (name == "VectorElement")
  {
    ret = new VectorElement(repr);
  }
  /*
  else if (name == "MixedElement")
  {
    ret = new MixedElement(repr);
  }
  else if (name == "EnrichedElement")
  {
    ret = new EnrichedElement(repr);
  }
  else if (name == "RestrictedElement")
  {
    ret = new RestrictedElement(repr);
  }
  else if (name == "TensorElement")
  {
    ret = new TensorElement(repr);
  }
  */
  else
  {
    error("Unknown type of ufl::FiniteElementBase: '" + name + "'");
  }

  return ret;
}

//-----------------------------------------------------------------------------
FiniteElementBase::FiniteElementBase(std::string const& name,
                                     Family::Type const& family,
                                     Cell const& cell, dolfin::uint const degree,
                                     QuadratureScheme quad_scheme,
                                     ValueArray value_shape) :
    Class(name),
    family_(family),
    cell_(cell),
    degree_(degree),
    quad_scheme_(quad_scheme),
    value_shape_(value_shape)
{
}

//-----------------------------------------------------------------------------
FiniteElementBase::FiniteElementBase(std::string const& name, repr_t repr) :
    Class(name, repr),
    family_(arg(0)),
    cell_(arg(1)),
    degree_(arg(2)),
    quad_scheme_(),
    value_shape_()
{
}

//-----------------------------------------------------------------------------
FiniteElementBase::~FiniteElementBase()
{
}

//-----------------------------------------------------------------------------
Family const FiniteElementBase::family() const
{
  return family_;
}

//-----------------------------------------------------------------------------
Cell const FiniteElementBase::cell() const
{
  return cell_;
}

//-----------------------------------------------------------------------------
dolfin::uint const FiniteElementBase::degree() const
{
  return degree_;
}

//-----------------------------------------------------------------------------
QuadratureScheme const FiniteElementBase::quadrature_scheme() const
{
  return quad_scheme_;
}

//-----------------------------------------------------------------------------
ValueArray const FiniteElementBase::value_shape() const
{
  return value_shape_;
}

//-----------------------------------------------------------------------------
bool FiniteElementBase::component_is_valid(ValueArray const& i) const
{
  dolfin::uint r = value_shape_.size();
  bool range_ok = true;
  for (size_t idx = 0; idx < value_shape_.size(); ++idx)
  {
    range_ok = range_ok && (i[idx] < value_shape_[idx]);
  }
  return (i.size() == r && range_ok);
}

//-----------------------------------------------------------------------------
void FiniteElementBase::check_component(ValueArray const& i) const
{
  if (!component_is_valid(i))
  {
    error("Requested component is invalid");
  }
}

//-----------------------------------------------------------------------------
Cell const FiniteElementBase::get_cell(FiniteElementBaseList const& elements)
{
  FiniteElementBaseList::const_iterator it = elements.begin();
  Cell ret = (*it)->cell();
  for (++it; it != elements.end(); ++it)
  {
    if (ret.repr() != (*it)->cell().repr())
    {
      error("All subelements of mixed element should have the same cell.");
    }
  }
  return ret;
}

//-----------------------------------------------------------------------------
dolfin::uint const FiniteElementBase::get_degree_max(
    FiniteElementBaseList const& elements)
{
  dolfin::uint ret = 0;
  for (FiniteElementBaseList::const_iterator it = elements.begin();
      it != elements.end(); ++it)
  {
    ret = std::max((*it)->degree(), ret);
  }
  return ret;
}

//-----------------------------------------------------------------------------
void FiniteElementBase::display() const
{
  Class::display();
  std::cout << std::setw(24) << "family" << " = " << this->family().str()
      << std::endl;
  std::cout << std::setw(24) << "cell" << " = " << this->cell().str()
      << std::endl;
  std::cout << std::setw(24) << "degree" << " = " << this->degree()
      << std::endl;
  std::cout << std::setw(24) << "quadrature_scheme" << " = "
      << this->quadrature_scheme().str() << std::endl;
  std::cout << std::endl;
}

}

