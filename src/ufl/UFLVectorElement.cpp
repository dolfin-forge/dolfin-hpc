// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/ufl/UFLVectorElement.h>

namespace ufl
{

using dolfin::error;

//-----------------------------------------------------------------------------
VectorElement::VectorElement(Family::Type family, Cell const& cell,
                             dolfin::uint const degree, dolfin::uint const dim) :
    FiniteElementSpace("VectorElement"),
    family_(Family::Vector),
    sub_element_(family, cell, degree),
    dim_(dim),
    value_shape_(ValueArray(1, dim_) + sub_element_.value_shape()),
    symmetry_(),
    sub_elements_(dim, &sub_element_),
    repr_(*this, sub_element_.family(), cell, sub_element_.degree(), dim_,
          sub_element_.quadrature_scheme())
{
  createStr();
}

//-----------------------------------------------------------------------------
VectorElement::VectorElement(repr_t const& repr) :
    FiniteElementSpace("VectorElement", repr),
    family_(Family::Vector),
    sub_element_(Family(arg(0)).type(), Cell(arg(1)),
                 type<dolfin::uint>(arg(2))),
    dim_(arg(3)),
    value_shape_(ValueArray(1, dim_) + sub_element_.value_shape()),
    symmetry_(),
    sub_elements_(dim_, &sub_element_),
    repr_(repr)
{
  createStr();
}

//-----------------------------------------------------------------------------
VectorElement::~VectorElement()
{
}

//-----------------------------------------------------------------------------
Family const& VectorElement::family() const
{
  return sub_element_.family();
}

//-----------------------------------------------------------------------------
Family::Type VectorElement::metatype() const
{
  return Family::Vector;
}

//-----------------------------------------------------------------------------
Cell const& VectorElement::cell() const
{
  return sub_element_.cell();
}

//-----------------------------------------------------------------------------
type<dolfin::uint> const& VectorElement::degree() const
{
  return sub_element_.degree();
}

//-----------------------------------------------------------------------------
ValueArray const& VectorElement::value_shape() const
{
  return value_shape_;
}

//-----------------------------------------------------------------------------
bool VectorElement::is_cellwise_constant() const
{
  bool ret = true;
  for (List::const_iterator it = sub_elements_.begin();
      it != sub_elements_.end(); ++it)
  {
    ret |= (*it)->is_cellwise_constant();
  }
  return ret;
}

//-----------------------------------------------------------------------------
std::map<dolfin::uint, dolfin::uint> const& VectorElement::symmetry() const
{
  return symmetry_;
}

//-----------------------------------------------------------------------------
std::pair<ValueArray, ValueArray> VectorElement::extract_subelement_component(
    ValueArray const&) const
{
  return std::pair<ValueArray, ValueArray>();
}

//-----------------------------------------------------------------------------
std::pair<dolfin::uint, FiniteElementSpace const *> VectorElement::extract_component(
    ValueArray const& i) const
{
  return sub_element_.extract_component(i);
}

//-----------------------------------------------------------------------------
dolfin::uint VectorElement::num_sub_elements() const
{
  return sub_elements_.size();
}

//-----------------------------------------------------------------------------
FiniteElementSpace::List const& VectorElement::sub_elements() const
{
  return sub_elements_;
}

//-----------------------------------------------------------------------------
Object::repr_t const& VectorElement::repr() const
{
  return repr_;
}

//-----------------------------------------------------------------------------
std::string const& VectorElement::str() const
{
  return str_;
}

//-----------------------------------------------------------------------------
void VectorElement::createStr()
{
  //TODO: Check mixed finite element definition
  if (sub_elements_.size() < 2)
  {
    error("A vector element should contain more than one subelement");
  }

  std::stringstream ssstr;
  ssstr << "<" << sub_element_.family().short_name()
        << " vector element of degree " << this->degree() << " on a "
        << cell().str() << ": " << sub_elements_.size() << " x "
        << sub_element_.str() << ">";
  str_ = ssstr.str();
}

}

