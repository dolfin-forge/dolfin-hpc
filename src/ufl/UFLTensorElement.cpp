// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#include <dolfin/ufl/UFLTensorElement.h>

namespace ufl
{

//-----------------------------------------------------------------------------
TensorElement::TensorElement(Family::Type family, Cell const& cell,
                             dolfin::uint const degree, dolfin::uint const dim) :
    FiniteElementBase("TensorElement"),
    family_(Family::Tensor),
    sub_element_(family, cell, degree),
    value_shape_(
        ValueArray(2, dim) + sub_element_.value_shape()),
    symmetry_(),
    sub_elements_(dim*dim, &sub_element_)
{
  createReprStr();
}

//-----------------------------------------------------------------------------
TensorElement::TensorElement(repr_t const& repr) :
    FiniteElementBase("TensorElement", repr),
    family_(Family::Tensor),
    sub_element_(Family(arg(0)).type(), Cell(arg(1)),
                 type<dolfin::uint>(arg(2))),
    value_shape_(
        ValueArray(2, type<dolfin::uint>(arg(3)))
            + sub_element_.value_shape()),
    symmetry_(),
    sub_elements_(type<dolfin::uint>(arg(3))*type<dolfin::uint>(arg(3)), &sub_element_)
{
  createReprStr();
}

//-----------------------------------------------------------------------------
TensorElement::~TensorElement()
{
}

//-----------------------------------------------------------------------------
Family const& TensorElement::family() const
{
  return sub_element_.family();
}

//-----------------------------------------------------------------------------
Family::Type TensorElement::metatype() const
{
  return Family::Tensor;
}

//-----------------------------------------------------------------------------
Cell const& TensorElement::cell() const
{
  return sub_element_.cell();
}

//-----------------------------------------------------------------------------
type<dolfin::uint> const& TensorElement::degree() const
{
  return sub_element_.degree();
}

//-----------------------------------------------------------------------------
ValueArray const& TensorElement::value_shape() const
{
  return value_shape_;
}

//-----------------------------------------------------------------------------
bool const TensorElement::is_cellwise_constant() const
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
std::map<dolfin::uint, dolfin::uint> const TensorElement::symmetry() const
{
  return symmetry_;
}

//-----------------------------------------------------------------------------
std::pair<ValueArray, ValueArray> const TensorElement::extract_subelement_component(
    ValueArray const& i) const
{
  return std::pair<ValueArray, ValueArray>();
}

//-----------------------------------------------------------------------------
std::pair<dolfin::uint, FiniteElementBase const *> const TensorElement::extract_component(
    ValueArray const& i) const
{
  return sub_element_.extract_component(i);
}

//-----------------------------------------------------------------------------
dolfin::uint const TensorElement::num_sub_elements() const
{
  return sub_elements_.size();
}

//-----------------------------------------------------------------------------
FiniteElementBase::List const& TensorElement::sub_elements() const
{
  return sub_elements_;
}

//-----------------------------------------------------------------------------
Object::repr_t const& TensorElement::repr() const
{
  return repr_;
}

//-----------------------------------------------------------------------------
std::string const& TensorElement::str() const
{
  return str_;
}

//-----------------------------------------------------------------------------
void TensorElement::createReprStr()
{
  // Create string representation
  std::stringstream ssrepr;
  ssrepr << "TensorElement(" << this->family().repr() << ", " << cell().repr()
      << ", " << this->degree() << ", " << quadrature_scheme().repr() << ")";
  repr_ = ssrepr.str();

  std::stringstream ssstr;
  ssstr << "<" << this->family().short_name() << " vector element of degree "
      << this->degree() << " on a " << cell().str() << ": "
      << sub_elements_.size() << " x " << sub_element_.str() << ">";
  str_ = ssstr.str();
}

}

