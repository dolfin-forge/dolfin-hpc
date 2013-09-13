// Copyright (C) 2013 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2013-09-12
// Last changed: 2013-09-12

#include <dolfin/fem/FiniteElement.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
FiniteElement::FiniteElement(std::string const& signature) :
    ufc_finite_element_(ElementLibrary::create_finite_element(signature)),
    finite_element_local_(true),
    degree_(0)

{

}

//-----------------------------------------------------------------------------
FiniteElement::FiniteElement(Mesh& mesh, Form& form, uint i) :
    ufc_finite_element_(NULL),
    finite_element_local_(true),
    degree_(0)
{
  // Check argument
  uint const num_arguments = form.form().rank()
      + form.form().num_coefficients();
  if (i >= num_arguments)
  {
    error("Illegal function index %d. Form only has %d arguments.", i,
        num_arguments);
  }

  // Create finite element
  ufc_finite_element_ = form.form().create_finite_element(i);

  // Update dof maps
  form.updateDofMaps(mesh);
}

//-----------------------------------------------------------------------------
FiniteElement::FiniteElement(ufc::finite_element& finite_element,
                             bool const finite_element_local) :
    ufc_finite_element_(&finite_element),
    finite_element_local_(finite_element_local),
    degree_(0)
{

}

FiniteElement::~FiniteElement()

{
  if(finite_element_local_)
    delete ufc_finite_element_;
}

//-----------------------------------------------------------------------------
uint const FiniteElement::degree() const
{
  return degree_;
}

//-----------------------------------------------------------------------------
void FiniteElement::init()
{
  dolfin_assert(ufc_finite_element_);
  // Lookup signature to get degree until part of UFC 2.x interface
  // for FiniteElement and VectorElement the token offset is 3
  // if  the tokem offset is 3
  char * sign = new char[128];
  std::strcpy(sign, ufc_finite_element_->signature());
  char * tok = std::strtok(sign, ",");
  size_t offset = 3;
  size_t t = 0;
  while (tok != NULL)
  {
    ++t;
    tok = std::strtok(NULL, ",");
    if (offset == t)
    {
      std::stringstream d;
      d << tok;
      d >> degree_;
    }
  }
}

//-----------------------------------------------------------------------------
ufc::finite_element* FiniteElement::create_sub_element(Array<uint> const& sub_system) const
{
  // Recursively extract sub element
  ufc::finite_element* sub_finite_element = create_sub_element(*ufc_finite_element_,
      sub_system);
  message(2, "Extracted finite element for sub system: %s",
      sub_finite_element->signature());

  return sub_finite_element;
}
//-----------------------------------------------------------------------------
ufc::finite_element* FiniteElement::create_sub_element(
    const ufc::finite_element& finite_element, Array<uint> const& sub_system)
{
  // Check if there are any sub systems
  if (finite_element.num_sub_elements() == 0)
  {
    error("Unable to extract sub system (there are no sub systems).");
  }

  // Check that a sub system has been specified
  if (sub_system.size() == 0)
  {
    error("Unable to extract sub system (no sub system specified).");
  }

  // Check the number of available sub systems
  if (sub_system[0] >= finite_element.num_sub_elements())
  {
    error("Unable to extract sub system %d (only %d sub systems defined).",
        sub_system[0], finite_element.num_sub_elements());
  }

  // Create sub system
  ufc::finite_element* sub_element = finite_element.create_sub_element(
      sub_system[0]);

  // Return sub system if sub sub system should not be extracted
  if (sub_system.size() == 1)
  {
    return sub_element;
  }

  // Otherwise, recursively extract the sub sub system
  Array<uint> sub_sub_system;
  for (uint i = 1; i < sub_system.size(); i++)
  {
    sub_sub_system.push_back(sub_system[i]);
  }
  ufc::finite_element* sub_sub_element = create_sub_element(*sub_element, sub_sub_system);
  delete sub_element;

  return sub_sub_element;
}
//-----------------------------------------------------------------------------

}

