// Copyright (C) 2013 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2013-09-12
// Last changed: 2013-09-12

#include <dolfin/fem/FiniteElement.h>

#include <dolfin/elements/ElementLibrary.h>
#include <dolfin/elements/Elements.h>
#include <dolfin/fem/Form.h>

#include <algorithm>
#include <iomanip>

namespace dolfin
{

//-----------------------------------------------------------------------------
FiniteElement::FiniteElement(ufc::finite_element const& element,
                             bool const owner) :
    ufc_finite_element_((owner ? &element : element.create())),
    sub_value_dims_(NULL)
{
  Initialize();
}

//-----------------------------------------------------------------------------
FiniteElement::FiniteElement(ufc::finite_element const& element, uint const i) :
    ufc_finite_element_(element.create_sub_element(i)),
    sub_value_dims_(NULL)
{
  Initialize();
}

//-----------------------------------------------------------------------------
FiniteElement::FiniteElement(ufc::finite_element const& element,
                             Array<uint> const& sub_system) :
    ufc_finite_element_(FiniteElement::create_sub_element(element, sub_system)),
    sub_value_dims_(NULL)
{
  Initialize();
}

//-----------------------------------------------------------------------------
FiniteElement::FiniteElement(CellType const& type, Form& form, uint const i) :
    ufc_finite_element_(NULL),
    sub_value_dims_(NULL)
{
  // Check argument
  uint const num_arguments = form.rank() + form.num_coefficients();
  if (i >= num_arguments)
  {
    error("Illegal function index %d. Form only has %d arguments.", i,
          num_arguments);
  }

  // Create finite element
  ufc_finite_element_ = form.create_finite_element(i);

  Initialize();
}

//-----------------------------------------------------------------------------
FiniteElement::FiniteElement(ufl::FiniteElementSpace const& element) :
    ufc_finite_element_(ElementLibrary::create_finite_element(element.repr())),
    sub_value_dims_(NULL)
{
  Initialize();
}

//-----------------------------------------------------------------------------
FiniteElement::FiniteElement(FiniteElement const& other) :
    ufc_finite_element_(other.create()),
    sub_value_dims_(NULL)
{
  Initialize();
}

//-----------------------------------------------------------------------------
FiniteElement::~FiniteElement()
{
  while (!flattened_.empty())
  {
    delete flattened_.back();
    flattened_.pop_back();
  }
  delete[] sub_value_dims_;
  delete[] sub_value_offs_;
  delete ufc_finite_element_;
  ufc_finite_element_ = NULL;
}

//-----------------------------------------------------------------------------
bool FiniteElement::operator ==(FiniteElement const& other) const
{
  return (std::strcmp(this->signature(), other.signature()) == 0);
}

//-----------------------------------------------------------------------------
bool FiniteElement::operator !=(FiniteElement const& other) const
{
  return !(*this == other);
}

//-----------------------------------------------------------------------------
void FiniteElement::Initialize()
{
  dolfin_assert(ufc_finite_element_);

  // Add sub value dimensions for mixed elements, packed by axis
  uint const max_dim = ufc_finite_element_->value_rank() + 1;
  sub_value_dims_ = new Array<uint> [max_dim];
  sub_value_offs_ = new Array<uint> [max_dim];
  uint nb_subs = this->num_sub_elements();
  if (nb_subs > 0)
  {
    uint * off = new uint[max_dim];
    std::fill_n(off, max_dim, 0);
    for (uint e = 0; e < nb_subs; ++e)
    {
      ufc::finite_element * sub_fe = ufc_finite_element_->create_sub_element(e);
      for (uint a = 0; a < max_dim; ++a)
      {
        sub_value_dims_[a].push_back(sub_fe->value_dimension(a));
        sub_value_offs_[a].push_back(off[a]);
        off[a] += sub_fe->value_dimension(a);
      }
      delete sub_fe;
    }
    delete [] off;
  }
  else
  {
    for (uint a = 0; a < max_dim; ++a)
    {
      sub_value_dims_[a].push_back(value_dimension(a));
      sub_value_offs_[a].push_back(0);
    }
  }
}

//-----------------------------------------------------------------------------
uint FiniteElement::value_size() const
{
  uint size = 1;
  for (uint i = 0; i < ufc_finite_element_->value_rank(); ++i)
  {
    size *= ufc_finite_element_->value_dimension(i);
  }
  return size;
}

//-----------------------------------------------------------------------------
ufc::finite_element*
FiniteElement::create_sub_element(Array<uint> const& sub_system) const
{
  return FiniteElement::create_sub_element(*ufc_finite_element_, sub_system);
}

//-----------------------------------------------------------------------------
ufc::finite_element*
FiniteElement::create_sub_element(const ufc::finite_element& finite_element,
                                  Array<uint> const& sub_system)
{
  // If the subsystem is empty return self
  if (sub_system.size() == 0)
  {
    //error("Unable to extract sub system (no sub system specified).");
    return finite_element.create();
  }

  // Check if there are any sub systems
  if (finite_element.num_sub_elements() == 0)
  {
    error("Unable to extract sub system (there are no sub systems).");
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
  ufc::finite_element* sub_sub_element = create_sub_element(*sub_element,
                                                            sub_sub_system);
  delete sub_element;

  return sub_sub_element;
}

//-----------------------------------------------------------------------------
Array<uint> const& FiniteElement::sub_value_dimensions(uint i) const
{
  return sub_value_dims_[i];
}

//-----------------------------------------------------------------------------
Array<uint> const& FiniteElement::sub_value_offsets(uint i) const
{
  return sub_value_offs_[i];
}

//-----------------------------------------------------------------------------
Array<ufc::finite_element const *> const& FiniteElement::flatten() const
{
  if (flattened_.empty())
  {
    flatten(ufc_finite_element_, flattened_);
  }
  return flattened_;
}

//-----------------------------------------------------------------------------
void FiniteElement::flatten(ufc::finite_element const * element,
                            Array<ufc::finite_element const *>& stack,
                            uint maxlevel)
{
  // Single root element or max level is set to zero, return immediately
  if (element->num_sub_elements() == 0 || maxlevel == 0)
  {
    stack.push_back(element->create());
    return;
  }
  // Go one level down
  for (uint s = 0; s < element->num_sub_elements(); ++s)
  {
    ufc::finite_element const * sub = element->create_sub_element(s);
    if (sub->num_sub_elements() == 0)
    {
      // Leaf element
      stack.push_back(sub);
    }
    else
    {
      // Branch
      FiniteElement::flatten(sub, stack, maxlevel - 1);
    }
  }
}

//-----------------------------------------------------------------------------
void FiniteElement::flatten(ufc::finite_element const * element,
                            Array<ufc::finite_element const *>& stack)
{
  // Single root element or max level is set to zero, return immediately
  if (element->num_sub_elements() == 0)
  {
    stack.push_back(element->create());
    return;
  }
  // Go one level down
  for (uint s = 0; s < element->num_sub_elements(); ++s)
  {
    ufc::finite_element const * sub = element->create_sub_element(s);
    if (sub->num_sub_elements() == 0)
    {
      // Leaf element
      stack.push_back(sub);
    }
    else
    {
      // Branch
      FiniteElement::flatten(sub, stack);
    }
  }
}

//-----------------------------------------------------------------------------
bool FiniteElement::is_vectorizable() const
{
  bool ret = true;
  Array<ufc::finite_element const*> const& flt = this->flatten();
  for (uint s = 1; s < flt.size(); ++s)
  {
    if (std::strcmp(flt[0]->signature(), flt[s]->signature()) != 0)
    {
      ret = false;
      break;
    }
  }
  return ret;
}

//-----------------------------------------------------------------------------
void FiniteElement::disp() const
{
  section("FiniteElement");
  begin("ufc::finite_element info");
  prm("Signature"             , ufc_finite_element_->signature());
  prm("Cell shape"            , ufc_finite_element_->cell_shape());
  prm("Topological dimension" , ufc_finite_element_->topological_dimension());
  prm("Geometric dimension"   , ufc_finite_element_->geometric_dimension());
  prm("Space dimension"       , ufc_finite_element_->space_dimension());
  prm("Value rank"            , ufc_finite_element_->value_rank());
  prm("Value dimension"       , ufc_finite_element_->value_dimension(0));
  prm("Nb of sub elements"    , ufc_finite_element_->num_sub_elements());
  end();
  end();
}

}

