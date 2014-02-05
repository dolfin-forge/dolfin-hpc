// Copyright (C) 2013 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2013-09-12
// Last changed: 2013-09-12

#include <dolfin/fem/FiniteElement.h>

#include <algorithm>
#include <iomanip>

namespace dolfin
{

//-----------------------------------------------------------------------------
FiniteElement::FiniteElement(std::string const& signature) :
    ufc_finite_element_(ElementLibrary::create_finite_element(signature)),
    finite_element_local_(true),
    sub_value_dims_(NULL),
    topo_dim_(0),
    geom_dim_(0)

{
  Initialize();
}

//-----------------------------------------------------------------------------
FiniteElement::FiniteElement(CellType& type, Form& form, uint const i) :
    ufc_finite_element_(NULL),
    finite_element_local_(true),
    sub_value_dims_(NULL),
    topo_dim_(0),
    geom_dim_(0)
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

  Initialize();
}

//-----------------------------------------------------------------------------
FiniteElement::FiniteElement(ufc::finite_element& finite_element,
                             bool const finite_element_local) :
    ufc_finite_element_(&finite_element),
    finite_element_local_(finite_element_local),
    sub_value_dims_(NULL),
    topo_dim_(0),
    geom_dim_(0)
{
  Initialize();
}

//-----------------------------------------------------------------------------
FiniteElement::~FiniteElement()
{
  if (finite_element_local_)
    delete ufc_finite_element_;

  delete[] sub_value_dims_;
  delete[] sub_value_offs_;
}

//-----------------------------------------------------------------------------
void FiniteElement::Initialize()
{
  dolfin_assert(ufc_finite_element_);

  //Topological_dimension and geometric_dimension are generated with UFC 2.x
  //but should be inferred from the cell shape if UFL is not enabled.

#if UFC_VERSION_MAJOR >= 2
  topo_dim_ = ufc_finite_element_->topological_dimension();
  geom_dim_ = ufc_finite_element_->geometric_dimension();
#else
  switch (ufc_finite_element_->cell_shape())
  {
    case ufc::interval:
    topo_dim_ = 1;
    geom_dim_ = 1;
    break;
    case ufc::triangle:
    topo_dim_ = 2;
    geom_dim_ = 2;
    break;
    case ufc::tetrahedron:
    topo_dim_ = 3;
    geom_dim_ = 3;
    break;
    default:
    error("Unknown cell type.");
    break;
  }
#endif

  // Add sub value dimensions for mixed elements, packed by axis
  sub_value_dims_ = new Array<uint> [geom_dim_];
  sub_value_offs_ = new Array<uint> [geom_dim_];
  uint nb_subs = this->num_sub_elements();
  if (nb_subs > 0)
  {
    uint * off = new uint[geom_dim_];
    std::fill_n(off, geom_dim_, 0);
    for (uint e = 0; e < nb_subs; ++e)
    {
      ufc::finite_element * sub_fe = ufc_finite_element_->create_sub_element(e);
      for (uint a = 0; a < geom_dim_; ++a)
      {
        sub_value_dims_[a].push_back(sub_fe->value_dimension(a));
        sub_value_offs_[a].push_back(off[a]);
        off[a] += sub_fe->value_dimension(a);
      }
      delete sub_fe;
    }
  }
  else
  {
    for (uint a = 0; a < geom_dim_; ++a)
    {
      sub_value_dims_[a].push_back(value_dimension(a));
      sub_value_offs_[a].push_back(0);
    }
  }
}

//-----------------------------------------------------------------------------
ufc::finite_element*
FiniteElement::create_sub_element(Array<uint> const& sub_system) const
{
  // Recursively extract sub element
  ufc::finite_element* sub_finite_element = create_sub_element(
      *ufc_finite_element_, sub_system);
  message(2, "Extracted finite element for sub system: %s",
          sub_finite_element->signature());

  return sub_finite_element;
}

//-----------------------------------------------------------------------------
ufc::finite_element*
FiniteElement::create_sub_element(const ufc::finite_element& finite_element,
                                  Array<uint> const& sub_system)
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

}

