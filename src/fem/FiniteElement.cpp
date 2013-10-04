// Copyright (C) 2013 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2013-09-12
// Last changed: 2013-09-12

#include <dolfin/elements/ElementLibrary.h>
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
    geom_dim_(0),
    degree_(0)

{
  init();
}

//-----------------------------------------------------------------------------
FiniteElement::FiniteElement(Mesh& mesh, Form& form, uint i) :
    ufc_finite_element_(NULL),
    finite_element_local_(true),
    sub_value_dims_(NULL),
    topo_dim_(0),
    geom_dim_(0),
    degree_(-1)
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

  init();
}

//-----------------------------------------------------------------------------
FiniteElement::FiniteElement(ufc::finite_element& finite_element,
                             bool const finite_element_local) :
    ufc_finite_element_(&finite_element),
    finite_element_local_(finite_element_local),
    sub_value_dims_(NULL),
    topo_dim_(0),
    geom_dim_(0),
    degree_(0)
{
  init();
}

FiniteElement::~FiniteElement()

{
  if (finite_element_local_)
    delete ufc_finite_element_;

  if (sub_value_dims_)
    delete[] sub_value_dims_;

  if (sub_value_offs_)
    delete[] sub_value_offs_;
}

//-----------------------------------------------------------------------------
void FiniteElement::init()
{
  dolfin_assert(ufc_finite_element_);

  // Set attributes
  FE::attributes const attr = ElementLibrary::get_attributes(ufc_finite_element_->signature());
  type_ = attr.type;
  family_ = attr.family;
  degree_ = attr.degree;

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
Array<uint> const&
FiniteElement::sub_value_dimensions(uint i) const
{
  return sub_value_dims_[i];
}
//-----------------------------------------------------------------------------
Array<uint> const&
FiniteElement::sub_value_offsets(uint i) const
{
  return sub_value_offs_[i];
}
//-----------------------------------------------------------------------------
void FiniteElement::info() const
{
  std::stringstream msg;
  uint const padding = 24;
  msg << std::endl;
  msg << std::setw(padding) << "signature = " << signature() << std::endl;
  std::string shape;
  switch (ufc_finite_element_->cell_shape())
  {
    case ufc::interval:
      shape = "";
      break;
    case ufc::triangle:
      shape = "triangle";
      break;
    case ufc::tetrahedron:
      shape = "tetrahedron";
      break;
    default:
      shape = "unknown";
      break;
  }
  msg << std::setw(padding) << "cell_shape = " << shape << std::endl;
  msg << std::setw(padding) << "topological_dimension = "
      << topological_dimension() << std::endl;
  msg << std::setw(padding) << "geometric_dimension = " << geometric_dimension()
      << std::endl;
  msg << std::setw(padding) << "space_dimension = " << space_dimension()
      << std::endl;
  msg << std::setw(padding) << "value_rank = " << value_rank() << std::endl;
  msg << std::setw(padding) << "value_dimension = " << std::endl;
  for (uint a = 0; a < geom_dim_; ++a)
  {
    msg << std::setw(padding + 1) << " [" << "axis " << a << "] : "
        << value_dimension(a) << std::endl;
  }
  msg << std::setw(padding) << "num_sub_elements = " << num_sub_elements()
      << std::endl;
  msg << std::setw(padding) << "sub_value_dimensions = " << std::endl;
  for (uint a = 0; a < geom_dim_; ++a)
  {
    msg << std::setw(padding + 1) << " [" << "axis " << a << "] : "
        << std::endl;
    for (uint sub = 0; sub < sub_value_dimensions(a).size(); ++sub)
    {
      msg << std::setw(padding + 12) << " [" << "element " << sub << "] : "
          << sub_value_dimensions(a)[sub] << " ( +" << sub_value_offsets(a)[sub]
          << ")" << std::endl;
    }
  }
  msg << std::setw(padding) << "degree = " << degree_ << std::endl;
  std::cout << msg.str();
}

}

