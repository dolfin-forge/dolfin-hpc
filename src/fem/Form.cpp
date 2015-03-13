// Copyright (C) 2007 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Aurélien Larcher, 2014.
//
// First added:  2007-12-10
// Last changed: 2014-02-26

#include <dolfin/fem/Form.h>

#include <dolfin/fem/CoefficientMap.h>
#include <dolfin/fem/FiniteElement.h>
#include <dolfin/fem/FiniteElementSpace.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
Form::Form(Mesh& mesh) :
    mesh_(mesh),
    dof_map_set_(*this, mesh)
{
}

//-----------------------------------------------------------------------------
Form::~Form()
{
}

//-----------------------------------------------------------------------------
void Form::update_dofmaps() const
{
  if (dof_map_set_.size() == 0)
  {
    dof_map_set_.update(*this, mesh_);
  }
}

//-----------------------------------------------------------------------------
DofMapSet& Form::dofmaps() const
{
  this->update_dofmaps();
  return dof_map_set_;
}

//-----------------------------------------------------------------------------
uint Form::coefficient_index(std::string const& name) const
{
  error("Not implemented without UFL support: \n"
          "uint Form::coefficient_index(const std::string& name) const");
  return 0;
}

//-----------------------------------------------------------------------------
uint Form::coefficient_number(std::string const& name) const
{
  return this->coefficient_index(name);
}

//-----------------------------------------------------------------------------
std::string Form::coefficient_name(uint i) const
{
  error("Not implemented without UFL support: \n"
        "std::string Form::coefficient_name(dolfin::uint i) const");
  return "";
}

//----------------------------------------------------------------------------
FiniteElementSpace * Form::create_space(uint i) const
{
  ufc::finite_element * test_f = this->form().create_finite_element(i);
  ufc::dofmap * test_d = this->form().create_dofmap(i);
  // For an argument the mesh is the one passed to the form and for coefficient
  // the mesh passed to the function.
  Mesh& mesh = (
      i < this->rank() ? this->mesh() : this->coefficients()[i - this->rank()]->mesh());
  return new FiniteElementSpace(mesh, *test_f, *test_d, true);
}

//----------------------------------------------------------------------------
FiniteElementSpace * Form::create_coefficient_space(
    std::string const& name) const
{
  return this->create_space(this->rank() + this->coefficient_index(name));
}

//-----------------------------------------------------------------------------
bool Form::check_coefficients(Array<Function*> const& coefficients) const
{
  // Check that we get the correct number of coefficients
  if (coefficients.size() != this->num_coefficients())
  {
    error("Incorrect number of coefficients: %d given but %d required.",
          coefficients.size(), this->num_coefficients());
  }

  // Check that all coefficients have valid value dimensions
  for (uint i = 0; i < coefficients.size(); ++i)
  {
    message(1, "Checking coefficient %d:", i);
    if (coefficients[i] == NULL)
    {
      error("Got NULL Function as coefficient %d labeled as '%s'.", i,
            this->coefficient_name(i).c_str());
    }

    if (coefficients[i]->type() == Function::empty)
    {
      error("Got empty Function as coefficient %d labeled as '%s'.", i,
            this->coefficient_name(i).c_str());
    }

    ufc::finite_element * fe = this->create_finite_element(i + this->rank());
    uint coef_rank = coefficients[i]->rank();
    uint fe_rank = fe->value_rank();
    message(1, "Coefficient rank: expected  = %d, provided = %d, ", fe_rank,
            coef_rank);
    if (fe_rank != coef_rank)
    {
      error("Invalid value rank of Function '%s' with index %d:\n"
            "Got %d but expecting %d.\n"
            "You may need to provide the rank of a user defined Function.",
            this->coefficient_name(i).c_str(), i, coef_rank, fe_rank);
    }

    for (uint j = 0; j < coef_rank; ++j)
    {
      uint dim = coefficients[i]->dim(j);
      uint fe_dim = fe->value_dimension(j);
      if (dim != fe_dim)
      {
        error(
            "Invalid value dimension %d of Function '%s' with index %d:\n"
            "got %d but expecting %d.\n"
            "You may need to provide the dimension of a user defined Function.",
            j, this->coefficient_name(i).c_str(), i, dim, fe_dim);
      }
    }
    delete fe;
  }

  // Check that the cell dimension matches the mesh dimension
  if (this->rank() + this->num_coefficients() > 0)
  {
    ufc::finite_element* element = this->create_finite_element(0);
    dolfin_assert(element);
    CellType::Type celltype = mesh().type().cellType();
    ufc::shape shape = element->cell_shape();

    if (celltype == CellType::interval && shape != ufc::interval)
    {
      error("Mesh cell type (intervals) does not match cell type of form.");
    }
    if (celltype == CellType::triangle && shape != ufc::triangle)
    {
      error("Mesh cell type (triangles) does not match cell type of form.");
    }
    if (celltype == CellType::tetrahedron && shape != ufc::tetrahedron)
    {
      error("Mesh cell type (tetrahedra) does not match cell type of form.");
    }
    delete element;
  }
  return true;
}

//-----------------------------------------------------------------------------
bool Form::check_index(uint i) const
{
  // Check argument
  uint const num_arguments = form().rank() + form().num_coefficients();
  if (i >= num_arguments)
  {
    error("Illegal function index %d. Form only has %d arguments.", i,
          num_arguments);
  }
  return true;
}

//-----------------------------------------------------------------------------
void Form::assign_coefficients(CoefficientMap const& coefficient_map,
                               Array<dolfin::Function *>& form_coefficients)
{
  form_coefficients.clear();
  for (uint i = 0; i < this->num_coefficients(); ++i)
  {
    std::string name = this->coefficient_name(i);
    if(coefficient_map.has(name))
    {
      form_coefficients.push_back(coefficient_map.get(name));
    }
    else
    {
      error("Missing coefficient named '%s' in CoefficientMap.", name.c_str());
    }
  }
}

}
