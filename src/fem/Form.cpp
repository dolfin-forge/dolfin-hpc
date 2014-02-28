// Copyright (C) 2007 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Aurélien Larcher, 2014.
//
// First added:  2007-12-10
// Last changed: 2014-02-26

#include <dolfin/fem/Form.h>

#include <dolfin/fem/FiniteElement.h>

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
uint Form::coefficient_number(std::string const& name) const
{
  error("Not implemented without UFL support: \n"
        "uint Form::coefficient_number(const std::string& name) const");
  return 0;
}

//-----------------------------------------------------------------------------
std::string Form::coefficient_name(uint const i) const
{
  error("Not implemented without UFL support: \n"
        "std::string Form::coefficient_name(dolfin::uint i) const");
  return "";
}

//-----------------------------------------------------------------------------
bool Form::check(Array<Function*> const& coefficients) const
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
    if (coefficients[i] == NULL)
    {
      error("Got NULL Function as coefficient %d.", i);
    }

    ufc::finite_element* fe = this->create_finite_element(i + this->rank());
    uint r = coefficients[i]->rank();
    uint fe_r = fe->value_rank();
    if (fe_r != r)
    {
      error("Invalid value rank of Function %d, got %d but expecting %d.",
            "You may need to provide the rank of a user defined Function.", i,
            r, fe_r);
    }

    for (uint j = 0; j < r; ++j)
    {
      uint dim = coefficients[i]->dim(j);
      uint fe_dim = fe->value_dimension(j);
      if (dim != fe_dim)
      {
        error(
            "Invalid value dimension %d of Function %d, got %d but expecting %d.",
            "You may need to provide the dimension of a user defined Function.",
            j, i, dim, fe_dim);
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

}
