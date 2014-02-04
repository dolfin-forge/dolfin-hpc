// Copyright (C) 2007 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2007-12-10
// Last changed:

#include <dolfin/fem/Form.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
Form::Form() :
    dof_map_set_(NULL)
{
}

//-----------------------------------------------------------------------------
Form::~Form()
{
  delete dof_map_set_;
}

//-----------------------------------------------------------------------------
void Form::update_dofmaps(Mesh& mesh) const
{
  if (!dof_map_set_)
  {
    // Create dof maps
    dof_map_set_ = new DofMapSet(*this, mesh);
  }
  else if (mesh.hash() != dof_map_set_->mesh().hash())
  {
    error("Attempt to reinitialize the form's dof map set with a different mesh"
          " topology.");
  }
}

//-----------------------------------------------------------------------------
DofMapSet& Form::dofmaps() const
{
  if (!dof_map_set_)
    error("Degree of freedom maps for Form have not been created.");

  return *dof_map_set_;
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

}
