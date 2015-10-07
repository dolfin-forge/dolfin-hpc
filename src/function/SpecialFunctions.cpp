// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Kristian B. Oelgaard, 2007, 2008.
// Modified by Martin Sandve Alnes, 2008.
// Modified by Garth N. Wells, 2008.
// Modified by Niclas Jansson, 2008-2013.
//
// First added:  2008-07-17
// Last changed: 2013-10-22

#include <dolfin/common/constants.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/fem/Form.h>
#include <dolfin/fem/UFC.h>
#include <dolfin/function/SpecialFunctions.h>
#include <dolfin/main/MPI.h>


using namespace dolfin;

//-----------------------------------------------------------------------------
MeshSize::MeshSize(Mesh& mesh) :
  Function(mesh)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void MeshSize::eval(real * values, real const * x) const
{
  values[0] = cell().diameter();
}
//-----------------------------------------------------------------------------
dolfin::uint MeshSize::rank() const
{
  return 0;
}
//-----------------------------------------------------------------------------
dolfin::uint MeshSize::dim(uint i) const
{
  return 1;
}
//-----------------------------------------------------------------------------
real MeshSize::min() const
{
  CellIterator c(mesh());
  real hmin = c->diameter();
  for (; !c.end(); ++c)
    hmin = std::min(hmin, c->diameter());

#ifdef HAVE_MPI
  // Compute the global minimum
  if (mesh().is_distributed())
  {
    real hmin_tmp = hmin;
    MPI_Allreduce(&hmin_tmp, &hmin, 1, MPI_DOUBLE, MPI_MIN, MPI::DOLFIN_COMM);
  }
#endif

  return hmin;
}
//-----------------------------------------------------------------------------
real MeshSize::max() const
{
  CellIterator c(mesh());
  real hmax = c->diameter();
  for (; !c.end(); ++c)
    hmax = std::max(hmax, c->diameter());

#ifdef HAVE_MPI
  // Compute the global maximum
  if (mesh().is_distributed())
  {
    real hmax_tmp = hmax;
    MPI_Allreduce(&hmax_tmp, &hmax, 1, MPI_DOUBLE, MPI_MAX, MPI::DOLFIN_COMM);
  }
#endif

  return hmax;
}
//-----------------------------------------------------------------------------
InvMeshSize::InvMeshSize(Mesh& mesh) :
  Function(mesh)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void InvMeshSize::eval(real * values, real const * x) const
{
  values[0] = 1.0 / cell().diameter();
}
//-----------------------------------------------------------------------------
dolfin::uint InvMeshSize::rank() const
{
  return 0;
}
//-----------------------------------------------------------------------------
dolfin::uint InvMeshSize::dim(uint i) const
{
  return 1;
}
//-----------------------------------------------------------------------------
CellVolume::CellVolume(Mesh& mesh) :
  Function(mesh)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void CellVolume::eval(real * values, real const * x) const
{
  values[0] = cell().volume();
}
//-----------------------------------------------------------------------------
dolfin::uint CellVolume::rank() const
{
  return 0;
}
//-----------------------------------------------------------------------------
dolfin::uint CellVolume::dim(uint i) const
{
  return 1;
}
//-----------------------------------------------------------------------------
real CellVolume::min() const
{
  CellIterator c(mesh());
  real hmin = c->volume();
  for (; !c.end(); ++c)
    hmin = std::min(hmin, c->volume());

#ifdef HAVE_MPI
  // Compute the global minimum
  if (mesh().is_distributed())
  {
    real hmin_tmp = hmin;
    MPI_Allreduce(&hmin_tmp, &hmin, 1, MPI_DOUBLE, MPI_MIN, MPI::DOLFIN_COMM);
  }
#endif

  return hmin;
}
//-----------------------------------------------------------------------------
real CellVolume::max() const
{
  CellIterator c(mesh());
  real hmax = c->volume();
  for (; !c.end(); ++c)
    hmax = std::max(hmax, c->volume());

#ifdef HAVE_MPI
  // Compute the global maximum
  if (mesh().is_distributed())
  {
    real hmax_tmp = hmax;
    MPI_Allreduce(&hmax_tmp, &hmax, 1, MPI_DOUBLE, MPI_MAX, MPI::DOLFIN_COMM);
  }
#endif

  return hmax;
}
//-----------------------------------------------------------------------------
InvCellVolume::InvCellVolume(Mesh& mesh) :
  Function(mesh)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void InvCellVolume::eval(real * values, real const * x) const
{
  values[0] = 1.0 / cell().volume();
}
//-----------------------------------------------------------------------------
dolfin::uint InvCellVolume::rank() const
{
  return 0;
}
//-----------------------------------------------------------------------------
dolfin::uint InvCellVolume::dim(uint i) const
{
  return 1;
}
//-----------------------------------------------------------------------------
AvgMeshSize::AvgMeshSize(Mesh& mesh) :
  Function(mesh)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void AvgMeshSize::eval(real * values, real const * x) const
{
  // If there is no facet (assembling on interior), return cell diameter
  if (facet() < 0)
    values[0] = cell().diameter();
  else
  {
    // Create facet from the global facet number
    Facet facet0(mesh(),
                 cell().entities(cell().mesh().topology().dim() - 1)[facet()]);

    // If there are two cells connected to the facet
    if (facet0.numEntities(cell().mesh().topology().dim()) == 2)
    {
      // Create the two connected cells and return the average of their diameter
      Cell cell0(mesh(), facet0.entities(cell().mesh().topology().dim())[0]);
      Cell cell1(mesh(), facet0.entities(cell().mesh().topology().dim())[1]);

      values[0] = (cell0.diameter() + cell1.diameter()) / 2.0;
    }
    // Else there is only one cell connected to the facet and the average is
    // the cell diameter
    else
      values[0] = cell().diameter();
  }
}
//-----------------------------------------------------------------------------
dolfin::uint AvgMeshSize::rank() const
{
  return 0;
}
//-----------------------------------------------------------------------------
dolfin::uint AvgMeshSize::dim(uint i) const
{
  return 1;
}
//-----------------------------------------------------------------------------
FacetNormal::FacetNormal(Mesh& mesh) :
  Function(mesh)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void FacetNormal::eval(real * values, real const * x) const
{
  if (facet() >= 0)
  {
    Point n = cell().normal(facet());
    for (uint i = 0; i < cell().dim(); ++i)
    {
      values[i] = n[i];
    }
  }
  else
  {
    for (uint i = 0; i < cell().dim(); ++i)
    {
      values[i] = 0.0;
    }
  }
}
//-----------------------------------------------------------------------------
dolfin::uint FacetNormal::rank() const
{
  return 1;
}
//-----------------------------------------------------------------------------
dolfin::uint FacetNormal::dim(uint i) const
{
  if (i > 0)
    error("Invalid dimension %d in FacetNormal::dim.", i);
  return mesh().geometry().dim();
}
//-----------------------------------------------------------------------------
FacetArea::FacetArea(Mesh& mesh) :
  Function(mesh)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void FacetArea::eval(real * values, real const * x) const
{
  if (facet() >= 0)
    values[0] = cell().facetArea(facet());
  else
    values[0] = 0.0;
}
//-----------------------------------------------------------------------------
dolfin::uint FacetArea::rank() const
{
  return 0;
}
//-----------------------------------------------------------------------------
dolfin::uint FacetArea::dim(uint i) const
{
  return 1;
}
//-----------------------------------------------------------------------------
InvFacetArea::InvFacetArea(Mesh& mesh) :
  Function(mesh)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void InvFacetArea::eval(real * values, real const * x) const
{
  if (facet() >= 0)
    values[0] = 1.0 / cell().facetArea(facet());
  else
    values[0] = 0.0;
}
//-----------------------------------------------------------------------------
dolfin::uint InvFacetArea::rank() const
{
  return 0;
}
//-----------------------------------------------------------------------------
dolfin::uint InvFacetArea::dim(uint i) const
{
  return 1;
}
//-----------------------------------------------------------------------------
OutflowFacet::OutflowFacet(Mesh& mesh, Form& form) :
  Function(mesh),
  mesh(mesh),
  form(form)
{
  // Some simple sanity checks on form
  if (!(form.rank() == 0 && form.num_coefficients() == 2))
  {
    error("Invalid form: rank = %d, number of coefficients = %d."
          "Must be rank 0 form with 2 coefficients.",
          form.rank(), form.num_coefficients());
  }
  if (!(form.num_cell_integrals() == 0
      && form.num_exterior_facet_integrals() == 1
      && form.num_interior_facet_integrals() == 0))
  {
    error("Invalid form: Must have exactly 1 exterior facet integral");
  }

  form.update_dofmaps();
  ufc = new UFC(form, mesh, form.dofmaps());
}
//-----------------------------------------------------------------------------
OutflowFacet::~OutflowFacet()
{
  delete ufc;
}
//-----------------------------------------------------------------------------
void OutflowFacet::eval(real * values, real const * x) const
{
  // If there is no facet (assembling on interior), return 0.0
  if (facet() < 0)
  {
    values[0] = 0.0;
  }
  else
  {
    // Copy cell, cannot call interpolate with const cell()
    Cell cell0(cell());
    ufc->update(cell0, mesh.distdata());

    // Interpolate coefficients on cell and current facet
    for (dolfin::uint i = 0; i < form.coefficients().size(); i++)
    {
      form.coefficients()[i]->interpolate(ufc->w[i], ufc->cell,
                                          *ufc->coefficient_elements[i],
                                          cell0, facet());
    }

    // Get exterior facet integral (we need to be able to tabulate ALL facets
    // of a given cell)
    ufc::exterior_facet_integral* integral = ufc->exterior_facet_integrals[0];

    // Call tabulate_tensor on exterior facet integral,
    // dot(velocity, facet_normal)
    integral->tabulate_tensor(ufc->A, ufc->w, ufc->cell, facet());
  }

  // If dot product is positive, the current facet is an outflow facet
  if (ufc->A[0] > DOLFIN_EPS)
  {
    values[0] = 1.0;
  }
  else
  {
    values[0] = 0.0;
  }
}
//-----------------------------------------------------------------------------
dolfin::uint OutflowFacet::rank() const
{
  return 0;
}
//-----------------------------------------------------------------------------
dolfin::uint OutflowFacet::dim(uint i) const
{
  return 1;
}
//-----------------------------------------------------------------------------

