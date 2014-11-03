// Copyright (C) 2014 Aurélien Larcher
// Licensed under the GNU LGPL Version 2.1.
//
// Refactoring of classes from UNICORN.
//
// First added:  2014-06-12
// Last changed: 2014-06-12

#ifndef DOLFIN_MESH_SMOOTHER_H
#define DOLFIN_MESH_SMOOTHER_H

#include <dolfin/mesh/MeshDependent.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/SubDomain.h>

#if HAVE_SUNPERF_H
#include <sunperf.h>
#elif HAVE_SCSL_CBLAS_H
#include <cmplrs/cblas.h>
#elif HAVE_GSL_CBLAS_H
extern "C"
{
#include <gsl_cblas.h>
}
#elif HAVE_CBLAS_H
extern "C"
{
#include <cblas.h>
}
#endif

#define RM(row,col,nrow) ((row) + ((nrow)*(col)))

namespace dolfin
{

class Cell;
class Mesh;
class Vector;

//-----------------------------------------------------------------------------
class MeshSmoother : public MeshDependent
{
public:

  //--- INTERFACE -------------------------------------------------------------

  virtual void smooth(MeshFunction<bool>& smoothed_cells,
                      MeshFunction<bool>& masked_vertices,
                      MeshFunction<real>& h0, GenericVector& node_values,
                      GenericVector& motion, bool reset) = 0;

  //---------------------------------------------------------------------------

  static void maph0(Mesh& mesh, Mesh& sub, MeshFunction<int>& cell_map,
                    MeshFunction<real>& h0, MeshFunction<real>& subh0);

  ///
  static bool onBoundary(Cell& cell);

  ///
  static void worstElement(Mesh& mesh, int& index,
                           MeshFunction<bool>& masked_cells);

  ///
  static void elementNhood(Mesh& mesh, Cell& element,
                           MeshFunction<bool>& elements, int depth);

  ///
  static void submesh(Mesh& mesh, Mesh& sub, MeshFunction<bool>& smoothed_cells,
                      MeshFunction<int>& old2new_vertex,
                      MeshFunction<int>& old2new_cell);

protected:

  // Sub domain for Dirichlet boundary condition
  class DirichletBoundary : public SubDomain
  {

  public:

    bool inside(const real* x, bool on_boundary) const
    {
      return on_boundary;
    }
  };

  /// Constructor
  MeshSmoother(Mesh& mesh);

  virtual ~MeshSmoother();

};

//-----------------------------------------------------------------------------
}

#endif /* DOLFIN_MESH_SMOOTHER */
