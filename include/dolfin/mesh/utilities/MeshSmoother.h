// Copyright (C) 2014 Aurélien Larcher
// Licensed under the GNU LGPL Version 2.1.

// Refactoring of classes from UNICORN.

#ifndef DOLFIN_MESH_SMOOTHER_H
#define DOLFIN_MESH_SMOOTHER_H

#include <dolfin/mesh/MeshDependent.h>
#include <dolfin/mesh/MeshValues.h>
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

#define RM( row, col, nrow ) ( ( row ) + ( ( nrow ) * ( col ) ) )

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

  virtual void smooth( MeshValues< bool, Cell > &   smoothed_cells,
                       MeshValues< bool, Vertex > & masked_vertices,
                       MeshValues< real, Cell > &   h0,
                       GenericVector &              node_values,
                       GenericVector &              motion,
                       bool                         reset ) = 0;

  //---------------------------------------------------------------------------

  static void maph0( MeshValues< int, Cell > &  cell_map,
                     MeshValues< real, Cell > & h0,
                     MeshValues< real, Cell > & subh0 );

  ///
  static auto onBoundary( Cell & cell ) -> bool;

  ///
  static void worstElement( int &                      index,
                            MeshValues< bool, Cell > & masked_cells );

  ///
  static void elementNhood( Cell &                     element,
                            MeshValues< bool, Cell > & elements,
                            int                        depth );

  ///
  static void submesh( Mesh &                      sub,
                       MeshValues< bool, Cell > &  smoothed_cells,
                       MeshValues< int, Vertex > & old2new_vertex,
                       MeshValues< int, Cell > &   old2new_cell );

protected:
  /// Constructor
  MeshSmoother( Mesh & mesh );

  ~MeshSmoother() override = default;
};

//-----------------------------------------------------------------------------
}

#endif /* DOLFIN_MESH_SMOOTHER */
