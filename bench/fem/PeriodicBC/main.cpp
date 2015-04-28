#include "PoissonP1.h"

#define DEBUG 1

#include <dolfin/math/basic.h>
#include <dolfin/fem/Assembler.h>
#include <dolfin/fem/FiniteElementSpace.h>
#include <dolfin/fem/PeriodicDofsMapping.h>
#include <dolfin/fem/ScratchSpace.h>
#include <dolfin/fem/SparsityPatternBuilder.h>
#include <dolfin/fem/SubSystem.h>
#include <dolfin/fem/UFC.h>
#include <dolfin/la/Matrix.h>
#include <dolfin/la/PETScMatrix.h>
#include <dolfin/la/SparsityPattern.h>
#include <dolfin/la/Vector.h>
#include <dolfin/main/MPI.h>
#include <dolfin/mesh/BoundaryMesh.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/IntersectionDetector.h>
#include <dolfin/mesh/MappedManifold.h>
#include <dolfin/mesh/PeriodicSubDomain.h>
#include <dolfin/mesh/UnitCube.h>
#include <dolfin/ufl/UFLFiniteElement.h>

#include <mpi.h>

using namespace dolfin;

real const XMIN = 0.0;
real const XMAX = 1.0;
real const YMIN = 0.0;
real const YMAX = 1.0;
real const ZMIN = 0.0;
real const ZMAX = 1.0;

//-----------------------------------------------------------------------------
class LeftToRight : public PeriodicSubDomain
{

public:

  /// Defines subdomain G which is the left side of the unit square
  bool inside(real const * x, bool const on_boundary) const
  {
    return on_boundary && abscmp(x[0], XMIN, 1.0e-6);
  }

  /// Defines the mapping from the right side (H) to the left side (G)
  void map(real const * xH, real* xG) const
  {
    xG[0] = xH[0] - XMAX;
    xG[1] = xH[1];
  }

};

//-----------------------------------------------------------------------------
class PeriodicSquare : public PeriodicSubDomain
{

public:

  /// Defines subdomain G which is the left side of the unit square
  bool inside(real const * x, bool const on_boundary) const
  {
    return on_boundary && abscmp(x[0], XMIN, 1.0e-6);
  }

  /// Defines the mapping from the right side (H) to the left side (G)
  void map(real const * xH, real* xG) const
  {
    if(close(xH[0], 1.0))
    {
      xG[0] = xH[0] - XMAX;
    }
    if(close(xH[1], 1.0))
    {
      xG[1] = xH[1] - YMAX;
    }
  }

};

//-----------------------------------------------------------------------------
class PeriodicCube : public PeriodicSubDomain
{

public:

  /// Defines subdomain G which is the left side of the unit square
  bool inside(real const * x, bool const on_boundary) const
  {
    return on_boundary && abscmp(x[0], XMIN, 1.0e-6);
  }

  /// Defines the mapping from the right side (H) to the left side (G)
  void map(real const * xH, real* xG) const
  {
    xG[0] = xH[0] - XMAX;
    xG[1] = xH[1] - YMAX;
    xG[2] = xH[2] - ZMAX;
  }

};

//-----------------------------------------------------------------------------
struct MappedDof
{
  uint global_index;
  uint facet_index;
  real coordinates[EuclideanSpace::MAX_DIMENSION];

  MappedDof() :
      global_index(0),
      facet_index(0),
      coordinates()
  {
  }

};

int main(int argc, char *argv[])
{
  {
    Mesh mesh("cube10R.xml.gz");
    LeftToRight subdomain;
    mesh.add_periodic_constraint(subdomain);
    ufl::FiniteElement cg1_3d(ufl::Family::CG, mesh.type(), 1);
    FiniteElementSpace spaceU(mesh, cg1_3d);
    PeriodicDofsMapping const& pdm = spaceU.dofmap().periodic_mapping();
    pdm.disp();
  }

  if(false)
  {
    Mesh mesh("square100R.xml.gz");
    LeftToRight subdomain;
    mesh.add_periodic_constraint(subdomain);

  //
  bool do_test0 = true;

  //
  if(do_test0)
  {
    Matrix Aas;
    as.assemble(Aas, a, true);
    Array<uint> columns;
    Array<real> values;
    Mat matAas = reinterpret_cast<PETScMatrix *>(Aas.instance())->mat();

    for (uint i = 0; i < pdm.num_Gdofs(); ++i)
    {
      uint row = pdm.get_Gindices()[i];
      if (!dmU.is_ghost(row))
      {
        Aas.getrow(row, columns, values);
        message("Row %d with %d entries.", row, columns.size());
      }
    }
  }

  return 0;
  //
  /*


  //
  Mat matAin = reinterpret_cast<PETScMatrix *>(Ain.instance())->mat();
  real * block = new real[pdm.max_local_dimension()];
  for (uint i = 0; i < A.size(0); ++i)
  {
   const int *cols = 0;
   const double *vals = 0;
   int ncols = 0;
   MatGetRow(matA, i, &ncols, &cols, &vals);
   message("%8d : %8d", i, ncols);
   MatRestoreRow(matA, i, &ncols, &cols, &vals);
  }
  delete[] block;

  return 0;
  */

  //A.apply();
  //spattern.disp();
  /*
  if (dolfin::MPI::numProcesses() > 1)
  {
    uint p = dolfin::MPI::processNumber();

    uint local_size = spattern.numLocalRows(p);
    uint* d_nzrow = new uint[local_size];
    uint* o_nzrow = new uint[local_size];
    spattern.numNonZeroPerRow(p, d_nzrow, o_nzrow);
    uint spattern_size_0 = spattern.size(0);
    uint spattern_size_1 = spattern.size(1);
    //const_cast<SparsityPattern&>(spattern).clear();
    //
    delete[] d_nzrow;
    delete[] o_nzrow;
  }
  else
  {
    uint* nzrow = new uint[spattern.size(0)];
      spattern.numNonZeroPerRow(nzrow);
    std::memset(block, 1, pdm.max_local_dimension()*sizeof(real) );
    A.ident(pdm.num_Gdofs(), pdm.get_Gindices());
    A.apply();


    //A.zero();
    message("Reassemble");
    //as.assemble(A, a, false);
    for (uint i = 0; i < A.size(0); ++i)
    {
      const int *cols = 0;
      const double *vals = 0;
      int ncols = 0;
      MatGetRow(matA, i, &ncols, &cols, &vals);
      message("%8d : %8d", i, ncols);
      MatRestoreRow(matA, i, &ncols, &cols, &vals);
    }

    delete[] nzrow;

  }
  */

  //

  //return 0;
}

/*
 // Below is PeriodicBC

 //--- Apply homogeneous Dirichlet BC for all G and I dofs -----------------
 message("Apply homogeneous Dirichlet BC for all G and I dofs");
 uint * GIrows = new uint[Gdofs.size() + Idofs.size()];
 real * GIzeros = new real[Gdofs.size() + Idofs.size()];
 uint GIii = 0;
 for (_set<uint>::const_iterator it = Gdofs.begin(); it != Gdofs.end(); ++it)
 {
 GIrows[GIii] = *it;
 GIzeros[GIii] = 0.0;
 ++GIii;
 }
 for (_set<uint>::const_iterator it = Idofs.begin(); it != Idofs.end(); ++it)
 {
 GIrows[GIii] = *it;
 GIzeros[GIii] = 0.0;
 ++GIii;
 }
 // Apply at I dofs
 A.ident(GIii, GIrows);
 b.set(GIzeros, GIii, GIrows);
 delete[] GIzeros;
 delete[] GIrows;

 //--- Local expansion coefficients ---------------------------------------
 uint const ldofsH_size = ldofsH.size();
 uint const ldofsG_size = ldofsG.size();
 uint coefsii = 0;
 for(uint ii = 0; ii < ldofsH_size; ++ii)
 {
 A.set(&coefsG[ii], 1, &ldofsH[ii*num_facet_dofsU], num_facet_dofsU, &ldofsG[ii*num_facet_dofsU]);
 }
 b.set(&coefsG[0], ldofsG.size(), &ldofsG[0]);

 //--- Compute expansion coefficients for G dofs --------------------------
 message("Compute expansions coefficients for G dofs in terms of H dofs");
 real * Ufacetcount = new real[Gdofs2facets_.size()];
 uint * Ufacetdofs = new uint[Gdofs2facets_.size()];
 uint ii = 0;
 for (_map<uint, Array<uint> >::const_iterator it = Gdofs2facets_.begin();
 it != Gdofs2facets_.end(); ++it)
 {
 uint const dof = it->first;
 Array<uint> const& facets = it->second;
 Ufacetcount[ii] = facets.size();
 Ufacetdofs[ii] = dof;
 ++ii;
 //for (Array<uint>::const_iterator f = facets.begin(); f != facets.end();
 //     ++f)
 //{
 //  Facet facet(mesh, *f);
 //  Cell(mesh, facet.entities(tdim)[0]);

 //}
 }
 FacetCount.vector().set(Ufacetcount, Gdofs2facets_.size(), Ufacetdofs);
 FacetCount.vector().apply();
 File facetcount("facetcount.pvd");
 facetcount << FacetCount;
 delete[] Ufacetcount;
 delete[] Ufacetdofs;

 //
 delete[] facetdofs_value;
 delete[] Gdofsx;
 delete[] Gdofsi;

 // Apply changes
 A.apply();
 b.apply();

 File fB("RHS.pvd");
 fB << B;

 }

 return 0;
 }
 */
