#include "PoissonP1.h"

#define DEBUG 1

#include <dolfin/common/Test.h>
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
#include <dolfin/math/basic.h>
#include <dolfin/main/MPI.h>
#include <dolfin/mesh/BoundaryMesh.h>
#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/IntersectionDetector.h>
#include <dolfin/mesh/MappedManifold.h>
#include <dolfin/mesh/PeriodicSubDomain.h>
#include <dolfin/mesh/UnitCube.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/ufl/UFLFiniteElement.h>

#include <mpi.h>

#include <cstring>

using namespace dolfin;

real const XMIN = 0.0;
real const XMAX = 1.0;
real const YMIN = 0.0;
real const YMAX = 1.0;
real const ZMIN = 0.0;
real const ZMAX = 1.0;

//-----------------------------------------------------------------------------
class LeftToRight2D : public PeriodicSubDomain
{

public:

  /// Defines subdomain G which is the left side of the unit square
  bool inside(real const * x, bool const on_boundary) const
  {
    return on_boundary && close(x[0], XMIN);
  }

  /// Defines the mapping from the right side (H) to the left side (G)
  void map(real const * xH, real* xG) const
  {
    xG[0] = xH[0] - XMAX;
    xG[1] = xH[1];
  }

};

//-----------------------------------------------------------------------------
class LeftToRight3D : public PeriodicSubDomain
{

public:

  /// Defines subdomain G which is the left side of the unit square
  bool inside(real const * x, bool const on_boundary) const
  {
    return on_boundary && close(x[0], XMIN);
  }

  /// Defines the mapping from the right side (H) to the left side (G)
  void map(real const * xH, real* xG) const
  {
    xG[0] = xH[0] - XMAX;
    xG[1] = xH[1];
    xG[2] = xH[2];
  }

};

//-----------------------------------------------------------------------------
class PeriodicSquareX : public PeriodicSubDomain
{

public:

  /// Defines subdomain G which is the left side of the unit square
  bool inside(real const * x, bool const on_boundary) const
  {
    bool on_faceX0 = (close(x[0], XMIN)
                      && x[1] >= (YMIN - BMARG) && x[1] <= (YMAX + BMARG));
    return on_boundary && (on_faceX0);
  }

  /// Defines the mapping from the right side (H) to the left side (G)
  void map(real const * xH, real* xG) const
  {
    xG[0] = xH[0] - XMAX;
    xG[1] = xH[1];
  }
};

//-----------------------------------------------------------------------------
class PeriodicSquareY : public PeriodicSubDomain
{

public:

  /// Defines subdomain G which is the left side of the unit square
  bool inside(real const * x, bool const on_boundary) const
  {
    bool on_faceY0 = (close(x[1], YMIN)
                      && x[0] >= (XMIN - BMARG) && x[0] <= (XMAX + BMARG));
    return on_boundary && (on_faceY0);
  }

  /// Defines the mapping from the right side (H) to the left side (G)
  void map(real const * xH, real* xG) const
  {
    xG[0] = xH[0];
    xG[1] = xH[1] - YMAX;
  }
};

//-----------------------------------------------------------------------------
class PeriodicCubeX : public PeriodicSubDomain
{

public:

  ///
  bool inside(real const * x, bool const on_boundary) const
  {
    bool on_faceX0 = (close(x[0], XMIN)
        && x[1] >= (YMIN - BMARG) && x[1] <= (YMAX + BMARG)
        && x[2] >= (ZMIN - BMARG) && x[2] <= (ZMAX + BMARG));
    return on_boundary && (on_faceX0);
  }

  /// Defines the mapping from the right side (H) to the left side (G)
  void map(real const * xH, real* xG) const
  {
      xG[0] = xH[0] - XMAX;
      xG[1] = xH[1];
      xG[2] = xH[2];
  }
};

//-----------------------------------------------------------------------------
class PeriodicCubeY : public PeriodicSubDomain
{

public:

  ///
  bool inside(real const * x, bool const on_boundary) const
  {
    bool on_faceY0 = (close(x[1], YMIN) && x[2] >= (ZMIN - BMARG)
        && x[2] <= (ZMAX + BMARG) && x[0] >= (XMIN - BMARG)
        && x[0] <= (XMAX + BMARG));
    return on_boundary && (on_faceY0);
  }

  /// Defines the mapping from the right side (H) to the left side (G)
  void map(real const * xH, real* xG) const
  {
      xG[0] = xH[0];
      xG[1] = xH[1] - YMAX;
      xG[2] = xH[2];
  }
};

//-----------------------------------------------------------------------------
class PeriodicCubeZ : public PeriodicSubDomain
{

public:

  ///
  bool inside(real const * x, bool const on_boundary) const
  {
    bool on_faceZ0 = (close(x[2], ZMIN) && x[0] >= (XMIN - BMARG)
        && x[0] <= (XMAX + BMARG) && x[1] >= (YMIN - BMARG)
        && x[1] <= (YMAX + BMARG));
    return on_boundary && (on_faceZ0);
  }

  /// Defines the mapping from the right side (H) to the left side (G)
  void map(real const * xH, real* xG) const
  {
    xG[0] = xH[0];
    xG[1] = xH[1];
    xG[2] = xH[2] - ZMAX;
  }
};

//-----------------------------------------------------------------------------
void write_facets(std::string name, Mesh& mesh, std::string set,
                  _set<uint> const& S)
{
  message("Number of " + set + " facets = %d", S.size());
  std::stringstream ss;
  ss << name << "_" + set + "facets" << dolfin::MPI::size() << "P.pvd";
  MeshValues<bool, Vertex> mf(mesh);
  mesh.init(mesh.type().facet_dim(), 0);
  for(_set<uint>::const_iterator it = S.begin(); it != S.end(); ++it)
  {
    Facet f(mesh, *it);
    for(VertexIterator v(f); ! v.end(); ++v)
    {
      mf.set(*v, true);
    }
  }
  File f(ss.str());
  f << mf;
}

//-----------------------------------------------------------------------------
void write_mesh(std::string name, Mesh& mesh)
{
  std::stringstream ss0;
  ss0 << name << "_mesh_" << dolfin::MPI::size() << "P.pvd";
  File f0(ss0.str());
  f0 << mesh;

  for(uint i = 0; i < mesh.periodic_mappings().size(); ++i)
  {
    std::stringstream ssi;
    ssi << i;
    write_facets(name, mesh, "G"+ssi.str(), mesh.periodic_mappings()[i]->Gfacets());
    write_facets(name, mesh, "H"+ssi.str(), mesh.periodic_mappings()[i]->Hfacets());
    write_facets(name, mesh, "I"+ssi.str(), mesh.periodic_mappings()[i]->Ifacets());

    std::stringstream ss1;
    ss1 << name << "_Gmesh" << ssi.str() << "_" << dolfin::MPI::size()
        << "P_" << dolfin::MPI::rank() << ".pvd";
    File f1(ss1.str());
    f1 << *mesh.periodic_mappings()[i];
  }
}

//-----------------------------------------------------------------------------
void write_Gdofs(std::string name, FiniteElementSpace const& space)
{
  message("Write Gdofs for test %s", name.c_str());
  Function G(space);
  PeriodicDofsMapping const& pdm = G.space().dofmap().periodic_mapping();
  //pdm.disp();
  real * blockG = new real[pdm.num_Gdofs()];
  message("Number of G dofs = %d", pdm.num_Gdofs());
  uint count = 0;
  uint Gdof = 0;
  uint * Hdofs = new uint[pdm.max_local_dimension()];
  for(uint ii = 0; ii < pdm.num_Gdofs(); ++ii)
  {
    pdm.tabulate_dofs(ii, &Gdof, Hdofs, count);
    blockG[ii] = count;
  }
  delete[] Hdofs;
  G.vector() = 0.0;
  G.vector().set(blockG, pdm.num_Gdofs(), pdm.get_Gindices());
  G.sync();
  delete[] blockG;
  std::stringstream ss;
  ss << name << "_G_" << dolfin::MPI::size() << "P.pvd";
  File fG(ss.str());
  fG << G;
  message("Done");
}

//-----------------------------------------------------------------------------
bool onEntity(real* coordinates, MeshEntity& entity)
{
  // Check if the coordinates are on the same line as the line segment
  if ( entity.dim() == 1 )
  {
    // Create points
    Point p(coordinates[0], coordinates[1]);
    Point v0 = Vertex(entity.mesh(), entity.entities(0)[0]).point();
    Point v1 = Vertex(entity.mesh(), entity.entities(0)[1]).point();

    // Create vectors
    Point v01 = v1 - v0;
    Point vp0 = v0 - p;
    Point vp1 = v1 - p;

    // Check if the length of the sum of the two line segments vp0 and vp1 is
    // equal to the total length of the entity
    if ( std::abs(v01.norm() - vp0.norm() - vp1.norm()) < DOLFIN_EPS )
    {
      return true;
    }
    else
    {
      return false;
    }
  }
  // Check if the coordinates are in the same plane as the triangular entity
  else if ( entity.dim() == 2 )
  {
    // Create points
    Point p(coordinates[0], coordinates[1], coordinates[2]);
    Point v0 = Vertex(entity.mesh(), entity.entities(0)[0]).point();
    Point v1 = Vertex(entity.mesh(), entity.entities(0)[1]).point();
    Point v2 = Vertex(entity.mesh(), entity.entities(0)[2]).point();

    // Create vectors
    Point v01 = v1 - v0;
    Point v02 = v2 - v0;
    Point vp0 = v0 - p;
    Point vp1 = v1 - p;
    Point vp2 = v2 - p;

    // Check if the sum of the area of the sub triangles is equal to the total
    // area of the entity
    if ( std::abs(v01.cross(v02).norm() - vp0.cross(vp1).norm() - vp1.cross(vp2).norm()
        - vp2.cross(vp0).norm()) < DOLFIN_EPS )
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  error("Unable to determine if given point is on entity (not implemented for given facet dimension).");

  return false;
}

//-----------------------------------------------------------------------------
void check_matching_facets(MappedManifold& mm, Point& p)
{
  Array<uint> matching;
  mm.intersector().overlap(p, matching);
  message("Matching with point of plane : %d", matching.size());
  std::set<uint> matching_set;
  matching_set.insert(matching.begin(), matching.end());
  message("Unique   with point of plane : %d", matching_set.size());
  for(std::set<uint>::const_iterator it = matching_set.begin();
      it != matching_set.end(); ++it)
  {
    Cell f(mm, *it);
    Point fp = f.midpoint();
    message("Facet midpoint : %8f, %8f, %8f", fp[0], fp[1], fp[2]);
    if(onEntity(&p[0], f))
    {
      message("On entity");
    }
    else
    {
      message("Not on entity");
    }
  }
}

//-----------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  Test t(argc, argv);
  bool test_id = false;
  bool test_2d = true;
  bool test_3d = true;

  if (test_id)
  {
    if (dolfin::MPI::size() == 1)
    {
      Mesh mesh("../../data/meshes/squareN128R.xml.gz");

      PeriodicSquareX subdomainX;
      MappedManifold mmX(mesh, subdomainX);
      Point pX(0.0, 0.5, 0.0);
      check_matching_facets(mmX, pX);

      PeriodicSquareY subdomainY;
      MappedManifold mmY(mesh, subdomainY);
      Point pY(0.5, 0.0, 0.0);
      check_matching_facets(mmY, pY);
    }

    if (dolfin::MPI::size() == 1)
    {
      Mesh mesh("../../data/meshes/cubeN32R.xml.gz");

      PeriodicCubeX subdomainX;
      MappedManifold mmX(mesh, subdomainX);
      Point pX(0.0, 0.5, 0.5);
      check_matching_facets(mmX, pX);

      PeriodicCubeY subdomainY;
      MappedManifold mmY(mesh, subdomainY);
      Point pY(0.5, 0.0, 0.5);
      check_matching_facets(mmY, pY);

      PeriodicCubeZ subdomainZ;
      MappedManifold mmZ(mesh, subdomainZ);
      Point pZ(0.5, 0.5, 0.0);
      check_matching_facets(mmZ, pZ);
    }
  }

  if (test_2d)
  {
    std::string name("LeftToRight2D");
    message("Square : " + name);
    Mesh mesh("../../data/meshes/squareN128R.xml.gz");
    LeftToRight2D subdomain;
    mesh.add_periodic_constraint(subdomain);
    ufl::FiniteElement cg1_2d(ufl::Family::CG, mesh.type(), 1);
    FiniteElementSpace spaceU(mesh, cg1_2d);
    PeriodicDofsMapping const& pdm = spaceU.dofmap().periodic_mapping();

    LeftToRight2D bc_subdomain;

    //
    uint Gdof = 0;
    Point G;
    real ** Hcoords = new real *[pdm.max_local_dimension()];
    for(uint c = 0; c < pdm.max_local_dimension(); ++c)
    {
      Hcoords[c] = new real[Space::MAX_DIMENSION];
    }
    uint count = 0;
    for(uint i = 0; i < pdm.num_Gdofs(); ++i)
    {
      pdm.tabulate_coordinates(i, &Gdof, &G[0], Hcoords, count);
      if(bc_subdomain.inside(&G[0], true))
      {

      }
    }

    //
    for(uint c = 0; c < pdm.max_local_dimension(); ++c)
    {
      delete[] Hcoords[c];
    }
    delete[] Hcoords;

    //
    write_Gdofs(name, spaceU);
    write_mesh(name, mesh);
  }

  return 0;

  if (test_2d)
  {
    std::string name("PeriodicSquare");
    message("Square : " + name);
    Mesh mesh("../../data/meshes/squareN128R.xml.gz");
    PeriodicSquareX subdomainX;
    mesh.add_periodic_constraint(subdomainX);
    PeriodicSquareY subdomainY;
    mesh.add_periodic_constraint(subdomainY);
    ufl::FiniteElement cg1_2d(ufl::Family::CG, mesh.type(), 1);
    FiniteElementSpace spaceU(mesh, cg1_2d);
    PeriodicDofsMapping const& pdm = spaceU.dofmap().periodic_mapping();
    //pdm.disp();

    //
    write_Gdofs(name, spaceU);
    write_mesh(name, mesh);
  }

  if (test_3d)
  {
    std::string name("LeftToRight3D");
    message("Cube : " + name);
    Mesh mesh("../../data/meshes/cubeN32R.xml.gz");
    LeftToRight3D subdomain;
    mesh.add_periodic_constraint(subdomain);
    ufl::FiniteElement cg1_3d(ufl::Family::CG, mesh.type(), 1);
    FiniteElementSpace spaceU(mesh, cg1_3d);

    //
    write_Gdofs(name, spaceU);
    write_mesh(name, mesh);
  }

  if (test_3d)
  {
    std::string name("PeriodicCubeX");
    message("Cube : " + name);
    Mesh mesh("../../data/meshes/cubeN32R.xml.gz");
    PeriodicCubeX subdomainX;
    mesh.add_periodic_constraint(subdomainX);
    ufl::FiniteElement cg1_3d(ufl::Family::CG, mesh.type(), 1);
    FiniteElementSpace spaceU(mesh, cg1_3d);

    //
    write_Gdofs(name, spaceU);
    write_mesh(name, mesh);
  }

  if (test_3d)
  {
    std::string name("PeriodicCubeY");
    message("Cube : " + name);
    Mesh mesh("../../data/meshes/cubeN32R.xml.gz");
    PeriodicCubeY subdomainY;
    mesh.add_periodic_constraint(subdomainY);
    ufl::FiniteElement cg1_3d(ufl::Family::CG, mesh.type(), 1);
    FiniteElementSpace spaceU(mesh, cg1_3d);

    //
    write_Gdofs(name, spaceU);
    write_mesh(name, mesh);
  }

  if (test_3d)
  {
    std::string name("PeriodicCubeZ");
    message("Cube : " + name);
    Mesh mesh("../../data/meshes/cubeN32R.xml.gz");
    PeriodicCubeZ subdomainZ;
    mesh.add_periodic_constraint(subdomainZ);
    ufl::FiniteElement cg1_3d(ufl::Family::CG, mesh.type(), 1);
    FiniteElementSpace spaceU(mesh, cg1_3d);

    //
    write_Gdofs(name, spaceU);
    write_mesh(name, mesh);
  }

  if (test_3d)
  {
    std::string name("PeriodicCubeXYZ");
    message("Cube : " + name);
    Mesh mesh("../../data/meshes/cubeN32R.xml.gz");
    PeriodicCubeX subdomainX;
    mesh.add_periodic_constraint(subdomainX);
    PeriodicCubeY subdomainY;
    mesh.add_periodic_constraint(subdomainY);
    PeriodicCubeZ subdomainZ;
    mesh.add_periodic_constraint(subdomainZ);
    ufl::FiniteElement cg1_3d(ufl::Family::CG, mesh.type(), 1);
    FiniteElementSpace spaceU(mesh, cg1_3d);

    //
    write_Gdofs(name, spaceU);
    write_mesh(name, mesh);
  }

  if (test_3d)
  {
    std::string name("PeriodicCubeYZX");
    message("Cube : " + name);
    Mesh mesh("../../data/meshes/cubeN32R.xml.gz");
    PeriodicCubeY subdomainY;
    mesh.add_periodic_constraint(subdomainY);
    PeriodicCubeZ subdomainZ;
    mesh.add_periodic_constraint(subdomainZ);
    PeriodicCubeX subdomainX;
    mesh.add_periodic_constraint(subdomainX);
    ufl::FiniteElement cg1_3d(ufl::Family::CG, mesh.type(), 1);
    FiniteElementSpace spaceU(mesh, cg1_3d);

    //
    write_Gdofs(name, spaceU);
    write_mesh(name, mesh);
  }

  if (test_3d)
  {
    std::string name("PeriodicCubeZXY");
    message("Cube : " + name);
    Mesh mesh("../../data/meshes/cubeN32R.xml.gz");
    PeriodicCubeZ subdomainZ;
    mesh.add_periodic_constraint(subdomainZ);
    PeriodicCubeX subdomainX;
    mesh.add_periodic_constraint(subdomainX);
    PeriodicCubeY subdomainY;
    mesh.add_periodic_constraint(subdomainY);
    ufl::FiniteElement cg1_3d(ufl::Family::CG, mesh.type(), 1);
    FiniteElementSpace spaceU(mesh, cg1_3d);

    //
    write_Gdofs(name, spaceU);
    write_mesh(name, mesh);
  }

  message("Finished tests");

  return 0;
}
