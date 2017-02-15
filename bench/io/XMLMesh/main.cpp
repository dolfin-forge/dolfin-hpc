
#include <dolfin/io/File.h>
#include <dolfin/main/init.h>
#include <dolfin/main/MPI.h>
#include <dolfin/mesh/BoundaryMesh.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshFunction.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
int main(int argc, char * argv[])
{
  dolfin_init(argc, argv);
  //---
  logm.verbose(1);
  logm.file();
//std::string filename("square");
  std::string filename("aneurysm");
  std::stringstream ss;
  ss << filename << dolfin::MPI::size();

  ///
  Mesh mesh(filename+".xml");
  MeshFunction<uint> p(mesh, mesh.topology().dim());
  p = dolfin::MPI::rank();
  File f(ss.str()+"_partitions.pvd");
  f << p;

//  mesh.init();

  BoundaryMesh& eb = mesh.exterior_boundary();
  File feb(filename+"_extboundary.pvd");
  feb << eb;
  BoundaryMesh& ib = mesh.interior_boundary();
  File fib(filename+"_intboundary.pvd");
  fib << ib;

  //---
  dolfin_finalize();
  return 0;
}
//-----------------------------------------------------------------------------
