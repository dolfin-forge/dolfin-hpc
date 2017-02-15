#include <dolfin.h>

#include "Poisson2D.h"
#include "PoissonP1.h"

using namespace dolfin;

int main(int argc, char** argv)
{
  dolfin_init(argc, argv);

  //
  {
    Mesh mesh("../../data/meshes/squareN100R.xml.gz");
    std::stringstream sg;
    sg << "_G" << dolfin::MPI::groupNumber();
    std::stringstream ss;
    ss << "mesh";
    ss << sg.str();
    message("Group %2d : Rank %2d\n", dolfin::MPI::groupNumber(),
            dolfin::MPI::rank());

    for (uint i = 1; i < dolfin::MPI::groupNumber(); ++i)
    {
      mesh.refine();
    }
    message("Mesh in MPI Group %d has %d cells", dolfin::MPI::groupNumber(),
            mesh.num_global_cells());
    File fmesh(ss.str() + ".pvd");
    fmesh << mesh;

    //
    // Create functions
    Source f(mesh);
    Flux g(mesh);

    // Create boundary condition
    Function u0(mesh, 0.0);
    DirichletBoundary boundary;
    DirichletBC bc(u0, mesh, boundary);

    // Define PDE
    PoissonBilinearForm a(mesh);
    PoissonLinearForm L(f, g);
    dolfin_set("PDE linear solver", "iterative");
    LinearPDE pde(a, L, mesh, bc);

    // Solve PDE
    Function u(mesh);

    pde.solve(u);

    // Save solution to file
    File file("poisson"+sg.str()+".pvd");
    file << u;

    //
  }

  dolfin_finalize();
  return 0;
}

