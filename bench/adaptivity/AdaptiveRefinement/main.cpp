#include <dolfin.h>

using namespace dolfin;

//
class CheckerBoard : public Value<CheckerBoard>
{
public:

  CheckerBoard(uint dim, uint& level, uint cycle, bool cyclic_flip) :
      dim_(dim),
      level_(level),
      cycle_(cycle),
      cyclic_flip_(cyclic_flip)
  {
  }

  void eval(real * values, const real * x) const
  {
    values[0] = 1.0;
    real p = std::pow(2., level_ % cycle_) * DOLFIN_PI;
    for (uint d = 0; d < dim_; ++d)
    {
      values[0] *= std::sin(p * x[d]);
    }
    // Flip every cycle except level 0 of each cycle
    real ith = level_;
    if (cyclic_flip_)
    {
      ith = (
          (level_ % cycle_ == 0) ? 0 : std::floor(real(level_) / real(cycle_)));
    }
    values[0] = std::max(std::pow(-1.0, ith) * values[0], 0.0);
  }

private:
  uint dim_;
  uint& level_;
  uint const cycle_;
  bool cyclic_flip_;
};

int main(int argc, char *argv[])
{
  dolfin_init(argc, argv);
  Mesh mesh("../../data/meshes/squareN100R.xml.gz");
  //---------------------------------------------------------------------------
  bool const save_file = false;

  //TEST:
  message("Checkerboard refinement:");
  {
    uint const LMAX = 8;
    uint level = 0;
    CheckerBoard cb(mesh.topology().dim(), level, 4, true);
    Analytic<CheckerBoard> em(mesh, cb);
    ufl::FiniteElement cg1(ufl::Family::CG, mesh.type(), 1);

    real th = 0.5;

    //
    dolfin_add("output_format", "vtk");
    dolfin_add("adapt_algorithm", "rivara");
    dolfin_add("adapt_project", false);
    for (; level < LMAX; ++level)
    {
      Function up0(mesh, cg1);
      up0.interpolate(em);

//      if (save_file)
//      {
//        std::stringstream ss0;
//        ss0 << "em" << level << ".pvd";
//        File f0(ss0.str());
//        f0 << em;
//      }

      MeshValues<bool, Cell> mrkr(mesh);
      real v;
      for (CellIterator cell(mesh); !cell.end(); ++cell)
      {
        cb.eval(&v, &cell->midpoint()[0]);
        mrkr.set(*cell, (v >= th));
      }

      if (save_file)
      {
        std::stringstream ss1;
        ss1 << "mrkr" << level << ".pvd";
        File f1(ss1.str());
        f1 << mrkr;
      }

      if (dolfin_get("adapt_project"))
      {
        dolfin_set("Load balancer redistribute", false);
        Array<Function *> to_project;
        to_project.push_back(&up0);
        AdaptiveRefinement::refine_and_project(mesh, to_project, mrkr);
        dolfin_set("adapt_projected", true);
        to_project.clear();

        std::stringstream ss;
        ss << "cb_project" << level << ".bin";
        File up0_file(ss.str());
        up0_file << up0;
      }
      else
      {
        AdaptiveRefinement::refine(mesh, mrkr);
      }

      if (save_file)
      {
        MeshValues<uint, Cell> rank(mesh);
        rank = dolfin::MPI::rank();
        std::stringstream ss2;
        ss2 << "rank" << level << ".pvd";
        File f2(ss2.str());
        f2 << rank;
      }
    }
  }

  dolfin_finalize();
  return 0;
}

