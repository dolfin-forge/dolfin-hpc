#include <dolfin.h>

#include "IthCoordinate.h"
#include "Coordinates.h"

#include <dolfin/function/FunctionDecomposition.h>

using namespace dolfin;

int main(int argc, char** argv)
{
  int ret = 0;

  //---------------------------------------------------------------------------
  Test t(argc, argv);

  Mesh sqR("square100R.xml.gz");

  message("Function decomposition");
  {
    typedef std::pair<ufl::Family::Type, uint> SpaceItem;
    Array<SpaceItem> spaces;
    spaces.push_back(SpaceItem(ufl::Family::DG, 0));
    spaces.push_back(SpaceItem(ufl::Family::CG, 1));
    spaces.push_back(SpaceItem(ufl::Family::CG, 2));

    Mesh& mesh = sqR;
    uint const gdim = mesh.geometry().dim();
    Coordinates X(mesh);
    Array<IthCoordinate *> Xi;
    for (uint d = 0; d < gdim; ++d)
    {
      Xi.push_back(new IthCoordinate(mesh, d));
    }

    //
    for (Array<SpaceItem>::const_iterator it = spaces.begin();
        it != spaces.end(); ++it)
    {
      ufl::FiniteElement Shi(it->first, mesh.type(), it->second);
      FiniteElementSpace Vh(mesh, Shi);
      Function Fi(Vh);

      ufl::VectorElement Shd(it->first, mesh.type(), it->second, gdim);
      FiniteElementSpace Wh(mesh, Shd);
      Function F(Wh);

      F.interpolate(X);
      Array<Function *> Si = FunctionDecomposition::compute(F);
      for (uint i = 0; i < Si.size(); ++i)
      {
        Fi.interpolate(*Xi[i]);
        Fi.vector() -= Si[i]->vector();
        Fi.sync_ghosts();
        real l2err = Fi.vector().norm(dolfin::l2);
        message("F[%d] - S[%d] l2 error = %16f", i, i, l2err);
      }

      // Cleanup
      while (!Si.empty())
      {
        delete Si.back();
        Si.pop_back();
      }
      Si.clear();
    }

    // Cleanup
    while (!Xi.empty())
    {
      delete Xi.back();
      Xi.pop_back();
    }
    Xi.clear();
  }

  message("That's all folks !");

  return ret;
}
