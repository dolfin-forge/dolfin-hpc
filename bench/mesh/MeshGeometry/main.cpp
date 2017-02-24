#include <dolfin.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
int main(int argc, char** argv)
{
  Test T(argc, argv);
  uint const N = std::pow(2, 24);
  //---
  T.begin("Constructor with zero size");
  {
    /*
     * Zero initialization should not fail
     */

    MeshGeometry mg;
    for (uint dim = 1; dim <= Space::MAX_DIMENSION; ++dim)
    {
      mg.init(dim, 0);
      mg.clear();
    }
  }
  T.end();
  //---
  T.begin("Constructor and re-initializations");
  {
    /*
     * Initialize with different geometric dimensions and sizes
     */

    MeshGeometry mg;
    Point x;
    for (uint dim = 1; dim <= Space::MAX_DIMENSION; ++dim)
    {
      uint const size = std::rand() % (N / dim);
      message("dim = %-10u; size = %-10u", dim, size);
      mg.init(dim, size);
      //
      for (uint n = 0; n < mg.size(); ++n)
      {
        for (uint i = 0; i < dim; ++i)
        {
          x[i] = std::rand();
        }

        /*
         * Test that saved values match the original
         */
        for (uint i = 0; i < dim; ++i)
        {
          mg.set(n, x);
        }
        if (mg.point(n).distance(x) > DOLFIN_EPS)
        {
          error("Point values do not match");
        }

        /*
         * Test that saved values match the original
         */
        for (uint i = 0; i < dim; ++i)
        {
          mg.set(n, i, x[i]);
        }
        if (mg.point(n).distance(x) > DOLFIN_EPS)
        {
          error("Point values do not match");
        }
      }
      //
      mg.clear();
    }
  }
  T.end();
  //---
  return 0;
}
//-----------------------------------------------------------------------------

