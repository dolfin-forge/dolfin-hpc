#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/fem/UFCCellIterator.h>
#include <dolfin/mesh/algorithm.h>

using namespace dolfin;

/*
 * These template functions provide a convenient way for defining operators on
 * mesh entities or even a ufc::cell using functors.
 */

//-----------------------------------------------------------------------------

struct Volume
{
  real operator()(Cell& cell)
  {
    return cell.volume();
  }
};

//-----------------------------------------------------------------------------

struct DistGlobalIndex
{
  uint operator()(UFCCell& cell)
  {
    uint min = cell.entity_indices[0][0];
    uint max = cell.entity_indices[0][0];
    for (uint i = 1; i < cell.num_vertices; ++i)
    {
      min = std::min(min, cell.entity_indices[0][i]);
      max = std::max(max, cell.entity_indices[0][i]);
    }
    return (max - min);
  }
};

//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_algorithm )
  {
    {
      UnitInterval mesh(42);
      real value = 1.0;
      foreach<CellIterator, std::min<real> >(mesh, Volume(), value);
      foreach<CellIterator, std::max<real> >(mesh, Volume(), value);
    }
    //---
    {
      UnitInterval mesh(42);
      uint value = mesh.global_size(0);
      foreach<UFCCellIterator, std::min<uint> >(mesh, DistGlobalIndex(), value);
      foreach<UFCCellIterator, std::min<uint> >(mesh, DistGlobalIndex(), value);
    }
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
