#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/mesh/algorithm.h>
#include <dolfin/fem/UFCCell.h>

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
  size_t operator()(Cell& cell_)
  {
    UFCCell cell( cell_ );
    size_t min = cell.entity_indices[0][0];
    size_t max = cell.entity_indices[0][0];
    for (size_t i = 1; i < cell.num_vertices; ++i)
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
      size_t value = mesh.global_size(0);
      foreach<CellIterator, std::min >(mesh, DistGlobalIndex(), value);
      foreach<CellIterator, std::min >(mesh, DistGlobalIndex(), value);
    }
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
