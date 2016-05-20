#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_CHECK

#include <dolfin/common/Test.h>
#include <dolfin/fem/UFCCellIterator.h>
#include <dolfin/mesh/algorithm.h>

#include <check.h>

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
START_TEST( test_algorithm )
  {
    int init_failed = 0;
    begin("test_algorithm");
    //---
    Test T;
    T.begin("Foreach min/max volume");
    {
      UnitInterval mesh(42);
      real value = 1.0;
      foreach<CellIterator>(mesh, Volume(), std::min<real>, value);
      message("%g", value);
      foreach<CellIterator>(mesh, Volume(), std::max<real>, value);
      message("%g", value);
    }
    T.end();
    //---
    T.begin("Foreach min/max distance global vertex index");
    {
      UnitInterval mesh(42);
      uint value = mesh.global_size(0);
      foreach<UFCCellIterator>(mesh, DistGlobalIndex(), std::min<uint>, value);
      message("%u", value);
      foreach<UFCCellIterator>(mesh, DistGlobalIndex(), std::max<uint>, value);
      message("%u", value);
    }
    T.end();
    //---
    end();
    skip();
    fail_unless( init_failed == 0 );
  }END_TEST
//-----------------------------------------------------------------------------

#endif
