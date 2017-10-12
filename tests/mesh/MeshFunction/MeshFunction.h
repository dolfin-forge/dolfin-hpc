#include <dolfin/common/Check.h>

#ifdef HAVE_CHECK

#include <dolfin/mesh/MeshValues.h>
#include <dolfin/mesh/PointCell.h>
#include <dolfin/mesh/IntervalCell.h>
#include <dolfin/mesh/TriangleCell.h>
#include <dolfin/mesh/TetrahedronCell.h>
#include <dolfin/mesh/QuadrilateralCell.h>
#include <dolfin/mesh/HexahedronCell.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
template<class CellType, typename T, class Entity>
void check_reference_cell()
{
  CellType cellt;
  Mesh     cellm;
  cellt.create_reference_cell(cellm);
  MeshValues<T, Entity> M(cellm);

  ck_assert(M.size() == cellm.size(M.dim()));
  ck_assert(M.dim() ==  entity_dimension<Entity>(cellm));

  // Value accessor
  {
    for (typename Entity::iterator it(cellm); !it.end(); ++it)
    {
      M(it->index()) = it->index();
    }
    for (typename Entity::iterator it(cellm); !it.end(); ++it)
    {
      ck_assert(M(it->index()) == it->index());
    }
  }
  // Array accessor
  {
    for (typename Entity::iterator it(cellm); !it.end(); ++it)
    {
      M[it->index()][0] = it->index();
    }
    for (typename Entity::iterator it(cellm); !it.end(); ++it)
    {
      ck_assert(M(it->index()) == it->index());
    }
  }

  // Display
  M.disp();

  // Copy constructor
  MeshValues<T, Entity> N(M);
  ck_assert(M == N);

}
//-----------------------------------------------------------------------------
START_TEST( test_MeshFunction )
  {
    int init_failed = 0;
    begin("test_MeshFunction");
    //---
    check_reference_cell<PointCell, uint, Vertex>();
    check_reference_cell<PointCell, uint, Cell>();
    //
    check_reference_cell<IntervalCell, uint, Vertex>();
    check_reference_cell<IntervalCell, uint, Cell>();
    //
    check_reference_cell<TriangleCell, uint, Vertex>();
    check_reference_cell<TriangleCell, uint, Edge>();
    check_reference_cell<TriangleCell, uint, Cell>();
    //
    check_reference_cell<TetrahedronCell, uint, Vertex>();
    check_reference_cell<TetrahedronCell, uint, Edge>();
    check_reference_cell<TetrahedronCell, uint, Face>();
    check_reference_cell<TetrahedronCell, uint, Cell>();
    //
    check_reference_cell<QuadrilateralCell, uint, Vertex>();
    check_reference_cell<QuadrilateralCell, uint, Edge>();
    check_reference_cell<QuadrilateralCell, uint, Cell>();
    //
    check_reference_cell<HexahedronCell, uint, Vertex>();
    check_reference_cell<HexahedronCell, uint, Edge>();
    check_reference_cell<HexahedronCell, uint, Face>();
    check_reference_cell<HexahedronCell, uint, Cell>();
    //---
    end();
    fail_unless( init_failed == 0 );
  }END_TEST
//-----------------------------------------------------------------------------

#endif
