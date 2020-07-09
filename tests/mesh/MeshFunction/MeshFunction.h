#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/mesh/EuclideanSpace.h>
#include <dolfin/mesh/MeshValues.h>
#include <dolfin/mesh/CellTypes.h>
#include <dolfin/mesh/CellIterator.h>
#include <dolfin/mesh/EdgeIterator.h>
#include <dolfin/mesh/FaceIterator.h>
#include <dolfin/mesh/FacetIterator.h>
#include <dolfin/mesh/VertexIterator.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
template<class CellType, typename T, class Entity>
void check_reference_cell()
{
  CellType cellt;
  Mesh     cellm( cellt, EuclideanSpace( cellt.dim() ) );
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

  // Copy constructor
  MeshValues<T, Entity> N(M);
  ck_assert(M == N);
}
//-----------------------------------------------------------------------------
template<typename T, typename V>
void check_conversion_vertices_x0()
{
  QuadrilateralCell cellt;
  Mesh              cellm( cellt, EuclideanSpace( cellt.dim() ) );
  cellt.create_reference_cell(cellm);
  cellm.geometry() -= Point(0.25, 0.25);
  cellm.geometry() *= 2.0;
  MeshValues<T, Vertex> M(cellm);
  for (Vertex::iterator v(M.mesh()); !v.end(); ++v)
  {
    M(*v) = static_cast<T>(v->x()[0]);
  }

  MeshValues<V, Vertex> N(M);
  ck_assert(M.dim() == N.dim());
  ck_assert(M.size() == N.size());

  MeshValues<V, Vertex> O(cellm);
  O = M;
  ck_assert(N == O);

  MeshValues<T, Vertex> P(O);
  ck_assert(M.dim() == P.dim());
  ck_assert(M.size() == P.size());
}
//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_MeshFunction )
  {
    // Check mesh functions on reference cells
    {
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
    }
    // Check conversions
    {
      // bool -> int
      check_conversion_vertices_x0<bool, int>();
      // bool -> uint
      check_conversion_vertices_x0<bool, uint>();
      // bool -> float
      check_conversion_vertices_x0<bool, float>();
      // bool -> real
      check_conversion_vertices_x0<bool, real>();

      // int -> bool
      check_conversion_vertices_x0<int, bool>();
      // int -> uint
      check_conversion_vertices_x0<int, uint>();
      // int -> float
      check_conversion_vertices_x0<int, float>();
      // int -> real
      check_conversion_vertices_x0<int, real>();

      // uint -> bool
      check_conversion_vertices_x0<uint, bool>();
      // uint -> int
      check_conversion_vertices_x0<uint, int>();
      // uint -> float
      check_conversion_vertices_x0<uint, float>();
      // uint -> real
      check_conversion_vertices_x0<uint, real>();

      // float -> bool
      check_conversion_vertices_x0<float, bool>();
      // float -> int
      check_conversion_vertices_x0<float, int>();
      // float -> uint
      check_conversion_vertices_x0<float, uint>();
      // float -> real
      check_conversion_vertices_x0<float, real>();

      // real -> bool
      check_conversion_vertices_x0<real, bool>();
      // real -> int
      check_conversion_vertices_x0<real, int>();
      // real -> uint
      check_conversion_vertices_x0<real, uint>();
      // real -> float
      check_conversion_vertices_x0<real, float>();
    }
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
