#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/mesh/CellTypes.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshEditor.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
template<class T>
void check_reference_cell_refinement()
{
  T cell;
  Mesh rc0;
  cell.create_reference_cell(rc0);
  Mesh rc1;
  MeshEditor me(rc1, cell);
  me.init_cells(cell.RefinementPattern::num_refined_cells(rc0));
  me.init_vertices(cell.RefinementPattern::num_refined_vertices(rc0));
  Cell mc0(rc0, 0);
  uint nci = 0;
  cell.refine_cell(mc0, me, nci);
  me.close();
}

//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_PointCell )
  {
    PointCell cell;
    // Topology
    ck_assert_int_eq(cell.dim(), 0);
    ck_assert_int_eq(cell.num_entities(0), 1);
    ck_assert_int_eq(cell.num_vertices(0), 1);
    ck_assert_int_eq(cell.num_entities(0, 0), 1);
    // Refinement pattern
    ck_assert_int_eq(cell.num_refined_cells(), 1);
    ck_assert_int_eq(cell.num_refined_vertices(0), 1);
    // Cell refinement
    check_reference_cell_refinement<PointCell>();
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_IntervalCell )
  {
    IntervalCell cell;
    // Topology
    ck_assert_int_eq(cell.dim(), 1);
    ck_assert_int_eq(cell.num_entities(0), 2);
    ck_assert_int_eq(cell.num_entities(1), 1);
    ck_assert_int_eq(cell.num_vertices(0), 1);
    ck_assert_int_eq(cell.num_vertices(1), 2);
    ck_assert_int_eq(cell.num_entities(0, 0), 1);
    ck_assert_int_eq(cell.num_entities(0, 1), 0);
    ck_assert_int_eq(cell.num_entities(1, 0), 2);
    ck_assert_int_eq(cell.num_entities(1, 1), 1);
    // Refinement pattern
    ck_assert_int_eq(cell.num_refined_cells(), 2);
    ck_assert_int_eq(cell.num_refined_vertices(0), 1);
    ck_assert_int_eq(cell.num_refined_vertices(1), 1);
    // Cell refinement
    check_reference_cell_refinement<IntervalCell>();
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_TriangleCell )
  {
    TriangleCell cell;
    // Topology
    ck_assert_int_eq(cell.dim(), 2);
    ck_assert_int_eq(cell.num_entities(0), 3);
    ck_assert_int_eq(cell.num_entities(1), 3);
    ck_assert_int_eq(cell.num_entities(2), 1);
    ck_assert_int_eq(cell.num_vertices(0), 1);
    ck_assert_int_eq(cell.num_vertices(1), 2);
    ck_assert_int_eq(cell.num_vertices(2), 3);
    ck_assert_int_eq(cell.num_entities(0, 0), 1);
    ck_assert_int_eq(cell.num_entities(0, 1), 0);
    ck_assert_int_eq(cell.num_entities(0, 2), 0);
    ck_assert_int_eq(cell.num_entities(1, 0), 2);
    ck_assert_int_eq(cell.num_entities(1, 1), 1);
    ck_assert_int_eq(cell.num_entities(1, 2), 0);
    ck_assert_int_eq(cell.num_entities(2, 0), 3);
    ck_assert_int_eq(cell.num_entities(2, 1), 3);
    ck_assert_int_eq(cell.num_entities(2, 2), 1);
    // Refinement pattern
    ck_assert_int_eq(cell.num_refined_cells(), 4);
    ck_assert_int_eq(cell.num_refined_vertices(0), 1);
    ck_assert_int_eq(cell.num_refined_vertices(1), 1);
    ck_assert_int_eq(cell.num_refined_vertices(2), 0);
    // Cell refinement
    check_reference_cell_refinement<TriangleCell>();
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_TetrahedronCell )
  {
    TetrahedronCell cell;
    // Topology
    ck_assert_int_eq(cell.dim(), 3);
    ck_assert_int_eq(cell.num_entities(0), 4);
    ck_assert_int_eq(cell.num_entities(1), 6);
    ck_assert_int_eq(cell.num_entities(2), 4);
    ck_assert_int_eq(cell.num_entities(3), 1);
    ck_assert_int_eq(cell.num_vertices(0), 1);
    ck_assert_int_eq(cell.num_vertices(1), 2);
    ck_assert_int_eq(cell.num_vertices(2), 3);
    ck_assert_int_eq(cell.num_vertices(3), 4);
    ck_assert_int_eq(cell.num_entities(0, 0), 1);
    ck_assert_int_eq(cell.num_entities(0, 1), 0);
    ck_assert_int_eq(cell.num_entities(0, 2), 0);
    ck_assert_int_eq(cell.num_entities(0, 3), 0);
    ck_assert_int_eq(cell.num_entities(1, 0), 2);
    ck_assert_int_eq(cell.num_entities(1, 1), 1);
    ck_assert_int_eq(cell.num_entities(1, 2), 0);
    ck_assert_int_eq(cell.num_entities(1, 3), 0);
    ck_assert_int_eq(cell.num_entities(2, 0), 3);
    ck_assert_int_eq(cell.num_entities(2, 1), 3);
    ck_assert_int_eq(cell.num_entities(2, 2), 1);
    ck_assert_int_eq(cell.num_entities(2, 3), 0);
    ck_assert_int_eq(cell.num_entities(3, 0), 4);
    ck_assert_int_eq(cell.num_entities(3, 1), 6);
    ck_assert_int_eq(cell.num_entities(3, 2), 4);
    ck_assert_int_eq(cell.num_entities(3, 3), 1);
    // Refinement pattern
    ck_assert_int_eq(cell.num_refined_cells(), 8);
    ck_assert_int_eq(cell.num_refined_vertices(0), 1);
    ck_assert_int_eq(cell.num_refined_vertices(1), 1);
    ck_assert_int_eq(cell.num_refined_vertices(2), 0);
    ck_assert_int_eq(cell.num_refined_vertices(3), 0);
    // Cell refinement
    check_reference_cell_refinement<TetrahedronCell>();
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_QuadrilateralCell )
  {
    QuadrilateralCell cell;
    // Topology
    ck_assert_int_eq(cell.dim(), 2);
    ck_assert_int_eq(cell.num_entities(0), 4);
    ck_assert_int_eq(cell.num_entities(1), 4);
    ck_assert_int_eq(cell.num_entities(2), 1);
    ck_assert_int_eq(cell.num_vertices(0), 1);
    ck_assert_int_eq(cell.num_vertices(1), 2);
    ck_assert_int_eq(cell.num_vertices(2), 4);
    ck_assert_int_eq(cell.num_entities(0, 0), 1);
    ck_assert_int_eq(cell.num_entities(0, 1), 0);
    ck_assert_int_eq(cell.num_entities(0, 2), 0);
    ck_assert_int_eq(cell.num_entities(1, 0), 2);
    ck_assert_int_eq(cell.num_entities(1, 1), 1);
    ck_assert_int_eq(cell.num_entities(1, 2), 0);
    ck_assert_int_eq(cell.num_entities(2, 0), 4);
    ck_assert_int_eq(cell.num_entities(2, 1), 4);
    ck_assert_int_eq(cell.num_entities(2, 2), 1);
    // Refinement pattern
    ck_assert_int_eq(cell.num_refined_cells(), 4);
    ck_assert_int_eq(cell.num_refined_vertices(0), 1);
    ck_assert_int_eq(cell.num_refined_vertices(1), 1);
    ck_assert_int_eq(cell.num_refined_vertices(2), 1);
    // Cell refinement
    check_reference_cell_refinement<QuadrilateralCell>();
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_HexahedronCell )
  {
    HexahedronCell cell;
    // Topology
    ck_assert_int_eq(cell.dim(), 3);
    ck_assert_int_eq(cell.num_entities(0), 8);
    ck_assert_int_eq(cell.num_entities(1), 12);
    ck_assert_int_eq(cell.num_entities(2), 6);
    ck_assert_int_eq(cell.num_entities(3), 1);
    ck_assert_int_eq(cell.num_vertices(0), 1);
    ck_assert_int_eq(cell.num_vertices(1), 2);
    ck_assert_int_eq(cell.num_vertices(2), 4);
    ck_assert_int_eq(cell.num_vertices(3), 8);
    ck_assert_int_eq(cell.num_entities(0, 0), 1);
    ck_assert_int_eq(cell.num_entities(0, 1), 0);
    ck_assert_int_eq(cell.num_entities(0, 2), 0);
    ck_assert_int_eq(cell.num_entities(0, 3), 0);
    ck_assert_int_eq(cell.num_entities(1, 0), 2);
    ck_assert_int_eq(cell.num_entities(1, 1), 1);
    ck_assert_int_eq(cell.num_entities(1, 2), 0);
    ck_assert_int_eq(cell.num_entities(1, 3), 0);
    ck_assert_int_eq(cell.num_entities(2, 0), 4);
    ck_assert_int_eq(cell.num_entities(2, 1), 4);
    ck_assert_int_eq(cell.num_entities(2, 2), 1);
    ck_assert_int_eq(cell.num_entities(2, 3), 0);
    ck_assert_int_eq(cell.num_entities(3, 0), 8);
    ck_assert_int_eq(cell.num_entities(3, 1), 12);
    ck_assert_int_eq(cell.num_entities(3, 2), 6);
    ck_assert_int_eq(cell.num_entities(3, 3), 1);
    // Refinement pattern
    ck_assert_int_eq(cell.num_refined_cells(), 8);
    ck_assert_int_eq(cell.num_refined_vertices(0), 1);
    ck_assert_int_eq(cell.num_refined_vertices(1), 1);
    ck_assert_int_eq(cell.num_refined_vertices(2), 1);
    ck_assert_int_eq(cell.num_refined_vertices(3), 1);
    // Cell refinement
    check_reference_cell_refinement<HexahedronCell>();
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
