#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_CHECK

#include <dolfin/common/Test.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/CellType.h>
#include <dolfin/mesh/PointCell.h>
#include <dolfin/mesh/IntervalCell.h>
#include <dolfin/mesh/TriangleCell.h>
#include <dolfin/mesh/TetrahedronCell.h>
#include <dolfin/mesh/QuadrilateralCell.h>
#include <dolfin/mesh/HexahedronCell.h>

#include <sstream>

#include <check.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
void check_reference_cell(CellType& cell, Mesh& refcell)
{
  Cell c(refcell, 0);
  begin("Properties");
  message("Volume       : %+e", cell.volume(c));
  message("Diameter     : %+e", cell.diameter(c));
  message("Circumradius : %+e", cell.circumradius(c));
  end();
  // Initialize all connectivities
  for (dolfin::uint i = 0; i <= refcell.topology().dim(); ++i)
  {
    for (dolfin::uint j = i; j <= refcell.topology().dim(); ++j)
    {
      refcell.init(i, j);
    }
  }
  //
  for (CellIterator c(refcell); !c.end(); ++c)
  {
    cell.check(*c);
  }
  dolfin::uint nc0 = refcell.num_cells();

  dolfin::uint const N = std::pow(2.0, (int) (5 - cell.dim()));
  for (dolfin::uint l = 1; l <= N; ++l)
  {
    refcell.refine();
    ck_assert_int_eq(refcell.num_cells(), nc0 * cell.num_refined_cells());
    nc0 = refcell.num_cells();
    MeshFunction<dolfin::uint> vi(refcell, 0);
    for (VertexIterator v(refcell); !v.end(); ++v)
    {
      vi.set(*v, v->index());
    }
  }
}
//-----------------------------------------------------------------------------
START_TEST( test_PointCell )
  {
    int init_failed = 0;
    begin("test_PointCell");
    //---
    PointCell cell;
    ck_assert_int_eq(cell.dim(), 0);
    ck_assert_int_eq(cell.num_entities(0), 1);
    ck_assert_int_eq(cell.num_vertices(0), 1);
    cell.disp();
    //
    CellType * ct0 = CellType::create(CellType::point);
    ck_assert_int_eq(cell.cellType(), ct0->cellType());
    delete ct0;
    //
    Mesh refcell = cell.create_reference_cell();
    // Most member functions are undefined
    //check_reference_cell(cell, refcell);
    //---
    end();
    skip();
    fail_unless( init_failed == 0 );
  }END_TEST
//-----------------------------------------------------------------------------
START_TEST( test_IntervalCell )
  {
    int init_failed = 0;
    begin("test_IntervalCell");
    //---
    IntervalCell cell;
    ck_assert_int_eq(cell.dim(), 1);
    ck_assert_int_eq(cell.num_entities(0), 2);
    ck_assert_int_eq(cell.num_entities(1), 1);
    ck_assert_int_eq(cell.num_vertices(0), 1);
    ck_assert_int_eq(cell.num_vertices(1), 2);
    cell.disp();
    //
    CellType * ct0 = CellType::create(CellType::interval);
    ck_assert_int_eq(cell.cellType(), ct0->cellType());
    delete ct0;
    // UFC convention
    Mesh refcell = cell.create_reference_cell();
    check_reference_cell(cell, refcell);
    //---
    end();
    fail_unless( init_failed == 0 );
  }END_TEST
//-----------------------------------------------------------------------------
START_TEST( test_TriangleCell )
  {
    int init_failed = 0;
    begin("test_TriangleCell");
    //---
    TriangleCell cell;
    ck_assert_int_eq(cell.dim(), 2);
    ck_assert_int_eq(cell.num_entities(0), 3);
    ck_assert_int_eq(cell.num_entities(1), 3);
    ck_assert_int_eq(cell.num_entities(2), 1);
    ck_assert_int_eq(cell.num_vertices(0), 1);
    ck_assert_int_eq(cell.num_vertices(1), 2);
    ck_assert_int_eq(cell.num_vertices(2), 3);
    cell.disp();
    //
    CellType * ct0 = CellType::create(CellType::triangle);
    ck_assert_int_eq(cell.cellType(), ct0->cellType());
    delete ct0;
    // UFC convention
    Mesh refcell = cell.create_reference_cell();
    check_reference_cell(cell, refcell);
    //---
    end();
    fail_unless( init_failed == 0 );
  }END_TEST
//-----------------------------------------------------------------------------
START_TEST( test_TetrahedronCell )
  {
    int init_failed = 0;
    begin("test_TetrahedronCell");
    //---
    TetrahedronCell cell;
    ck_assert_int_eq(cell.dim(), 3);
    ck_assert_int_eq(cell.num_entities(0), 4);
    ck_assert_int_eq(cell.num_entities(1), 6);
    ck_assert_int_eq(cell.num_entities(2), 4);
    ck_assert_int_eq(cell.num_entities(3), 1);
    ck_assert_int_eq(cell.num_vertices(0), 1);
    ck_assert_int_eq(cell.num_vertices(1), 2);
    ck_assert_int_eq(cell.num_vertices(2), 3);
    ck_assert_int_eq(cell.num_vertices(3), 4);
    cell.disp();
    //
    CellType * ct0 = CellType::create(CellType::tetrahedron);
    ck_assert_int_eq(cell.cellType(), ct0->cellType());
    delete ct0;
    // UFC convention
    Mesh refcell = cell.create_reference_cell();
    check_reference_cell(cell, refcell);
    //---
    end();
    fail_unless( init_failed == 0 );
  }END_TEST
//-----------------------------------------------------------------------------
START_TEST( test_QuadrilateralCell )
  {
    int init_failed = 0;
    begin("test_QuadrilateralCell");
    //---
    QuadrilateralCell cell;
    ck_assert_int_eq(cell.dim(), 2);
    ck_assert_int_eq(cell.num_entities(0), 4);
    ck_assert_int_eq(cell.num_entities(1), 4);
    ck_assert_int_eq(cell.num_entities(2), 1);
    ck_assert_int_eq(cell.num_vertices(0), 1);
    ck_assert_int_eq(cell.num_vertices(1), 2);
    ck_assert_int_eq(cell.num_vertices(2), 4);
    cell.disp();
    //
    CellType * ct0 = CellType::create(CellType::quadrilateral);
    ck_assert_int_eq(cell.cellType(), ct0->cellType());
    delete ct0;
    // UFC convention
    Mesh refcell = cell.create_reference_cell();
    Cell c(refcell, 0);
    refcell.init(1, 0);
    refcell.init(cell.dim(), 1);
    cell.check(c);
    //---
    end();
    fail_unless( init_failed == 0 );
  }END_TEST
//-----------------------------------------------------------------------------
START_TEST( test_HexahedronCell )
  {
    int init_failed = 0;
    begin("test_HexahedronCell");
    //---
    HexahedronCell cell;
    ck_assert_int_eq(cell.dim(), 3);
    ck_assert_int_eq(cell.num_entities(0), 8);
    ck_assert_int_eq(cell.num_entities(1), 12);
    ck_assert_int_eq(cell.num_entities(2), 6);
    ck_assert_int_eq(cell.num_entities(3), 1);
    ck_assert_int_eq(cell.num_vertices(0), 1);
    ck_assert_int_eq(cell.num_vertices(1), 2);
    ck_assert_int_eq(cell.num_vertices(2), 4);
    ck_assert_int_eq(cell.num_vertices(3), 8);
    cell.disp();
    // Check cell type enum
    CellType * ct0 = CellType::create(CellType::hexahedron);
    delete ct0;
    // UFC convention
    Mesh refcell = cell.create_reference_cell();
    Cell c(refcell, 0);
    refcell.init(1, 0);
    refcell.init(2, 0);
    refcell.init(cell.dim(), 1);
    cell.check(c);
    //---
    end();
    fail_unless( init_failed == 0 );
  }END_TEST
//-----------------------------------------------------------------------------

#endif
