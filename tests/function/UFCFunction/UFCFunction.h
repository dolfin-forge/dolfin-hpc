#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_CHECK

#include <dolfin/common/Test.h>
#include <dolfin/function/UFCFunction.h>
#include <dolfin/function/impl/CellDiameter.h>
#include <dolfin/fem/UFCCellIterator.h>
#include <dolfin/mesh/Vertex.h>

#include <check.h>

using namespace dolfin;

//-----------------------------------------------------------------------------

struct MidpointDistance : public ValueSpace<>
{
  ///
  void evaluate(real* values, const real* x, const UFCCell& cell) const
  {
    values[0] = cell.cell->midpoint().dist(Point(cell.geometric_dimension, x));
  }

};

//-----------------------------------------------------------------------------
START_TEST( test_UFCFunction )
{
  int init_failed = 0;
  Test T;
  //---
  T.begin("test_UFCFunction : MidpointDistance");
  {
    Array<CellType *> cells = CellType::create_all();
    for (Array<CellType *>::const_iterator it = cells.begin();
          it != cells.end(); ++it)
    {
      Mesh refcell;
      (*it)->create_reference_cell(refcell);
      begin("Reference cell: %s", (*it)->str().c_str());
      UFCFunction<MidpointDistance> f(refcell);
      real value;
      for (UFCCellIterator c(refcell); !c.end(); ++c)
      {
        for (VertexIterator v(c.cell()); !v.end(); ++v)
        {
          f.evaluate(&value, v->x(), *c);
          message("Distance from vertex %d to midpoint : %g", v.pos(), value);
        }
      }
      end();
    }
    cells.free();
  }
  T.end();
  //---
  T.begin("test_UFCFunction : CellDiameter");
  {
    Array<CellType *> cells = CellType::create_all();
    for (Array<CellType *>::const_iterator it = cells.begin();
        it != cells.end(); ++it)
    {
      Mesh refcell;
      (*it)->create_reference_cell(refcell);
      begin("Reference cell: %s", (*it)->str().c_str());
      UFCFunction<CellDiameter> f(refcell);
      real value;
      for (UFCCellIterator c(refcell); !c.end(); ++c)
      {
        for (VertexIterator v(c.cell()); !v.end(); ++v)
        {
          f.evaluate(&value, v->x(), *c);
          message("Cell diameter : %g", v.pos(), value);
        }
      }
      end();
    }
    cells.free();
  }
  T.end();
  //---
  ck_assert( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------

#endif
