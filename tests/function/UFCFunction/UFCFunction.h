#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/function/UFCFunction.h>
#include <dolfin/function/Operators.h>
#include <dolfin/fem/UFCCellIterator.h>
#include <dolfin/mesh/Vertex.h>

using namespace dolfin;

//-----------------------------------------------------------------------------

struct MidpointDistance : public ValueSpace<>
{
  void evaluate(real* values, const real* x, const UFCCell& cell) const
  {
    values[0] = (*cell).midpoint().dist(Point(cell.geometric_dimension, x));
  }
};

//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_UFCFunction )
  {
    {
      Array<CellType *> cells = CellType::create_all();
      for (Array<CellType *>::const_iterator it = cells.begin();
          it != cells.end(); ++it)
      {
        CellType * cell = CellType::create(**it);
        Mesh refcell( *cell, EuclideanSpace( cell->dim() ) );
        (*it)->create_reference_cell(refcell);
        UFCFunction<MidpointDistance> f(refcell);
        real value;
        for (UFCCellIterator c(refcell); !c.end(); ++c)
        {
          for (VertexIterator v(c.cell()); !v.end(); ++v)
          {
            f.evaluate(&value, v->x(), *c);
          }
        }
        delete cell;
      }
      free( cells );
    }
    //---
    {
      Array<CellType *> cells = CellType::create_all();
      for (Array<CellType *>::const_iterator it = cells.begin();
          it != cells.end(); ++it)
      {
        CellType * cell = CellType::create(**it);
        Mesh refcell( *cell, EuclideanSpace( cell->dim() ) );
        (*it)->create_reference_cell(refcell);
        UFCFunction<Diameter<Cell> > f(refcell);
        real value;
        for (UFCCellIterator c(refcell); !c.end(); ++c)
        {
          for (VertexIterator v(c.cell()); !v.end(); ++v)
          {
            f.evaluate(&value, v->x(), *c);
          }
        }
        delete cell;
      }
      free( cells );
    }
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
