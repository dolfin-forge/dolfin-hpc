#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/mesh/MeshData.h>
#include <dolfin/mesh/UnitSquare.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_MeshData_add )
  {
    UnitSquare m(4, 4);
    MeshData d(m);
    MeshValues<bool, Cell>   bC(m); d.add(bC); ck_assert(d.count(bC));
    MeshValues<bool, Vertex> bV(m); d.add(bV); ck_assert(d.count(bV));
    MeshValues<uint, Cell>   uC(m); d.add(uC); ck_assert(d.count(uC));
    MeshValues<uint, Vertex> uV(m); d.add(uV); ck_assert(d.count(uV));
    MeshValues<real, Cell>   rC(m); d.add(rC); ck_assert(d.count(rC));
    MeshValues<real, Vertex> rV(m); d.add(rV); ck_assert(d.count(rV));
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_MeshData_insert )
  {
    UnitSquare m(4, 4);
    MeshData d(m);
    MeshValues<bool, Cell>   bC(m);
    d.insert(bC); ck_assert(d.count(bC) == 1);
    d.insert(bC); ck_assert(d.count(bC) == 1);
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_MeshData_iterator )
  {
    UnitSquare m(4, 4);
    {
      MeshData d(m);
      ck_assert(d.size() == 0);
      uint ii = 0;
      for (MeshData::iterator<bool, Facet> it(d); it.valid(); ++it, ++ii)
      {
        begin("%u", it.pos()); it->disp(); endblock();
      }
      ck_assert(ii == 0);
    }
    //---
    {
      MeshData d(m);
      MeshValues<uint, Cell>   f0(m); d.add(f0); f0 = 0;
      MeshValues<uint, Cell>   f1(m); d.add(f1); f1 = 1;
      ck_assert(d.size() == 2);
      uint ii = d.size<uint, Cell>(); ck_assert(ii == 2);
      for (MeshData::iterator<uint, Cell> it(d); it.valid(); ++it, --ii)
      {
        ck_assert(it->dim()  == m.topology().dim());
        ck_assert(it->size() == f0.size());
        begin("%u", it.pos()); it->disp(); endblock();
        for (Cell::iterator c(m); !c.end(); ++c)
        {
          ck_assert((*it)(*c) == it.pos());
        }
      }
      ck_assert(ii == 0);
    }
    //---
    {
      MeshData d(m);
      MeshValues<real, Vertex>   f0(m); d.add(f0); f0 = 0.;
      MeshValues<real, Vertex>   f1(m); d.add(f1); f1 = 1.;
      MeshValues<real, Vertex>   f2(m); d.add(f2); f2 = 2.;
      ck_assert(d.size() == 3);
      uint ii = d.size<real, Vertex>(); ck_assert(ii == 3);
      for (MeshData::iterator<real, Vertex> it(d); it.valid(); ++it, --ii)
      {
        ck_assert(it->dim() == 0);
        ck_assert(it->size() == f0.size());
        begin("%u", it.pos()); it->disp(); endblock();
        for (Vertex::iterator v(m); !v.end(); ++v)
        {
          ck_assert((*it)(*v) == it.pos());
        }
      }
      ck_assert(ii == 0);
    }
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
