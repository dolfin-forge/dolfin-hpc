#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/fem/DofMap.h>
#include <dolfin/fem/DofNumbering.h>
#include <dolfin/mesh/CellType.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/ufl/UFLFiniteElementSpace.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_DofNumbering )
  {
    CellType * cell = NULL;
    ufc::dofmap * dofmap = NULL;
    ufl::ElementList const& L = ElementLibrary::elements();
    for (ufl::FiniteElementSpace const * e = L.first(); L.valid(); e = L.next())
    {
      cell = CellType::create(e->cell());
      dofmap = ElementLibrary::create_dof_map(DofMap::make_signature(e->repr()));
      //---
      Mesh refcell;
      cell->create_reference_cell(refcell);
      DofNumbering * numbering = DofNumbering::create(refcell, *dofmap);
      delete numbering;
      //---
      delete dofmap;
      delete cell;
    }
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
