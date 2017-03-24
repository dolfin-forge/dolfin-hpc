#include <dolfin/config/dolfin_config.h>
#include <dolfin/common/Test.h>
#include <dolfin/log/log.h>

#ifdef HAVE_CHECK

#include <dolfin/fem/DofMap.h>
#include <dolfin/fem/DofNumbering.h>
#include <dolfin/mesh/CellType.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/ufl/UFLFiniteElementSpace.h>

#include <check.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
START_TEST( test_DofNumbering )
{
  int init_failed = 0;
  Test T;
  //---
  T.begin("test_DofNumbering");
  {
    ufl::ElementList const& list = ElementLibrary::elements();
    for (ufl::FiniteElementSpace const * it = list.first(); list.valid();
        it = list.next())
    {
      begin("%s", it->str().c_str());
      //---
      CellType * cell = CellType::create(it->cell());
      ufc::dofmap * ufc_dofmap = ElementLibrary::create_dof_map(DofMap::make_signature(it->repr()));

      Mesh refcell;
      cell->create_reference_cell(refcell);
      DofNumbering * numbering = DofNumbering::create(refcell, *ufc_dofmap);
      numbering->disp();
      delete numbering;

      delete ufc_dofmap;
      delete cell;
      //---
      end();
    }
  }
  T.end();
  //---
  ck_assert( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------

#endif
