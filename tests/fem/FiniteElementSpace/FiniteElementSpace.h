#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/elements/Elements.h>
#include <dolfin/elements/ElementLibrary.h>
#include <dolfin/fem/DofMap.h>
#include <dolfin/fem/FiniteElement.h>
#include <dolfin/fem/FiniteElementSpace.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/EuclideanSpace.h>

using namespace dolfin;

using ufl::Cell;
using ufl::Domain;
using ufl::Family;
using ufl::Object;
using ufl::Space;

//----------------------------------------------------------------------------
DOLFIN_START_TEST( test_FiniteElementSpace )
{
  uint const deg_max = 2;
  std::vector<Family::Type> v;
  v.push_back(Family::DG);
  v.push_back(Family::CG);

  for (std::vector<Family::Type>::const_iterator it = v.begin(); it != v.end();
       ++it)
  {
    Family f(*it);
    uint d_min = f.degree_min();
    uint d_max = std::min(f.degree_max(), std::max(d_min, deg_max));
    ufl::Domain::Set domains = f.domains();

    for (ufl::Domain::Set::const_iterator dom_it = domains.begin();
         dom_it != domains.end(); ++dom_it)
    {
      Domain dom(*dom_it);
      ufl::Cell cell(dom);
      // uint const dim = cell.topological_dimension();

      CellType * ctype = CellType::create(cell);
      Mesh refcell( *ctype, EuclideanSpace( ctype->dim() ) );
      ctype->create_reference_cell(refcell);
      for (uint d = d_min; d <= d_max; ++d)
      {
        ufl::FiniteElement uflfem(*it, cell, d);
        message(uflfem.repr());
        //---
        ufl::FiniteElementSpace * factuflfem =
            ufl::FiniteElementSpace::create(uflfem.repr());
        delete factuflfem;
        //---
        dolfin::FiniteElement fem(uflfem);
        ck_assert(fem.signature() == uflfem.repr());
        //---
        ufl::FiniteElement uflfemd(Object::repr_t(fem.signature()));
        ck_assert(uflfem.repr() == uflfemd.repr());
        //---
        ufc::dofmap * ufcdm =
            ElementLibrary::create_dof_map(
                DofMap::make_signature(fem.signature()));
        DofMap dm(refcell, *ufcdm, true);
        ck_assert(dm.signature() == DofMap::make_signature(uflfem.repr()));
        //---
        FiniteElementSpace femspace(refcell, fem, dm, false);
      }
      delete ctype;
    }
  }
}
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
