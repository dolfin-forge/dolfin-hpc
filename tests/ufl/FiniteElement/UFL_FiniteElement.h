#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/mesh/Mesh.h>
#include <dolfin/ufl/UFLFamily.h>
#include <dolfin/ufl/UFLFiniteElement.h>
#include <dolfin/ufl/UFLVectorElement.h>

using namespace dolfin;

using ufl::Domain;
using ufl::Family;
using ufl::Object;
using ufl::Space;

//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_UFL_FiniteElement )
  {
    dolfin::uint const deg_max = 2;
    std::vector<Family::Type> v;
    v.push_back(Family::DG);
    v.push_back(Family::CG);
    
    for (std::vector<Family::Type>::const_iterator it = v.begin();
        it != v.end(); ++it)
    {
      Family f(*it);
      dolfin::uint d_min = f.degree_min();
      dolfin::uint d_max = std::min(f.degree_max(), std::max(d_min, deg_max));
      ufl::Domain::Set domains = f.domains();

      for (ufl::Domain::Set::const_iterator dom_it = domains.begin();
          dom_it != domains.end(); ++dom_it)
      {
        Domain dom(*dom_it);
        ufl::Cell cell(dom);
        for (dolfin::uint d = d_min; d <= d_max; ++d)
        {
          ufl::FiniteElement uflfem(*it, cell, d);
          uflfem.display();

          ufl::FiniteElementSpace * factuflfem =
              ufl::FiniteElementSpace::create(uflfem.repr());
          factuflfem->display();
          if (*factuflfem != uflfem)
          {
            init_failed += 1;
          }
          delete factuflfem;
        }
      }
    }
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------
  
#endif
