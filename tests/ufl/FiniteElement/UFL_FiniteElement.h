#ifdef HAVE_CHECK

#include <check.h>

#include <dolfin/config/dolfin_config.h>

#include <dolfin/log/log.h>
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
START_TEST( test_UFL_FiniteElement )
{
  int init_failed = 0;
  
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
	begin("Creating UFLFiniteElement:");
	ufl::FiniteElement uflfem(*it, cell, d);
	message(uflfem.repr());
	uflfem.display();
	end();
	skip();
	
	begin("Creating UFLFiniteElement from factory function:");
	ufl::FiniteElementBase * factuflfem = ufl::FiniteElementBase::create(uflfem.repr());
	message(factuflfem->repr());
	factuflfem->display();
	if(*factuflfem != uflfem)
	{
	  init_failed += 1;
	}
	delete factuflfem;
	end();
	skip();
	
      }
    }
  }
  
  fail_unless( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------

#endif
