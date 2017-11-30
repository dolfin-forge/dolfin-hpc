#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/mesh/UnitSquare.h>
#include <dolfin/ufl/UFLFiniteElement.h>
#include <dolfin/ufl/UFLVectorElement.h>
#include <dolfin/ufl/UFLMixedElement.h>

//-----------------------------------------------------------------------------
START_TEST( test_UFL_MixedElement )
{

  dolfin::UnitSquare mesh(4,4);

  ufl::VectorElement Uspace(ufl::Family::CG, mesh.type(), 2,
                            mesh.geometry().dim());
  Uspace.display();

  ufl::FiniteElement Pspace(ufl::Family::CG, mesh.type(), 1);
  Pspace.display();

  ufl::FiniteElementSpace::List spaces;
  spaces.push_back(&Uspace);
  spaces.push_back(&Pspace);

  ufl::MixedElement UPspace(spaces);
  UPspace.display();

  ufl::MixedElement UPspaceFormRepr(UPspace.repr());
  UPspaceFormRepr.display();

  UPspaceFormRepr.sub_elements()[0]->display();
  UPspaceFormRepr.sub_elements()[1]->display();

  ck_assert_msg(UPspace == UPspaceFormRepr,
		"Representation string differ");

}END_TEST
//-----------------------------------------------------------------------------

#endif
