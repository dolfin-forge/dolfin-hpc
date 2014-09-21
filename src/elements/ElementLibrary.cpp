
#include <dolfin/elements/ElementLibrary.h>

#include "element_library.inc"

namespace dolfin
{

//-----------------------------------------------------------------------------
ufc::finite_element* ElementLibrary::create_finite_element(std::string const signature)
{
  return create_finite_element(signature.c_str());
}

//-----------------------------------------------------------------------------
ufc::dofmap* ElementLibrary::create_dof_map(std::string const signature)
{
  return create_dof_map(signature.c_str());
}



}
