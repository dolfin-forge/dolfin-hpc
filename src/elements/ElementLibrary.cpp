#include <dolfin/config/dolfin_config.h>
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

#if ENABLE_UFL

//-----------------------------------------------------------------------------
ElementLibrary::ElementsTable const init_ElementsTable()
{
  ElementLibrary::ElementsTable ret;
  ret.insert(ElementLibrary::ElementsItem( FE::LAGRANGE1DP1S, FE::attributes(FE::FunctionSpace::FiniteElement, FE::Family::Lagrange, FE::Cell::interval, 1, 1, 1)));
  ret.insert(ElementLibrary::ElementsItem( FE::LAGRANGE1DP2S, FE::attributes(FE::FunctionSpace::FiniteElement, FE::Family::Lagrange, FE::Cell::interval, 1, 2, 1)));
  ret.insert(ElementLibrary::ElementsItem( FE::LAGRANGE2DP1S, FE::attributes(FE::FunctionSpace::FiniteElement, FE::Family::Lagrange, FE::Cell::triangle, 2, 1, 1)));
  ret.insert(ElementLibrary::ElementsItem( FE::LAGRANGE2DP2S, FE::attributes(FE::FunctionSpace::FiniteElement, FE::Family::Lagrange, FE::Cell::triangle, 2, 2, 1)));
  ret.insert(ElementLibrary::ElementsItem( FE::LAGRANGE3DP1S, FE::attributes(FE::FunctionSpace::FiniteElement, FE::Family::Lagrange, FE::Cell::tetrahedron, 3, 1, 1)));
  ret.insert(ElementLibrary::ElementsItem( FE::LAGRANGE3DP2S, FE::attributes(FE::FunctionSpace::FiniteElement, FE::Family::Lagrange, FE::Cell::tetrahedron, 3, 2, 1)));
  ret.insert(ElementLibrary::ElementsItem( FE::LAGRANGE2DP1V, FE::attributes(FE::FunctionSpace::VectorElement, FE::Family::Lagrange, FE::Cell::triangle, 2, 1, 2)));
  ret.insert(ElementLibrary::ElementsItem( FE::LAGRANGE2DP2V, FE::attributes(FE::FunctionSpace::VectorElement, FE::Family::Lagrange, FE::Cell::triangle, 2, 2, 2)));
  ret.insert(ElementLibrary::ElementsItem( FE::LAGRANGE3DP1V, FE::attributes(FE::FunctionSpace::VectorElement, FE::Family::Lagrange, FE::Cell::tetrahedron, 3, 1, 3)));
  ret.insert(ElementLibrary::ElementsItem( FE::LAGRANGE3DP2V, FE::attributes(FE::FunctionSpace::VectorElement, FE::Family::Lagrange, FE::Cell::tetrahedron, 3, 2, 3)));
  ret.insert(ElementLibrary::ElementsItem( FE::DG1DP0S, FE::attributes(FE::FunctionSpace::FiniteElement, FE::Family::DG, FE::Cell::interval, 1, 0, 1)));
  ret.insert(ElementLibrary::ElementsItem( FE::DG1DP1S, FE::attributes(FE::FunctionSpace::FiniteElement, FE::Family::DG, FE::Cell::interval, 1, 1, 1)));
  ret.insert(ElementLibrary::ElementsItem( FE::DG1DP2S, FE::attributes(FE::FunctionSpace::FiniteElement, FE::Family::DG, FE::Cell::interval, 1, 2, 1)));
  ret.insert(ElementLibrary::ElementsItem( FE::DG2DP0S, FE::attributes(FE::FunctionSpace::FiniteElement, FE::Family::DG, FE::Cell::triangle, 2, 0, 1)));
  ret.insert(ElementLibrary::ElementsItem( FE::DG2DP1S, FE::attributes(FE::FunctionSpace::FiniteElement, FE::Family::DG, FE::Cell::triangle, 2, 1, 1)));
  ret.insert(ElementLibrary::ElementsItem( FE::DG2DP2S, FE::attributes(FE::FunctionSpace::FiniteElement, FE::Family::DG, FE::Cell::triangle, 2, 2, 1)));
  ret.insert(ElementLibrary::ElementsItem( FE::DG3DP0S, FE::attributes(FE::FunctionSpace::FiniteElement, FE::Family::DG, FE::Cell::tetrahedron, 3, 0, 1)));
  ret.insert(ElementLibrary::ElementsItem( FE::DG3DP1S, FE::attributes(FE::FunctionSpace::FiniteElement, FE::Family::DG, FE::Cell::tetrahedron, 3, 1, 1)));
  ret.insert(ElementLibrary::ElementsItem( FE::DG3DP2S, FE::attributes(FE::FunctionSpace::FiniteElement, FE::Family::DG, FE::Cell::tetrahedron, 3, 2, 1)));
  ret.insert(ElementLibrary::ElementsItem( FE::DG2DP0V, FE::attributes(FE::FunctionSpace::VectorElement, FE::Family::DG, FE::Cell::triangle, 2, 0, 2)));
  ret.insert(ElementLibrary::ElementsItem( FE::DG2DP1V, FE::attributes(FE::FunctionSpace::VectorElement, FE::Family::DG, FE::Cell::triangle, 2, 1, 2)));
  ret.insert(ElementLibrary::ElementsItem( FE::DG2DP2V, FE::attributes(FE::FunctionSpace::VectorElement, FE::Family::DG, FE::Cell::triangle, 2, 2, 2)));
  ret.insert(ElementLibrary::ElementsItem( FE::DG3DP0V, FE::attributes(FE::FunctionSpace::VectorElement, FE::Family::DG, FE::Cell::tetrahedron, 3, 0, 3)));
  ret.insert(ElementLibrary::ElementsItem( FE::DG3DP1V, FE::attributes(FE::FunctionSpace::VectorElement, FE::Family::DG, FE::Cell::tetrahedron, 3, 1, 3)));
  ret.insert(ElementLibrary::ElementsItem( FE::DG3DP2V, FE::attributes(FE::FunctionSpace::VectorElement, FE::Family::DG, FE::Cell::tetrahedron, 3, 2, 3)));
  ret.insert(ElementLibrary::ElementsItem( FE::BDM2DP1, FE::attributes(FE::FunctionSpace::FiniteElement, FE::Family::BDM, FE::Cell::triangle, 2, 1, 2)));
  return ret;
}

//-----------------------------------------------------------------------------
ElementLibrary::ElementsTable const ElementLibrary::Elements = init_ElementsTable();

//-----------------------------------------------------------------------------
FE::attributes const ElementLibrary::get_attributes(const char* signature)
{
  return get_attributes(std::string(signature));
}

//-----------------------------------------------------------------------------
FE::attributes const ElementLibrary::get_attributes(std::string const signature)
{
  ElementsTable::const_iterator it = Elements.find(signature);
  if (it == Elements.end())
  {
    return FE::get_attributes(signature);
  }
  else
  {
    return it->second;
  }
}

#endif

}
