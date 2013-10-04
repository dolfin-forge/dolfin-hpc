#include <dolfin/config/dolfin_config.h>
#include <dolfin/elements/ElementLibrary.h>

#include "element_library.inc"

namespace dolfin
{

#if ENABLE_UFL

//-----------------------------------------------------------------------------
Array<std::string> const init_types()
{
  Array<std::string> ret;
  ret.push_back(FE::FINITE_ELEMENT);
  ret.push_back(FE::VECTOR_ELEMENT);
  ret.push_back(FE::MIXED_ELEMENT);
  return ret;
}

Array<std::string> const ElementLibrary::Types = init_types();

//-----------------------------------------------------------------------------
Array<std::string> const init_families()
{
  Array<std::string> ret;
  ret.push_back(FE::LAGRANGE);
  ret.push_back(FE::DG);
  ret.push_back(FE::BDM);
  return ret;
}

Array<std::string> const ElementLibrary::Families = init_families();

//-----------------------------------------------------------------------------
ElementLibrary::ElementsTable const init_elements()
{
  ElementLibrary::ElementsTable ret;
  ret.insert(ElementLibrary::ElementsItem( FE::LAGRANGE1DP1S, FE::attributes(FE::FINITE_ELEMENT, FE::LAGRANGE, ufc::interval, 1, 1, 1)));
  ret.insert(ElementLibrary::ElementsItem( FE::LAGRANGE1DP2S, FE::attributes(FE::FINITE_ELEMENT, FE::LAGRANGE, ufc::interval, 1, 2, 1)));
  ret.insert(ElementLibrary::ElementsItem( FE::LAGRANGE2DP1S, FE::attributes(FE::FINITE_ELEMENT, FE::LAGRANGE, ufc::triangle, 2, 1, 1)));
  ret.insert(ElementLibrary::ElementsItem( FE::LAGRANGE2DP2S, FE::attributes(FE::FINITE_ELEMENT, FE::LAGRANGE, ufc::triangle, 2, 2, 1)));
  ret.insert(ElementLibrary::ElementsItem( FE::LAGRANGE3DP1S, FE::attributes(FE::FINITE_ELEMENT, FE::LAGRANGE, ufc::tetrahedron, 3, 1, 1)));
  ret.insert(ElementLibrary::ElementsItem( FE::LAGRANGE3DP2S, FE::attributes(FE::FINITE_ELEMENT, FE::LAGRANGE, ufc::tetrahedron, 3, 2, 1)));
  ret.insert(ElementLibrary::ElementsItem( FE::LAGRANGE2DP1V, FE::attributes(FE::VECTOR_ELEMENT, FE::LAGRANGE, ufc::triangle, 2, 1, 2)));
  ret.insert(ElementLibrary::ElementsItem( FE::LAGRANGE2DP2V, FE::attributes(FE::VECTOR_ELEMENT, FE::LAGRANGE, ufc::triangle, 2, 2, 2)));
  ret.insert(ElementLibrary::ElementsItem( FE::LAGRANGE3DP1V, FE::attributes(FE::VECTOR_ELEMENT, FE::LAGRANGE, ufc::tetrahedron, 3, 1, 3)));
  ret.insert(ElementLibrary::ElementsItem( FE::LAGRANGE3DP2V, FE::attributes(FE::VECTOR_ELEMENT, FE::LAGRANGE, ufc::tetrahedron, 3, 2, 3)));
  ret.insert(ElementLibrary::ElementsItem( FE::DG1DP0S, FE::attributes(FE::FINITE_ELEMENT, FE::DG, ufc::interval, 1, 0, 1)));
  ret.insert(ElementLibrary::ElementsItem( FE::DG1DP1S, FE::attributes(FE::FINITE_ELEMENT, FE::DG, ufc::interval, 1, 1, 1)));
  ret.insert(ElementLibrary::ElementsItem( FE::DG1DP2S, FE::attributes(FE::FINITE_ELEMENT, FE::DG, ufc::interval, 1, 2, 1)));
  ret.insert(ElementLibrary::ElementsItem( FE::DG2DP0S, FE::attributes(FE::FINITE_ELEMENT, FE::DG, ufc::triangle, 2, 0, 1)));
  ret.insert(ElementLibrary::ElementsItem( FE::DG2DP1S, FE::attributes(FE::FINITE_ELEMENT, FE::DG, ufc::triangle, 2, 1, 1)));
  ret.insert(ElementLibrary::ElementsItem( FE::DG2DP2S, FE::attributes(FE::FINITE_ELEMENT, FE::DG, ufc::triangle, 2, 2, 1)));
  ret.insert(ElementLibrary::ElementsItem( FE::DG3DP0S, FE::attributes(FE::FINITE_ELEMENT, FE::DG, ufc::tetrahedron, 3, 0, 1)));
  ret.insert(ElementLibrary::ElementsItem( FE::DG3DP1S, FE::attributes(FE::FINITE_ELEMENT, FE::DG, ufc::tetrahedron, 3, 1, 1)));
  ret.insert(ElementLibrary::ElementsItem( FE::DG3DP2S, FE::attributes(FE::FINITE_ELEMENT, FE::DG, ufc::tetrahedron, 3, 2, 1)));
  ret.insert(ElementLibrary::ElementsItem( FE::DG2DP0V, FE::attributes(FE::VECTOR_ELEMENT, FE::DG, ufc::triangle, 2, 0, 2)));
  ret.insert(ElementLibrary::ElementsItem( FE::DG2DP1V, FE::attributes(FE::VECTOR_ELEMENT, FE::DG, ufc::triangle, 2, 1, 2)));
  ret.insert(ElementLibrary::ElementsItem( FE::DG2DP2V, FE::attributes(FE::VECTOR_ELEMENT, FE::DG, ufc::triangle, 2, 2, 2)));
  ret.insert(ElementLibrary::ElementsItem( FE::DG3DP0V, FE::attributes(FE::VECTOR_ELEMENT, FE::DG, ufc::tetrahedron, 3, 0, 3)));
  ret.insert(ElementLibrary::ElementsItem( FE::DG3DP1V, FE::attributes(FE::VECTOR_ELEMENT, FE::DG, ufc::tetrahedron, 3, 1, 3)));
  ret.insert(ElementLibrary::ElementsItem( FE::DG3DP2V, FE::attributes(FE::VECTOR_ELEMENT, FE::DG, ufc::tetrahedron, 3, 2, 3)));
  ret.insert(ElementLibrary::ElementsItem( FE::BDM2DP1, FE::attributes(FE::FINITE_ELEMENT, FE::BDM, ufc::triangle, 2, 1, 2)));
  return ret;
}

//-----------------------------------------------------------------------------
ElementLibrary::ElementsTable const ElementLibrary::Elements = init_elements();

#endif

}
