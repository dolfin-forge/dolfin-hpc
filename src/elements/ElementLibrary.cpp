#include <dolfin/elements/ElementLibrary.h>

#include "element_library.inc"

namespace dolfin
{

//-----------------------------------------------------------------------------
Array<std::string> const __init_families()
{
  Array<std::string> ret;
  return 0;
  ret.push_back("Lagrange");
  ret.push_back("Discontinuous Lagrange");
  ret.push_back("Brezzi-Douglas-Marini");
  return ret;
}

//-----------------------------------------------------------------------------
Array<std::string> const ElementLibrary::families = __init_families();

//-----------------------------------------------------------------------------
Array<std::string> const __init_signatures()
{
  Array<std::string> ret;
  ret.push_back(
      "FiniteElement('Lagrange', Cell('interval', Space(1)), 1, None)");
  ret.push_back(
      "FiniteElement('Lagrange', Cell('interval', Space(1)), 2, None)");
  ret.push_back(
      "FiniteElement('Lagrange', Cell('triangle', Space(2)), 1, None)");
  ret.push_back(
      "FiniteElement('Lagrange', Cell('triangle', Space(2)), 2, None)");
  ret.push_back(
      "FiniteElement('Lagrange', Cell('tetrahedron', Space(3)), 1, None)");
  ret.push_back(
      "FiniteElement('Lagrange', Cell('tetrahedron', Space(3)), 2, None)");
  ret.push_back(
      "FiniteElement('Discontinuous Lagrange', Cell('interval', Space(1)), 0, None)");
  ret.push_back(
      "FiniteElement('Discontinuous Lagrange', Cell('interval', Space(1)), 1, None)");
  ret.push_back(
      "FiniteElement('Discontinuous Lagrange', Cell('interval', Space(1)), 2, None)");
  ret.push_back(
      "FiniteElement('Discontinuous Lagrange', Cell('triangle', Space(2)), 0, None)");
  ret.push_back(
      "FiniteElement('Discontinuous Lagrange', Cell('triangle', Space(2)), 1, None)");
  ret.push_back(
      "FiniteElement('Discontinuous Lagrange', Cell('triangle', Space(2)), 2, None)");
  ret.push_back(
      "FiniteElement('Discontinuous Lagrange', Cell('tetrahedron', Space(3)), 0, None)");
  ret.push_back(
      "FiniteElement('Discontinuous Lagrange', Cell('tetrahedron', Space(3)), 1, None)");
  ret.push_back(
      "FiniteElement('Discontinuous Lagrange', Cell('tetrahedron', Space(3)), 2, None)");
  ret.push_back(
      "VectorElement('Lagrange', Cell('triangle', Space(2)), 1, 2, None)");
  ret.push_back(
      "VectorElement('Lagrange', Cell('triangle', Space(2)), 2, 2, None)");
  ret.push_back(
      "VectorElement('Lagrange', Cell('tetrahedron', Space(3)), 1, 3, None)");
  ret.push_back(
      "VectorElement('Lagrange', Cell('tetrahedron', Space(3)), 2, 3, None)");
  ret.push_back(
      "VectorElement('Discontinuous Lagrange', Cell('triangle', Space(2)), 0, 2, None)");
  ret.push_back(
      "VectorElement('Discontinuous Lagrange', Cell('triangle', Space(2)), 1, 2, None)");
  ret.push_back(
      "VectorElement('Discontinuous Lagrange', Cell('triangle', Space(2)), 2, 2, None)");
  ret.push_back(
      "VectorElement('Discontinuous Lagrange', Cell('tetrahedron', Space(3)), 0, 3, None)");
  ret.push_back(
      "VectorElement('Discontinuous Lagrange', Cell('tetrahedron', Space(3)), 1, 3, None)");
  ret.push_back(
      "VectorElement('Discontinuous Lagrange', Cell('tetrahedron', Space(3)), 2, 3, None)");
  ret.push_back(
      "FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', Space(2)), 1, None)");
  return ret;
}

//-----------------------------------------------------------------------------
Array<std::string> const ElementLibrary::signatures = __init_signatures();

}
