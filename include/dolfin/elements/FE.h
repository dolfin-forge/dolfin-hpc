// Copyright (C) 2013 Aurélien Larcher
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2013-10-07 (merged from branch larcher)
// Last changed: 2013-10-07

#ifndef FE_H_
#define FE_H_

#include <dolfin/config/dolfin_config.h>
#include <dolfin/common/Array.h>

using dolfin::uint;
using dolfin::Array;

#include <ufc.h>

#include <map>

namespace FE
{

static char const INTERVAL     [] =  "interval";
static char const TRIANGLE     [] =  "triangle";
static char const QUADRILATERAL[] =  "quadrilateral";
static char const TETRAHEDRON  [] =  "tetrahedron";
static char const HEXAHEDRON   [] =  "hexahedron";

#if ENABLE_UFL
//UFC2.1
//-----------------------------------------------------------------------------
static char const FINITE_ELEMENT     [] =  "FiniteElement";
static char const VECTOR_ELEMENT     [] =  "VectorElement";
static char const MIXED_ELEMENT      [] =  "MixedElement";

//-----------------------------------------------------------------------------
static char const LAGRANGE     [] =  "Lagrange";
static char const LAGRANGE1DP1S[] =  "FiniteElement('Lagrange', Cell('interval', Space(1)), 1, None)";
static char const LAGRANGE1DP2S[] =  "FiniteElement('Lagrange', Cell('interval', Space(1)), 2, None)";
static char const LAGRANGE2DP1S[] =  "FiniteElement('Lagrange', Cell('triangle', Space(2)), 1, None)";
static char const LAGRANGE2DP2S[] =  "FiniteElement('Lagrange', Cell('triangle', Space(2)), 2, None)";
static char const LAGRANGE3DP1S[] =  "FiniteElement('Lagrange', Cell('tetrahedron', Space(3)), 1, None)";
static char const LAGRANGE3DP2S[] =  "FiniteElement('Lagrange', Cell('tetrahedron', Space(3)), 2, None)";
static char const LAGRANGE2DP1V[] =  "VectorElement('Lagrange', Cell('triangle', Space(2)), 1, 2, None)";
static char const LAGRANGE2DP2V[] =  "VectorElement('Lagrange', Cell('triangle', Space(2)), 2, 2, None)";
static char const LAGRANGE3DP1V[] =  "VectorElement('Lagrange', Cell('tetrahedron', Space(3)), 1, 3, None)";
static char const LAGRANGE3DP2V[] =  "VectorElement('Lagrange', Cell('tetrahedron', Space(3)), 2, 3, None)";

//-----------------------------------------------------------------------------
static char const DG           [] =  "Discontinuous Lagrange";
static char const DG1DP0S      [] =  "FiniteElement('Discontinuous Lagrange', Cell('interval', Space(1)), 0, None)";
static char const DG1DP1S      [] =  "FiniteElement('Discontinuous Lagrange', Cell('interval', Space(1)), 1, None)";
static char const DG1DP2S      [] =  "FiniteElement('Discontinuous Lagrange', Cell('interval', Space(1)), 2, None)";
static char const DG2DP0S      [] =  "FiniteElement('Discontinuous Lagrange', Cell('triangle', Space(2)), 0, None)";
static char const DG2DP1S      [] =  "FiniteElement('Discontinuous Lagrange', Cell('triangle', Space(2)), 1, None)";
static char const DG2DP2S      [] =  "FiniteElement('Discontinuous Lagrange', Cell('triangle', Space(2)), 2, None)";
static char const DG3DP0S      [] =  "FiniteElement('Discontinuous Lagrange', Cell('tetrahedron', Space(3)), 0, None)";
static char const DG3DP1S      [] =  "FiniteElement('Discontinuous Lagrange', Cell('tetrahedron', Space(3)), 1, None)";
static char const DG3DP2S      [] =  "FiniteElement('Discontinuous Lagrange', Cell('tetrahedron', Space(3)), 2, None)";
static char const DG2DP0V      [] =  "VectorElement('Discontinuous Lagrange', Cell('triangle', Space(2)), 0, 2, None)";
static char const DG2DP1V      [] =  "VectorElement('Discontinuous Lagrange', Cell('triangle', Space(2)), 1, 2, None)";
static char const DG2DP2V      [] =  "VectorElement('Discontinuous Lagrange', Cell('triangle', Space(2)), 2, 2, None)";
static char const DG3DP0V      [] =  "VectorElement('Discontinuous Lagrange', Cell('tetrahedron', Space(3)), 0, 3, None)";
static char const DG3DP1V      [] =  "VectorElement('Discontinuous Lagrange', Cell('tetrahedron', Space(3)), 1, 3, None)";
static char const DG3DP2V      [] =  "VectorElement('Discontinuous Lagrange', Cell('tetrahedron', Space(3)), 2, 3, None)";

//-----------------------------------------------------------------------------
static char const BDM          [] =  "Brezzi-Douglas-Marini";
static char const BDM2DP1      [] =  "FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', Space(2)), 1, None)";

#else
// UFC1.1
//-----------------------------------------------------------------------------
static char const FINITE_ELEMENT     [] =  "FiniteElement";
static char const VECTOR_ELEMENT     [] =  "VectorElement";
static char const MIXED_ELEMENT      [] =  "MixedElement";

//-----------------------------------------------------------------------------
static char const LAGRANGE     [] =  "Lagrange";
static char const LAGRANGE1DP1S[] =  "Lagrange finite element of degree 1 on a interval";
static char const LAGRANGE1DP2S[] =  "Lagrange finite element of degree 2 on a interval";
static char const LAGRANGE2DP1S[] =  "Lagrange finite element of degree 1 on a triangle";
static char const LAGRANGE2DP2S[] =  "Lagrange finite element of degree 2 on a triangle";
static char const LAGRANGE3DP1S[] =  "Lagrange finite element of degree 1 on a tetrahedron";
static char const LAGRANGE3DP2S[] =  "Lagrange finite element of degree 2 on a tetrahedron";
static char const LAGRANGE2DP1V[] =  "Mixed finite element: [Lagrange finite element of degree 1 on a triangle, Lagrange finite element of degree 1 on a triangle]";
static char const LAGRANGE2DP2V[] =  "Mixed finite element: [Lagrange finite element of degree 2 on a triangle, Lagrange finite element of degree 2 on a triangle]";
static char const LAGRANGE3DP1V[] =  "Mixed finite element: [Lagrange finite element of degree 1 on a tetrahedron, Lagrange finite element of degree 1 on a tetrahedron, Lagrange finite element of degree 1 on a tetrahedron]";
static char const LAGRANGE3DP2V[] =  "Mixed finite element: [Lagrange finite element of degree 2 on a tetrahedron, Lagrange finite element of degree 2 on a tetrahedron, Lagrange finite element of degree 2 on a tetrahedron]";

//-----------------------------------------------------------------------------
static char const DG           [] =  "Discontinuous Lagrange";
static char const DG1DP0S      [] =  "Discontinuous Lagrange finite element of degree 0 on a interval";
static char const DG1DP1S      [] =  "Discontinuous Lagrange finite element of degree 1 on a interval";
static char const DG1DP2S      [] =  "Discontinuous Lagrange finite element of degree 2 on a interval";
static char const DG2DP0S      [] =  "Discontinuous Lagrange finite element of degree 0 on a triangle";
static char const DG2DP1S      [] =  "Discontinuous Lagrange finite element of degree 1 on a triangle";
static char const DG2DP2S      [] =  "Discontinuous Lagrange finite element of degree 2 on a triangle";
static char const DG3DP0S      [] =  "Discontinuous Lagrange finite element of degree 0 on a tetrahedron";
static char const DG3DP1S      [] =  "Discontinuous Lagrange finite element of degree 1 on a tetrahedron";
static char const DG3DP2S      [] =  "Discontinuous Lagrange finite element of degree 0 on a tetrahedron";
static char const DG2DP0V      [] =  "Mixed finite element: [Discontinuous Lagrange finite element of degree 0 on a triangle, Discontinuous Lagrange finite element of degree 0 on a triangle]";
static char const DG2DP1V      [] =  "Mixed finite element: [Discontinuous Lagrange finite element of degree 1 on a triangle, Discontinuous Lagrange finite element of degree 1 on a triangle]";
static char const DG2DP2V      [] =  "Mixed finite element: [Discontinuous Lagrange finite element of degree 2 on a triangle, Discontinuous Lagrange finite element of degree 2 on a triangle]";
static char const DG3DP0V      [] =  "Mixed finite element: [Discontinuous Lagrange finite element of degree 0 on a tetrahedron, Discontinuous Lagrange finite element of degree 0 on a tetrahedron, Discontinuous Lagrange finite element of degree 0 on a tetrahedron]";
static char const DG3DP1V      [] =  "Mixed finite element: [Discontinuous Lagrange finite element of degree 1 on a tetrahedron, Discontinuous Lagrange finite element of degree 1 on a tetrahedron, Discontinuous Lagrange finite element of degree 1 on a tetrahedron]";
static char const DG3DP2V      [] =  "Mixed finite element: [Discontinuous Lagrange finite element of degree 2 on a tetrahedron, Discontinuous Lagrange finite element of degree 2 on a tetrahedron, Discontinuous Lagrange finite element of degree 2 on a tetrahedron]";

//-----------------------------------------------------------------------------
static char const BDM          [] =  "Brezzi-Douglas-Marini";
static char const BDM2DP1      [] =  "Brezzi-Douglas-Marini finite element of degree 1 on a triangle";

#endif

//-----------------------------------------------------------------------------
class FunctionSpace
{

public:

  //---------------------------------------------------------------------------
  enum Type { FiniteElement, VectorElement, MixedElement };
  //---------------------------------------------------------------------------
  /// List of finite element spaces
  static Array<std::string> const Names;
  //---------------------------------------------------------------------------
  static std::string const type2string(Type const& t);
  static Type const string2type(std::string const& s);

private:

  static std::map<Type, std::string> const TypesOfFunctionSpaces;
  static std::map<std::string, Type> const NamesOfFunctionSpaces;

};

inline std::string const FunctionSpace::type2string(FunctionSpace::Type const& t)
{
  return TypesOfFunctionSpaces.find(t)->second;
}

inline FunctionSpace::Type const FunctionSpace::string2type(std::string const& s)
{
  return NamesOfFunctionSpaces.find(s)->second;
}

//-----------------------------------------------------------------------------
class Family
{

public:

  //---------------------------------------------------------------------------
  enum Type { Lagrange, DG, BDM };
  //---------------------------------------------------------------------------
  /// List of finite element families
  static Array<std::string> const Names;
  //---------------------------------------------------------------------------
  static std::string const type2string(Type const& t);
  static Type const string2type(std::string const& s);

private:

  static std::map<Type, std::string> const TypesOfFamilies;
  static std::map<std::string, Type> const NamesOfFamilies;

};

inline std::string const Family::type2string(Family::Type const& t)
{
  return TypesOfFamilies.find(t)->second;
}

inline Family::Type const Family::string2type(std::string const& s)
{
  return NamesOfFamilies.find(s)->second;
}

//-----------------------------------------------------------------------------
class Cell
{

public:

  //---------------------------------------------------------------------------
  //typedef ufc::shape Type;
  enum Type {interval, triangle, quadrilateral, tetrahedron, hexahedron};
  //---------------------------------------------------------------------------
  /// List of cell types
  static Array<std::string> const Names;
  //---------------------------------------------------------------------------
  static std::string const type2string(Type const& t);
  static Type const string2type(std::string const& s);

private:

  static std::map<Type, std::string> const TypesOfCells;
  static std::map<std::string, Type> const NamesOfCells;

};

inline std::string const Cell::type2string(Cell::Type const& t)
{
  return TypesOfCells.find(t)->second;
}

inline Cell::Type const Cell::string2type(std::string const& s)
{
  return NamesOfCells.find(s)->second;
}

//-----------------------------------------------------------------------------
class attributes {

public:

      FunctionSpace::Type const type;
      Family::Type const family;
      Cell::Type const shape;
      uint const space_dim;
      uint const degree;
      uint const value_dim;

      attributes(FunctionSpace::Type tp, Family::Type fm, Cell::Type sh, uint sp, uint dg, uint vl) :
        type(tp), family(fm), shape(sh), space_dim(sp), degree(dg), value_dim(vl)
      {

      }

      attributes(std::string const& tp, std::string const& fm, std::string const& sh, uint sp, uint dg, uint vl) :
        type(FunctionSpace::string2type(tp)), family(Family::string2type(fm)), shape(Cell::string2type(sh)), space_dim(sp), degree(dg), value_dim(vl)
      {

      }

      void display()
      {
        std::cout << std::endl;
      }

};

//-----------------------------------------------------------------------------
std::string const get_signature(FE::FunctionSpace::Type const space,
                                FE::Family::Type const family,
                                FE::Cell::Type const cell,
                                uint const space_dim,
                                uint const degree,
                                uint const value_dim = 0);

//-----------------------------------------------------------------------------
attributes const get_attributes(std::string const signature);

}

#endif /* FE_H_ */
