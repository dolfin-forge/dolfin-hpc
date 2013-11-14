// Copyright (C) 2013 Aurélien Larcher
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2013-10-07 (merged from branch larcher)
// Last changed: 2013-10-07

#ifndef FE_H_
#define FE_H_

#include <dolfin/config/dolfin_config.h>
#include <dolfin/common/Array.h>

using dolfin::Array;
using dolfin::uint;

#include <ufc.h>

#include <map>

namespace FE
{

//-----------------------------------------------------------------------------
struct attributes {
      char const * type;
      char const * family;
      ufc::shape shape;
      uint space;
      uint degree;
      uint value;

      attributes(char const * tp, char const * fm, ufc::shape sh, uint sp, uint dg, uint vl) :
        type(tp),family(fm), shape(sh), space(sp), degree(dg), value(vl)
      {

      }
};

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

}

#endif /* FE_H_ */
