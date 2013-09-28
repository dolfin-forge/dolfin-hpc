#ifndef FE_H_
#define FE_H_

#include <dolfin/common/Array.h>

using dolfin::Array;

namespace FE
{

//-----------------------------------------------------------------------------
char const * signature(std::string family, unsigned int rank, unsigned int nsdim, unsigned int degree);

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

static char const BDM          [] =  "Brezzi-Douglas-Marini";
static char const BDM1DP2      [] =  "FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', Space(2)), 1, None)";

//-----------------------------------------------------------------------------
Array<std::string> const init_families();
static Array<std::string> const Families = init_families();

//-----------------------------------------------------------------------------
Array<std::string> const init_elements();
static Array<std::string> const Elements = init_elements();

}

#endif /* FE_H_ */
