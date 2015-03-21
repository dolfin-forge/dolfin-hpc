//
//
//

#include <dolfin/mesh/MeshQualityFunction.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
MeshQualityFunction::MeshQualityFunction(Mesh& mesh, uint p) :
    Function(mesh),
    p_(p),
    mqual_(mesh)
{
}

//-----------------------------------------------------------------------------
MeshQualityFunction::~MeshQualityFunction()
{
}

//-----------------------------------------------------------------------------
uint MeshQualityFunction::rank() const
{
  return 0;
}

//-----------------------------------------------------------------------------
uint MeshQualityFunction::dim(uint i) const
{
  return 1;
}

//-----------------------------------------------------------------------------
void MeshQualityFunction::eval(real* values, const real* x) const
{
  //FIXME: should we have to cast?
  Cell& c = const_cast<Cell&>(cell());
  real const qK = mqual_.mean_ratio(c);
  values[0] = 1.0 / std::pow(qK, static_cast<real>(p_));
}

}

