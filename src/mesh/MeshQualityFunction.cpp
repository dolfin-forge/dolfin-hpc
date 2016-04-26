//
//
//

#include <dolfin/mesh/MeshQualityFunction.h>

#include <dolfin/fem/UFCCell.h>

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
void MeshQualityFunction::evaluate(real* values, const real* x,
                                   const ufc::cell& cell) const
{
  UFCCell const& ufc_cell = static_cast<UFCCell const&>(cell);
  real const qK = mqual_.mean_ratio(*(ufc_cell.cell));
  values[0] = 1.0 / std::pow(qK, static_cast<real>(p_));
}

//-----------------------------------------------------------------------------

}

