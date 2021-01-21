// This class was imported from UNICORN and did not possess any copyright.
// It was provided with a second argument 'p' to avoid the same mistakes as in
// the UNICORN implementation where p could have the value 2 or 4 (or even more)
// at the same time in different parts of the solver.

#ifndef MESH_QUALITY_FUNCTION_H
#define MESH_QUALITY_FUNCTION_H

#include <dolfin/fem/UFCCell.h>
#include <dolfin/function/Function.h>
#include <dolfin/mesh/MeshQuality.h>

namespace dolfin
{

/**
 *  DOCUMENTATION:
 *
 *  @class
 *
 *  @brief
 *
 *  Cell quality: q_K
 *
 *  The mesh quality function is defined as 1./q_K^p
 *
 *  with p = 2 by default (the bigger, the sharper).
 *
 *
 */


class MeshQualityFunction : public Function
{

public:

  // The default value for p ought to be 2
  MeshQualityFunction(Mesh& mesh, uint p);

  //
  ~MeshQualityFunction() override = default;

  //
  auto rank() const -> uint override;

  //
  auto dim(uint i) const -> uint override;

  //
  void evaluate(real* values, const real* x, const ufc::cell& cell) const override;

private:

  uint const p_;
  MeshQuality mqual_;
};

//-----------------------------------------------------------------------------
inline auto MeshQualityFunction::rank() const -> uint
{
  return 0;
}

//-----------------------------------------------------------------------------
inline auto MeshQualityFunction::dim( uint ) const -> uint
{
  return 1;
}

//-----------------------------------------------------------------------------
inline void MeshQualityFunction::evaluate( real * values,
                                           const real *,
                                           const ufc::cell & cell ) const
{
  UFCCell const & ufc_cell = static_cast< UFCCell const & >( cell );
  real const      qK       = mqual_.mean_ratio( *ufc_cell );
  values[0]                = 1.0 / std::pow( qK, static_cast< real >( p_ ) );
}

//-----------------------------------------------------------------------------

}

#endif /* MESH_QUALITY_FUNCTION_H */
