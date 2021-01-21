
#ifndef DOLFIN_MESH_QUALITY_H
#define DOLFIN_MESH_QUALITY_H

#include <dolfin/mesh/MeshDependent.h>

#include <dolfin/common/constants.h>
#include <dolfin/mesh/EquiAffineMapping.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshValues.h>

namespace dolfin
{

class MeshQuality : public MeshDependent
{

public:

  /// Constructor
  MeshQuality(Mesh& mesh);

  ///
  auto is_inverted(uint& first) -> bool;

  /// Cell quality based on mean-ratio metric
  auto mean_ratio(Cell const& cell) const -> real;

  ///
  void compute();

  ///
  void disp();

  //--- PUBLIC ATTRIBUTES -----------------------------------------------------

  real mu_min;
  real mu_max;
  real mu_avg;
  real h_min;
  real h_max;
  real h_avg;
  real vol_min;
  real vol_max;

private:

  ///
  static auto reduceMinReal(real val) -> real;
  static auto reduceMaxReal(real val) -> real;
  static auto reduceAvgReal(real val) -> real;

  MeshValues<int, Cell> orientation_;
  mutable EquiAffineMapping mapping_;
  Point bbox_min_;
  Point bbox_max_;

};

}

#endif
