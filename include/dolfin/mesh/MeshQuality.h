
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
  bool is_inverted(uint& first);

  /// Cell quality based on mean-ratio metric
  real mean_ratio(Cell const& cell) const;

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
  static real reduceMinReal(real val);
  static real reduceMaxReal(real val);
  static real reduceAvgReal(real val);

  MeshValues<int, Cell> orientation_;
  mutable EquiAffineMapping mapping_;
  Point bbox_min_;
  Point bbox_max_;

};

}

#endif
