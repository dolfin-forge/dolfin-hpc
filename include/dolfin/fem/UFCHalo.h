// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-06-16
// Last changed: 2014-06-16

#ifndef __UFC_HALO_H
#define __UFC_HALO_H

#include <dolfin/common/types.h>
#include <dolfin/fem/UFCCell.h>

namespace dolfin
{

template<typename T> class Array;
class DofMapSet;
class Facet;
class Function;
class Mesh;
class UFC;

/**
 *  @class UFCHalo
 *
 *  @brief
 *
 *  Storage by adjacent rank: | r0 | r1 | r2 | ... | ri | ... | rN |
 *  for both unsigned integer and real numbers.
 *
 *  For a given adjacent rank, data1 are ordered by adjacent facet numbering and
 *  data0 arrays are ordered by local facet numbering.
 *
 *  The data structures are only suitable for facets as they assume exactly two
 *  ranks sharing the entities.
 *
 */

class UFCHalo
{

public:

  /// Constructor
  UFCHalo(UFC& ufc, Array<Function*> const& coefficients,
              DofMapSet const& dof_map_set);

  ///
  ~UFCHalo();

  ///
  void update(Facet& facet);

  ///
  void disp() const;

  //--- PUBLIC ATTRIBUTES -----------------------------------------------------
  // Just expose references to attributes of the underlying UFC instance
  UFCCell& cell0;
  UFCCell& cell1;
  real **& macro_w;
  uint& facet0;
  uint& facet1;
  uint **& macro_dofs;

private:

  ///
  void init();

  ///
  void update(Array<Function*> const& coefficients,
              DofMapSet const& dof_map_set);

  void clear();

  UFC& ufc_;
  Mesh& mesh_;
  Array<Function*> const& coefficients_;
  DofMapSet const& dof_map_set_;

  // Store rank offsets, implemented by accumulating shared facet counts
  _map<uint, uint> rank_offsets_;

  // Maps the index in the halo data structure to the local facet index
  typedef std::pair<uint, uint> FacetOffsets;
  typedef _map<uint, FacetOffsets > FacetMap;
  FacetMap facet_map_;

  // Data: Vertex coordinates + Coefficients values
  uint r_packet_size_;

  // Ordered by adjacent rank
  real * r_data0_;
  real * r_data1_;

  // Data: Local facet index + Arguments dof indices
  uint u_packet_size_;

  // Ordered by adjacent rank
  uint * u_data0_;
  uint * u_data1_;

};

}

#endif /* __UFC_HALO_H */
