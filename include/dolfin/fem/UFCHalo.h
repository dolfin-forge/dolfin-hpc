// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

// This class is the first implementation allowing assembly of halo facets
// contributions like jumps.

#ifndef __DOLFIN_UFC_HALO_H
#define __DOLFIN_UFC_HALO_H

#include <dolfin/common/types.h>
#include <dolfin/fem/UFCCell.h>

namespace dolfin
{

class DofMapSet;
class Facet;
class Coefficient;
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
  UFCHalo( UFC &                                ufc,
           std::vector< Coefficient * > const & coefficients,
           DofMapSet const &                    dof_map_set );

  ///
  ~UFCHalo();

  ///
  void update( Facet & facet );

  ///
  void disp() const;

  //--- PUBLIC ATTRIBUTES -----------------------------------------------------
  // Just expose references to attributes of the underlying UFC instance
  UFCCell &  cell0;
  UFCCell &  cell1;
  real **&   macro_w;
  size_t &   facet0;
  size_t &   facet1;
  size_t **& macro_dofs;

private:
  ///
  void init();

  ///
  void update( std::vector< Coefficient * > const & coefficients,
               DofMapSet const &                    dof_map_set );

  void clear();

  UFC &                                ufc_;
  Mesh &                               mesh_;
  std::vector< Coefficient * > const & coefficients_;
  DofMapSet const &                    dof_map_set_;

  // Store rank offsets, implemented by accumulating shared facet counts
  _map< size_t, size_t > rank_offsets_;

  // Maps the index in the halo data structure to the local facet index
  typedef std::pair< size_t, size_t >  FacetOffsets;
  typedef _map< size_t, FacetOffsets > FacetMap;
  FacetMap                             facet_map_;

  // Data: Vertex coordinates + Coefficients values
  size_t r_packet_size_;

  // Ordered by adjacent rank
  real * r_data0_;
  real * r_data1_;

  // Data: Local facet index + Arguments dof indices
  size_t u_packet_size_;

  // Ordered by adjacent rank
  size_t * u_data0_;
  size_t * u_data1_;
};

} /* namespace dolfin */

#endif /* __DOLFIN_UFC_HALO_H */
