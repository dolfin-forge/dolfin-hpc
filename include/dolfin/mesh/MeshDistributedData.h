// Copyright (C) 2008 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_MESH_DISTRIBUTED_DATA_H
#define __DOLFIN_MESH_DISTRIBUTED_DATA_H

#include <dolfin/common/DistributedData.h>
#include <dolfin/common/types.h>

namespace dolfin
{

class Mesh;
class MeshEntity;

/**
 *  @class  MeshDistributedData
 *
 *  @brief  Implements data structures for distributed mesh topology using an
 *          array of distributed data for each topological dimension.
 */

class MeshDistributedData
{

public:
  /// Constructor
  MeshDistributedData( uint dim );

  /// Copy constructor
  MeshDistributedData( MeshDistributedData const & other );

  /// Destructor
  ~MeshDistributedData();

  /// Assignment
  auto operator=( MeshDistributedData const & other ) -> MeshDistributedData &;

  /// Equality
  auto operator==( MeshDistributedData const & other ) const -> bool;

  /// Non-equality
  auto operator!=( MeshDistributedData const & other ) const -> bool;

  /// Access to distributed data
  inline auto operator[]( uint dim ) -> DistributedData &
  {
    dolfin_assert( dim <= dim_ );
    return data_[dim];
  }

  /// Access to distributed data (const)
  inline auto operator[]( uint dim ) const -> DistributedData const &
  {
    dolfin_assert( dim <= dim_ );
    return data_[dim];
  }

  /// Return topological dimension of the distributed data
  auto dim() const -> uint;

  /// Display basic information
  void disp() const;

private:
  /// Clear
  void clear();

  // Topological dimensions
  uint dim_;

  // Distributed data
  std::vector< DistributedData > data_;
};

} /* namespace dolfin */

#endif /* __DOLFIN_MESH_DISTRIBUTED_DATA_H */
