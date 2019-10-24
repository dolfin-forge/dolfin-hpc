// Copyright (C) 2008 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_MESH_DISTRIBUTED_DATA_H
#define __DOLFIN_MESH_DISTRIBUTED_DATA_H

#include <dolfin/common/types.h>
#include <dolfin/common/Array.h>
#include <dolfin/common/DistributedData.h>

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
  MeshDistributedData(uint dim);

  /// Copy constructor
  MeshDistributedData(MeshDistributedData const& other);

  /// Destructor
  ~MeshDistributedData();

  /// Assignment
  MeshDistributedData& operator=(MeshDistributedData const& other);

  /// Equality
  bool operator==(MeshDistributedData const& other) const;

  /// Non-equality
  bool operator!=(MeshDistributedData const& other) const;

  /// Access to distributed data
  inline DistributedData& operator[](uint dim)
  {
    dolfin_assert(dim <= dim_);
    return data_[dim];
  }

  /// Access to distributed data (const)
  inline DistributedData const& operator[](uint dim) const
  {
    dolfin_assert(dim <= dim_);
    return data_[dim];
  }

  /// Return topological dimension of the distributed data
  uint dim() const;

  /// Display basic information
  void disp() const;

private:

  /// Clear
  void clear();

  // Topological dimensions
  uint dim_;

  // Distributed data
  DistributedData * data_;

};

} /* namespace dolfin */

#endif /* __DOLFIN_MESH_DISTRIBUTED_DATA_H */
