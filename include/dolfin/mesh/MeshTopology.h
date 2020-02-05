// Copyright (C) 2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

// This class has been extensively tortured since 2015 due to flawed design and
// non-robust implementation.
// Additionally the topology computation was rewritten to respect abstraction
// layers and avoid storing cell - cell neighbours connectivities.
// The updated interface gives access to the distributed data to allow writing
// generic serial/parallel code with no requirement of parallel guards for
// sections reserved to distributed meshes.
// Additionally it benefits from several inconsistency fixes merged from mesh
// distributed data.

#ifndef __DOLFIN_MESH_TOPOLOGY_H
#define __DOLFIN_MESH_TOPOLOGY_H

#include <dolfin/common/Tokenized.h>
#include <dolfin/common/Clonable.h>

#include <dolfin/log/log.h>
#include <dolfin/common/types.h>
#include <dolfin/common/Array.h>
#include <dolfin/common/Distributed.h>
#include <dolfin/mesh/CellType.h>
#include <dolfin/mesh/Connectivity.h>
#include <dolfin/mesh/MeshRenumber.h>
#include <dolfin/mesh/MeshDistributedData.h>

namespace dolfin
{

class Mesh;
class Connectivity;

/**
 *
 *  @class  MeshTopology
 *
 *  @brief  MeshTopology stores the topology of a mesh, consisting of mesh
 *          entities and connectivity (incidence relations of mesh entities).
 *          Note that the mesh entities do not need to be stored, only the
 *          number of entities and the connectivity. Any numbering scheme for
 *          the mesh entities is stored separately in a MeshFunction over the
 *          entities.
 *          A mesh entity e may be identified globally as a pair e = (dim, i),
 *          where dim is the topological dimension and i is the index of the
 *          entity within that topological dimension.
 *
 */

class MeshTopology: public Clonable<MeshTopology>, public Distributed<MeshTopology>
{
  // Set max connectivity matrix dimension
  static uint const CMAX = 4;

  // Save some limbo at MeshEntity construction until classes are rewritten
  friend class MeshEntity;

public:

  /// Create mesh topology for given cell type
  MeshTopology(CellType const& type, Comm& comm, bool frozen);

  /// Copy constructor
  MeshTopology(MeshTopology const& other);

  /// Destructor
  ~MeshTopology();

  /// Swap instances
  friend void swap( MeshTopology& a, MeshTopology& b );

  /// Assignment (Disabled)
  MeshTopology & operator=( const MeshTopology & other );

  /// Equality
  bool operator==(MeshTopology const& other) const;

  /// Non-equality
  bool operator!=(MeshTopology const& other) const;

  /// Set topology entities for given dimension
  /// Optionally specify the global number of entities for a distributed mesh.
  /// If the topology is not distributed, any value different than zero or the
  /// number of local entities will trigger an error.
  void init(uint dim, uint nlocal, uint nglobal = 0);

  /// Finalize: check, reorder, renumber
  void finalize();

  ///
  CellType const& type() const;

  ///
  CellType const& type(uint i) const;

  /// Remap local entities of given dimension
  void remap(uint d0,  Array<uint> const& mapping);

  //--- Connectivity ----------------------------------------------------------

  /// Return connectivity for given pair of topological dimensions
  Connectivity& operator()(uint d0, uint d1);

  /// Return connectivity for given pair of topological dimensions
  Connectivity const& operator()(uint d0, uint d1) const;

  /// Return topological dimension
  uint dim() const;

  /// Return number of entities in the local topology for given dimension
  uint size(uint dim) const;

  /// Return pointer to connectivity for given pair
  Connectivity * connectivity(uint d0, uint d1 = 0);

  /// Return pointer to connectivity for given pair (const)
  Connectivity const * connectivity(uint d0, uint d1 = 0) const;

  //--- Distributed data ------------------------------------------------------

  /// Return mesh distribution data if the topology is distributed
  MeshDistributedData& distdata();

  /// Return mesh distribution data if the topology is distributed (const)
  MeshDistributedData const& distdata() const;

  /// Return number of entities in the global topology for given dimension
  uint global_size(uint dim) const;

  /// Return offset of global indices on current rank
  uint offset(uint dim) const;

  /// Return number of given entities
  uint num_owned(uint dim) const;

  /// Return number of given entities
  uint num_shared(uint dim) const;

  /// Return number of given entities
  uint num_ghost(uint dim) const;

  //---------------------------------------------------------------------------

  /// Display data
  void disp() const;

  //--- TOKENIZED -------------------------------------------------------------

  /// Return token identifying the internal state of mesh topology
  int token() const;

private:

  /// Update token value
  void update_token();

  //---------------------------------------------------------------------------

  /// Compute connectivity for given pair of topological dimensions
  Connectivity const * compute(uint d0, uint d1) const;

  /// Compute entities for given topological dimension
  Connectivity const * entities(uint di) const;

  /// Compute transpose for given pair of topological dimensions
  Connectivity const * transpose(uint d0, uint d1) const;

  /// Compute connectivity for given triple of topological dimensions
  Connectivity const * intersection(uint d0, uint di, uint d1) const;

public:

  /// Force renumbering of mesh topology entities
  /// @todo public for the moment but this just legacy of bad design
  void renumber() const;

private:
  /// Force reordering of mesh topology connectivities
  void reorder() const;

  ///
  CellType const * type_;

  /// Topological dimension
  uint dim_;

  // Topology cannot be modified
  bool frozen_;

  /// Connectivity for pairs of topological dimensions
  mutable Connectivity * C_[CMAX][CMAX];

  /// Distributed mesh topology data
  MeshDistributedData distdata_;

  //
  int timestamp_;

};

} /* namespace dolfin */

#endif /* __DOLFIN_MESH_TOPOLOGY_H */
