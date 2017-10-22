// Copyright (C) 2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Aurelien Larcher, 2015-2016.
//
// This class has been extensively tortured since 2015 due to flawed design and
// non-robust implementation.
// Additionally the topology computation was rewritten to respect abstraction
// layers and avoid storing cell - cell neighbours connectivities.
// The updated interface gives access to the distributed data to allow writing
// generic serial/parallel code with no requirement of parallel guards for
// sections reserved to distributed meshes.
// Additionally it benefits from several inconsistency fixes merged from mesh
// distributed data.
//
// First added:  2006-05-08
// Last changed: 2007-11-30

#ifndef __DOLFIN_MESH_TOPOLOGY_H
#define __DOLFIN_MESH_TOPOLOGY_H

#include <dolfin/common/Tokenized.h>
#include <dolfin/common/Clonable.h>

#include <dolfin/log/log.h>
#include <dolfin/common/types.h>
#include <dolfin/common/Array.h>
#include <dolfin/mesh/CellType.h>
#include <dolfin/mesh/MeshConnectivity.h>
#include <dolfin/mesh/MeshRenumber.h>

namespace dolfin
{

class MeshConnectivity;
class MeshDistributedData;

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

class MeshTopology: public Clonable<MeshTopology>
{
  // Save some limbo at MeshEntity construction until classes are rewritten
  friend class MeshEntity;

public:

  /// Create mesh topology
  MeshTopology();

  /// Create mesh topology for given cell type
  MeshTopology(CellType const& type, bool frozen = false);

  /// Copy constructor
  MeshTopology(MeshTopology const& other);

  /// Destructor
  ~MeshTopology();

  /// Assignment
  MeshTopology const& operator=(MeshTopology const& other);

  /// Swap instances
  void swap(MeshTopology& other);

  /// Equality
  bool operator==(MeshTopology const& other) const;

  /// Non-equality
  bool operator!=(MeshTopology const& other) const;

  /// Initialize topology for given cell type
  void init(CellType const& type, bool frozen = false);

  /// Set topology entities for given dimension
  /// Optionally specify the global number of entities for a distributed mesh.
  /// If the topology is not distributed, any value different than zero or the
  /// number of local entities will trigger an error.
  void init(uint dim, uint nlocal, uint nglobal = 0);

  /// Clear all data
  void clear();

  /// Finalize: check, reorder, renumber
  void finalize();

  ///
  CellType const& type(uint i) const;

  /// Remap local entities of given dimension
  void remap(uint dim,  Array<uint> const& mapping);

  //--- Connectivity ----------------------------------------------------------

  /// Return connectivity for given pair of topological dimensions
  MeshConnectivity& operator()(uint d0, uint d1);

  /// Return connectivity for given pair of topological dimensions
  MeshConnectivity const& operator()(uint d0, uint d1) const;

  /// Return topological dimension
  uint dim() const;

  /// Return number of entities in the local topology for given dimension
  uint size(uint dim) const;

  /// Return if connectivity for given pair is computed
  bool is_computed(uint d0, uint d1) const;

  /// Return if entities exist
  bool entities_exist(uint dim) const;

  //--- Distributed data ------------------------------------------------------

  /// Set the topology as distributed
  void set_distributed();

  /// Return if the topology is distributed
  bool is_distributed() const;

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
  void compute_connectivity(uint d0, uint d1) const;

public:

  /// Force renumbering of mesh topology entities
  ///FIXME: public for the moment but this just legacy of bad design
  void renumber() const;

private:

  /// Force reordering of mesh topology connectivities
  void reorder() const;

  ///
  CellType const * type_;

  /// Topological dimension
  uint dim_;

  /// Number of mesh vertices
  uint num_vertices_;
  bool ini_vertices_;

  /// Connectivity for pairs of topological dimensions
  MeshConnectivity ** connectivity_;

  /// Distributed mesh topology data
  MeshDistributedData * distdata_;

  //
  bool frozen_;

  //
  int timestamp_;

};

} /* namespace dolfin */

#endif /* __DOLFIN_MESH_TOPOLOGY_H */
