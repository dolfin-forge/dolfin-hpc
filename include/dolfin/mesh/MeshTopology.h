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

#include <dolfin/log/log.h>
#include <dolfin/common/types.h>
#include <dolfin/common/Array.h>
#include <dolfin/main/MPI.h>
#include <dolfin/mesh/MeshConnectivity.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/MeshRenumber.h>

namespace dolfin
{

class Mesh;
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

class MeshTopology
{

public:

  /// Create empty mesh topology
  MeshTopology(Mesh& mesh);

  /// Copy constructor
  MeshTopology(MeshTopology const& other);

  /// Destructor
  ~MeshTopology();

  /// Assignment
  MeshTopology const& operator=(MeshTopology const& other);

  /// Equality
  bool operator==(MeshTopology const& other) const;

  /// Non-equality
  bool operator!=(MeshTopology const& other) const;

  /// Initialize topology of given maximum dimension and optionally set a flag
  /// to not distribute the topology
  void init(uint dim, bool distribute = true);

  /// Initialize topology entities for given maximum dimension
  /// Optionally specify the global number of entities for a distributed mesh.
  /// If the topology is not distributed, any value different than zero or the
  /// number of local entities will trigger an error.
  void init(uint dim, uint num_local, uint num_global = 0);

  /// Clear all data
  void clear();

  /// Finalize: check, reorder, renumber
  void finalize();

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

  /// Return if the topology is distributed
  bool is_distributed() const;

  /// Return mesh distribution data if the topology is distributed
  MeshDistributedData& distdata();

  /// Return mesh distribution data if the topology is distributed (const)
  MeshDistributedData const& distdata() const;

  /// Return number of entities in the global topology for given dimension
  uint global_size(uint dim) const;

  /// Return number of given entities
  uint num_owned(uint dim) const;

  /// Return number of given entities
  uint num_shared(uint dim) const;

  /// Return number of given entities
  uint num_ghost(uint dim) const;

  /// Return global index of mesh entity
  uint get_global(MeshEntity const& entity) const;

  /// Return local index of mesh entity: local-to-global then global-to-local
  /// Can be used for checking bijectivity of mappings
  uint get_local(MeshEntity const& entity) const;

  /// Return if the given mesh entity is shared
  bool is_owned(MeshEntity const& entity) const;

  /// Return if the given mesh entity is shared
  bool is_shared(MeshEntity const& entity) const;

  /// Return if the given mesh entity is ghosted
  bool is_ghost(MeshEntity const& entity) const;

  /// Return owner of the entity
  uint get_owner(MeshEntity const& entity) const;

  //---------------------------------------------------------------------------

  /// Display data
  void disp() const;

  //--- CHECK ROUTINES --------------------------------------------------------

  /// Check
  void check() const;

  //--- TOKENIZED -------------------------------------------------------------

  /// Return token identifying the internal state of mesh topology
  int token() const;

private:

  /// Update token value
  void update_token();

  //---------------------------------------------------------------------------

  /// Compute connectivity for given pair of topological dimensions
  void compute_connectivity(uint d0, uint d1) const;

  /// Force renumbering of mesh topology entities
  void renumber() const;

  /// Force reordering of mesh topology connectivities
  void reorder() const;

  ///
  Mesh& mesh_;

  /// Topological dimension
  uint dim_;

  /// Number of mesh vertices
  uint num_vertices_;

  /// Connectivity for pairs of topological dimensions
  MeshConnectivity ** connectivity_;

  /// Distributed mesh topology data
  MeshDistributedData * distdata_;

  //
  int timestamp_;

};

} /* namespace dolfin */

#endif /* __DOLFIN_MESH_TOPOLOGY_H */
