// Copyright (C) 2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
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

  /// Initialize topology of given maximum dimension
  void init(uint dim, uint num_vertices);

  /// Clear all data
  void clear();

  /// Finalize
  void finalize();

  /// Renumber mesh topology entities
  void renumber();

  /// Remap local entities of given dimension
  void remap(uint dim,  Array<uint> const& map);

  /// Return connectivity for given pair of topological dimensions
  MeshConnectivity& operator()(uint d0, uint d1);

  /// Return connectivity for given pair of topological dimensions
  MeshConnectivity const& operator()(uint d0, uint d1) const;

  /// Return mesh distribution data
  MeshDistributedData& distdata();

  /// Return mesh distribution data (const version)
  MeshDistributedData const& distdata() const;

  //---------------------------------------------------------------------------

  /// Return topological dimension
  uint dim() const;

  /// Return number of entities in the local topology for given dimension
  uint size(uint dim) const;

  /// Return number of entities in the global topology for given dimension
  uint global_size(uint dim) const;

  /// Return number of given entities
  uint num_owned(uint dim) const;

  /// Return number of given entities
  uint num_shared(uint dim) const;

  /// Return number of given entities
  uint num_ghosts(uint dim) const;

  /// Return if connectivity for given pair is computed
  bool is_computed(uint d0, uint d1) const;

  /// Return if entities exist
  bool entities_exist(uint d) const;

  ///
  bool is_distributed() const;

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

  ///
  void order();

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

  /// Return true iff topology is ordered according to the UFC numbering
  mutable bool is_ordered_;

  //
  int timestamp_;

  //
  uint renumbering_count_;

};

} /* namespace dolfin */

#endif /* __DOLFIN_MESH_TOPOLOGY_H */
