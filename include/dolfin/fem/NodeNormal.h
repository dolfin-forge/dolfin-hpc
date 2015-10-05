// Copyright (C) 2007 Murtazo Nazarov
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2009.
// Modified by Aurélien Larcher, 2012-14. (rewrite, extension to any element)
//
// First added:  2007-05-01
// Last changed: 2014-05-22

#ifndef __DOLFIN_NODENORMAL_H
#define __DOLFIN_NODENORMAL_H

#include <dolfin/fem/BoundaryNormal.h>
#include <dolfin/mesh/Mesh.h>

#include <dolfin/common/constants.h>
#include <dolfin/common/Array.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/BoundaryMesh.h>

#include <map>


namespace dolfin
{

class SubDomain;

class NodeNormal : public BoundaryNormal
{

public:

  enum Type
  {
    none, unit, facet
  };

  /// Create normal, tangents for the boundary of mesh
  NodeNormal(Mesh& mesh, Type w = none, real alpha = 1.57);

  /// Create normal, tangents for the boundary of mesh for given subdomain
  NodeNormal(Mesh& mesh, SubDomain const& subdomain, Type w = none,
             real alpha = 1.57);

  /// Destructor
  ~NodeNormal();

  /// Compute the orthogonal basis
  void compute();

  /// Returns the node type
  uint node_type(uint node_id) const;

  /// Assignment
  NodeNormal& operator=(NodeNormal& node_normal);

private:

  /// Cleanup
  void Clear();

  /// Compute boundary normal basis
  void ComputeP1(Mesh& mesh, Array<Function>& basis);
  void ComputePk(Mesh& mesh, Array<Function>& basis);

  //--- ATTRIBUTES ------------------------------------------------------------

  Mesh& mesh_;

  SubDomain const * const subdomain_;
  bool const no_subdomain_;

  /// Maximum absolute angle between two neighbouring facets to be discriminated
  /// as belonging to different hyperplanes.
  real const alpha_max_;

  /// Type of weight used for computing the node normal from facet normals.
  Type const type_;

  ///
  _map<uint, uint> node_type_;

  struct FacetData
  {
    uint global_index;
    real weight;
    Point normal;
    _set<uint> nodes;

    void disp() const
    {
      cout << "FacetData" << endl;
      cout << "---------" << endl;
      // Begin indentation
      begin("");
      cout << "global_index : " << global_index << endl;
      cout << "nodes        : " << (uint) nodes.size() << endl;
      cout << "weight       : " << weight << endl;
      cout << "normal       : " << normal.x() << ", "
      << normal.y() << ", "
      << normal.z() << endl;
      // End indentation
      end();
      cout << endl;
    }
  };

  struct NodeData
  {
    uint node_type;
    Array<uint> dofs;
    Array<uint> adjs;
    Array<FacetData *> facets;

    void disp() const
    {
      cout << "NodeData" << endl;
      cout << "--------" << endl;
      // Begin indentation
      begin("");
      cout << "node_type    : " << node_type << endl;
      cout << "dofs         : " << (uint) dofs.size() << endl;
      cout << "adjs         : " << (uint) adjs.size() << endl;
      cout << "facets       : " << (uint) facets.size() << endl;
      // End indentation
      end();
      cout << endl;
    }
  };

};

}
#endif

