// Copyright (C) 2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-05-08
// Last changed: 2006-11-01

#include <dolfin/log/dolfin_log.h>
#include <dolfin/mesh/MeshConnectivity.h>
#include <dolfin/mesh/MeshTopology.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
MeshTopology::MeshTopology() :
    dim_(0),
    num_entities_(NULL),
    connectivity_(NULL),
    distdata_(*this),
    ordered_(false),
    timestamp_(std::time(NULL)),
    renumbering_count_(0)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
MeshTopology::MeshTopology(MeshTopology const& topology) :
    dim_(0),
    num_entities_(NULL),
    connectivity_(NULL),
    distdata_(*this),
    ordered_(false),
    timestamp_(std::time(NULL)),
    renumbering_count_(0)
{
  *this = topology;
}
//-----------------------------------------------------------------------------
MeshTopology::~MeshTopology()
{
  clear();
}
//-----------------------------------------------------------------------------
MeshTopology const& MeshTopology::operator=(MeshTopology const& topology)
{
  // Clear old data if any
  clear();

  // Allocate data
  dim_ = topology.dim_;
  num_entities_ = new uint[dim_ + 1];
  connectivity_ = new MeshConnectivity*[dim_ + 1];
  for (uint d = 0; d <= dim_; d++)
  {
    connectivity_[d] = new MeshConnectivity[dim_ + 1];
  }

  // Copy data
  if (dim_ > 0)
  {
    for (uint d = 0; d <= dim_; ++d)
    {
      num_entities_[d] = topology.num_entities_[d];
    }
    for (uint d0 = 0; d0 <= dim_; ++d0)
    {
      for (uint d1 = 0; d1 <= dim_; ++d1)
      {
        connectivity_[d0][d1] = topology.connectivity_[d0][d1];
      }
    }
  }
  distdata_ = topology.distdata_;
  ordered_ = topology.ordered_;
  timestamp_ = topology.timestamp_;
  renumbering_count_ = topology.renumbering_count_;

  return *this;
}
//-----------------------------------------------------------------------------
void MeshTopology::clear()
{
  // Clear parallel data structures
  distdata_.clear();

  // Delete number of mesh entities
  delete[] num_entities_;
  num_entities_ = NULL;

  // Delete mesh connectivity
  if (connectivity_)
  {
    for (uint d = 0; d <= dim_; ++d)
    {
      delete[] connectivity_[d];
    }
    delete[] connectivity_;
  }
  connectivity_ = NULL;
  ordered_ = false;
  timestamp_ = 0;
  renumbering_count_ = 0;

  // Reset dimension
  dim_ = 0;
}
//-----------------------------------------------------------------------------
void MeshTopology::init(uint dim)
{
  // Clear old data if any
  clear();

  timestamp_ = std::time(NULL); // Reset token

  // Initialize number of mesh entities
  num_entities_ = new uint[dim + 1];
  for (uint d = 0; d <= dim; d++)
  {
    num_entities_[d] = 0;
  }

  // Initialize mesh connectivity
  connectivity_ = new MeshConnectivity*[dim + 1];
  for (uint d = 0; d <= dim; d++)
  {
    connectivity_[d] = new MeshConnectivity[dim + 1];
  }

  // Save dimension
  dim_ = dim;
  distdata_.init(dim_);
}
//-----------------------------------------------------------------------------
void MeshTopology::init(uint dim, uint size)
{
  dolfin_assert(num_entities_);dolfin_assert(dim <= dim_);

  num_entities_[dim] = size;
}
//-----------------------------------------------------------------------------
uint MeshTopology::compute_entities(Mesh& mesh, uint dim) const
{
  return TopologyComputation::computeEntities(mesh, dim);
}
//-----------------------------------------------------------------------------
void MeshTopology::compute_connectivity(Mesh& mesh, uint d0, uint d1) const
{
  TopologyComputation::computeConnectivity(mesh, d0, d1);
}
//-----------------------------------------------------------------------------
int MeshTopology::token() const
{
  return timestamp_ + size(0) + size(dim_); // FIXME
}
//-----------------------------------------------------------------------------
void MeshTopology::disp() const
{
  // Begin indentation
  cout << "Mesh topology" << endl;
  begin("-------------");
  cout << endl;

  // Check if empty
  if (dim_ == 0)
  {
    cout << "empty" << endl << endl;
    end();
    return;
  }

  // Display topological dimension
  cout << "Topological dimension: " << dim_ << endl << endl;

  // Display number of entities for each topological dimension
  cout << "Number of entities:" << endl;
  begin("");
  for (uint d = 0; d <= dim_; ++d)
  {
    cout << "dim = " << d << ": " << num_entities_[d] << endl;
  }
  end();
  cout << endl;

  // Display matrix of connectivities
  cout << "Connectivity:" << endl;
  begin("");
  cout << " ";
  for (uint d1 = 0; d1 <= dim_; ++d1)
  {
    cout << " " << d1;
  }
  cout << endl;
  for (uint d0 = 0; d0 <= dim_; ++d0)
  {
    cout << d0;
    for (uint d1 = 0; d1 <= dim_; ++d1)
    {
      if (connectivity_[d0][d1].size() > 0)
      {
        cout << " x";
      }
      else
      {
        cout << " -";
      }
    }
    cout << endl;
  }
  cout << endl;
  end();

  // Display connectivity for each topological dimension
  for (uint d0 = 0; d0 <= dim_; ++d0)
  {
    for (uint d1 = 0; d1 <= dim_; ++d1)
    {
      if (connectivity_[d0][d1].size() == 0)
      {
        continue;
      }
      cout << "Connectivity " << d0 << " -- " << d1 << ":" << endl;
      begin("");
      connectivity_[d0][d1].disp();
      end();
      cout << endl;
    }
  }

  // End indentation
  end();
}
//-----------------------------------------------------------------------------

}


