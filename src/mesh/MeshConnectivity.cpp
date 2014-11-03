// Copyright (C) 2006-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-05-09
// Last changed: 2007-03-01

#include <dolfin/log/dolfin_log.h>
#include <dolfin/mesh/MeshConnectivity.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
MeshConnectivity::MeshConnectivity() :
    size_(0),
    num_entities_(0),
    connections_(NULL),
    offsets_(NULL)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
MeshConnectivity::MeshConnectivity(MeshConnectivity const& other) :
    size_(0),
    num_entities_(0),
    connections_(NULL),
    offsets_(NULL)
{
  *this = other;
}
//-----------------------------------------------------------------------------
MeshConnectivity::~MeshConnectivity()
{
  clear();
}
//-----------------------------------------------------------------------------
MeshConnectivity const& MeshConnectivity::operator=(
    MeshConnectivity const& other)
{
  // Clear old data if any
  clear();

  // Allocate data
  size_ = other.size_;
  num_entities_ = other.num_entities_;
  connections_ = new uint[size_];
  offsets_ = new uint[num_entities_ + 1];

  // Copy data
  for (uint i = 0; i < size_; ++i)
  {
    connections_[i] = other.connections_[i];
  }
  if (num_entities_ > 0)
  {
    for (uint e = 0; e <= num_entities_; ++e)
    {
      offsets_[e] = other.offsets_[e];
    }
  }

  return *this;
}
//-----------------------------------------------------------------------------
void MeshConnectivity::clear()
{
  size_ = 0;
  num_entities_ = 0;
  delete[] connections_;
  connections_ = NULL;
  delete[] offsets_;
  offsets_ = NULL;
}
//-----------------------------------------------------------------------------
void MeshConnectivity::init(uint num_entities, uint num_connections)
{
  // Clear old data if any
  clear();

  // Compute the total size
  size_ = num_entities * num_connections;
  this->num_entities_ = num_entities;

  // Allocate data
  connections_ = new uint[size_];
  offsets_ = new uint[num_entities + 1];

  // Initialize data
  for (uint i = 0; i < size_; ++i)
  {
    connections_[i] = 0;
  }
  for (uint e = 0; e <= num_entities; ++e)
  {
    offsets_[e] = e * num_connections;
  }
}
//-----------------------------------------------------------------------------
void MeshConnectivity::init(Array<uint> const& num_connections)
{
  // Clear old data if any
  clear();

  // Initialize offsets and compute total size
  num_entities_ = num_connections.size();
  offsets_ = new uint[num_entities_ + 1];
  size_ = 0;
  for (uint e = 0; e < num_entities_; ++e)
  {
    offsets_[e] = size_;
    size_ += num_connections[e];
  }
  offsets_[num_entities_] = size_;

  // Initialize connections
  connections_ = new uint[size_];
  for (uint i = 0; i < size_; ++i)
  {
    connections_[i] = 0;
  }
}
//-----------------------------------------------------------------------------
void MeshConnectivity::set(uint entity, uint connection, uint pos)
{
  dolfin_assert(entity < num_entities_);
  dolfin_assert(pos < offsets_[entity + 1] - offsets_[entity]);

  connections_[offsets_[entity] + pos] = connection;
}
//-----------------------------------------------------------------------------
void MeshConnectivity::set(uint entity, const Array<uint>& connections)
{
  dolfin_assert(entity < num_entities_);
  dolfin_assert(connections.size() == offsets_[entity + 1] - offsets_[entity]);

  // Copy data
  for (uint i = 0; i < connections.size(); ++i)
  {
    this->connections_[offsets_[entity] + i] = connections[i];
  }
}
//-----------------------------------------------------------------------------
void MeshConnectivity::set(uint entity, uint const * connections)
{
  dolfin_assert(entity < num_entities_);
  dolfin_assert(connections_);

  // Copy data
  uint const num_connections = offsets_[entity + 1] - offsets_[entity];
  for (uint i = 0; i < num_connections; ++i)
  {
    this->connections_[offsets_[entity] + i] = connections[i];
  }
}
//-----------------------------------------------------------------------------
void MeshConnectivity::set(const Array<Array<uint> >& connections)
{
  // Clear old data if any
  clear();

  // Initialize offsets and compute total size
  num_entities_ = connections.size();
  offsets_ = new uint[num_entities_ + 1];
  size_ = 0;
  for (uint e = 0; e < num_entities_; ++e)
  {
    offsets_[e] = size_;
    size_ += connections[e].size();
  }
  offsets_[num_entities_] = size_;

  // Initialize connections
  this->connections_ = new uint[size_];
  for (uint e = 0; e < num_entities_; ++e)
  {
    for (uint i = 0; i < connections[e].size(); ++i)
    {
      this->connections_[offsets_[e] + i] = connections[e][i];
    }
  }
}
//-----------------------------------------------------------------------------
void MeshConnectivity::disp() const
{
  // Begin indentation
  cout << "MeshConnectivity" << endl;
  begin("----------------");
  cout << endl;

  // Check if there are any connections
  if (size_ == 0)
  {
    cout << "empty" << endl;
  }
  else
  {
    // Display all connections
    for (uint e = 0; e < num_entities_; ++e)
    {
      cout << e << ":";
      for (uint i = offsets_[e]; i < offsets_[e + 1]; ++i)
      {
        cout << " " << connections_[i];
      }
      cout << endl;
    }
  }
  // End indentation
  end();
}
//-----------------------------------------------------------------------------

}

