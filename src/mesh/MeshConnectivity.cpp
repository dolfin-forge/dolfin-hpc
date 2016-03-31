// Copyright (C) 2006-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-05-09
// Last changed: 2007-03-01

#include <dolfin/mesh/MeshConnectivity.h>

#include <dolfin/common/Array.h>
#include <dolfin/log/LogStream.h>

#include <cstring>

namespace dolfin
{

//-----------------------------------------------------------------------------
MeshConnectivity::MeshConnectivity() :
    is_initialized_(false),
    num_entities_(0),
    size_(0),
    connections_(NULL),
    offsets_(NULL)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
MeshConnectivity::MeshConnectivity(MeshConnectivity const& other) :
    is_initialized_(false),
    num_entities_(0),
    size_(0),
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
  is_initialized_ = other.is_initialized_;
  num_entities_ = other.num_entities_;
  size_ = other.size_;
  if(size_ > 0)
  {
    connections_ = new uint[size_];
  }
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
bool MeshConnectivity::operator==(MeshConnectivity const& other) const
{
  if (this == &other)
  {
    return true;
  }
  //
  if(is_initialized_ != other.is_initialized_)
  {
    return false;
  }
  //
  if (size_ != other.size_)
  {
    return false;
  }
  //
  if (num_entities_ != other.num_entities_)
  {
    return false;
  }
  //
  for (uint i = 0; i < size_; ++i)
  {
    if (connections_[i] != other.connections_[i])
    {
      return false;
    }
  }
  //
  for (uint e = 0; e <= num_entities_; ++e)
  {
    if (offsets_[e] != other.offsets_[e])
    {
      return false;
    }
  }
  //
  return true;
}
//-----------------------------------------------------------------------------
bool MeshConnectivity::operator!=(MeshConnectivity const& other) const
{
  return !(*this == other);
}
//-----------------------------------------------------------------------------
void MeshConnectivity::init(uint num_entities, uint num_connections)
{
  clear();
  //
  is_initialized_= true;
  num_entities_ = num_entities;
  size_ = num_entities * num_connections;
  offsets_ = new uint[num_entities + 1];
  for (uint e = 0; e <= num_entities; ++e)
  {
    offsets_[e] = e * num_connections;
  }
  if (size_ > 0)
  {
    connections_ = new uint[size_];
    for (uint i = 0; i < size_; ++i)
    {
      connections_[i] = 0;
    }
  }
}
//-----------------------------------------------------------------------------
void MeshConnectivity::init(Array<uint> const& num_connections)
{
  clear();
  //
  is_initialized_= true;
  num_entities_ = num_connections.size();
  size_ = 0;
  offsets_ = new uint[num_entities_ + 1];
  for (uint e = 0; e < num_entities_; ++e)
  {
    offsets_[e] = size_;
    size_ += num_connections[e];
  }
  offsets_[num_entities_] = size_;
  if (size_ > 0)
  {
    connections_ = new uint[size_];
    for (uint i = 0; i < size_; ++i)
    {
      connections_[i] = 0;
    }
  }
}
//-----------------------------------------------------------------------------
void MeshConnectivity::clear()
{
  is_initialized_= false;
  num_entities_ = 0;
  size_ = 0;
  delete[] connections_;
  connections_ = NULL;
  delete[] offsets_;
  offsets_ = NULL;
}
//-----------------------------------------------------------------------------
uint * MeshConnectivity::operator()()
{
  return connections_;
}
//-----------------------------------------------------------------------------
uint const * MeshConnectivity::operator()() const
{
  return connections_;
}
//-----------------------------------------------------------------------------
bool MeshConnectivity::is_initialized() const
{
  return is_initialized_;
}
//-----------------------------------------------------------------------------
uint MeshConnectivity::num_entities() const
{
  return num_entities_;
}
//-----------------------------------------------------------------------------
uint MeshConnectivity::size() const
{
  return size_;
}
//-----------------------------------------------------------------------------
void MeshConnectivity::set(uint entity, uint connection, uint pos)
{
  dolfin_assert(entity < num_entities_);
  dolfin_assert(pos < offsets_[entity + 1] - offsets_[entity]);

  connections_[offsets_[entity] + pos] = connection;
}
//-----------------------------------------------------------------------------
void MeshConnectivity::set(uint entity, Array<uint> const& connections)
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
void MeshConnectivity::set(Array<Array<uint> > const& connections)
{
  // Clear old data if any
  clear();

  // Initialize offsets and compute total size
  num_entities_ = connections.size();
  size_ = 0;
  offsets_ = new uint[num_entities_ + 1];
  for (uint e = 0; e < num_entities_; ++e)
  {
    offsets_[e] = size_;
    size_ += connections[e].size();
  }
  offsets_[num_entities_] = size_;

  // Initialize connections
  if(size_ > 0)
  {
    this->connections_ = new uint[size_];
    for (uint e = 0; e < num_entities_; ++e)
    {
      for (uint i = 0; i < connections[e].size(); ++i)
      {
        this->connections_[offsets_[e] + i] = connections[e][i];
      }
    }
  }
}

//-----------------------------------------------------------------------------
void MeshConnectivity::remap_left(Array<uint> const& map)
{
  message(1, "MeshConnectivity : remap left");
  if (num_entities_ > 0)
  {
    if(map.size() != num_entities_)
    {
      error("MeshConnectivity : remap_left mapping has invalid size");
    }
    // Set sizes in place of offsets to depict connectivities layout
    uint * onew = new uint[num_entities_];
    onew[0] = 0;
    for(uint e = 0; e < num_entities_; ++e)
    {
      dolfin_assert(map[e] < num_entities_);
      onew[e + 1] = offsets_[map[e] + 1] - offsets_[map[e]];
    }
    delete[] offsets_;
    offsets_ = onew;
    uint * cnew = new uint[size_];
    for (uint e = 0; e < num_entities_; ++e)
    {
      std::memcpy(&cnew[onew[e]], &connections_[offsets_[map[e]]],
                  offsets_[e + 1] * sizeof(uint));
      offsets_[e + 1] += offsets_[e];
    }
    delete[] connections_;
    connections_ = cnew;
  }
}
//-----------------------------------------------------------------------------
void MeshConnectivity::remap_right(Array<uint> const& map)
{
  message(1, "MeshConnectivity : remap right");
  for (uint i = 0; i < size_; ++i)
  {
    connections_[i] = map[connections_[i]];
  }
}
//-----------------------------------------------------------------------------
void MeshConnectivity::disp() const
{
  section("MeshConnectivity");
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
  end();
}
//-----------------------------------------------------------------------------
void MeshConnectivity::check() const
{
  message("MeshConnectivity : check");

  /**
   *  CHECK:
   *
   *  Check connectivity size, number of entities and verify that connected
   *  entities are not listed twice.
   *
   */

  if(size_ < num_entities_)
  {
    message("MeshConnectivity : connectivities size %u < %u number of entities",
            size_, num_entities_);
  }
  //
  for(uint e0 = 0; e0 < this->num_entities(); ++e0)
  {
    _set<uint> ce;
    for(uint e1 = 0; e1 <this->size(e0); ++e1)
    {
      uint ec = (*this)(e0)[e1];
      if(ce.count(ec) > 0)
      {
        error("Entity %u appears twice in connectivities for %u", ec, e0);
      }
      ce.insert(ec);
    }
  }
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */

