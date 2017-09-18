// Copyright (C) 2006-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-05-09
// Last changed: 2007-03-01

#include <dolfin/mesh/MeshConnectivity.h>

#include <dolfin/common/Array.h>
#include <dolfin/log/LogStream.h>

#include <algorithm>
#include <cstring>

namespace dolfin
{

//-----------------------------------------------------------------------------
MeshConnectivity::MeshConnectivity() :
    is_initialized_(false),
    num_entities_(0),
    size_(0),
    min_connections_(0),
    max_connections_(0),
    offsets_(NULL),
    connections_(NULL)
{
}
//-----------------------------------------------------------------------------
MeshConnectivity::MeshConnectivity(MeshConnectivity const& other) :
    is_initialized_(false),
    num_entities_(0),
    size_(0),
    min_connections_(0),
    max_connections_(0),
    offsets_(NULL),
    connections_(NULL)
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
  if(&other != this)
  {
    clear();
    //
    is_initialized_ = other.is_initialized_;
    num_entities_ = other.num_entities_;
    size_ = other.size_;
    if(other.offsets_ != NULL)
    {
      offsets_ = new uint[num_entities_ + 1];
      for (uint e = 0; e <= num_entities_; ++e)
      {
        offsets_[e] = other.offsets_[e];
      }
    }
    if (other.connections_ != NULL)
    {
      connections_ = new uint[size_];
      for (uint i = 0; i < size_; ++i)
      {
        connections_[i] = other.connections_[i];
      }
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
  if (!objptrcmp(offsets_, other.offsets_))
  {
    return false;
  }
  //
  if (!cmp<uint>(size_, connections_, other.connections_))
  {
    return false;
  }
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
  init(NULL, num_entities, num_connections);
}
//-----------------------------------------------------------------------------
void MeshConnectivity::init(uint * connectivity, uint num_entities,
                            uint num_connections)
{
  clear();
  //
  is_initialized_= true;
  num_entities_ = num_entities;
  size_ = num_entities * num_connections;
  min_connections_ = num_connections;
  max_connections_ = num_connections;
  offsets_ = new uint[num_entities + 1];
  for (uint e = 0; e <= num_entities; ++e)
  {
    offsets_[e] = e * num_connections;
  }
  if (connectivity != NULL)
  {
    connections_ = connectivity;
  }
  else if (size_ > 0)
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
  init(NULL, num_connections);
}
//-----------------------------------------------------------------------------
void MeshConnectivity::init(uint * connectivity, Array<uint> const& num_connections)
{
  clear();
  //
  is_initialized_= true;
  num_entities_ = num_connections.size();
  size_ = 0;
  min_connections_ = 0;
  max_connections_ = 0;
  offsets_ = new uint[num_entities_ + 1];
  offsets_[0] = 0;
  if (num_entities_ > 0)
  {
    size_ = num_connections[0];
    min_connections_ = num_connections[0];
    max_connections_ = num_connections[0];
    for (uint e = 1; e < num_entities_; ++e)
    {
      offsets_[e] = size_;
      size_ += num_connections[e];
      min_connections_ = std::min(min_connections_, num_connections[e]);
      max_connections_ = std::max(max_connections_, num_connections[e]);
    }
    offsets_[num_entities_] = size_;

    if (connectivity != NULL)
    {
      connections_ = connectivity;
    }
    else
    {
      connections_ = new uint[size_];
      std::fill_n(connections_, size_, 0);
    }
  }
  else if (connectivity != NULL)
  {
    error("MeshConnectivity : assigning empty connectivity with non-zero pointer");
  }
}
//-----------------------------------------------------------------------------
void MeshConnectivity::clear()
{
  is_initialized_= false;
  num_entities_ = 0;
  size_ = 0;
  min_connections_ = 0;
  max_connections_ = 0;
  delete[] offsets_;
  offsets_ = NULL;
  delete[] connections_;
  connections_ = NULL;
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
uint MeshConnectivity::min_connections() const
{
  return min_connections_;
}
//-----------------------------------------------------------------------------
uint MeshConnectivity::max_connections() const
{
  return max_connections_;
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
  for (uint i = 0; i < connections.size(); ++i)
  {
    connections_[offsets_[entity] + i] = connections[i];
  }
}
//-----------------------------------------------------------------------------
void MeshConnectivity::set(uint entity, uint const * connections)
{
  dolfin_assert(entity < num_entities_);
  dolfin_assert(connections_);
  for (uint i = 0; i < offsets_[entity + 1] - offsets_[entity]; ++i)
  {
    connections_[offsets_[entity] + i] = connections[i];
  }
}
//-----------------------------------------------------------------------------
void MeshConnectivity::set(Array<uint> const& connectivity)
{
  if (connectivity.size() != size_)
  {
    error("MeshConnectivity : provided connectivity size %u does no match %u",
          connectivity.size(), size_);
  }
  if (connections_ == NULL)
  {
    error("MeshConnectivity : connectivity is not initialized");
  }
  std::copy(connectivity.begin(), connectivity.end(), connections_);
}
//-----------------------------------------------------------------------------
void MeshConnectivity::set(Array<Array<uint> > const& connections)
{
  clear();
  //
  is_initialized_= true;
  num_entities_ = connections.size();
  size_ = 0;
  min_connections_ = 0;
  max_connections_ = 0;
  offsets_ = new uint[num_entities_ + 1];
  offsets_[0] = 0;
  if (num_entities_ > 0)
  {
    size_ = connections[0].size();
    min_connections_ = size_;
    max_connections_ = size_;
    for (uint e = 1; e < num_entities_; ++e)
    {
      offsets_[e] = size_;
      uint const s = connections[e].size();
      size_ += s;
      min_connections_ = std::min(min_connections_, s);
      max_connections_ = std::max(max_connections_, s);
    }
    offsets_[num_entities_] = size_;
    connections_ = new uint[size_];
    for (uint e = 0; e < num_entities_; ++e)
    {
      for (uint i = 0; i < connections[e].size(); ++i)
      {
        connections_[offsets_[e] + i] = connections[e][i];
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
    uint * ocpy = new uint[num_entities_ + 1];
    std::swap(offsets_, ocpy);
    uint * ccpy = new uint[size_];
    std::swap(connections_, ccpy);
    offsets_[0] = 0;
    for (uint e = 0; e < num_entities_; ++e)
    {
      dolfin_assert(map[e] < num_entities_);
      uint const ii = map[e];
      offsets_[e + 1] = offsets_[e] + ocpy[ii + 1] - ocpy[ii];
      std::copy(ccpy + ocpy[ii], ccpy + ocpy[ii + 1], connections_+ offsets_[e]);
    }
    delete[] ocpy;
    delete[] ccpy;
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
  endblock();
}
//-----------------------------------------------------------------------------
MeshConnectivity const& MeshConnectivity::operator>>(Array<uint>& A) const
{
  A.assign(connections_, connections_ + size_);
  // Set stride if the graph is regular
  if(min_connections_ == max_connections_) A %= min_connections_;
  return *this;
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

