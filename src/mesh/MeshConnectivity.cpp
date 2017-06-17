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
    order_(0),
    s_(0),
    min_degree_(0),
    max_degree_(0),
    offsets_(NULL),
    connections_(NULL)
{
}
//-----------------------------------------------------------------------------
MeshConnectivity::MeshConnectivity(uint order, uint degree, uint * graph) :
  is_initialized_(false),
  order_(0),
  s_(0),
  min_degree_(0),
  max_degree_(0),
  offsets_(NULL),
  connections_(NULL)
{
  init(order, degree, graph);
}
//-----------------------------------------------------------------------------
MeshConnectivity::MeshConnectivity(Array<uint> const& valency, uint * graph) :
  is_initialized_(false),
  order_(0),
  s_(0),
  min_degree_(0),
  max_degree_(0),
  offsets_(NULL),
  connections_(NULL)
{
  init(valency, graph);
}
//-----------------------------------------------------------------------------
MeshConnectivity::MeshConnectivity(MeshConnectivity const& other) :
    is_initialized_(false),
    order_(0),
    s_(0),
    min_degree_(0),
    max_degree_(0),
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
    order_ = other.order_;
    s_ = other.s_;
    if(other.offsets_ != NULL)
    {
      offsets_ = new uint[order_ + 1];
      for (uint e = 0; e <= order_; ++e)
      {
        offsets_[e] = other.offsets_[e];
      }
    }
    if (other.connections_ != NULL)
    {
      connections_ = new uint[s_];
      for (uint i = 0; i < s_; ++i)
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
  if (s_ != other.s_)
  {
    return false;
  }
  //
  if (order_ != other.order_)
  {
    return false;
  }
  //
  if (!objptrcmp(offsets_, other.offsets_))
  {
    return false;
  }
  //
  if (!cmp<uint>(s_, connections_, other.connections_))
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
void MeshConnectivity::init(uint order, uint degree)
{
  init(order, degree, NULL);
}
//-----------------------------------------------------------------------------
void MeshConnectivity::init(uint order, uint degree, uint * graph)
{
  clear();
  //
  is_initialized_= true;
  order_ = order;
  s_ = order * degree;
  min_degree_ = degree;
  max_degree_ = degree;
  offsets_ = new uint[order + 1];
  for (uint e = 0; e <= order; ++e)
  {
    offsets_[e] = e * degree;
  }
  if (graph != NULL)
  {
    connections_ = graph;
  }
  else if (s_ > 0)
  {
    connections_ = new uint[s_];
    for (uint i = 0; i < s_; ++i)
    {
      connections_[i] = 0;
    }
  }
}
//-----------------------------------------------------------------------------
void MeshConnectivity::init(Array<uint> const& valency)
{
  init(valency, NULL);
}
//-----------------------------------------------------------------------------
void MeshConnectivity::init(Array<uint> const& valency, uint * graph)
{
  clear();
  //
  is_initialized_= true;
  order_ = valency.size();
  s_ = 0;
  min_degree_ = 0;
  max_degree_ = 0;
  offsets_ = new uint[order_ + 1];
  offsets_[0] = 0;
  if (order_ > 0)
  {
    s_ = valency[0];
    min_degree_ = valency[0];
    max_degree_ = valency[0];
    for (uint e = 1; e < order_; ++e)
    {
      offsets_[e] = s_;
      s_ += valency[e];
      min_degree_ = std::min(min_degree_, valency[e]);
      max_degree_ = std::max(max_degree_, valency[e]);
    }
    offsets_[order_] = s_;

    if (graph != NULL)
    {
      connections_ = graph;
    }
    else
    {
      connections_ = new uint[s_];
      std::fill_n(connections_, s_, 0);
    }
  }
  else if (graph != NULL)
  {
    error("MeshConnectivity : assigning empty connectivity with non-zero pointer");
  }
}
//-----------------------------------------------------------------------------
void MeshConnectivity::clear()
{
  is_initialized_= false;
  order_ = 0;
  s_ = 0;
  min_degree_ = 0;
  max_degree_ = 0;
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
uint MeshConnectivity::order() const
{
  return order_;
}
//-----------------------------------------------------------------------------
uint MeshConnectivity::entries() const
{
  return s_;
}
//-----------------------------------------------------------------------------
uint MeshConnectivity::min_degree() const
{
  return min_degree_;
}
//-----------------------------------------------------------------------------
uint MeshConnectivity::max_degree() const
{
  return max_degree_;
}
//-----------------------------------------------------------------------------
void MeshConnectivity::set(uint entity, uint const * connections)
{
  dolfin_assert(entity < order_);
  dolfin_assert(connections_);
  for (uint i = 0; i < offsets_[entity + 1] - offsets_[entity]; ++i)
  {
    connections_[offsets_[entity] + i] = connections[i];
  }
}
//-----------------------------------------------------------------------------
void MeshConnectivity::set(Array<uint> const& connectivity)
{
  if (connectivity.size() != s_)
  {
    error("MeshConnectivity : provided connectivity size %u does no match %u",
          connectivity.size(), s_);
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
  order_ = connections.size();
  s_ = 0;
  min_degree_ = 0;
  max_degree_ = 0;
  offsets_ = new uint[order_ + 1];
  offsets_[0] = 0;
  if (order_ > 0)
  {
    s_ = connections[0].size();
    min_degree_ = s_;
    max_degree_ = s_;
    for (uint e = 1; e < order_; ++e)
    {
      offsets_[e] = s_;
      uint const s = connections[e].size();
      s_ += s;
      min_degree_ = std::min(min_degree_, s);
      max_degree_ = std::max(max_degree_, s);
    }
    offsets_[order_] = s_;
    connections_ = new uint[s_];
    for (uint e = 0; e < order_; ++e)
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
  if (order_ > 0)
  {
    if(map.size() != order_)
    {
      error("MeshConnectivity : remap_left mapping has invalid size");
    }
    // Set sizes in place of offsets to depict connectivities layout
    uint * ocpy = new uint[order_ + 1];
    std::swap(offsets_, ocpy);
    uint * ccpy = new uint[s_];
    std::swap(connections_, ccpy);
    offsets_[0] = 0;
    for (uint e = 0; e < order_; ++e)
    {
      dolfin_assert(map[e] < order_);
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
  for (uint i = 0; i < s_; ++i)
  {
    connections_[i] = map[connections_[i]];
  }
}
//-----------------------------------------------------------------------------
void MeshConnectivity::disp() const
{
  section("MeshConnectivity");
  if (s_ == 0)
  {
    cout << "empty" << endl;
  }
  else
  {
    // Display all connections
    for (uint e = 0; e < order_; ++e)
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
  A.assign(connections_, connections_ + s_);
  // Set stride if the graph is regular
  if(min_degree_ == max_degree_) A %= min_degree_;
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

  if(s_ < order_)
  {
    message("MeshConnectivity : connectivities size %u < %u number of entities",
            s_, order_);
  }
  //
  for(uint e0 = 0; e0 < this->order(); ++e0)
  {
    _set<uint> ce;
    for(uint e1 = 0; e1 <this->degree(e0); ++e1)
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

