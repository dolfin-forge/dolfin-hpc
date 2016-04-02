// Copyright (C) 2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-05-08
// Last changed: 2014-11-03

#include <dolfin/mesh/MeshTopology.h>

#include <dolfin/log/log.h>
#include <dolfin/math/basic.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/EntityKey.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshConnectivity.h>

#include <ctime>

namespace dolfin
{

//-----------------------------------------------------------------------------
MeshTopology::MeshTopology(Mesh& mesh) :
    mesh_(mesh),
    dim_(0),
    num_vertices_(0),
    connectivity_(NULL),
    distdata_(NULL),
    is_ordered_(false),
    is_numbered_(false),
    timestamp_(0)
{
}
//-----------------------------------------------------------------------------
MeshTopology::MeshTopology(MeshTopology const& other) :
    mesh_(other.mesh_),
    dim_(0),
    num_vertices_(0),
    connectivity_(NULL),
    distdata_(NULL),
    is_ordered_(false),
    is_numbered_(false),
    timestamp_(0)
{
  *this = other;
}
//-----------------------------------------------------------------------------
MeshTopology::~MeshTopology()
{
  clear();
}
//-----------------------------------------------------------------------------
MeshTopology const& MeshTopology::operator=(MeshTopology const& other)
{
  clear();

  dim_ = other.dim_;
  num_vertices_ = other.num_vertices_;
  connectivity_ = new MeshConnectivity*[dim_ + 1];
  for (uint d0 = 0; d0 <= dim_; ++d0)
  {
    connectivity_[d0] = new MeshConnectivity[dim_ + 1];
    for (uint d1 = 0; d1 <= dim_; ++d1)
    {
      connectivity_[d0][d1] = other.connectivity_[d0][d1];
    }
  }
  if (other.distdata_ != NULL)
  {
    distdata_ = new MeshDistributedData(*other.distdata_);
  }
  is_ordered_ = other.is_ordered_;
  is_numbered_ = other.is_numbered_;
  timestamp_ = other.timestamp_;

  return *this;
}
//-----------------------------------------------------------------------------
bool MeshTopology::operator==(MeshTopology const& other) const
{
  if (this == &other)
  {
    return true;
  }
  //
  if (dim_ != other.dim_)
  {
    return false;
  }
  //
  for (uint i = 0; i <= dim_; ++i)
  {
    if (this->size(i) != other.size(i))
    {
      return false;
    }
  }
  //
  if (connectivity_)
  {
    for (uint d0 = 0; d0 <= dim_; ++d0)
    {
      for (uint d1 = 0; d1 <= dim_; ++d1)
      {
        if (connectivity_[d0][d1] != other.connectivity_[d0][d1])
        {
          return false;
        }
      }
    }
  }
  //
  if (!objptrcmp(distdata_, other.distdata_))
  {
    return false;
  }
  //
  return true;
}
//-----------------------------------------------------------------------------
bool MeshTopology::operator!=(MeshTopology const& other) const
{
  return !(*this == other);
}
//-----------------------------------------------------------------------------
void MeshTopology::init(uint dim, uint num_vertices)
{
  if (connectivity_ != NULL)
  {
    error("MeshTopology : clear instance before reinitializing");
  }
  if (dim == 0 && num_vertices > 1)
  {
    error("MeshTopology : more than one vertex for a mesh of dimension zero");
  }

  dim_ = dim;
  num_vertices_ = num_vertices;
  connectivity_ = new MeshConnectivity*[dim + 1];
  for (uint d = 0; d <= dim; ++d)
  {
    connectivity_[d] = new MeshConnectivity[dim + 1];
  }
  if (MPI::numProcesses() > 1)
  {
    distdata_ = new MeshDistributedData(dim_);
  }

  //
  update_token();
}
//-----------------------------------------------------------------------------
void MeshTopology::clear()
{
  // Clear parallel data structures
  delete distdata_;
  distdata_ = NULL;

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
  is_ordered_ = false;
  is_numbered_ = false;
  timestamp_ = 0;
  num_vertices_ = 0;
  dim_ = 0;
}
//-----------------------------------------------------------------------------
void MeshTopology::finalize()
{
  if(!connectivity_[dim_][0].is_initialized())
  {
    error("MeshTopology : cell -> vertices connectivity does not exist");
  }
  if(!is_ordered_)
  {
    reorder();
  }
  if(!is_numbered_)
  {
    renumber();
  }
}
//-----------------------------------------------------------------------------
void MeshTopology::reorder() const
{
  message(1, "MeshTopology : order");
  CellType const& cell_type = mesh_.type();
  for (CellIterator cell(mesh_); !cell.end(); ++cell)
  {
    cell_type.order_entities(*cell);
  }
  is_ordered_ = true;
}
//-----------------------------------------------------------------------------
void MeshTopology::renumber() const
{
  if (distdata_ == NULL)
  {
    return;
  }

  message(1, "MeshTopology : renumber");
  if(!MeshRenumber::renumber(mesh_))
  {
    error("Triggered mesh renumbering for nothing");
  }
  is_numbered_ = true;
}
//-----------------------------------------------------------------------------
void MeshTopology::remap(uint dim, Array<uint> const& map)
{
  if (connectivity_)
  {
    uint d0 = dim;
    for (uint d1 = 0; d1 <= dim_; ++d1)
    {
      if (connectivity_[d0][d1].size() > 0)
      {
        connectivity_[d0][d1].remap_left(map);
        update_token();
        is_ordered_ = false;
      }
      if (connectivity_[d1][d0].size() > 0)
      {
        connectivity_[d1][d0].remap_right(map);
        update_token();
        is_ordered_ = false;
      }
    }

    if (distdata_ != NULL)
    {
      error("MeshTopology : remapping entities in parallel is not allowed.");
    }
  }
}
//-----------------------------------------------------------------------------
MeshConnectivity& MeshTopology::operator()(uint d0, uint d1)
{
  dolfin_assert(d0 <= dim_ && d1 <= dim_);
  if(!connectivity_[d0][d1].is_initialized())
  {
    compute_connectivity(d0, d1);
  }
  return connectivity_[d0][d1];
}
//-----------------------------------------------------------------------------
MeshConnectivity const& MeshTopology::operator()(uint d0, uint d1) const
{
  dolfin_assert(d0 <= dim_ && d1 <= dim_);
  if(!connectivity_[d0][d1].is_initialized())
  {
    compute_connectivity(d0, d1);
  }
  return connectivity_[d0][d1];
}
//-----------------------------------------------------------------------------
MeshDistributedData& MeshTopology::distdata()
{
  if (distdata_ == NULL)
  {
    error("MeshDistributedData : returning distributed data of serial mesh");
  }
  return *distdata_;
}
//-----------------------------------------------------------------------------
MeshDistributedData const& MeshTopology::distdata() const
{
  if (distdata_ == NULL)
  {
    error("MeshDistributedData : returning distributed data of serial mesh");
  }
  return *distdata_;
}
//-----------------------------------------------------------------------------
uint MeshTopology::dim() const
{
  return dim_;
}
//-----------------------------------------------------------------------------
uint MeshTopology::size(uint dim) const
{
  dolfin_assert(dim <= dim_);
  return (dim == 0 ? num_vertices_ : (*this)(dim, 0).num_entities());
}
//-----------------------------------------------------------------------------
uint MeshTopology::global_size(uint dim) const
{
  return (distdata_ ? distdata_->num_global(dim) : this->size(dim));
}
//-----------------------------------------------------------------------------
uint MeshTopology::num_owned(uint dim) const
{
  return (distdata_ ? distdata_->num_owned(dim) : this->size(dim));
}
//-----------------------------------------------------------------------------
uint MeshTopology::num_shared(uint dim) const
{
  return (distdata_ ? distdata_->num_shared(dim) : 0);
}
//-----------------------------------------------------------------------------
uint MeshTopology::num_ghosts(uint dim) const
{
  return (distdata_ ? distdata_->num_ghost(dim) : 0);
}
//-----------------------------------------------------------------------------
bool MeshTopology::is_computed(uint d0, uint d1) const
{
  return connectivity_[d0][d1].is_initialized();
}
//-----------------------------------------------------------------------------
bool MeshTopology::entities_exist(uint d) const
{
  return (d == 0 ? num_vertices_ > 0 : connectivity_[d][0].is_initialized());
}
//-----------------------------------------------------------------------------
void MeshTopology::compute_connectivity(uint d0, uint d1) const
{
  if ((d0 == dim_ && d1 == 0) || connectivity_[d0][d1].is_initialized())
  {
    /*
     *  Return if connectivity exists or if cell -> vertices connectivity, which
     *  is supposed to be provided, is the one called.
     *
     */

    return;
  }
  else if ((d0 == dim_ || d1 == dim_) && !connectivity_[dim_][0].is_initialized())
  {
    /*
     *  For these connectivities, cell -> vertices connectivities should exist
     *
     */

    error("MeshTopology : trying to initialize cell-based connectivity but "
          "the mesh was not provided with cells - vertices connectivity.");
  }
  else if ((d0 > 0 && d1 == 0) || (d0 == dim_ && d1 > 0 && d1 < dim_))
  {
    /*
     *  Compute entities to obtain (tdim, d) and (d, 0) connectivities
     *
     */

    uint const dim = (d1 == 0 ? d0 : d1);
    MeshConnectivity const& cv = (*this)(dim_, 0);
    dolfin_assert(cv.is_initialized());
    MeshConnectivity& ce = connectivity_[dim_][dim];
    MeshConnectivity& ev = connectivity_[dim][0];

    // Initialize local array of entities
    CellType const& cell_type = mesh_.type();
    uint const m = cell_type.num_entities(dim);
    uint const n = cell_type.num_vertices(dim);
    dolfin_assert(this->size(0) > 0);
    Array<uint> * vertex_entities = new Array<uint> [this->size(0)];
    uint ** entities = new uint*[m];
    for (uint e = 0; e < m; ++e)
    {
      entities[e] = new uint[n];
      for (uint v = 0; v < n; ++v)
      {
        entities[e][v] = 0;
      }
    }

    // Collect entities and create cell -> entities connectivities
    EntityKey key(n);
    Array<EntityKey *> entities_list;
    ce.init(this->size(dim_), m);
    for (uint c = 0; c < cv.num_entities(); ++c)
    {
      cell_type.create_entities(entities, dim, cv(c));
      for (uint e = 0; e < m; ++e)
      {
        key.set(entities[e], entities_list.size());
        uint v0 = entities[e][0];
        dolfin_assert(v0 < this->size(0));
        if (vertex_entities[v0].size() == 0)
        {
          entities_list.push_back(new EntityKey(n, entities[e], key.idx));
        }
        else
        {
          // If the first vertex does not contain entity, append new key
          for (uint i = 0; i < vertex_entities[v0].size(); ++i)
          {
            if (*entities_list[vertex_entities[v0][i]] == key)
            {
              key.idx = entities_list[vertex_entities[v0][i]]->idx;
              goto add_entity;
            }
          }
          entities_list.push_back(new EntityKey(n, entities[e], key.idx));
        }

        add_entity: for (uint v = 0; v < n; ++v)
        {
          vertex_entities[entities[e][v]].push_back(key.idx);
        }
        ce.set(c, key.idx, e);
      }
    }

    // Create entity -> vertices connectivities from collected entities
    ev.init(entities_list.size(), n);
    for (uint e = 0; e < entities_list.size(); ++e)
    {
      ev.set(e, entities_list[e]->indices);
    }

    // Cleanup
    entities_list.free();
    for (uint e = 0; e < m; ++e)
    {
      delete[] entities[e];
    }
    delete[] entities;
    delete[] vertex_entities;
  }
  else if (d0 < d1)
  {
    /*
     *  Compute connectivities from transpose.
     *
     */

    MeshConnectivity& c01 = connectivity_[d0][d1];
    compute_connectivity(d1, d0);
    MeshConnectivity const& c10 = connectivity_[d1][d0];

    // Compute from transpose
    Array<uint> conn(this->size(d0), 0);
    dolfin_assert(c10.num_entities() > 0);
    for (uint e1 = 0; e1 < c10.num_entities(); ++e1)
    {
      for (uint e0 = 0; e0 < c10.size(e1); ++e0)
      {
        conn[c10(e1)[e0]]++;
      }
    }
    c01.init(conn);
    //
    conn = 0;
    for (uint e1 = 0; e1 < c10.num_entities(); ++e1)
    {
      for (uint e0 = 0; e0 < c10.size(e1); ++e0)
      {
        c01.set(c10(e1)[e0], e1, conn[c10(e1)[e0]]++);
      }
    }
    dolfin_assert(c01.num_entities() == this->size(d0));
  }
  else if (d0 == d1)
  {
    /*
     *  Compute neighbours for given dimension
     *
     */

    MeshConnectivity& c01 = connectivity_[d0][d1];
    uint const di = (d0 == 0 ? dim_ : 0);

    // Compute connectivity d0 - d - d1 and take intersection
    compute_connectivity(d0, di);
    MeshConnectivity const& c0d = connectivity_[d0][di];
    compute_connectivity(di, d0);
    MeshConnectivity const& cd0 = connectivity_[di][d0];
    std::set<uint> entities;
    Array<uint> conn(this->size(d0), 0);
    for (uint e0 = 0; e0 < c0d.num_entities(); ++e0)
    {
      entities.clear();
      for (uint i = 0; i < c0d.size(e0); ++i)
      {
        uint const e = c0d(e0)[i];
        for (uint j = 0; j < cd0.size(e); ++j)
        {
          uint const e1 = cd0(e)[j];
          // An entity is not a neighbor to itself
          if (e0 != e1)
          {
            entities.insert(e1);
          }
        }
      }
      conn[e0] = entities.size();
    }
    c01.init(conn);
    for (uint e0 = 0; e0 < c0d.num_entities(); ++e0)
    {
      entities.clear();
      for (uint i = 0; i < c0d.size(e0); ++i)
      {
        uint const e = c0d(e0)[i];
        for (uint j = 0; j < cd0.size(e); ++j)
        {
          uint const e1 = cd0(e)[j];
          // An entity is not a neighbor to itself
          if (e0 != e1)
          {
            entities.insert(e1);
          }
        }
      }
      uint pos = 0;
      for (std::set<uint>::iterator it = entities.begin(); it != entities.end();
           ++it, ++pos)
      {
        c01.set(e0, *it, pos);
      }
    }
    dolfin_assert(c01.num_entities() == this->size(d0));
  }
  else
  {
    /*
     *  Compute connectivities between edges/faces and cell/vertices
     *
     */

    MeshConnectivity& c01 = connectivity_[d0][d1];

    // Former code was a special case taking intersection with d = 0
    // Compute connectivity d0 - d - d1 and take intersection
    compute_connectivity(d0, 0);
    MeshConnectivity const& c0v = connectivity_[d0][0];
    compute_connectivity(d1, 0);
    MeshConnectivity const& c1v = connectivity_[d1][0];
    compute_connectivity(0, d1);
    MeshConnectivity const& cv1 = connectivity_[0][d1];
    std::set<uint> entities;
    Array<uint> conn(this->size(d0), 0);
    for (uint e0 = 0; e0 < c0v.num_entities(); ++e0)
    {
      entities.clear();
      for (uint i = 0; i < c0v.size(e0); ++i)
      {
        uint const e = c0v(e0)[i];
        for (uint j = 0; j < cv1.size(e); ++j)
        {
          uint const e1 = cv1(e)[j];
          if (contains(c0v(e0), c0v.size(e0), c1v(e1), c1v.size(e1)))
          {
            entities.insert(e1);
          }
        }
      }
      conn[e0] = entities.size();
    }
    c01.init(conn);
    for (uint e0 = 0; e0 < c0v.num_entities(); ++e0)
    {
      entities.clear();
      for (uint i = 0; i < c0v.size(e0); ++i)
      {
        uint const e = c0v(e0)[i];
        for (uint j = 0; j < cv1.size(e); ++j)
        {
          uint const e1 = cv1(e)[j];
          if (contains(c0v(e0), c0v.size(e0), c1v(e1), c1v.size(e1)))
          {
            entities.insert(e1);
          }
        }
      }
      // Add the connected entities
      uint pos = 0;
      for (std::set<uint>::iterator it = entities.begin(); it != entities.end();
           ++it, ++pos)
      {
        c01.set(e0, *it, pos);
      }
    }
    dolfin_assert(c01.num_entities() == this->size(d0));
  }

  /*
   * If the created connectivity needs ordering then reset the flag to trigger
   * reordering
   */

  if(mesh_.type().connectivity_needs_ordering(d0, d1))
  {
    is_ordered_ = false;
  }

  message(1, "MeshTopology : computed connectivity (%u, %u)", d0, d1);
}
//-----------------------------------------------------------------------------
bool MeshTopology::is_distributed() const
{
  return (distdata_ != NULL);
}
//-----------------------------------------------------------------------------
void MeshTopology::disp() const
{
  section("MeshTopology");
  //---
  cout << "Dimension: " << dim_ << endl;
  skip();
  begin("Number of entities:");
  if (num_vertices_ == 0)
  {
    cout << "empty" << endl;
  }
  else
  {
    cout << "0: " << num_vertices_ << endl;
    for (uint d = 1; d <= dim_; ++d)
    {
      if (connectivity_[d][0].is_initialized())
      {
        cout << d << ": " << connectivity_[d][0].num_entities() << endl;
      }
      else
      {
        cout << d << ": " << "not computed" << endl;
      }
    }
  }
  end();
  skip();
  begin("Connectivity:");
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
  //---
  end();
  skip();
}
//-----------------------------------------------------------------------------
void MeshTopology::update_token()
{
  timestamp_ = std::time(NULL);
}
//-----------------------------------------------------------------------------
int MeshTopology::token() const
{
  return timestamp_ ^ size(0) ^ size(dim_);
}
//-----------------------------------------------------------------------------
void MeshTopology::check() const
{
  message("MeshTopology: check");

  /**
   *  CHECK:
   *
   *  Mesh entities connectivities should follow the convention provided by the
   *  cell type.
   *
   */

  message("MeshTopology: check connectivities");

  for (uint d0 = 0; d0 <= dim_; ++d0)
  {
    for (uint d1 = 0; d1 <= dim_; ++d1)
    {
      if (this->is_computed(d0, d1))
      {
        message("(%u,%u)", d0, d1);
        connectivity_[d0][d1].check();
      }
    }
  }

  message("MeshTopology: check ordering of entities on cells");
  Mesh& mesh = const_cast<Mesh&>(mesh_);
  for (CellIterator c(mesh); !c.end(); ++c)
  {
    mesh.type().check(*c);
  }

  if (this->is_distributed())
  {
    message("MeshTopology: check distributed data");
    //distdata_->check();
  }
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
