// Copyright (C) 2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Aurelien Larcher, 2016.
// Full rewrite.
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
    ini_vertices_(false),
    connectivity_(NULL),
    distdata_(NULL),
    timestamp_(0)
{
}
//-----------------------------------------------------------------------------
MeshTopology::MeshTopology(MeshTopology const& other) :
    mesh_(other.mesh_),
    dim_(0),
    num_vertices_(0),
    ini_vertices_(false),
    connectivity_(NULL),
    distdata_(NULL),
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
  ini_vertices_ = other.ini_vertices_;
  if (other.connectivity_ != NULL)
  {
    connectivity_ = new MeshConnectivity*[dim_ + 1];
    for (uint d0 = 0; d0 <= dim_; ++d0)
    {
      connectivity_[d0] = new MeshConnectivity[dim_ + 1];
      for (uint d1 = 0; d1 <= dim_; ++d1)
      {
        connectivity_[d0][d1] = other.connectivity_[d0][d1];
      }
    }
  }
  if (other.distdata_ != NULL)
  {
    distdata_ = new MeshDistributedData(*other.distdata_);
  }
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
  if (connectivity_ != NULL)
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
void MeshTopology::init(uint dim, bool distribute /* = true */)
{
  if (connectivity_ != NULL)
  {
    error("MeshTopology : clear instance before reinitializing");
  }
  dim_ = dim;
  connectivity_ = new MeshConnectivity*[dim + 1];
  for (uint d = 0; d <= dim; ++d)
  {
    connectivity_[d] = new MeshConnectivity[dim + 1];
  }
  if (distribute)
  {
    distdata_ = new MeshDistributedData(dim_);
  }
  //
  update_token();
}
//-----------------------------------------------------------------------------
void MeshTopology::init(uint dim, uint num_local, uint num_global /* = 0 */)
{
  init(NULL, dim, num_local, num_global);
}
//-----------------------------------------------------------------------------
void MeshTopology::init(uint * connectivity, uint dim, uint num_local,
                        uint num_global /* = 0 */)
{
  if (connectivity_ == NULL)
  {
    error("MeshTopology : initializing entities of dimension %u but topology "
          "dimension is initialized", dim);
  }
  // NOTE: cases have to be split that way to handle point meshes for which the
  //       cell dimension is equal to the vertex dimension
  // Vertices
  if(dim == 0)
  {
    num_vertices_ = num_local;
    ini_vertices_ = true;
  }
  // Edges/Faces
  if (dim < 0 && dim < dim_)
  {
    if(num_local > 0 && num_vertices_ == 0)
    {
      error("MeshTopology : initializing non-zero number of entities of "
            "dimension %u but topology containes zero vertices", dim);
    }
    // Well Well Well *erm* *erm* *erm* OOP gone wrong
    connectivity_[dim_][0].init(connectivity, num_local,
                                mesh_.type().num_vertices(dim));
  }
  // Cells
  if (dim_ == dim)
  {
    if(num_local > 0 && num_vertices_ == 0)
    {
      error("MeshTopology : initializing non-zero number of cells but topology "
            "contains zero vertices");
    }
    // Well Well Well *erm* *erm* *erm* OOP gone wrong
    connectivity_[dim_][0].init(connectivity, num_local,
                                mesh_.type().num_vertices(dim));
  }
  // Overflow
  if (dim_ < dim)
  {
    error("MeshTopology : initializing entities of dimension %u but topology"
          "dimension is %u", dim, dim_);
  }
  // Set size of distributed data
  if (distdata_ != NULL)
  {
    if (num_global > 0 && num_global < num_local)
    {
      error("MeshTopology : number of global entities lower than number of "
            " local entities %u < %u", num_local, num_global);
    }
    this->distdata()[dim].set_size(num_local, num_global);
  }
  else
  {
    // In serial, require that the number of global entities is not initialized
    // of is equal to the number of local entities
    if ((num_global > 0) && (num_local != num_global))
    {
      error("MeshTopology : invalid number of global entities set in serial");
    }
  }

}
//-----------------------------------------------------------------------------
void MeshTopology::clear()
{
  // Clear parallel data structures
  delete distdata_;
  distdata_ = NULL;

  // Delete mesh connectivity
  if (connectivity_ != NULL)
  {
    for (uint d = 0; d <= dim_; ++d)
    {
      delete[] connectivity_[d];
    }
    delete[] connectivity_;
  }
  connectivity_ = NULL;
  timestamp_ = 0;
  num_vertices_ = 0;
  ini_vertices_ = false;
  dim_ = 0;
}
//-----------------------------------------------------------------------------
void MeshTopology::finalize()
{
  if(!connectivity_[dim_][0].is_initialized())
  {
    warning("MeshTopology : cell -> vertices connectivity does not exist");
  }

  // Reorder cells according to UFC convention
  reorder();

  // Finalize distributed data
  if (distdata_ != NULL)
  {
    for (uint d = 0; d <= dim_; ++d)
    {
      if (this->entities_exist(d))
      {
        DistributedData& ddata = this->distdata()[d];
        if (!ddata.is_finalized())
        {
          ddata.finalize();
        }
        if (this->size(d) != ddata.local_size())
        {
          error("MeshEditor : vertex size mismatch between topology '%u' and "
                "distributed data '%u'", this->size(d), ddata.local_size());
        }
      }
    }
  }
  // Do not renumber automatically !
  // This would cause issues for boundary meshes and some mesh algorithms.
}
//-----------------------------------------------------------------------------
CellType const& MeshTopology::type(uint i) const
{
  return mesh_.type();
}
//-----------------------------------------------------------------------------
void MeshTopology::remap(uint dim, Array<uint> const& mapping)
{
  if (connectivity_ != NULL)
  {
    uint d0 = dim;
    for (uint d1 = 0; d1 <= dim_; ++d1)
    {
      if (connectivity_[d0][d1].size() > 0)
      {
        connectivity_[d0][d1].remap_left(mapping);
        update_token();
      }
      if (connectivity_[d1][d0].size() > 0)
      {
        connectivity_[d1][d0].remap_right(mapping);
        update_token();
      }
    }
    reorder();
    // Remap distributed data
    if (distdata_ != NULL)
    {
      (*distdata_)[d0].remap_numbering(mapping);
    }
  }
}
//-----------------------------------------------------------------------------
MeshConnectivity& MeshTopology::operator()(uint d0, uint d1)
{
  dolfin_assert(d0 <= dim_ && d1 <= dim_);
  if(connectivity_ != NULL && !connectivity_[d0][d1].is_initialized())
  {
    compute_connectivity(d0, d1);
  }
  return connectivity_[d0][d1];
}
//-----------------------------------------------------------------------------
MeshConnectivity const& MeshTopology::operator()(uint d0, uint d1) const
{
  dolfin_assert(d0 <= dim_ && d1 <= dim_);
  if(connectivity_ != NULL && !connectivity_[d0][d1].is_initialized())
  {
    compute_connectivity(d0, d1);
  }
  return connectivity_[d0][d1];
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
bool MeshTopology::is_computed(uint d0, uint d1) const
{
  dolfin_assert(d0 <= dim_);
  dolfin_assert(d1 <= dim_);
  return connectivity_[d0][d1].is_initialized();
}
//-----------------------------------------------------------------------------
bool MeshTopology::entities_exist(uint dim) const
{
  dolfin_assert(dim <= dim_);
  return (dim == 0 ?
            (ini_vertices_ == true) : connectivity_[dim][0].is_initialized());
}
//-----------------------------------------------------------------------------
bool MeshTopology::is_distributed() const
{
  return (distdata_ != NULL);
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
uint MeshTopology::global_size(uint dim) const
{
  return (distdata_ ? (*distdata_)[dim].global_size() : this->size(dim));
}
//-----------------------------------------------------------------------------
uint MeshTopology::offset(uint dim) const
{
  return (distdata_ ? (*distdata_)[dim].offset() : 0);
}
//-----------------------------------------------------------------------------
uint MeshTopology::num_owned(uint dim) const
{
  return (distdata_ ? (*distdata_)[dim].num_owned() : this->size(dim));
}
//-----------------------------------------------------------------------------
uint MeshTopology::num_shared(uint dim) const
{
  return (distdata_ ? (*distdata_)[dim].num_shared() : 0);
}
//-----------------------------------------------------------------------------
uint MeshTopology::num_ghost(uint dim) const
{
  return (distdata_ ? (*distdata_)[dim].num_ghost() : 0);
}
//-----------------------------------------------------------------------------
void MeshTopology::compute_connectivity(uint d0, uint d1) const
{
  if (connectivity_ == NULL)
  {
    error("MeshTopology : connectivity array does not exist");
  }

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
    Array<uint> * vertex_entities = (this->size(0) == 0 ?
                                      NULL : new Array<uint> [this->size(0)]);
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

    // New entities have been computed: trigger renumbering
    renumber();
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
   *
   */

  if(mesh_.type().connectivity_needs_ordering(d0, d1))
  {
    reorder();
  }

  message(1, "MeshTopology : computed connectivity (%u, %u)", d0, d1);
}
//-----------------------------------------------------------------------------
void MeshTopology::reorder() const
{
  //FIXME: this test ensures that boundary meshes are not reordered
  if (mesh_.geometry().dim() == dim_)
  {
    message(1, "MeshTopology : order");
    CellType const& cell_type = mesh_.type();
    MeshTopology& topology = const_cast<MeshTopology&>(*this);
    uint const num_cells = this->size(dim_);
    for (uint i = 0; i < num_cells; ++i)
    {
      cell_type.order_entities(topology, i);
    }
  }
}
//-----------------------------------------------------------------------------
void MeshTopology::renumber() const
{
  if (distdata_ == NULL)
  {
    return;
  }

  message(1, "MeshTopology : renumber");

  MeshTopology& topology = const_cast<MeshTopology&>(*this);
  if(!MeshRenumber::renumber(topology))
  {
    warning("MeshTopology: triggered mesh renumbering for nothing");
  }
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
      dolfin_assert(connectivity_ != NULL);
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
  if (connectivity_ == NULL)
  {
    cout << "empty" << endl;
  }
  else
  {
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
  }
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

} /* namespace dolfin */
