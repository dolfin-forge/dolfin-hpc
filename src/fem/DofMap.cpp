// Copyright (C) 2007-2008 Anders Logg and Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Martin Alnes, 2008
// Modified by Niclas Jansson, 2009
// Modified by Aurélien Larcher, 2013  (Pk bug, extension and partial rewrite)
//
// First added:  2007-03-01
// Last changed: 2009-11-01

#include <dolfin/fem/DofMap.h>

#include <dolfin/config/dolfin_config.h>
#include <dolfin/common/types.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/fem/DofMapCache.h>
#include <dolfin/fem/UFCCell.h>
#include <dolfin/fem/SubSystem.h>
#include <dolfin/common/Array.h>
#include <dolfin/elements/ElementLibrary.h>
#include <dolfin/fem/UFC.h>
#include <dolfin/main/MPI.h>

#include <dolfin/mesh/BoundaryMesh.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/MeshData.h>
#include <dolfin/mesh/MeshFunction.h>
#include <cstring>
#include <cstdlib>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

namespace dolfin
{

//-----------------------------------------------------------------------------
#if UFC_VERSION_MAJOR == 1
std::string const DofMap::SIGN_PREFIX = "FFC dof map for ";
#elif UFC_VERSION_MAJOR == 2
std::string const DofMap::SIGN_PREFIX = "FFC dofmap for ";
#endif

//-----------------------------------------------------------------------------
DofMap::DofMap(Mesh& mesh, ufc::form const& form, uint const i) :
    MeshDependent(mesh),
    ufc_mesh_(mesh),
    type_(DofMap::ufc_default),
    offset_(0),
    ufc_dofmap_(form.create_dofmap(i)),
    hash_(make_hash(mesh, *ufc_dofmap_)),
    local_size_(0),
    vertex_map_(NULL),
    pretabulated_dofmap_(NULL),
    pretabulated_dofmap_size_(0)
{
  init();
}

//-----------------------------------------------------------------------------
DofMap::DofMap(Mesh& mesh, ufc::dofmap& dofmap, bool const owner) :
    MeshDependent(mesh),
    ufc_mesh_(mesh),
    type_(DofMap::ufc_default),
    offset_(0),
    ufc_dofmap_((owner ? &dofmap : dofmap.create())),
    hash_(make_hash(mesh, *ufc_dofmap_)),
    local_size_(0),
    vertex_map_(NULL),
    pretabulated_dofmap_(NULL),
    pretabulated_dofmap_size_(0)
{
  init();
}

//-----------------------------------------------------------------------------
DofMap::DofMap(DofMap const& dofmap, uint i) :
    MeshDependent(dofmap.mesh()),
    ufc_mesh_(dofmap.mesh()),
    type_(DofMap::ufc_default),
    offset_(0),
    ufc_dofmap_(dofmap.create_sub_dofmap(i)),
    hash_(make_hash(mesh(), *ufc_dofmap_)),
    local_size_(0),
    vertex_map_(NULL),
    pretabulated_dofmap_(NULL),
    pretabulated_dofmap_size_(0)
{
  message(1, "Extracted dof map for subspace: %s", ufc_dofmap_->signature());

  init();
}

//-----------------------------------------------------------------------------
DofMap::DofMap(DofMap const& dofmap, Array<uint> const& sub_system,
               uint& offset) :
    MeshDependent(dofmap.mesh()),
    ufc_mesh_(dofmap.mesh()),
    type_(DofMap::ufc_default),
    offset_(0),
    ufc_dofmap_(dofmap.create_sub_dofmap(sub_system, offset_)),
    hash_(make_hash(mesh(), *ufc_dofmap_)),
    local_size_(0),
    vertex_map_(NULL),
    pretabulated_dofmap_(NULL),
    pretabulated_dofmap_size_(0)
{
  // Check that dof map has not be re-ordered
  offset = offset_;
  message(1, "Extracted dof map for sub system: %s", ufc_dofmap_->signature());
  message(1, "Offset for sub system: %d", offset);

  // Reset offset
  offset_ = 0;
  init();
}

//-----------------------------------------------------------------------------
DofMap::~DofMap()
{
  delete[] pretabulated_dofmap_;
  delete[] vertex_map_;
  while (!flattened_.empty())
  {
    delete flattened_.back();
    flattened_.pop_back();
  }
  delete ufc_dofmap_;
}

//-----------------------------------------------------------------------------
bool DofMap::operator ==(DofMap const& other) const
{
  return (this->hash() == other.hash());
}

//-----------------------------------------------------------------------------
bool DofMap::operator !=(DofMap const& other) const
{
  return !(*this == other);
}

//-----------------------------------------------------------------------------
DofMap& DofMap::acquire(Mesh& mesh, Form const& form, uint const i)
{
  return DofMapCache::instance().acquire(mesh, form, i);
}

//-----------------------------------------------------------------------------
DofMap& DofMap::acquire(Mesh& mesh, ufc::dofmap& dofmap, bool owner)
{
  return DofMapCache::instance().acquire(mesh, dofmap, owner);
}

//-----------------------------------------------------------------------------
void DofMap::release(DofMap& dofmap)
{
  return DofMapCache::instance().release(dofmap);
}

//-----------------------------------------------------------------------------
ufc::dofmap* DofMap::create_sub_dofmap(Array<uint> const& sub_system) const
{
  // Reset offset
  uint local_offset = 0;

  // Recursively extract sub dof map
  ufc::dofmap* sub_dofmap = DofMap::create_sub_dofmap(*ufc_dofmap_, sub_system,
                                                      local_offset);
  message(1, "Extracted ufc dof map for sub system: %s",
          sub_dofmap->signature());
  message(1, "Local offset for sub system: %d", local_offset);

  return sub_dofmap;
}

//-----------------------------------------------------------------------------
ufc::dofmap* DofMap::create_sub_dofmap(Array<uint> const& sub_system,
                                       uint& local_offset) const
{
  // Reset offset
  local_offset = 0;

  // Recursively extract sub dof map
  ufc::dofmap* sub_dofmap = DofMap::create_sub_dofmap(*ufc_dofmap_, sub_system,
                                                      local_offset);
  message(1, "Extracted ufc dof map for sub system: %s",
          sub_dofmap->signature());
  message(1, "Local offset for sub system: %d", local_offset);

  return sub_dofmap;
}

//-----------------------------------------------------------------------------
ufc::dofmap* DofMap::create_sub_dofmap(ufc::dofmap const& dofmap,
                                       Array<uint> const& sub_system,
                                       uint& local_offset)
{
  // Check that a sub system has been specified
  if (sub_system.size() == 0)
  {
    //error("Unable to extract sub system (no sub system specified).");
    return dofmap.create();
  }

  // Check if there are any sub systems
  if (dofmap.num_sub_dofmaps() == 0)
  {
    error("Unable to extract sub system (there are no sub systems).");
  }

  // Check the number of available sub systems
  if (sub_system[0] >= dofmap.num_sub_dofmaps())
  {
    error("Unable to extract sub system %d (only %d sub systems defined).",
          sub_system[0], dofmap.num_sub_dofmaps());
  }

  // Add to offset if necessary
  for (uint i = 0; i < sub_system[0]; i++)
  {
    ufc::dofmap * ufc_sub_dofmap = dofmap.create_sub_dofmap(i);
    local_offset += ufc_sub_dofmap->local_dimension();
    delete ufc_sub_dofmap;
  }

  // Create sub system
  ufc::dofmap* sub_dofmap = dofmap.create_sub_dofmap(sub_system[0]);

  // Return sub system if sub sub system should not be extracted
  if (sub_system.size() == 1) return sub_dofmap;

  // Otherwise, recursively extract the sub sub system
  Array<uint> sub_sub_system;
  for (uint i = 1; i < sub_system.size(); i++)
  {
    sub_sub_system.push_back(sub_system[i]);
  }
  ufc::dofmap* sub_sub_dofmap = DofMap::create_sub_dofmap(*sub_dofmap,
                                                          sub_sub_system,
                                                          local_offset);
  delete sub_dofmap;

  return sub_sub_dofmap;
}

//-----------------------------------------------------------------------------
ufc::dofmap* DofMap::create_sub_dofmap(UFCMesh& ufc_mesh,
                                       ufc::dofmap const& dofmap,
                                       Array<uint> const& sub_system,
                                       uint& global_offset)
{
  // Check that a sub system has been specified
  if (sub_system.size() == 0)
  {
    //error("Unable to extract sub system (no sub system specified).");
    return dofmap.create();
  }

  // Check if there are any sub systems
  if (dofmap.num_sub_dofmaps() == 0)
  {
    error("Unable to extract sub system (there are no sub systems).");
  }

  // Check the number of available sub systems
  if (sub_system[0] >= dofmap.num_sub_dofmaps())
  {
    error("Unable to extract sub system %d (only %d sub systems defined).",
          sub_system[0], dofmap.num_sub_dofmaps());
  }

  // Add to offset if necessary
  for (uint i = 0; i < sub_system[0]; i++)
  {
    ufc::dofmap * ufc_sub_dofmap = dofmap.create_sub_dofmap(i);
    //Avoid creating a DOLFIN dofmap, just calling the static init_ufc function
    DofMap::initUFC(ufc_mesh, *ufc_sub_dofmap);
    global_offset += ufc_sub_dofmap->global_dimension();
    delete ufc_sub_dofmap;
  }

  // Create sub system
  ufc::dofmap* sub_dofmap = dofmap.create_sub_dofmap(sub_system[0]);

  // Return sub system if sub sub system should not be extracted
  if (sub_system.size() == 1) return sub_dofmap;

  // Otherwise, recursively extract the sub sub system
  Array<uint> sub_sub_system;
  for (uint i = 1; i < sub_system.size(); i++)
  {
    sub_sub_system.push_back(sub_system[i]);
  }
  ufc::dofmap* sub_sub_dofmap = DofMap::create_sub_dofmap(ufc_mesh, *sub_dofmap,
                                                          sub_sub_system,
                                                          global_offset);
  delete sub_dofmap;

  return sub_sub_dofmap;
}

//-----------------------------------------------------------------------------
void DofMap::initUFC(UFCMesh& ufc_mesh, ufc::dofmap& ufc_dofmap)
{
  Mesh * mesh = const_cast<Mesh *>(ufc_mesh.mesh);

  if (ufc_dofmap.geometric_dimension() != mesh->topology().dim())
  {
    dolfin::error("The geometric dimension of the dof map is not equal to the "
                  "topological dimension of the mesh.");
  }

  // Order vertices, so entities will be created correctly according to convention
  mesh->order();

  // Initialize mesh entities used by dof map
  for (uint d = 0; d <= mesh->topology().dim(); d++)
  {
    if (ufc_dofmap.needs_mesh_entities(d))
    {
      mesh->init(d);
    }
  }

  // Initialize UFC mesh data (must be done after entities are created)
  ufc_mesh.init(*mesh);

  // Initialize UFC dof map
  const bool init_cells = ufc_dofmap.init_mesh(ufc_mesh);
  if (init_cells)
  {
    CellIterator cell(*mesh);
    UFCCell ufc_cell(*cell);
    for (; !cell.end(); ++cell)
    {
      ufc_cell.update(*cell, mesh->distdata());
      ufc_dofmap.init_cell(ufc_mesh, ufc_cell);
    }
    ufc_dofmap.init_cell_finalize();
  }
}

//-----------------------------------------------------------------------------
void DofMap::init()
{
  //dolfin_debug("Initializing dof map...");

  // Initialize UFC data structures
  DofMap::initUFC(ufc_mesh_, *ufc_dofmap_);

  // Build the DOLFIN dofmap
  build();

  // Information for mixed elements
  uint const nb_sub = this->num_sub_dofmaps();
  if (nb_sub > 0)
  {
    // Set offsets and local dimensions
    uint off = 0;
    for (uint i = 0; i < nb_sub; ++i)
    {
      ufc::dofmap * subdm = this->create_sub_dofmap(i);
      sub_dofmaps_dims_.push_back(subdm->local_dimension());
      sub_dofmaps_offs_.push_back(off);
      off += subdm->local_dimension();
      delete subdm;
    }
  }
  else
  {
    sub_dofmaps_dims_.push_back(this->local_dimension());
    sub_dofmaps_offs_.push_back(0);
  }
}
//-----------------------------------------------------------------------------
void DofMap::tabulate_dofs(uint* dofs, UFCCell const& ufc_cell,
                           uint cell_index) const
{
  tabulate_dofs(dofs, ufc_cell, *ufc_cell.cell);
}

//-----------------------------------------------------------------------------
void DofMap::tabulate_dofs(uint* dofs, ufc::cell const& ufc_cell,
                           Cell const& cell) const
{
  dolfin_assert(dofs != NULL);
  // Either lookup pretabulated values (if build() has been called)
  // or ask the ufc::dofmap to tabulate the values
  switch (type_)
    {
    case real_space:
      {
        uint const rank = MPI::processNumber();
        for (uint i = 0; i < local_dimension(); ++i)
        {
          dofs[i] = rank * local_dimension() + i;
        }
      }
      break;
    case scalar_p1:
      {
        for (uint i = 0; i < local_dimension(); ++i)
        {
          dofs[i] = ufc_cell.entity_indices[0][i];
          dolfin_assert(mesh().distdata().has_global(dofs[i], 0));
        }
      }
      break;
    case scalar_dg0:
      {
        dofs[0] = ufc_cell.index;
        dolfin_assert(
            mesh().distdata().has_global(dofs[0], mesh().topology().dim()));
      }
      break;
    case vector_p1:
      {
        uint const num_cellverts = cell.numEntities(0);
        for (uint k = 0; k < num_leaf_spaces_; ++k)
        {
          for (uint i = 0; i < num_cellverts; ++i)
          {
            dofs[i + k * num_cellverts] = vertex_map_[cell.entities(0)[i]] + k;
          }
        }
      }
      break;
    case vector_dg0:
      {
        for (uint i = 0; i < num_leaf_spaces_; ++i)
        {
          dofs[i] = (cell.index() + i * mesh().numCells()) + offset_;
        }
      }
      break;
    case generic:
      {
        uint const local_dim = ufc_dofmap_->local_dimension();
        std::memcpy(dofs, &pretabulated_dofmap_[local_dim * cell.index()],
                    sizeof(uint) * local_dim);
      }
      break;
    case ufc_default:
      {
        ufc_dofmap_->tabulate_dofs(dofs, ufc_mesh_, ufc_cell);
      }
      break;
    default:
      error("Unknown dofmap type.");
      break;
    }
}

//-----------------------------------------------------------------------------
Array<ufc::dofmap const *> const& DofMap::flatten() const
{
  if (flattened_.empty())
  {
    flatten(ufc_dofmap_, flattened_);
  }
  return flattened_;
}

//-----------------------------------------------------------------------------
void DofMap::flatten(ufc::dofmap const * dofmap,
                     Array<ufc::dofmap const *>& stack, uint maxlevel) const
{
  // Single root element or max level is set to zero, return immediately
  if (dofmap->num_sub_dofmaps() == 0 || maxlevel == 0)
  {
    stack.push_back(dofmap->create());
    return;
  }
  // Go one level down
  for (uint s = 0; s < dofmap->num_sub_dofmaps(); ++s)
  {
    ufc::dofmap const * sub = dofmap->create_sub_dofmap(s);
    if (sub->num_sub_dofmaps() == 0)
    {
      // Leaf dofmap
      stack.push_back(sub);
    }
    else
    {
      // Branch
      flatten(sub, stack, maxlevel - 1);
    }
  }
}

//-----------------------------------------------------------------------------
Array<uint> const& DofMap::sub_dofmaps_dimensions() const
{
  return sub_dofmaps_dims_;
}

//-----------------------------------------------------------------------------
Array<uint> const& DofMap::sub_dofmaps_offsets() const
{
  return sub_dofmaps_offs_;
}

//-----------------------------------------------------------------------------
bool DofMap::renumbered() const
{
  return (pretabulated_dofmap_ || type_ > -1 || vertex_map_);
}

//-----------------------------------------------------------------------------
uint DofMap::local_size() const
{
  return local_size_;
}

//-----------------------------------------------------------------------------
uint const * DofMap::dofsmapping() const
{
  if (pretabulated_dofmap_ == NULL)
  {
    pretabulateAllDofs();
  }
  return pretabulated_dofmap_;
}

//-----------------------------------------------------------------------------
uint DofMap::dofsmapping_size() const
{
  return pretabulated_dofmap_size_;
}

//--------------------------------------------------------------------------
void DofMap::pretabulateAllDofs() const
{
  delete[] pretabulated_dofmap_;
  pretabulated_dofmap_ = new uint[pretabulated_dofmap_size_];
  uint *ip = &pretabulated_dofmap_[0];
  uint const local_dim = this->local_dimension();
  MeshDistributedData& distdata = mesh().distdata();
  CellIterator cell(mesh());
  UFCCell ufc_cell(*cell);
  for (; !cell.end(); ++cell)
  {
    // cell indices for real valued function
    ufc_cell.update(*cell, distdata);
    this->tabulate_dofs(ip, ufc_cell, *cell);
    ip += local_dim;
  }
  type_ = generic;
}

//-----------------------------------------------------------------------------
void DofMap::build()
{

  // Cleanup dofmap structures and set dimension of pretabulated array
  delete[] pretabulated_dofmap_;
  pretabulated_dofmap_ = NULL;
  pretabulated_dofmap_size_ = this->local_dimension() * mesh().numCells();
  map_.clear();

  if (mesh().is_distributed())
  {
#ifdef HAVE_MPI
    Mesh& thismesh = this->mesh();
    thismesh.renumber();

    uint pe_size = MPI::numProcesses();
    uint rank = MPI::processNumber();

    // Determine type of dofmap numbering and build
    Array<ufc::dofmap const *> const& flt = this->flatten();
    num_leaf_spaces_ = flt.size();
    bool can_vectorize = (num_leaf_spaces_ == 1) ? false : true;
    char const * dm0 = flt[0]->signature();
    for (uint i = 1; i < num_leaf_spaces_; ++i)
    {
      if (std::strcmp(dm0, flt[i]->signature()) != 0)
      {
        can_vectorize = false;
        break;
      }
    }

    // Build
    if (ufc_dofmap_->global_dimension() == ufc_dofmap_->local_dimension())
    {
      type_ = real_space;
      local_size_ = ufc_dofmap_->local_dimension();
    }
    else if (ufc_dofmap_->global_dimension() == thismesh.global_numVertices())
    {
      // Scalar Lagrange P1
      type_ = scalar_p1;
      local_size_ = thismesh.distdata().num_owned(0);

      //DEBUG: Create ghosts list
      MeshDistributedData& distdata = thismesh.distdata();
      for (MeshGhostIterator git(distdata, 0); !git.end(); ++git)
      {
        ghosts_.insert(distdata.get_vertex_global(git.index()));
      }
    }
    else if (ufc_dofmap_->global_dimension() == thismesh.global_numCells())
    {
      // Scalar Discontinuous Lagrange P0
      type_ = scalar_dg0;
      local_size_ = thismesh.numCells();

      // No ghosted dofs
    }
    else if (can_vectorize
        && (ufc_dofmap_->global_dimension()
            == num_leaf_spaces_ * thismesh.global_numVertices()))
    {
      // Vector Lagrange P1
      type_ = vector_p1;
      uint vdim = num_leaf_spaces_;
      uint num_local = thismesh.distdata().num_owned(0);

      uint num_dofs = vdim * num_local;
      uint offset = 0;

#if ( MPI_VERSION > 1 )
      MPI_Exscan(&num_dofs, &offset, 1, MPI_UNSIGNED, MPI_SUM,
                 MPI::DOLFIN_COMM);
#else
      MPI_Scan(&num_dofs, &offset, 1, MPI_UNSIGNED, MPI_SUM, MPI::DOLFIN_COMM);
      offset -= num_dofs;
#endif
      _map<uint, uint> v_offset;

      for (VertexIterator v(thismesh); !v.end(); ++v)
      {
        if (!v->is_ghost())
        {
          v_offset[thismesh.distdata().get_global(*v)] = offset;
          offset += vdim;
        }
      }

      Array<uint> *ghost_buff = new Array<uint> [pe_size];
      for (MeshGhostIterator iter(thismesh.distdata(), 0); !iter.end(); ++iter)
      {
        ghost_buff[iter.owner()].push_back(
            thismesh.distdata().get_vertex_global(iter.index()));
      }

      MPI_Status status;
      Array<uint> send_buff;
      uint src, dest;
      uint recv_size = thismesh.distdata().num_ghost(0);
      int recv_count, recv_size_gh, send_size;

      for (uint i = 0; i < pe_size; i++)
      {
        send_size = ghost_buff[i].size();
        MPI_Reduce(&send_size, &recv_size_gh, 1, MPI_INT, MPI_SUM, i,
                   MPI::DOLFIN_COMM);
      }

      uint *recv_ghost = new uint[recv_size_gh];
      uint *recv_buff = new uint[recv_size];

      for (uint j = 1; j < pe_size; j++)
      {
        src = (rank - j + pe_size) % pe_size;
        dest = (rank + j) % pe_size;

        MPI_Sendrecv(&ghost_buff[dest][0], ghost_buff[dest].size(),
                     MPI_UNSIGNED, dest, 1, recv_ghost, recv_size_gh,
                     MPI_UNSIGNED, src, 1, MPI::DOLFIN_COMM, &status);
        MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);

        for (int k = 0; k < recv_count; k++)
          send_buff.push_back(v_offset[recv_ghost[k]]);

        MPI_Sendrecv(&send_buff[0], send_buff.size(), MPI_UNSIGNED, src, 2,
                     recv_buff, recv_size, MPI_UNSIGNED, dest, 2,
                     MPI::DOLFIN_COMM, &status);
        MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);

        for (int j = 0; j < recv_count; j++)
          v_offset[ghost_buff[dest][j]] = recv_buff[j];

        send_buff.clear();
      }

      delete[] recv_ghost;
      delete[] recv_buff;

      delete[] vertex_map_;
      vertex_map_ = NULL;

      vertex_map_ = new uint[thismesh.numVertices()];
      for (VertexIterator v(thismesh); !v.end(); ++v)
      {
        vertex_map_[v->index()] =
            v_offset[thismesh.distdata().get_vertex_global(v->index())];
      }

      v_offset.clear();

      local_size_ = vdim * num_local;

      for (uint i = 0; i < pe_size; i++)
      {
        ghost_buff[i].clear();
      }
      delete[] ghost_buff;

      //DEBUG: Create ghosts list
      MeshDistributedData& distdata = thismesh.distdata();
      for (MeshGhostIterator git(distdata, 0); !git.end(); ++git)
      {
        for (uint k = 0; k < ufc_dofmap_->num_sub_dofmaps(); ++k)
        {
          uint gdof = vertex_map_[distdata.get_vertex_global(git.index())] + k;
          ghosts_.insert(gdof);
        }
      }
    }
    else if (can_vectorize
        && (ufc_dofmap_->global_dimension()
            == num_leaf_spaces_ * thismesh.global_numCells()))
    {
      // Vector Discontinuous Lagrange P0
      type_ = vector_dg0;
      uint vdim = num_leaf_spaces_;
      uint num_dofs = vdim * thismesh.numCells();
      uint offset = 0;

#if ( MPI_VERSION > 1 )
      MPI_Exscan(&num_dofs, &offset, 1, MPI_UNSIGNED, MPI_SUM,
                 MPI::DOLFIN_COMM);
#else
      MPI_Scan(&num_dofs, &offset, 1, MPI_UNSIGNED, MPI_SUM, MPI::DOLFIN_COMM);
      offset -= num_dofs;
#endif
      offset_ = offset;
      local_size_ = vdim * thismesh.numCells();

      // No ghosted dofs
    }
    else
    {
      type_ = generic;

      _set<uint> owned_dofs;
      _set<uint> shared_dofs;  // shared dofs that are not ghosted (!)
      _set<uint> ghost_dofs;
      _map<uint, std::vector<uint> > dof2index;

      //--- Attribute ownership -----------------------------------------------
      // HPC distribution cannot be used with VectorElements
      bool const hpc_distribution = (ufc_dofmap_->num_sub_dofmaps() == 0);
      if (hpc_distribution)
      {
        distributeByVote(ufc_mesh_, ufc_dofmap_, owned_dofs, shared_dofs,
                         ghost_dofs, dof2index);
      }
      else
      {
        distributeByEntities(ufc_mesh_, ufc_dofmap_, owned_dofs, shared_dofs,
                             ghost_dofs, dof2index);
      }
      dolfin_assert(owned_dofs.size() > shared_dofs.size());
#if DEBUG
      message("Owned  dofs : %d", owned_dofs.size());
      message("Shared dofs : %d", shared_dofs.size());
      message("Ghost  dofs : %d", ghost_dofs.size());
#endif

      //--- Global renumbering ------------------------------------------------

      // Compute offset for owned and non shared dofs
      uint range = owned_dofs.size();
      local_size_ = range;
      uint offset = 0;
#if ( MPI_VERSION > 1 )
      MPI_Exscan(&range, &offset, 1, MPI_UNSIGNED, MPI_SUM, MPI::DOLFIN_COMM);
#else
      MPI_Scan(&range, &offset, 1, MPI_UNSIGNED, MPI_SUM, MPI::DOLFIN_COMM);
      offset -= range;
#endif
      uint const offset_ = offset;

      // Compute renumbering for local and owned shared dofs
      std::vector<uint> send_buffer;
      pretabulated_dofmap_ = new uint[pretabulated_dofmap_size_];
      for (_set<uint>::iterator it = owned_dofs.begin();
      it != owned_dofs.end(); ++it, ++offset)
      {
        for(std::vector<uint>::iterator di = dof2index[*it].begin();
        di != dof2index[*it].end(); ++di)
        {
          pretabulated_dofmap_[*di] = offset;
        }

        if (shared_dofs.find(*it) != shared_dofs.end())
        {
          send_buffer.push_back(*it);
          send_buffer.push_back(offset);
        }
      }

      // Exchange new dof numbers for shared dofs
      MPI_Status status;
      int recv_count;
      uint pe_size = MPI::numProcesses();
      uint rank = MPI::processNumber();
      uint src, dest, max_recv;
      uint local_size = send_buffer.size();
      MPI_Allreduce(&local_size, &max_recv, 1, MPI_UNSIGNED, MPI_MAX,
                    MPI::DOLFIN_COMM);
      uint * recv_buffer = new uint[max_recv];
      for (uint k = 1; k < pe_size; ++k)
      {
        src = (rank - k + pe_size) % pe_size;
        dest = (rank + k) % pe_size;

        MPI_Sendrecv(&send_buffer[0], send_buffer.size(), MPI_UNSIGNED, dest, 1,
                     recv_buffer, max_recv, MPI_UNSIGNED, src, 1,
                     MPI::DOLFIN_COMM, &status);
        MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);

        for (int i = 0; i < recv_count; i += 2)
        {
          uint ghost_index = recv_buffer[i + 1];
          dolfin_assert(
              ghost_index < offset_ || ghost_index >= offset_ + local_size_);

          // Assign new dof number for shared dofs
          if (ghost_dofs.find(recv_buffer[i]) != ghost_dofs.end())
          {
            for (std::vector<uint>::iterator di =
                dof2index[recv_buffer[i]].begin();
                di != dof2index[recv_buffer[i]].end(); ++di)
            {
              pretabulated_dofmap_[*di] = ghost_index;
            }
            // Create ghost dofs list
            ghosts_.insert(ghost_index);
          }
        }
      }
      delete[] recv_buffer;

#if DEBUG
      uint const low = offset_;
      uint const high = offset_ + local_size_;
      for (_set<uint>::const_iterator git = ghosts_.begin();
      git != ghosts_.end(); ++git)
      {
        uint ii = *git;
        if( ii >= low && ii < high)
        {
          error("Ghost dof index %d found to be in ownership range."
          "Ownership range of dofs : %d, %d", ii, low, high);
        }
      }
#endif
      map_.clear();
    }

#if DEBUG
    if (type_ == real_space)
    {
      dolfin_assert(local_size_ == global_dimension());
    }
    else
    {
      // The sum of the local sizes should be the global size
      uint loc_s = local_size_;
      uint glb_s = 0;
      uint const expected_glob_s = this->global_dimension();
      MPI_Allreduce(&loc_s, &glb_s, 1, MPI_UNSIGNED, MPI_SUM, MPI::DOLFIN_COMM);
      if (glb_s != expected_glob_s)
      {
        error(
            "The sum of local dofmap sizes is not equal to the global dimension."
            "Sum: '%d' ; Global dimension : '%d'",
            glb_s, expected_glob_s);
      }
    }
#endif

#endif
  }
  else
  {
    type_ = ufc_default;
    offset_ = 0;
    local_size_ = global_dimension();
  }
}

//-----------------------------------------------------------------------------
std::map<uint, uint> DofMap::getMap() const
{
  return map_;
}

//-----------------------------------------------------------------------------
bool DofMap::is_ghost(uint i) const
{
  return (ghosts_.find(i) != ghosts_.end());
}

//-----------------------------------------------------------------------------
std::string const& DofMap::hash() const
{
  return hash_;
}

//-----------------------------------------------------------------------------
void DofMap::disp() const
{
  cout << "DofMap" << endl;
  cout << "------" << endl;

  // Begin indentation
  begin("");

  // Display UFC dofmap information
  cout << "ufc::dofmap info" << endl;
  cout << "-----------------" << endl;
  begin("");
  cout << "Signature            : " << ufc_dofmap_->signature() << endl;
  cout << "Global dimension     : " << ufc_dofmap_->global_dimension() << endl;
  cout << "Local dimension      : " << ufc_dofmap_->local_dimension() << endl;
  cout << "Geometric dimension  : " << ufc_dofmap_->geometric_dimension()
       << endl;
  cout << "Number of subdofmaps : " << ufc_dofmap_->num_sub_dofmaps() << endl;
  cout << "Number of facet dofs : " << ufc_dofmap_->num_facet_dofs() << endl;
  cout << endl;
  end();

#if DEBUG
  Mesh& dolfin_mesh = mesh();
  cout << "tabulate_dofs output" << endl;
  cout << "--------------------" << endl;
  begin("");
  {
    Mesh& dolfin_mesh = mesh();
    uint tdim = dolfin_mesh.topology().dim();
    uint num_dofs = ufc_dofmap_->local_dimension();
    uint* dofs = new uint[num_dofs];
    CellIterator cell(dolfin_mesh);
    UFCCell ufc_cell(*cell);
    for (; !cell.end(); ++cell)
    {
      ufc_cell.update(*cell, dolfin_mesh.distdata());

      ufc_dofmap_->tabulate_dofs(dofs, ufc_mesh_, ufc_cell);

      cout << "Cell " << ufc_cell.entity_indices[tdim][0] << ":  ";
      for (uint j = 0; j < num_dofs; ++j)
      {
        cout << dofs[j];
        if (j < num_dofs - 1) cout << ", ";
      }
      cout << endl;
    }
    delete[] dofs;
    cout << endl;
  }
  end();

  cout << "tabulate_coordinates output" << endl;
  cout << "---------------------------" << endl;
  begin("");
  {
    uint tdim = dolfin_mesh.topology().dim();
    uint gdim = ufc_dofmap_->geometric_dimension();
    uint num_dofs = ufc_dofmap_->local_dimension();
    double** coordinates = new double*[num_dofs];
    for (uint d = 0; d < num_dofs; ++d)
    {
      coordinates[d] = new double[gdim];
    }
    CellIterator cell(dolfin_mesh);
    UFCCell ufc_cell(*cell);
    for (; !cell.end(); ++cell)
    {
      ufc_cell.update(*cell, dolfin_mesh.distdata());

      ufc_dofmap_->tabulate_coordinates(coordinates, ufc_cell);

      cout << "Cell " << ufc_cell.entity_indices[tdim][0] << ":  ";
      for (uint j = 0; j < num_dofs; j++)
      {
        cout << "(";
        for (uint k = 0; k < gdim; ++k)
        {
          cout << coordinates[j][k];
          if (k < gdim - 1) cout << ", ";
        }
        cout << ")";
        if (j < num_dofs - 1) cout << ",  ";
      }
      cout << endl;
    }
    for (uint d = 0; d < num_dofs; ++d)
    {
      delete[] coordinates[d];
    }
    delete[] coordinates;
    cout << endl;
  }
  end();
#endif

  // End indentation
  end();
}
//-----------------------------------------------------------------------------
void DofMap::distributeByVote(UFCMesh& ufc_mesh, ufc::dofmap * ufc_dofmap,
                              _set<uint>& owned_dofs, _set<uint>& shared_dofs,
                              _set<uint>& ghost_dofs,
                              _map<uint, std::vector<uint> >& dof2index)
                              {
#ifdef HAVE_MPI
                              message("Distribute dofs by voting process.");

                              owned_dofs.clear();
                              shared_dofs.clear();
                              ghost_dofs.clear();
                              dof2index.clear();

                              //
                              Mesh& mesh = *const_cast<Mesh *>(ufc_mesh.mesh);
                              uint const local_dim = ufc_dofmap->local_dimension();
  uint * dofs = new uint[local_dim];
  uint const nb_facet_dofs = ufc_dofmap->num_facet_dofs();
  uint * facet_dofs = new uint[nb_facet_dofs];
  _map<uint, uint> dof_vote;

  Cell c_tmp(mesh, 1);
  UFCCell ufc_cell(c_tmp);

  // Initialize random number generator differently on each process
  // FIXME: if the ghosts are randomly generated for a given mesh then
  // two instances for the same mesh and same numbering have different
  // ghosts which leads to several issues:
  // - non matching ghosts with consequence of accessing the wrong dof
  // - different local size which leads to crash
  // We have to remove the randomness until a better solution is found.
  // srand((uint)time(0) + MPI::processNumber());
  srand(MPI::seed());

  // List shared dofs and assign vote
  std::vector<uint> send_buffer;
  BoundaryMesh& interior_boundary = mesh.interior_boundary();
  for (CellIterator bc(interior_boundary); !bc.end(); ++bc)
  {
    Facet f(mesh, interior_boundary.facet_index(*bc));
    Cell c(mesh, f.entities(mesh.topology().dim())[0]);

    uint local_facet = c.index(f);

    ufc_cell.update(c, mesh.distdata());
    ufc_dofmap->tabulate_dofs(dofs, ufc_mesh, ufc_cell);
    ufc_dofmap->tabulate_facet_dofs(facet_dofs, local_facet);

    for (uint i = 0; i < nb_facet_dofs; i++)
    {
      // Assign an ownership vote for each "shared" dof
      uint dofidx = dofs[facet_dofs[i]];
      if (shared_dofs.find(dofidx) == shared_dofs.end())
      {
        shared_dofs.insert(dofidx);
        dof_vote[dofidx] = (uint) (rand() + (real) MPI::processNumber());
        send_buffer.push_back(dofidx);
        send_buffer.push_back(dof_vote[dofidx]);
      }
    }
  }

  // Decide ownership of "shared" dofs
  MPI_Status status;
  int recv_count;
  uint pe_size = MPI::numProcesses();
  uint rank = MPI::processNumber();
  uint src, dest, max_recv;
  uint local_size = send_buffer.size();
  MPI_Allreduce(&local_size, &max_recv, 1, MPI_UNSIGNED, MPI_MAX,
                MPI::DOLFIN_COMM);
  uint *recv_buffer = new uint[max_recv];
  for (uint k = 1; k < MPI::numProcesses(); ++k)
  {
    src = (rank - k + pe_size) % pe_size;
    dest = (rank + k) % pe_size;
    MPI_Sendrecv(&send_buffer[0], send_buffer.size(), MPI_UNSIGNED, dest, 1,
                 recv_buffer, max_recv, MPI_UNSIGNED, src, 1, MPI::DOLFIN_COMM,
                 &status);
    MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);

    for (int i = 0; i < recv_count; i += 2)
    {
      if (shared_dofs.find(recv_buffer[i]) != shared_dofs.end())
      {
        // Move dofs with higher ownership votes from shared to forbidden
        if (recv_buffer[i + 1] < dof_vote[recv_buffer[i]]
            || (recv_buffer[i + 1] == dof_vote[recv_buffer[i]] && (src < rank)))
        {
          ghost_dofs.insert(recv_buffer[i]);
          shared_dofs.erase(recv_buffer[i]);
        }
      }
    }
  }
  send_buffer.clear();

  // Mark all non forbidden dofs as owned by the processes
  for (CellIterator c(mesh); !c.end(); ++c)
  {
    ufc_cell.update(*c, mesh.distdata());
    ufc_dofmap->tabulate_dofs(dofs, ufc_mesh, ufc_cell);
    for (uint i = 0; i < local_dim; ++i)
    {
      if (ghost_dofs.find(dofs[i]) == ghost_dofs.end())
      {
        // Mark dof as owned
        owned_dofs.insert(dofs[i]);
      }

      // Create mapping from dof to dofmap offset
      dof2index[dofs[i]].push_back(c->index() * local_dim + i);
    }
  }

  // Cleanup
  delete[] recv_buffer;
  delete[] facet_dofs;
  delete[] dofs;

#endif
}
//-----------------------------------------------------------------------------
void DofMap::distributeByEntities(UFCMesh& ufc_mesh, ufc::dofmap * ufc_dofmap,
                                  _set<uint>& owned_dofs,
                                  _set<uint>& shared_dofs,
                                  _set<uint>& ghost_dofs,
                                  _map<uint, std::vector<uint> >& dof2index)
                                  {
                                    message("Distribute dofs by mesh entities ownership.");

  owned_dofs.clear();
  shared_dofs.clear();
  ghost_dofs.clear();
  dof2index.clear();

  //
  Mesh& mesh = *const_cast<Mesh *>(ufc_mesh.mesh);
  uint const tdim = mesh.topology().dim();
  uint const local_dim = ufc_dofmap->local_dimension();
  uint * dofs = new uint[local_dim];
  CellIterator cell(mesh);
  UFCCell ufc_cell(*cell);

  // Compute facets and facet - cell connectivity if not already computed
  mesh.init(tdim - 1);
  mesh.init(tdim - 1, tdim);
  mesh.order();

  // Cache number of dofs per mesh entity
  Array<uint> num_entity_dofs;
  uint max_num_entity_dof = 0;
  for (uint d = 0; d <= mesh.topology().dim(); ++d)
  {
    uint const ned = ufc_dofmap->num_entity_dofs(d);
    num_entity_dofs.push_back(ned);
    max_num_entity_dof = std::max(max_num_entity_dof, ned);
  }
  dolfin_assert(max_num_entity_dof > 0);
  uint * entity_local_dofs = new uint[max_num_entity_dof];
  uint * entity_dofs = new uint[max_num_entity_dof];

  // Loop, baby, loop !
  for (; !cell.end(); ++cell)
  {
    ufc_cell.update(*cell, mesh.distdata());
    ufc_dofmap->tabulate_dofs(dofs, ufc_mesh, ufc_cell);

    // Create mapping from dof to dofmap offset
    for (uint i = 0; i < local_dim; ++i)
    {
      dof2index[dofs[i]].push_back(cell->index() * local_dim + i);
    }

    // Add dofs restricted to the cell as owned
    uint const num_celldofs = num_entity_dofs[tdim];
    ufc_dofmap->tabulate_entity_dofs(entity_local_dofs, tdim, 0);
    for (uint dof = 0; dof < num_celldofs; ++dof)
    {
      owned_dofs.insert(dofs[entity_local_dofs[dof]]);
    }

    // Decide ownership
    for (uint d = 0; d < mesh.topology().dim(); ++d)
    {
      uint const num_dofs = num_entity_dofs[d];
      for (MeshEntityIterator m(*cell, d); !m.end(); ++m)
      {
        // Get the dof indices for the entity
        ufc_dofmap->tabulate_entity_dofs(entity_local_dofs, d, m.pos());
        for (uint dof = 0; dof < num_dofs; ++dof)
        {
          entity_dofs[dof] = dofs[entity_local_dofs[dof]];
        }

        //
        if (mesh.distdata().is_shared(*m))
        {
          if (mesh.distdata().is_ghost(*m))
          {
            ghost_dofs.insert(entity_dofs, entity_dofs + num_dofs);
          }
          else
          {
            shared_dofs.insert(entity_dofs, entity_dofs + num_dofs);
            owned_dofs.insert(entity_dofs, entity_dofs + num_dofs);
          }
        }
        else
        {
          owned_dofs.insert(entity_dofs, entity_dofs + num_dofs);
        }
      }
    }
  }
  delete[] entity_dofs;
  delete[] entity_local_dofs;
  delete[] dofs;
}
//-----------------------------------------------------------------------------
    }
