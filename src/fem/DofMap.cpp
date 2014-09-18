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
  message(0, "Extracted dof map for subspace: %s", ufc_dofmap_->signature());

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
  message(0, "Extracted dof map for sub system: %s", ufc_dofmap_->signature());
  message(0, "Offset for sub system: %d", offset);

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
ufc::dofmap* DofMap::create_sub_dofmap(Array<uint> const& sub_system,
                                       uint& offset) const
{
  // Reset offset
  offset = 0;

  // Recursively extract sub dof map
  ufc::dofmap* sub_dofmap = DofMap::create_sub_dofmap(ufc_mesh_, *ufc_dofmap_,
                                                      sub_system, offset);
  message(0, "Extracted ufc dof map for sub system: %s",
          sub_dofmap->signature());
  message(0, "Offset for sub system: %d", offset);

  return sub_dofmap;
}

//-----------------------------------------------------------------------------
ufc::dofmap* DofMap::create_sub_dofmap(UFCMesh& ufc_mesh,
                                       ufc::dofmap const& dofmap,
                                       Array<uint> const& sub_system,
                                       uint& offset)
{
  // Check if there are any sub systems
  if (dofmap.num_sub_dofmaps() == 0)
  {
    error("Unable to extract sub system (there are no sub systems).");
  }

  // Check that a sub system has been specified
  if (sub_system.size() == 0)
  {
    error("Unable to extract sub system (no sub system specified).");
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
    DofMap::init_ufc(ufc_mesh, *ufc_sub_dofmap);
    offset += ufc_sub_dofmap->global_dimension();
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
                                                          offset);
  delete sub_dofmap;

  return sub_sub_dofmap;
}

//-----------------------------------------------------------------------------
void DofMap::init_ufc(UFCMesh& ufc_mesh, ufc::dofmap& ufc_dofmap)
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
  DofMap::init_ufc(ufc_mesh_, *ufc_dofmap_);

  // Build the DOLFIN dofmap
  build();

  // Set offsets and local dimensions for mixed elements
  uint const nb_sub = this->num_sub_dofmaps();
  if (nb_sub > 0)
  {
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

  pretabulated_dofmap_size_ = this->local_dimension() * mesh().numCells();
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
    case scalar_p1:
      {
        for (uint i = 0; i < local_dimension(); ++i)
        {
          dofs[i] = ufc_cell.entity_indices[0][i];
          dolfin_assert(mesh().distdata().have_global(dofs[i], 0));
        }
      }
      break;
    case scalar_dg0:
      {
        dofs[0] = ufc_cell.index;
        dolfin_assert(mesh().distdata().have_global(dofs[0], mesh().topology().dim()));
      }
      break;
    case vector_p1:
      {
        uint const num_cellverts = cell.numEntities(0);
        for (uint k = 0; k < ufc_dofmap_->num_sub_dofmaps(); ++k)
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
        for (uint i = 0; i < ufc_dofmap_->num_sub_dofmaps(); ++i)
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
                     Array<ufc::dofmap const *>& stack) const
{
  // Single root dofmap
  if (dofmap->num_sub_dofmaps() == 0)
  {
    stack.push_back(dofmap->create());
    return;
  }
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
      flatten(sub, stack);
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
    pretabulate_all_dofs();
  }
  return pretabulated_dofmap_;
}

//-----------------------------------------------------------------------------
uint DofMap::dofsmapping_size() const
{
  return pretabulated_dofmap_size_;
}

//--------------------------------------------------------------------------
void DofMap::pretabulate_all_dofs() const
{
  CellIterator ref_cell(mesh());
  UFCCell ufc_cell(*ref_cell);
  uint const local_dim = this->local_dimension();

  pretabulated_dofmap_ = new uint[pretabulated_dofmap_size_];
  uint *ip = &pretabulated_dofmap_[0];
  MeshDistributedData& distdata = mesh().distdata();
  for (CellIterator cell(mesh()); !cell.end(); ++cell)
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

  delete[] pretabulated_dofmap_;
  pretabulated_dofmap_ = NULL;

  map.clear();

  if (mesh().is_distributed())
  {
#ifdef HAVE_MPI
    Mesh& dolfin_mesh = mesh();
    uint *dofs = new uint[local_dimension()];

    uint pe_size = MPI::numProcesses();
    uint rank = MPI::processNumber();

    //
    dolfin_mesh.renumber();

    if (ufc_dofmap_->global_dimension()
        == dolfin_mesh.distdata().global_numVertices())
    {
      // Scalar Lagrange P1
      type_ = scalar_p1;
      local_size_ = dolfin_mesh.numVertices()
          - dolfin_mesh.distdata().num_ghost(0);
    }
    else if (ufc_dofmap_->global_dimension()
        == dolfin_mesh.distdata().global_numCells())
    {
      // Scalar Discontinuous Lagrange P0
      type_ = scalar_dg0;
      local_size_ = dolfin_mesh.numCells();
    }
    else if (ufc_dofmap_->global_dimension()
        == ufc_dofmap_->num_sub_dofmaps()
            * dolfin_mesh.distdata().global_numVertices())
    {
      // Vector Lagrange P1
      type_ = vector_p1;
      uint gdim = ufc_dofmap_->num_sub_dofmaps();
      uint num_local = dolfin_mesh.numVertices()
          - dolfin_mesh.distdata().num_ghost(0);

      uint num_dofs = gdim * num_local;
      uint offset = 0;

#if ( MPI_VERSION > 1 )
      MPI_Exscan(&num_dofs, &offset, 1, MPI_UNSIGNED, MPI_SUM,
                 MPI::DOLFIN_COMM);
#else
      MPI_Scan(&num_dofs, &offset, 1, MPI_UNSIGNED, MPI_SUM, MPI::DOLFIN_COMM);
      offset -= num_dofs;
#endif
      _map<uint, uint> v_offset;

      for (VertexIterator v(dolfin_mesh); !v.end(); ++v)
      {
        if (!dolfin_mesh.distdata().is_ghost(v->index(), 0))
        {
          v_offset[dolfin_mesh.distdata().get_global(*v)] = offset;
          offset += gdim;
        }
      }

      Array<uint> *ghost_buff = new Array<uint> [pe_size];
      for (MeshGhostIterator iter(dolfin_mesh.distdata(), 0); !iter.end();
          ++iter)
        ghost_buff[iter.owner()].push_back(
            dolfin_mesh.distdata().get_vertex_global(iter.index()));

      MPI_Status status;
      Array<uint> send_buff;
      uint src, dest;
      uint recv_size = dolfin_mesh.distdata().num_ghost(0);
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

      vertex_map_ = new uint[dolfin_mesh.numVertices()];
      for (VertexIterator v(dolfin_mesh); !v.end(); ++v)
        vertex_map_[v->index()] =
            v_offset[dolfin_mesh.distdata().get_vertex_global(v->index())];

      v_offset.clear();

      local_size_ = gdim * num_local;

      for (uint i = 0; i < pe_size; i++)
        ghost_buff[i].clear();
      delete[] ghost_buff;

    }
    else if (ufc_dofmap_->global_dimension()
        == ufc_dofmap_->num_sub_dofmaps()
            * dolfin_mesh.distdata().global_numCells())
    {
      // Vector Discontinuous Lagrange P0
      type_ = vector_dg0;
      uint gdim = ufc_dofmap_->num_sub_dofmaps();
      uint num_dofs = gdim * dolfin_mesh.numCells();
      uint offset = 0;

#if ( MPI_VERSION > 1 )
      MPI_Exscan(&num_dofs, &offset, 1, MPI_UNSIGNED, MPI_SUM,
                 MPI::DOLFIN_COMM);
#else
      MPI_Scan(&num_dofs, &offset, 1, MPI_UNSIGNED, MPI_SUM, MPI::DOLFIN_COMM);
      offset -= num_dofs;
#endif
      offset_ = offset;
      local_size_ = gdim * dolfin_mesh.numCells();
    }
    else
    {
      type_ = generic;
      BoundaryMesh interior_boundary(dolfin_mesh, BoundaryMesh::interior);
      MeshFunction<uint>* cell_map = interior_boundary.data().meshFunction(
          "cell map");

      std::vector<uint> send_buffer;
      _set<uint> shared_dofs, forbidden_dofs, owned_dofs;
      _map<uint, uint> dof_vote;
      _map<uint, std::vector<uint> > dof2index;

      uint n = local_dimension();
      pretabulated_dofmap_ = new uint[n * mesh().numCells()];
      uint *facet_dofs = new uint[num_facet_dofs()];

      Cell c_tmp(dolfin_mesh, 1);
      UFCCell ufc_cell(c_tmp);

      // Initialize random number generator differently on each process
      // FIXME: if the ghosts are randomly generated for a given mesh then
      // two instances for the same mesh and same numbering have different
      // ghosts which leads to several issues:
      // - non matching ghosts with consequence of accessing the wrong dof
      // - different local size which leads to crash
      // We have to remove the randomness until a better solution is found.
      // srand((uint)time(0) + MPI::processNumber());
      srand(MPI::processNumber());

      // Decide ownership of shared dofs
      for (CellIterator bc(interior_boundary); !bc.end(); ++bc)
      {
        Facet f(dolfin_mesh, cell_map->get(*bc));
        Cell c(dolfin_mesh, f.entities(dolfin_mesh.topology().dim())[0]);

        uint local_facet = c.index(f);

        ufc_cell.update(c, dolfin_mesh.distdata());
        ufc_dofmap_->tabulate_dofs(dofs, ufc_mesh_, ufc_cell);
        ufc_dofmap_->tabulate_facet_dofs(facet_dofs, local_facet);

        for (uint i = 0; i < num_facet_dofs(); i++)
        {
          // Assign an ownership vote for each "shared" dof
          uint dofidx = dofs[facet_dofs[i]];
          if (shared_dofs.find(dofidx) == shared_dofs.end())
          {
            shared_dofs.insert(dofidx);
            dof_vote[dofidx] = (uint) rand() + (uint) MPI::processNumber();
            send_buffer.push_back(dofidx);
            send_buffer.push_back(dof_vote[dofidx]);
          }
        }
      }

      // Decide ownership of "shared" dofs
      MPI_Status status;
      int recv_count;
      uint src, dest, max_recv;
      uint num_proc = MPI::numProcesses();
      uint proc_num = MPI::processNumber();
      uint local_size = send_buffer.size();
      MPI_Allreduce(&local_size, &max_recv, 1, MPI_UNSIGNED, MPI_MAX,
                    MPI::DOLFIN_COMM);
      uint *recv_buffer = new uint[max_recv];
      for (uint k = 1; k < MPI::numProcesses(); ++k)
      {
        src = (proc_num - k + num_proc) % num_proc;
        dest = (proc_num + k) % num_proc;
        MPI_Sendrecv(&send_buffer[0], send_buffer.size(), MPI_UNSIGNED, dest, 1,
                     recv_buffer, max_recv, MPI_UNSIGNED, src, 1,
                     MPI::DOLFIN_COMM, &status);
        MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);

        for (int i = 0; i < recv_count; i += 2)
        {
          if (shared_dofs.find(recv_buffer[i]) != shared_dofs.end())
          {
            // Move dofs with higher ownership votes from shared to forbidden
            if (recv_buffer[i + 1] < dof_vote[recv_buffer[i]])
            {
              forbidden_dofs.insert(recv_buffer[i]);
              shared_dofs.erase(recv_buffer[i]);
            }
          }
        }
      }

      send_buffer.clear();

      // Mark all non forbidden dofs as owned by the processes
      for (CellIterator c(dolfin_mesh); !c.end(); ++c)
      {
        ufc_cell.update(*c, dolfin_mesh.distdata());
        ufc_dofmap_->tabulate_dofs(dofs, ufc_mesh_, ufc_cell);
        for (uint i = 0; i < n; i++)
        {
          if (forbidden_dofs.find(dofs[i]) == forbidden_dofs.end())
          {
            // Mark dof as owned
            owned_dofs.insert(dofs[i]);
          }

          // Create mapping from dof to dofmap offset
          dof2index[dofs[i]].push_back(c->index() * n + i);
        }
      }

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

      // Compute renumbering for local and owned shared dofs
      for (_set<uint>::iterator it = owned_dofs.begin();
      it != owned_dofs.end(); ++it, offset++)
      {
        for(std::vector<uint>::iterator di = dof2index[*it].begin();
        di != dof2index[*it].end(); ++di)
        pretabulated_dofmap_[*di] = offset;

        if (shared_dofs.find(*it) != shared_dofs.end())
        {
          send_buffer.push_back(*it);
          send_buffer.push_back(offset);
        }
      }

      // Exchange new dof numbers for shared dofs
      delete[] recv_buffer;
      local_size = send_buffer.size();
      MPI_Allreduce(&local_size, &max_recv, 1, MPI_UNSIGNED, MPI_MAX,
                    MPI::DOLFIN_COMM);
      recv_buffer = new uint[max_recv];
      for (uint k = 1; k < MPI::numProcesses(); ++k)
      {
        src = (proc_num - k + num_proc) % num_proc;
        dest = (proc_num + k) % num_proc;

        MPI_Sendrecv(&send_buffer[0], send_buffer.size(), MPI_UNSIGNED, dest, 1,
                     recv_buffer, max_recv, MPI_UNSIGNED, src, 1,
                     MPI::DOLFIN_COMM, &status);
        MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);

        for (int i = 0; i < recv_count; i += 2)
        {
          // Assign new dof number for shared dofs
          if (forbidden_dofs.find(recv_buffer[i]) != forbidden_dofs.end())
          {
            for (std::vector<uint>::iterator di =
                dof2index[recv_buffer[i]].begin();
                di != dof2index[recv_buffer[i]].end(); ++di)
              pretabulated_dofmap_[*di] = recv_buffer[i + 1];
          }
        }
      }
      delete[] recv_buffer;
      map.clear();
      delete[] facet_dofs;
    }
    delete[] dofs;

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
std::map<dolfin::uint, dolfin::uint> DofMap::getMap() const
{
  return map;
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
        if (j < num_dofs - 1)
        cout << ", ";
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
          if (k < gdim - 1)
          cout << ", ";
        }
        cout << ")";
        if (j < num_dofs - 1)
        cout << ",  ";
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

}
