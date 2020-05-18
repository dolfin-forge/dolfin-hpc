// Copyright (C) 2007-2008 Anders Logg and Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/fem/DofMap.h>

#include <dolfin/config/dolfin_config.h>
#include <dolfin/common/types.h>
#include <dolfin/fem/DofNumbering.h>
#include <dolfin/fem/DofMapCache.h>
#include <dolfin/fem/PeriodicDofsMapping.h>
#include <dolfin/fem/UFCCell.h>
#include <dolfin/mesh/BoundaryMesh.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/main/MPI.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
DofMap::DofMap(Mesh& mesh, ufc::form const& form, uint const i) :
    MeshDependent(mesh),
    offset_(0),
    ufc_dofmap_(form.create_dofmap(i)),
    numbering_(DofNumbering::create(mesh, *ufc_dofmap_)),
    hash_(make_hash(mesh, *ufc_dofmap_)),
    periodic_dofmap_(NULL)
{
  init();
}

//-----------------------------------------------------------------------------
DofMap::DofMap(Mesh& mesh, ufc::dofmap& dofmap, bool const owner) :
    MeshDependent(mesh),
    offset_(0),
    ufc_dofmap_((owner ? &dofmap : dofmap.create())),
    numbering_(DofNumbering::create(mesh, *ufc_dofmap_)),
    hash_(make_hash(mesh, *ufc_dofmap_)),
    periodic_dofmap_(NULL)
{
  init();
}

//-----------------------------------------------------------------------------
DofMap::DofMap(DofMap const& dofmap, uint i) :
    MeshDependent(dofmap.mesh()),
    offset_(0),
    ufc_dofmap_(dofmap.create_sub_dofmap(i)),
    numbering_(DofNumbering::create(dofmap.mesh(), *ufc_dofmap_)),
    hash_(make_hash(mesh(), *ufc_dofmap_)),
    periodic_dofmap_(NULL)
{
  message(1, "DofMap: Extracted dof map for subspace: %s",
          ufc_dofmap_->signature());

  init();
}

//-----------------------------------------------------------------------------
DofMap::DofMap(DofMap const& dofmap, Array<uint> const& subsystem, uint& offset) :
    MeshDependent(dofmap.mesh()),
    offset_(0),
    ufc_dofmap_(dofmap.create_sub_dofmap(subsystem, offset_)),
    numbering_(DofNumbering::create(dofmap.mesh(), *ufc_dofmap_)),
    hash_(make_hash(mesh(), *ufc_dofmap_)),
    periodic_dofmap_(NULL)
{
  // Check that dof map has not be re-ordered
  offset = offset_;
  message(1, "DofMap: Extracted dof map for sub system: %s",
          ufc_dofmap_->signature());
  message(1, "DofMap: Offset for sub system: %d", offset);

  // Reset offset
  offset_ = 0;
  init();
}

//-----------------------------------------------------------------------------
DofMap::~DofMap()
{
  delete periodic_dofmap_;
  flattened_.free();
  delete numbering_;
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
  message(1, "DofMap: Extracted ufc dof map for sub system: %s",
          sub_dofmap->signature());
  message(1, "DofMap: Local offset for sub system: %d", local_offset);

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
  message(1, "DofMap: Extracted ufc dof map for sub system: %s",
          sub_dofmap->signature());
  message(1, "DofMap: Local offset for sub system: %d", local_offset);

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
void DofMap::init()
{
  // Build the DOLFIN dofmap
  message(1, "DofMap: init dofmap for signature:\n %s", this->signature());
  numbering_->build();
  message(1, "DofMap: offset = %u; size = %u", numbering_->offset(),
          numbering_->size());

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
void DofMap::tabulate_dofs(uint* dofs, ufc::cell const& ufc_cell, Cell const& cell) const
{
  dolfin_assert(dofs != NULL);
  numbering_->tabulate_dofs(dofs, ufc_cell, cell);
}

//-----------------------------------------------------------------------------
void DofMap::tabulate_dofs(uint* dofs, UFCCell const& ufc_cell) const
{
  dolfin_assert(dofs != NULL);
  numbering_->tabulate_dofs(dofs, ufc_cell, *ufc_cell);
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
                     Array<ufc::dofmap const *>& stack, uint maxlevel)
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
      DofMap::flatten(sub, stack, maxlevel - 1);
    }
  }
}

//-----------------------------------------------------------------------------
void DofMap::flatten(ufc::dofmap const * dofmap,
                     Array<ufc::dofmap const *>& stack)
{
  // Single root element or max level is set to zero, return immediately
  if (dofmap->num_sub_dofmaps() == 0)
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
      DofMap::flatten(sub, stack);
    }
  }
}

//-----------------------------------------------------------------------------
bool DofMap::is_vectorizable() const
{
  return DofMap::can_vectorize(this->flatten());
}

//-----------------------------------------------------------------------------
bool DofMap::can_vectorize(Array<ufc::dofmap const *> flattened)
{
  bool ret = true;
  for (uint s = 1; s < flattened.size(); ++s)
  {
    if (std::strcmp(flattened[0]->signature(), flattened[s]->signature()) != 0)
    {
      ret = false;
      break;
    }
  }
  return ret;
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
  // FIXME:
  return (true);
}

//-----------------------------------------------------------------------------
uint DofMap::local_size() const
{
  return numbering_->size();
}

//-----------------------------------------------------------------------------
uint const * DofMap::dofsmapping() const
{
  return numbering_->block();
}

//-----------------------------------------------------------------------------
uint DofMap::dofsmapping_size() const
{
  return numbering_->block_size();
}

//--------------------------------------------------------------------------
PeriodicDofsMapping const& DofMap::periodic_mapping() const
{
  if (periodic_dofmap_ == NULL)
  {
    periodic_dofmap_ = new PeriodicDofsMapping(*this);
  }
  return *periodic_dofmap_;
}

//-----------------------------------------------------------------------------
_ordered_map<uint, uint> DofMap::get_map() const
{
  return map_;
}

//-----------------------------------------------------------------------------
bool DofMap::is_shared(uint index) const
{
  return numbering_->is_shared(index);
}

//-----------------------------------------------------------------------------
bool DofMap::is_ghost(uint index) const
{
  return numbering_->is_ghost(index);
}

//-----------------------------------------------------------------------------
std::string const& DofMap::hash() const
{
  return hash_;
}

//-----------------------------------------------------------------------------
void DofMap::disp() const
{
  section("DofMap");
  begin("ufc::dofmap info");
  prm("Signature"           , ufc_dofmap_->signature());
  prm("Global dimension"    , ufc_dofmap_->global_dimension());
  prm("Local dimension"     , ufc_dofmap_->local_dimension());
  prm("Geometric dimension" , ufc_dofmap_->geometric_dimension());
  prm("Number of subdofmaps", ufc_dofmap_->num_sub_dofmaps());
  prm("Number of facet dofs", ufc_dofmap_->num_facet_dofs());
  end();
  end();
}
//-----------------------------------------------------------------------------
bool DofMap::check(bool throw_error)
{
  bool ret = true;

  message("Check dofs distribution");
  message("signature   : %s", this->signature());
  bool distributed_by_entities_ = true;
  message("by entities : %d", distributed_by_entities_);
  //FIXME:
  // Check if a dof is owned twice, send to everyone to make sure
  BoundaryMesh& boundary = mesh().interior_boundary();
  Mesh& mesh = this->mesh();
  MeshDistributedData& distdata = mesh.distdata();
  uint const tdim = mesh.topology_dimension();

  if (this->num_facet_dofs() == 0)
  {
    return true;
  }

  //
  _ordered_set<uint> shared_owned;
  typedef _ordered_map<uint, std::pair<uint, uint> > EntitiesDofMap;
  EntitiesDofMap shared_owned_entities;
  Cell c0(mesh, 0);
  UFCCell ufc_cell(c0);
  uint * cell_dofs = new uint[this->local_dimension()];
  uint * num_entity_dofs = new uint[tdim];
  for (uint d = 0; d < tdim; ++d)
  {
    mesh.init(tdim, d);
    mesh.init(mesh.type().facet_dim(), d);
    num_entity_dofs[d] = ufc_dofmap_->num_entity_dofs(d);
  }
  uint const num_facet_dofs = this->num_facet_dofs();
  uint * loc_entity_dofs = new uint[this->num_facet_dofs()];
  for (CellIterator bcell(boundary); !bcell.end(); ++bcell)
  {
    Facet f(boundary.mesh(), boundary.facet_index(*bcell));
    dolfin_assert(f.is_shared());

    // Tabulate cell dofs
    Cell cell(mesh, f.entities(tdim)[0]);
    ufc_cell.update(cell);
    this->tabulate_dofs(cell_dofs, ufc_cell);

    // Check that all shared facet dofs are shared
    uint const local_facet = cell.index(f);
    ufc_dofmap_->tabulate_facet_dofs(loc_entity_dofs, local_facet);
    _ordered_set<uint> shared_facet_dofs;
    _ordered_set<uint> local_facet_dofs;
    for (uint dof = 0; dof < num_facet_dofs; ++dof)
    {
      uint gdof = cell_dofs[loc_entity_dofs[dof]];
      if (!is_shared(gdof))
      {
        error("Dof with local_index %u on shared facet is not shared",
              loc_entity_dofs[dof]);
      }
      shared_facet_dofs.insert(gdof);
      local_facet_dofs.insert(loc_entity_dofs[dof]);
    }

    // Get boundary facet dofs
    uint const facet_dim = mesh.type().facet_dim();
    ufc_dofmap_->tabulate_entity_dofs(loc_entity_dofs, facet_dim, local_facet);
    for (uint dof = 0; dof < num_entity_dofs[facet_dim]; ++dof)
    {
      dolfin_assert(is_shared(cell_dofs[loc_entity_dofs[dof]]));
      if (!is_ghost(cell_dofs[loc_entity_dofs[dof]]))
      {
        if (distributed_by_entities_)
        {
          if (distdata[f.dim()].is_ghost(f.index()))
          {
            error("Non-ghosted dof's facet is ghosted");
          }
          shared_owned_entities[cell_dofs[loc_entity_dofs[dof]]] = std::pair<
              uint, uint>(distdata[f.dim()].get_global(f.index()), f.dim());
        }
        shared_owned.insert(cell_dofs[loc_entity_dofs[dof]]);
      }
    }

    // Check lower dimensional entities of boundary facet
    for (uint d = 0; d < (facet_dim); ++d)
    {
      for (MeshEntityIterator m(f, d); !m.end(); ++m)
      {
        dolfin_assert(m->is_shared());
        // Get the dof indices for the entity
        uint local_index = cell.index(*m);
        MeshEntity e(mesh, d, cell.entities(d)[local_index]);
        if (!distdata[e.dim()].is_shared(e.index()))
        {
          error("Entity of dim %u in shared facet is not shared", d);
        }
        ufc_dofmap_->tabulate_entity_dofs(loc_entity_dofs, d, local_index);

        for (uint dof = 0; dof < num_entity_dofs[d]; ++dof)
        {
          if (local_facet_dofs.count(loc_entity_dofs[dof]) == 0)
          {
            error("Dof index is not a facet dof:\n"
                  "Entity of dim %u with local index %u",
                  d, cell.index(*m));
          }
          uint gdof = cell_dofs[loc_entity_dofs[dof]];
          if (shared_facet_dofs.count(gdof) == 0)
          {
            error("Dof index is not a shared facet dof:\n"
                  "Entity of dim %u with local index %u",
                  d, cell.index(*m));
          }
          dolfin_assert(is_shared(gdof));
          if (!is_ghost(gdof))
          {
            if (distributed_by_entities_)
            {
              if (distdata[m->dim()].is_ghost(m->index()))
              {
                error("Non-ghosted dof's entity of dim %u is ghosted", d);
              }
              shared_owned_entities[gdof] = std::pair<uint, uint>(
                  distdata[m->dim()].get_global(m->index()), m->dim());
            }
            shared_owned.insert(gdof);
          }
        }
      }
    }
  }
  delete[] num_entity_dofs;
  delete[] loc_entity_dofs;
  delete[] cell_dofs;

#if HAVE_MPI

  // Simple version
  uint rank = dolfin::MPI::rank();
  uint pe_size = dolfin::MPI::size();

  Array<uint> sendbuf;
  if (distributed_by_entities_)
  {
    for (EntitiesDofMap::const_iterator it = shared_owned_entities.begin();
        it != shared_owned_entities.end(); ++it)
    {
      sendbuf.push_back(it->second.first);
      sendbuf.push_back(it->second.second);
      sendbuf.push_back(it->first);
    }
  }
  else
  {
#ifdef __SUNPRO_CC
    for (_ordered_set<uint>::iterator it = shared_owned.begin();
        it != shared_owned.end(); ++it)
    {
      sendbuf.push_back(*it);
    }
#else
    sendbuf.insert(sendbuf.end(), shared_owned.begin(), shared_owned.end());
#endif
  }

  int recv_count = sendbuf.size();
  int maxrev_count = 0;
  MPI::all_reduce<MPI::max>( recv_count, maxrev_count );
  dolfin_assert(maxrev_count > 0);
  uint * recvbuf = new uint[maxrev_count];

  _set<uint> owned_more;
  for (uint p = 1; p < pe_size; ++p)
  {
    int src = (rank - p + pe_size) % pe_size;
    int dest = (rank + p) % pe_size;

    //
    recv_count = MPI::sendrecv( &sendbuf[0], sendbuf.size(), dest,
                                recvbuf, maxrev_count, src, 0);
    for (uint i = 0; i < uint(recv_count);)
    {
      if (distributed_by_entities_)
      {
        uint const idx = recvbuf[i++];
        uint const dim = recvbuf[i++];
        if (distdata[dim].has_global(idx))
        {
          uint loc_id = distdata[dim].get_local(idx);
          if (!distdata[dim].is_ghost(loc_id))
          {
            if (throw_error)
            {
              error("Entity on rank %u is owned by rank %u", rank, src);
            }
            else
            {
              warning("Entity on rank %u is owned by rank %u", rank, src);
            }
          }
        }
      }
      uint const dof = recvbuf[i++];
      if (is_shared(dof) && !is_ghost(dof))
      {
        owned_more.insert(dof);
        if (throw_error)
        {
          error("Degree of freedom on rank %u is also owned by rank %u", rank,
                src);
        }
        else
        {
          warning("Degree of freedom on rank %u is also owned by rank %u", rank,
                  src);
        }
      }
    }
  }
  ret &= owned_more.empty();
  delete[] recvbuf;
#else
  MAYBE_UNUSED(throw_error);
#endif
  return ret;
}

//-----------------------------------------------------------------------------
std::string const DofMap::make_signature(std::string const& finite_element)
{
  return "FFC dofmap for " + finite_element;
}

//-----------------------------------------------------------------------------
std::string const DofMap::make_hash(Mesh& mesh, ufc::dofmap const& ufc_dofmap)
{
  return std::string(ufc_dofmap.signature()) + mesh.hash();
}

//-----------------------------------------------------------------------------

}
