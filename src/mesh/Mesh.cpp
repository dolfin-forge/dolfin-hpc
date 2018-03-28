// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Johan Hoffman, 2007.
// Modified by Garth N. Wells 2007.
// Modified by Balthasar Reuter, 2013.
// Modified by Aurélien Larcher, 2013.
//
// First added:  2006-05-09
// Last changed: 2013-03-22

#include <dolfin/mesh/Mesh.h>

#include <dolfin/io/File.h>
#include <dolfin/main/PE.h>
#include <dolfin/mesh/BoundaryMesh.h>
#include <dolfin/mesh/IntersectionDetector.h>
#include <dolfin/mesh/MappedManifold.h>
#include <dolfin/mesh/MeshData.h>
#include <dolfin/mesh/MeshPartition.h>
#include <dolfin/mesh/MPIMeshCommunicator.h>
#include <dolfin/mesh/Space.h>
#include <dolfin/mesh/UniformRefinement.h>
#include <dolfin/parameter/parameters.h>

#include <fstream>
#include <sstream>

namespace dolfin
{

#define DOLFIN_DEFAULT_MESH_NAME  "mesh"
#define DOLFIN_DEFAULT_MESH_LABEL "DOLFIN mesh"

//-----------------------------------------------------------------------------
Mesh::Mesh() :
    Variable(DOLFIN_DEFAULT_MESH_NAME, DOLFIN_DEFAULT_MESH_LABEL),
    topology_(NULL),
    geometry_(NULL),
    exterior_boundary_(NULL),
    interior_boundary_(NULL),
    intersection_detector_(NULL),
    timestamp_(time(0))
{
  // Do nothing
}
//-----------------------------------------------------------------------------
Mesh::Mesh(CellType const& ctype, Space const& space) :
    Variable(DOLFIN_DEFAULT_MESH_NAME, DOLFIN_DEFAULT_MESH_LABEL),
    topology_(new MeshTopology(ctype, DOLFIN_COMM_SELF, !this->reordering())),
    geometry_(new MeshGeometry(space)),
    exterior_boundary_(NULL),
    interior_boundary_(NULL),
    intersection_detector_(NULL),
    timestamp_(time(0))
{
}
//-----------------------------------------------------------------------------
Mesh::Mesh(CellType const& ctype, Space const& space, Comm& comm) :
    Variable(DOLFIN_DEFAULT_MESH_NAME, DOLFIN_DEFAULT_MESH_LABEL),
    topology_(new MeshTopology(ctype, comm, !this->reordering())),
    geometry_(new MeshGeometry(space)),
    exterior_boundary_(NULL),
    interior_boundary_(NULL),
    intersection_detector_(NULL),
    timestamp_(time(0))
{
}
//-----------------------------------------------------------------------------
Mesh::Mesh(Mesh const& other) :
    Variable(other.name(), other.label()),
    topology_(cloneptr(other.topology_)),
    geometry_(cloneptr(other.geometry_)),
    exterior_boundary_(copyptr(other.exterior_boundary_)),
    interior_boundary_(copyptr(other.interior_boundary_)),
    intersection_detector_(copyptr(other.intersection_detector_)),
    timestamp_(other.timestamp_)
{
  for(Array<MappedManifold *>::iterator it = other.periodic_mappings_.begin();
      it != other.periodic_mappings_.end(); ++it)
  {
    periodic_mappings_.push_back(new MappedManifold(*this, (*it)->subdomain()));
  }
}
//-----------------------------------------------------------------------------
Mesh::Mesh(std::string const& filename) :
    Variable(DOLFIN_DEFAULT_MESH_NAME, DOLFIN_DEFAULT_MESH_LABEL),
    topology_(NULL),
    geometry_(NULL),
    exterior_boundary_(NULL),
    interior_boundary_(NULL),
    intersection_detector_(NULL),
    timestamp_(time(0))
{
  File file(filename);
  file >> *this;
  this->distribute();
}
//-----------------------------------------------------------------------------
Mesh::~Mesh()
{
  timestamp_ = 0;
  delete topology_;
  topology_ = NULL;
  delete geometry_;
  geometry_ = NULL;
  delete exterior_boundary_;
  exterior_boundary_ = NULL;
  delete interior_boundary_;
  interior_boundary_ = NULL;
  delete intersection_detector_;
  intersection_detector_ = NULL;
  while(!periodic_mappings_.empty())
  {
    delete periodic_mappings_.back();
    periodic_mappings_.pop_back();
  }
}
//-----------------------------------------------------------------------------
void Mesh::swap(Mesh& other)
{
  if (this == &other) return;
  std::swap(topology_             , other.topology_);
  std::swap(geometry_             , other.geometry_);
  std::swap(exterior_boundary_    , other.exterior_boundary_);
  std::swap(interior_boundary_    , other.interior_boundary_);
  std::swap(intersection_detector_, other.intersection_detector_);
  std::swap(periodic_mappings_    , other.periodic_mappings_);
  std::swap(timestamp_            , other.timestamp_);
}
//-----------------------------------------------------------------------------
bool Mesh::operator ==(Mesh const& other) const
{
  if (!objptrcmp(topology_, other.topology_))
  {
    return false;
  }
  if (!objptrcmp(geometry_, other.geometry_))
  {
    return false;
  }
  return true;
}
//-----------------------------------------------------------------------------
bool Mesh::operator !=(Mesh const& other) const
{
  return !(*this == other);
}
//-----------------------------------------------------------------------------
bool Mesh::empty() const
{
  return (topology_ == NULL && geometry_ == NULL);
}
//-----------------------------------------------------------------------------
CellType const& Mesh::type() const
{
  dolfin_assert(topology_);
  return topology_->type();
}
//-----------------------------------------------------------------------------
MeshTopology& Mesh::topology()
{
  dolfin_assert(topology_);
  return *topology_;
}
//-----------------------------------------------------------------------------
MeshTopology const& Mesh::topology() const
{
  dolfin_assert(topology_);
  return *topology_;
}
//-----------------------------------------------------------------------------
uint Mesh::topology_dimension() const
{
  dolfin_assert(topology_);
  return topology_->dim();
}
//-----------------------------------------------------------------------------
uint Mesh::size(uint dim) const
{
  dolfin_assert(topology_);
  return topology_->size(dim);
}
//-----------------------------------------------------------------------------
uint Mesh::num_vertices() const
{
  dolfin_assert(topology_);
  return topology_->size(0);
}
//-----------------------------------------------------------------------------
uint Mesh::num_cells() const
{
  return topology_->size(topology_->dim());
}
//-----------------------------------------------------------------------------
uint* Mesh::cells()
{
  dolfin_assert(topology_);
  return (*topology_)(topology_->dim(), 0)();
}
//-----------------------------------------------------------------------------
uint const * Mesh::cells() const
{
  dolfin_assert(topology_);
  return (*topology_)(topology_->dim(), 0)();
}
//-----------------------------------------------------------------------------
void Mesh::init(uint dim) const
{
  dolfin_assert(topology_);
  topology_->size(dim);
}
//-----------------------------------------------------------------------------
void Mesh::init(uint d0, uint d1) const
{
  dolfin_assert(topology_);
  (*topology_)(d0, d1).order();
}
//-----------------------------------------------------------------------------
void Mesh::init() const
{
  for (uint d0 = 0; d0 <= topology_dimension(); ++d0)
  {
    for (uint d1 = 0; d1 <= topology_dimension(); ++d1)
    {
      init(d0, d1);
    }
  }
}
//-----------------------------------------------------------------------------
BoundaryMesh& Mesh::exterior_boundary()
{
  ///FIXME: Improve hash logic to regenerate boundary at topology change
  if (exterior_boundary_ == NULL || exterior_boundary_->invalid_mesh_topology())
  {
    if(exterior_boundary_)
    {
      warning("Recomputing mesh exterior boundary");
    }
    delete exterior_boundary_;
    exterior_boundary_ = new BoundaryMesh(*this, BoundaryMesh::exterior);
  }
  return *exterior_boundary_;
}
//-----------------------------------------------------------------------------
BoundaryMesh& Mesh::interior_boundary()
{
  ///FIXME: Improve hash logic to regenerate boundary at topology change
  if (interior_boundary_ == NULL || interior_boundary_->invalid_mesh_topology())
  {
    if(interior_boundary_)
    {
      warning("Recomputing mesh interior boundary");
    }
    delete interior_boundary_;
    interior_boundary_ = new BoundaryMesh(*this, BoundaryMesh::interior);
  }
  return *interior_boundary_;
}
//-----------------------------------------------------------------------------
bool Mesh::serial_io() const
{
  return (PE::size() == 1) || dolfin_get("Mesh read in serial");
}
//-----------------------------------------------------------------------------
bool Mesh::parallel_io() const
{
  return !this->serial_io();
}
//-----------------------------------------------------------------------------
bool Mesh::is_distributed() const
{
  return topology().is_distributed();
}
//-----------------------------------------------------------------------------
MeshDistributedData& Mesh::distdata()
{
  return topology().distdata();
}
//-----------------------------------------------------------------------------
MeshDistributedData const& Mesh::distdata() const
{
  return topology().distdata();
}
//-----------------------------------------------------------------------------
uint Mesh::global_size(uint dim) const
{
  dolfin_assert(topology_);
  return topology_->global_size(dim);
}
//-----------------------------------------------------------------------------
uint Mesh::num_global_vertices() const
{
  dolfin_assert(topology_);
  return topology_->global_size(0);
}
//-----------------------------------------------------------------------------
uint Mesh::num_global_cells() const
{
  dolfin_assert(topology_);
  return topology_->global_size(topology_->dim());
}
//-----------------------------------------------------------------------------
Space const& Mesh::space() const
{
  dolfin_assert(geometry_);
  return geometry_->space();
}
//-----------------------------------------------------------------------------
MeshGeometry& Mesh::geometry()
{
  dolfin_assert(geometry_);
  return *geometry_;
}
//-----------------------------------------------------------------------------
MeshGeometry const& Mesh::geometry() const
{
  dolfin_assert(geometry_);
  return *geometry_;
}
//-----------------------------------------------------------------------------
uint Mesh::geometry_dimension() const
{
  dolfin_assert(geometry_);
  return geometry_->dim();
}
//-----------------------------------------------------------------------------
IntersectionDetector& Mesh::intersector()
{
  ///FIXME: Improve hash logic to regenerate detector at topology change
  if (intersection_detector_ == NULL)
  {
    if(intersection_detector_)
    {
      warning("Recreating mesh intersection detector");
    }
    delete intersection_detector_;
    intersection_detector_ = new IntersectionDetector(*this);
  }
  return *intersection_detector_;
}
//-----------------------------------------------------------------------------
void Mesh::partition(MeshValues<uint, Cell>& partitions)
{
  MeshPartition::partition(partitions);
}
//-----------------------------------------------------------------------------
void Mesh::partition(MeshValues<uint, Cell>& partitions, MeshValues<uint, Cell>& weight)
{
  MeshPartition::partition(partitions, weight);
}
//-----------------------------------------------------------------------------
void Mesh::partition_geom(MeshValues<uint, Vertex>& partitions)
{
  MeshPartition::partition_geom(partitions);
}
//-----------------------------------------------------------------------------
void Mesh::distribute()
{
  if (this->parallel_io())
  {
    // COMMENT: At this point the distributed data cannot be empty as the file
    //          format is supposed to fill it.
    if(!topology().is_distributed())
    {
      error("The topology of a mesh read in parallel should be distributed.");
    }
    MeshValues<uint, Cell> partitions(*this);
    partition(partitions);
    distribute(partitions);
    //FIXME: following the legacy behaviour entities are always renumbered
    topology().renumber();
  }
}
//-----------------------------------------------------------------------------
void Mesh::distribute(MeshValues<uint, Cell>& distribution)
{
  MPIMeshCommunicator::distribute(distribution);
}
//-----------------------------------------------------------------------------
void Mesh::distribute(MeshValues<uint, Vertex>& distribution)
{
  MPIMeshCommunicator::distribute(distribution);
}
//-----------------------------------------------------------------------------
void Mesh::distribute(MeshValues<uint, Cell>& distribution, MeshData& data)
{
  MPIMeshCommunicator::distribute(distribution, &data);
}
//-----------------------------------------------------------------------------
void Mesh::refine()
{
  UniformRefinement R; R(*this);
}
//-----------------------------------------------------------------------------
bool Mesh::has_periodic_constraint() const
{
  return (!periodic_mappings_.empty());
}
//-----------------------------------------------------------------------------
void Mesh::add_periodic_constraint(PeriodicSubDomain const& periodic)
{
  periodic_mappings_.push_back(new MappedManifold(*this, periodic));
}
//-----------------------------------------------------------------------------
Array<MappedManifold *> const& Mesh::periodic_mappings() const
{
  for(Array<MappedManifold *>::iterator it = periodic_mappings_.begin();
      it != periodic_mappings_.end(); ++it)
  {
    if((*it)->invalid_mesh())
    {
      PeriodicSubDomain const * p = &(*it)->subdomain();
      delete (*it);
      warning("Recreating mesh periodic mapping");
      Mesh& mesh = const_cast<Mesh&>(*this);
      (*it) = new MappedManifold(mesh, *p);
    }
  }
  return periodic_mappings_;
}
//-----------------------------------------------------------------------------
std::string const Mesh::hash() const
{
  std::stringstream ss;
  ss << "Mesh@" << this << ":"
      << (topology_ ? topology_->type().description() : "empty")
      << ":time" << timestamp_
      << ":T" << (topology_ ? topology_->token() : 0)
      << ":G" << geometry_->token();
  return ss.str();
}
//-----------------------------------------------------------------------------
void Mesh::disp() const
{
  section("Mesh");
  if (this->empty())
  {
    message("Empty");
  }
  else
  {
    section("Topology");
    message("distributed : %u", this->is_distributed());
    if(this->is_distributed())
    {
      message("cells    : local = %12u ; global = %12u",
              this->num_cells(), this->num_global_cells());
      message("vertices : local = %12u ; global = %12u",
              this->num_vertices(), this->num_global_vertices());
    }
    else
    {
      message("cells    : %12u", this->num_cells());
      message("vertices : %12u", this->num_vertices());
    }
    end();
  }
  end();
}
//-----------------------------------------------------------------------------
void Mesh::check() const
{
  MPIMeshCommunicator::check(const_cast<Mesh&>(*this));
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
