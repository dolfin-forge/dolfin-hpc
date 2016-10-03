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
#include <dolfin/mesh/BoundaryMesh.h>
#include <dolfin/mesh/IntersectionDetector.h>
#include <dolfin/mesh/LocalMeshRefinement.h>
#include <dolfin/mesh/MappedManifold.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/MeshPartition.h>
#include <dolfin/mesh/MPIMeshCommunicator.h>
#include <dolfin/mesh/UniformMeshRefinement.h>
#include <dolfin/parameter/parameters.h>

#include <fstream>
#include <sstream>

namespace dolfin
{

//-----------------------------------------------------------------------------
Mesh::Mesh() :
    Variable("mesh", "DOLFIN mesh"),
    cell_type_(NULL),
    topology_(*this),
    geometry_(),
    exterior_boundary_(NULL),
    interior_boundary_(NULL),
    intersection_detector_(NULL),
    timestamp_(time(0))
{
  // Do nothing
}
//-----------------------------------------------------------------------------
Mesh::Mesh(Mesh const& mesh) :
    Variable("mesh", "DOLFIN mesh"),
    cell_type_(NULL),
    topology_(*this),
    geometry_(),
    exterior_boundary_(NULL),
    interior_boundary_(NULL),
    intersection_detector_(NULL),
    timestamp_(time(0))
{
  *this = mesh;
}
//-----------------------------------------------------------------------------
Mesh::Mesh(std::string const& filename) :
    Variable("mesh", "DOLFIN mesh"),
    cell_type_(NULL),
    topology_(*this),
    geometry_(),
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
  clear();
}
//-----------------------------------------------------------------------------
Mesh const& Mesh::operator=(Mesh const& other)
{
  clear();

  rename(other.name(), other.label());

  if (other.cell_type_)
  {
    cell_type_ = other.cell_type_->clone();
  }

  topology_ = other.topology_;
  geometry_ = other.geometry_;
  timestamp_ = other.timestamp_;

  for(Array<MappedManifold *>::iterator it = other.periodic_mappings_.begin();
      it != periodic_mappings_.end(); ++it)
  {
    PeriodicSubDomain const& p = (*it)->subdomain();
    periodic_mappings_ .push_back(new MappedManifold(*this, p));
  }

  return *this;
}
//-----------------------------------------------------------------------------
bool Mesh::operator ==(Mesh const& other) const
{
  return this->hash() == other.hash();
}
//-----------------------------------------------------------------------------
bool Mesh::operator !=(Mesh const& other) const
{
  return this->hash() != other.hash();
}
//-----------------------------------------------------------------------------
void Mesh::clear()
{
  timestamp_ = 0;
  delete cell_type_;
  cell_type_ = NULL;
  topology_.clear();
  geometry_.clear();
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
CellType& Mesh::type()
{
  dolfin_assert(cell_type_);
  return *cell_type_;
}
//-----------------------------------------------------------------------------
CellType const& Mesh::type() const
{
  dolfin_assert(cell_type_);
  return *cell_type_;
}
//-----------------------------------------------------------------------------
MeshTopology& Mesh::topology()
{
  return topology_;
}
//-----------------------------------------------------------------------------
MeshTopology const& Mesh::topology() const
{
  return topology_;
}
//-----------------------------------------------------------------------------
uint Mesh::size(uint dim) const
{
  return topology_.size(dim);
}
//-----------------------------------------------------------------------------
uint Mesh::num_cells() const
{
  return topology_.size(topology_.dim());
}
//-----------------------------------------------------------------------------
uint* Mesh::cells()
{
  return topology_(topology_.dim(), 0)();
}
//-----------------------------------------------------------------------------
uint const * Mesh::cells() const
{
  return topology_(topology_.dim(), 0)();
}
//-----------------------------------------------------------------------------
void Mesh::init(uint dim) const
{
  uint size = topology_.size(dim);
}
//-----------------------------------------------------------------------------
void Mesh::init(uint d0, uint d1) const
{
  uint size = topology_(d0, d1).size();
}
//-----------------------------------------------------------------------------
void Mesh::init() const
{
  for (uint d0 = 0; d0 <= topology().dim(); ++d0)
  {
    for (uint d1 = 0; d1 <= topology().dim(); ++d1)
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
  return (MPI::numProcesses() == 1) || dolfin_get("Mesh read in serial");
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
  return topology_.global_size(dim);
}
//-----------------------------------------------------------------------------
uint Mesh::num_global_cells() const
{
  return topology_.global_size(topology_.dim());
}
//-----------------------------------------------------------------------------
MeshGeometry& Mesh::geometry()
{
  return geometry_;
}
//-----------------------------------------------------------------------------
MeshGeometry const& Mesh::geometry() const
{
  return geometry_;
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
void Mesh::partition(MeshFunction<uint>& partitions)
{
  MeshPartition::partition(*this, partitions);
}
//-----------------------------------------------------------------------------
void Mesh::partition(MeshFunction<uint>& partitions, MeshFunction<uint>& weight)
{
  MeshPartition::partition(*this, partitions, weight);
}
//-----------------------------------------------------------------------------
void Mesh::partition_geom(MeshFunction<uint>& partitions)
{
  MeshPartition::partition_geom(*this, partitions);
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
    MeshFunction<uint> partitions;
    partition(partitions);
    distribute(partitions);
    //FIXME: following the legacy behaviour entities are always renumbered
    topology().renumber();
  }
}
//-----------------------------------------------------------------------------
void Mesh::distribute(MeshFunction<uint>& distribution)
{
  MPIMeshCommunicator::distribute(*this, distribution);
}
//-----------------------------------------------------------------------------
void Mesh::refine()
{
  message("No cells marked for refinement, assuming uniform mesh refinement.");
  UniformMeshRefinement::refine(*this);
}
//-----------------------------------------------------------------------------
void Mesh::refine(MeshFunction<bool>& cell_markers, bool refine_boundary,
                  bool load_balance)
{
  LocalMeshRefinement::refineMeshByEdgeBisection(*this, cell_markers,
                                                 refine_boundary, load_balance);
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
  ss << "Mesh@" << this << ":" << this->type().description() << ":time"
      << timestamp_ << ":T" << this->topology().token() << ":G"
      << this->geometry().token();
  return ss.str();
}
//-----------------------------------------------------------------------------
void Mesh::disp() const
{
  section("Mesh");
  topology_.disp();
  geometry_.disp();
  begin("Cell type");
  if (cell_type_)
  {
    cout << cell_type_->description() << endl;
  }
  else
  {
    cout << "undefined" << endl;
  }
  end();
  end();
}
//-----------------------------------------------------------------------------
void Mesh::check() const
{
  message("Check: mesh");
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
