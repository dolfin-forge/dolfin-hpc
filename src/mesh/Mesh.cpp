// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/mesh/Mesh.h>

#include <dolfin/io/File.h>
#include <dolfin/mesh/BoundaryMesh.h>
#include <dolfin/mesh/EuclideanSpace.h>
#include <dolfin/mesh/IntersectionDetector.h>
#include <dolfin/mesh/MPIMeshCommunicator.h>
#include <dolfin/mesh/MappedManifold.h>
#include <dolfin/mesh/MeshData.h>
#include <dolfin/mesh/MeshPartition.h>
#include <dolfin/mesh/Space.h>
#include <dolfin/mesh/TetrahedronCell.h>
#include <dolfin/mesh/UniformRefinement.h>
#include <dolfin/mesh/VertexIterator.h>

#include <sstream>

namespace dolfin
{

#define DOLFIN_DEFAULT_MESH_NAME  "mesh"
#define DOLFIN_DEFAULT_MESH_LABEL "DOLFIN mesh"

//-----------------------------------------------------------------------------

Mesh::Mesh()
  : Mesh( TetrahedronCell(), EuclideanSpace( 3 ), DOLFIN_COMM )
{
}

//-----------------------------------------------------------------------------

Mesh::Mesh( CellType const & ctype, Space const & space, Comm & comm )
  : Variable( DOLFIN_DEFAULT_MESH_NAME, DOLFIN_DEFAULT_MESH_LABEL )
  , topology_( ctype, comm, !this->reordering() )
  , geometry_( space )
  , exterior_boundary_( nullptr )
  , interior_boundary_( nullptr )
  , intersection_detector_( nullptr )
  , extent_min_( DOLFIN_REAL_MAX, DOLFIN_REAL_MAX, DOLFIN_REAL_MAX )
  , extent_max_( DOLFIN_REAL_MIN, DOLFIN_REAL_MIN, DOLFIN_REAL_MIN )
  , timestamp_( time( nullptr ) )
{
}

//-----------------------------------------------------------------------------


Mesh::Mesh( Mesh const & other )
  : Variable( other.name(), other.label() )
  , topology_( other.topology_ )
  , geometry_( other.geometry_ )
  , exterior_boundary_( copyptr( other.exterior_boundary_ ) )
  , interior_boundary_( copyptr( other.interior_boundary_ ) )
  , intersection_detector_( copyptr( other.intersection_detector_ ) )
  , extent_min_( other.extent_min_ )
  , extent_max_( other.extent_max_ )
  , timestamp_( other.timestamp_ )
{
  for ( MappedManifold * m : other.periodic_mappings_ )
  {
    periodic_mappings_.push_back( new MappedManifold( *this, m->subdomain() ) );
  }
}

//-----------------------------------------------------------------------------

Mesh::Mesh( std::string const & filename )
  : Mesh()
{
  File( filename ) >> *this;
  this->distribute();
}

//-----------------------------------------------------------------------------

Mesh::~Mesh()
{
  timestamp_ = 0;

  delete exterior_boundary_;
  exterior_boundary_ = nullptr;

  delete interior_boundary_;
  interior_boundary_ = nullptr;

  delete intersection_detector_;
  intersection_detector_ = nullptr;

  destruct( periodic_mappings_ );
}

//-----------------------------------------------------------------------------

void swap( Mesh& a, Mesh& b )
{
  using std::swap;

  swap( a.topology_,              b.topology_              );
  swap( a.geometry_,              b.geometry_              );
  swap( a.exterior_boundary_,     b.exterior_boundary_     );
  swap( a.interior_boundary_,     b.interior_boundary_     );
  swap( a.intersection_detector_, b.intersection_detector_ );
  swap( a.periodic_mappings_,     b.periodic_mappings_     );
  swap( a.extent_min_,            b.extent_min_            );
  swap( a.extent_max_,            b.extent_max_            );
  swap( a.timestamp_,             b.timestamp_             );
}

//-----------------------------------------------------------------------------

Mesh const & Mesh::operator=( Mesh const & other )
{
  Mesh tmp(other);
  swap( *this, tmp );

  return *this;
}

//-----------------------------------------------------------------------------

bool Mesh::operator ==(Mesh const& other) const
{
  if ( not ( topology_ == other.topology_ ) )
  {
    return false;
  }

  if ( not ( geometry_ == other.geometry_ ) )
  {
    return false;
  }

  if ( not ( exterior_boundary_ == other.exterior_boundary_ ) )
  {
    return false;
  }

  if ( not ( interior_boundary_ == other.interior_boundary_ ) )
  {
    return false;
  }

  if ( not ( intersection_detector_ == other.intersection_detector_ ) )
  {
    return false;
  }

  return true;
}

//-----------------------------------------------------------------------------

void Mesh::extent_update()
{
  Point min_( +DOLFIN_REAL_MAX, +DOLFIN_REAL_MAX, +DOLFIN_REAL_MAX );
  Point max_( -DOLFIN_REAL_MAX, -DOLFIN_REAL_MAX, -DOLFIN_REAL_MAX );

  for ( VertexIterator v( *this ); not v.end(); ++v )
  {
    Point const & p = v->point();
    min_[0] = std::min( min_[0], p[0] );
    min_[1] = std::min( min_[1], p[1] );
    min_[2] = std::min( min_[2], p[2] );

    max_[0] = std::max( max_[0], p[0] );
    max_[1] = std::max( max_[1], p[1] );
    max_[2] = std::max( max_[2], p[2] );
  }

#ifdef HAVE_MPI
  MPI::Communicator & comm = topology().comm();
  MPI::all_reduce< MPI::min >( &min_[0], &extent_min_[0], 3, comm );
  MPI::all_reduce< MPI::max >( &max_[0], &extent_max_[0], 3, comm );
#else
  extent_min_ = min_;
  extent_max_ = max_;
#endif
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
  /// @todo Improve hash logic to regenerate boundary at topology change
  if (exterior_boundary_ == nullptr || exterior_boundary_->invalid_mesh_topology())
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
  /// @todo Improve hash logic to regenerate boundary at topology change
  if (interior_boundary_ == nullptr || interior_boundary_->invalid_mesh_topology())
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


IntersectionDetector& Mesh::intersector()
{
  /// @todo Improve hash logic to regenerate detector at topology change
  if (intersection_detector_ == nullptr)
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
    if(!topology().distributed())
    {
      error("The topology of a mesh read in parallel should be distributed.");
    }
    MeshValues<uint, Cell> partitions(*this);
    partition(partitions);
    distribute(partitions);
    /// @todo following the legacy behaviour entities are always renumbered
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
  UniformRefinement::refine(*this);
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
  ss << "Mesh@" << this << ":";
  ss << topology_.type().description();
  ss << ":time" << timestamp_;
  ss << ":T" << topology_.token();
  ss << ":G" << geometry_.token();
  return ss.str();
}

//-----------------------------------------------------------------------------

void Mesh::disp() const
{
  section( "Mesh" );
  section( "Topology" );
  message( "distributed : %u", topology_.distributed() );
  if ( this->is_distributed() )
  {
    message( "cells    : local = %12u ; global = %12u",
             this->num_cells(), this->num_global_cells() );
    message( "vertices : local = %12u ; global = %12u",
             this->num_vertices(), this->num_global_vertices() );
  }
  else
  {
    message( "cells    : %12u", this->num_cells() );
    message( "vertices : %12u", this->num_vertices() );
  }
  end();
  section( "Extent" );
  message( "Min: [%f, %f, %f]", extent_min_[0], extent_min_[1], extent_min_[2] );
  message( "Max: [%f, %f, %f]", extent_max_[0], extent_max_[1], extent_max_[2] );
  end();
}

//-----------------------------------------------------------------------------

void Mesh::check() const
{
  MPIMeshCommunicator::check(const_cast<Mesh&>(*this));
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */
