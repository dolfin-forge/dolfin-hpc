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
#include <dolfin/mesh/UniformRefinement.h>
#include <dolfin/mesh/celltypes/TetrahedronCell.h>
#include <dolfin/mesh/entities/iterators/VertexIterator.h>

#include <sstream>

namespace dolfin
{

#define DOLFIN_DEFAULT_MESH_NAME "mesh"
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

auto swap( Mesh & a, Mesh & b ) -> void
{
  using std::swap;

  swap( a.topology_, b.topology_ );
  swap( a.geometry_, b.geometry_ );
  swap( a.exterior_boundary_, b.exterior_boundary_ );
  swap( a.interior_boundary_, b.interior_boundary_ );
  swap( a.intersection_detector_, b.intersection_detector_ );
  swap( a.periodic_mappings_, b.periodic_mappings_ );
  swap( a.extent_min_, b.extent_min_ );
  swap( a.extent_max_, b.extent_max_ );
  swap( a.timestamp_, b.timestamp_ );
}

//-----------------------------------------------------------------------------

auto Mesh::operator=( Mesh const & other ) -> Mesh const &
{
  Mesh tmp( other );
  swap( *this, tmp );

  return *this;
}

//-----------------------------------------------------------------------------

auto Mesh::operator==( Mesh const & other ) const -> bool
{
  if ( not( topology_ == other.topology_ ) )
  {
    return false;
  }

  if ( not( geometry_ == other.geometry_ ) )
  {
    return false;
  }

  if ( not( exterior_boundary_ == other.exterior_boundary_ ) )
  {
    return false;
  }

  if ( not( interior_boundary_ == other.interior_boundary_ ) )
  {
    return false;
  }

  if ( not( intersection_detector_ == other.intersection_detector_ ) )
  {
    return false;
  }

  return true;
}

//-----------------------------------------------------------------------------

auto Mesh::extent_update() -> void
{
  Point min_( +DOLFIN_REAL_MAX, +DOLFIN_REAL_MAX, +DOLFIN_REAL_MAX );
  Point max_( -DOLFIN_REAL_MAX, -DOLFIN_REAL_MAX, -DOLFIN_REAL_MAX );

  for ( VertexIterator v( *this ); not v.end(); ++v )
  {
    real const * x = v->x();
    for ( size_t i = 0; i < geometry_.dim(); ++i )
    {
      min_[i] = std::min( min_[i], x[i] );
      max_[i] = std::max( max_[i], x[i] );
    }
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

auto Mesh::init() const -> void
{
  for ( size_t d0 = 0; d0 <= topology_dimension(); ++d0 )
  {
    for ( size_t d1 = 0; d1 <= topology_dimension(); ++d1 )
    {
      init( d0, d1 );
    }
  }
}

//-----------------------------------------------------------------------------

auto Mesh::exterior_boundary() -> BoundaryMesh &
{
  /// @todo Improve hash logic to regenerate boundary at topology change
  if ( exterior_boundary_ == nullptr
       || exterior_boundary_->invalid_mesh_topology() )
  {
    if ( exterior_boundary_ )
    {
      warning( "Recomputing mesh exterior boundary" );
    }
    delete exterior_boundary_;
    exterior_boundary_ = new BoundaryMesh( *this, BoundaryMesh::exterior );
  }
  return *exterior_boundary_;
}

//-----------------------------------------------------------------------------

auto Mesh::interior_boundary() -> BoundaryMesh &
{
  /// @todo Improve hash logic to regenerate boundary at topology change
  if ( interior_boundary_ == nullptr
       || interior_boundary_->invalid_mesh_topology() )
  {
    if ( interior_boundary_ )
    {
      warning( "Recomputing mesh interior boundary" );
    }
    delete interior_boundary_;
    interior_boundary_ = new BoundaryMesh( *this, BoundaryMesh::interior );
  }
  return *interior_boundary_;
}

//-----------------------------------------------------------------------------

auto Mesh::intersector() -> IntersectionDetector &
{
  /// @todo Improve hash logic to regenerate detector at topology change
  if ( intersection_detector_ == nullptr )
  {
    if ( intersection_detector_ )
    {
      warning( "Recreating mesh intersection detector" );
    }
    delete intersection_detector_;
    intersection_detector_ = new IntersectionDetector( *this );
  }
  return *intersection_detector_;
}

//-----------------------------------------------------------------------------

auto Mesh::partition( MeshValues< size_t, Cell > & partitions ) -> void
{
  MeshPartition::partition( partitions );
}

//-----------------------------------------------------------------------------

auto Mesh::partition( MeshValues< size_t, Cell > & partitions,
                      MeshValues< size_t, Cell > & weight ) -> void
{
  MeshPartition::partition( partitions, weight );
}

//-----------------------------------------------------------------------------

auto Mesh::partition_geom( MeshValues< size_t, Vertex > & partitions ) -> void
{
  MeshPartition::partition_geom( partitions );
}

//-----------------------------------------------------------------------------

auto Mesh::distribute() -> void
{
  if ( this->parallel_io() )
  {
    // COMMENT: At this point the distributed data cannot be empty as the file
    //          format is supposed to fill it.
    if ( !topology().distributed() )
    {
      error( "The topology of a mesh read in parallel should be distributed." );
    }
    MeshValues< size_t, Cell > partitions( *this );
    partition( partitions );
    distribute( partitions );
    /// @todo following the legacy behaviour entities are always renumbered
    topology().renumber();
  }
}

//-----------------------------------------------------------------------------

auto Mesh::distribute( MeshValues< size_t, Cell > & distribution ) -> void
{
  MPIMeshCommunicator::distribute( distribution );
}

//-----------------------------------------------------------------------------

auto Mesh::distribute( MeshValues< size_t, Vertex > & distribution ) -> void
{
  MPIMeshCommunicator::distribute( distribution );
}

//-----------------------------------------------------------------------------

auto Mesh::distribute( MeshValues< size_t, Cell > & distribution,
                       MeshData &                   data ) -> void
{
  MPIMeshCommunicator::distribute( distribution, &data );
}

//-----------------------------------------------------------------------------

auto Mesh::refine() -> void
{
  UniformRefinement::refine( *this );
}

//-----------------------------------------------------------------------------

auto Mesh::has_periodic_constraint() const -> bool
{
  return ( !periodic_mappings_.empty() );
}

//-----------------------------------------------------------------------------

auto Mesh::add_periodic_constraint( PeriodicSubDomain const & periodic ) -> void
{
  periodic_mappings_.push_back( new MappedManifold( *this, periodic ) );
}

//-----------------------------------------------------------------------------

auto Mesh::periodic_mappings() const -> std::vector< MappedManifold * > const &
{
  for ( std::vector< MappedManifold * >::iterator it =
          periodic_mappings_.begin();
        it != periodic_mappings_.end();
        ++it )
  {
    if ( ( *it )->invalid_mesh() )
    {
      PeriodicSubDomain const * p = &( *it )->subdomain();
      delete ( *it );
      warning( "Recreating mesh periodic mapping" );
      Mesh & mesh = const_cast< Mesh & >( *this );
      ( *it )     = new MappedManifold( mesh, *p );
    }
  }
  return periodic_mappings_;
}

//-----------------------------------------------------------------------------

auto Mesh::hash() const -> std::string const
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

auto Mesh::disp() const -> void
{
  section( "Mesh" );
  section( "Topology" );
  message( "distributed : %u", topology_.distributed() );
  if ( this->is_distributed() )
  {
    message( "cells    : local = %12u ; global = %12u",
             this->num_cells(),
             this->num_global_cells() );
    message( "vertices : local = %12u ; global = %12u",
             this->num_vertices(),
             this->num_global_vertices() );
  }
  else
  {
    message( "cells    : %12u", this->num_cells() );
    message( "vertices : %12u", this->num_vertices() );
  }
  end();
  section( "Extent" );
  message(
    "Min: [%f, %f, %f]", extent_min_[0], extent_min_[1], extent_min_[2] );
  message(
    "Max: [%f, %f, %f]", extent_max_[0], extent_max_[1], extent_max_[2] );
  end();
}

//-----------------------------------------------------------------------------

auto Mesh::check() const -> void
{
  MPIMeshCommunicator::check( const_cast< Mesh & >( *this ) );
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */
