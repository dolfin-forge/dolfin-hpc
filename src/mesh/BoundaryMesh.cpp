// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/mesh/BoundaryMesh.h>

#include <dolfin/log/log.h>
#include <dolfin/main/PE.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/SubDomain.h>
#include <dolfin/mesh/entities/Facet.h>
#include <dolfin/mesh/entities/iterators/CellIterator.h>
#include <dolfin/mesh/entities/iterators/FacetIterator.h>
#include <dolfin/mesh/entities/iterators/VertexIterator.h>

#include <iostream>

namespace dolfin
{

//-----------------------------------------------------------------------------

BoundaryMesh::BoundaryMesh( Mesh & mesh, BoundaryMesh::Type type )
  : Mesh( mesh.type(), mesh.space(), mesh.topology().comm() )
  , MeshDependent( static_cast< Mesh & >( mesh ) )
  , type_( type )
  , boundary_of_boundary_( false )
  , cell_map_()
  , vertex_map_()
  , subdomain_( nullptr )
{
  init( mesh, type );
}

//-----------------------------------------------------------------------------

BoundaryMesh::BoundaryMesh( BoundaryMesh & mesh, BoundaryMesh::Type type )
  : Mesh( mesh.type(), mesh.space(), mesh.topology().comm() )
  , MeshDependent( static_cast< Mesh & >( mesh ) )
  , type_( type )
  , boundary_of_boundary_( true )
  , cell_map_()
  , vertex_map_()
  , subdomain_( nullptr )
{
  init( mesh, type );
}

//-----------------------------------------------------------------------------

BoundaryMesh::BoundaryMesh( Mesh &             mesh,
                            SubDomain const &  subdomain,
                            BoundaryMesh::Type type )
  : Mesh( mesh.type(), mesh.space(), mesh.topology().comm() )
  , MeshDependent( mesh )
  , type_( type )
  , boundary_of_boundary_( false )
  , cell_map_()
  , vertex_map_()
  , subdomain_( &subdomain )
{
  init( mesh, type );
}

//-----------------------------------------------------------------------------

BoundaryMesh::BoundaryMesh( BoundaryMesh const & other )
  : Mesh( other )
  , MeshDependent( other.mesh() )
  , type_( other.type_ )
  , boundary_of_boundary_( other.boundary_of_boundary_ )
  , cell_map_( other.cell_map_ )
  , vertex_map_( other.vertex_map_ )
  , subdomain_( other.subdomain_ )
{
}

//-----------------------------------------------------------------------------

void BoundaryMesh::init( Mesh & mesh, BoundaryMesh::Type type )
{
  switch ( type )
  {
    case BoundaryMesh::exterior:
      // Exterior boundary i.e facets at the domain boundary
      compute( mesh, true, false );
      break;
    case BoundaryMesh::interior:
      // Interior boundary i.e facets between processors
      compute( mesh, false, true );
      break;
    case BoundaryMesh::full:
      // Full boundary including facets between processors
      compute( mesh, true, true );
      break;
    default:
      error( "Unknown boundary mesh type." );
      break;
  }
}

//-----------------------------------------------------------------------------

BoundaryMesh::BoundaryMesh( BoundaryMesh &    boundary,
                            SubDomain const & subdomain,
                            bool              inside )
  : Mesh( boundary.type(), boundary.space(), boundary.topology().comm() )
  , MeshDependent( boundary.mesh() )
  , type_( boundary.type_ )
  , boundary_of_boundary_( true )
  , cell_map_()
  , vertex_map_()
  , subdomain_( &subdomain )
{
  Mesh &       mesh = boundary.mesh();
  size_t const gdim = mesh.geometry_dimension();
  size_t const tdim = mesh.topology_dimension();

  if ( tdim == 1 )
  {
    vertex_map_.clear();
    for ( VertexIterator v( boundary ); not v.end(); ++v )
    {
      if ( subdomain_->inside( v->x(), not v->is_shared() ) )
      {
        vertex_map_.push_back( boundary.vertex_map_[v->index()] );
      }
    }
    cell_map_ = vertex_map_;

    // Create boundary vertices and cells
    MeshEditor editor( *this, mesh.type().facetType(), gdim );
    editor.init_vertices( vertex_map_.size() );
    editor.init_cells( cell_map_.size() );
    MeshGeometry const & geom = mesh.geometry();
    for ( size_t i = 0; i < vertex_map_.size(); ++i )
    {
      editor.add_vertex( i, geom.x( vertex_map_[i] ) );
      editor.add_cell( i, &vertex_map_[i] );
    }
    // If the mesh is distributed, set global numbering and copy the ownership
    if ( mesh.is_distributed() )
    {
      this->distdata()[0].assign( mesh.distdata()[0], vertex_map_ );
    }
    editor.close();
  }
  else
  {

    size_t const          num_verts = boundary.vertex_map_.size();
    std::vector< size_t > boundary_vertices( num_verts, num_verts );
    for ( CellIterator c( boundary ); not c.end(); ++c )
    {
      bool const on_boundary = not c->is_shared();
      if ( ( inside and subdomain_->enclosed( *c, on_boundary ) )
           or ( not inside and subdomain_->overlap( *c, on_boundary ) ) )
      {
        cell_map_.push_back( c->index() );
        for ( VertexIterator v( *c ); not v.end(); ++v )
        {
          size_t const vertex_index = v->index();
          if ( boundary_vertices[v->index()] == num_verts )
          {
            boundary_vertices[vertex_index] = vertex_map_.size();
            vertex_map_.push_back( boundary.vertex_map_[vertex_index] );
          }
        }
      }
    }

    // Create boundary vertices and cells
    MeshEditor editor( *this, mesh.type().facetType(), gdim );
    editor.init_vertices( vertex_map_.size() );
    for ( size_t i = 0; i < vertex_map_.size(); ++i )
    {
      editor.add_vertex( i, mesh.geometry().x( vertex_map_[i] ) );
    }
    editor.init_cells( cell_map_.size() );
    //
    size_t const          d                  = mesh.type().facet_dim();
    size_t const          num_facet_vertices = mesh.type().num_vertices( d );
    std::vector< size_t > facet_vertices( num_facet_vertices );
    for ( size_t i = 0; i < cell_map_.size(); ++i )
    {
      Cell c( boundary, cell_map_[i] );
      cell_map_[i] = boundary.cell_map_[cell_map_[i]];
      for ( VertexIterator v( c ); not v.end(); ++v )
      {
        facet_vertices[v.pos()] = boundary_vertices[v->index()];
      }
      editor.add_cell( i, facet_vertices.data() );
    }

    // If the mesh is distributed, set global numbering and copy the ownership
    if ( mesh.is_distributed() )
    {
      this->distdata()[0].assign( mesh.distdata()[0], vertex_map_ );
      this->distdata()[d].assign( mesh.distdata()[d], cell_map_ );
    }
    editor.close();
  }
}

//-----------------------------------------------------------------------------

BoundaryMesh::~BoundaryMesh()
{
  // Do nothing
}

//-----------------------------------------------------------------------------

void BoundaryMesh::compute( Mesh & mesh, bool exterior, bool interior )
{
  size_t const gdim = mesh.geometry_dimension();
  size_t const tdim = mesh.topology_dimension();

  // If the boundary is full then no need to compute the facet map
  bool const full = interior and exterior;

  message( 1, "BoundaryMesh : compute %s boundary",
              ( full ? "full" : ( interior ? "interior" : "exterior" ) ) );

  if ( tdim <= 1 )
  {
    vertex_map_.clear();
    for ( VertexIterator v( mesh ); not v.end(); ++v )
    {
      // Boundary facets are connected to exactly one cell
      bool const bndr = ( v->num_entities( tdim ) == 1 );
      // Interior facets are shared while exterior facets are not
      bool const shrd = v->is_shared();
      if ( bndr
           and ( full or ( interior and shrd ) or ( exterior and not shrd ) ) )
      {
        if ( ( subdomain_ == nullptr )
             or ( subdomain_->inside( v->x(), bndr and not shrd ) ) )
        {
          vertex_map_.push_back( v->index() );
        }
      }
    }
    cell_map_ = vertex_map_;

    // Create boundary vertices and cells
    MeshEditor editor( *this, mesh.type().facetType(), gdim );
    editor.init_vertices( vertex_map_.size() );
    editor.init_cells( cell_map_.size() );
    MeshGeometry const & geom = mesh.geometry();
    for ( size_t i = 0; i < vertex_map_.size(); ++i )
    {
      editor.add_vertex( i, geom.x( vertex_map_[i] ) );
      editor.add_cell( i, &vertex_map_[i] );
    }
    // If the mesh is distributed, set global numbering and copy the ownership
    if ( mesh.is_distributed() )
    {
      this->distdata()[0].assign( mesh.distdata()[0], vertex_map_ );
    }
    editor.close();
  }
  else
  {
    size_t const pe_size   = PE::size();
    size_t const num_verts = mesh.size( 0 );
    cell_map_.clear();
    vertex_map_.clear();

    std::vector< std::vector< size_t > > shared_vertices( pe_size );
    std::vector< size_t > boundary_vertices( num_verts, num_verts );

    for ( FacetIterator f( mesh ); not f.end(); ++f )
    {
      // Boundary facets are connected to exactly one cell
      bool const bndr = ( f->num_entities( tdim ) == 1 );
      // Interior facets are shared while exterior facets are not
      bool const shrd = f->is_shared();

#if DEBUG
      if ( shrd and not f->has_all_vertices_shared() )
      {
        error( "Shared facets have not all vertices shared" );
      }
#endif
      if ( bndr
           and ( full or ( interior and shrd ) or ( exterior and not shrd ) ) )
      {
        if ( subdomain_ != nullptr )
        {
          bool inside = false;
          for ( VertexIterator v( *f ); not v.end(); ++v )
          {
            if ( subdomain_->inside( v->x(), bndr and not shrd ) )
            {
              inside = true;
            }
          }
          if ( not inside )
          {
            continue;
          }
        }
        cell_map_.push_back( f->index() );
        for ( VertexIterator v( *f ); not v.end(); ++v )
        {
          size_t const vertex_index = v->index();
          if ( boundary_vertices[vertex_index] == num_verts )
          {
            boundary_vertices[vertex_index] = vertex_map_.size();
            vertex_map_.push_back( vertex_index );
            // Collect shared vertices per adjacent ranks
            if ( v->is_shared() )
            {
              for ( size_t const adj : mesh.distdata()[0].get_shared_adj( vertex_index ) )
              {
                dolfin_assert( mesh.distdata()[0].get_adj_ranks().count( adj ) );
                shared_vertices[adj].push_back( v->global_index() );
              }
            }
          }
        }
      }
    }

    // The mesh distribution is not constrained: non-owner of exterior facets
    // may own lower-dimensional entities on the exterior boundary.
    if ( mesh.is_distributed() )
    {
#if HAVE_MPI
      DistributedData const & distdata = mesh.distdata()[0];
      // _set< size_t > const &  vadjs    = distdata.get_adj_ranks();
      size_t                  recvmax  = max_array_size( shared_vertices );
      MPI::all_reduce_in_place< MPI::max >( recvmax );
      std::vector< std::vector< size_t > > recvbuf( pe_size,
                                                    std::vector< size_t >( recvmax ) );

      std::vector< MPI_Request > sendreq( pe_size );
      std::vector< MPI_Request > recvreq( pe_size );
      std::vector< MPI_Status >  status( pe_size );

      _set< size_t > added_vertices( vertex_map_.begin(), vertex_map_.end() );

      for ( size_t i = 0; i < pe_size; ++i )
      {
        MPI::check_error( MPI_Irecv( recvbuf[i].data(), recvbuf[i].size(),
                                     MPI_type< size_t >::value,
                                     i, 0, mesh.topology().comm(),
                                     &recvreq[i] ) );
      }

      // FIXME safety
      MPI_Barrier( mesh.topology().comm() );

      for ( size_t i = 0; i < pe_size; ++i )
      {
        MPI::check_error( MPI_Isend( shared_vertices[i].data(),
                                     shared_vertices[i].size(),
                                     MPI_type< size_t >::value,
                                     i, 0, mesh.topology().comm(),
                                     &sendreq[i] ) );
      }

      // FIXME safety
      MPI_Barrier( mesh.topology().comm() );

      MPI::check_error( MPI_Waitall( pe_size, sendreq.data(), status.data() ) );


      for ( size_t i = 0; i < pe_size; ++i )
      {
        int recvcount = 0;
        MPI::check_error( MPI_Wait( &recvreq[i], &status[i] ) );
        MPI::check_error( MPI_Get_count( &status[i], MPI_type< size_t >::value,
                                         &recvcount ) );

        for ( int k = 0; k < recvcount; ++k )
        {
          dolfin_assert( distdata.has_global( recvbuf[i][k] ) );
          size_t const local_index = distdata.get_local( recvbuf[i][k] );
          if ( added_vertices.count( local_index ) == 0 )
          {
            vertex_map_.push_back( local_index );
            added_vertices.insert( local_index );
          }
        }
      }

      // // If one adjacent rank has no boundary cell but the current rank has
      // // ghost entities owned by this rank then any algorithm based on mesh
      // // facets is bound to fail.
      // // Naive implementation for now, just send all the ghost vertices
      // for ( size_t const adj : vadjs )
      // {
      //   MPI::check_error( MPI_Send( &shared_vertices[adj][0],
      //                               shared_vertices[adj].size(),
      //                               MPI_type< size_t >::value,
      //                               adj, 0, MPI::DOLFIN_COMM ) );
      // }

      // for ( size_t const adj : vadjs )
      // {
      //   MPI::check_error( MPI_Recv( recvbuf.data(), recvbuf.size(),
      //                               MPI_type< size_t >::value,
      //                               adj, 0, MPI::DOLFIN_COMM, &status ) );
      //   MPI::check_error( MPI_Get_count( &status, MPI_type< size_t >::value,
      //                                    &recvcount ) );

      //   for ( int k = 0; k < recvcount; ++k )
      //   {
      //     dolfin_assert( distdata.has_global( recvbuf[k] ) );
      //     size_t const local_index = distdata.get_local( recvbuf[k] );
      //     if ( added_vertices.count( local_index ) == 0 )
      //     {
      //       vertex_map_.push_back( local_index );
      //       added_vertices.insert( local_index );
      //     }
      //   }
      // }

#endif /* HAVE_MPI */
    }

    // Create boundary vertices and cells
    MeshEditor editor( *this, mesh.type().facetType(), gdim );
    editor.init_vertices( vertex_map_.size() );
    for ( size_t i = 0; i < vertex_map_.size(); ++i )
    {
      editor.add_vertex( i, mesh.geometry().x( vertex_map_[i] ) );
    }
    editor.init_cells( cell_map_.size() );
    size_t const d                  = mesh.type().facet_dim();
    size_t const num_facet_vertices = mesh.type().num_vertices( d );
    std::vector< size_t > facet_vertices( num_facet_vertices );
    CellType const &      celltype = mesh.type();

    for ( size_t i = 0; i < cell_map_.size(); ++i )
    {
      Facet facet( mesh, cell_map_[i] );
      for ( size_t v = 0; v < num_facet_vertices; ++v )
      {
        facet_vertices[v] = boundary_vertices[facet.entities( 0 )[v]];
      }
      // Reorder vertices so facet is right-oriented w.r.t. facet normal
      celltype.order_facet( facet_vertices.data(), facet );
      editor.add_cell( i, facet_vertices.data() );
    }

    // If the mesh is distributed, set global numbering and copy the ownership
    if ( mesh.is_distributed() )
    {
      this->distdata()[0].assign( mesh.distdata()[0], vertex_map_ );
      this->distdata()[d].assign( mesh.distdata()[d], cell_map_ );
    }
    editor.close();
  }

  message( 1, "BoundaryMesh : number of cells = %u, number of vertices %u",
              cell_map_.size(), vertex_map_.size() );
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */
