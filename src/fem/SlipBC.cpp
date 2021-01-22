// Copyright (C) 2007 Murtazo Nazarov
// Licensed under the GNU LGPL Version 2.1.

// Existing code for Dirichlet BC is used

#include <dolfin/fem/SlipBC.h>

#include <dolfin/fem/BilinearForm.h>
#include <dolfin/fem/NodeNormal.h>
#include <dolfin/fem/ScratchSpace.h>
#include <dolfin/fem/UFC.h>
#include <dolfin/la/PETScMatrix.h>
#include <dolfin/main/MPI.h>
#include <dolfin/mesh/SubDomain.h>
#include <dolfin/mesh/entities/Facet.h>
#include <dolfin/mesh/entities/Vertex.h>
#include <dolfin/mesh/entities/iterators/CellIterator.h>
#include <dolfin/mesh/entities/iterators/VertexIterator.h>
#include <dolfin/parameter/parameters.h>

#include <cmath>
#include <string>

#if ( __sgi )
#define fmax( a, b ) ( a > b ? a : b );
#endif

namespace dolfin
{

//-----------------------------------------------------------------------------
SlipBC::SlipBC( Mesh & mesh, SubDomain const & sub_domain )
  : BoundaryCondition( "SlipBC", mesh, sub_domain )
  , mesh( mesh )
  , node_normal( new NodeNormal( mesh ) )
  , node_normal_local( true )
  , As( nullptr )
  , As_local( true )
{
  // Do nothing
}
//-----------------------------------------------------------------------------
SlipBC::SlipBC( Mesh &            mesh,
                SubDomain const & sub_domain,
                NodeNormal &      normals )
  : BoundaryCondition( "SlipBC", mesh, sub_domain )
  , mesh( mesh )
  , node_normal( &normals )
  , node_normal_local( false )
  , As( nullptr )
  , As_local( true )
{
  // Do nothing
}
//-----------------------------------------------------------------------------
SlipBC::SlipBC( Mesh &            mesh,
                SubDomain const & sub_domain,
                SubSystem const & sub_system )
  : BoundaryCondition( "SlipBC", mesh, sub_domain, sub_system )
  , mesh( mesh )
  , node_normal( new NodeNormal( mesh ) )
  , node_normal_local( true )
  , As( nullptr )
  , As_local( true )
{
  // Do nothing
}
//-----------------------------------------------------------------------------
SlipBC::~SlipBC()
{
  if ( node_normal_local )
    delete node_normal;

  if ( As_local )
    delete As;
}
//-----------------------------------------------------------------------------
void SlipBC::apply( GenericMatrix &      A,
                    GenericVector &      b,
                    BilinearForm const & form )
{
  message( "Applying SlipBC boundary conditions to linear system." );

  FiniteElementSpace const & fullspace = form.test_space();

  // If the subsystem is not empty then the scratch space corresponds to the
  // given subspace, otherwise it is the space itself.
  ScratchSpace scratch( fullspace, sub_system() );

  // Optimize loops if the space is linear Lagrange: the space dimension is
  // the number of vertices times the number of components
  // bool const is_P1 = (scratch.space_dimension
  //                       == mesh.type().num_entities(0) * scratch.size);

  std::string const la_backend =
    dolfin_get< std::string >( "linear algebra backend" );

  if ( As == nullptr || As->size( 0 ) != A.size( 0 )
       || As->size( 1 ) != A.size( 1 ) )
  {
    // Create data structure for local assembly data
    if ( la_backend == "JANPACK" )
    {
      As       = reinterpret_cast< Matrix * >( &A );
      As_local = false;
    }
    else
    {
      delete As;
      As = new Matrix();
#if HAVE_PETSC
      ( *( As->instance() ) ).down_cast< PETScMatrix >().dup( A );
#endif
    }

    // Initialize ghosts for rhs vector using the full space dofmap
    {
      _ordered_set< size_t > rows;
      DofMap const &         fulldofmap = fullspace.dofmap();
      std::vector< size_t >  celldofs( fulldofmap.num_element_dofs() );

      for ( CellIterator c( mesh ); !c.end(); ++c )
      {
        scratch.cell.update( *c );
        fulldofmap.tabulate_dofs( celldofs.data(), scratch.cell );
        rows.insert( celldofs.begin(), celldofs.end() );
      }

      _ordered_map< size_t, size_t > mapping;
      b.init_ghosted( rows.size(), rows, mapping );
    }

    // Initialize normal field on given space and compute at the boundary
    if ( node_normal_local ) /// @todo add test for uninitialized external
                             /// NodeNormal objects
    {
      if ( sub_system().depth() == 0 )
      {
        node_normal->init( fullspace );
      }
      else
      {
        FiniteElementSpace subspace( fullspace, sub_system() );
        node_normal->init( subspace );
      }
      node_normal->compute();
    }

    // Initialize local data structures
    for ( size_t i = 0; i < mesh.type().dim(); ++i )
    {
      a[i].reserve( A.size( 0 ) );
      a_slip_row[i].reserve( A.size( 0 ) );
      a_col_indices[i].reserve( A.size( 0 ) );
      row_indices.insert( i );
    }
  }

  // Copy global stiffness matrix into temporary one (if needed)
  if ( la_backend != "JANPACK" )
  {
    *( As->instance() ) = A;
  }

  applySlipBC( A, b, form, scratch );

  // Apply changes in the temporary matrix
  As->apply();

  // Apply changes in the stiffness matrix and load vector (if needed)
  if ( la_backend != "JANPACK" )
  {
    A = *( As->instance() );
  }
  b.apply();
}
//-----------------------------------------------------------------------------
void SlipBC::apply( GenericMatrix &,
                    GenericVector &,
                    GenericVector const &,
                    BilinearForm const & )
{
  error( "SlipBC not implemented for non linear systems" );
}

//-----------------------------------------------------------------------------
void SlipBC::applySlipBC_P1( GenericMatrix &      A,
                             GenericVector &      b,
                             BilinearForm const & form,
                             ScratchSpace &       scratch )
{
  BoundaryMesh & boundary = mesh.exterior_boundary();
  if ( boundary.num_cells() )
  {
    // Used to get the global dof indices
    DofMap const &        Udofmap = form.test_space().dofmap();
    std::vector< size_t > fulldofs( Udofmap.num_element_dofs() );

    // If any, to get the node id
    bool const     same_space = ( sub_system().depth() == 0 );
    DofMap const & Ndofmap    = node_normal->basis()[0].space().dofmap();

    std::vector< size_t > node_Udofs;
    std::vector< size_t > node_Ndofs;
    size_t const          gdim = mesh.type().dim();
    size_t const          vdim =
      mesh.type().num_entities( 0 ); // number of nodes = vertices

    for ( VertexIterator v( boundary ); !v.end(); ++v )
    {
      Vertex vertex( mesh, boundary.vertex_index( *v ) );

      // Skip vertices not inside the sub domain
      if ( !this->sub_domain().enclosed( vertex, true ) )
      {
        continue;
      }

      size_t node = vertex.index();
      if ( vertex.is_owned() )
      {
        Cell cell( mesh, ( vertex.entities( gdim ) )[0] );
        scratch.cell.update( cell );

        // Find the vertex position in the cell
        std::vector< size_t > const & cvi = cell.entities( 0 );
        size_t                        ci  = 0;
        for ( ci = 0; ci < cell.num_entities( 0 ); ++ci )
        {
          if ( cvi[ci] == node )
            break;
        }
        size_t const vindex = ci;

        // Get the component dofs of U for the given node
        size_t voff = 0;
        Udofmap.tabulate_dofs( fulldofs.data(), scratch.cell );
        for ( size_t i = 0; i < scratch.size; ++i, voff += vdim )
        {
          node_Udofs.push_back( fulldofs[scratch.offset + voff + vindex] );
        }

        // The node id is the global index of the 1st comp. of the node normal
        // which is the same as the fullspace dof is there is no subsystem
        if ( same_space )
        {
          size_t node_id = node_Udofs[0];
          applyNodeBC( A, b, mesh, node_id, node_Udofs, node_Udofs );
          node_Udofs.clear();
        }
        else
        {
          // Get the component dofs of N for the given node
          size_t voff = 0;
          Ndofmap.tabulate_dofs( scratch.dofs, scratch.cell );
          for ( size_t i = 0; i < scratch.size; ++i, voff += vdim )
          {
            node_Ndofs.push_back( scratch.dofs[voff + vindex] );
          }
          size_t node_id = node_Ndofs[0];
          applyNodeBC( A, b, mesh, node_id, node_Udofs, node_Ndofs );
          node_Udofs.clear();
          node_Ndofs.clear();
        }
      }
    }
  }
}

//-----------------------------------------------------------------------------
void SlipBC::applySlipBC( GenericMatrix &      A,
                          GenericVector &      b,
                          BilinearForm const & form,
                          ScratchSpace &       scratch )
{
  // Be careful for now
  size_t const gdim = mesh.geometry_dimension();
  dolfin_assert( scratch.size == gdim );

  BoundaryMesh & boundary = mesh.exterior_boundary();
  if ( boundary.num_cells() != 0 )
  {
    DofMap const &        Udofmap = form.test_space().dofmap();
    std::vector< size_t > fulldofs( Udofmap.num_element_dofs() );

    // If any, to get the node id defined as the dof of the first normal comp.
    bool const     same_space = ( sub_system().depth() == 0 );
    DofMap const & Ndofmap    = node_normal->basis()[0].space().dofmap();

    //
    std::vector< size_t > node_Udofs;
    std::vector< size_t > node_Ndofs;
    size_t const          num_facet_dofs  = scratch.dof_map->num_facet_dofs();
    size_t const          num_facet_nodes = num_facet_dofs / scratch.size;
    size_t const          vdim =
      scratch.local_dimension / scratch.size; // component offset
    _set< size_t > visited_nodes;
    for ( CellIterator c( boundary ); !c.end(); ++c )
    {
      Facet facet( mesh, boundary.facet_index( *c ) );

      // Skip facets outside the sub domain
      if ( !this->sub_domain().overlap( facet, true ) )
      {
        continue;
      }

      //
      Cell         cell( mesh, facet.entities( gdim )[0] );
      size_t const local_facet = cell.index( facet );
      scratch.cell.update( cell );
      scratch.finite_element->tabulate_dof_coordinates(
        scratch.coordinates.data(), scratch.cell.coordinates.data() );
      // Attention, local-to-local mapping.
      scratch.dof_map->tabulate_facet_dofs( scratch.facet_dofs, local_facet );

      // Tabulate full space and nodenormal space if needed
      Udofmap.tabulate_dofs( fulldofs.data(), scratch.cell );
      if ( !same_space )
      {
        Ndofmap.tabulate_dofs( scratch.dofs, scratch.cell );
      }

      // Apply slipbc on each non ghosted facet dof node of the subspace
      for ( size_t ni = 0; ni < num_facet_nodes; ++ni )
      {
        // Get the cell local index of the first dof of node ni
        size_t const ni_celldof0 = scratch.facet_dofs[ni];
        // Get the global index of the first dof of node ni
        size_t const ni_Udof0 = fulldofs[scratch.offset + ni_celldof0];

        // Skip the node if it has already been encountered
        if ( visited_nodes.count( ni_Udof0 ) > 0 )
        {
          continue;
        }

        visited_nodes.insert( ni_Udof0 ); // mark as visited

        // Skip the node if its component dofs are ghosted
        if ( Udofmap.is_ghost( ni_Udof0 ) )
        {
          continue;
        }

        // Skip the node if the subdomain is defined geometrically
        if ( !sub_domain().inside(
               &scratch.coordinates[ni_celldof0 * Space::MAX_DIMENSION],
               true ) )
        {
          continue;
        }

        // Get the component dofs of U for the given node
        size_t voff = 0;
        for ( size_t i = 0; i < scratch.size; ++i, voff += vdim )
        {
          node_Udofs.push_back( fulldofs[scratch.offset + voff + ni_celldof0] );
          dolfin_assert( not Udofmap.is_ghost( node_Udofs.back() ) );
        }

        // The node id is the global index of the 1st comp. of the node normal
        // which is the same as the fullspace dof is there is no subsystem
        if ( same_space )
        {
          applyNodeBC( A, b, mesh, node_Udofs[0], node_Udofs, node_Udofs );
          node_Udofs.clear();
        }
        else
        {
          // Get the component dofs of N for the given node
          size_t voff = 0;
          for ( size_t i = 0; i < scratch.size; ++i, voff += vdim )
          {
            node_Ndofs.push_back( scratch.dofs[voff + ni_celldof0] );
          }

          applyNodeBC( A, b, mesh, node_Ndofs[0], node_Udofs, node_Ndofs );
          node_Udofs.clear();
          node_Ndofs.clear();
        }
      }
    }
    visited_nodes.clear();
  }
}

//-----------------------------------------------------------------------------
void SlipBC::applyNodeBC( GenericMatrix &               A,
                          GenericVector &               b,
                          Mesh const &                  mesh,
                          size_t const                  node_id,
                          std::vector< size_t > const & Udofs,
                          std::vector< size_t > const & Ndofs )
{
  // Naive reimplementation -- Aurélien
  // The node type defines the number of discriminated surfaces at the node.
  // Therefore it is the number of constrained directions up to the topological
  // dimension
  size_t const gdim   = mesh.topology_dimension();
  size_t const n_type = std::min( node_normal->node_type( node_id ), gdim );
  dolfin_assert( n_type > 0 );

  // Initialize set of row indices for reordering
  _ordered_set< size_t >           row_idx = row_indices;
  _ordered_set< size_t >::iterator it      = row_idx.begin();

  //--- Fill data structures for each space coordinate ---
  std::vector< Function > & basis_functions = node_normal->basis();
  for ( size_t i = 0; i < gdim; ++i )
  {
    // Copy non-zero entries from the stiffness matrix into local row
    A.getrow( Udofs[i], a_col_indices[i], a[i] );

    // Reset rhs for slip
    size_t nb_cols = a_col_indices[i].size();
    a_slip_row[i].resize( nb_cols );
    std::fill( a_slip_row[i].begin(), a_slip_row[i].end(), 0.0 );

    // Reset lhs for slip
    l_slip[i] = 0.0;

    // Zero the row in the matrix copy
    As->set( &a_slip_row[i][0], 1, &Udofs[i], nb_cols, &a_col_indices[i][0] );

    // Fill component i-th basis vector
    real( &v )[Space::MAX_DIMENSION] = basis_[i];
    basis_functions[i].vector().get( &v[0], gdim, &Ndofs[0] );

    // Determine maximum component (to be simplified)
    it     = row_idx.begin();
    max[i] = ( *it );
    for ( ; it != row_idx.end(); ++it )
    {
      if ( std::fabs( v[*it] ) > std::fabs( v[max[i]] ) )
      {
        max[i] = *it;
      }
    }
    row_idx.erase( max[i] );
  }

  // Integer n_type represents the number of constrained directions
  if ( n_type < gdim ) // At least one free direction
  {
    //--- For each constrained direction
    for ( size_t i = 0; i < n_type; ++i )
    {
      // Set row to the global index
      row[i] = Udofs[max[i]];

      // Update the LHS row with the vector components
      As->set( &basis_[i][0], 1, &row[i], gdim, &Udofs[0] );

      // Reset rhs for slip
      l_slip[i] = 0.0;
    }
    //--- Copy RHS to local vector (here since a surface node is most probable)
    b.get( &l[0], gdim, &Udofs[0] );
    //--- TODO: Row swapping of free directions assumes same row structure
    if ( gdim == 3 && n_type == 1 )
    {
      // Find position of the diagonal columns in the vectors
      size_t const maxn_1     = max[1];
      size_t       diag_idx_1 = 0;
      while ( a_col_indices[maxn_1][diag_idx_1] != Udofs[maxn_1] )
      {
        ++diag_idx_1;
      }
      real row_d1 = std::fabs( a[0][diag_idx_1] * basis_[1][0]
                               + a[1][diag_idx_1] * basis_[1][1]
                               + a[2][diag_idx_1] * basis_[1][2] )
                    - std::fabs( a[0][diag_idx_1] * basis_[2][0]
                                 + a[1][diag_idx_1] * basis_[2][1]
                                 + a[2][diag_idx_1] * basis_[2][2] );

      size_t const maxn_2     = max[2];
      size_t       diag_idx_2 = 0;
      while ( a_col_indices[maxn_2][diag_idx_2] != Udofs[maxn_2] )
      {
        ++diag_idx_2;
      }
      real row_d2 = std::fabs( a[0][diag_idx_2] * basis_[2][0]
                               + a[1][diag_idx_2] * basis_[2][1]
                               + a[2][diag_idx_2] * basis_[2][2] );
      row_d2 -= std::fabs( a[0][diag_idx_2] * basis_[1][0]
                           + a[1][diag_idx_2] * basis_[1][1]
                           + a[2][diag_idx_2] * basis_[1][2] );

      if ( row_d1 < 0.0 || row_d2 < 0.0 )
      {
        max[1] = maxn_2;
        max[2] = maxn_1;
      }
    }
    //--- For each free direction
    for ( size_t i = n_type; i < gdim; ++i )
    {
      // Set row to the global index
      row[i] = Udofs[max[i]];

      // Project equation on the tangential vector: a[j].tau_i
      // Beware, we assume here that the component contributions were
      // inserted in the same order !
      size_t const nb_cols = a_col_indices[max[i]].size();
      for ( size_t k = 0; k < gdim; ++k )
      {
        for ( size_t j = 0; j < nb_cols; ++j )
        {
          a_slip_row[i][j] += a[k][j] * basis_[i][k];
        }
        l_slip[i] += l[k] * basis_[i][k];
      }
      // Update the LHS row with the projection of the equation on tau_i
      As->setrow( row[i], a_col_indices[max[i]], a_slip_row[i] );
    }
    //--- Apply local equation RHS to the copy of the matrix
    b.set( &l_slip[0], gdim, &row[0] );
  }
  else // All directions are constrained
  {
    for ( size_t i = 0; i < gdim; ++i )
    {
      // Find position of the diagonal in the vectors
      size_t diag_idx = 0;
      while ( a_col_indices[i][diag_idx] != Udofs[i] )
      {
        ++diag_idx;
      }

      // Scale the diagonal entry and update the LHS diagonal
      real diag_val = std::fabs( a[i][diag_idx] );
      As->set( &diag_val, 1, &Udofs[i], 1, &Udofs[i] );
    }
    // Apply local equation RHS to the copy of the matrix
    b.set( &l_slip[0], gdim, &Udofs[0] );
  }

  //--- DEBUGGING
  // real entries[3] = { 1.0, 2.0, 3.0 };
  // b.set(&entries[0], gdim, &Udofs[0]);
  //
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
