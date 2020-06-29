// Copyright (C) 2010 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/adaptivity/AdaptiveRefinement.h>

#include <dolfin/common/system.h>
#include <dolfin/config/dolfin_config.h>
#include <dolfin/fem/DofMap.h>
#include <dolfin/fem/FiniteElementSpace.h>
#include <dolfin/fem/UFCCell.h>
#include <dolfin/function/Function.h>
#include <dolfin/function/FunctionDecomposition.h>
#include <dolfin/io/BinaryFile.h>
#include <dolfin/la/Vector.h>
#include <dolfin/mesh/LoadBalancer.h>
#include <dolfin/mesh/MPIMeshCommunicator.h>
#include <dolfin/mesh/MeshData.h>
#include <dolfin/mesh/RivaraRefinement.h>
#include <dolfin/parameter/parameters.h>

#include <algorithm>
#include <fstream>
#include <limits>

namespace dolfin
{

namespace AdaptiveRefinement
{

#ifdef DOLFIN_HAVE_MPI

void redistribute_func( Mesh & mesh, Function const & f,
                        Array< real > & values_, Array< uint > & rows_,
                        MeshValues< uint, Cell > const & distribution );

/// Project function on new mesh i.e. interpolation on non-matching meshes.
void project( Mesh & new_mesh, Array<Function> & f_post, Function & projected );

//-----------------------------------------------------------------------------
void refine(Mesh& mesh, MeshValues<bool, Cell>& cell_marker)
{
  message("Adaptive refinement");

  uint const numcellsbefore = mesh.num_global_cells();
  uint const numvertsbefore = mesh.global_size(0);
  message("  - cells    before: %d", numcellsbefore);
  message("  - vertices before: %d", numvertsbefore);

  std::string marked_filename("marked");
  std::string const marked_format = dolfin_get<std::string>("output_format");
  if (marked_format == "vtk")
  {
    marked_filename += ".pvd";
  }
  else if (marked_format == "binary")
  {
    marked_filename += ".bin";
  }

  File( marked_filename ) << cell_marker;

  if ( dolfin_get<std::string>("adapt_algorithm") == "rivara")
  {
    RivaraRefinement::refine(mesh, cell_marker);
  }
  else
  {
    dolfin::error("Unknown refinement algorithm");
  }

  uint const numcellsafter = mesh.num_global_cells();
  dolfin_assert(numcellsafter >= numcellsbefore);

  uint const numvertsafter = mesh.global_size(0);
  dolfin_assert(numvertsafter >= numvertsbefore);

  message("  - cells    after: %d", numcellsafter);
  message("  - vertices after: %d", numvertsafter);

  skip();
}
//-----------------------------------------------------------------------------
void refine_and_project( Mesh& mesh,
                         FunctionMapping const& functions,
                         MeshValues<bool, Cell>& cell_marker)
{
  dolfin_set("Load balancer redistribute", false);
  message("Adaptive refinement (with projection)");
  message("  - cells    before: %d", mesh.num_global_cells());
  message("  - vertices before: %d", mesh.global_size(0));

  std::string const refine_type = dolfin_get<std::string>("adapt_algorithm");
  if (refine_type == "simple")
  {
    LoadBalancer::balance(mesh, cell_marker);
  }
  else if (refine_type == "rivara")
  {
    LoadBalancer::balance(mesh, cell_marker, LoadBalancer::LEPP);
  }
  else
  {
    dolfin::error("Unknown refinement algorithm");
  }

  std::string marked_filename( "marked" );
  std::string const marked_format = dolfin_get<std::string>("output_format");
  if (marked_format == "vtk")
  {
    marked_filename += ".pvd";
  }
  else if (marked_format == "binary")
  {
    marked_filename += ".bin";
  }

  File( marked_filename ) << cell_marker;

  MeshValues<uint, Cell>& partitions = LoadBalancer::partitions(mesh);

  Array< Array<Function *> > coarse( functions.size() );
  Array< Array< Array< real > > > x_values( functions.size() );
  Array< Array< Array< uint > > > x_rows( functions.size() );

  for (uint f = 0; f < functions.size(); ++f )
  {
    Function * func = functions[f].second;
    uint const num_sub = func->space().element().num_sub_elements();

    if (num_sub == 0)
    {
      continue;
    }

    // make sure data is synchronized
    func->vector().apply();

    // decompose function
    coarse[f] = FunctionDecomposition::compute(*func);

    dolfin_assert( num_sub <= coarse[f].size() );

    // resize data redistribution arrays
    x_values[f].resize( num_sub );
    x_rows[f].resize( num_sub );

    // redistribute decomposed functions
    for (uint i = 0; i < num_sub; ++i)
    {
      AdaptiveRefinement::redistribute_func(mesh, *coarse[f][i], x_values[f][i],
                                            x_rows[f][i], partitions);
    }
  }

  // distribute the mesh and cell_marker according to the obtained partitioning
  MeshData D(mesh);
  D.add(cell_marker);
  MPIMeshCommunicator::distribute( partitions, &D );

  // global renumbering is not necessary
  mesh.distdata()[0].valid_numbering = true;

  // refine the mesh and keep a copy of the original / coarse mesh,
  // needed in the coarse functions
  Mesh new_mesh = mesh;
  RivaraRefinement::refine(new_mesh, cell_marker, 0.0, 0.0, 0.0, false);
  new_mesh.topology().renumber();

  for (uint f = 0; f < functions.size(); ++f )
  {
    FiniteElementSpace const& space = functions[f].second->space();
    uint const num_sub = space.element().num_sub_elements();

    Array<Function> post;
    /// @todo Invalid for scalar functions due to the zero subspace assumption
    for (uint i = 0; i < num_sub; ++i)
    {
      post.push_back( Function( FiniteElementSpace( coarse[f][i]->space() ) ) );

      dolfin_assert( x_values[f][i].size() == x_rows[f][i].size() );

      post.back().vector().set(x_values[f][i].data(),
                               x_values[f][i].size(),
                               x_rows[f][i].data());
      post.back().sync();
    }

    // project functions back on the mesh
    FiniteElementSpace projected_space(new_mesh, space);
    Function proj(projected_space);
    AdaptiveRefinement::project(new_mesh, post, proj);

#ifdef ENABLE_MPIIO
    std::string const filename = "projected_" + functions[f].first + ".bin";
#else
    std::stringstream pf;
    pf << "projected_" << functions[f].first << "_" << MPI::rank() << ".bin";
    std::string const filename = pf.str();
#endif

    File( filename ) << proj.vector();
  }

  // cleanup coarse functions
	for ( uint i = 0; i < functions.size(); ++i )
	{
		while ( !coarse[i].empty() )
		{
			delete coarse[i].back();
			coarse[i].pop_back();
		}
	}

	swap( mesh, new_mesh );
  mesh.topology().renumber();
  LoadBalancer::clear(mesh);

  message("  - cells    after: %d", mesh.num_global_cells());
  message("  - vertices after: %d", mesh.global_size(0));
}
//-----------------------------------------------------------------------------
void redistribute_func( Mesh& mesh, Function const& f,
                        Array<real> & values_,
                        Array<uint> & rows_,
                        MeshValues<uint, Cell> const& distribution )
{
  message( 1, "Redistributing Function: %p", &f );

  uint const pe_rank = MPI::rank();
  uint const pe_size = MPI::size();

  Array< Array< real > > send_buffer( pe_size );
  Array< Array< uint > > send_buffer_indices( pe_size );

  MeshValues<bool, Vertex> marked(mesh, false);

  Array<std::pair<uint, real> > recv_data;
  uint local_size = 0;

  for (CellIterator c(mesh); !c.end(); ++c)
  {
    uint target_proc = distribution(*c);

    /// @todo Only P1 friendly.
    for (VertexIterator v(*c); !v.end(); ++v)
    {
      real value = 0.0;
      uint global_index = v->global_index();
      f.vector().get(&value, 1, &global_index);

      if ( not v->is_ghost() and not marked(*v) )
      {
        if ( target_proc == pe_rank )
        {
          recv_data.push_back( std::make_pair(global_index, value) );
        }
        else
        {
          send_buffer[target_proc].push_back(value);
          send_buffer_indices[target_proc].push_back(global_index);
        }

        marked(*v) = true;
      }
    }

    local_size = std::max(local_size,
                          static_cast<uint>( send_buffer[target_proc].size() ));
  }

  uint recv_size = 0;
  MPI::all_reduce<MPI::max>( local_size, recv_size );

  Array< real > recv_buffer( recv_size, 0.0 );
  Array< uint > recv_buffer_indices( recv_size, 0 );

  for (uint j = 1; j < pe_size; ++j)
  {
    int src = (pe_rank - j + pe_size) % pe_size;
    int dst = (pe_rank + j) % pe_size;

#if defined( DEBUG )
    int recv_count_v =
#endif
    MPI::sendrecv( send_buffer[dst], dst, recv_buffer, src, 1 );

    int recv_count_i = MPI::sendrecv( send_buffer_indices[dst], dst,
                                      recv_buffer_indices, src, 2 );

    dolfin_assert( recv_count_v == recv_count_i );

    for (int i = 0; i < recv_count_i; ++i)
    {
      recv_data.push_back( std::make_pair( recv_buffer_indices[i],
                                           recv_buffer[i] ) );
    }
  }

  rows_.resize( recv_data.size() );
  values_.resize( recv_data.size() );

  for ( uint i = 0; i < recv_data.size(); ++i )
  {
    rows_[i] = recv_data[i].first;
    values_[i] = recv_data[i].second;
  }
}
//-----------------------------------------------------------------------------
void project( Mesh& new_mesh, Array<Function>& f_post, Function& projected )
{
  message( 1, "Projecting %u functions to function: %p",
           f_post.size(), &projected );

  FiniteElementSpace const & space = projected.space();

  Array< real > vv ( projected.vector().local_size() );
  Array< uint > indices( projected.vector().local_size() );
  Array< uint > local_indices( projected.space().dofmap().local_dimension() );
  uint i = 0;
  MeshValues<bool, Vertex> processed( new_mesh, false );

  projected.vector().zero();
  projected.sync();

  real gts_tol = dolfin_get<real>("GTS Tolerance");
  real geom_tol = dolfin_get<real>("Geometrical Tolerance Tetrahedron");

  dolfin_set("GTS Tolerance", 1e-10);
  dolfin_set("Geometrical Tolerance Tetrahedron", 1e-8);

  if ( ( space.family() != ufl::Family::CG ) || ( space.degree() != 1 ) )
  {
    error("AdaptiveRefinement::project only implemented for P1");
  }

  UFCCell ufccell( *CellIterator( new_mesh ) );
  for (CellIterator c(new_mesh); !c.end(); ++c)
  {
    ufccell.update(*c);
    space.dofmap().tabulate_dofs(local_indices.data(), ufccell, *c);

    /// @todo Only P1 friendly.
    for (VertexIterator v(*c); !v.end(); ++v)
    {
      uint ci   = 0;
      Array<uint> const & cvi = c->entities(0);
      for (; ci < c->num_entities(0); ci++)
      {
        if (cvi[ci] == v->index())
        {
          break;
        }
      }

      if (v->is_ghost() || processed(*v))
      {
        continue;
      }
      processed(*v) = true;

      real x[3] = { v->x()[0], v->x()[1], v->x()[2] };
      real test_value = 0.0;
      f_post[0].eval( &test_value, &x[0] );
      if (test_value == std::numeric_limits<real>::infinity())
      {
        for (EdgeIterator e(*v); !e.end(); ++e)
        {
          Array<uint> const & edge_v = e->entities(0);
          uint const index = ( edge_v[0] != v->index() ) ? 0 : 1;
          Vertex v_e( new_mesh, edge_v[index] );

          x[0] = v_e.x()[0];
          x[1] = v_e.x()[1];
          x[2] = v_e.x()[2];
          f_post[0].eval( &test_value, &x[0] );

          if (test_value != std::numeric_limits<real>::infinity())
          {
            break;
          }
        }

        if (test_value == std::numeric_limits<real>::infinity())
        {
          error("Couldn't find any suitable projection point");
        }
      }

      vv[i] = test_value;
      indices[i++] = local_indices[ci];

      for ( uint d = 1; d < f_post.size(); ++d )
      {
        f_post[d].eval( &test_value, &x[0] );
        if (test_value != std::numeric_limits<real>::infinity())
        {
          vv[i] = test_value;
          indices[i++] = local_indices[ci + d * c->num_entities(0)];
        }
      }
    }
  }

  dolfin_set("GTS Tolerance", gts_tol);
  dolfin_set("Geometrical Tolerance Tetrahedron", geom_tol);

  if (i > 0)
  {
    projected.vector().set(vv.data(), i, indices.data());
    projected.vector().apply();
  }

  projected.sync();
}
//-----------------------------------------------------------------------------
#else // no MPI
//-----------------------------------------------------------------------------
void refine( Mesh& /* unused */, MeshValues<bool, Cell>& /* unused */ )
{
  warning("Adaptive Refinement only implemented for MPI");
}
//-----------------------------------------------------------------------------
void refine_and_project( Mesh& /* unused */,
                         Array<Function *> const& /* unused */,
                         MeshValues<bool, Cell>& /* unused */)
{
  warning("Adaptive Refinement only implemented for MPI");
}
#endif

} // namespace AdaptiveRefinement

} // namespace dolfin

