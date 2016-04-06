// Copyright (C) 2008 Johan Jansson
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2009-2010.
// Modified by Balthasar Reuter, 2013.
//

#include <dolfin/common/constants.h>
#include <dolfin/config/dolfin_config.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/main/MPI.h>
#include <dolfin/math/basic.h>

#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/LoadBalancer.h>
#include <dolfin/mesh/RivaraRefinement.h>
#include <dolfin/mesh/MeshFunctionConverter.h>
#include <dolfin/parameter/parameters.h>
#include <dolfin/mesh/DMesh.h>
#include <dolfin/mesh/DVertex.h>
#include <dolfin/mesh/DCell.h>

#ifdef HAVE_LIBGEOM
#include <Geometry.h>
#endif

#ifdef HAVE_MPI
#include <mpi.h>
#endif

namespace dolfin
{

//-----------------------------------------------------------------------------
void RivaraRefinement::refine(Mesh& mesh,
                              MeshFunction<bool>& cell_marker,
                              real tf, real tb, real ts, bool balance)
{
  message("Refining simplicial mesh by recursive Rivara bisection without boundary smoothing.");

  // Start Loadbalancer
  if(MPI::numProcesses() > 1 && balance)
  {
    begin("Load balancing");
    // Tune loadbalancer using machine specific parameters, if available
    if( tf > 0.0 && tb > 0.0 && ts > 0.0)
    {
      LoadBalancer::balance(mesh, cell_marker, tf, tb, ts, LoadBalancer::LEPP);
    }
    else
    {
      LoadBalancer::balance(mesh, cell_marker, LoadBalancer::LEPP);
    }
    end();
  }

  DMesh dmesh;
  dmesh.imp(mesh);

  std::vector<bool> dmarked(mesh.num_cells());
  for (CellIterator ci(mesh); !ci.end(); ++ci)
  {
    if(cell_marker.get(ci->index()) == true)
    {
      dmarked[ci->index()] = true;
    }
    else
    {
      dmarked[ci->index()] = false;
    }
  }

  dmesh.bisectMarked(dmarked);

  Mesh omesh;
  dmesh.exp(omesh);

  mesh = omesh;
}
//-----------------------------------------------------------------------------
#ifdef HAVE_LIBGEOM
//-----------------------------------------------------------------------------
void RivaraRefinement::refine(Mesh& mesh,
                              MeshFunction<bool>& cell_marker,
                              libgeom::Geometry& geom,
                              MeshFunction<int>& patch_id_list,
                              MeshFunction<float>& bnd_u,
                              MeshFunction<float>& bnd_v,
                              real tf, real tb, real ts, bool balance)
{
  message("Refining simplicial mesh by recursive Rivara bisection with boundary smoothing");


  // Start Loadbalancer
  if(MPI::numProcesses() > 1 && balance)
  {
    MeshFunction<real> patch_id_double(mesh);
    MeshFunction<real> u_double(mesh);
    MeshFunction<real> v_double(mesh);

    MeshFunction<real> new_patch_id_double(mesh);
    MeshFunction<real> new_u_double(mesh);
    MeshFunction<real> new_v_double(mesh);

    MeshFunctionConverter::cast(patch_id_list, patch_id_double);
    MeshFunctionConverter::cast(bnd_u, u_double);
    MeshFunctionConverter::cast(bnd_v, v_double);

    Array< std::pair< MeshFunction<real> *, MeshFunction<real> * > > vertex_functions;

    vertex_functions.push_back( std::make_pair(&patch_id_double, &new_patch_id_double) );
    vertex_functions.push_back( std::make_pair(&u_double, &new_u_double) );
    vertex_functions.push_back( std::make_pair(&v_double, &new_v_double) );

    begin("Load balancing");
    // Tune loadbalancer using machine specific parameters, if available
    if( tf > 0.0 && tb > 0.0 && ts > 0.0)
    {
      LoadBalancer::balance(mesh, cell_marker, vertex_functions, tf, tb, ts, LoadBalancer::LEPP);
    }
    else
    {
      LoadBalancer::balance(mesh, cell_marker, vertex_functions, LoadBalancer::LEPP);
    }
    end();

    if(new_patch_id_double.size() > 0)
    {
      MeshFunctionConverter::cast(new_patch_id_double, patch_id_list);
    }
    if(new_u_double.size() > 0 )
    {
      MeshFunctionConverter::cast(new_u_double, bnd_u);
    }
    if(new_v_double.size() > 0 )
    {
      MeshFunctionConverter::cast(new_v_double, bnd_v);
    }
  }

  if (MPI::numProcesses() > 1) mesh.renumber();

  DMesh dmesh;
  dmesh.imp(mesh, patch_id_list, bnd_u , bnd_v);

  std::vector<bool> dmarked(mesh.num_cells());
  for (CellIterator ci(mesh); !ci.end(); ++ci)
  {
    if(cell_marker.get(*ci) == true)
    {
      dmarked[ci->index()] = true;
    }
    else
    {
      dmarked[ci->index()] = false;
    }
  }

  dmesh.bisectMarked(dmarked, geom);


  Mesh omesh;
  MeshFunction<int> patch_id_omesh(omesh);
  MeshFunction<float> u_omesh(omesh);
  MeshFunction<float> v_omesh(omesh);
  dmesh.exp(omesh, patch_id_omesh, u_omesh, v_omesh);

  mesh = omesh;
  bnd_u = MeshFunction<float>(mesh);
  bnd_v = MeshFunction<float>(mesh);
  patch_id_list = MeshFunction<int>(mesh);

  bnd_u.init(0);
  bnd_v.init(0);
  patch_id_list.init(0);

  for (VertexIterator v(mesh); !v.end(); ++v)
  {
    bnd_u.set(v->index(), u_omesh.get(v->index()));
    bnd_v.set(v->index(), v_omesh.get(v->index()));
    patch_id_list.set(v->index(), patch_id_omesh.get(v->index()));
  }

  MPI_Barrier(dolfin::MPI::DOLFIN_COMM);

}
//-----------------------------------------------------------------------------
#endif // HAVE_LIBGEOM
//-----------------------------------------------------------------------------

}

