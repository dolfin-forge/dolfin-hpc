// Copyright (C) 2008 Johan Jansson
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/mesh/utilities/RivaraRefinement.h>

#include <dolfin/mesh/DMesh.h>
#include <dolfin/mesh/LoadBalancer.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
void RivaraRefinement::refine( Mesh& mesh, MeshValues<bool, Cell>& cell_marker,
                               real tf, real tb, real ts, bool balance )
{
  message("Refining simplicial mesh by recursive Rivara bisection without boundary smoothing.");

  // Start Loadbalancer
  if(MPI::size() > 1 && balance)
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

  DMesh dmesh(mesh);
  dmesh.bisectMarked(cell_marker);

  Mesh omesh(mesh.type(), mesh.space(), mesh.topology().comm());
  dmesh.exp(omesh);

  swap( mesh, omesh );
  mesh.topology().renumber();
}
//-----------------------------------------------------------------------------

}

