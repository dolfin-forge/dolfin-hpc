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
#include <dolfin/parameter/parameters.h>
#include <dolfin/mesh/DMesh.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
void RivaraRefinement::refine(Mesh& mesh,
                              MeshValues<bool, Cell>& cell_marker,
                              real tf, real tb, real ts, bool balance)
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
  message("bisectMarked");
  dmesh.bisectMarked(cell_marker);

  Mesh omesh;
  dmesh.exp(omesh);

  mesh.swap(omesh);
  mesh.topology().renumber();
}
//-----------------------------------------------------------------------------

}

