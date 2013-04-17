// Copyright (C) 2008 Johan Jansson
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2009-2010.
//

#include <dolfin/main/MPI.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/LoadBalancer.h>
#include <dolfin/mesh/RivaraRefinement.h>
#include <dolfin/mesh/DMesh.h>
#include <dolfin/mesh/DVertex.h>
#include <dolfin/mesh/DCell.h>

#include <vector>
#include <list>

#ifdef HAVE_MPI
#include <mpi.h>
#endif 

using namespace dolfin;
//-----------------------------------------------------------------------------
void RivaraRefinement::refine(Mesh& mesh, 
			      MeshFunction<bool>& cell_marker,
			      real tf, real tb, real ts, bool balance)
{
  message("Refining simplicial mesh by recursive Rivara bisection.");

  // Start Loadbalancer
  if(MPI::numProcesses() > 1 && balance) {
    begin("Load balancing");
    // Tune loadbalancer using machine specific parameters, if available
    if( tf > 0.0 && tb > 0.0 && ts > 0.0)
      LoadBalancer::balance(mesh, cell_marker, tf, tb, ts, LoadBalancer::LEPP);
    else
      LoadBalancer::balance(mesh, cell_marker, LoadBalancer::LEPP);
    end();
  }

  if (MPI::numProcesses() > 1) mesh.renumber();

  DMesh dmesh;
  dmesh.imp(mesh);
  
  std::vector<bool> dmarked(mesh.numCells());
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
  
  dmesh.bisectMarked(dmarked);


  // Remove deleted cells from global list
  for(std::list<DCell* >::iterator it = dmesh.cells.begin();
      it != dmesh.cells.end(); )
  {
    
    DCell* dc = *it;
    
    if(dc->deleted)
      it = dmesh.cells.erase(it);
    else
      it++;
  }
  
  Mesh omesh;    
  dmesh.exp(omesh);
  
  mesh = omesh;
}
//-----------------------------------------------------------------------------
