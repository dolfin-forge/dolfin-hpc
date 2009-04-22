// Copyright (C) 2008 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2008-03-03
// Last changed: 2008-07-17


#include "LoadBalancer.h"
#include "MeshFunction.h"
#include "Cell.h"
#include "Edge.h"
#include "Vertex.h"
#include <dolfin/common/Array.h>
#include <dolfin/main/MPI.h>

#include <mpi.h>

using namespace dolfin;
//-----------------------------------------------------------------------------
void LoadBalancer::balance(Mesh& mesh, MeshFunction<bool>& cell_marker)
{  
  balance(mesh, cell_marker, 0.0, 0.0, 0.0);
}
//-----------------------------------------------------------------------------
void LoadBalancer::balance(Mesh& mesh, MeshFunction<bool>& cell_marker,
			   real tf, real tb, real ts)
{  
  // Construct weight function 
  MeshFunction<uint> weight;
  uint w_local, w_sum, w_max;
  real w_avg;
  weight_function(mesh, cell_marker, weight, &w_local);
  
  // Preliminary evalution of load imbalance
  MPI_Allreduce(&w_local, &w_max, 1, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&w_local, &w_sum, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
  
  w_avg = (real) w_sum / (real) MPI::numProcesses();
  
  const real  threshold = 1.04;
  real imbalance = (real) w_max  / (real) w_avg  ;
  if( threshold > imbalance ) {
    message("Load imbalance %0.2f percent, below threshold.",
	    (imbalance - 1.0) * 100);
    return;
  }
  else
    message("Repartitioning %0.2f percent load imbalance.",
	    (imbalance - 1.0) * 100);

  // Repartition mesh
  MeshFunction<uint> partitions;
  mesh.partition(partitions, weight);

  // Calculate process reassignment
  uint max_sendrecv;
  process_reassignment(partitions, &max_sendrecv);
  
  // Determine if the computation gains from repartitioning/reassignment
  if( tf > 0.0 && tb > 0.0 && ts > 0.0 ) 
    if(!computational_gain(mesh, weight, partitions, max_sendrecv, tf, tb, ts))
      return;

  // Distribute mesh according to new partition function
  MeshFunction<bool> new_cell_marker;
  mesh.distribute(partitions, cell_marker, new_cell_marker);
  cell_marker.init(mesh, mesh.topology().dim());  
  for(CellIterator c(mesh); !c.end(); ++c)
    cell_marker.set(*c, new_cell_marker.get(*c));

  /*
   * REMOVE
   * Debug, calculate load imbalance after repartitioning
   */
  weight_function(mesh, cell_marker, weight, &w_local);
  
  MPI_Allreduce(&w_local, &w_max, 1, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&w_local, &w_sum, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
  
  w_avg = (real) w_sum / (real) MPI::numProcesses();
  real new_imbalance = (real)w_max  / (real)w_avg  ;
  message("%0.2f percent load imbalance after repartitioning.", 
	  (new_imbalance - 1.0) * 100);

}
//-----------------------------------------------------------------------------
void LoadBalancer::weight_function(Mesh& mesh, 
				   MeshFunction<bool>& cell_marker,
				   MeshFunction<uint>& weight,
				   uint* w_sum)
{
  // Calculated weights for the dual graph with a simulated mesh refinement

  weight.init(mesh, mesh.topology().dim()); 
  MeshFunction<bool> used_cell(mesh, mesh.topology().dim());
  MeshFunction<bool> used_edge(mesh, 1);
  used_cell = false;
  used_edge = false;
  real max, l;
  uint index = 0;
  weight = 1;
  //  *w_sum = mesh.numCells();
  for(CellIterator c(mesh); !c.end(); ++c) {
    if( cell_marker.get(*c) && !used_cell.get(*c)) {
      max = 0.0;
      for(EdgeIterator e(*c); !e.end(); ++e) {
	if(!used_edge.get(*e)){
	  l = e->length();
	  if(max < l) {
	    max = l;
	    index = e->index();	    
	  }
	}
      }
      if(max == 0.0)
	continue; 
      Edge le(mesh, index);
      for(CellIterator nc(le); !nc.end(); ++nc) {
	if(!used_cell.get(*nc)) {
	  //	  *w_sum++;
	  weight.set(*nc, weight.get(*nc) + 1);
	  used_cell.set(*nc, true);
	  for(EdgeIterator e(*nc); !e.end(); ++e)
	    used_edge.set(*e, true);
	}
      }
    }
  }

  *w_sum = 0;
  for(CellIterator c(mesh); !c.end(); ++c)
    *w_sum += weight.get(*c);
  

}
//-----------------------------------------------------------------------------
void LoadBalancer::process_reassignment(MeshFunction<uint>& partitions,
					uint* max_sendrecv)
{
  uint pe_size = MPI::numProcesses();
  uint rank = MPI::processNumber();

  // Construct my row in the similarity matrix
  uint *sim_row = new uint[pe_size];
  uint *sorted_indices = new uint[pe_size];
  memset(sim_row, 0, pe_size * sizeof(uint));
  for(uint i = 0; i<partitions.size(); i++)
    sim_row[ partitions.get(i) ] ++;

  // Gather similarity matrix on root node
  uint m = pe_size * pe_size ;
  uint *SiM = new uint[m]; // FIXME, this should only be needon on rank == 0
  //  MPI_Gather(sim_row, pe_size, MPI_UNSIGNED, 
  //	     SiM, pe_size, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

  MPI_Allgather(sim_row, pe_size, MPI_UNSIGNED,
		SiM, pe_size, MPI_UNSIGNED, MPI_COMM_WORLD);

  uint *sorted = new uint[m];
  pradixsort_matrix(&sorted_indices[0],&SiM[0], pe_size);
  MPI_Gather(sorted_indices, pe_size, MPI_UNSIGNED, 
	     sorted, pe_size, MPI_UNSIGNED, 0 , MPI_COMM_WORLD);


  uint *map  = new uint[pe_size];
  
  if(rank == 0) {
    
    uint num_assigned = 0;
    uint idx, idy;
    bool *unassigned_x  = new bool[pe_size];
    bool *unassigned_y  = new bool[pe_size];
    for(uint i = 0; i<pe_size; i++){
      unassigned_y[i] = true;
      unassigned_x[i] = true;
    }
      
    uint tail = m;
    // Greedy approximation for the MWBG problem
    while(num_assigned < pe_size) {
      for(uint i = tail; i > 0; i--) {
	idx = (sorted[i-1]) / pe_size;
	idy = (sorted[i-1]) % pe_size;
	if(unassigned_x[idx] && unassigned_y[idy]) {
	  unassigned_y[idy] = false;
	  unassigned_x[idx] = false;
	  map[idx] = idy;
	  tail = i;
	  num_assigned++;
	  break;
	}
      }
    }

    delete[] unassigned_x;
    delete[] unassigned_y;
  }

  MPI_Bcast(map, pe_size, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

  
  // Reassign processors
  for(uint i = 0; i < partitions.size(); i++)  
    partitions.set(i, map[ partitions.get(i) ]);

  // Calculate maximum number to send from processor
  *max_sendrecv = 0;
  for(uint i = 0; i < pe_size; i++)
    if( map[i] != rank )
      *max_sendrecv = std::max(*max_sendrecv, sim_row[i]);

  delete[] map;  
  delete[] sim_row;
  delete[] sorted_indices;
  delete[] sorted;
  delete[] SiM; // FIXME, this should only be needon on rank == 0

}
//-----------------------------------------------------------------------------
bool LoadBalancer::computational_gain(Mesh& mesh,
				      MeshFunction<uint>& weight,
				      MeshFunction<uint>& partitions,
				      uint max_sendrecv,
				      real tf, real tb, real ts)
{
  uint pe_size = MPI::numProcesses();

  uint w_old = 0;
  uint *tmp_w = new uint[pe_size];
  memset(tmp_w, 0, pe_size * sizeof(uint));
  for(CellIterator c(mesh); !c.end(); ++c) {
    w_old += weight.get(*c);      
    tmp_w[ partitions.get(*c) ] += weight.get(*c);      
  }

  uint w_new;
  for(uint i = 0; i < pe_size; i++)
    MPI_Reduce(&tmp_w[i], &w_new, 1, MPI_UNSIGNED, MPI_SUM, i, MPI_COMM_WORLD);

  uint w_oldmax, w_newmax;
  MPI_Allreduce(&w_old, &w_oldmax, 1, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&w_new, &w_newmax, 1, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);

  uint tmp = max_sendrecv; // No MPI aliasing on BG/L
  MPI_Allreduce(&tmp, &max_sendrecv, 1, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);
  message("**** %d - %d = %d  w_nmax / w_omax = %f maxsr = %d ****", w_oldmax, w_newmax, w_oldmax - w_newmax, (real) w_newmax/ (real) w_oldmax, max_sendrecv);

  real ndims =  (real) mesh.type().numVertices(mesh.topology().dim());
  real comp_cost = tf * ((ndims * ndims) * ( (real) (w_oldmax - w_newmax)));
  //  comp_cost +=  ts + tb *((real) (w_oldmax - w_newmax));
  real comm_cost = ts + tb * (ndims * (real) max_sendrecv);

  message("comp_cost %g comm_cost %g", comp_cost, comm_cost);

  delete[] tmp_w;
  return (comp_cost > comm_cost);
}
//-----------------------------------------------------------------------------
void LoadBalancer::radixsort_matrix(uint* res, uint* Matrix, uint m,bool desc)
{

  uint *index = new uint[m];
  uint *tmp = new uint[m];

  for(uint i = 0; i < m; i++)
    index[i] = i ;

  uint *count = new uint[256];
  uint *map = new uint[256];
  for(uint i = 0; i < 4; i++) {
    memcpy(tmp, index, m * sizeof(uint));
    memset(count, 0, 256 * sizeof(uint));
    for(uint j = 0; j < m; j++) 
      count[ ((Matrix[tmp[j]]) >> (8 * i)) & 0xff ]++;
        
    map[0] = 0;
    for(uint j = 1; j < 256; j++)
      map[j] = map[j-1] + count[j-1];
      
    for(uint j = 0; j < m; j++)
      index[ map[ ((Matrix[tmp[j]]) >> (8 * i)) & 0xff ]++ ] = tmp[j];        
  }

  if(desc) {
    uint j = m - 1;
    for(uint i = 0; i < m; i++)
      res[i] = index[j--];
  }
  else
    memcpy(res, index, m * sizeof(uint));

  delete[] map;
  delete[] count;
  delete[] tmp;
  delete[] index;
}
//-----------------------------------------------------------------------------
void LoadBalancer::pradixsort_matrix(uint* res, uint* Matrix, uint m)
{

  uint *index = new uint[m];
  uint *tmp = new uint[m];

  uint rank = MPI::processNumber();
  uint pe_size = MPI::numProcesses();

  uint *recvbuff = new uint[2 * pe_size];
  Array<uint> *sendbuff = new Array<uint>[pe_size];
  
  MPI_Status status;

  for(uint i = 0; i < m; i++) 
    index[i] = i + rank * pe_size;

  uint *count = new uint[256];
  uint *glb_count = new uint[256];
  uint *map = new uint[256];
  uint offset, glb_scan;
  for(uint i = 0; i < 4; i++) {
    memcpy(tmp, index, m * sizeof(uint));
    memset(count, 0, 256 * sizeof(uint));
    for(uint j = 0; j < m; j++) 
      count[ ((Matrix[tmp[j]]) >> (8 * i)) & 0xff ]++;
    

    MPI_Allreduce(count, glb_count, 256, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);

    offset = 0;
    for(uint j = 0; j < 256; j++) {
      MPI_Scan(&count[j], &glb_scan, 1, MPI_UNSIGNED, 
		 MPI_SUM, MPI_COMM_WORLD);
      map[j] = glb_scan + offset;
      offset += glb_count[j];      
    }

    for(uint j = 0; j < 256; j++) 
      map[j]--;

    for(uint j = m; j > 0 ; j--) {
      uint glb_index = map[ ((Matrix[tmp[j-1]]) >> (8 * i)) & 0xff ]--;
      uint target = floor( (double) glb_index / (double) pe_size);
      if(target != rank) {
	sendbuff[target].push_back(glb_index);
	sendbuff[target].push_back(tmp[j-1]);	
      }
      else
	index[glb_index%pe_size] = tmp[j-1];
    }
    uint src, dest;
    int recv_size;
    
    for(uint j = 1; j < pe_size; j++){

      src = (rank - j + pe_size) % pe_size;
      dest = (rank + j) % pe_size;

      MPI_Sendrecv(&sendbuff[dest][0], sendbuff[dest].size(), MPI_UNSIGNED,
		   dest, 0, recvbuff, 2*pe_size, MPI_UNSIGNED, src,
		   0, MPI_COMM_WORLD, &status);
      MPI_Get_count(&status,MPI_UNSIGNED,&recv_size);
      for(int k = 0; k < recv_size; k +=2){
	index[recvbuff[k]%pe_size]  = recvbuff[k+1];
      }
      sendbuff[dest].clear();
    }
  }

  memcpy(res, index, m * sizeof(uint));

  for(uint i = 0; i < pe_size; i++)
    sendbuff[i].clear();
  delete[] sendbuff;

  delete[] map;
  delete[] glb_count;
  delete[] count;
  delete[] recvbuff;
  delete[] tmp;
  delete[] index;


}
//-----------------------------------------------------------------------------


