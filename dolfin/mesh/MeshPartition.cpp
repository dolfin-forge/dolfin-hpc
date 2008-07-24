// Copyright (C) 2007 Magnus Vikstrom.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Anders Logg, 2008.
// Modified by Niclas Jansson, 2008.
//
// First added:  2007-04-03
// Last changed: 2008-02-11

#include <dolfin/graph/Graph.h>
#include <dolfin/graph/GraphPartition.h>
#include "MeshPartition.h"
#include "MeshFunction.h"
#include "MeshRenumber.h"
#include <dolfin/parameter/parameters.h>

#include "Vertex.h"
#include "Cell.h"

#include <parmetis.h>
#include <mpi.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
void MeshPartition::partition(Mesh& mesh,
                              MeshFunction<uint>& partitions,
                              uint num_partitions)
{
  partitions.init(mesh, mesh.topology().dim());
  Graph graph(mesh);
  GraphPartition::partition(graph, num_partitions, partitions.values());

  bool report_edge_cut = dolfin_get("report edge cut");
  if(report_edge_cut)
    GraphPartition::edgecut(graph, num_partitions, partitions.values());
}
//-----------------------------------------------------------------------------
void MeshPartition::partition(Mesh& mesh, MeshFunction<uint>& partitions)
{
  partitionCommonMetis(mesh, partitions, 0);
}
//-----------------------------------------------------------------------------
void MeshPartition::partition(Mesh& mesh, MeshFunction<uint>& partitions,
			      MeshFunction<uint>& weight)
{
  partitionCommonMetis(mesh, partitions, &weight);
}
//-----------------------------------------------------------------------------
void MeshPartition::partitionCommonMetis(Mesh& mesh, 
					 MeshFunction<uint>& partitions,
					 MeshFunction<uint>* weight)
{

  // Metis assumes vertices numbered from process 0 
  MeshRenumber::renumber_vertices(mesh);

  float ubvec = 1.05; /* magic magic */  
  int numflag = 0;    // C-style numbering
  int edgecut = 0;    
  int wgtflag, ncon;
  if( weight ) {
    wgtflag = 2;    // Weights on vertices only
    ncon = 1;       // One weight per vertex
  }
  else {
    wgtflag = 0;    // Turn off graph weights 
    ncon = 0;       // No weights on vertices
  }

  // Duplicate MPI communicator
  MPI_Comm comm;
  MPI_Comm_dup(MPI_COMM_WORLD, &comm);

  // Get information about the PE
  int size, rank;  
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);


  idxtype *elmdist = new idxtype[size + 1];
  int ncells = mesh.numCells();
  elmdist[rank] = ncells;
  MPI_Allgather(&elmdist[rank], 1, MPI_INT, elmdist, 
		1, MPI_INT, MPI_COMM_WORLD);

  idxtype *elmwgt = NULL;
  if( weight ) {
    elmwgt = new idxtype[ncells];
    for(CellIterator c(mesh); !c.end(); ++c) 
      elmwgt[c->index()] = static_cast<idxtype>(weight->get(*c));
  }

  int sum_elm = elmdist[0];  
  int tmp_elm;
  elmdist[0] = 0;
  for(int i=1;i<size+1;i++){    
    tmp_elm = elmdist[i];
    elmdist[i] = sum_elm;
    sum_elm = tmp_elm + sum_elm;
  }

  int nvertices = mesh.type().numVertices(mesh.topology().dim());
  int ncnodes = nvertices - 1 ;

  idxtype *eptr = new idxtype[ncells + 1];
  eptr[0] = 0;
  for(uint i=1;i < (mesh.numCells() + 1);i++)
    eptr[i] = eptr[i-1] + nvertices;

  int *eind =  new idxtype[nvertices * ncells];  
  int i = 0;
  for(CellIterator c(mesh); !c.end(); ++c)
    for(VertexIterator v(*c); !v.end(); ++v)
      eind[i++] = mesh.distdata().get_global(*v);

  idxtype *part = new idxtype[ncells];

  float *tpwgts = new float[size];
  for(i=0; i<size; i++)
    tpwgts[i] =  1.0/(float)(size); /* magic things happens */

  // default options
  int options[3] = {1, 0, 15};


  ParMETIS_V3_PartMeshKway(elmdist, eptr, eind, elmwgt, &wgtflag,&numflag,
                           &ncon,&ncnodes,&size, tpwgts, &ubvec,
                           options, &edgecut, part,&comm);

  delete[] eind;
  delete[] elmdist;
  delete[] tpwgts;
  delete[] eptr;
  if(weight)
    delete[] elmwgt;

  // Create partition function
  partitions.init(mesh, mesh.topology().dim());
  partitions = size;
  for(CellIterator cell(mesh); !cell.end(); ++cell)
    partitions.set(*cell, (uint) part[ cell->index() ]);

}
//-----------------------------------------------------------------------------
void MeshPartition::partition_geom(Mesh& mesh, MeshFunction<uint>& partitions)
{
  // Duplicate MPI communicator
  MPI_Comm comm; 
  MPI_Comm_dup(MPI_COMM_WORLD, &comm);

  int size, rank;
  // Get information about the PE
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  // Gather number of locally stored vertices for each processor
  idxtype *vtxdist = new idxtype[size+1];  
  vtxdist[rank] = static_cast<idxtype> (mesh.numVertices());


  MPI_Allgather(&vtxdist[rank], 1, MPI_INT, vtxdist, 1, 
		MPI_INT, MPI_COMM_WORLD);

  int i,tmp;
  int sum = vtxdist[0];  
  vtxdist[0] = 0;
  for(i=1;i<size+1;i++){    
    tmp = vtxdist[i];
    vtxdist[i] = sum;
    sum = tmp + sum;
  }

  idxtype *part = new idxtype[mesh.numVertices()];
  int gdim =  static_cast<int>( mesh.geometry().dim() );
  float *xdy = new float[gdim * mesh.numVertices()];
  // FIXME cast (double *) -> (float *)
  i = 0;
  for(VertexIterator vertex(mesh); !vertex.end(); ++vertex) {
    xdy[i] = static_cast<float>(vertex->point().x());
    xdy[i+1] = static_cast<float>(vertex->point().y());
    if(gdim > 2)
      xdy[i+2] = static_cast<float>(vertex->point().z());
    i +=gdim;
  }

  ParMETIS_V3_PartGeom(vtxdist,&gdim,xdy,part,&comm);

  // Create meshfunction from partitions
  partitions.init(mesh,0);
  for(VertexIterator vertex(mesh); !vertex.end(); ++vertex)
    partitions.set(*vertex, static_cast<uint>( part[vertex->index()]) );

  delete[] xdy;
  delete[] part;
  delete[] vtxdist;
}
//-----------------------------------------------------------------------------
