// Copyright (C) 2007 Magnus Vikstrom.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Anders Logg, 2008.
// Modified by Niclas Jansson, 2008-2009.
// Modified by Aurélien Larcher 2012.
//
// First added:  2007-04-03
// Last changed: 2012-11-27

#include <dolfin/config/dolfin_config.h>
#include <dolfin/graph/Graph.h>
#include <dolfin/graph/GraphPartition.h>
#include <dolfin/mesh/MeshPartition.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/MeshRenumber.h>
#include <dolfin/parameter/parameters.h>

#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Cell.h>

#ifdef HAVE_MPI
#include <mpi.h>
#include <parmetis.h>
#endif

using namespace dolfin;

//-----------------------------------------------------------------------------
void MeshPartition::partition(Mesh& mesh, MeshFunction<uint>& partitions,
			      uint num_partitions)
{
  partitions.init(mesh, mesh.topology().dim());
  Graph graph(mesh);
  GraphPartition::partition(graph, num_partitions, partitions.values());
  
  bool report_edge_cut = dolfin_get("report edge cut");
  if (report_edge_cut)
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
#ifdef HAVE_MPI
//-----------------------------------------------------------------------------
void MeshPartition::partitionCommonMetis(Mesh& mesh,
					 MeshFunction<uint>& partitions,
					 MeshFunction<uint>* weight)
{

  // Metis assumes vertices numbered from process 0
  MeshRenumber::renumber_vertices(mesh);
  
#if PARMETIS_MAJOR_VERSION > 3
  real_t ubvec = 1.05;
  idx_t numflag = 0;    // C-style numbering
  idx_t edgecut = 0;
  idx_t wgtflag, ncon;
#else
  float ubvec = 1.05;
  int numflag = 0;   
  int edgecut = 0;
  int wgtflag, ncon;
#endif
  if (weight)
  {
    wgtflag = 2;    // Weights on vertices only
    ncon = 1;       // One weight per vertex
  }
  else
  {
    wgtflag = 0;    // Turn off graph weights
#if PARMETIS_MAJOR_VERSION > 3
    ncon = 1;       // No weights on vertices
#else
    ncon = 0;       // No weights on vertices
#endif
  }
  
  // Duplicate MPI communicator
  MPI_Comm comm;
  MPI_Comm_dup(MPI::DOLFIN_COMM, &comm);

  // Get information about the PE
  int size, rank;
  MPI_Comm_size(MPI::DOLFIN_COMM, &size);
  MPI_Comm_rank(MPI::DOLFIN_COMM, &rank);


#if PARMETIS_MAJOR_VERSION > 3
  idx_t *elmdist = new idx_t[size + 1];
#else
  idxtype *elmdist = new idxtype[size + 1];
#endif

  int ncells = mesh.numCells();

  if (ncells == 0)
  {
    dolfin::error("One mesh partition contains zero cells.");
  }
  
  elmdist[rank] = ncells;
  MPI_Allgather(&ncells, 1, MPI_INT, elmdist, 1, MPI_INT, MPI::DOLFIN_COMM);

#if PARMETIS_MAJOR_VERSION > 3
  idx_t *elmwgt = NULL;
#else
  idxtype *elmwgt = NULL;
#endif

  if( weight ) {
#if PARMETIS_MAJOR_VERSION > 3
    elmwgt = new idx_t[ncells];
#else
    elmwgt = new idxtype[ncells];
#endif
    for(CellIterator c(mesh); !c.end(); ++c) 
#if PARMETIS_MAJOR_VERSION > 3
      elmwgt[c->index()] = static_cast<idx_t>(weight->get(*c));
#else
      elmwgt[c->index()] = static_cast<idxtype>(weight->get(*c));
#endif
  }
  
  int sum_elm = elmdist[0];
  int tmp_elm;
  elmdist[0] = 0;
  for (int i = 1; i < size + 1; i++)
  {
    tmp_elm = elmdist[i];
    elmdist[i] = sum_elm;
    sum_elm = tmp_elm + sum_elm;
  }
  
  int nvertices = mesh.type().numVertices(mesh.topology().dim());
  int ncnodes = nvertices - 1;
  
#if PARMETIS_MAJOR_VERSION > 3
  idx_t *eptr = new idx_t[ncells + 1];
#else
  idxtype *eptr = new idxtype[ncells + 1];
#endif

  eptr[0] = 0;
  for(uint i=1;i < (mesh.numCells() + 1);i++)
    eptr[i] = eptr[i-1] + nvertices;

#if PARMETIS_MAJOR_VERSION > 3
  int *eind =  new idx_t[nvertices * ncells];  
#else
  int *eind =  new idxtype[nvertices * ncells];  
#endif

  int i = 0;
  for (CellIterator c(mesh); !c.end(); ++c)
    for (VertexIterator v(*c); !v.end(); ++v)
      eind[i++] = mesh.distdata().get_global(*v);


#if PARMETIS_MAJOR_VERSION > 3
  idx_t *part = new idx_t[ncells];
  real_t *tpwgts = new real_t[size];
#else
  idxtype *part = new idxtype[ncells];
  float *tpwgts = new float[size];
#endif


  for (i = 0; i < size; i++)
#if PARMETIS_MAJOR_VERSION > 3
    tpwgts[i] = 1.0 / (real_t) (size);
#else
    tpwgts[i] = 1.0 / (float) (size);
#endif

  // default options
#if PARMETIS_MAJOR_VERSION > 3
  idx_t options[3] = { 1, 0, 15 };
#else
  int options[3] = { 1, 0, 15 };
#endif
  
  ParMETIS_V3_PartMeshKway(elmdist, eptr, eind, elmwgt, &wgtflag, &numflag,
			   &ncon, &ncnodes, &size, tpwgts, &ubvec, options, &edgecut, part,
			   &comm);
  
  delete[] eind;
  delete[] elmdist;
  delete[] tpwgts;
  delete[] eptr;
  if (weight)
    delete[] elmwgt;

  // Create partition function
  partitions.init(mesh, mesh.topology().dim());
  partitions = size;
  for (CellIterator cell(mesh); !cell.end(); ++cell)
    partitions.set(*cell, (uint) part[cell->index()]);

  delete[] part;
  MPI_Comm_free(&comm);
}
//-----------------------------------------------------------------------------
void MeshPartition::partition_geom(Mesh& mesh, MeshFunction<uint>& partitions)
{
  // Duplicate MPI communicator
  MPI_Comm comm;
  MPI_Comm_dup(MPI::DOLFIN_COMM, &comm);
  
  int size, rank;
  // Get information about the PE
  MPI_Comm_size(MPI::DOLFIN_COMM, &size);
  MPI_Comm_rank(MPI::DOLFIN_COMM, &rank);

  // Gather number of locally stored vertices for each processor
#if PARMETIS_MAJOR_VERSION > 3
  idx_t *vtxdist = new idx_t[size+1];  
  vtxdist[rank] = static_cast<idx_t> (mesh.numVertices());
  idx_t local_vertices = vtxdist[rank];
#else
  idxtype *vtxdist = new idxtype[size+1];  
  vtxdist[rank] = static_cast<idxtype> (mesh.numVertices());
  idxtype local_vertices = vtxdist[rank];
#endif

  MPI_Allgather(&local_vertices, 1, MPI_INT, vtxdist, 1, MPI_INT,
		MPI::DOLFIN_COMM);

  int i;
#if PARMETIS_MAJOR_VERSION > 3
  idx_t tmp;
  idx_t sum = vtxdist[0];
#else
  int tmp;
  int sum = vtxdist[0];
#endif
  vtxdist[0] = 0;
  for (i = 1; i < size + 1; i++)
  {
    tmp = vtxdist[i];
    vtxdist[i] = sum;
    sum = tmp + sum;
  }

#if PARMETIS_MAJOR_VERSION > 3
  idx_t *part = new idx_t[mesh.numVertices()];
  idx_t gdim =  static_cast<idx_t>( mesh.geometry().dim() );
  real_t *xdy = new real_t[gdim * mesh.numVertices()];
#else
  idxtype *part = new idxtype[mesh.numVertices()];
  int gdim =  static_cast<int>( mesh.geometry().dim() );
  float *xdy = new float[gdim * mesh.numVertices()];
#endif


  i = 0;
  for (VertexIterator vertex(mesh); !vertex.end(); ++vertex)
  {
#if PARMETIS_MAJOR_VERSION > 3
    xdy[i] = static_cast<real_t>(vertex->point().x());
    xdy[i + 1] = static_cast<real_t>(vertex->point().y());
    if (gdim > 2)
      xdy[i + 2] = static_cast<real_t>(vertex->point().z());
#else
    xdy[i] = static_cast<float>(vertex->point().x());
    xdy[i + 1] = static_cast<float>(vertex->point().y());
    if (gdim > 2)
      xdy[i + 2] = static_cast<float>(vertex->point().z());
#endif
    i += gdim;
  }

  ParMETIS_V3_PartGeom(vtxdist, &gdim, xdy, part, &comm);

  // Create meshfunction from partitions
  partitions.init(mesh, 0);
  for (VertexIterator vertex(mesh); !vertex.end(); ++vertex)
    partitions.set(*vertex, static_cast<uint>(part[vertex->index()]));

  delete[] xdy;
  delete[] part;
  delete[] vtxdist;
  MPI_Comm_free(&comm);
}
//-----------------------------------------------------------------------------
#else
//-----------------------------------------------------------------------------
void MeshPartition::partitionCommonMetis(Mesh& mesh,
					 MeshFunction<uint>& partitions,
					 MeshFunction<uint>* weight)
{
  error("ParMetis needs MPI");
}
//-----------------------------------------------------------------------------
void MeshPartition::partition_geom(Mesh& mesh, MeshFunction<uint>& partitions)
{
  error("ParMetis needs MPI");
}
//-----------------------------------------------------------------------------
#endif

