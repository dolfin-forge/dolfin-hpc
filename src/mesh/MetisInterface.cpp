// Copyright (C) 2015 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2015-01-30
// Last changed: 2015-01-30

#include <dolfin/config/dolfin_config.h>
#include <dolfin/common/timing.h>
#include <dolfin/math/basic.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/MeshRenumber.h>
#include <dolfin/mesh/MetisInterface.h>
#include <dolfin/parameter/parameters.h>

#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Cell.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#ifdef HAVE_PARMETIS
#include <parmetis.h>

#if PARMETIS_MAJOR_VERSION > 3
#define pm_idx_t  idx_t
#define pm_real_t real_t
#define pm_ncon   1
#else
#define pm_idx_t  idxtype
#define pm_real_t float
#define pm_ncon   0
#endif

#endif

#include <cstring>

namespace dolfin
{

#ifdef HAVE_PARMETIS
//-----------------------------------------------------------------------------
void MetisInterface::partitionCommonMetis(Mesh& mesh, 
                                          MeshFunction<uint>& partitions,
                                          MeshFunction<uint>* weight)
{
  if(!mesh.is_distributed())
  {
    return;
  }

  mesh.distdata()[0].renumber_global();
  tic();

  pm_real_t ubvec = 1.05;
  pm_idx_t numflag = 0;    // C-style numbering
  pm_idx_t edgecut = 0;
  pm_idx_t wgtflag, ncon;
  if (weight)
  {
    wgtflag = 2;    // Weights on vertices only
    ncon = 1;       // One weight per vertex
  }
  else
  {
    wgtflag = 0;    // Turn off graph weights
    ncon = pm_ncon; // No weights on vertices
  }

  // Duplicate MPI communicator
  MPI_Comm comm;
  MPI_Comm_dup(MPI::DOLFIN_COMM, &comm);

  // Get information about the PE
  int size, rank;
  MPI_Comm_size(MPI::DOLFIN_COMM, &size);
  MPI_Comm_rank(MPI::DOLFIN_COMM, &rank);

  pm_idx_t *elmdist = new pm_idx_t[size + 1];

  uint const tdim = mesh.topology().dim();
  int ncells = mesh.num_cells();

  if (ncells == 0)
  {
    dolfin::warning("MetisInterface : mesh partition contains zero cells.");
  }
  
  elmdist[rank] = ncells;
  MPI_Allgather(&ncells, 1, MPI_INT, elmdist, 1, MPI_INT, MPI::DOLFIN_COMM);

  pm_idx_t *elmwgt = NULL;

  if( weight )
  {
    elmwgt = new pm_idx_t[ncells];
    for(CellIterator c(mesh); !c.end(); ++c)
    {
      elmwgt[c->index()] = static_cast<pm_idx_t>(weight->get(*c));
    }
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

  int nvertices = mesh.type().num_vertices(tdim);
  int ncnodes = nvertices - 1;

  pm_idx_t *eptr = new pm_idx_t[ncells + 1];

  eptr[0] = 0;
  for(uint i=1;i < (mesh.num_cells() + 1);i++)
  {
    eptr[i] = eptr[i-1] + nvertices;
  }

  pm_idx_t *eind =  new pm_idx_t[nvertices * ncells];

  int i = 0;
  for (CellIterator c(mesh); !c.end(); ++c)
  {
    for (VertexIterator v(*c); !v.end(); ++v)
    {
      eind[i++] = v->global_index();
    }
  }

  pm_idx_t *part = new pm_idx_t[ncells];
  pm_real_t *tpwgts = new pm_real_t[size];

  for (i = 0; i < size; ++i)
  {
    tpwgts[i] = 1.0 / (pm_real_t) (size);
  }

  // default options
  pm_idx_t options[3] = { 1, 0, 15 };

  ParMETIS_V3_PartMeshKway(elmdist, eptr, eind, elmwgt, &wgtflag,&numflag,
                           &ncon,&ncnodes,&size, tpwgts, &ubvec, options,
                           &edgecut, part, &comm);

  delete[] eind;
  delete[] elmdist;
  delete[] tpwgts;
  delete[] eptr;
  delete[] elmwgt;

  // Create partition function
  partitions.init(mesh, tdim);
  partitions = size;
  uint lreassigned = 0;
  for (CellIterator cell(mesh); !cell.end(); ++cell)
  {
    partitions.set(*cell, (uint) part[cell->index()]);
    if(part[cell->index()] != rank)
    {
      ++lreassigned;
    }
  }
  uint greassigned = 0;
  MPI::numGlobalSum(lreassigned, greassigned);
  message(1, "METISPartMeshKway reassigned local = %2.2f %%, global = %2.2f %%",
          percent(lreassigned, mesh.size(tdim)),
          percent(greassigned, mesh.distdata()[tdim].global_size()));
  tocd();

  delete[] part;
  MPI_Comm_free(&comm);
}
//-----------------------------------------------------------------------------
void MetisInterface::partitionGeomMetis(Mesh& mesh, 
                                        MeshFunction<uint>& partitions)
{
  if(!mesh.is_distributed())
  {
    return;
  }

  tic();

  // Duplicate MPI communicator
  MPI_Comm comm;
  MPI_Comm_dup(MPI::DOLFIN_COMM, &comm);

  int size, rank;
  // Get information about the PE
  MPI_Comm_size(MPI::DOLFIN_COMM, &size);
  MPI_Comm_rank(MPI::DOLFIN_COMM, &rank);

  // Gather number of locally stored vertices for each processor
  pm_idx_t *vtxdist = new pm_idx_t[size+1];
  vtxdist[rank] = static_cast<pm_idx_t> (mesh.size(0));
  pm_idx_t local_vertices = vtxdist[rank];

  MPI_Allgather(&local_vertices, 1, MPI_INT, vtxdist, 1, 
		MPI_INT, MPI::DOLFIN_COMM);

  int i;
  pm_idx_t tmp;
  pm_idx_t sum = vtxdist[0];
  vtxdist[0] = 0;
  for (i = 1; i < size + 1; i++)
  {
    tmp = vtxdist[i];
    vtxdist[i] = sum;
    sum = tmp + sum;
  }

  pm_idx_t *part = new pm_idx_t[mesh.size(0)];
  pm_idx_t gdim =  static_cast<pm_idx_t>( mesh.geometry().dim() );
  pm_real_t *xdy = new pm_real_t[gdim * mesh.size(0)];

  i = 0;
  for (VertexIterator vertex(mesh); !vertex.end(); ++vertex, i += gdim)
  {
    for (int d = 0; d < gdim; ++d)
    {
      xdy[i + d] = static_cast<pm_real_t>(vertex->x()[d]);
    }
  }

  ParMETIS_V3_PartGeom(vtxdist, &gdim, xdy, part, &comm);

  // Create meshfunction from partitions
  partitions.init(mesh, 0);
  uint lreassigned = 0;
  for (VertexIterator vertex(mesh); !vertex.end(); ++vertex)
  {
    partitions.set(*vertex, static_cast<uint>(part[vertex->index()]));
    if(part[vertex->index()] != rank)
    {
      ++lreassigned;
    }
  }
  uint greassigned = 0;
  MPI::numGlobalSum(lreassigned, greassigned);
  message(1, "METISPartGeom reassigned local = %2.2f %%, global = %2.2f %%",
          percent(lreassigned, mesh.size(0)),
          percent(greassigned, mesh.distdata()[0].global_size()));
  tocd();

  delete[] xdy;
  delete[] part;
  delete[] vtxdist;
  MPI_Comm_free(&comm);
}
//-----------------------------------------------------------------------------
#else
//-----------------------------------------------------------------------------
void MetisInterface::partitionCommonMetis(Mesh& mesh, 
                                          MeshFunction<uint>& partitions,
                                          MeshFunction<uint>* weight)
{
  error("DOLFIN needs to be built with ParMetis support");
}
//-----------------------------------------------------------------------------
void MetisInterface::partitionGeomMetis(Mesh& mesh, 
                                        MeshFunction<uint>& partitions)
{
  error("DOLFIN needs to be built with ParMetis support");
}
//-----------------------------------------------------------------------------
#endif

} /* namespace dolfin */
