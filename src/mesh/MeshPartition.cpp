// Copyright (C) 2007 Magnus Vikstrom.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Anders Logg, 2008.
// Modified by Niclas Jansson, 2008-2015.
// Modified by Aurélien Larcher 2012.
//
// First added:  2007-04-03
// Last changed: 2015-01-05

#include <dolfin/config/dolfin_config.h>
#include <dolfin/graph/Graph.h>
#include <dolfin/graph/GraphPartition.h>
#include <dolfin/mesh/MeshPartition.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/MeshRenumber.h>
#include <dolfin/mesh/MetisInterface.h>
#include <dolfin/mesh/ZoltanInterface.h>
#include <dolfin/parameter/parameters.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
void MeshPartition::partition(Mesh& mesh, MeshFunction<uint>& partitions,
			      uint num_partitions)
{
  partitions.init(mesh, mesh.topology().dim());
  Graph graph(mesh);
  GraphPartition::partition(graph, num_partitions, &partitions[0]);
  
  bool report_edge_cut = dolfin_get("report edge cut");
  if (report_edge_cut)
    GraphPartition::edgecut(graph, num_partitions, &partitions[0]);
}
//-----------------------------------------------------------------------------
#ifdef HAVE_MPI
void MeshPartition::partition(Mesh& mesh, MeshFunction<uint>& partitions)
{
  const std::string method = dolfin_get("Mesh partitioner");

  if (method == "parmetis")
    MetisInterface::partitionCommonMetis(mesh, partitions, 0);
  else if (method == "zoltan")
    ZoltanInterface::partitionCommonZoltan(mesh, partitions, 0);
  else
    error("Unknown mesh partitioner");
}
//-----------------------------------------------------------------------------
void MeshPartition::partition(Mesh& mesh, MeshFunction<uint>& partitions,
			      MeshFunction<uint>& weight)
{
  const std::string method = dolfin_get("Mesh partitioner");

  if (method == "parmetis")
    MetisInterface::partitionCommonMetis(mesh, partitions, &weight);
  else if (method == "zoltan")
    ZoltanInterface::partitionCommonZoltan(mesh, partitions, &weight);
  else
    error("Unknown mesh partitioner");
}
//-----------------------------------------------------------------------------
void MeshPartition::partition_geom(Mesh& mesh, MeshFunction<uint>& partitions)
{
  const std::string method = dolfin_get("Mesh partitioner");

  if (method == "parmetis")
    MetisInterface::partitionGeomMetis(mesh, partitions);
  else if (method == "zoltan")
    ZoltanInterface::partitionGeomZoltan(mesh, partitions);
  else
    error("Unknown mesh partitioner");
}
//-----------------------------------------------------------------------------
#else
void MeshPartition::partition(Mesh& mesh, MeshFunction<uint>& partitions)
{
  error("Mesh partitioning requires MPI");
}
//-----------------------------------------------------------------------------
void MeshPartition::partition(Mesh& mesh, MeshFunction<uint>& partitions,
			      MeshFunction<uint>& weight)
{
  error("Mesh partitioning requires MPI");
}
//-----------------------------------------------------------------------------
void MeshPartition::partition_geom(Mesh& mesh, MeshFunction<uint>& partitions)
{
  error("Geometric mesh partitioning requires MPI");
}
#endif
//-----------------------------------------------------------------------------


