// Copyright (C) 2007 Magnus Vikstrom.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/mesh/MeshPartition.h>

#include <dolfin/config/dolfin_config.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshRenumber.h>
#include <dolfin/mesh/MetisInterface.h>
#include <dolfin/mesh/ZoltanInterface.h>
#include <dolfin/parameter/parameters.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
#ifdef DOLFIN_HAVE_MPI
void MeshPartition::partition(MeshValues<uint, Cell>& partitions)
{
  const std::string method = dolfin_get<std::string>("Mesh partitioner");

  if (method == "parmetis")
    MetisInterface::partitionCommonMetis(partitions.mesh(), partitions, 0);
  else if (method == "zoltan")
    ZoltanInterface::partitionCommonZoltan(partitions.mesh(), partitions, 0);
  else
    error("Unknown mesh partitioner");
}
//-----------------------------------------------------------------------------
void MeshPartition::partition(MeshValues<uint, Cell>& partitions,
                              MeshValues<uint, Cell>& weight)
{
  const std::string method = dolfin_get<std::string>("Mesh partitioner");

  if (method == "parmetis")
    MetisInterface::partitionCommonMetis(partitions.mesh(), partitions, &weight);
  else if (method == "zoltan")
    ZoltanInterface::partitionCommonZoltan(partitions.mesh(), partitions, &weight);
  else
    error("Unknown mesh partitioner");
}
//-----------------------------------------------------------------------------
void MeshPartition::partition_geom(MeshValues<uint, Vertex>& partitions)
{
  const std::string method = dolfin_get<std::string>("Mesh partitioner");

  if (method == "parmetis")
    MetisInterface::partitionGeomMetis(partitions.mesh(), partitions);
  else if (method == "zoltan")
    ZoltanInterface::partitionGeomZoltan(partitions.mesh(), partitions);
  else
    error("Unknown mesh partitioner");
}
//-----------------------------------------------------------------------------
#else
void MeshPartition::partition(MeshValues<uint, Cell>&)
{
  error("Mesh partitioning requires MPI");
}
//-----------------------------------------------------------------------------
void MeshPartition::partition(MeshValues<uint, Cell>&, MeshValues<uint, Cell>&)
{
  error("Mesh partitioning requires MPI");
}
//-----------------------------------------------------------------------------
void MeshPartition::partition_geom(MeshValues<uint, Vertex>&)
{
  error("Geometric mesh partitioning requires MPI");
}
#endif
//-----------------------------------------------------------------------------

} /* namespace dolfin */
