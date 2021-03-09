// Copyright (C) 2015 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_METIS_INTERFACE_H
#define __DOLFIN_METIS_INTERFACE_H

#include <dolfin/common/types.h>

#include <dolfin/mesh/MeshValues.h>

namespace dolfin
{

/**
 *  @namespace  MetisInterface
 *
 *  @brief   This namespace provides an interface to ParMETIS
 *
 */

namespace MetisInterface
{

void partitionCommonMetis( Mesh &                       mesh,
                           MeshValues< size_t, Cell > & partitions,
                           MeshValues< size_t, Cell > * weight );

void partitionGeomMetis( Mesh &                         mesh,
                         MeshValues< size_t, Vertex > & partitions );

} /* namespace MetisInterface */

} /* namespace dolfin */

#endif /* __DOLFIN_METIS_INTERFACE_H */
