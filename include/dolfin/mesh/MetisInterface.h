// Copyright (C) 2015 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2015-01-30
// Last changed: 2015-01-30

#ifndef __DOLFIN_METIS_INTERFACE_H
#define __DOLFIN_METIS_INTERFACE_H

#include <dolfin/common/types.h>

namespace dolfin
{

/**
 *  @class  MetisInterface
 *
 *  @brief   This class provides an interface to ParMETIS
 *
 */

class MetisInterface
{

public:
  
  static void partitionCommonMetis(Mesh& mesh,
                                   MeshValues<uint, Cell>& partitions,
                                   MeshValues<uint, Cell> * weight);

  static void partitionGeomMetis(Mesh& mesh,
                                 MeshValues<uint, Vertex> & partitions);

};

} /* namespace dolfin */

#endif /* __DOLFIN_METIS_INTERFACE_H */
