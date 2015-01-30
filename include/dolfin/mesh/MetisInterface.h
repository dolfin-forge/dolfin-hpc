// Copyright (C) 2015 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2015-01-30
// Last changed: 2015-01-30

#ifndef __METIS_INTERFACE_H
#define __METIS_INTERFACE_H

#include <dolfin/common/types.h>

namespace dolfin
{
  /// This class provides an interface to ParMETIS 
  
  class MetisInterface
  {
  public:
    
    static void partitionCommonMetis(Mesh& mesh,
				     MeshFunction<uint>& partitions,
				     MeshFunction<uint>* weight);
   
    static void partitionGeomMetis(Mesh& mesh, 
				   MeshFunction<uint>& partitions);
  };

}

#endif
