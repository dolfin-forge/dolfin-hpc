// Copyright (C) 2015 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2015-01-30
// Last changed: 2015-01-30

#ifndef __ZOLTAN_INTERFACE_H
#define __ZOLTAN_INTERFACE_H

#include <dolfin/common/types.h>

#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_ZOLTAN
#include <zoltan_cpp.h>
#endif


namespace dolfin
{
  /// This class provides an interface to Zoltan
  
  class ZoltanInterface
  {
  public:

    static void partitionCommonZoltan(Mesh& mesh,
				      MeshFunction<uint>& partitions,
				      MeshFunction<uint>* weight);
    
    static void partitionGeomZoltan(Mesh& mesh, 
				    MeshFunction<uint>& partitions);
    
#ifdef HAVE_ZOLTAN
    
  private:
    
    /// Zoltan callbacks
    
    // Generall functions
    static int partitionZoltanNumCells(void *data, int *ierr);
    static int partitionZoltanNumVertices(void *data, int *ierr);
    
    static void partitionZoltanCellList(void *data, 
					int num_gid_entries, 
					int num_lid_entries, 
					ZOLTAN_ID_PTR global_ids, 
					ZOLTAN_ID_PTR local_ids, 
					int wgt_dim, 
					float *obj_wgts, 
					int *ierr);
    
    static void partitionZoltanVertexList(void *data, 
					  int num_gid_entries, 
					  int num_lid_entries, 
					  ZOLTAN_ID_PTR global_ids, 
					  ZOLTAN_ID_PTR local_ids, 
					  int wgt_dim, 
					  float *obj_wgts, 
					  int *ierr);
    
    // Geometry based functions
    static int partitionZoltanNumGeom(void *data, int *ierr);
    
    static void partitionZoltanGeomCoords(void *data, 
					  int num_gid_entries,
					  int num_lid_entries,
					  int num_obj, 
					  ZOLTAN_ID_PTR global_ids, 
					  ZOLTAN_ID_PTR local_ids, 
					  int num_dim, 
					  double *geom_vec, 
					  int *ierr);
    

    // Internal common partition function
    static void partitionZoltanInternal(Mesh& mesh, 
					MeshFunction<uint>& partitions, Zoltan *zz);
    
#endif
    
  };

}

#endif
