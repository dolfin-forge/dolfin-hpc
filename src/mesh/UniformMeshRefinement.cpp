// Copyright (C) 2006-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2008.
//
// First added:  2006-06-08
// Last changed: 2008-07-07

#include <dolfin/math/dolfin_math.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshTopology.h>
#include <dolfin/mesh/MeshGeometry.h>
#include <dolfin/mesh/MeshConnectivity.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/UniformMeshRefinement.h>

#include <dolfin/mesh/RefinementManager.h>
#include <dolfin/main/MPI.h>

#include <algorithm>


using namespace dolfin;

//-----------------------------------------------------------------------------
void UniformMeshRefinement::refine(Mesh& mesh)
{
  // Only know how to refine simplicial meshes
  refineSimplex(mesh);
}
//-----------------------------------------------------------------------------
void UniformMeshRefinement::refine(Mesh& mesh, libgeom::Geometry& geom, MeshFunction<int>& patch_id_list, MeshFunction<float>& bnd_u, MeshFunction<float>& bnd_v )
{
	if(mesh.geometry().dim()==3)
	{
  	// Only know how to refine simplicial meshes
  	refineSimplex(mesh, geom, patch_id_list, bnd_u,  bnd_v  );
	}
	else if(mesh.geometry().dim()==2)
	{
		refineSimplex(mesh, geom, patch_id_list, bnd_u);
	}
	else{
		std::cout<<"This is not implemented for the dimension of you geometry."<<std::endl;
	}
}
//-----------------------------------------------------------------------------
void UniformMeshRefinement::refine(Mesh& mesh, libgeom::Geometry& geom, MeshFunction<int>& patch_id_list, MeshFunction<float>& bnd_u )
{
  // Only know how to refine simplicial meshes
  refineSimplex(mesh, geom, patch_id_list, bnd_u  );
}
//-----------------------------------------------------------------------------
void UniformMeshRefinement::refineSimplex(Mesh& mesh)
{
  message(1, "Refining simplicial mesh uniformly.");
	//std::cout<<"Start of uniform refinement"<<std::endl;
  //mesh.disp();

  // Generate cell - edge connectivity if not generated
  mesh.init(mesh.topology().dim(), 1);
	//std::cout<<"After init(mesh.topology().dim(), 1)"<<std::endl;
  //mesh.disp();

  // Generate edge - vertex connectivity if not generated
  mesh.init(1, 0);
	//std::cout<<"After init(1, 0)"<<std::endl; 
	//mesh.disp();

  // Get cell type
  const CellType& cell_type = mesh.type();

  // Create a new refinement manager  FIXME:remove for the serial case
  RefinementManager refman(mesh);

  // Create new mesh and open for editing
  Mesh refined_mesh;
  MeshEditor editor;
  editor.open(refined_mesh, cell_type.cellType(),
	      mesh.topology().dim(), mesh.geometry().dim());
  
  // Get size of mesh
  const uint num_vertices = mesh.size(0);
  const uint num_edges = mesh.size(1);
  const uint num_cells = mesh.size(mesh.topology().dim());

  // Specify number of vertices and cells
  editor.initVertices(num_vertices + num_edges);
  editor.initCells(ipow(2, mesh.topology().dim())*num_cells);

  uint* edge_vert;
  Array<uint> shared_edge;
    
  uint vertex = 0;
  if(MPI::numProcesses() > 1){
    // Add old vertices
    for (VertexIterator v(mesh); !v.end(); ++v) {
      refined_mesh.distdata().set_map(vertex, mesh.distdata().get_global(*v), 0);
      
      if( mesh.distdata().is_ghost(v->index(), 0) ) {
	refined_mesh.distdata().set_ghost(vertex, 0);
	refined_mesh.distdata().set_ghost_owner(vertex,
						mesh.distdata().get_owner(*v), 0);
      }
      else if(mesh.distdata().is_shared(v->index(), 0))
	refined_mesh.distdata().set_shared(vertex, 0);

      editor.addVertex(vertex++, v->point());
    }
    
    for (EdgeIterator e(mesh); !e.end(); ++e) {
      edge_vert = e->entities(0);
      // If the edge is shared and lies between processes
      // process new vertex inside refinement manager
      if( refman.on_boundary(*e) ){
	
	// Add the new vertex inside the refinement manager
	refman.addVertex(edge_vert, vertex, refined_mesh);
	
	// Buffer edge information for mapping phase
	shared_edge.push_back(edge_vert[0]);
	shared_edge.push_back(edge_vert[1]);
	shared_edge.push_back(vertex);
      } 
      else
        refman.addVertex(vertex, refined_mesh);

      editor.addVertex(vertex++, e->midpoint());
    }

  }
  else {
    
    // Add old vertices
    for (VertexIterator v(mesh); !v.end(); ++v)
      editor.addVertex(vertex++, v->point());
    
    // Add new vertices
    for (EdgeIterator e(mesh); !e.end(); ++e)
      editor.addVertex(vertex++, e->midpoint());
    
  }

  // Add cells
  uint current_cell = 0;
  for (CellIterator c(mesh); !c.end(); ++c)
    cell_type.refineCell(*c, editor, current_cell);

  editor.close();
  
  // Map global numbers to unassigned shared vertices
  if(MPI::numProcesses() > 1)
    refman.map_new_vertices(shared_edge, mesh, refined_mesh);

  // Overwrite old mesh with refined mesh
  mesh = refined_mesh;
  mesh.distdata().invalid_numbering();
  mesh.renumber();

}

// overloaded method to include mesh boundary smoothing according to the geometry for a surface
//TODO the parallel case is untouched therefor no smoothing possible
void UniformMeshRefinement::refineSimplex(Mesh& mesh, libgeom::Geometry& geom, MeshFunction<int>& patch_id_list, MeshFunction<float>& bnd_u, MeshFunction<float>& bnd_v  )
{
  message(1, "Refining simplicial mesh uniformly with added boundary smoothing.");

  // Generate cell - edge connectivity if not generated
  mesh.init(mesh.topology().dim(), 1);
  

  // Generate edge - vertex connectivity if not generated
  mesh.init(1, 0);
  // Get cell type
  const CellType& cell_type = mesh.type();

  // Create a new refinement manager  FIXME:remove for the serial case
  RefinementManager refman(mesh);

  // Create new mesh and open for editing
  Mesh refined_mesh;
  MeshEditor editor;
  editor.open(refined_mesh, cell_type.cellType(),
	mesh.topology().dim(), mesh.geometry().dim());


  // Get size of mesh
  const uint num_vertices = mesh.size(0);
  const uint num_edges = mesh.size(1);
  const uint num_cells = mesh.size(mesh.topology().dim());

  // Specify number of vertices and cells
  editor.initVertices(num_vertices + num_edges);
  editor.initCells(ipow(2, mesh.topology().dim())*num_cells);

	//create the Meshfunctions for the refined mesh
	MeshFunction<int> refined_patch_id_list(refined_mesh);
	MeshFunction<float> refined_bnd_u(refined_mesh);
	MeshFunction<float> refined_bnd_v(refined_mesh);

	//initialize the new Meshfunctions
	refined_patch_id_list.init(0);
	refined_bnd_u.init(0);
	refined_bnd_v.init(0);
	
  uint* edge_vert;
  Array<uint> shared_edge;
    
  uint vertex = 0;
//parallel version without boundary smoothing so far TODO include boundary smoothing
  if(MPI::numProcesses() > 1){
    // Add old vertices
    for (VertexIterator v(mesh); !v.end(); ++v) {
      refined_mesh.distdata().set_map(vertex, mesh.distdata().get_global(*v), 0);
      
      if( mesh.distdata().is_ghost(v->index(), 0) ) {
				refined_mesh.distdata().set_ghost(vertex, 0);
				refined_mesh.distdata().set_ghost_owner(vertex,
						mesh.distdata().get_owner(*v), 0);
      }
      else if(mesh.distdata().is_shared(v->index(), 0))
				refined_mesh.distdata().set_shared(vertex, 0);

      editor.addVertex(vertex, v->point());

			refined_patch_id_list.set(	vertex, patch_id_list.get(v->index()));
			refined_bnd_u.set(					vertex,	bnd_u.get(v->index()));
			refined_bnd_v.set(					vertex, bnd_v.get(v->index()));
			vertex++;
    }
    
    for (EdgeIterator e(mesh); !e.end(); ++e) {
      edge_vert = e->entities(0);
			int patch_id_vertex0 = patch_id_list.get(edge_vert[0]);
			int patch_id_vertex1 = patch_id_list.get(edge_vert[1]);
			//don´t in general add the mitpoint. add midpoint if only one or no vertex of the edge is on the boundary. In the case that all the vertices ly on the boundary one needs to adjust the point so that it is also on the geometry.

			if(patch_id_vertex0 < 0 || patch_id_vertex1 < 0)
			{
      	// If the edge is shared and lies between processes
      	// process new vertex inside refinement manager
      	if( refman.on_boundary(*e) ){
	
					// Add the new vertex inside the refinement manager
					refman.addVertex(edge_vert, vertex, refined_mesh);
	
					// Buffer edge information for mapping phase
					shared_edge.push_back(edge_vert[0]);
					shared_edge.push_back(edge_vert[1]);
					shared_edge.push_back(vertex);
      	} 
     		else
      		refman.addVertex(vertex, refined_mesh);

     		editor.addVertex(vertex, e->midpoint());

				refined_patch_id_list.set( 	vertex, -1);
				refined_bnd_u.set(					vertex,	-1);
				refined_bnd_v.set(					vertex, -1);
				vertex++;
			}
			else{
				// If the edge is shared and lies between processes
      	// process new vertex inside refinement manager
      	if( refman.on_boundary(*e) )
				{
	
					// Add the new vertex inside the refinement manager
					refman.addVertex(edge_vert, vertex, refined_mesh);
	
					// Buffer edge information for mapping phase
					shared_edge.push_back(edge_vert[0]);
					shared_edge.push_back(edge_vert[1]);
					shared_edge.push_back(vertex);
      	} 
     	 	else
      	  refman.addVertex(vertex, refined_mesh);

				//TODO test with multiple patches
				if(patch_id_vertex0 == patch_id_vertex1)
				{
					float u,v;
					std::vector<float> pt;
					geom.get_midpoint(patch_id_vertex0, bnd_u.get(edge_vert[0]),  bnd_v.get(edge_vert[0]),  bnd_u.get(edge_vert[1]),  bnd_v.get(edge_vert[1]), u , v , pt);
					
					editor.addVertex(vertex, Point(pt[0],pt[1],pt[2]));
					refined_patch_id_list.set( 	vertex, patch_id_vertex0);
					refined_bnd_u.set(					vertex,	u);
					refined_bnd_v.set(					vertex, v);
					vertex++;
				}
				else
				{	
					//Vertex point1(mesh,edge_vert[0]);
					libgeom::Point3D midpoint(e->midpoint().x(), e->midpoint().y(), e->midpoint().z());
					libgeom::Point3D r0,r1;
					float u0, u1, v0, v1;
					real dist0 = geom.find_closest_point(midpoint, r0, u0, v0, patch_id_vertex0); 
					real dist1 = geom.find_closest_point(midpoint, r1, u1, v1, patch_id_vertex1);
					if(dist0 < dist1)
					{
						editor.addVertex(vertex, Point( r0.x(), r0.y(), r0.z() ) );
						refined_patch_id_list.set( 	vertex, patch_id_vertex0);
						refined_bnd_u.set(					vertex,	u0);
						refined_bnd_v.set(					vertex, v0);
						vertex++;

					}
					else
					{
						editor.addVertex(vertex, Point( r1.x(), r1.y(), r1.z() ) );
						refined_patch_id_list.set( 	vertex, patch_id_vertex1);
						refined_bnd_u.set(					vertex,	u1);
						refined_bnd_v.set(					vertex, v1);
						vertex++;
					}

				}
			}
    }

  }//serial version
  else {
    // Add old vertices
    for (VertexIterator v(mesh); !v.end(); ++v)
		{
      editor.addVertex(vertex, v->point());
			refined_patch_id_list.set(	vertex, patch_id_list.get(v->index()));
			refined_bnd_u.set(					vertex,	bnd_u.get(v->index()));
			refined_bnd_v.set(					vertex, bnd_v.get(v->index()));
			vertex++;
    }
    // Add new vertices
		// right now I assume single patch geometries... TODO fix to multiple patches
    for (EdgeIterator e(mesh); !e.end(); ++e)
		{    
			//don´t in general add the mitpoint. add midpoint if only one or no vertex of the edge is on the boundary. In the case that all the vertices ly on the boundary one needs to adjust the point so that it is also on the geometry.
			edge_vert = e->entities(0);
			
			int patch_id_vertex0 = patch_id_list.get(edge_vert[0]);
			int patch_id_vertex1 = patch_id_list.get(edge_vert[1]);
			if(patch_id_vertex0 < 0 || patch_id_vertex1 < 0)
			{
				editor.addVertex(vertex, e->midpoint());
				refined_patch_id_list.set( 	vertex, -1);
				refined_bnd_u.set(					vertex,	-1);
				refined_bnd_v.set(					vertex, -1);
				vertex++;
			}
			else{
				//TODO test with multiple patches
				if(patch_id_vertex0 == patch_id_vertex1)
				{
					float u,v;
					std::vector<float> pt;
					geom.get_midpoint(patch_id_vertex0, bnd_u.get(edge_vert[0]),  bnd_v.get(edge_vert[0]),  bnd_u.get(edge_vert[1]),  bnd_v.get(edge_vert[1]), u , v , pt);
					
					editor.addVertex(vertex, Point(pt[0],pt[1],pt[2]));
					refined_patch_id_list.set( 	vertex, patch_id_vertex0);
					refined_bnd_u.set(					vertex,	u);
					refined_bnd_v.set(					vertex, v);
					vertex++;
				}
				else
				{	
					//Vertex point1(mesh,edge_vert[0]);
					libgeom::Point3D midpoint(e->midpoint().x(), e->midpoint().y(), e->midpoint().z());
					libgeom::Point3D r0,r1;
					float u0, u1, v0, v1;
					real dist0 = geom.find_closest_point(midpoint, r0, u0, v0, patch_id_vertex0); 
					real dist1 = geom.find_closest_point(midpoint, r1, u1, v1, patch_id_vertex1);
					if(dist0 < dist1)
					{
						editor.addVertex(vertex, Point( r0.x(), r0.y(), r0.z() ) );
						refined_patch_id_list.set( 	vertex, patch_id_vertex0);
						refined_bnd_u.set(					vertex,	u0);
						refined_bnd_v.set(					vertex, v0);
						vertex++;

					}
					else
					{
						editor.addVertex(vertex, Point( r1.x(), r1.y(), r1.z() ) );
						refined_patch_id_list.set( 	vertex, patch_id_vertex1);
						refined_bnd_u.set(					vertex,	u1);
						refined_bnd_v.set(					vertex, v1);
						vertex++;
					}

				}
			}
    }
		
  }

  // Add cells
  uint current_cell = 0;
  for (CellIterator c(mesh); !c.end(); ++c)
    cell_type.refineCell(*c, editor, current_cell);

  editor.close();
  
  // Map global numbers to unassigned shared vertices
  if(MPI::numProcesses() > 1)
    refman.map_new_vertices(shared_edge, mesh, refined_mesh);

  // Overwrite old mesh with refined mesh
	
	std::swap(mesh, refined_mesh);
	bnd_u = MeshFunction<float>(mesh);
	bnd_u.init(0);
	bnd_v = MeshFunction<float>(mesh);
	bnd_v.init(0);
	patch_id_list =MeshFunction<int>(mesh);
	patch_id_list.init(0);
	
//TODO easiest but not fastest way. 
	for (VertexIterator v(mesh); !v.end(); ++v)
	{
		bnd_u.set(v->index(), refined_bnd_u.get(v->index()));
		bnd_v.set(v->index(), refined_bnd_v.get(v->index()));
		patch_id_list.set(v->index(), refined_patch_id_list.get(v->index()));
	}


  mesh.distdata().invalid_numbering();
  mesh.renumber();

}
//-----------------------------------------------------------------------------
// overloaded method to include mesh boundary smoothing according to the geometry for a surface
//TODO the parallel case is untouched therefor no smoothing possible
//TODO include libgeom geometry.h
void UniformMeshRefinement::refineSimplex(Mesh& mesh, libgeom::Geometry& geom, MeshFunction<int>& patch_id_list, MeshFunction<float>& bnd_u  )
{
  message(1, "Refining simplicial mesh uniformly with added boundary smoothing.");

  // Generate cell - edge connectivity if not generated
  mesh.init(mesh.topology().dim(), 1);
  

  // Generate edge - vertex connectivity if not generated
  mesh.init(1, 0);
  // Get cell type
  const CellType& cell_type = mesh.type();

  // Create a new refinement manager  FIXME:remove for the serial case
  RefinementManager refman(mesh);

  // Create new mesh and open for editing
  Mesh refined_mesh;
  MeshEditor editor;
  editor.open(refined_mesh, cell_type.cellType(),
	      mesh.topology().dim(), mesh.geometry().dim());


  // Get size of mesh
  const uint num_vertices = mesh.size(0);
  const uint num_edges = mesh.size(1);
  const uint num_cells = mesh.size(mesh.topology().dim());

  // Specify number of vertices and cells
  editor.initVertices(num_vertices + num_edges);
  editor.initCells(ipow(2, mesh.topology().dim())*num_cells);

	//create the Meshfunctions for the refined mesh
	MeshFunction<int> refined_patch_id_list(refined_mesh);
	MeshFunction<float> refined_bnd_u(refined_mesh);

	//initialize the new Meshfunctions
	refined_patch_id_list.init(0);
	refined_bnd_u.init(0);
	
  uint* edge_vert;
  Array<uint> shared_edge;
    
  uint vertex = 0;
//parallel version 
  if(MPI::numProcesses() > 1)
	{
    // Add old vertices
    for (VertexIterator v(mesh); !v.end(); ++v) 
		{
      refined_mesh.distdata().set_map(vertex, mesh.distdata().get_global(*v), 0);
      
      if( mesh.distdata().is_ghost(v->index(), 0) ) 
			{
				refined_mesh.distdata().set_ghost(vertex, 0);
				refined_mesh.distdata().set_ghost_owner(vertex, mesh.distdata().get_owner(*v), 0);
      }
      else if(mesh.distdata().is_shared(v->index(), 0))
			{
				refined_mesh.distdata().set_shared(vertex, 0);
			}
      editor.addVertex(vertex, v->point());

			refined_patch_id_list.set(	vertex, patch_id_list.get(v->index()));
			refined_bnd_u.set(					vertex,	bnd_u.get(v->index()));
			vertex++;
    }
    
    for (EdgeIterator e(mesh); !e.end(); ++e) 
		{
      edge_vert = e->entities(0);
			int patch_id_vertex0 = patch_id_list.get(edge_vert[0]);
			int patch_id_vertex1 = patch_id_list.get(edge_vert[1]);
			//don´t in general add the mitpoint. add midpoint if only one or no vertex of the edge is on the boundary. In the case that all the vertices ly on the boundary one needs to adjust the point so that it is also on the geometry.

			if(patch_id_vertex0 < 0 || patch_id_vertex1 < 0)
			{
      	// If the edge is shared and lies between processes
      	// process new vertex inside refinement manager
      	if( refman.on_boundary(*e) ){
	
					// Add the new vertex inside the refinement manager
					refman.addVertex(edge_vert, vertex, refined_mesh);
	
					// Buffer edge information for mapping phase
					shared_edge.push_back(edge_vert[0]);
					shared_edge.push_back(edge_vert[1]);
					shared_edge.push_back(vertex);
      	} 
     		else
      		refman.addVertex(vertex, refined_mesh);

     		editor.addVertex(vertex, e->midpoint());

				refined_patch_id_list.set( 	vertex, -1);
				refined_bnd_u.set(					vertex,	-1);
				vertex++;
			}
			else
			{
				// If the edge is shared and lies between processes
      	// process new vertex inside refinement manager
      	if( refman.on_boundary(*e) )
				{	
					// Add the new vertex inside the refinement manager
					refman.addVertex(edge_vert, vertex, refined_mesh);
	
					// Buffer edge information for mapping phase
					shared_edge.push_back(edge_vert[0]);
					shared_edge.push_back(edge_vert[1]);
					shared_edge.push_back(vertex);
      	} 
     	 	else
				{
      	  refman.addVertex(vertex, refined_mesh);
				}
     	 	//editor.addVertex(vertex++, e->midpoint());

				//this is the point where I assume there is only one patch I only use patch id of vertex 0 TODO fix it to multiple patches
				
				/*float u;
				std::vector<float> pt;
				geom.get_midpoint(patch_id_vertex0, bnd_u.get(edge_vert[0]),  bnd_u.get(edge_vert[1]),   u , pt);

				editor.addVertex(vertex, Point(pt[0],pt[1],pt[2]));
				refined_patch_id_list.set( 	vertex, patch_id_vertex0);
				refined_bnd_u.set(					vertex,	u);
				vertex++;*/
				if(patch_id_vertex0 == patch_id_vertex1)
				{
					float u;
					std::vector<float> pt;
					geom.get_midpoint(patch_id_vertex0, bnd_u.get(edge_vert[0]),  bnd_u.get(edge_vert[1]), u , pt);
					
					editor.addVertex(vertex, Point(pt[0],pt[1],pt[2]));
					refined_patch_id_list.set( 	vertex, patch_id_vertex0);
					refined_bnd_u.set(					vertex,	u);
					vertex++;
				}
				else
				{	
					//Vertex point1(mesh,edge_vert[0]);
					libgeom::Point3D midpoint(e->midpoint().x(), e->midpoint().y(), e->midpoint().z());
					libgeom::Point3D r0,r1;
					libgeom::REAL u0, u1;
					real dist0 = geom.find_closest_point(midpoint, r0, u0, patch_id_vertex0); 
					real dist1 = geom.find_closest_point(midpoint, r1, u1, patch_id_vertex1);
					if(dist0 < dist1)
					{
						editor.addVertex(vertex, Point( r0.x(), r0.y(), r0.z() ) );
						refined_patch_id_list.set( 	vertex, patch_id_vertex0);
						refined_bnd_u.set(					vertex,	u0);
						vertex++;

					}
					else
					{
						editor.addVertex(vertex, Point( r1.x(), r1.y(), r1.z() ) );
						refined_patch_id_list.set( 	vertex, patch_id_vertex1);
						refined_bnd_u.set(					vertex,	u1);
						vertex++;
					}

				}
			}
    }

  }//serial version
  else {
    // Add old vertices
    for (VertexIterator v(mesh); !v.end(); ++v)
		{
      editor.addVertex(vertex, v->point());
			refined_patch_id_list.set(	vertex, patch_id_list.get(v->index()));
			refined_bnd_u.set(					vertex,	bnd_u.get(v->index()));
			vertex++;
    }
    // Add new vertices
		// right now I assume single patch geometries... TODO fix to multiple patches
    for (EdgeIterator e(mesh); !e.end(); ++e)
		{    
			//don´t in general add the mitpoint. add midpoint if only one or no vertex of the edge is on the boundary. In the case that all the vertices ly on the boundary one needs to adjust the point so that it is also on the geometry.
			edge_vert = e->entities(0);
			
			int patch_id_vertex0 = patch_id_list.get(edge_vert[0]);
			int patch_id_vertex1 = patch_id_list.get(edge_vert[1]);
			if(patch_id_vertex0 < 0 || patch_id_vertex1 < 0)
			{
				editor.addVertex(vertex, e->midpoint());
				refined_patch_id_list.set( 	vertex, -1);
				refined_bnd_u.set(					vertex,	-1);
				vertex++;
			}
			else{
				//this is the point where I assume there is only one patch I only use patch id of vertex 0 TODO fix it to multiple patches
				/*float u;
				std::vector<float> pt;
				geom.get_midpoint(patch_id_vertex0, bnd_u.get(vertices[0]),  bnd_u.get(vertices[1]), u  , pt);*/

				if(patch_id_vertex0 == patch_id_vertex1)
				{
					float u;
					std::vector<float> pt;
					geom.get_midpoint(patch_id_vertex0, bnd_u.get(edge_vert[0]),  bnd_u.get(edge_vert[1]), u , pt);
					
					editor.addVertex(vertex, Point(pt[0],pt[1],pt[2]));
					refined_patch_id_list.set( 	vertex, patch_id_vertex0);
					refined_bnd_u.set(					vertex,	u);
					vertex++;
				}
				else
				{	
					//Vertex point1(mesh,edge_vert[0]);
					libgeom::Point3D midpoint(e->midpoint().x(), e->midpoint().y(), e->midpoint().z());
					libgeom::Point3D r0,r1;
					float u0, u1;
					real dist0 = geom.find_closest_point(midpoint, r0, u0, patch_id_vertex0); 
					real dist1 = geom.find_closest_point(midpoint, r1, u1, patch_id_vertex1);
					if(dist0 < dist1)
					{
						editor.addVertex(vertex, Point( r0.x(), r0.y(), r0.z() ) );
						refined_patch_id_list.set( 	vertex, patch_id_vertex0);
						refined_bnd_u.set(					vertex,	u0);
						vertex++;

					}
					else
					{
						editor.addVertex(vertex, Point( r1.x(), r1.y(), r1.z() ) );
						refined_patch_id_list.set( 	vertex, patch_id_vertex1);
						refined_bnd_u.set(					vertex,	u1);
						vertex++;
					}

				}
				//std::cout<<"u: "<<bnd_u.get(vertices[0])<<" "<<bnd_u.get(vertices[1])<<std::endl;
				//std::cout<<"midpoint"<<e->midpoint().x()<<" "<<e->midpoint().y()<<" "<<e->midpoint().z()<<std::endl;
				//std::cout<<"new midpoint"<<pt[0]<<" "<<pt[1]<<" "<<pt[2]<<std::endl;
				//editor.addVertex(vertex, Point(pt[0],pt[1],pt[2]));
				//refined_patch_id_list.set( 	vertex, patch_id_vertex0);
				//refined_bnd_u.set(					vertex,	u);
				//vertex++;
			}
    }
		
  }

  // Add cells
  uint current_cell = 0;
  for (CellIterator c(mesh); !c.end(); ++c)
    cell_type.refineCell(*c, editor, current_cell);

  editor.close();
  
  // Map global numbers to unassigned shared vertices
  if(MPI::numProcesses() > 1)
    refman.map_new_vertices(shared_edge, mesh, refined_mesh);

  // Overwrite old mesh with refined mesh
	
	std::swap(mesh, refined_mesh);
	bnd_u = MeshFunction<float>(mesh);
	bnd_u.init(0);
	patch_id_list =MeshFunction<int>(mesh);
	patch_id_list.init(0);
	
//TODO easiest but not fastest way. 
	for (VertexIterator v(mesh); !v.end(); ++v)
	{
		bnd_u.set(v->index(), refined_bnd_u.get(v->index()));
		patch_id_list.set(v->index(), refined_patch_id_list.get(v->index()));
	}

	//mesh.disp();
  mesh.distdata().invalid_numbering();
  mesh.renumber();

}
//-----------------------------------------------------------------------------
