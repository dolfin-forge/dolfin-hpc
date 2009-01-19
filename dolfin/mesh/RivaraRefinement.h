// Copyright (C) 2008 Johan Jansson
// Licensed under the GNU LGPL Version 2.1.
//

#ifndef __RIVARAREFINEMENT_H
#define __RIVARAREFINEMENT_H

#include <list>
#include <vector>

#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Edge.h>

namespace dolfin
{
  class RivaraRefinement
  {
  public:
    
    /// Refine simplicial mesh locally by recursive edge bisection 
    static void refine(Mesh& mesh, 
		       MeshFunction<bool>& cell_marker,
		       MeshFunction<uint>& cell_map);
  };

  class DVertex;
  class DCell;
  class DMesh;
  class DEdge;
  
  class DVertex
  {
  public:
    DVertex();

    int id;

    std::list<DCell *> cells;
    Point p;
  };
    
  class DCell
  {
  public:
    DCell();

    int id;
    int parent_id;

    std::vector<DVertex *> vertices;

    bool deleted;
  };
      
  class DMesh
  {
  public:
    DMesh();

    void addVertex(DVertex* v);
    void addCell(DCell* c,
		 std::vector<DVertex*> vs, int parent_id);
    void removeCell(DCell* c);

    void imp(Mesh& mesh);
    void exp(Mesh& mesh, std::vector<int>& new2old_arr);
    void number();

    void bisect(DCell* dcell, DVertex* hangv,
		DVertex* hv0, DVertex* hv1);

    void bisectMarked(std::vector<bool> marked_ids);

    DCell* opposite(DCell* dcell, DVertex* v1, DVertex* v2);


    std::list<DVertex *> vertices;
    std::list<DCell *> cells;

    CellType* cell_type;
    uint d;
  };


}

#endif
