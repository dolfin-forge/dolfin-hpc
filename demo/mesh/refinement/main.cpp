// Copyright (C) 2006-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Johan Hoffman 2006.
// Modified by Niclas Jansson 2013. 
//
// First added:  2006-10-26
// Last changed: 2007-06-31

#include <dolfin.h>

using namespace dolfin;

int main()
{
  File *f_mesh;
  for (unsigned int j = 0; j < 2; j++) 
  {
    // Create mesh of unit square
    UnitSquare mesh(5, 5);
    
    if (j == 0)
      f_mesh = new File("simple.pvd");
    else
      f_mesh = new File("rivara.pvd");
    (*f_mesh) << mesh;
    
    // Uniform refinement
    mesh.refine();
    (*f_mesh) << mesh;
    
    // Refine mesh close to x = (0.5, 0.5)
    Point p(0.5, 0.5);
    for (unsigned int i = 0; i < 5; i++)
    {
      // Mark cells for refinement
      MeshFunction<bool> cell_markers(mesh, mesh.topology().dim());
      cell_markers = false;
      for (CellIterator c(mesh); !c.end(); ++c)
      {
	if (c->midpoint().distance(p) < 0.1)
	  cell_markers(*c) = true;
      }
      
      // Refine mesh    
      if (j == 0)
	mesh.refine(cell_markers); 
      else
	RivaraRefinement::refine(mesh, cell_markers);
      
      // Smooth mesh
      mesh.smooth();
      
      (*f_mesh) << mesh;
    }
    delete f_mesh;
  } 
  return 0;
}
