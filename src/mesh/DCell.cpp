// Copyright (C) 2008 Johan Jansson
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2009-2010.
// Modified by Balthasar Reuter, 2013
//

#include <dolfin/mesh/DCell.h>
#include <dolfin/mesh/DVertex.h>

using namespace dolfin;

//------------------------------------------------------------------------------
DCell::DCell() : id(0), parent_id(0), vertices(0), deleted(false), nref(0)
{
}
//-----------------------------------------------------------------------------
bool DCell::has_edge(DVertex *v1, DVertex *v2)
{
  uint found = 0;
  for ( std::vector<DVertex*>::iterator it = vertices.begin() ;
        it != vertices.end(); ++it )
    if ( *it == v1 || *it == v2 )
      found++;
  return (found == 2);
}
//------------------------------------------------------------------------------
real DCell::volume() const
{
  switch(vertices.size())
  {
    case 2: // interval
      {
        // volume of interval is the length, i. e. distance of endpoints
        return vertices[0]->p.distance(vertices[1]->p);
      }
      break;
    case 3: // triangle
      {
        // formula from http://mathworld.wolfram.com/TriangleArea.html
        Point p01 = vertices[1]->p - vertices[0]->p;
        Point p02 = vertices[2]->p - vertices[0]->p;
        Point n = p01.cross(p02);
        return 0.5 * n.norm();
      }
      break;
    case 4: // tetrahedron
      {
        // formula from http://mathworld.wolfram.com/Tetrahedron.html
        Point p01 = vertices[1]->p - vertices[0]->p;
        Point p02 = vertices[2]->p - vertices[0]->p;
        Point p03 = vertices[3]->p - vertices[0]->p;
        Point n = p02.cross(p03);
        real val = p01.dot(n);
        return (1.0/6.0) * std::abs(val);
      }
      break;
    default:
      error("Volume computation for given CellType not defined / implemented!");
      break;
  }
  return 0.0;
}
//------------------------------------------------------------------------------
real DCell::diameter() const
{
  switch(vertices.size())
  {
    case 2: //interval
      {
        // diameter of interval is the length, i. e. distance of endpoints
        return vertices[0]->p.distance(vertices[1]->p);
      }
      break;
    case 3: //triangle
      {
        // formula from TriangleCell::diameter
        real a = vertices[1]->p.distance(vertices[2]->p);
        real b = vertices[0]->p.distance(vertices[2]->p);
        real c = vertices[0]->p.distance(vertices[1]->p);
        return 0.5 * a * b * c / volume();
      }
      break;
    case 4: // tetrahedron
      {
        // formula from TetrahedronCell::diameter
        real a = vertices[1]->p.distance(vertices[2]->p);
        real b = vertices[0]->p.distance(vertices[2]->p);
        real c = vertices[0]->p.distance(vertices[1]->p);
        real aa = vertices[0]->p.distance(vertices[3]->p);
        real bb = vertices[1]->p.distance(vertices[3]->p);
        real cc = vertices[2]->p.distance(vertices[3]->p);
        real la = a*aa;
        real lb = b*bb;
        real lc = c*cc;
        real s = 0.5*(la+lb+lc);
        real area = sqrt(s*(s-la)*(s-lb)*(s-lc));
        return area / ( 3.0 * volume() );
      }
      break;
    default:
      error("Diameter computation for given CellType not defined / implemented!");
      break;
  }
  return 0.0;
}
//------------------------------------------------------------------------------
uint DCell::orientation() const
{
  switch(vertices.size())
  {
    case 2: //interval
      {
        Point p01 = vertices[1]->p - vertices[0]->p;
        Point n(-p01.y(), p01.x());
        return ( n.dot(p01) < 0.0 ? 1 : 0 );
      }
      break;
    case 3: //triangle
      {
        Point p01 = vertices[1]->p - vertices[0]->p;
        Point p02 = vertices[2]->p - vertices[0]->p;
        Point n(-p01.y(), p01.x());
        return ( n.dot(p02) < 0.0 ? 1 : 0 );
      }
      break;
    case 4: // tetrahedron
      {
        Point p01 = vertices[1]->p - vertices[0]->p;
        Point p02 = vertices[2]->p - vertices[0]->p;
        Point p03 = vertices[3]->p - vertices[0]->p;
        Point n = p01.cross(p02);
        return ( n.dot(p03) < 0.0 ? 1 : 0 );
      }
      break;
    default:
      error("Orientation computation for given CellType not defined / implemented!");
      break;
  }
  return 0;
}
//------------------------------------------------------------------------------
