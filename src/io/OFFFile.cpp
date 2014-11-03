// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/config/dolfin_config.h>

#include <dolfin/io/OFFFile.h>
#include <dolfin/mesh/CellType.h>
#include <dolfin/mesh/MeshEditor.h>

#include <algorithm>
#include <iostream>
#include <fstream>
#include <set>
#include <sstream>

namespace dolfin
{

//-----------------------------------------------------------------------------
OFFFile::OFFFile(std::string const filename) :
    GenericFile(filename)
{
  type = "OFF";
}
//-----------------------------------------------------------------------------
OFFFile::~OFFFile()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void OFFFile::operator>>(Mesh& mesh)
{
  std::ifstream off("2O2BWXQL_001_timewise_ip4msct.off");
  if (!off.is_open())
  {
    dolfin::error("Unable to open file.");
  }

  std::string line;

  // First line is OFF1 file type
  get_next_line(off, line);
  line.erase(remove_if(line.begin(), line.end(), isspace), line.end());
  if (line != "OFF1")
  {
    dolfin::error("File format is not Geomview Object File Version 1.");
  }

  // Get number of entities
  // 0: vertex
  // 1: faces
  // 2: lines (ignored)
  std::vector<uint> num_entities;
  get_next_line(off, line);
  split_line(line, num_entities);

  // Create mesh editor
  // Since this is a surface mesh the topology dimension is 2
  // while the geometric dimension is 3
  MeshEditor editor(mesh, CellType::triangle, 2, 3);
  editor.initVertices(num_entities[0]);
  editor.initCells(num_entities[1]);

  // Get vertices
  for (uint v = 0; v < num_entities[0]; ++v)
  {
    get_next_line(off, line);
    std::vector<real> coords;
    split_line(line, coords);
    editor.addVertex(v, &coords[0]);
  }

  // Get face connectivities
  for (uint c = 0; c < num_entities[1]; ++c)
  {
    get_next_line(off, line);
    std::vector<uint> v;

    // A line contains the type + connectivities
    split_line(line, v);

    switch(v[0])
    {
      case 3: // Only treat triangles
        editor.addCell(c, &v[1]);
        break;
      default:
        break;
    }
  }

  // Close mesh editor
  editor.close();

  // Close filestream
  off.close();
}
//-----------------------------------------------------------------------------
void OFFFile::get_next_line(std::ifstream& filestream, std::string& line)
{
  // Skip leading comments
   while (true)
   {
     std::getline(filestream, line);
     if (line[0] != '#')
     {
       break;
     }
   }
}

}
