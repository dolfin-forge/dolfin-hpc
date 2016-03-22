// Copyright (C) 2003-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2008.
//
// First added:  2003-10-21
// Last changed: 2008-05-21

#ifndef __DOLFIN_XML_MESH_H
#define __DOLFIN_XML_MESH_H

#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/MeshFunction.h>
#include "XMLObject.h"

namespace dolfin
{

class Mesh;

class XMLMesh : public XMLObject
{
public:
  
  XMLMesh(Mesh& mesh);
  ~XMLMesh();

  void startElement(const xmlChar* name, const xmlChar** attrs);
  void endElement(const xmlChar* name);

  void open(std::string filename);
  bool close();

private:
  
  enum ParserState
  {
    ROOT,
    IN_MESH,
    IN_VERTICES,
    IN_CELLS
  };

  void beginMesh(const xmlChar* name, const xmlChar** attrs);
  void readVertices(const xmlChar* name, const xmlChar** attrs);
  void readCells(const xmlChar* name, const xmlChar** attrs);
  void readVertex(const xmlChar* name, const xmlChar** attrs);
  void readCell(const xmlChar* name, const xmlChar** attrs);
  void endMesh();

  Mesh& mesh_;
  ParserState state_;
  MeshEditor * editor_;

  bool parallel_;
  uint cell_count_;
  uint vertex_offset_;
  uint vertex_range_end_;
  uint cell_offset_;
  uint cell_range_end_;

  uint * local_vertices_;
  uint * shared_vertices_;
  uint num_local_vertices_;
  uint num_local_cells_;

  uint * vertex_owner_;
  Array<uint> cell_buffer_;
  _set<uint> nonlocal_vertices_;

};

}

#endif
