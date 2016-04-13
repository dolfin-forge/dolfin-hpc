// Copyright (C) 2003-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2008.
// Modified by Niclas Jansson, 2015-2016.
//
// This class was cleaned up and merged with the parallel implementation which
// got some fixes in the process.
//
// First added:  2003-10-21
// Last changed: 2008-05-21

#ifndef __DOLFIN_XML_MESH_H
#define __DOLFIN_XML_MESH_H

#include <dolfin/io/XMLObject.h>

#include <dolfin/math/LinearDistribution.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/MeshFunction.h>

namespace dolfin
{

class CellType;
class Mesh;

class XMLMesh : public XMLObject
{

public:
  
  ///
  XMLMesh(Mesh& mesh);

  ///
  ~XMLMesh();

  ///
  void startElement(const xmlChar* name, const xmlChar** attrs);

  ///
  void endElement(const xmlChar* name);

  void open(std::string const& filename);

  ///
  bool close();

private:
  
  enum ParserState
  {
    ROOT,
    IN_MESH,
    IN_VERTICES,
    IN_CELLS
  };

  ///
  void clear();

  /// Parsing function
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
  CellType * cell_type_;
  LinearDistribution * vertex_dist_;
  LinearDistribution * cell_dist_;
  uint cell_count_;

  uint * vertex_owner_;
  Array<uint> cell_buffer_;
  _set<uint> nonlocal_vertices_;

};

}

#endif
