// Copyright (C) 2003-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2008.
//
// First added:  2003-10-21
// Last changed: 2008-05-21

#ifndef __DOLFIN_PXML_MESH_H
#define __DOLFIN_PXML_MESH_H

#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/MeshFunction.h>
#include "XMLObject.h"

namespace dolfin
{

class Mesh;

class PXMLMesh : public XMLObject
{
public:
  
  PXMLMesh(Mesh& mesh);
  ~PXMLMesh();

  void startElement(const xmlChar* name, const xmlChar** attrs);
  void endElement(const xmlChar* name);

  void open(std::string filename);
  bool close();

private:
  
  enum ParserState
  {
    OUTSIDE,
    INSIDE_MESH,
    INSIDE_VERTICES,
    INSIDE_CELLS,
    INSIDE_DATA,
    INSIDE_MESH_FUNCTION,
    INSIDE_ARRAY,
    DONE
  };

  void readMesh(const xmlChar* name, const xmlChar** attrs);
  void readVertices(const xmlChar* name, const xmlChar** attrs);
  void readCells(const xmlChar* name, const xmlChar** attrs);
  void readVertex(const xmlChar* name, const xmlChar** attrs);
  void readInterval(const xmlChar* name, const xmlChar** attrs);
  void readTriangle(const xmlChar* name, const xmlChar** attrs);
  void readTetrahedron(const xmlChar* name, const xmlChar** attrs);
  void readMeshFunction(const xmlChar* name, const xmlChar** attrs);
  void readArray(const xmlChar* name, const xmlChar** attrs);
  void readMeshEntity(const xmlChar* name, const xmlChar** attrs);
  void readArrayElement(const xmlChar* name, const xmlChar** attrs);

  void closeMesh();

  Mesh& mesh_;
  ParserState state_;
  MeshEditor * editor_;
  MeshFunction<uint> * f_;
  Array<uint> * a_;

  uint numParsedVertices_;
  uint numParsedCells_;
  uint startIndex_vert_;
  uint endIndex_vert_;
  uint startIndex_cell_;
  uint endIndex_cell_;

  uint * local_vertices_;
  uint * shared_vertices_;
  uint numLocalVertices_;
  uint numLocalCells_;

  _map<uint,bool> own_vertex_;
  _map<uint,bool> used_vertex_;
  Array<uint> cell_buffer_;
  _map<uint,bool> shared_buffer_;

};

}

#endif
