/*
 * Copyright (C) 2003-2008 Anders Logg.
 * Licensed under the GNU LGPL Version 2.1.
 *
 * Based on the old XMLFile, XMLObject and XMLMesh classes in DOLFIN HPC
 * Refactored by Niclas Jansson, 2020
 */

#include <DOLFINxml.h>
#include <iostream>
#include <cstring>

void DOLFINxml::load_mesh(std::string& filename) {
  std::cout << "Reading DOLFIN XML \n";
  
  xmlSAXHandler sax;
  memset( &sax, 0, sizeof(sax));

  sax.startDocument = sax_start_document;
  sax.endDocument = sax_end_document;
  sax.startElement = sax_start_element;
  sax.endElement = sax_end_element;

  if (xmlSAXUserParseFile(&sax, (void *) this, filename.c_str()) < 0)
    throw std::runtime_error( "Failed to read \"" + filename + "\"" );
}

void DOLFINxml::start_element(const xmlChar *name, const xmlChar **attrs) {
  switch(state) {
  case ROOT:
    if (xmlStrcasecmp(name, (xmlChar *) "mesh") == 0) {
      const char *ctype = parse<const char *>(name, attrs, "celltype");
      if (strcmp(ctype, "triangle") == 0) cell_type = 2;
      else if (strcmp(ctype, "tetrahedron") == 0) cell_type = 3;
      else
	throw std::runtime_error("Invalid cell type");
      gdim = parse<uint32_t>(name, attrs, "dim");
      state = MESH;
    }
    break;
  case MESH:
    if (xmlStrcasecmp(name, (xmlChar *) "vertices") == 0) {
      num_vertices = parse<uint32_t>(name, attrs, "size");
      vertices = new double[num_vertices * gdim];
      vp = &vertices[0];
      state = VERTICES;
    }
    else if (xmlStrcasecmp(name, (xmlChar *) "cells") == 0) {
      num_cells = parse<uint32_t>(name, attrs, "size");
      cells = new uint32_t[num_cells * ( cell_type + 1)];
      cp = &cells[0];
      state = CELLS;
    }
    break;
  case VERTICES:
    if (xmlStrcasecmp(name, (xmlChar *) "vertex") == 0) {
      *(vp++) = parse<double>(name, attrs, "x");
      *(vp++) = parse<double>(name, attrs, "y");
      if (gdim == 3)
	*(vp++) = parse<double>(name, attrs, "z");
      parsed_vertices++;
    }
    break;
  case CELLS:
    static const char * const vertex_attr[8] =
      { "v0", "v1", "v2", "v3", "v4", "v5", "v6", "v7" };
    for (int i = 0; i < (cell_type + 1); i++) 
      *(cp++) = parse<uint32_t>(name, attrs, vertex_attr[i]);
    parsed_cells++;
    break;
  default:
    break;
  }
}

void DOLFINxml::end_element(const xmlChar *name) {
  switch(state) {
  case MESH:
    if (xmlStrcasecmp(name, (xmlChar *) "mesh") == 0) {
      state = ROOT;
    }
    break;
  case VERTICES:
    if (xmlStrcasecmp(name, (xmlChar *) "vertices") == 0) {
      if (parsed_vertices != num_vertices)
	throw std::runtime_error("Couldn't find all vertices");
      std::cout << "Found " << parsed_vertices << " vertices\n";
      state = MESH;
    }
    break;
  case CELLS:
    if (xmlStrcasecmp(name, (xmlChar *) "cells") == 0) {
      if (parsed_vertices != num_vertices)
	throw std::runtime_error("Couldn't find all cells");
      std::cout << "Found " << parsed_cells << " cells\n";
      state = MESH;

    }
    break;
  default:
    break;
  }
}

/*
 * SAX callback functions
 */
void sax_start_document(void *ctx) { /* do nothing */ }

void sax_end_document(void *ctx) { /* do nothing */ }

void sax_start_element(void *ctx, const xmlChar *name, const xmlChar **attrs) {
  ((DOLFINxml *) ctx)->start_element(name, attrs);
}

void sax_end_element(void *ctx, const xmlChar *name) {
  ((DOLFINxml *) ctx)->end_element(name);
}
