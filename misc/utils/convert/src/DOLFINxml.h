/*
 * Copyright (C) 2003-2008 Anders Logg.
 * Licensed under the GNU LGPL Version 2.1.
 *
 * Based on the old XMLFile, XMLObject and XMLMesh classes in DOLFIN HPC
 * Refactored by Niclas Jansson, 2020
 */

#ifndef DOLFIN_CONVERT_XML_MESH_H
#define DOLFIN_CONVERT_XML_MESH_H

#include <Mesh.h>
#include <libxml/parser.h>
#include <stdexcept>

class DOLFINxml : public Mesh
{
public:
  DOLFINxml() : Mesh(),
		state(ROOT),
		vp(NULL),
		cp(NULL),
		parsed_vertices(0),
		parsed_cells(0) {}

  void load_mesh(std::string& filename);

  void start_element(const xmlChar *name, const xmlChar **atrs);
  void end_element(const xmlChar *name);

private:

  enum sax_state {ROOT, MESH, VERTICES, CELLS};
  sax_state state;

  double *vp;
  uint32_t *cp;

  uint32_t parsed_vertices;
  uint32_t parsed_cells;

  template <class T> T read(const xmlChar * s) {return T();};

  template <>
   const char * read(const xmlChar *s) {
    return reinterpret_cast<const char *>(s);
  }

  template <>
  double read(const xmlChar *s) {
    return strtod(reinterpret_cast<const char *>(s), NULL);
  }

  template <>
  uint32_t read(const xmlChar *s) {
    return (uint32_t) atoi(reinterpret_cast<const char *>(s));
  }

  template <class T>
  T parse(const xmlChar* name, const xmlChar** attrs, const char * label) {
     if (attrs == NULL) {
       throw std::runtime_error("Missing attribute");
     }
     for (uint32_t i = 0; attrs[i] != NULL; ++i) {
       if (xmlStrcasecmp(attrs[i], (xmlChar *) label) == 0) {
	 if (attrs[i + 1] == NULL) {
	   throw std::runtime_error("Missing value for an attribute");
	 }
	 return read<T>(attrs[i + 1]);
       }
     }
     throw std::runtime_error("Missing attribute");
     return T();
  }
};

/*
 * SAX callback functions
 */
void sax_start_document(void *ctx);
void sax_end_document(void *ctx);
void sax_start_element(void *ctx, const xmlChar *name, const xmlChar **attrs);
void sax_end_element(void *ctx, const xmlChar *name);

#endif
