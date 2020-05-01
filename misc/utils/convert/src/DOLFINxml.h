/*
 * ----------------------------------------------------------------------------
 * "THE BEER-WARE LICENSE" (Revision 42):
 * <njansson@kth.se> wrote this file. As long as you retain this notice you
 * can do whatever you want with this stuff. If we meet some day, and you think
 * this stuff is worth it, you can buy me a beer in return Niclas Jansson
 * ----------------------------------------------------------------------------
 */

#ifndef DOLFIN_CONVERT_XML_MESH_H
#define DOLFIN_CONVERT_XML_MESH_H

#include <cstdlib>
#include <cstring>
#include <Mesh.h>
#include <libxml/xmlreader.h>

struct DOLFINxml : public Mesh
{

  DOLFINxml() : Mesh() {};
		

  void load_mesh(std::string& filename) {
    /* TODO rewrite as a SAX parser */
    xmlTextReaderPtr xml_reader = xmlNewTextReaderFilename(filename.c_str());
    if (xml_reader == NULL)
      throw std::runtime_error( "Failed to open file \"" + filename + "\"" );    

    for (int i = 0; i < 3; i++)
      xmlTextReaderRead(xml_reader);    
    parse_header(xml_reader);

    for (int i = 0; i < 2; i++)
      xmlTextReaderRead(xml_reader);
    parse_vertices(xml_reader);
    
    for (int i = 0; i < 2; i++)
      xmlTextReaderRead(xml_reader);
    parse_cells(xml_reader);

    xmlFreeTextReader(xml_reader);
  }
  
  void parse_header(xmlTextReaderPtr xml_reader) {


    
    gdim = atoi((const char *)xmlTextReaderGetAttribute(xml_reader,
							(const xmlChar *)"dim"));
    const char *ctype = (const char *)
      xmlTextReaderGetAttribute(xml_reader, (const xmlChar *) "celltype");
    
    if (strcmp(ctype, "triangle") == 0) cell_type = 2;
    else if (strcmp(ctype, "tetrahedron") == 0) cell_type = 3;

  }

  void parse_vertices(xmlTextReaderPtr xml_reader) {

    num_vertices = atoi((const char *)
			xmlTextReaderGetAttribute(xml_reader,
						  (const xmlChar *) "size"));
    vertices = new double[num_vertices * gdim];
    double *dp = &vertices[0];
    xmlTextReaderRead(xml_reader);
    xmlTextReaderRead(xml_reader);
    for (int i = 0 ; i < num_vertices;
	 xmlTextReaderRead(xml_reader),
	   xmlTextReaderRead(xml_reader), i++) {
      *(dp++) = atof((const char *)
		     xmlTextReaderGetAttribute(xml_reader,
					       (const xmlChar *) "x"));
      *(dp++) = atof((const char *)
		     xmlTextReaderGetAttribute(xml_reader,
					       (const xmlChar *) "y"));
      if ( gdim == 3)
	*(dp++) = atof((const char *)
		       xmlTextReaderGetAttribute(xml_reader,
						 (const xmlChar *) "z"));      
    }
    std::cout << "Found " << num_vertices << " vertices\n";
  }

  void parse_cells(xmlTextReaderPtr xml_reader) {

    num_cells = atoi((const char *)
		     xmlTextReaderGetAttribute(xml_reader,
					       (const xmlChar *) "size"));
    cells = new uint32_t[num_cells * ( cell_type + 1)];
    uint32_t *dp = &cells[0];
    xmlTextReaderRead(xml_reader);    
    xmlTextReaderRead(xml_reader);    
    for (int i = 0 ; i < num_cells;
	 xmlTextReaderRead(xml_reader),
	   xmlTextReaderRead(xml_reader), i++) {
      *(dp++) = atoi((const char *)
		     xmlTextReaderGetAttribute(xml_reader,
					       (const xmlChar *) "v0"));
      *(dp++) = atoi((const char *)
		     xmlTextReaderGetAttribute(xml_reader,
					       (const xmlChar *) "v1"));
      *(dp++) = atoi((const char *)
		     xmlTextReaderGetAttribute(xml_reader,
					       (const xmlChar *) "v2"));
      if (cell_type == 3)
	*(dp++) = atoi((const char *)
		       xmlTextReaderGetAttribute(xml_reader,
						 (const xmlChar *)"v3"));           
    }
    std::cout << "Found " << num_cells << " cells\n";
	
  }
};

#endif
