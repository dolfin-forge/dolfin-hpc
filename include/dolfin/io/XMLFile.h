// Copyright (C) 2003-2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

// This class is still an abstraction nightmare and the code could be sexier...

#ifndef __DOLFIN_XML_FILE_H
#define __DOLFIN_XML_FILE_H

#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_XML
#include <libxml/parser.h>
#endif

#include <dolfin/common/types.h>
#include <dolfin/la/Vector.h>
#include <dolfin/la/GenericMatrix.h>
#include <dolfin/io/GenericFile.h>
#include <dolfin/io/XMLMeshFunction.h>

namespace dolfin
{

  class Mesh;
  template <typename T> class MeshFunction;
  class ParameterList;
  class XMLObject;

  class XMLFile : public GenericFile
  {
  public:

    XMLFile(const std::string filename);
    ~XMLFile();

    // Input
    void operator>> (Mesh& mesh);

    // Output
    void operator<< (Mesh& mesh);

    virtual void operator>> (MeshFunction<int>& meshfunction);
    virtual void operator>> (MeshFunction<uint>& meshfunction);
    virtual void operator>> (MeshFunction<real>& meshfunction);
    virtual void operator>> (MeshFunction<bool>& meshfunction);

    // Friends
    #ifdef HAVE_XML
    friend void sax_start_element (void *ctx, const xmlChar *name, const xmlChar **attrs);
    friend void sax_end_element   (void *ctx, const xmlChar *name);
    #endif

  private:

    void parseFile();
    void parseSAX();

    FILE* openFile();
    void  closeFile(FILE* fp);

    template<typename T>
    void read_meshfunction( MeshFunction<T> & meshfunction );

    // Implementation for specific class (output)
    XMLObject* xmlObject;

    // True if header is written (need to close)
    bool header_written;

    // Most recent position in file
    long mark;

  };

  // Callback functions for the SAX interface
#ifdef HAVE_XML
  void sax_start_document (void *ctx);
  void sax_end_document   (void *ctx);
  void sax_start_element  (void *ctx, const xmlChar *name, const xmlChar **attrs);
  void sax_end_element    (void *ctx, const xmlChar *name);

  void sax_warning     (void *ctx, const char *msg, ...);
  void sax_error       (void *ctx, const char *msg, ...);
  void sax_fatal_error (void *ctx, const char *msg, ...);
#endif

//-----------------------------------------------------------------------------
template<typename T>
void XMLFile::read_meshfunction( MeshFunction<T> & meshfunction )
{
  message(1, "Reading meshfunction from file %s.", filename.c_str());

  delete xmlObject;
  xmlObject = new XMLMeshFunction<T>(meshfunction);
  parseFile();
}

}
#endif
