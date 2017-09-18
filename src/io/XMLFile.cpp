// Copyright (C) 2002-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Erik Svensson, 2003.
// Modified by Garth N. Wells, 2006.
// Modified by Ola Skavhaug, 2006.
// Modified by Magnus Vikstrom, 2007.
// Modified by Niclas Jansson, 2008.
// Modified by Aurélien Larcher, 2013.
//
// First added:  2002-12-03
// Last changed: 2013-09-13

#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_XML

#include <stdarg.h>

#include <dolfin/log/log.h>
#include <dolfin/common/constants.h>
#include <dolfin/common/Array.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/main/MPI.h>

#include <dolfin/io/XMLObject.h>
#include <dolfin/io/XMLMesh.h>
#include <dolfin/io/XMLFile.h>

#include <iomanip>

namespace dolfin
{

//-----------------------------------------------------------------------------
XMLFile::XMLFile(const std::string filename) :
    GenericFile(filename),
    header_written(false),
    mark(0)
{
  type = "XML";
  xmlObject = NULL;
}
//-----------------------------------------------------------------------------
XMLFile::~XMLFile()
{
  delete xmlObject;
}
//-----------------------------------------------------------------------------
void XMLFile::operator>>(Mesh& mesh)
{
  message(1, "Reading mesh from file %s.", filename.c_str());

  delete xmlObject;
  xmlObject = new XMLMesh(mesh);
  parseFile();
}
//-----------------------------------------------------------------------------
void XMLFile::operator<<(Mesh& mesh)
{

  if(MPI::size() == 1)
  {
    this->write();

    // Open file
    FILE *fp = openFile();

    // Get cell type
    CellType::Type cell_type = mesh.type().cellType();

    // Write mesh in XML format
    fprintf(fp, "<mesh celltype=\"%s\" dim=\"%u\">\n",
            mesh.type().str().c_str(), mesh.geometry().dim());

    fprintf(fp, "<vertices size=\"%u\">\n", mesh.size(0));

    switch (mesh.geometry().dim())
      {
      case 1:
        for (VertexIterator v(mesh); !v.end(); ++v)
        {
          fprintf(fp,
                  "<vertex index=\"%u\" x=\"%g\"/>\n", v->index(),
                  v->x()[0]);
        }
        break;
      case 2:
        for (VertexIterator v(mesh); !v.end(); ++v)
        {
          fprintf(fp,
                  "<vertex index=\"%u\" x=\"%g\" y=\"%g\"/>\n",
                  v->index(), v->x()[0], v->x()[1]);
        }
        break;
      case 3:
        for (VertexIterator v(mesh); !v.end(); ++v)
        {
          fprintf(fp,
                  "<vertex index=\"%u\" x=\"%g\" y=\"%g\" z=\"%g\" />\n",
                  v->index(), v->x()[0], v->x()[1], v->x()[2]);
        }
        break;
      default:
        error("The XML mesh file format only supports 1D, 2D and 3D meshes.");
        break;
      }

    fprintf(fp, "</vertices>\n");
    fprintf(fp, "<cells size=\"%u\">\n", mesh.num_cells());

    switch (cell_type)
      {
      case CellType::interval:
        for (CellIterator c(mesh); !c.end(); ++c)
        {
          uint* vertices = c->entities(0);
          dolfin_assert(vertices);
          fprintf(fp,
                  "<interval index=\"%u\" v0=\"%u\" v1=\"%u\"/>\n",
                  c->index(), vertices[0], vertices[1]);
        }
        break;
      case CellType::triangle:
        for (CellIterator c(mesh); !c.end(); ++c)
        {
          uint* vertices = c->entities(0);
          dolfin_assert(vertices);
          fprintf(
              fp,
              "<triangle index=\"%u\" v0=\"%u\" v1=\"%u\" v2=\"%u\"/>\n",
              c->index(), vertices[0], vertices[1], vertices[2]);
        }
        break;
      case CellType::tetrahedron:
        for (CellIterator c(mesh); !c.end(); ++c)
        {
          uint* vertices = c->entities(0);
          dolfin_assert(vertices);
          fprintf(
              fp,
              "<tetrahedron index=\"%u\" v0=\"%u\" v1=\"%u\" v2=\"%u\" v3=\"%u\"/>\n",
              c->index(), vertices[0], vertices[1], vertices[2], vertices[3]);
        }
        break;
      case CellType::quadrilateral:
        for (CellIterator c(mesh); !c.end(); ++c)
        {
          uint* vertices = c->entities(0);
          dolfin_assert(vertices);
          fprintf(
              fp,
              "<quadrilateral index=\"%u\" v0=\"%u\" v1=\"%u\" v2=\"%u\" v3=\"%u\"/>\n",
              c->index(), vertices[0], vertices[1], vertices[2], vertices[3]);
        }
        break;
      case CellType::hexahedron:
        for (CellIterator c(mesh); !c.end(); ++c)
        {
          uint* vertices = c->entities(0);
          dolfin_assert(vertices);
          fprintf(
              fp,
              "<hexahedron index=\"%u\" v0=\"%u\" v1=\"%u\" v2=\"%u\" v3=\"%u\" v4=\"%u\" v5=\"%u\" v6=\"%u\" v7=\"%u\"/>\n",
              c->index(), vertices[0], vertices[1], vertices[2], vertices[3],
                          vertices[4], vertices[5], vertices[6], vertices[7]);
        }
        break;
      default:
        error("Unknown cell type: %u.", cell_type);
        break;
      }

    fprintf(fp, "</cells>\n");
    fprintf(fp, "</mesh>\n");

    // Close file
    closeFile(fp);
  }
  else
  {
    if(!mesh.is_distributed())
    {
      error("Cannot write serial mesh in parallel");
    }

#ifdef ENABLE_MPIIO

    //-------------------------------------------------------------------------
    // Delete if file exists
    MPI_File_delete((char *) filename.c_str(), MPI_INFO_NULL );

    // Open file
    MPI_File fh;
    MPI_Offset curr_offset = 0;
    MPI_File_open(dolfin::MPI::DOLFIN_COMM, (char *) filename.c_str(),
                  MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
    MPI_File_set_view(fh, 0, MPI_CHAR, MPI_CHAR, (char *) "native", MPI_INFO_NULL );

    // Write DOLFIN XML format header
    std::string const hdr =
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
        "<dolfin xmlns:dolfin=\"http://www.fenics.org/dolfin/\">\n";
    MPI_File_write_all(fh, (void *) hdr.c_str(), hdr.size(), MPI_CHAR,
                       MPI_STATUS_IGNORE);
    curr_offset += hdr.size();
    //-------------------------------------------------------------------------
    // Mesh properties
    uint const tdim = mesh.topology().dim();
    uint const gdim = mesh.geometry().dim();
    CellType::Type cell_type = mesh.type().cellType();
    std::string const cell_str(mesh.type().str());
    MeshDistributedData const& distdata = mesh.distdata();
    uint const num_owned_vertices = distdata[0].num_owned();
    uint const num_global_vertices = distdata[0].global_size();
    uint const num_owned_cells = distdata[tdim].num_owned();
    uint const num_global_cells = distdata[tdim].global_size();

    // Open mesh and vertices
    std::stringstream s0;
    s0
    << "<mesh celltype=\"" << cell_str << "\" " << "dim=\"" << gdim << "\">\n"
    << "<vertices size=\"" << num_global_vertices << "\">\n";
    MPI_File_write_all(fh, (void *) s0.str().c_str(), s0.str().size(), 
		       MPI_CHAR, MPI_STATUS_IGNORE);
    curr_offset += s0.str().size();

    // Determine maximum size of line to allocate buffer
    uint const max_digits = 24; // %g takes the shortest.
    std::stringstream svline;
    svline << "<vertex index=\"" << num_global_vertices << "\"";
    for(uint i = 0; i < gdim; ++i)
    {
      // Conservative
      svline << " x=\"" << std::setw(max_digits) << "0" << "\"";
    }
    svline << "/>\n";

    // Fill buffer
    uint vsize = 0;
    uint vline = 0;
    uint vline_max = svline.str().size();
    uint const vsize_max = vline_max * num_owned_vertices;
    char * vbuffer = new char[vsize_max];
    switch (gdim)
      {
      case 1:
        for (VertexIterator v(mesh); !v.end(); ++v)
        {
          if(!v->is_ghost())
          {
            vline = sprintf(vbuffer + vsize,
              "<vertex index=\"%u\" x=\"%g\"/>\n",
              v->global_index(), v->x()[0]);
            dolfin_assert(vline <= vline_max);
            vsize += vline;
          }
        }
        break;
      case 2:
        for (VertexIterator v(mesh); !v.end(); ++v)
        {
          if(!v->is_ghost())
          {
            vline = sprintf(vbuffer + vsize,
              "<vertex index=\"%u\" x=\"%g\" y=\"%g\"/>\n",
              v->global_index(), v->x()[0], v->x()[1]);
            dolfin_assert(vline <= vline_max);
            vsize += vline;
          }
        }
        break;
      case 3:
        for (VertexIterator v(mesh); !v.end(); ++v)
        {
          if(!v->is_ghost())
          {
            vline = sprintf(vbuffer + vsize,
              "<vertex index=\"%u\" x=\"%g\" y=\"%g\" z=\"%g\" />\n",
              v->global_index(), v->x()[0], v->x()[1], v->x()[2]);
            dolfin_assert(vline <= vline_max);
            vsize += vline;
          }
        }
        break;
      default:
        error("The XML mesh file format only supports 1D, 2D and 3D meshes.");
        break;
      }

    // Write buffer
    if(vsize > vsize_max)
    {
      error("Buffer overflow writing vertices: %d > %d", vsize, vsize_max);
    }
    uint voffset = 0;
#if ( MPI_VERSION > 1 )
    MPI_Exscan(&vsize, &voffset, 1, MPI_UNSIGNED, MPI_SUM,
               dolfin::MPI::DOLFIN_COMM);
#else
    MPI_Scan(&vsize, &voffset, 1, MPI_UNSIGNED, MPI_SUM,
	     dolfin::MPI::DOLFIN_COMM);
    voffset -= vsize;
#endif
    MPI_File_write_at_all(fh, curr_offset + voffset, vbuffer, vsize, MPI_CHAR,
                          MPI_STATUS_IGNORE);
    MPI_Allreduce(&vsize, &voffset, 1, MPI_UNSIGNED, MPI_SUM,
                  dolfin::MPI::DOLFIN_COMM);
    curr_offset += voffset;
    delete [] vbuffer;

    // Close vertices and open cells
    std::stringstream s1;
    s1
    << "</vertices>\n"
    << "<cells size=\"" << num_global_cells << "\">\n";
    MPI_File_write_at_all(fh, curr_offset, (void *) s1.str().c_str(), 
			  s1.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
    curr_offset += s1.str().size();

    // Determine maximum size of line to allocate buffer
    std::stringstream scline;
    scline << "<" << cell_str << " index=\""<< num_global_cells << "\"";
    for (uint i = 0; i < mesh.type().num_entities(0); ++i)
    {
      scline << "<" << cell_str
             << " v" << i << "=\""<< num_global_vertices << "\"";
    }
    scline << "/>\n";

    // Fill buffer
    uint csize = 0;
    uint cline = 0;
    uint const cline_max = scline.str().size();
    uint const csize_max = cline_max * num_owned_cells;
    char * cbuffer = new char[csize_max];
    switch (cell_type)
      {
      case CellType::interval:
        for (CellIterator c(mesh); !c.end(); ++c)
        {
          uint* vertices = c->entities(0);
          dolfin_assert(vertices);
          cline = sprintf(cbuffer + csize,
            "<interval index=\"%u\" v0=\"%u\" v1=\"%u\"/>\n",
            c->global_index(),
            distdata[0].get_global(vertices[0]),
            distdata[0].get_global(vertices[1]));
          dolfin_assert(cline <= cline_max);
          csize += cline;
        }
        break;
      case CellType::triangle:
        for (CellIterator c(mesh); !c.end(); ++c)
        {
          uint* vertices = c->entities(0);
          dolfin_assert(vertices);
          cline = sprintf(cbuffer + csize,
            "<triangle index=\"%u\" v0=\"%u\" v1=\"%u\" v2=\"%u\"/>\n",
            c->global_index(),
            distdata[0].get_global(vertices[0]),
            distdata[0].get_global(vertices[1]),
            distdata[0].get_global(vertices[2]));
          dolfin_assert(cline <= cline_max);
          csize += cline;
        }
        break;
      case CellType::tetrahedron:
        for (CellIterator c(mesh); !c.end(); ++c)
        {
          uint* vertices = c->entities(0);
          dolfin_assert(vertices);
          cline = sprintf(cbuffer + csize,
            "<tetrahedron index=\"%u\" v0=\"%u\" v1=\"%u\" v2=\"%u\" v3=\"%u\"/>\n",
            c->global_index(),
            distdata[0].get_global(vertices[0]),
            distdata[0].get_global(vertices[1]),
            distdata[0].get_global(vertices[2]),
            distdata[0].get_global(vertices[3]));
          dolfin_assert(cline <= cline_max);
          csize += cline;
        }
        break;

      case CellType::quadrilateral:
        for (CellIterator c(mesh); !c.end(); ++c)
        {
          uint* vertices = c->entities(0);
          dolfin_assert(vertices);
          cline = sprintf(cbuffer + csize,
            "<quadrilateral index=\"%u\" v0=\"%u\" v1=\"%u\" v2=\"%u\" v3=\"%u\"/>\n",
            c->global_index(),
            distdata[0].get_global(vertices[0]),
            distdata[0].get_global(vertices[1]),
            distdata[0].get_global(vertices[2]),
            distdata[0].get_global(vertices[3]));
          dolfin_assert(cline <= cline_max);
          csize += cline;
        }
        break;
      case CellType::hexahedron:
        for (CellIterator c(mesh); !c.end(); ++c)
        {
          uint* vertices = c->entities(0);
          dolfin_assert(vertices);
          cline = sprintf(cbuffer + csize,
            "<hexahedron index=\"%u\" v0=\"%u\" v1=\"%u\" v2=\"%u\" v3=\"%u\" v4=\"%u\" v5=\"%u\" v6=\"%u\" v7=\"%u\"/>\n",
            c->global_index(),
            distdata[0].get_global(vertices[0]),
            distdata[0].get_global(vertices[1]),
            distdata[0].get_global(vertices[2]),
            distdata[0].get_global(vertices[3]),
            distdata[0].get_global(vertices[4]),
            distdata[0].get_global(vertices[5]),
            distdata[0].get_global(vertices[6]),
            distdata[0].get_global(vertices[7]));
          dolfin_assert(cline <= cline_max);
          csize += cline;
        }
        break;
      default:
        error("Unknown cell type: %u.", cell_type);
        break;
      }

    // Write buffer
    if(csize > csize_max)
    {
      error("Buffer overflow writing cells: %d > %d", csize, csize_max);
    }
    uint coffset = 0;
#if ( MPI_VERSION > 1 )
    MPI_Exscan(&csize, &coffset, 1, MPI_UNSIGNED, MPI_SUM,
               dolfin::MPI::DOLFIN_COMM);
#else
    MPI_Scan(&csize, &coffset, 1, MPI_UNSIGNED, MPI_SUM,
	     dolfin::MPI::DOLFIN_COMM);
    coffset -= csize;
#endif
    MPI_File_write_at_all(fh, curr_offset + coffset, cbuffer, csize, MPI_CHAR,
                          MPI_STATUS_IGNORE);
    MPI_Allreduce(&csize, &coffset, 1, MPI_UNSIGNED, MPI_SUM,
                  dolfin::MPI::DOLFIN_COMM);
    curr_offset += coffset;
    delete [] cbuffer;

    // Close cells and mesh
    std::string const s2("</cells>\n</mesh>\n");
    MPI_File_write_at_all(fh, curr_offset, (void *) s2.c_str(), s2.size(), 
			  MPI_CHAR, MPI_STATUS_IGNORE);
    curr_offset += s2.size();

    //-------------------------------------------------------------------------
    // Write DOLFIN XML format footer
    std::string const ftr("</dolfin>\n");
    MPI_File_write_at_all(fh, curr_offset, (void *) ftr.c_str(), ftr.size(), 
			  MPI_CHAR, MPI_STATUS_IGNORE);

    // Close file
    MPI_File_close(&fh);

#else
    parallel_write_not_impl("Mesh");
#endif
  }

  message(1, "Saved mesh to file %s in DOLFIN XML format.", filename.c_str());
}
//-----------------------------------------------------------------------------
FILE* XMLFile::openFile()
{

  // Open file
  FILE *fp = fopen(filename.c_str(), "r+");
  if(fp == NULL)
  {
    error("XMLFile : impossible to open file '%s'", filename.c_str());
  }

  // Step to position before previously written footer
  //printf("Stepping to position: %ld\n", mark);
  fseek(fp, mark, SEEK_SET);
  fflush(fp);

  // Write DOLFIN XML format header
  if (!header_written)
  {
    fprintf(fp, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n\n");
    fprintf(fp, "<dolfin xmlns:dolfin=\"http://www.fenics.org/dolfin/\">\n");

    header_written = true;
  }

  return fp;
}
//-----------------------------------------------------------------------------
void XMLFile::closeFile(FILE* fp)
{
  // Get position in file before writing footer
  mark = ftell(fp);
  //printf("Position in file before writing footer: %ld\n", mark);

  // Write DOLFIN XML format footer
  if (header_written) fprintf(fp, "</dolfin>\n");

  // Close file
  fclose(fp);
}
//-----------------------------------------------------------------------------
void XMLFile::parseFile()
{
  // Notify that file is being opened
  xmlObject->open(filename);

  // Parse file using the SAX interface
  parseSAX();

  // Notify that file is being closed
  if (!xmlObject->close()) error("Unable to find data in XML file.");
}
//-----------------------------------------------------------------------------
void XMLFile::parseSAX()
{
  // Set up the sax handler. Note that it is important that we initialise
  // all (24) fields, even the ones we don't use!
  xmlSAXHandler sax =
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

  // Set up handlers for parser events
  sax.startDocument = sax_start_document;
  sax.endDocument = sax_end_document;
  sax.startElement = sax_start_element;
  sax.endElement = sax_end_element;
  sax.warning = sax_warning;
  sax.error = sax_error;
  sax.fatalError = sax_fatal_error;

  // Parse file
  xmlSAXUserParseFile(&sax, (void *) xmlObject, filename.c_str());
}
//-----------------------------------------------------------------------------
// Callback functions for the SAX interface
//-----------------------------------------------------------------------------
void sax_start_document(void *ctx)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void sax_end_document(void *ctx)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void sax_start_element(void *ctx, const xmlChar *name,
                               const xmlChar **attrs)
{
  ((XMLObject *) ctx)->startElement(name, attrs);
}
//-----------------------------------------------------------------------------
void sax_end_element(void *ctx, const xmlChar *name)
{
  ((XMLObject *) ctx)->endElement(name);
}
//-----------------------------------------------------------------------------
void sax_warning(void *ctx, const char *msg, ...)
{
  va_list args;
  va_start(args, msg);
  char buffer[DOLFIN_LINELENGTH];
  vsnprintf(buffer, DOLFIN_LINELENGTH, msg, args);
  warning("Incomplete XML data: " + std::string(buffer));
  va_end(args);
}
//-----------------------------------------------------------------------------
void sax_error(void *ctx, const char *msg, ...)
{
  va_list args;
  va_start(args, msg);
  char buffer[DOLFIN_LINELENGTH];
  vsnprintf(buffer, DOLFIN_LINELENGTH, msg, args);
  error("Illegal XML data: " + std::string(buffer));
  va_end(args);
}
//-----------------------------------------------------------------------------
void sax_fatal_error(void *ctx, const char *msg, ...)
{
  va_list args;
  va_start(args, msg);
  char buffer[DOLFIN_LINELENGTH];
  vsnprintf(buffer, DOLFIN_LINELENGTH, msg, args);
  error("Illegal XML data: " + std::string(buffer));
  va_end(args);
}
//-----------------------------------------------------------------------------
}

#endif



