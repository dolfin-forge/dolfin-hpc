// Copyright (C) 2005-2007 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Anders Logg 2005-2006.
// Modified by Kristian Oelgaard 2006.
// Modified by Niclas Jansson 2008-2009.
// Modified by Kaspar Mueller 2013.
// Modified by Aurélien Larcher 2014.
//
// First added:  2005-07-05
// Last changed: 2014-02-08

#include <dolfin/io/VTKFile.h>

#include <dolfin/fem/FiniteElementSpace.h>
#include <dolfin/function/Function.h>
#include <dolfin/io/Encoder.h>
#include <dolfin/la/Vector.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/Vertex.h>

#include <stdint.h>

namespace dolfin
{

//----------------------------------------------------------------------------
VTKFile::VTKFile(const std::string filename) : GenericFile(filename), _t(0)
{
  type = "VTK";
}
//----------------------------------------------------------------------------
VTKFile::VTKFile(const std::string filename, real& t) :
  GenericFile(filename), _t(&t)
{
  type = "VTK";
}
//----------------------------------------------------------------------------
VTKFile::~VTKFile()
{
  // Do nothing
}
//----------------------------------------------------------------------------
void VTKFile::operator<<(Mesh& mesh)
{
  // Update vtu file name and clear file
  vtuNameUpdate(counter);

  // Only the root updates the pvd file
  if(MPI::processNumber() == 0)
  {

    if (MPI::numProcesses() > 1)
    {
      // Update pvtu file name and clear file
      pvtuNameUpdate(counter);

      // Write pvtu file
      pvtuFileWrite();
    }

    // Write pvd file
    pvdFileWrite(counter);
  }

  // Write headers
  VTKHeaderOpen(mesh);

  // Write mesh
  MeshWrite(mesh);

  // Close headers
  VTKHeaderClose();

  // Increase the number of times we have saved the mesh
  ++counter;

//  dolfin_debug(message(1, "Saved mesh %s (%s) to file %s in VTK format.",
//          mesh.name().c_str(), mesh.label().c_str(), filename.c_str());)
}
//----------------------------------------------------------------------------
void VTKFile::operator<<(MeshFunction<int>& meshfunction)
{
  MeshFunctionWrite(meshfunction);
}
//----------------------------------------------------------------------------
void VTKFile::operator<<(MeshFunction<unsigned int>& meshfunction)
{
  MeshFunctionWrite(meshfunction);
}
//----------------------------------------------------------------------------
void VTKFile::operator<<(MeshFunction<double>& meshfunction)
{
  MeshFunctionWrite(meshfunction);
}
//----------------------------------------------------------------------------
void VTKFile::operator<<(MeshFunction<bool>& meshfunction)
{
  MeshFunctionWrite(meshfunction);
}
//----------------------------------------------------------------------------
void VTKFile::operator<<(Function& u)
{
  std::pair<Function*, std::string> f(&u, u.name());
  std::vector<std::pair<Function*, std::string> > tmp;
  tmp.push_back(f);
  write_dataset(tmp);
}
//----------------------------------------------------------------------------
void VTKFile::operator<<(std::vector<std::pair<Function*, std::string> >& f)
{
  write_dataset(f);
}
//----------------------------------------------------------------------------
void VTKFile::write_dataset(std::vector<std::pair<Function*, std::string> >& f)
{
  dolfin_assert(f.size() > 0);

  // Update vtu file name and clear file
  vtuNameUpdate(counter);

  // Write pvd file

  // Only the root updates the pvd file
  if(MPI::processNumber() == 0)
  {

    if(MPI::numProcesses() > 1)
    {
      // Update pvtu file name and clear file
      pvtuNameUpdate(counter);

      // Write pvtu file
      pvtuFileWrite_func(f);
    }

    // Write pvd file
    pvdFileWrite(counter);
  }

  Mesh& mesh = f[0].first->mesh();

  // Write headers
  VTKHeaderOpen(mesh);

  // Write Mesh
  MeshWrite(mesh);

  // Write results
  ResultsWrite(f);

  // Close headers
  VTKHeaderClose();

  // Increase the number of times we have saved the function
  ++counter;

}
//----------------------------------------------------------------------------
void VTKFile::MeshWrite(Mesh& mesh) const
{
  // Open file
  FILE* fp = fopen(vtu_filename.c_str(), "a");

  uint const num_mesh_verts = mesh.numVertices();
  uint const num_mesh_cells = mesh.numCells();
  uint const num_cell_verts = mesh.type().numEntities(0);
  uint const cell_verts_block_size = num_cell_verts * num_mesh_cells;

  // Write vertex positions
  fprintf(fp, "<Points>  \n");
  fprintf(fp, "<DataArray  type=\"Float32\"  NumberOfComponents=\"3\"  format=\"binary\"> \n");
  std::vector<float> v_data(3*num_mesh_verts);
  std::vector<float>::iterator entry = v_data.begin();
  for (VertexIterator v(mesh); !v.end(); ++v)
  {
    Point p = v->point();
    *entry++ = p.x();
    *entry++ = p.y();
    *entry++ = p.z();
  }

  // Create encoded stream
  std::stringstream base64_stream;
  encode_stream(base64_stream, v_data);
  fprintf(fp, "%s\n", base64_stream.str().c_str());
  fprintf(fp, "</DataArray>  \n");
  fprintf(fp, "</Points>  \n");

  // Write cell connectivity
  fprintf(fp, "<Cells>  \n");
  fprintf(fp, "<DataArray  type=\"Int32\"  Name=\"connectivity\"  format=\"binary\"> \n");
  std::vector<uint32_t> c_data(cell_verts_block_size);
  std::vector<uint32_t>::iterator c_entry = c_data.begin();
  for (CellIterator c(mesh); !c.end(); ++c)
  {
    for (VertexIterator v(*c); !v.end(); ++v)
    {
      *c_entry++ = v->index();
    }
  }

  // Create encoded stream
  std::stringstream base64_c_stream;
  encode_stream(base64_c_stream, c_data);
  fprintf(fp, "%s\n", base64_c_stream.str().c_str());
  fprintf(fp, "</DataArray> \n");

  // Write offset into connectivity array for the end of each cell
  fprintf(fp, "<DataArray  type=\"Int32\"  Name=\"offsets\"  format=\"binary\">  \n");
  std::vector<uint32_t>::iterator cc_entry = c_data.begin();

  for (uint offsets = 1; offsets <= num_mesh_cells; ++offsets)
  {
    *cc_entry++ = offsets * num_cell_verts;
  }

  std::stringstream base64_cc_stream;
  encode_stream(base64_cc_stream, c_data);
  fprintf(fp, "%s\n",  base64_cc_stream.str().c_str());
  fprintf(fp, "</DataArray> \n");

  //Write cell type
  fprintf(fp, "<DataArray  type=\"UInt8\"  Name=\"types\"  format=\"binary\">  \n");

  uint8_t typeval = 0;
  switch(mesh.type().cellType())
  {
    case CellType::tetrahedron:
      typeval = uint8_t(10);
      break;
    case CellType::triangle:
      typeval = uint8_t(5);
      break;
    case CellType::interval:
      typeval = uint8_t(3);
      break;
    default:
      error("Provided mesh type is not supported");
      break;
  }
  std::vector<uint8_t> t_data(num_mesh_cells, typeval);

  // Create encoded stream
  std::stringstream base64_t_stream;
  encode_stream(base64_t_stream, t_data);
  fprintf(fp, "%s\n", base64_t_stream.str().c_str());
  fprintf(fp, "</DataArray> \n");
  fprintf(fp, "</Cells> \n");

  // Close file
  fclose(fp);
}
//----------------------------------------------------------------------------
void VTKFile::ResultsWrite(std::vector<std::pair<Function*, std::string> > f) const
{
  // Open file
  FILE *fp = fopen(vtu_filename.c_str(), "a");

  // Assume same mesh for all data arrays
  Mesh& mesh = f[0].first->mesh();

  //--- Interpolation to vertices for general functions-----------------------
  fprintf(fp, "<PointData> \n");
  for (std::vector<std::pair<Function*, std::string> >::iterator it = f.begin();
       it != f.end(); it++)
  {
    Function* u = it->first;

    // Check type of function space
    if(u->type() == Function::discrete && u->space().is_cellwise_constant())
    {
      // These are cell based function and will be written as CellData
      continue;
    }

    //
    uint const value_rank = u->rank();
    uint const value_dim = u->dim(0);
    if ( value_dim > 3 )
    {
      error("Cannot handle VTK file with number of components > 3.");
    }

    // Write function data at mesh vertices
    std::string& name = it->second;
    switch(value_rank)
    {
      case 0:
        fprintf(fp, "<DataArray  type=\"Float32\"  Name=\"%s\"  format=\"binary\">   ", name.c_str());
        break;
      case 1:
        fprintf(fp, "<DataArray  type=\"Float32\"  Name=\"%s\"  NumberOfComponents=\"3\" format=\"binary\">  ", name.c_str());
        break;
      default:
        error("Only scalar and vectors functions can be saved in VTK format.");
        break;
    }

    // Allocate memory for function values at vertices
    uint const num_verts = mesh.numVertices();
    uint size = 1;
    for (uint i = 0; i < value_rank; i++)
    {
      size *= u->dim(i);
    }
    real* values = new real[size*num_verts];

    // Get function values at vertices
    u->interpolate_vertex_values(values);

    // Copy values
    // Should be value_dim^value_rank but 1^0 = 1 and (2|3)^1 = (2|3) so ...
    std::vector<float> data;

    switch(value_dim)
    {
    case 1:
      data.resize(num_verts);
      {
        std::vector<float>::iterator entry = data.begin();
        for (VertexIterator vertex(mesh); !vertex.end(); ++vertex)
        {
          *entry++ = values[ vertex->index() ];
        }
      }
      break;
    case 2:
      data.resize(3*num_verts);
      {
        std::vector<float>::iterator entry = data.begin();
        for (VertexIterator vertex(mesh); !vertex.end(); ++vertex)
        {
          *entry++ = values[ vertex->index() ];
          *entry++ = values[ vertex->index() + num_verts];
          *entry++ = 0.0;
        }
      }
      break;
    case 3:
    default:
      data.resize(3*num_verts);
      {
        std::vector<float>::iterator entry = data.begin();
        for (VertexIterator vertex(mesh); !vertex.end(); ++vertex)
        {
          *entry++ = values[ vertex->index() ];
          *entry++ = values[ vertex->index() + num_verts];
          *entry++ = values[ vertex->index() + 2*num_verts];
        }
      }
      break;
    }

    //
    delete [] values;

    // Create encoded stream
    std::stringstream base64_stream;
    encode_stream(base64_stream, data);
    fprintf(fp, "%s\n", base64_stream.str().c_str());
    fprintf(fp, "</DataArray> \n");

  }
  fprintf(fp, "</PointData> \n");

  //--- Cellwise functions----------------------------------------------------
  fprintf(fp, "<CellData> \n");
  for (std::vector<std::pair<Function*, std::string> >::iterator it = f.begin();
       it != f.end(); it++)
  {
    Function* u = it->first;

    // Check type of function space
    if(u->type() != Function::discrete || !u->space().is_cellwise_constant())
    {
      // These are not cell based functions and will be written as PointData
      continue;
    }

    //
    uint const value_rank = u->rank();
    uint const value_dim = u->dim(0);
    if ( value_dim > 3 )
    {
      error("Cannot handle VTK file with number of components > 3.");
    }

    //
    DofMap const& dm_u   = u->dofmap();
    uint nbdofspercell_u    = dm_u.local_dimension();

    // Write function data in cell
    std::string& name = it->second;
    std::vector<float> data;
    switch(value_rank)
    {
      case 0:
        fprintf(fp, "<DataArray  type=\"Float32\"  Name=\"%s\"  format=\"binary\">   ", name.c_str());
        break;
      case 1:
        fprintf(fp, "<DataArray  type=\"Float32\"  Name=\"%s\"  NumberOfComponents=\"3\" format=\"binary\">  ", name.c_str());
        break;
      default:
        error("Only scalar and vectors functions can be saved in VTK format.");
        break;
    }

    // Copy values
    // Should be value_dim^value_rank but 1^0 = 1 and (2|3)^1 = (2|3) so ...
    real * values = new real[dm_u.dofsmapping_size()];
    u->get_block(values);

    switch(value_dim)
    {
    case 1:
      data.resize(mesh.numCells());
      {
        std::vector<float>::iterator entry = data.begin();
        uint ii = 0;
        for (CellIterator cell(mesh); !cell.end(); ++cell, ++ii)
        {
          *entry++ = values[ cell->index() ];
        }
      }
      break;
    case 2:
      data.resize(3*mesh.numCells());
      {
        std::vector<float>::iterator entry = data.begin();
        uint ii = 0;
        for (CellIterator cell(mesh); !cell.end(); ++cell, ++ii)
        {
          *entry++ = values[ cell->index()*nbdofspercell_u ];
          *entry++ = values[ cell->index()*nbdofspercell_u + 1 ];
          *entry++ = 0.0;
        }
      }
      break;
    case 3:
    default:
      data.resize(3*mesh.numCells());
      {
        std::vector<float>::iterator entry = data.begin();
        uint ii = 0;
        for (CellIterator cell(mesh); !cell.end(); ++cell, ++ii)
        {
          *entry++ = values[ cell->index()*nbdofspercell_u ];
          *entry++ = values[ cell->index()*nbdofspercell_u + 1 ];
          *entry++ = values[ cell->index()*nbdofspercell_u + 2 ];
        }
      }
      break;
    }

    //
    delete [] values;

    // Create encoded stream
    std::stringstream base64_stream;
    encode_stream(base64_stream, data);
    fprintf(fp, "%s\n", base64_stream.str().c_str());
    fprintf(fp, "</DataArray> \n");
  }
  fprintf(fp, "</CellData> \n");


  // Close file
  fclose(fp);

}
//----------------------------------------------------------------------------
void VTKFile::pvdFileWrite(uint num)
{
  std::fstream pvdFile;

  if( num == 0)
  {
    // Open pvd file
    pvdFile.open(filename.c_str(), std::ios::out|std::ios::trunc);
    // Write header
    pvdFile << "<?xml version=\"1.0\"?> " << std::endl;
    pvdFile << "<VTKFile type=\"Collection\" version=\"0.1\" > " << std::endl;
    pvdFile << "<Collection> " << std::endl;
  }
  else
  {
    // Open pvd file
    pvdFile.open(filename.c_str(),  std::ios::out|std::ios::in);
    pvdFile.seekp(mark);
  }

  // Remove directory path from name for pvd file
  std::string fname;
  if(MPI::numProcesses() > 1)
  {
    fname.assign(pvtu_filename, filename.find_last_of("/") + 1, pvtu_filename.size());
  }
  else
  {
    fname.assign(vtu_filename, filename.find_last_of("/") + 1, vtu_filename.size());
  }

  // Data file name
  if(_t)
  {
    pvdFile << "<DataSet timestep=\"" << *_t << "\" part=\"0\"" << " file=\""
            <<  fname <<  "\"/>" << std::endl;
  }
  else
  {
    pvdFile << "<DataSet timestep=\"" << num << "\" part=\"0\"" << " file=\""
            <<  fname <<  "\"/>" << std::endl;
  }
  mark = pvdFile.tellp();

  // Close headers
  pvdFile << "</Collection> " << std::endl;
  pvdFile << "</VTKFile> " << std::endl;

  // Close file
  pvdFile.close();

}
//----------------------------------------------------------------------------
void VTKFile::pvtuFileWrite(bool mesh_function)
{
  std::fstream pvtuFile;


  // Open pvtu file
  pvtuFile.open(pvtu_filename.c_str(), std::ios::out|std::ios::trunc);
  // Write header
  pvtuFile << "<?xml version=\"1.0\"?> " << std::endl;
  pvtuFile << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\">" << std::endl;
  pvtuFile << "<PUnstructuredGrid GhostLevel=\"0\">" << std::endl;
  pvtuFile << "<PCellData>" << std::endl;
  pvtuFile << "<PDataArray  type=\"UInt32\"  Name=\"connectivity\" />" << std::endl;
  pvtuFile << "<PDataArray  type=\"UInt32\"  Name=\"offsets\" />" << std::endl;
  pvtuFile << "<PDataArray  type=\"UInt8\"  Name=\"types\" />"  << std::endl;
  pvtuFile<<"</PCellData>" << std::endl;

  pvtuFile << "<PPoints>" <<std::endl;
  pvtuFile << "<PDataArray  type=\"Float32\"  NumberOfComponents=\"3\" />" << std::endl;
  pvtuFile << "</PPoints>" << std::endl;

  if (mesh_function)
  {
    pvtuFile << "<PCellData Scalars=\"U\">" << std::endl;
    pvtuFile << "<PDataArray  type=\"Float32\"  Name=\"U\" />" << std::endl;
    pvtuFile<<"</PCellData>" << std::endl;
  }

  std::string fname;
  // Remove rank from vtu filename ( <rank>.vtu)
  fname.assign(vtu_filename, filename.find_last_of("/") + 1, vtu_filename.size() - 5 );
  for(uint i=0; i< MPI::numProcesses(); i++)
  {
    pvtuFile << "<Piece Source=\"" << fname << i << ".vtu\"/>" << std::endl;
  }


  pvtuFile << "</PUnstructuredGrid>" << std::endl;
  pvtuFile << "</VTKFile>" << std::endl;
  pvtuFile.close();

}//----------------------------------------------------------------------------
void VTKFile::pvtuFileWrite_func(std::vector<std::pair<Function*, std::string> > f)
{
  std::fstream pvtuFile;


  // Open pvtu file
  pvtuFile.open(pvtu_filename.c_str(), std::ios::out|std::ios::trunc);
  // Write header
  pvtuFile << "<?xml version=\"1.0\"?> " << std::endl;
  pvtuFile << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\">" << std::endl;
  pvtuFile << "<PUnstructuredGrid GhostLevel=\"0\">" << std::endl;

  pvtuFile << "<PPointData>" << std::endl;
  for (std::vector<std::pair<Function*, std::string> >::iterator it = f.begin();
       it != f.end(); it++)
  {
    Function* u = it->first;
    std::string& name = it->second;

    // Check type of function space
    if(u->type() == Function::discrete && u->space().is_cellwise_constant())
    {
      // These are cell based function and will be written as CellData
      continue;
    }

    if(u->rank() == 0)
    {
      pvtuFile << "<PDataArray  type=\"Float32\"  Name=\"" << name << "\" />"
               << std::endl;
    }
    else
    {
      // Get number of components
      pvtuFile << "<PDataArray  type=\"Float32\"  Name=\"" << name
               << "\"  NumberOfComponents=\"3\" />" << std::endl;
    }
  }
  pvtuFile << "</PPointData>" << std::endl;

  pvtuFile << "<PCellData>" << std::endl;
  for (std::vector<std::pair<Function*, std::string> >::iterator it = f.begin();
       it != f.end(); it++)
  {
    Function* u = it->first;
    std::string& name = it->second;

    // Check type of function space
    if(!u->type() == Function::discrete || !u->space().is_cellwise_constant())
    {
      // These are not cell based functions and will be written as PointData
      continue;
    }

    if(u->rank() == 0)
    {
      pvtuFile << "<PDataArray  type=\"Float32\"  Name=\"" << name << "\" />"
               << std::endl;
    }
    else
    {
      // Get number of components
      pvtuFile << "<PDataArray  type=\"Float32\"  Name=\"" << name
               << "\"  NumberOfComponents=\"3\" />" << std::endl;
    }
  }
  pvtuFile << "</PCellData>" << std::endl;


  pvtuFile << "<PCellData>" << std::endl;
  pvtuFile << "<PDataArray  type=\"UInt32\"  Name=\"connectivity\"  />" << std::endl;
  pvtuFile << "<PDataArray  type=\"UInt32\"  Name=\"offsets\" />" << std::endl;
  pvtuFile << "<PDataArray  type=\"UInt8\"  Name=\"types\" />"  << std::endl;
  pvtuFile<<"</PCellData>" << std::endl;

  pvtuFile << "<PPoints>" <<std::endl;
  pvtuFile << "<PDataArray  type=\"Float32\"  NumberOfComponents=\"3\" />" << std::endl;
  pvtuFile << "</PPoints>" << std::endl;

  std::string fname;

  // Remove rank from vtu filename ( <rank>.vtu)
  fname.assign(vtu_filename, filename.find_last_of("/") + 1, vtu_filename.size() - 5 );
  for(uint i=0; i< MPI::numProcesses(); i++)
    pvtuFile << "<Piece Source=\"" << fname << i << ".vtu\"/>" << std::endl;


  pvtuFile << "</PUnstructuredGrid>" << std::endl;
  pvtuFile << "</VTKFile>" << std::endl;
  pvtuFile.close();

}
//----------------------------------------------------------------------------
void VTKFile::VTKHeaderOpen(Mesh& mesh) const
{
  // Open file
  FILE *fp = fopen(vtu_filename.c_str(), "a");

    // Figure out endianness of machine
  std::string endianness = "";
#if defined HAVE_BIG_ENDIAN
  endianness = "byte_order=\"BigEndian\"";;
#else
  endianness = "byte_order=\"LittleEndian\"";
#endif

  std::string compressor = "";
#ifdef HAVE_LIBZ
  compressor = "compressor=\"vtkZLibDataCompressor\"";
#endif

  // Write headers
  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\"  version=\"0.1\" %s %s>\n",
          endianness.c_str(), compressor.c_str());
  fprintf(fp, "<UnstructuredGrid>  \n");
  fprintf(fp, "<Piece  NumberOfPoints=\" %8u\"  NumberOfCells=\" %8u\">  \n",
          mesh.numVertices(), mesh.numCells());

  // Close file
  fclose(fp);
}
//----------------------------------------------------------------------------
void VTKFile::VTKHeaderClose() const
{
  // Open file
  FILE *fp = fopen(vtu_filename.c_str(), "a");

  // Close headers
  fprintf(fp, "</Piece> \n </UnstructuredGrid> \n </VTKFile>");

  // Close file
  fclose(fp);
}
//----------------------------------------------------------------------------
void VTKFile::vtuNameUpdate(const int counter)
{
  std::string filestart, extension;
  std::ostringstream fileid, newfilename;

  fileid.fill('0');
  fileid.width(6);

  filestart.assign(filename, 0, filename.find("."));
  extension.assign(filename, filename.find("."), filename.size());

  fileid << counter;
  newfilename << filestart << fileid.str() << "_" << MPI::processNumber() <<".vtu";
  vtu_filename = newfilename.str();

  // Make sure file is empty
  FILE* fp = fopen(vtu_filename.c_str(), "w");
  fclose(fp);
}
//----------------------------------------------------------------------------
void VTKFile::pvtuNameUpdate(const int counter)
{
  std::string filestart, extension;
  std::ostringstream fileid, newfilename;

  fileid.fill('0');
  fileid.width(6);

  filestart.assign(filename, 0, filename.find("."));
  extension.assign(filename, filename.find("."), filename.size());

  fileid << counter;
  newfilename << filestart << fileid.str() << ".pvtu";

  pvtu_filename = newfilename.str();

  // Make sure file is empty
  FILE* fp = fopen(pvtu_filename.c_str(), "w");
  fclose(fp);
}
//----------------------------------------------------------------------------
template<class T>
void VTKFile::MeshFunctionWrite(T& meshfunction)
{
  // Update vtu file name and clear file
  vtuNameUpdate(counter);

  // Write pvd file
  if(MPI::processNumber() == 0)
  {
    pvtuNameUpdate(counter);
    pvdFileWrite(counter);
    pvtuFileWrite(true);

  }

  Mesh& mesh = meshfunction.mesh();

  if( meshfunction.dim() != mesh.topology().dim() )
  {
    error("VTK output of mesh functions is implemented for cell-based functions only.");
  }

  // Write headers
  VTKHeaderOpen(mesh);

  // Write mesh
  MeshWrite(mesh);

  std::vector<float> data;
  data.resize(mesh.numCells());
  std::vector<float>::iterator entry = data.begin();

  for (CellIterator cell(mesh); !cell.end(); ++cell)
    *entry++ = static_cast<float>(meshfunction.get(cell->index()));

  // Open file
  FILE *fp = fopen(vtu_filename.c_str(), "a");
  fprintf(fp,"<CellData  Scalars=\"U\">\n" );
  fprintf(fp,"<DataArray  type=\"Float32\"  Name=\"U\"  format=\"binary\">\n");

  // Create encoded stream
  std::stringstream base64_stream;
  encode_stream(base64_stream, data);
  fprintf(fp, "%s\n", base64_stream.str().c_str());

  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</CellData>\n");

  // Close file
  fclose(fp);

  // Close headers
  VTKHeaderClose();

  // Increase the number of times we have saved the mesh function
  ++counter;

  message("saved mesh function %d times.", counter);
  message("Saved mesh function ( %s ) to file %s in VTK format.",
          mesh.label().c_str(), filename.c_str());
}
//-----------------------------------------------------------------------------
template<typename T>
void VTKFile::encode_stream(std::stringstream& stream,
                            const std::vector<T>& data) const
{

#ifdef HAVE_LIBZ
  encode_inline_compressed_base64(stream, data);
#else
  warning("zlib must be configured to enable compressed VTK output. "
          "Using uncompressed base64 encoding instead.");
  encode_inline_base64(stream, data);
#endif

}
//----------------------------------------------------------------------------
template<typename T>
void VTKFile::encode_inline_base64(std::stringstream& stream,
                                   const std::vector<T>& data) const
{
  const uint32_t size = data.size()*sizeof(T);
  Encoder::encode_base64(&size, 1, stream);
  Encoder::encode_base64(data, stream);
}
//----------------------------------------------------------------------------
#ifdef HAVE_LIBZ
template<typename T>
void VTKFile::encode_inline_compressed_base64(std::stringstream& stream,
                                              const std::vector<T>& data) const
{
  uint32_t header[4];
  header[0] = 1;
  header[1] = data.size()*sizeof(T);
  header[2] = 0;

  // Compress data
  std::pair<unsigned char*, unsigned long> compressed_data = Encoder::compress_data(data);

  // Length of compressed data
  header[3] = compressed_data.second;

  // Encode header
  Encoder::encode_base64(&header[0], 4, stream);

  // Encode data
  Encoder::encode_base64(compressed_data.first, compressed_data.second, stream);

  // Free compress data
  delete[] compressed_data.first;
}
#endif
//----------------------------------------------------------------------------

}

