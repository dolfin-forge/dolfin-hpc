// Copyright (C) 2005-2007 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Anders Logg 2005-2006.
// Modified by Kristian Oelgaard 2006.
// Modified by Niclas Jansson 2008-2009.
//
// First added:  2005-07-05
// Last changed: 2009-11-01

#include <boost/cstdint.hpp>
#include <boost/detail/endian.hpp>

#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/function/Function.h>
#include <dolfin/la/Vector.h>
#include "Encoder.h"
#include "PVTKFile.h"


using namespace dolfin;

//----------------------------------------------------------------------------
PVTKFile::PVTKFile(const std::string filename) : GenericFile(filename), _t(0)
{
  type = "VTK";
}
//----------------------------------------------------------------------------
PVTKFile::PVTKFile(const std::string filename, real& t) : 
  GenericFile(filename), _t(&t)
{
  type = "VTK";
}
//----------------------------------------------------------------------------
PVTKFile::~PVTKFile()
{
  // Do nothing
}
//----------------------------------------------------------------------------
void PVTKFile::operator<<(Mesh& mesh)
{
  // Update vtu file name and clear file
  vtuNameUpdate(counter);

  // Only the root updates the pvd file
  if(MPI::processNumber() == 0) {    
    // Update pvtu file name and clear file

    pvtuNameUpdate(counter);
    // Write pvd file
    pvdFileWrite(counter);
    
    // Write pvtu file
    pvtuFileWrite();
  }

  // Write headers
  VTKHeaderOpen(mesh);

  // Write mesh
  MeshWrite(mesh);
  
  // Close headers
  VTKHeaderClose();

  // Increase the number of times we have saved the mesh
  counter++;

  message(1, "Saved mesh %s (%s) to file %s in VTK format.",
          mesh.name().c_str(), mesh.label().c_str(), filename.c_str());
}
//----------------------------------------------------------------------------
void PVTKFile::operator<<(MeshFunction<int>& meshfunction)
{
  MeshFunctionWrite(meshfunction);
}
//----------------------------------------------------------------------------
void PVTKFile::operator<<(MeshFunction<unsigned int>& meshfunction)
{
  MeshFunctionWrite(meshfunction);
}
//----------------------------------------------------------------------------
void PVTKFile::operator<<(MeshFunction<double>& meshfunction)
{
  MeshFunctionWrite(meshfunction);
}
//----------------------------------------------------------------------------
void PVTKFile::operator<<(Function& u)
{
  std::pair<Function*, std::string> f(&u, "U");
  std::vector<std::pair<Function*, std::string> > tmp;
  tmp.push_back(f);
  write_dataset(tmp);
}
//----------------------------------------------------------------------------
void PVTKFile::operator<<(std::vector<std::pair<Function*, std::string> >& f) 
{
  write_dataset(f);
}
//----------------------------------------------------------------------------
void PVTKFile::write_dataset(std::vector<std::pair<Function*, std::string> >& f)
{ 
  dolfin_assert(f.size() > 0);
    
  // Update vtu file name and clear file
  vtuNameUpdate(counter);

  // Write pvd file

  // Only the root updates the pvd file
  if(MPI::processNumber() == 0) {
    
    // Update pvtu file name and clear file
    pvtuNameUpdate(counter);

    // Write pvd file
    pvdFileWrite(counter);

    // Write pvtu file
    pvtuFileWrite_func(f);
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
  counter++;
  /*
  cout << "Saved function " << u.name() << " (" << u.label()
       << ") to file " << filename << " in VTK format." << endl;
  */
}
//----------------------------------------------------------------------------
void PVTKFile::MeshWrite(Mesh& mesh) const
{
  // Open file
  FILE* fp = fopen(vtu_filename.c_str(), "a");

  // Write vertex positions
  fprintf(fp, "<Points>  \n");
  fprintf(fp, "<DataArray  type=\"Float32\"  NumberOfComponents=\"3\"  format=\"binary\"> \n");
  std::vector<float> data;
  data.resize(3*mesh.numVertices());
  std::vector<float>::iterator entry = data.begin();
  for (VertexIterator v(mesh); !v.end(); ++v)
  {
    Point p = v->point();
    *entry++ = p.x();     
    *entry++ = p.y();     
    *entry++ = p.z();        
  }
  
  // Create encoded stream
  std::stringstream base64_stream;
  encode_stream(base64_stream, data);
  fprintf(fp, "%s\n", base64_stream.str().c_str());
  fprintf(fp, "</DataArray>  \n");
  fprintf(fp, "</Points>  \n");
  
  // Write cell connectivity
  fprintf(fp, "<Cells>  \n");
  fprintf(fp, "<DataArray  type=\"Int32\"  Name=\"connectivity\"  format=\"binary\"> \n");
  std::vector<boost::uint32_t> c_data;
  c_data.resize(mesh.numCells() * mesh.type().numEntities(0));
  std::vector<boost::uint32_t>::iterator c_entry = c_data.begin(); 
  for (CellIterator c(mesh); !c.end(); ++c)
  {
    for (VertexIterator v(*c); !v.end(); ++v)
      *c_entry++ = v->index();
  }  
  
  // Create encoded stream
  std::stringstream base64_c_stream;
  encode_stream(base64_c_stream, c_data);
  fprintf(fp, "%s\n", base64_c_stream.str().c_str());
  fprintf(fp, "</DataArray> \n");

  // Write offset into connectivity array for the end of each cell
  fprintf(fp, "<DataArray  type=\"Int32\"  Name=\"offsets\"  format=\"binary\">  \n");
  std::vector<boost::uint32_t>::iterator cc_entry = c_data.begin(); 
  for (uint offsets = 1; offsets <= mesh.numCells(); offsets++)
    *cc_entry++ = offsets*mesh.type().numEntities(0);

  std::stringstream base64_cc_stream;
  encode_stream(base64_cc_stream, c_data);
  fprintf(fp, "%s\n",  base64_cc_stream.str().c_str());
  fprintf(fp, "</DataArray> \n");
  
  //Write cell type
  fprintf(fp, "<DataArray  type=\"UInt8\"  Name=\"types\"  format=\"binary\">  \n");
  std::vector<boost::uint8_t> t_data;
  t_data.resize(mesh.numCells());
  std::vector<boost::uint8_t>::iterator t_entry = t_data.begin();
  for (uint types = 1; types <= mesh.numCells(); types++)
  {
    if (mesh.type().cellType() == CellType::tetrahedron )
      *t_entry++ = boost::uint8_t(10);
    if (mesh.type().cellType() == CellType::triangle )
      *t_entry++ = boost::uint8_t(5);
    if (mesh.type().cellType() == CellType::interval )
      *t_entry++ = boost::uint8_t(3);
  }

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
void PVTKFile::ResultsWrite(std::vector<std::pair<Function*, std::string> > f) const
{
  // Open file
  FILE *fp = fopen(vtu_filename.c_str(), "a");
  
  // Assume same mesh for all data arrays
  Mesh& mesh = f[0].first->mesh();

  fprintf(fp, "<PointData> \n");  
  for (std::vector<std::pair<Function*, std::string> >::iterator it = f.begin();
       it != f.end(); it++) 
  {
    Function* u = it->first;
    std::string& name = it->second;
    const uint rank = u->rank();
    if(rank > 1)
      error("Only scalar and vectors functions can be saved in VTK format.");
    
    // Get number of components
    const uint dim = u->dim(0);
    
  // Allocate memory for function values at vertices
    uint size = mesh.numVertices();
    for (uint i = 0; i < u->rank(); i++)
      size *= u->dim(i);
    real* values = new real[size];
    
    // Get function values at vertices
    u->interpolate(values);
    
    // Write function data at mesh vertices
    if ( rank == 0 )
    {
      //      fprintf(fp, "<PointData  Scalars=\"%s\"> \n", name.c_str());
      fprintf(fp, "<DataArray  type=\"Float32\"  Name=\"%s\"  format=\"binary\">	 ", name.c_str());
    }
    else
    {
      //      fprintf(fp, "<PointData  Vectors=\"%s\"> \n", name.c_str());
      fprintf(fp, "<DataArray  type=\"Float32\"  Name=\"%s\"  NumberOfComponents=\"3\" format=\"binary\">	 ",
	      name.c_str());	
    }
    
    if ( dim > 3 )
      warning("Cannot handle VTK file with number of components > 3. Writing first three components only");
    
    std::vector<float> data;
    data.resize(size);
    std::vector<float>::iterator entry = data.begin();
    
    for (VertexIterator vertex(mesh); !vertex.end(); ++vertex)
    {    
      if ( rank == 0 )
	*entry++ = values[ vertex->index() ] ;
      else if ( u->dim(0) == 2 ) 
      {
	*entry++ = values[ vertex->index() ];
	*entry++ = values[ vertex->index() + mesh.numVertices()];
	*entry++ = 0.0;
      }
      else
      {
	*entry++ = values[ vertex->index() ];
	*entry++ = values[ vertex->index() + mesh.numVertices()];
	*entry++ = values[ vertex->index() + 2*mesh.numVertices()];
      }    
    }	 
    
    // Create encoded stream
    std::stringstream base64_stream;
    encode_stream(base64_stream, data);
    fprintf(fp, "%s\n", base64_stream.str().c_str());
    
    fprintf(fp, "</DataArray> \n");

    delete [] values;
  }
  fprintf(fp, "</PointData> \n");
  
  // Close file
  fclose(fp);
 
}
//----------------------------------------------------------------------------
void PVTKFile::pvdFileWrite(uint num)
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
  fname.assign(pvtu_filename, filename.find_last_of("/") + 1, pvtu_filename.size()); 
  
  // Data file name 
  if(_t)
    pvdFile << "<DataSet timestep=\"" << *_t << "\" part=\"0\"" << " file=\"" <<  fname <<  "\"/>" << std::endl; 
  else
    pvdFile << "<DataSet timestep=\"" << num << "\" part=\"0\"" << " file=\"" <<  fname <<  "\"/>" << std::endl; 
  mark = pvdFile.tellp();
  
  // Close headers
  pvdFile << "</Collection> " << std::endl;
  pvdFile << "</VTKFile> " << std::endl;
  
  // Close file
  pvdFile.close();  

}
//----------------------------------------------------------------------------
void PVTKFile::pvtuFileWrite(bool mesh_function)
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
    pvtuFile << "<Piece Source=\"" << fname << i << ".vtu\"/>" << std::endl; 
  
  
  pvtuFile << "</PUnstructuredGrid>" << std::endl;
  pvtuFile << "</VTKFile>" << std::endl;
  pvtuFile.close();
    
}//----------------------------------------------------------------------------
void PVTKFile::pvtuFileWrite_func(std::vector<std::pair<Function*, std::string> > f)
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
    
    if(u->rank() == 0) {
      //      pvtuFile << "<PPointData Scalars=\"" << name << "\">" << std::endl;    
      pvtuFile << "<PDataArray  type=\"Float32\"  Name=\"" << name << "\" />" << std::endl;
    }
    else {
      //      pvtuFile << "<PPointData Vectors=\"" << name << "\">" << std::endl;    
      pvtuFile << "<PDataArray  type=\"Float32\"  Name=\"" << name << "\"  NumberOfComponents=\"3\" />" << std::endl;
    }
  }    
  pvtuFile << "</PPointData>" << std::endl;
    

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
void PVTKFile::VTKHeaderOpen(Mesh& mesh) const
{
  // Open file
  FILE *fp = fopen(vtu_filename.c_str(), "a");
  
    // Figure out endianness of machine
  std::string endianness = "";
#if defined BOOST_LITTLE_ENDIAN
  endianness = "byte_order=\"LittleEndian\"";
#elif defined BOOST_BIG_ENDIAN
  endianness = "byte_order=\"BigEndian\"";;
#else
  error("Unable to determine the endianness of the machine for VTK binary output.");
#endif
  
  std::string compressor = "";
#ifdef HAS_ZLIB
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
void PVTKFile::VTKHeaderClose() const
{
  // Open file
  FILE *fp = fopen(vtu_filename.c_str(), "a");
  
  // Close headers
  fprintf(fp, "</Piece> \n </UnstructuredGrid> \n </VTKFile>"); 	
  
  // Close file
  fclose(fp);
}
//----------------------------------------------------------------------------
void PVTKFile::vtuNameUpdate(const int counter) 
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
void PVTKFile::pvtuNameUpdate(const int counter)
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
void PVTKFile::MeshFunctionWrite(T& meshfunction) 
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
    error("VTK output of mesh functions is implemenetd for cell-based functions only.");    

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
  counter++;

  cout << "saved mesh function " << counter << " times." << endl;

  cout << "Saved mesh function " << mesh.name() << " (" << mesh.label()
       << ") to file " << filename << " in VTK format." << endl;
}    
//-----------------------------------------------------------------------------
template<typename T>
void PVTKFile::encode_stream(std::stringstream& stream, 
                            const std::vector<T>& data) const
{

#ifdef HAS_ZLIB
  encode_inline_compressed_base64(stream, data);
#else
  warning("zlib must be configured to enable compressed VTK output. Using uncompressed base64 encoding instead.");
  encode_inline_base64(stream, data);
#endif
  
}
//----------------------------------------------------------------------------
template<typename T>
void PVTKFile::encode_inline_base64(std::stringstream& stream, 
                                   const std::vector<T>& data) const
{
  const boost::uint32_t size = data.size()*sizeof(T);
  Encoder::encode_base64(&size, 1, stream);
  Encoder::encode_base64(data, stream);
}
//----------------------------------------------------------------------------
#ifdef HAS_ZLIB
template<typename T>
void PVTKFile::encode_inline_compressed_base64(std::stringstream& stream, 
                                              const std::vector<T>& data) const
{
  boost::uint32_t header[4];
  header[0] = 1;
  header[1] = data.size()*sizeof(T);
  header[2] = 0;

  // Compress data
  std::pair<boost::shared_array<unsigned char>, unsigned long> compressed_data = Encoder::compress_data(data);

  // Length of compressed data
  header[3] = compressed_data.second;

  // Encode header
  Encoder::encode_base64(&header[0], 4, stream);

  // Encode data
  Encoder::encode_base64(compressed_data.first.get(), compressed_data.second, stream);
}
#endif
//----------------------------------------------------------------------------
