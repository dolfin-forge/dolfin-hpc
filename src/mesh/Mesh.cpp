// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Johan Hoffman, 2007.
// Modified by Garth N. Wells 2007.
//
// First added:  2006-05-09
// Last changed: 2008-07-16

#include <sstream>
#include <fstream>

#include <dolfin/io/File.h>
#include <dolfin/mesh/ALE.h>
#include <dolfin/mesh/UniformMeshRefinement.h>
#include <dolfin/mesh/LocalMeshRefinement.h>
#include <dolfin/mesh/LocalMeshCoarsening.h>
#include <dolfin/mesh/TopologyComputation.h>

#include <dolfin/mesh/MeshFunction.h>

#include <dolfin/mesh/MeshSmoothing.h>
#include <dolfin/mesh/MeshOrdering.h>
#include <dolfin/mesh/MeshPartition.h>
#include <dolfin/mesh/BoundaryMesh.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/MPIMeshCommunicator.h>
#include <dolfin/mesh/MeshData.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshRenumber.h>
#include <dolfin/parameter/parameters.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
Mesh::Mesh()
  : Variable("mesh", "DOLFIN mesh"), _data(0), _cell_type(0), _ordered(false)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
Mesh::Mesh(const Mesh& mesh)
  : Variable("mesh", "DOLFIN mesh"), _data(0), _cell_type(0), _ordered(false)
{
  *this = mesh;
}
//-----------------------------------------------------------------------------
Mesh::Mesh(std::string filename)
  : Variable("mesh", "DOLFIN mesh"), _data(0), _cell_type(0), _ordered(false)
{
  File file(filename);
  file >> *this;
  
  if( MPI::numProcesses() > 1 && !dolfin_get("Mesh read in serial")) {
    MeshFunction<uint> partitions;
    partition(partitions);
    distribute(partitions);
    renumber();
  }
}
//-----------------------------------------------------------------------------
/*Mesh::Mesh(std::string filename, std::string geomDataFilename)
  : Variable("mesh", "DOLFIN mesh"), _data(0), _cell_type(0), _ordered(false)
{
  File file(filename);
  file >> *this;
	std::cout<<this->str()<<std::endl;
	//create meshfunctions for the geometry data
	MeshFunction<int>* u_ = 0;
	u_ = this->data().createMeshFunctionInt("u");
	dolfin_assert(u_);
	u_->init(0);
	MeshFunction<real>* v_ = 0;
	v_ = this->data().createMeshFunctionReal("v");
	dolfin_assert(v_);
	v_->init(0);
	MeshFunction<real>* patch_id_ = 0;
	patch_id_ = this->data().createMeshFunctionReal("patch_id");
	dolfin_assert(patch_id_);
	patch_id_->init(0);

	real u_value, v_value;
	int patch_id_value;

	std::ifstream ifile;
  ifile.open(geomDataFilename.c_str());	
	
  int i = 0;
  while(ifile >> patch_id_value >> u_value >> v_value)
  {
	//	if(i<15)
	//		std::cout << patch_id_temp << " " << x << " " << y << std::endl;
		patch_id_->set(i,patch_id_value);
    u_->set(i, u_value);
    v_->set(i, v_value);
    i++;
  }
  
  ifile.close();	

  if( MPI::numProcesses() > 1 && !dolfin_get("Mesh read in serial")) {
    MeshFunction<uint> partitions;
    partition(partitions);
    distribute(partitions);
    renumber();
  }
}*/
//-----------------------------------------------------------------------------
Mesh::~Mesh()
{
  clear();
}
//-----------------------------------------------------------------------------
const Mesh& Mesh::operator=(const Mesh& mesh)
{
  clear();

  _topology = mesh._topology;
  _geometry = mesh._geometry;
  _distdata = mesh._distdata;
  
  if (mesh._cell_type)
    _cell_type = CellType::create(mesh._cell_type->cellType());
  
  rename(mesh.name(), mesh.label());

  return *this;
}
//-----------------------------------------------------------------------------
MeshData& Mesh::data()
{
  if (!_data)
    _data = new MeshData(*this);

  return *_data;
}
//-----------------------------------------------------------------------------
dolfin::uint Mesh::init(uint dim)
{
  return TopologyComputation::computeEntities(*this, dim);
}
//-----------------------------------------------------------------------------
void Mesh::init(uint d0, uint d1)
{
  TopologyComputation::computeConnectivity(*this, d0, d1);
}
//-----------------------------------------------------------------------------
void Mesh::init()
{
  // Compute all entities
  for (uint d = 0; d <= topology().dim(); d++)
    init(d);

  // Compute all connectivity
  for (uint d0 = 0; d0 <= topology().dim(); d0++)
    for (uint d1 = 0; d1 <= topology().dim(); d1++)
      init(d0, d1);
}
//-----------------------------------------------------------------------------
void Mesh::clear()
{
  _topology.clear();
  _geometry.clear();
  delete _cell_type;
  _cell_type = 0;
  delete _data;
  _data = 0;
}
//-----------------------------------------------------------------------------
void Mesh::order()
{
  if (_ordered)
    message(1, "Mesh has already been ordered, no need to reorder entities.");
  else
    MeshOrdering::order(*this);
}
//-----------------------------------------------------------------------------
bool Mesh::ordered() const
{
  return _ordered;
}
//-----------------------------------------------------------------------------
void Mesh::refine()
{
  message("No cells marked for refinement, assuming uniform mesh refinement.");
  UniformMeshRefinement::refine(*this);
}
//-----------------------------------------------------------------------------
void Mesh::refine(libgeom::Geometry& geom, MeshFunction<int>& patch_id_list, MeshFunction<float>& bnd_u, MeshFunction<float>& bnd_v)
{
  message("No cells marked for refinement, assuming uniform mesh refinement with geometry informations.");
  UniformMeshRefinement::refine(*this, geom, patch_id_list, bnd_u, bnd_v);
}
//-----------------------------------------------------------------------------
void Mesh::refine(libgeom::Geometry& geom, MeshFunction<int>& patch_id_list, MeshFunction<float>& bnd_u)
{
  message("No cells marked for refinement, assuming uniform mesh refinement with geometry informations.");
  UniformMeshRefinement::refine(*this, geom, patch_id_list, bnd_u);
}
//-----------------------------------------------------------------------------
void Mesh::refine(MeshFunction<bool>& cell_markers, bool refine_boundary,
		  bool load_balance)
{
  LocalMeshRefinement::refineMeshByEdgeBisection(*this, cell_markers,
                                                 refine_boundary, load_balance);
}
//-----------------------------------------------------------------------------
void Mesh::coarsen()
{
  // FIXME: Move implementation to separate class and just call function here

  message("No cells marked for coarsening, assuming uniform mesh coarsening.");
  MeshFunction<bool> cell_marker(*this);
  cell_marker.init(this->topology().dim());
  for (CellIterator c(*this); !c.end(); ++c)
    cell_marker.set(c->index(),true);

  LocalMeshCoarsening::coarsenMeshByEdgeCollapse(*this,cell_marker);
}
//-----------------------------------------------------------------------------
void Mesh::coarsen(MeshFunction<bool>& cell_markers, bool coarsen_boundary)
{
  LocalMeshCoarsening::coarsenMeshByEdgeCollapse(*this, cell_markers,
                                                 coarsen_boundary);
}
//-----------------------------------------------------------------------------
void Mesh::move(Mesh& boundary, ALEType method)
{
  ALE::move(*this, boundary, method);
}
//-----------------------------------------------------------------------------
void Mesh::smooth() 
{
  MeshSmoothing::smooth(*this);
}
//-----------------------------------------------------------------------------
void Mesh::partition(MeshFunction<uint>& partitions)
{
  MeshPartition::partition(*this, partitions);
}
//-----------------------------------------------------------------------------
void Mesh::partition(MeshFunction<uint>& partitions,MeshFunction<uint>& weight)
{
  MeshPartition::partition(*this, partitions, weight);
}
//-----------------------------------------------------------------------------
void Mesh::partition(MeshFunction<uint>& partitions, uint num_partitions)
{
  // Partition mesh
  MeshPartition::partition(*this, partitions, num_partitions);
}
//-----------------------------------------------------------------------------
void Mesh::partition_geom(MeshFunction<uint>& partitions)
{
  MeshPartition::partition_geom(*this, partitions);
}
//-----------------------------------------------------------------------------
void Mesh::distribute(MeshFunction<uint>& distribution)
{
  MPIMeshCommunicator::distribute(*this, distribution);
}
//-----------------------------------------------------------------------------
void Mesh::distribute(MeshFunction<uint>& distribution, 
		      MeshFunction<bool>& cell_markers,
		      MeshFunction<bool>& new_cell_markers)
{
  MPIMeshCommunicator::distribute(*this, distribution, 
				  cell_markers, new_cell_markers);
}
//-----------------------------------------------------------------------------
void Mesh::renumber()
{
  MeshRenumber::renumber(*this);
}
//-----------------------------------------------------------------------------
void Mesh::disp() const
{
  // Begin indentation
  cout << "Mesh data" << endl;
  begin("---------");
  cout << endl;

  // Display topology and geometry
  _topology.disp();
  _geometry.disp();

  // Display cell type
  cout << "Cell type" << endl;
  cout << "---------" << endl;
  begin("");
  if (_cell_type)
    cout << _cell_type->description() << endl;
  else
    cout << "undefined" << endl;
  end();

  // Display mesh data
  if (_data)
  {
    cout << endl;
    _data->disp();
  }
  
  // End indentation
  end();
}
//-----------------------------------------------------------------------------
std::string Mesh::str() const
{
  std::ostringstream stream;
  stream << "[Mesh of topological dimension "
         << topology().dim()
         << " with "
         << numVertices()
         << " and "
         << numCells()
         << " cells]";
  return stream.str();
}
//-----------------------------------------------------------------------------
dolfin::LogStream& dolfin::operator<< (LogStream& stream, const Mesh& mesh)
{
  stream << mesh.str();
  return stream;
}
//-----------------------------------------------------------------------------
