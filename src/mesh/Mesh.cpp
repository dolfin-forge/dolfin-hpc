// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Johan Hoffman, 2007.
// Modified by Garth N. Wells 2007.
// Modified by Balthasar Reuter, 2013.
// Modified by Aurélien Larcher, 2013.
//
// First added:  2006-05-09
// Last changed: 2013-03-22


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
#include <dolfin/mesh/MeshPartition.h>
#include <dolfin/mesh/BoundaryMesh.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/IntersectionDetector.h>
#include <dolfin/mesh/MPIMeshCommunicator.h>
#include <dolfin/mesh/MeshData.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/parameter/parameters.h>


#ifdef HAVE_LIBGEOM
#include <Geometry.h>
#endif


namespace dolfin {

//-----------------------------------------------------------------------------
Mesh::Mesh() :
  Variable("mesh", "DOLFIN mesh"),
  _is_distributed(false),
  _topology(),
  _geometry(),
  _data(0),
  _cell_type(0),
  _exterior_boundary(NULL),
  _interior_boundary(NULL),
  _intersection_detector(NULL),
  _timestamp(time(0))
{
  // Do nothing
}
//-----------------------------------------------------------------------------
Mesh::Mesh(Mesh const& mesh) :
  Variable("mesh", "DOLFIN mesh"),
  _is_distributed(mesh._is_distributed),
  _topology(),
  _geometry(),
  _data(0),
  _cell_type(0),
  _exterior_boundary(NULL),
  _interior_boundary(NULL),
  _intersection_detector(NULL),
  _timestamp(time(0))
{
  *this = mesh;
}
//-----------------------------------------------------------------------------
Mesh::Mesh(std::string filename) :
  Variable("mesh", "DOLFIN mesh"),
  _is_distributed(false),
  _topology(),
  _geometry(),
  _data(0),
  _cell_type(0),
  _exterior_boundary(NULL),
  _interior_boundary(NULL),
  _intersection_detector(NULL),
  _timestamp(time(0))
{
  File file(filename);
  file >> *this;

  if (MPI::numProcesses() > 1 && !dolfin_get("Mesh read in serial"))
  {
      MeshFunction<uint> partitions;
      partition(partitions);
      distribute(partitions);
      renumber();
      _is_distributed = true;
  }
}
//-----------------------------------------------------------------------------
Mesh::~Mesh()
{
  clear();
}
//-----------------------------------------------------------------------------
const Mesh& Mesh::operator=(const Mesh& mesh)
{
  clear();

  _is_distributed = mesh._is_distributed;
  _topology = mesh._topology;
  _geometry = mesh._geometry;
  _timestamp = mesh._timestamp;

  if (mesh._cell_type)
  {
    _cell_type = CellType::create(mesh._cell_type->cellType());
  }

  rename(mesh.name(), mesh.label());

  return *this;
}
//-----------------------------------------------------------------------------
MeshData& Mesh::data()
{
  if (!_data)
  {
    _data = new MeshData(*this);
  }
  return *_data;
}
//-----------------------------------------------------------------------------
BoundaryMesh& Mesh::exterior_boundary()
{
  ///FIXME: Improve hash logic to regenerate bounday at topology change
  if(_exterior_boundary == NULL
      || (_exterior_boundary->mesh_hash() != this->hash()))
  {
    delete _exterior_boundary;
    _exterior_boundary = new BoundaryMesh(*this, BoundaryMesh::exterior);
  }
  return *_exterior_boundary;
}
//-----------------------------------------------------------------------------
BoundaryMesh& Mesh::interior_boundary()
{
  ///FIXME: Improve hash logic to regenerate bounday at topology change
  if(_interior_boundary == NULL
      || (_interior_boundary->mesh_hash() != this->hash()))
  {
    delete _interior_boundary;
    _interior_boundary = new BoundaryMesh(*this, BoundaryMesh::interior);
  }
  return *_interior_boundary;
}
//-----------------------------------------------------------------------------
IntersectionDetector& Mesh::intersector()
{
  ///FIXME: Improve hash logic to regenerate detector at topology change
  if(_intersection_detector == NULL)
  {
    delete _intersection_detector;
    _intersection_detector = new IntersectionDetector(*this);
  }
  return *_intersection_detector;
}
//-----------------------------------------------------------------------------
dolfin::uint Mesh::init(uint dim) const
{
  Mesh* mesh = const_cast<Mesh*>(this);
  return _topology.compute_entities(*mesh, dim);
}
//-----------------------------------------------------------------------------
void Mesh::init(uint d0, uint d1) const
{
  Mesh* mesh = const_cast<Mesh*>(this);
  _topology.compute_connectivity(*mesh, d0, d1);
}
//-----------------------------------------------------------------------------
void Mesh::init() const
{
  // Compute all entities
  for (uint d = 0; d <= topology().dim(); ++d)
  {
    init(d);
  }

  // Compute all connectivity
  for (uint d0 = 0; d0 <= topology().dim(); ++d0)
  {
    for (uint d1 = 0; d1 <= topology().dim(); ++d1)
    {
      init(d0, d1);
    }
  }
}
//-----------------------------------------------------------------------------
void Mesh::clear()
{
  _timestamp = 0;
  _topology.clear();
  _geometry.clear();
  delete _cell_type;
  _cell_type = NULL;
  delete _data;
  _data = NULL;
  delete _exterior_boundary;
  _exterior_boundary = NULL;
  delete _interior_boundary;
  _interior_boundary = NULL;
  delete _intersection_detector;
  _intersection_detector = NULL;
}
//-----------------------------------------------------------------------------
void Mesh::order()
{
  if (_topology.is_ordered())
  {
    message(1, "Mesh has already been ordered, no need to reorder entities.");
  }
  else
  {
    _topology.order(*this);
  }
}
//-----------------------------------------------------------------------------
bool Mesh::ordered() const
{
  return _topology.is_ordered();
}
//-----------------------------------------------------------------------------
void Mesh::refine()
{
  message("No cells marked for refinement, assuming uniform mesh refinement.");
  UniformMeshRefinement::refine(*this);
}
//-----------------------------------------------------------------------------
#ifdef HAVE_LIBGEOM
//-----------------------------------------------------------------------------
void Mesh::refine(libgeom::Geometry& geom,
                  MeshFunction<int>& patch_id_list,
                  MeshFunction<float>& bnd_u,
                  MeshFunction<float>& bnd_v)
{
  message("No cells marked for refinement, "
          "assuming uniform mesh refinement with geometry informations.");
  UniformMeshRefinement::refine(*this, geom, patch_id_list, bnd_u, bnd_v);
}
//-----------------------------------------------------------------------------
void Mesh::refine(libgeom::Geometry& geom,
                  MeshFunction<int>& patch_id_list, MeshFunction<float>& bnd_u)
{
  message("No cells marked for refinement, "
          "assuming uniform mesh refinement with geometry informations.");
  UniformMeshRefinement::refine(*this, geom, patch_id_list, bnd_u);
}
//-----------------------------------------------------------------------------
#endif // HAVE_LIBGEOM
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
  MeshFunction<bool> cell_marker(*this, this->topology().dim());
  cell_marker = true;

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
void Mesh::partition(MeshFunction<uint>& partitions, MeshFunction<uint>& weight)
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
  MPIMeshCommunicator::distribute(*this, distribution, cell_markers,
                                  new_cell_markers);
}
//-----------------------------------------------------------------------------
void Mesh::distribute(MeshFunction<uint>& distribution,
                      Array<std::pair<MeshFunction<uint> *,
                      MeshFunction<uint> *> >& cell_functions)
{
  MPIMeshCommunicator::distribute(*this, distribution, cell_functions);
}
//-----------------------------------------------------------------------------
void Mesh::distribute(MeshFunction<uint>& distribution,
                      Array<std::pair<MeshFunction<double> *,
                      MeshFunction<double> *> >& vertex_functions)
{
  MPIMeshCommunicator::distribute(*this, distribution, vertex_functions);
}
//-----------------------------------------------------------------------------
void Mesh::distribute(MeshFunction<uint>& distribution,
                      Array<std::pair<MeshFunction<uint> *,
                      MeshFunction<uint> *> >& cell_functions,
                      Array<std::pair<MeshFunction<double> *,
                      MeshFunction<double> *> >& vertex_functions)
{
  MPIMeshCommunicator::distribute(*this, distribution, cell_functions,
                                  vertex_functions);
}
//-----------------------------------------------------------------------------
void Mesh::renumber()
{
  MeshRenumber::renumber(*this);
//  _timestamp = time(0);
}

//-----------------------------------------------------------------------------
std::string const Mesh::hash() const
{
  std::stringstream ss;
  ss << "Mesh@" << this << ":" << this->type().description()
      << ":time" << _timestamp
      << ":T" << this->topology().token()
      << ":G" << this->geometry().token();
  return ss.str();
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
  stream << "[Mesh of topological dimension " << topology().dim() << " with "
      << numVertices() << " and " << numCells() << " cells]";
  return stream.str();
}
//-----------------------------------------------------------------------------
LogStream& operator<<(LogStream& stream, const Mesh& mesh)
{
  stream << mesh.str();
  return stream;
}
//-----------------------------------------------------------------------------

}

