// Copyright (C) 2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-05-08
// Last changed: 2006-11-01

#include <dolfin/log/dolfin_log.h>
#include <dolfin/mesh/MeshConnectivity.h>
#include <dolfin/mesh/MeshTopology.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
MeshTopology::MeshTopology() :
    _dim(0),
    _num_entities(NULL),
    _connectivity(NULL),
    _distdata(*this),
    _ordered(false)

{
  // Do nothing
}
//-----------------------------------------------------------------------------
MeshTopology::MeshTopology(const MeshTopology& topology) :
    _dim(0),
    _num_entities(NULL),
    _connectivity(NULL),
    _distdata(*this),
    _ordered(false)
{
  *this = topology;
}
//-----------------------------------------------------------------------------
MeshTopology::~MeshTopology()
{
  clear();
}
//-----------------------------------------------------------------------------
const MeshTopology& MeshTopology::operator=(const MeshTopology& topology)
{
  // Clear old data if any
  clear();

  // Allocate data
  _dim = topology._dim;
  _num_entities = new uint[_dim + 1];
  _connectivity = new MeshConnectivity*[_dim + 1];
  for (uint d = 0; d <= _dim; d++)
  {
    _connectivity[d] = new MeshConnectivity[_dim + 1];
  }

  // Copy data
  if (_dim > 0)
  {
    for (uint d = 0; d <= _dim; ++d)
    {
      _num_entities[d] = topology._num_entities[d];
    }
    for (uint d0 = 0; d0 <= _dim; ++d0)
    {
      for (uint d1 = 0; d1 <= _dim; ++d1)
      {
        _connectivity[d0][d1] = topology._connectivity[d0][d1];
      }
    }
  }
  _distdata = topology._distdata;

  return *this;
}
//-----------------------------------------------------------------------------
void MeshTopology::clear()
{
  // Clear parallel data structures
  _distdata.clear();

  // Delete number of mesh entities
  if (_num_entities) delete[] _num_entities;
  _num_entities = 0;

  // Delete mesh connectivity
  if (_connectivity)
  {
    for (uint d = 0; d <= _dim; ++d)
    {
      delete[] _connectivity[d];
    }
    delete[] _connectivity;
  }
  _connectivity = 0;

  // Reset dimension
  _dim = 0;
}
//-----------------------------------------------------------------------------
void MeshTopology::init(uint dim)
{
  // Clear old data if any
  clear();

  // Initialize number of mesh entities
  _num_entities = new uint[dim + 1];
  for (uint d = 0; d <= dim; d++)
  {
    _num_entities[d] = 0;
  }

  // Initialize mesh connectivity
  _connectivity = new MeshConnectivity*[dim + 1];
  for (uint d = 0; d <= dim; d++)
  {
    _connectivity[d] = new MeshConnectivity[dim + 1];
  }

  // Save dimension
  _dim = dim;
  _distdata.init(_dim);
}
//-----------------------------------------------------------------------------
void MeshTopology::init(uint dim, uint size)
{
  dolfin_assert(_num_entities);dolfin_assert(dim <= _dim);

  _num_entities[dim] = size;
}
//-----------------------------------------------------------------------------
void MeshTopology::disp() const
{
  // Begin indentation
  cout << "Mesh topology" << endl;
  begin("-------------");
  cout << endl;

  // Check if empty
  if (_dim == 0)
  {
    cout << "empty" << endl << endl;
    end();
    return;
  }

  // Display topological dimension
  cout << "Topological dimension: " << _dim << endl << endl;

  // Display number of entities for each topological dimension
  cout << "Number of entities:" << endl;
  begin("");
  for (uint d = 0; d <= _dim; ++d)
  {
    cout << "dim = " << d << ": " << _num_entities[d] << endl;
  }
  end();
  cout << endl;

  // Display matrix of connectivities
  cout << "Connectivity:" << endl;
  begin("");
  cout << " ";
  for (uint d1 = 0; d1 <= _dim; ++d1)
  {
    cout << " " << d1;
  }
  cout << endl;
  for (uint d0 = 0; d0 <= _dim; ++d0)
  {
    cout << d0;
    for (uint d1 = 0; d1 <= _dim; ++d1)
    {
      if (_connectivity[d0][d1].size() > 0)
      {
        cout << " x";
      }
      else
      {
        cout << " -";
      }
    }
    cout << endl;
  }
  cout << endl;
  end();

  // Display connectivity for each topological dimension
  for (uint d0 = 0; d0 <= _dim; ++d0)
  {
    for (uint d1 = 0; d1 <= _dim; ++d1)
    {
      if (_connectivity[d0][d1].size() == 0)
      {
        continue;
      }
      cout << "Connectivity " << d0 << " -- " << d1 << ":" << endl;
      begin("");
      _connectivity[d0][d1].disp();
      end();
      cout << endl;
    }
  }

  // End indentation
  end();
}
//-----------------------------------------------------------------------------
