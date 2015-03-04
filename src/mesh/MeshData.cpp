// Copyright (C) 2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2008-05-19
// Last changed: 2015-03-03

#include <dolfin/mesh/MeshData.h>

#include <dolfin/common/Array.h>
#include <dolfin/mesh/MeshFunction.h>

namespace dolfin
{

typedef std::map<std::string, MeshFunction<dolfin::uint>*>::iterator mf_iterator;
typedef std::map<std::string, MeshFunction<dolfin::uint>*>::const_iterator mf_const_iterator;

typedef std::map<std::string, Array<dolfin::uint>*>::iterator a_iterator;
typedef std::map<std::string, Array<dolfin::uint>*>::const_iterator a_const_iterator;

//-----------------------------------------------------------------------------
MeshData::MeshData(Mesh& mesh) :
    mesh_(mesh)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
MeshData::~MeshData()
{
  clear();
}
//-----------------------------------------------------------------------------
void MeshData::clear()
{
  for (mf_iterator it = meshfunctions_.begin(); it != meshfunctions_.end(); ++it)
  {
    delete it->second;
  }
  meshfunctions_.clear();

  for (a_iterator it = arrays_.begin(); it != arrays_.end(); ++it)
  {
    delete it->second;
  }
  arrays_.clear();
}
//-----------------------------------------------------------------------------
MeshFunction<dolfin::uint>* MeshData::createMeshFunction(std::string name)
{
  // Check if data already exists
  mf_iterator it = meshfunctions_.find(name);
  if (it != meshfunctions_.end())
  {
    warning("Mesh data named \"%s\" already exists.", name.c_str());
    return it->second;
  }

  // Create new data
  MeshFunction<uint>* f = new MeshFunction<uint>(mesh_);
  dolfin_assert(f);

  // Add to map
  meshfunctions_[name] = f;

  return f;
}
//-----------------------------------------------------------------------------
Array<dolfin::uint>* MeshData::createArray(std::string name, uint size)
{
  // Check if data already exists
  a_iterator it = arrays_.find(name);
  if (it != arrays_.end())
  {
    warning("Mesh data named \"%s\" already exists.", name.c_str());
    return it->second;
  }

  // Create new data
  Array<uint>* a = new Array<uint>(size);
  *a = 0;

  // Add to map
  arrays_[name] = a;

  return a;
}
//-----------------------------------------------------------------------------
MeshFunction<dolfin::uint>* MeshData::meshFunction(std::string name)
{
  MeshFunction<dolfin::uint>* ret = NULL;
  // Check if data exists
  mf_iterator it = meshfunctions_.find(name);
  if (it != meshfunctions_.end())
  {
    ret = it->second;
  }
  else
  {
    warning("Mesh data does not contain the requested mesh function.");
  }
  dolfin_assert(ret != NULL);
  return ret;
}
//-----------------------------------------------------------------------------
Array<dolfin::uint>* MeshData::array(std::string name)
{
  Array<dolfin::uint>* ret = NULL;
  // Check if data exists
  a_iterator it = arrays_.find(name);
  if (it != arrays_.end())
  {
    ret = it->second;
  }
  else
  {
    warning("Mesh data does not contain the requested array.");
  }
  dolfin_assert(ret != NULL);
  return ret;
}
//-----------------------------------------------------------------------------
void MeshData::disp() const
{
  // Begin indentation
  cout << "Auxiliary mesh data" << endl;
  begin("-------------------");
  cout << endl;

  for (mf_const_iterator it = meshfunctions_.begin(); it != meshfunctions_.end();
      ++it)
  {
    cout << "MeshFunction<uint> of size " << it->second->size()
        << " on entities of topological dimension " << it->second->dim()
        << ": \"" << it->first << "\"" << endl;
  }

  for (a_const_iterator it = arrays_.begin(); it != arrays_.end(); ++it)
    cout << "Array<uint> of size " << static_cast<uint>(it->second->size())
        << ": \"" << it->first << "\"" << endl;

  // End indentation
  end();
}
//-----------------------------------------------------------------------------

}
