// Copyright (C) 2008 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Aurélien Larcher, 2014-2016.
//
// First added:  2008-07-03
// Last changed: 2016-04-04

#include <dolfin/mesh/MeshDistributedData.h>

#include <dolfin/log/log.h>
#include <dolfin/log/LogStream.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
MeshDistributedData::MeshDistributedData(uint dim) :
    dim_(dim),
    data_(new DistributedData[dim_ + 1])
{
}
//-----------------------------------------------------------------------------
MeshDistributedData::MeshDistributedData(MeshDistributedData const& other) :
    dim_(other.dim_),
    data_(new DistributedData[dim_ + 1])
{
  for (uint i = 0; i <= dim_; ++i)
  {
    data_[i] = other.data_[i];
  }
}
//-----------------------------------------------------------------------------
MeshDistributedData::~MeshDistributedData()
{
  clear();
}
//-----------------------------------------------------------------------------
MeshDistributedData& MeshDistributedData::operator=(MeshDistributedData const& other)
{
  clear();
  dim_ = other.dim_;
  data_ = new DistributedData[dim_+1];
  for (uint i = 0; i <= dim_; ++i)
  {
    data_[i] = other.data_[i];
  }
  return *this;
}
//-----------------------------------------------------------------------------
bool MeshDistributedData::operator==(MeshDistributedData const& other) const
{
  if(this == &other)
  {
    return true;
  }
  //
  if (dim_ != other.dim_)
  {
    return false;
  }
  //
  for (uint i = 0; i <= dim_; ++i)
  {
    if(data_[i] != other.data_[i])
    {
      return false;
    }
  }
  //
  return true;
}
//-----------------------------------------------------------------------------
bool MeshDistributedData::operator!=(MeshDistributedData const& other) const
{
  return !(*this == other);
}
//-----------------------------------------------------------------------------
void MeshDistributedData::clear()
{
  delete[] data_;
  data_ = NULL;
  dim_ = 0;
}
//-----------------------------------------------------------------------------
uint MeshDistributedData::dim() const
{
  return dim_;
}
//-----------------------------------------------------------------------------
void MeshDistributedData::disp() const
{
  section("MeshDistributedData");
  cout << "Topological dimension     : " << (uint) dim_ << endl;
  for (uint d = 0; d <= dim_; ++d)
  {
    skip();
    section("");
    (*this)[d].disp();
    endblock();
  }
  end();
  skip();
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
