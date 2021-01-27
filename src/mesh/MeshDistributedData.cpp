// Copyright (C) 2008 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/mesh/MeshDistributedData.h>

#include <dolfin/log/log.h>
#include <dolfin/log/LogStream.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
MeshDistributedData::MeshDistributedData(uint dim) :
    dim_(dim),
    data_(std::vector<DistributedData>(dim_ + 1))
{
}
//-----------------------------------------------------------------------------
MeshDistributedData::MeshDistributedData(MeshDistributedData const& other) :
    dim_(other.dim_),
    data_(std::vector<DistributedData>(dim_ + 1))
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
auto MeshDistributedData::operator=(MeshDistributedData const& other) -> MeshDistributedData&
{
  clear();
  dim_ = other.dim_;
  data_.resize(dim_+1);
  for (uint i = 0; i <= dim_; ++i)
  {
    data_[i] = other.data_[i];
  }
  return *this;
}
//-----------------------------------------------------------------------------
auto MeshDistributedData::operator==(MeshDistributedData const& other) const -> bool
{
  if (dim_ != other.dim_)
  {
    return false;
  }

  for (uint i = 0; i <= dim_; ++i)
  {
    if(data_[i] != other.data_[i])
    {
      return false;
    }
  }

  return true;
}
//-----------------------------------------------------------------------------
auto MeshDistributedData::operator!=(MeshDistributedData const& other) const -> bool
{
  return !(*this == other);
}
//-----------------------------------------------------------------------------
void MeshDistributedData::clear()
{
  dim_ = 0;
}
//-----------------------------------------------------------------------------
auto MeshDistributedData::dim() const -> uint
{
  return dim_;
}
//-----------------------------------------------------------------------------
void MeshDistributedData::disp() const
{
  section("MeshDistributedData");
  cout << "Topological dimension     : " << (uint) dim_ << "\n";
  for (uint d = 0; d <= dim_; ++d)
  {
    skip();
    section("");
    (*this)[d].disp();
    end();
  }
  end();
  skip();
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
