/******************************************************************************
 * Copyright 2012, 2013 Aurélien Larcher
 *
 * Licensed under the EUPL, Version 1.1 only (the "Licence");
 * You may not use this work except in compliance with the Licence.
 * You may obtain a copy of the Licence at:
 *
 * http://ec.europa.eu/idabc/eupl5
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the Licence is distributed on an "AS IS" basis,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the Licence for the specific language governing permissions and
 * limitations under the Licence.
 ******************************************************************************/

#include <dolfin/function/Field.h>

#include <dolfin/log/log.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
Field::Field(Mesh& mesh) :
    storage_(1, Function(mesh)),
    bcs_(),
    Udt_(mesh)
{
}
//-----------------------------------------------------------------------------
Field::Field(FiniteElementSpace const& space) :
    storage_(1, Function(space)),
    bcs_(),
    Udt_(space)
{
}
//-----------------------------------------------------------------------------
Field::~Field()
{
}
//-----------------------------------------------------------------------------

Function& Field::dt()
{
  return Udt_;
}
//-----------------------------------------------------------------------------
uint Field::depth() const
{
  return storage_.size() - 1;
}
//-----------------------------------------------------------------------------
void Field::allocate(uint max_depth)
{
  for (uint i = this->depth(); i < max_depth; ++i)
  {
    // Use copy constructor
    storage_.push_back(Function(storage_.at(0)));
  }
}
//-----------------------------------------------------------------------------
void Field::init(FiniteElementSpace const& space)
{
  for (uint i = 0; i < storage_.size(); ++i)
  {
    if (storage_[i].empty())
    {
      storage_[i].init(space);
    }
  }
}
//-----------------------------------------------------------------------------
void Field::shift()
{
  for (uint i = this->depth(); i > 0; --i)
  {
    storage_[i].swap(storage_[i - 1]);
  }
}
//-----------------------------------------------------------------------------
Field& Field::operator<<(uint n)
{
  n = n % (this->depth() + 1);
  for (uint i = 0; i < n; ++i)
  {
    this->shift();
  }
  return *this;
}
//-----------------------------------------------------------------------------
Function& Field::operator[](uint i)
{
  if (i > this->depth())
  {
    error("Field: requested level '%d', maximum = '%d'", i, this->depth());
  }
  return storage_[i];
}
//-----------------------------------------------------------------------------
Function const& Field::operator[](uint i) const
{
  if (i > this->depth())
  {
    error("Field: requested level '%d', maximum = '%d'", i, this->depth());
  }
  return storage_[i];
}
//-----------------------------------------------------------------------------
SetOfBCs& Field::bcs()
{
  return bcs_;
}
//-----------------------------------------------------------------------------
SetOfBCs const& Field::bcs() const
{
  return bcs_;
}
//-----------------------------------------------------------------------------
void Field::disp() const
{
  section("Field");
  message("Storage size        : %d", storage_.size());
  message("Boundary conditions : %d", bcs_.size());
  begin("Functions           :");
  for (uint i = 0; i < storage_.size(); ++i)
  {
    if (storage_[i].empty())
    {
      message("%2d. Empty", i);
    }
    else
    {
      message("%2d. t = %e; min = %+e, max = %+e", storage_[i].time(), i,
              storage_[i].min(), storage_[i].max());
    }
  }
  end();
  end();
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
