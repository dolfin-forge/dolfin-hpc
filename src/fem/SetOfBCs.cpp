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

#include <dolfin/fem/SetOfBCs.h>

#include <dolfin/log/log.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
SetOfBCs::SetOfBCs()
{
}

//-----------------------------------------------------------------------------
SetOfBCs::~SetOfBCs()
{
}

//-----------------------------------------------------------------------------
void SetOfBCs::add(BoundaryCondition& bc)
{
  this->push_back(&bc);
}

//-----------------------------------------------------------------------------
void SetOfBCs::apply(GenericMatrix& AA, GenericVector& bb,
                     BilinearForm& form) const
{
  for (SetOfBCs::const_iterator it = this->begin(); it != this->end(); ++it)
  {
    message("Apply %s boundary condition", (*it)->type().c_str());
    (*it)->apply(AA, bb, form);
  }
}

//-----------------------------------------------------------------------------
void SetOfBCs::apply(GenericMatrix& AA, GenericVector& bb,
                     BilinearForm& form, SubSystem const sub) const
{
  for (SetOfBCs::const_iterator it = this->begin(); it != this->end(); ++it)
  {
    message("Apply %s boundary condition to subsystem %s",
            (*it)->type().c_str(), sub.str().c_str());
    (*it)->apply(AA, bb, form, sub);
  }
}

//-----------------------------------------------------------------------------
bool SetOfBCs::has(std::string const& type) const
{
  bool ret = false;

  for (SetOfBCs::const_iterator it = this->begin(); it != this->end(); ++it)
  {
    if ((*it)->type() == type)
    {
      ret = true;
      break;
    }
  }

  return ret;
}

//-----------------------------------------------------------------------------
Array<BoundaryCondition *> SetOfBCs::get(std::string const& type) const
{
  Array<BoundaryCondition *> ret;
  SetOfBCs::const_iterator it = this->begin();
  for (; it != this->end(); ++it)
  {
    if ((*it)->type() == type)
    {
      ret.push_back(*it);
      break;
    }
  }
  return ret;
}

//-----------------------------------------------------------------------------
uint SetOfBCs::size() const
{
  return Array<BoundaryCondition *>::size();
}

//-----------------------------------------------------------------------------
void SetOfBCs::disp() const
{
  section("SetOfBcs");
  if (this->size() > 0)
  {
    uint ii = 0;
    for (SetOfBCs::const_iterator it = this->begin(); it != this->end(); ++it)
    {
      message("%2d: type = %s (@%p)", ++ii, (*it)->type().c_str(), *it);
    }
  }
  else
  {
    message("Empty");
  }
  end();
}

//-----------------------------------------------------------------------------
void SetOfBCs::clear()
{
  while (!this->empty())
  {
    delete this->back();
    this->pop_back();
  }
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */
