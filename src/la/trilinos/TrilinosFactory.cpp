// Copyright (C) 2005-2006 Ola Skavhaug.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_TRILINOS

#include <dolfin/la/trilinos/TrilinosFactory.h>

namespace dolfin
{

namespace trilinos
{

//-----------------------------------------------------------------------------

auto Factory::createMatrix() const -> trilinos::Matrix *
{
  return new trilinos::Matrix();
}

//-----------------------------------------------------------------------------

auto Factory::createVector() const -> trilinos::Vector *
{
  return new trilinos::Vector();
}

//-----------------------------------------------------------------------------

auto Factory::createPattern() const -> SparsityPattern *
{
  return new SparsityPattern();
}

//-----------------------------------------------------------------------------

// Singleton instance
Factory Factory::factory;

//-----------------------------------------------------------------------------

} // namespace trilinos

} // namespace dolfin

#endif // HAVE_TRILINOS
