// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-28
// Last changed: 2014-01-28

#include <dolfin/ufl/UFLFamily.h>
#include <dolfin/ufl/UFLElementList.h>

namespace ufl
{

//-----------------------------------------------------------------------------
Family::Family(Family::Type const& t) :
    ufl::type<std::string>(ElementList::Supported().name(t),
                           "'" + ElementList::Supported().name(t) + "'"),
    type_(t)
{
}

//-----------------------------------------------------------------------------
Family::~Family()
{
}

} /* namespace icorne */
