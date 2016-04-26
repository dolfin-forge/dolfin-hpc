// Copyright (C) 2013 Aurélien Larcher
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2013-07-15 (merged from branch larcher)
// Last changed: 2013-07-15

#include <dolfin/fem/SlipFrictionBC.h>
#include <dolfin/fem/BoundaryNormal.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
SlipFrictionBC::SlipFrictionBC(Coefficient& beta, Mesh& mesh,
                               SubDomain const& sub_domain) :
    BoundaryCondition("SlipFriction", mesh, sub_domain),
    slipbc_(mesh, sub_domain),
    beta_(beta)
{
}

////-----------------------------------------------------------------------------
SlipFrictionBC::SlipFrictionBC(Coefficient& beta,
                               MeshFunction<uint>& sub_domains, uint sub_domain) :
    BoundaryCondition("SlipFriction", sub_domains, sub_domain),
    slipbc_(sub_domains, sub_domain),
    beta_(beta)
{
}

//-----------------------------------------------------------------------------
SlipFrictionBC::SlipFrictionBC(Coefficient& beta, Mesh& mesh,
                               SubDomain const& sub_domain,
                               SubSystem const& sub_system) :
    BoundaryCondition("SlipFriction", mesh, sub_domain, sub_system),
    slipbc_(mesh, sub_domain, sub_system),
    beta_(beta)
{
}

////-----------------------------------------------------------------------------
SlipFrictionBC::SlipFrictionBC(Coefficient& beta,
                               MeshFunction<uint>& sub_domains, uint sub_domain,
                               SubSystem const& sub_system) :
    BoundaryCondition("SlipFriction", sub_domains, sub_domain, sub_system),
    slipbc_(sub_domains, sub_domain, sub_system),
    beta_(beta)
{
}

//-----------------------------------------------------------------------------
SlipFrictionBC::~SlipFrictionBC()
{
}

//-----------------------------------------------------------------------------
BoundaryNormal& SlipFrictionBC::normal()
{
  return slipbc_.normal();
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */
