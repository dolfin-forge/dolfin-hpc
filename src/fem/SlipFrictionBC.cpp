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
SlipFrictionBC::SlipFrictionBC(Mesh& mesh, const SubDomain& sub_domain, real beta) :
    BoundaryCondition("SlipFriction"),
    slipbc_(mesh, sub_domain),
    beta_(beta),
    expr_(sub_domain, beta),
    Fbeta_(mesh, expr_)
{
}

//-----------------------------------------------------------------------------
SlipFrictionBC::SlipFrictionBC(BoundaryNormal& normal,
                               const SubDomain& sub_domain,
                               real beta) :
    BoundaryCondition("SlipFriction"),
    slipbc_(normal, sub_domain),
    beta_(beta),
    expr_(sub_domain, beta),
    Fbeta_(normal.mesh(), expr_)
{
}

////-----------------------------------------------------------------------------
//SlipFrictionBC::SlipFrictionBC(MeshFunction<uint>& sub_domains, uint sub_domain,
//                               real beta) :
//    BoundaryCondition("SlipFriction"),
//    slipbc_(sub_domains, sub_domain),
//    beta_(beta),
//    expr_(sub_domain, beta),
//    Fbeta_(sub_domains.mesh(), expr_)
//{
//}

//-----------------------------------------------------------------------------
SlipFrictionBC::SlipFrictionBC(Mesh& mesh, const SubDomain& sub_domain,
                               SubSystem const& sub_system, real beta) :
    BoundaryCondition("SlipFriction"),
    slipbc_(mesh, sub_domain, sub_system),
    beta_(beta),
    expr_(sub_domain, beta),
    Fbeta_(mesh, expr_)
{
}

////-----------------------------------------------------------------------------
//SlipFrictionBC::SlipFrictionBC(MeshFunction<uint>& sub_domains, uint sub_domain,
//                               SubSystem const& sub_system, real beta) :
//    BoundaryCondition("SlipFriction"),
//    slipbc_(sub_domains, sub_domain, sub_system),
//    beta_(beta),
//    expr_(sub_domain, beta),
//    Fbeta_(sub_domains.mesh(), expr_)
//{
//}

//-----------------------------------------------------------------------------
SlipFrictionBC::~SlipFrictionBC()
{
}

//-----------------------------------------------------------------------------
BoundaryNormal& SlipFrictionBC::normal()
{
  return slipbc_.normal();
}

} /* namespace dolfin */
