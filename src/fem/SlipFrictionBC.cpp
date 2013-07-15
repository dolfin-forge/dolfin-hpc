/*
 * SlipFriction.cpp
 *
 *  Created on: Jul 15, 2013
 *      Author: larcher
 */

#include <dolfin/fem/SlipFrictionBC.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
SlipFrictionBC::SlipFrictionBC(Mesh& mesh, SubDomain& sub_domain, real beta) :
    BoundaryCondition("SlipFriction"),
    slipbc_(mesh, sub_domain),
    beta_(beta),
    expr_(sub_domain, beta),
    Fbeta_(mesh, expr_)
{
}

//-----------------------------------------------------------------------------
SlipFrictionBC::SlipFrictionBC(Mesh& mesh, SubDomain& sub_domain,
                               NodeNormal& node_normal, real beta) :
    BoundaryCondition("SlipFriction"),
    slipbc_(mesh, sub_domain, node_normal),
    beta_(beta),
    expr_(sub_domain, beta),
    Fbeta_(mesh, expr_)
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
SlipFrictionBC::SlipFrictionBC(Mesh& mesh, SubDomain& sub_domain,
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

} /* namespace dolfin */
