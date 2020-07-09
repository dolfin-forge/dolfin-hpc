
#include <dolfin/mesh/MeshQualityFunction.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
MeshQualityFunction::MeshQualityFunction(Mesh& mesh, uint p) :
    Function(mesh),
    p_(p),
    mqual_(mesh)
{
}

//-----------------------------------------------------------------------------
MeshQualityFunction::~MeshQualityFunction()
{
}

//-----------------------------------------------------------------------------

}

