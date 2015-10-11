// Copyright (C) 2015 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//

#include <dolfin/mesh/StructuredGrid.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
StructuredGrid::StructuredGrid(CellType const& type, uint N) :
    Mesh(),
    bbox_(type.dim()),
    n_(N)
{
  init(type);
}

//-----------------------------------------------------------------------------
StructuredGrid::~StructuredGrid()
{
}

//-----------------------------------------------------------------------------
void StructuredGrid::init(CellType const& type)
{
  if (n_ == 0)
  {
    error("StructuredGrid::init : provided number of cells per axis is zero");
  }

  type.disp();

  //
  uint n[EuclideanSpace::MAX_DIMENSION];
  uint m[EuclideanSpace::MAX_DIMENSION] = { 0 };
  n[0] = (type.dim() > 0 ? n_ : 0);
  n[1] = (type.dim() > 1 ? n_ : 0);
  n[2] = (type.dim() > 2 ? n_ : 0);
  uint const dim = type.dim();
  MeshEditor editor(*this, type, dim);
  uint num_verts = 1;
  uint num_bricks = 1;
  for (uint i = 0; i < type.dim(); ++i)
  {
    num_verts *= (n[i] + 1);
    num_bricks *= n[i];
    m[i] = num_verts;
  }

  //
  uint i[EuclideanSpace::MAX_DIMENSION] = { 0 };
  real h[EuclideanSpace::MAX_DIMENSION] = { 0.0 };
  real x[EuclideanSpace::MAX_DIMENSION] = { 0.0 };

  // Create vertices
  uint vertex = 0;
  for (uint d = 0; d < type.dim(); ++d)
  {
    i[d] = 0;
    h[d] = 1.0 / static_cast<real>(n[d]);
    x[d] = bbox_[0][d];
  }
  editor.initVertices(num_verts);
  while (vertex < num_verts)
  {
    editor.addVertex(vertex++, x);
    x[0] = bbox_[0][0] + (vertex % m[0]) * h[0];
    for (uint d = 1; d < dim; ++d)
    {
      if (vertex % m[d - 1] == 0)
      {
        i[d] = (i[d] + 1) % (n[d] + 1);
        x[d] = bbox_[0][d] + i[d] * h[d];
      }
    }
  }

  // Create cells
  uint cell = 0;
  uint const v0 = 0;
  uint const v1 = 1;
  uint const v2 = (n[0] + 1);
  uint const v3 = (n[0] + 1) + 1;
  uint const v4 = v0 + (n[0] + 1) * (n[1] + 1);
  uint const v5 = v1 + (n[0] + 1) * (n[1] + 1);
  uint const v6 = v2 + (n[0] + 1) * (n[1] + 1);
  uint const v7 = v3 + (n[0] + 1) * (n[1] + 1);
  switch (type.cellType())
    {
    case CellType::interval:
      {
        uint cv[2] = { v0, v1 };
        editor.initCells(num_bricks);
        for (uint i = 0; i < n[0]; ++i)
        {
          editor.addCell(cell++, &cv[0]);
          ++cv[0];
          ++cv[1];
        }
      }
      break;
    case CellType::quadrilateral:
      {
        uint const co[4] = { v0, v1, v3, v2 };
        uint cv[4] = { 0 };
        editor.initCells(num_bricks);
        for (uint j = 0; j < n[1]; ++j)
        {
          for (uint i = 0; i < n[0]; ++i)
          {
            cv[0] = j * (n[0] + 1) + i;
            cv[1] = cv[0] + co[1];
            cv[2] = cv[0] + co[2];
            cv[3] = cv[0] + co[3];
            editor.addCell(cell++, &cv[0]);
          }
        }
      }
      break;
    case CellType::hexahedron:
      {
        uint const co[8] = { v0, v1, v3, v2, v4, v5, v7, v6 };
        uint cv[8] = { 0 };
        editor.initCells(num_bricks);
        for (uint k = 0; k < n[2]; ++k)
        {
          for (uint j = 0; j < n[1]; ++j)
          {
            for (uint i = 0; i < n[0]; ++i)
            {
              cv[0] = k * (n[1] + 1) * (n[0] + 1) + j * (n[0] + 1) + i;
              cv[1] = cv[0] + co[1];
              cv[2] = cv[0] + co[2];
              cv[3] = cv[0] + co[3];
              cv[4] = cv[0] + co[4];
              cv[5] = cv[0] + co[5];
              cv[6] = cv[0] + co[6];
              cv[7] = cv[0] + co[7];
              editor.addCell(cell++, &cv[0]);
            }
          }
        }
      }
      break;
    case CellType::triangle:
      {
        uint const co[6] = { v0, v1, v3, v0, v2, v3 };
        uint cv[6] = { 0 };
        editor.initCells(2 * num_bricks);
        for (uint j = 0; j < n[1]; ++j)
        {
          for (uint i = 0; i < n[0]; ++i)
          {
            cv[0] = j * (n[0] + 1) + i;
            cv[1] = cv[0] + co[1];
            cv[2] = cv[0] + co[2];
            cv[3] = cv[0] + co[3];
            cv[4] = cv[0] + co[4];
            cv[5] = cv[0] + co[5];
            editor.addCell(cell++, &cv[0]);
            editor.addCell(cell++, &cv[3]);
          }
        }
      }
      break;
    case CellType::tetrahedron:
      {
        uint const co[24] = { v0, v1, v3, v7,
                              v0, v1, v7, v5,
                              v0, v5, v7, v4,
                              v0, v3, v2, v7,
                              v0, v6, v4, v7,
                              v0, v2, v6, v7 };
        uint cv[24] = { 0 };
        editor.initCells(6 * num_bricks);
        for (uint k = 0; k < n[2]; ++k)
        {
          for (uint j = 0; j < n[1]; ++j)
          {
            for (uint i = 0; i < n[0]; ++i)
            {
              cv[0] = k * (n[1] + 1) * (n[0] + 1) + j * (n[0] + 1) + i;
              for (uint c = 1; c < 24; ++c)
              {
                cv[c] = cv[0] + co[c];
              }
              editor.addCell(cell++, &cv[0]);
              editor.addCell(cell++, &cv[4]);
              editor.addCell(cell++, &cv[8]);
              editor.addCell(cell++, &cv[12]);
              editor.addCell(cell++, &cv[16]);
              editor.addCell(cell++, &cv[20]);
            }
          }
        }
      }
      break;
    default:
      error("StructuredGrid::init : unsupported cell type");
      break;
    }

  //
  editor.close();
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */
