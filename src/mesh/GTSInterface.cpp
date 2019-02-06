// Copyright (C) 2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Johan Jansson 2006.
// Modified by Ola Skavhaug 2006.
// Modified by Dag Lindbo 2008.
//
// First added:  2006-06-21
// Last changed: 2006-12-01

#include <dolfin/config/dolfin_config.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/common/Array.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Cell.h>

#include <dolfin/mesh/GTSInterface.h>
#include <dolfin/parameter/parameters.h>

#ifdef HAVE_GTS
#include <gts.h>
#endif

namespace dolfin
{

//-----------------------------------------------------------------------------
GTSInterface::GTSInterface(Mesh& mesh) :
    mesh_(mesh),
    tree_(NULL)
{
  if (mesh.geometry_dimension() > 3)
  {
    error("Sorry, GTS interface not implemented for meshes of dimension %d.",
          mesh.geometry_dimension());
  }

  buildCellTree();
}
//-----------------------------------------------------------------------------
GTSInterface::~GTSInterface()
{
#ifdef HAVE_GTS
  gts_bb_tree_destroy(tree_, 1);
#endif
}
//-----------------------------------------------------------------------------
GtsBBox* GTSInterface::bboxCell(Cell& c) const
{
#ifdef HAVE_GTS
  GtsBBox* bbox;
  Point p;

  VertexIterator v(c);
  p = v->point();

  bbox = gts_bbox_new(gts_bbox_class(), (gpointer) c.index(), p[0], p[1], p[2],
                      p[0], p[1], p[2]);

  for (++v; !v.end(); ++v)
  {
    p = v->point();
    if (p[0] > bbox->x2) bbox->x2 = p[0];
    if (p[0] < bbox->x1) bbox->x1 = p[0];
    if (p[1] > bbox->y2) bbox->y2 = p[1];
    if (p[1] < bbox->y1) bbox->y1 = p[1];
    if (p[2] > bbox->z2) bbox->z2 = p[2];
    if (p[2] < bbox->z1) bbox->z1 = p[2];
  }
  return bbox;

#else
  error("Missing GTS");
  return 0;
#endif
}
//-----------------------------------------------------------------------------
GtsBBox* GTSInterface::bboxPoint(Point const& p) const
{
#ifdef HAVE_GTS

  GtsBBox* bbox;

  real btol = dolfin_get("GTS Tolerance");
  bbox = gts_bbox_new(gts_bbox_class(), (gpointer) NULL, p[0] - btol, p[1] - btol,
                      p[2] - btol, p[0] + btol, p[1] + btol, p[2] + btol);

  return bbox;

#else
  error("Missing GTS");
  return 0;
#endif
}
//-----------------------------------------------------------------------------
GtsBBox* GTSInterface::bboxPoint(Point const& p1, Point const& p2) const
{
#ifdef HAVE_GTS

  GtsBBox* bbox;

  real x1, x2;
  real y1, y2;
  real z1, z2;

  if (p1[0] < p2[0])
  {
    x1 = p1[0];
    x2 = p2[0];
  }
  else
  {
    x1 = p2[0];
    x2 = p1[0];
  }
  if (p1[1] < p2[1])
  {
    y1 = p1[1];
    y2 = p2[1];
  }
  else
  {
    y1 = p2[1];
    y2 = p1[1];
  }
  if (p1[2] < p2[2])
  {
    z1 = p1[2];
    z2 = p2[2];
  }
  else
  {
    z1 = p2[2];
    z2 = p1[2];
  }

  bbox = gts_bbox_new(gts_bbox_class(), (gpointer) NULL, x1, y1, z1, x2, y2, z2);

  return bbox;

#else
  error("Missing GTS");
  return 0;
#endif
}
//-----------------------------------------------------------------------------
void GTSInterface::buildCellTree()
{
#ifdef HAVE_GTS

  if (tree_)
  {
    error("GTS tree already initialized");
  }

  GSList* bboxes = NULL;

  for (CellIterator ci(mesh_); !ci.end(); ++ci)
  {
    Cell& c = *ci;
    bboxes = g_slist_prepend(bboxes, bboxCell(c));
  }

  // A null pointer is returned if the mesh is empty
  tree_ = gts_bb_tree_new(bboxes);
  g_slist_free(bboxes);

#else
  error("Missing GTS");
#endif
}
//-----------------------------------------------------------------------------
void GTSInterface::overlap(Cell& c, Array<uint>& cells) const
{
#ifdef HAVE_GTS
  GtsBBox* bbprobe;
  GtsBBox* bb;
  GSList* overlaps = 0, *overlaps_base;
  uint boundedcell;

  CellType const& type = mesh_.type();

  bbprobe = bboxCell(c);

  overlaps = gts_bb_tree_overlap(tree_, bbprobe);
  overlaps_base = overlaps;

  while (overlaps)
  {
    bb = (GtsBBox *) overlaps->data;
    boundedcell = (uint)(long) bb->bounded;

	Cell close(mesh_, boundedcell);
    if (type.intersects(c, close))
    {
      cells.push_back(boundedcell);
    }
    overlaps = overlaps->next;
  }

  g_slist_free(overlaps_base);
  gts_object_destroy(GTS_OBJECT(bbprobe));

#else
  error("Missing GTS");
#endif
}
//-----------------------------------------------------------------------------
void GTSInterface::overlap(Point const& p, Array<uint>& cells) const
{
#ifdef HAVE_GTS
  GtsBBox* bbprobe;
  GtsBBox* bb;
  GSList* overlaps = 0, *overlaps_base;
  uint boundedcell;

  CellType const& type = mesh_.type();

  bbprobe = bboxPoint(p);

  overlaps = gts_bb_tree_overlap(tree_, bbprobe);
  overlaps_base = overlaps;

  while (overlaps)
  {
    bb = (GtsBBox *) overlaps->data;
    boundedcell = (uint)(long)bb->bounded;

	Cell close(mesh_, boundedcell);

    if (type.intersects(close, p))
    {
      cells.push_back(boundedcell);
    }

    overlaps = overlaps->next;
  }

  g_slist_free(overlaps_base);
  gts_object_destroy(GTS_OBJECT(bbprobe));

#else
  error("Missing GTS");
#endif
}
//-----------------------------------------------------------------------------
void GTSInterface::overlap(Point const& p1, Point const& p2,
                           Array<uint>& cells) const
{
#ifdef HAVE_GTS
  GtsBBox* bbprobe;
  GtsBBox* bb;
  GSList* overlaps = 0, *overlaps_base;
  uint boundedcell;

  CellType const& type = mesh_.type();

  bbprobe = bboxPoint(p1, p2);

  overlaps = gts_bb_tree_overlap(tree_, bbprobe);
  overlaps_base = overlaps;

  while (overlaps)
  {
    bb = (GtsBBox *) overlaps->data;
    boundedcell = (uint) (long) bb->bounded;

    Cell close(mesh_, boundedcell);

    if (type.intersects(close, p1, p2)) cells.push_back(boundedcell);

    overlaps = overlaps->next;
  }
  g_slist_free(overlaps_base);
  gts_object_destroy(GTS_OBJECT(bbprobe));

#else
  error("Missing GTS");
#endif
}
//-----------------------------------------------------------------------------

}

