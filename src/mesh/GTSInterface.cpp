// Copyright (C) 2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/mesh/GTSInterface.h>

#include <dolfin/config/dolfin_config.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Point.h>
#include <dolfin/mesh/entities/Cell.h>
#include <dolfin/mesh/entities/Facet.h>
#include <dolfin/mesh/entities/Vertex.h>
#include <dolfin/mesh/entities/iterators/CellIterator.h>
#include <dolfin/mesh/entities/iterators/VertexIterator.h>
#include <dolfin/parameter/parameters.h>

#ifdef HAVE_GTS
#include <gts.h>
#endif

namespace dolfin
{

//-----------------------------------------------------------------------------

GTSInterface::GTSInterface( Mesh & mesh )
  : mesh_( mesh )
  , tree_( nullptr )
{
  if ( mesh_.geometry_dimension() > 3 )
  {
    error( "GTSInterface: not implemented for meshes of dimension %d.",
           mesh.geometry_dimension() );
  }

  buildCellTree();
}

//-----------------------------------------------------------------------------

GTSInterface::~GTSInterface()
{
#ifdef HAVE_GTS
  if ( tree_ != nullptr )
    gts_bb_tree_destroy( tree_, 1 );
#endif
}

//-----------------------------------------------------------------------------

auto GTSInterface::bboxCell( Cell & c ) const -> GtsBBox *
{
#ifdef HAVE_GTS
  VertexIterator v( c );
  real * p = v->x();

  GtsBBox *bbox = gts_bbox_new( gts_bbox_class(),
                                reinterpret_cast< gpointer >( c.index() ),
                                p[0], p[1], p[2],
                                p[0], p[1], p[2] );

  for ( ++v; !v.end(); ++v )
  {
    p = v->x();

    if ( p[0] > bbox->x2 )
      bbox->x2 = p[0];
    if ( p[0] < bbox->x1 )
      bbox->x1 = p[0];

    if ( p[1] > bbox->y2 )
      bbox->y2 = p[1];
    if ( p[1] < bbox->y1 )
      bbox->y1 = p[1];

    if ( p[2] > bbox->z2 )
      bbox->z2 = p[2];
    if ( p[2] < bbox->z1 )
      bbox->z1 = p[2];
  }
  return bbox;

#else
  error( "GTSInterface: Missing GTS" );
  return 0;
#endif
}

//-----------------------------------------------------------------------------

auto GTSInterface::bboxPoint( Point const & p ) const -> GtsBBox *
{
#ifdef HAVE_GTS

  real const btol = dolfin_get< real >( "GTS Tolerance" );
  GtsBBox * bbox  = gts_bbox_new( gts_bbox_class(),
                                  static_cast< gpointer >( nullptr ),
                                  p[0] - btol, p[1] - btol, p[2] - btol,
                                  p[0] + btol, p[1] + btol, p[2] + btol );

  return bbox;

#else
  error( "GTSInterface: Missing GTS" );
  return 0;
#endif
}

//-----------------------------------------------------------------------------

auto GTSInterface::bboxPoint( Point const & p1, Point const & p2 ) const
  -> GtsBBox *
{
#ifdef HAVE_GTS

  real x1 = p2[0];
  real y1 = p2[1];
  real z1 = p2[2];

  real x2 = p1[0];
  real y2 = p1[1];
  real z2 = p1[2];

  if ( p1[0] < p2[0] )
  {
    x1 = p1[0];
    x2 = p2[0];
  }

  if ( p1[1] < p2[1] )
  {
    y1 = p1[1];
    y2 = p2[1];
  }

  if ( p1[2] < p2[2] )
  {
    z1 = p1[2];
    z2 = p2[2];
  }

  GtsBBox * bbox = gts_bbox_new( gts_bbox_class(),
                                 static_cast< gpointer >( nullptr ),
                                 x1, y1, z1,
                                 x2, y2, z2 );

  return bbox;

#else
  error( "GTSInterface: Missing GTS" );
  return 0;
#endif
}

//-----------------------------------------------------------------------------

void GTSInterface::buildCellTree()
{
#ifdef HAVE_GTS

  if ( tree_ != nullptr )
  {
    error( "GTSInterface: tree already initialized" );
  }

  GSList * bboxes = nullptr;

  for ( CellIterator ci( mesh_ ); !ci.end(); ++ci )
  {
    Cell & c = *ci;
    bboxes   = g_slist_prepend( bboxes, bboxCell( c ) );
  }

  // A null pointer is returned if the mesh is empty
	if ( bboxes != nullptr )
  {
    tree_ = gts_bb_tree_new( bboxes );
    g_slist_free( bboxes );
  }

#else
  error( "GTSInterface: Missing GTS" );
#endif
}

//-----------------------------------------------------------------------------

void GTSInterface::overlap( Cell & c, std::vector< size_t > & cells ) const
{
#ifdef HAVE_GTS
	if ( tree_ == nullptr )
  {
    CellType const & type = mesh_.type();

    GtsBBox * bbprobe       = bboxCell( c );
    GSList *  overlaps      = gts_bb_tree_overlap( tree_, bbprobe );
    GSList *  overlaps_base = overlaps;

    while ( overlaps != nullptr )
    {
      GtsBBox * bb  = static_cast< GtsBBox * >( overlaps->data );
      size_t    idx = static_cast< size_t >( GPOINTER_TO_UINT( bb->bounded ) );

      Cell close( mesh_, idx );

      if ( type.intersects( c, close ) )
      {
        cells.push_back( idx );
      }

      overlaps = overlaps->next;
    }

    if ( overlaps_base != nullptr )
      g_slist_free( overlaps_base );

    if ( bbprobe != nullptr )
      gts_object_destroy( GTS_OBJECT( bbprobe ) );
  }

#else
  error( "GTSInterface: Missing GTS" );
#endif
}

//-----------------------------------------------------------------------------

void GTSInterface::overlap( Point const &           p,
                            std::vector< size_t > & cells ) const
{
#ifdef HAVE_GTS
	if ( tree_ != nullptr )
  {
    CellType const & type = mesh_.type();

    GtsBBox * bbprobe       = bboxPoint( p );
    GSList *  overlaps      = gts_bb_tree_overlap( tree_, bbprobe );
    GSList *  overlaps_base = overlaps;

    while ( overlaps != nullptr )
    {
      GtsBBox * bb  = static_cast< GtsBBox * >( overlaps->data );
      size_t    idx = static_cast< size_t >( GPOINTER_TO_UINT( bb->bounded ) );

      Cell close( mesh_, idx );

      if ( type.intersects( close, p ) )
      {
        cells.push_back( idx );
      }

      overlaps = overlaps->next;
    }

    if ( overlaps_base != nullptr )
      g_slist_free( overlaps_base );

    if ( bbprobe != nullptr )
      gts_object_destroy( GTS_OBJECT( bbprobe ) );
  }

#else
  error( "GTSInterface: Missing GTS" );
#endif
}

//-----------------------------------------------------------------------------

void GTSInterface::overlap( Point const &           p1,
                            Point const &           p2,
                            std::vector< size_t > & cells ) const
{
#ifdef HAVE_GTS
	if ( tree_ != nullptr )
  {
    CellType const & type = mesh_.type();

    GtsBBox * bbprobe       = bboxPoint( p1, p2 );
    GSList *  overlaps      = gts_bb_tree_overlap( tree_, bbprobe );
    GSList *  overlaps_base = overlaps;

    while ( overlaps )
    {
      GtsBBox * bb  = static_cast< GtsBBox * >( overlaps->data );
      size_t    idx = static_cast< size_t >( GPOINTER_TO_UINT( bb->bounded ) );

      Cell close( mesh_, idx );

      if ( type.intersects( close, p1, p2 ) )
        cells.push_back( idx );

      overlaps = overlaps->next;
    }

    if ( overlaps_base != nullptr )
      g_slist_free( overlaps_base );

    if ( bbprobe != nullptr )
      gts_object_destroy( GTS_OBJECT( bbprobe ) );
  }

#else
  error( "GTSInterface: Missing GTS" );
#endif
}

//-----------------------------------------------------------------------------

} // namespace dolfin
