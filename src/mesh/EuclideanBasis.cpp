
#include <dolfin/mesh/EuclideanBasis.h>

#include <dolfin/math/basic.h>

namespace dolfin
{

//-----------------------------------------------------------------------------

auto EuclideanBasis::compute( size_t                      gdim,
                              Point                       B[],
                              std::vector< real > const & N,
                              std::vector< real > const & W,
                              real                        cosalpha_max,
                              bool                        weighted ) -> size_t
{
  if ( N.size() != W.size() * gdim )
  {
    error( "EuclideanBasis : size mismatch for normals (%u) and weights (%u)",
           N.size(),
           W.size() );
  }

  //--- Determine vertex type by discriminating surfaces ----------------
  std::vector< Point >             nS;
  std::vector< real >              wS;
  _ordered_set< size_t >           Rnormals;
  size_t const                     num_facets = W.size();
  _ordered_set< size_t >::iterator it         = Rnormals.begin();
  real                             wSa        = 0.0;
  for ( size_t i = 0; i < num_facets; ++i )
  {
    Rnormals.insert( it, i );
    wSa += W[i];
    dolfin_assert( W[i] > 0.0 );
  }

  while ( Rnormals.size() > 0 )
  {
    // Initialize new surface normal and weight with reference
    it                  = Rnormals.begin();
    size_t const rfacet = ( *it );
    Point        nSx;
    real         wSx = W[rfacet];
    for ( size_t d = 0; d < gdim; ++d )
    {
      nSx[d] = wSx * N[gdim * rfacet + d];
    }

    _ordered_set< size_t > Unormals;
    Unormals.insert( Unormals.begin(), rfacet );
    // Loop through remaining normals indexes
    for ( ++it; it != Rnormals.end(); ++it )
    {
      size_t const cfacet = ( *it );
      dolfin_assert( W[cfacet] > 0 );

      // Compute the scalar product with the reference normal
      real cosalpha = 0.0;
#ifdef _CRAYC
#if _RELEASE_MAJOR == 8 && _RELEASE_MINOR < 3
#pragma _CRI novector
#endif
#endif
      for ( size_t d = 0; d < gdim; ++d )
      {
        cosalpha += N[gdim * rfacet + d] * N[gdim * cfacet + d];
      }
      if ( cosalpha > cosalpha_max )
      {
        // Add contribution to surface normal
        for ( size_t d = 0; d < gdim; ++d )
        {
          nSx[d] += W[cfacet] * N[gdim * cfacet + d];
        }
        wSx += W[cfacet];
        // Take into account that the normal is used
        Unormals.insert( Unormals.begin(), cfacet );
      }
    }

    // Add unit surface normal to the list of surface normals
    // message("Surface normal %i", vertex_type);
    nSx /= nSx.norm();
    nS.push_back( nSx );
    if ( weighted )
    {
      wSx /= wSa;
    }
    else
    {
      wSx = 1.0;
    }
    wS.push_back( wSx );

    // Next loop we add a new surface and discriminate again against
    // the remaining normals
    for ( it = Unormals.begin(); it != Unormals.end(); ++it )
    {
      Rnormals.erase( *it );
    }
  }
  dolfin_assert( wS.size() > 0 );

  //--- Compute vertex normal -------------------------------------------
  // The typical algorithm would be:
  // n_k    = sum_{i=1}^k nS_i
  // tau1_k = |n_{k-1}|^2 nS_k - (n_{k-1} . nS_k ) n_{k-1}
  // tau2_k = n_k ^ tau1_k
  // and such that in 2d tau2_k = ez = (0 , 0, 1)
  for ( size_t d = 0; d < gdim; ++d )
  {
    B[0][d] = 0.0;
    for ( size_t s = 0; s < nS.size(); ++s )
    {
      B[0][d] += wS[s] * nS[s][d];
    }
  }
  B[0] /= B[0].norm();

  switch ( gdim )
  {
    case 1:
      break;
    case 2:
      B[1][0] = -B[0][1];
      B[1][1] = +B[0][0];
      break;
    case 3:
      if ( nS.size() == 1 )
      {
        if ( std::fabs( B[0][0] ) >= 0.5 || std::fabs( B[0][1] ) >= 0.5 )
        {
          // t1 = rotation in (x, y)
          B[1][0] = -B[0][1];
          B[1][1] = +B[0][0];
          B[1][2] = +0.0;
          B[1] /= B[1].norm();

          // t2 = n ^ t1
          B[2][0] = -B[0][2] * B[1][1];
          B[2][1] = +B[0][2] * B[1][0];
          B[2][2] = +B[0][0] * B[1][1] - B[0][1] * B[1][0];
          // B[2] is de facto normalized
        }
        else
        {
          // t1 = rotation in (y, z)
          B[1][0] = +0.0;
          B[1][1] = -B[0][2];
          B[1][2] = +B[0][1];
          B[1] /= B[1].norm();

          // t2 = n ^ t1
          B[2][0] = +B[0][1] * B[1][2] - B[0][2] * B[1][1];
          B[2][1] = -B[0][0] * B[1][2];
          B[2][2] = +B[0][0] * B[1][1];
          // B[2] is de facto normalized
        }
      }
      else
      {
        size_t const k = nS.size() - 1;
        // t2 = n ^ nSk / || n ^ nSk ||
        B[2][0] = +B[0][1] * nS[k][2] - nS[k][1] * B[0][2];
        B[2][1] = +B[0][2] * nS[k][0] - nS[k][2] * B[0][0];
        B[2][2] = +B[0][0] * nS[k][1] - nS[k][0] * B[0][1];
        B[2] /= B[2].norm();

        // t1 = t2 ^ n
        B[1][0] = +B[2][1] * B[0][2] - B[2][2] * B[0][1];
        B[1][1] = +B[2][2] * B[0][0] - B[2][0] * B[0][2];
        B[1][2] = +B[2][0] * B[0][1] - B[2][1] * B[0][0];
        // B[1] is de facto normalized
      }
      break;
    default:
      error( "Unsupported geometric dimension." );
      break;
  }

#if DEBUG
  for ( size_t i = 0; i < gdim; ++i )
  {
    for ( size_t d = 0; d < gdim; ++d )
    {
      if ( B[i][d] != B[i][d] )
      {
        error( "Component %d of vector %d is Not-a-Number.", d, i );
      }
    }
    real en = B[i].norm();
    if ( !abscmp( en, 1.0, gdim * DOLFIN_EPS ) )
    {
      error( "Basis is not normal: ||e%u|| = %e", i, en );
    }
    size_t j  = ( i + 1 ) % gdim;
    real   sp = B[i].dot( B[j] );
    if ( not abscmp( sp, 0.0, 1e-13 ) )
    {
      error( "Basis is not orthogonal: e%u . e%u = %e", i, j, sp );
    }
  }
#endif

  return nS.size();
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */
