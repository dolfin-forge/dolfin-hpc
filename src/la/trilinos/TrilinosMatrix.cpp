// Copyright (C) 2021 Julian Hornich
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_TRILINOS

#include <dolfin/la/trilinos/TrilinosMatrix.h>

#include <dolfin/la/trilinos/TrilinosFactory.h>

// FIXME lets say, this is not... optimal?!
#include <Tpetra_Access.hpp>
#include <Tpetra_applyDirichletBoundaryCondition.hpp>
#include <Tpetra_Apply_Helpers.hpp>
#include <Tpetra_ApplyOp.hpp>
#include <Tpetra_Assembly_Helpers.hpp>
#include <Tpetra_BlockCrsMatrix_decl.hpp>
#include <Tpetra_BlockCrsMatrix_def.hpp>
#include <Tpetra_BlockCrsMatrix_fwd.hpp>
#include <Tpetra_BlockCrsMatrix_Helpers_decl.hpp>
#include <Tpetra_BlockCrsMatrix_Helpers_def.hpp>
#include <Tpetra_BlockCrsMatrix_Helpers.hpp>
#include <Tpetra_BlockCrsMatrix.hpp>
#include <Tpetra_BlockMultiVector_decl.hpp>
#include <Tpetra_BlockMultiVector_def.hpp>
#include <Tpetra_BlockMultiVector_fwd.hpp>
#include <Tpetra_BlockMultiVector.hpp>
#include <Tpetra_BlockVector_decl.hpp>
#include <Tpetra_BlockVector_def.hpp>
#include <Tpetra_BlockVector_fwd.hpp>
#include <Tpetra_BlockVector.hpp>
#include <Tpetra_BlockView.hpp>
#include <Tpetra_CombineMode.hpp>
#include <Tpetra_ComputeGatherMap.hpp>
#include <Tpetra_computeRowAndColumnOneNorms_decl.hpp>
#include <Tpetra_computeRowAndColumnOneNorms_def.hpp>
#include <Tpetra_computeRowAndColumnOneNorms.hpp>
#include <Tpetra_ConfigDefs.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_createDeepCopy_CrsMatrix_decl.hpp>
#include <Tpetra_createDeepCopy_CrsMatrix_def.hpp>
#include <Tpetra_createDeepCopy_CrsMatrix.hpp>
#include <Tpetra_createDeepCopy.hpp>
#include <Tpetra_CrsGraph_decl.hpp>
#include <Tpetra_CrsGraph_def.hpp>
#include <Tpetra_CrsGraph_fwd.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_CrsMatrix_decl.hpp>
#include <Tpetra_CrsMatrix_def.hpp>
#include <Tpetra_CrsMatrix_fwd.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_CrsMatrixMultiplyOp_fwd.hpp>
#include <Tpetra_CrsMatrixMultiplyOp.hpp>
#include <Tpetra_CrsMatrix_UQ_PCE_def.hpp>
#include <Tpetra_CrsMatrix_UQ_PCE.hpp>
#include <Tpetra_Details_allReduceView.hpp>
#include <Tpetra_Details_Behavior.hpp>
#include <Tpetra_Details_Blas.hpp>
#include <Tpetra_Details_castAwayConstDualView.hpp>
#include <Tpetra_Details_checkGlobalError.hpp>
#include <Tpetra_Details_checkLaunchBlocking.hpp>
#include <Tpetra_Details_checkView.hpp>
#include <Tpetra_Details_computeOffsets.hpp>
#include <Tpetra_Details_CooMatrix.hpp>
#include <Tpetra_Details_copyConvert.hpp>
#include <Tpetra_Details_copyOffsets.hpp>
#include <Tpetra_Details_createMirrorView.hpp>
#include <Tpetra_Details_crsMatrixAssembleElement.hpp>
#include <Tpetra_Details_CrsPadding.hpp>
#include <Tpetra_Details_crsUtils.hpp>
#include <Tpetra_Details_DefaultTypes.hpp>
#include <Tpetra_Details_determineLocalTriangularStructure.hpp>
#include <Tpetra_Details_DualViewUtil.hpp>
#include <Tpetra_Details_EquilibrationInfo.hpp>
#include <Tpetra_Details_extractBlockDiagonal.hpp>
#include <Tpetra_Details_extractMpiCommFromTeuchos.hpp>
#include <Tpetra_Details_fill.hpp>
#include <Tpetra_Details_fill_MP_Vector.hpp>
#include <Tpetra_Details_fill_UQ_PCE.hpp>
#include <Tpetra_Details_FixedHashTable_decl.hpp>
#include <Tpetra_Details_FixedHashTable_def.hpp>
#include <Tpetra_Details_FixedHashTable.hpp>
#include <Tpetra_Details_gathervPrint.hpp>
#include <Tpetra_Details_getDiagCopyWithoutOffsets_decl.hpp>
#include <Tpetra_Details_getDiagCopyWithoutOffsets_def.hpp>
#include <Tpetra_Details_getDiagCopyWithoutOffsets.hpp>
#include <Tpetra_Details_getEntryOnHost.hpp>
#include <Tpetra_Details_getGraphDiagOffsets_decl.hpp>
#include <Tpetra_Details_getGraphDiagOffsets_def.hpp>
#include <Tpetra_Details_getGraphDiagOffsets.hpp>
#include <Tpetra_Details_getNumDiags.hpp>
#include <Tpetra_Details_Hash.hpp>
#include <Tpetra_Details_iallreduce.hpp>
#include <Tpetra_Details_initializeKokkos.hpp>
#include <Tpetra_Details_isInterComm.hpp>
#include <Tpetra_Details_lclDot.hpp>
#include <Tpetra_Details_leftScaleLocalCrsMatrix.hpp>
#include <Tpetra_Details_localDeepCopyRowMatrix_decl.hpp>
#include <Tpetra_Details_localDeepCopyRowMatrix_def.hpp>
#include <Tpetra_Details_localDeepCopyRowMatrix.hpp>
#include <Tpetra_Details_LocalMap_fwd.hpp>
#include <Tpetra_Details_LocalMap.hpp>
#include <Tpetra_Details_localRowOffsets_decl.hpp>
#include <Tpetra_Details_localRowOffsets_def.hpp>
#include <Tpetra_Details_localRowOffsets.hpp>
#include <Tpetra_Details_makeColMap_decl.hpp>
#include <Tpetra_Details_makeColMap_def.hpp>
#include <Tpetra_Details_makeColMap.hpp>
#include <Tpetra_Details_makeOptimizedColMap.hpp>
#include <Tpetra_Details_makeValidVerboseStream.hpp>
#include <Tpetra_Details_Merge.hpp>
#include <Tpetra_Details_mpiIsInitialized.hpp>
#include <Tpetra_Details_MpiTypeTraits.hpp>
#include <Tpetra_Details_normImpl.hpp>
#include <Tpetra_Details_OrdinalTraits.hpp>
#include <Tpetra_Details_packCrsGraph_decl.hpp>
#include <Tpetra_Details_packCrsGraph_def.hpp>
#include <Tpetra_Details_packCrsGraph.hpp>
#include <Tpetra_Details_packCrsMatrix_decl.hpp>
#include <Tpetra_Details_packCrsMatrix_def.hpp>
#include <Tpetra_Details_packCrsMatrix.hpp>
#include <Tpetra_Details_PackTraits.hpp>
#include <Tpetra_Details_PackTriples.hpp>
#include <Tpetra_Details_printOnce.hpp>
#include <Tpetra_Details_Profiling.hpp>
#include <Tpetra_Details_radixSort.hpp>
#include <Tpetra_Details_ReadTriples.hpp>
#include <Tpetra_Details_reallocDualViewIfNeeded.hpp>
#include <Tpetra_Details_residual.hpp>
#include <Tpetra_Details_rightScaleLocalCrsMatrix.hpp>
#include <Tpetra_Details_ScalarViewTraits.hpp>
#include <Tpetra_Details_scaleBlockDiagonal.hpp>
#include <Tpetra_Details_shortSort.hpp>
#include <Tpetra_Details_StaticView.hpp>
#include <Tpetra_Details_Transfer_decl.hpp>
#include <Tpetra_Details_Transfer_def.hpp>
#include <Tpetra_Details_Transfer_fwd.hpp>
#include <Tpetra_Details_Transfer.hpp>
#include <Tpetra_Details_unpackCrsGraphAndCombine_decl.hpp>
#include <Tpetra_Details_unpackCrsGraphAndCombine_def.hpp>
#include <Tpetra_Details_unpackCrsGraphAndCombine.hpp>
#include <Tpetra_Details_unpackCrsMatrixAndCombine_decl.hpp>
#include <Tpetra_Details_unpackCrsMatrixAndCombine_def.hpp>
#include <Tpetra_Details_unpackCrsMatrixAndCombine.hpp>
#include <Tpetra_Directory_decl.hpp>
#include <Tpetra_Directory_def.hpp>
#include <Tpetra_Directory_fwd.hpp>
#include <Tpetra_Directory.hpp>
#include <Tpetra_DirectoryImpl_decl.hpp>
#include <Tpetra_DirectoryImpl_def.hpp>
#include <Tpetra_DirectoryImpl.hpp>
#include <Tpetra_DistObject_decl.hpp>
#include <Tpetra_DistObject_def.hpp>
#include <Tpetra_DistObject_fwd.hpp>
#include <Tpetra_DistObject.hpp>
#include <Tpetra_Distribution.hpp>
#include <Tpetra_DistributionLowerTriangularBlock.hpp>
#include <Tpetra_DistributionMM.hpp>
#include <Tpetra_Distributor.hpp>
#include <Tpetra_EpetraRowMatrix.hpp>
#include <Tpetra_ETIHelperMacros.h>
#include <Tpetra_Exceptions_decl.hpp>
#include <Tpetra_Exceptions_def.hpp>
#include <Tpetra_Exceptions.hpp>
#include <Tpetra_Export_decl.hpp>
#include <Tpetra_Export_def.hpp>
#include <Tpetra_Export_fwd.hpp>
#include <Tpetra_Export.hpp>
#include <Tpetra_FECrsGraph_decl.hpp>
#include <Tpetra_FECrsGraph_def.hpp>
#include <Tpetra_FECrsGraph_fwd.hpp>
#include <Tpetra_FECrsGraph.hpp>
#include <Tpetra_FECrsMatrix_decl.hpp>
#include <Tpetra_FECrsMatrix_def.hpp>
#include <Tpetra_FECrsMatrix.hpp>
#include <Tpetra_FEMultiVector_decl.hpp>
#include <Tpetra_FEMultiVector_def.hpp>
#include <Tpetra_FEMultiVector_fwd.hpp>
#include <Tpetra_FEMultiVector.hpp>
#include <Tpetra_for_each.hpp>
#include <Tpetra_for_each_MultiVector.hpp>
#include <Tpetra_HashTable_decl.hpp>
#include <Tpetra_HashTable_def.hpp>
#include <Tpetra_HashTable.hpp>
#include <Tpetra_iallreduce.hpp>
#include <Tpetra_idot.hpp>
#include <Tpetra_Import_decl.hpp>
#include <Tpetra_Import_def.hpp>
#include <Tpetra_ImportExportData_decl.hpp>
#include <Tpetra_ImportExportData_def.hpp>
#include <Tpetra_ImportExportData_fwd.hpp>
#include <Tpetra_ImportExportData.hpp>
#include <Tpetra_Import_fwd.hpp>
#include <Tpetra_Import.hpp>
#include <Tpetra_Import_Util2.hpp>
#include <Tpetra_Import_Util.hpp>
#include <Tpetra_KokkosRefactor_Details_MultiVectorDistObjectKernels.hpp>
#include <Tpetra_KokkosRefactor_Details_MultiVectorDistObjectKernels_MP_Vector.hpp>
#include <Tpetra_KokkosRefactor_Details_MultiVectorDistObjectKernels_UQ_PCE.hpp>
#include <Tpetra_KokkosRefactor_Details_MultiVectorLocalDeepCopy.hpp>
#include <Tpetra_KokkosRefactor_Details_MultiVectorLocalDeepCopy_MP_Vector.hpp>
#include <Tpetra_KokkosRefactor_Details_MultiVectorLocalDeepCopy_UQ_PCE.hpp>
#include <Tpetra_leftAndOrRightScaleCrsMatrix_decl.hpp>
#include <Tpetra_leftAndOrRightScaleCrsMatrix_def.hpp>
#include <Tpetra_leftAndOrRightScaleCrsMatrix.hpp>
#include <Tpetra_LocalCrsMatrixOperator_decl.hpp>
#include <Tpetra_LocalCrsMatrixOperator_def.hpp>
#include <Tpetra_LocalCrsMatrixOperator_fwd.hpp>
#include <Tpetra_LocalCrsMatrixOperator.hpp>
#include <Tpetra_LocalOperator_fwd.hpp>
#include <Tpetra_LocalOperator.hpp>
#include <Tpetra_Map_decl.hpp>
#include <Tpetra_Map_def.hpp>
#include <Tpetra_Map_fwd.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MatrixIO_decl.hpp>
#include <Tpetra_MatrixIO_def.hpp>
#include <Tpetra_MatrixIO.hpp>
#include <Tpetra_MultiVector_decl.hpp>
#include <Tpetra_MultiVector_def.hpp>
#include <Tpetra_MultiVector_fwd.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Operator_fwd.hpp>
#include <Tpetra_Operator.hpp>
#include <Tpetra_Packable_fwd.hpp>
#include <Tpetra_Packable.hpp>
#include <Tpetra_replaceDiagonalCrsMatrix_decl.hpp>
#include <Tpetra_replaceDiagonalCrsMatrix_def.hpp>
#include <Tpetra_replaceDiagonalCrsMatrix.hpp>
#include <Tpetra_RowGraph_decl.hpp>
#include <Tpetra_RowGraph_def.hpp>
#include <Tpetra_RowGraph_fwd.hpp>
#include <Tpetra_RowGraph.hpp>
#include <Tpetra_RowMatrix_decl.hpp>
#include <Tpetra_RowMatrix_def.hpp>
#include <Tpetra_RowMatrix_fwd.hpp>
#include <Tpetra_RowMatrix.hpp>
#include <Tpetra_RowMatrixTransposer_decl.hpp>
#include <Tpetra_RowMatrixTransposer_def.hpp>
#include <Tpetra_RowMatrixTransposer_fwd.hpp>
#include <Tpetra_RowMatrixTransposer.hpp>
#include <Tpetra_SrcDistObject.hpp>
#include <Tpetra_TieBreak_fwd.hpp>
#include <Tpetra_TieBreak.hpp>
#include <Tpetra_transform.hpp>
#include <Tpetra_transform_MultiVector.hpp>
#include <Tpetra_TsqrAdaptor.hpp>
#include <Tpetra_TsqrAdaptor_MP_Vector.hpp>
#include <Tpetra_Util.hpp>
#include <Tpetra_Vector_decl.hpp>
#include <Tpetra_Vector_def.hpp>
#include <Tpetra_Vector_fwd.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_Version.hpp>
#include <Tpetra_withLocalAccess.hpp>
#include <Tpetra_withLocalAccess_MultiVector.hpp>

namespace dolfin
{

namespace trilinos
{

//-----------------------------------------------------------------------------

Matrix::Matrix()
  : Variable( "A", "a sparse Matrix" )
	, mat_( nullptr )
{
  // do nothing
}

//-----------------------------------------------------------------------------

Matrix::Matrix( size_t M, size_t N, bool distributed )
  : Variable( "A", "a sparse Matrix" )
	, mat_( nullptr )
{
	init( M, N, distributed );
}

//-----------------------------------------------------------------------------

Matrix::Matrix( Matrix const & A )
  : Variable( "A", "a sparse Matrix" )
	, mat_( nullptr )
{
	// delegate to assignment operator
  *this = A;
}

//-----------------------------------------------------------------------------

Matrix::~Matrix()
{
}

//-----------------------------------------------------------------------------

auto Matrix::copy() const -> Matrix *
{
  return new trilinos::Matrix( *this );
}

//-----------------------------------------------------------------------------

auto Matrix::size( size_t dim ) const -> size_t
{
  // FIXME
  return DOLFIN_SIZE_T_MAX;
}

//-----------------------------------------------------------------------------

auto Matrix::zero() -> void
{
  dolfin_assert( not mat_.is_null() );

  if ( mat_->isFillComplete() )
    mat_->resumeFill();

  mat_->setAllToScalar( 0.0 );
}

//-----------------------------------------------------------------------------

auto Matrix::apply( FinalizeType finaltype ) -> void
{
  dolfin_assert( not mat_.is_null() );

  if (  finaltype == FINALIZE or finaltype == FLUSH )
  {
    mat_->fillComplete();
  }
  else
  {
    error( "trilinos::Matrix: Unknown apply mode" );
  }
}

//-----------------------------------------------------------------------------

auto Matrix::disp( size_t precision ) const -> void
{
  // FIXME
}

//-----------------------------------------------------------------------------

auto Matrix::init( const GenericSparsityPattern & sparsity_pattern ) -> void
{
  SparsityPattern const & spattern = reinterpret_cast< SparsityPattern const & >( sparsity_pattern );

  size_t const nLocalRows = spattern.size( 0 );
  size_t const nLocalCols = spattern.size( 1 );

  init( nLocalRows, nLocalCols, true );
}

//-----------------------------------------------------------------------------

auto Matrix::init( size_t M, size_t N ) -> void
{
  init( M, N, true );
}

//-----------------------------------------------------------------------------

auto Matrix::init( size_t M, size_t N, bool ) -> void
{
  size_t const pe_rank = MPI::rank();
  size_t const pe_size = MPI::size();

  size_t const nLocalRows = M;
  size_t const nLocalCols = N;

  std::vector< size_t > globalRows( pe_size, 0 );
  std::vector< size_t > globalCols( pe_size, 0 );

  MPI::all_gather( nLocalRows, globalRows );
  MPI::all_gather( nLocalCols, globalCols );

  size_t const nGlobalRows = std::accumulate( globalRows.begin(), globalRows.end(), 0);
  size_t const nGlobalCols = std::accumulate( globalCols.begin(), globalCols.end(), 0);

  // Create the row and column index lists on each processor
  std::vector< GO > rowIndices( nLocalRows, 0 );
  std::vector< GO > colIndices( nLocalCols, 0 );
  {
    GO rowStart = std::accumulate( globalRows.begin(), globalRows.begin() + pe_rank, 0 );
    std::iota(rowIndices.begin(), rowIndices.end(), rowStart);

    GO colStart = std::accumulate( globalCols.begin(), globalCols.begin() + pe_rank, 0 );
    std::iota(colIndices.begin(), colIndices.end(), colStart);
  }

  using Teuchos::MpiComm;
  using Teuchos::rcp;
  using Teuchos::RCP;

  RCP< MpiComm< int > const > comm( new MpiComm< int >( MpiComm< int >( MPI::DOLFIN_COMM ) ) );

  // Create the row map
  RCP< TPMap const > rowMap = rcp( new TPMap( nGlobalRows, rowIndices.data(),
                                              nLocalRows, indexBase, comm ) );

  // Create the column map
  RCP< TPMap const > colMap = rcp( new TPMap( nGlobalCols, colIndices.data(),
                                              nLocalCols, indexBase, comm ) );

  // Create a Tpetra sparse matrix whose rows have distribution
  // given by the row Map and column Map.
  RCP< TPMatrix > A( new TPMatrix( rowMap, colMap, 0 ) );
}

//-----------------------------------------------------------------------------

auto Matrix::get( real *         block,
                  size_t         m,
                  const size_t * rows,
                  size_t         n,
                  const size_t * cols ) const -> void
{
  // FIXME
  error( "trilinos::Matrix: get not implemented." );
}

//-----------------------------------------------------------------------------

auto Matrix::set( const real *   block,
                  size_t         m,
                  const size_t * rows,
                  size_t         n,
                  const size_t * cols ) -> void
{
  dolfin_assert( not mat_.is_null() );
  dolfin_assert( not mat_->isFillComplete() );

  // Tpetra View of column indices
  Teuchos::ArrayView< GO const > column_idx( cols, n );
  for ( size_t i = 0; i < m; ++i )
  {
    Teuchos::ArrayView< real const > data( block + i * n, n );
    mat_->replaceGlobalValues( rows[i], column_idx, data );
  }
}

//-----------------------------------------------------------------------------

auto Matrix::add( const real *   block,
                  size_t         m,
                  const size_t * rows,
                  size_t         n,
                  const size_t * cols ) -> void
{
  dolfin_assert( not mat_.is_null() );
  dolfin_assert( not mat_->isFillComplete() );

  // Tpetra View of column indices
  Teuchos::ArrayView< GO const > column_idx( cols, n );
  for ( size_t i = 0; i < m; ++i )
  {
    Teuchos::ArrayView< real const > data( block + i * n, n );
    mat_->sumIntoGlobalValues( rows[i], column_idx, data );
  }
}
//-----------------------------------------------------------------------------

auto Matrix::norm( std::string norm_type ) const -> real
{
  real norm = 0.0;

  if ( norm_type == "l1" )
  {
    error( "trilinos::Matrix: unimplemented L1 norm." );
  }
  else if ( norm_type == "linf" )
  {
    error( "trilinos::Matrix: unimplemented inf norm." );
  }
  else if ( norm_type == "frobenius" )
  {
    norm = mat_->getFrobeniusNorm();
  }
  else
  {
    error( "trilinos::Matrix: Unknown norm type." );
  }

  return norm;
}

//-----------------------------------------------------------------------------

auto Matrix::getrow( size_t                  row,
                     std::vector< size_t > & columns,
                     std::vector< real > &   values ) const -> void
{
  dolfin_assert( not mat_.is_null() );

  size_t const ncols = mat_->getNumEntriesInGlobalRow( row );
  if ( ncols == Teuchos::OrdinalTraits< size_t >::invalid() )
  {
    error( "trilinos::Matrix: Row %d not in range", row );
  }

  columns.resize( ncols );
  values.resize( ncols );

  Teuchos::ArrayView< GO >   _columns( columns );
  Teuchos::ArrayView< real > _values( values );

  size_t n = 0;
  mat_->getGlobalRowCopy( row, _columns, _values, n );

  dolfin_assert( n == ncols );
}

//-----------------------------------------------------------------------------

auto Matrix::setrow( size_t                        row,
                     const std::vector< size_t > & columns,
                     const std::vector< real > &   values ) -> void
{
  dolfin_assert( not mat_.is_null() );
  dolfin_assert( not mat_->isFillComplete() );

  if ( columns.size() != values.size() )
  {
    error( "trilinos::Matrix: Number of columns and values don't match" );
  }

  // Handle case n = 0
  if ( columns.size() != 0 )
  {
    // Tpetra View of column indices
    Teuchos::ArrayView< GO const> column_idx( columns );

    // Tpetra View of values
    Teuchos::ArrayView< real const > data( values );

    mat_->replaceGlobalValues( row, column_idx, data );
  }
}

//-----------------------------------------------------------------------------

auto Matrix::zero( size_t m, const size_t * rows ) -> void
{
  dolfin_assert( not mat_.is_null() );
  dolfin_assert( not mat_->isFillComplete() );

  for ( size_t i = 0; i < m; ++i )
  {
    size_t const ncols = mat_->getNumEntriesInGlobalRow( rows[i] );

    std::vector< GO >        colcols( ncols );
    Teuchos::ArrayView< GO > cols( colcols );

    std::vector< real >        coldata( ncols );
    Teuchos::ArrayView< real > data( coldata );

    size_t n = 0;
    mat_->getGlobalRowCopy( rows[i], cols, data, n );
    dolfin_assert( n == ncols );

    std::fill( coldata.begin(), coldata.end(), 0.0 );
    mat_->replaceGlobalValues( rows[i], cols, data );
  }
}
//-----------------------------------------------------------------------------

auto Matrix::ident( size_t m, const size_t * rows ) -> void
{
  dolfin_assert( not mat_.is_null() );

  if ( mat_->isFillComplete() )
    mat_->resumeFill();

  // Clear affected rows to zero
  zero( m, rows );

  // Get map of locally available columns
  Teuchos::RCP< TPMap const > colmap( mat_->getColMap() );

  real const one = 1;
  GO         col = 0;

  Teuchos::ArrayView< real const > data( &one, 1 );
  Teuchos::ArrayView< GO >         column_idx( &col, 1 );

  // Set diagonal entries where possible
  for ( size_t i = 0; i < m; ++i )
  {
    if ( colmap->isNodeGlobalElement( rows[i] ) )
    {
      col = rows[i];
      mat_->replaceGlobalValues( rows[i], column_idx, data );
    }
  }
}

//-----------------------------------------------------------------------------

auto Matrix::mult( const GenericVector & x,
                   GenericVector &       y,
                   bool                  transposed ) const -> void
{
  dolfin_assert( not mat_.is_null() );

  trilinos::Vector const & X = x.down_cast< trilinos::Vector const >();
  trilinos::Vector &       Y = y.down_cast< trilinos::Vector >();

  if ( not transposed )
  {

    if ( size( 1 ) != X.size() )
    {
      error( "trilinos::Matrix: Non-matching dimensions %d and %d for matrix-vector product",
             size( 1 ), X.size() );
    }

    // // Resize RHS if empty
    // if ( Y.size() == 0 )
    //   init_vector( Y, 0 );
    dolfin_assert( not Y.x_.is_null() );

    if ( size( 0 ) != Y.size() )
    {
      error( "trilinos::Matrix: Vector for matrix-vector result has wrong size" );
    }

    mat_->apply( *X.x_, *Y.x_ );
  }
  else // transposed
  {
    if ( size( 0 ) != X.size() )
    {
      error(
        "TpetraMatrix.cpp",
        "compute transpose matrix-vector product with Tpetra matrix",
        "Non-matching dimensions for transpose matrix-vector product" );
    }

    // // Resize RHS if empty
    // if ( Y.size() == 0 )
    //   init_vector( Y, 1 );
    dolfin_assert( not Y.x_.is_null() );

    if ( size( 1 ) != Y.size() )
    {
      error( "trilinos::Matrix: Vector for transpose matrix-vector result has wrong size" );
    }

    mat_->apply( *X.x_, *Y.x_ghosted_, Teuchos::TRANS );
  }
}

//-----------------------------------------------------------------------------

auto Matrix::operator*=( real a ) -> const Matrix &
{
  dolfin_assert( not mat_.is_null() );

  if ( mat_->isFillComplete() )
    mat_->resumeFill();

  mat_->scale( a );
  return *this;
}

//-----------------------------------------------------------------------------

auto Matrix::operator/=( real a ) -> const Matrix &
{
  dolfin_assert( not mat_.is_null() );

  if ( mat_->isFillComplete() )
    mat_->resumeFill();

  mat_->scale( 1.0 / a );
  return *this;
}

//-----------------------------------------------------------------------------

auto Matrix::operator=( GenericMatrix const & A ) -> const GenericMatrix &
{
  *this = A.down_cast< trilinos::Matrix const >();
  return *this;
}

//-----------------------------------------------------------------------------

auto Matrix::nz() const -> size_t
{
  // FIXME
  return DOLFIN_SIZE_T_MAX;
}

//-----------------------------------------------------------------------------

auto Matrix::factory() const -> LinearAlgebraFactory &
{
  return trilinos::Factory::instance();
}

//-----------------------------------------------------------------------------

auto Matrix::norm( const Norm type ) const -> real
{
  real norm = 0.0;

  switch ( type )
  {
    case l1:
      error( "trilinos::Matrix: unimplemented L1 norm." );
      break;
    case linf:
      error( "trilinos::Matrix: unimplemented inf norm." );
      break;
    case frobenius:
      norm = mat_->getFrobeniusNorm();
      break;
    default:
      error( "Unknown norm type." );
      break;
  }

  return norm;
}

//-----------------------------------------------------------------------------

auto Matrix::operator=( const Matrix & A ) -> const Matrix &
{
  mat_ = Teuchos::RCP< TPMatrix >( new TPMatrix( A.mat_->getRowMap(),
                                                 A.mat_->getColMap(),
                                                 0 ) );

  return *this;
}

//-----------------------------------------------------------------------------

auto Matrix::operator+=( const Matrix & A ) -> const Matrix &
{
  dolfin_assert( not mat_.is_null() );

  trilinos::Matrix const & AA = A.down_cast< trilinos::Matrix const >();
  dolfin_assert( not AA.mat_.is_null() );

  real const one = 1.0;

  mat_ = Teuchos::rcp_dynamic_cast< TPMatrix >( mat_->add( one, *AA.mat_, one,
                                                            Teuchos::null,
                                                            Teuchos::null,
                                                            Teuchos::null ) );

  return *this;
}

//-----------------------------------------------------------------------------

} // end namespace trilinos

} // end namespace dolfin

#endif // HAVE_TRILINOS
