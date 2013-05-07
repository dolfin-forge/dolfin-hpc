// Copyright (C) 2013 Balthasar Reuter.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2013-03-25
// Last changed: 2013-04-03

#ifndef __COARSENING_MANAGER_H
#define __COARSENING_MANAGER_H

#include <dolfin/common/types.h>
#include <dolfin/common/Array.h>
#include <dolfin/common/List.h>
#include <dolfin/mesh/MeshFunction.h>

#define ____USE_D_MESH____            // Activate usage of dynamic mesh

namespace dolfin
{
  class Mesh;
  class Vertex;

#ifdef ____USE_D_MESH____
  class DMesh;
  class DVertex;
  class DCell;
#endif

  /// Assists LocalMeshCoarsening by providing relevant information about mesh
  /// entities. 
  class CoarseningManager
  {
  public:
    class IndexMap;

    /// Create an empty instance, initialize it by calling init()
    CoarseningManager();

    /// Create a new instance and initialize it.
    explicit CoarseningManager(MeshFunction<bool>& cell_marker, 
                               bool coarsen_boundary = false);

    ~CoarseningManager();

    /// Initialize the coarsening manager: build independent set, list of cells
    /// that are marked for coarsening and boundary information
    void init(MeshFunction<bool>& cell_marker, bool coarsen_boundary = false);

    /// Check if a vertex is on any (interior or domain) boundary
    inline bool isBoundaryVertex(uint index)
    { return _int_bnd_vertices.at(index) || _bnd_vertices.at(index); }

    /// Check if a cell is on any (interior or domain) boundary
    inline bool isBoundaryCell(uint index)
    { return _int_bnd_cells.at(index) || _bnd_cells.at(index); }

    /// Check if a vertex is on the interior boundary
    inline bool isInteriorBoundaryVertex(uint index)
    { return _int_bnd_vertices.at(index); }

    /// Check if a cell is on the interior boundary
    inline bool isInteriorBoundaryCell(uint index)
    { return _int_bnd_cells.at(index); }

    /// Check if a vertex is on the domain boundary
    inline bool isDomainBoundaryVertex(uint index)
    { return _bnd_vertices.at(index); }

    /// Check if a cell is on the domain boundary
    inline bool isDomainBoundaryCell(uint index)
    { return _bnd_cells.at(index); }

    /// Check if a vertex is forbidden, i. e. part of the independent set
    inline bool isForbiddenVertex(uint index)
    { return _forbidden_vertices.at(index); }

#ifdef ____USE_D_MESH____
    bool checkDCellNumbering(uint max_index);

    /// Gives access to the dynamic mesh
    inline DMesh * dmesh()
    { return _dmesh; }

    /// Gives access to the dynamic mesh
    inline DMesh const * dmesh() const
    { return _dmesh; }
#endif

#ifndef ____USE_D_MESH____
    /// Gives access to the IndexMap for cells
    inline IndexMap& cell_map()
    { return _cell_map; }

    /// Gives access to the IndexMap for cells
    inline IndexMap const & cell_map() const
    { return _cell_map; }

    /// Gives access to the IndexMap for vertices
    inline IndexMap& vertex_map()
    { return _vertex_map; }

    /// Gives access to the IndexMap for vertices
    inline IndexMap const & vertex_map() const
    { return _vertex_map; }
#endif

#ifdef ____USE_D_MESH____
    /// Gives access to the list of cells for coarsening
    inline List<DCell *>& cells_to_coarsen()
    { return _cells_to_coarsen; }

    /// Gives access to the list of cells for coarsening
    inline List<DCell *> const & cells_to_coarsen() const
    { return _cells_to_coarsen; }
#else // ____USE_D_MESH____
    /// Gives access to the list of cells for coarsening
    inline List<uint>& cells_to_coarsen()
    { return _cells_to_coarsen; }

    /// Gives access to the list of cells for coarsening
    inline List<uint> const & cells_to_coarsen() const
    { return _cells_to_coarsen; }
#endif // ____USE_D_MESH____

    /// Gives access to the list of cells that have been attempted to be 
    /// coarsened but failed because of missing data from other processes
    inline List<uint>& cells_to_request()
    { return _cells_to_request; }

    /// Gives access to the list of cells that have been attempted to be 
    /// coarsened but failed because of missing data from other processes
    inline List<uint> const & cells_to_request() const
    { return _cells_to_request; }

    /// Gives access to the list of vertices that have been attempted to be 
    /// deleted but failed because of missing data from other processes
    inline List<uint>& vertices_to_request()
    { return _vertices_to_request; }

    /// Gives access to the list of vertices that have been attempted to be 
    /// deleted but failed because of missing data from other processes
    inline List<uint> const & vertices_to_request() const
    { return _vertices_to_request; }

    /// Migrates cells according to request list
    bool migrate(Mesh& mesh, bool repeat);

#ifndef ____USE_D_MESH____
    /// Update distributed data
    void updateDistdata(Mesh& mesh, Mesh& coarse_mesh);
#endif

  private:
    /// Helper class that implements a BiMap-like datastructure, that allows lookup in both
    /// directions in O(1) and requires 2n of memory. Other than a map it also allows multiple
    /// entries with same key on the right hand side, e. g. two pairs like (0,-1) and (5,-1) can
    /// be contained at the same time. It is based on Array as datastructures and therefore 
    /// provides all comfort that is also available for STL-vectors.
    ///
    /// *Templates*
    ///   T1
    ///     type of the left hand side values of the pair
    ///   T2
    ///     type of the right hand side values of the pair
    template<
      typename T1, 
      typename T2
    > class BiMap {
    public:
      typedef std::pair<T1,T2> value_type;

    private:
      typedef Array<T1> left_type;
      typedef Array<T2> right_type;

    public:
      /// Creates a new BiMap with specified size
      ///
      /// *Arguments*
      ///   size (uint)
      ///     Size of the new BiMap
      explicit BiMap(uint size = 0) : _left( left_type(size) ), _right( right_type(size) ) {}

      /// Initializes the BiMap to a new size without keeping the existing entries, 
      /// similar to creating a new BiMap
      ///
      /// *Arguments*
      ///   size (uint)
      ///     Size of the BiMap
      inline void init(uint size) 
      { 
        _left = left_type(size);
        _right = right_type(size);
      }

      /// Resizes the Bimap to the specified size, keeping as much values as possible.
      ///
      /// *Arguments*
      ///   size (uint)
      ///     New size of the BiMap
      inline void resize(uint size)
      {
        _left.resize(size);
        _right.resize(size);
      }

      /// Sets a pair to the specified values
      ///
      /// *Arguments*
      ///   value (std::pair<T1,T2> const &)
      ///     The new values
      inline void set(value_type const & value) { set(value.first, value.second); }

      /// Sets a pair to the specified values
      ///
      /// *Arguments*
      ///   left (T1 const &)
      ///     New value for left hand side
      ///   right (T2 const &)
      ///     New value for right hand side
      inline void set(T1 const & left, T2 const & right)
      {
        if ( left >= 0 )
          _left[left] = right;

        if ( right >= 0 )
          _right[right] = left;
      }

      /// Gives acces to left hand side datastructure, e. g. for element access.
      /// A constant variant of the method exists as well
      ///
      /// *Returns*
      ///   Array<T1>&
      ///     A reference to the left hand side datastructure.
      inline left_type &        left()        { return _left;  }
      inline left_type  const & left()  const { return _left;  }

      /// Gives acces to right hand side datastructure, e. g. for element access.
      /// A constant variant of the method exists as well
      ///
      /// *Returns*
      ///   Array<T2>&
      ///     A reference to the right hand side datastructure.
      inline right_type &       right()       { return _right; }
      inline right_type const & right() const { return _right; }

    private:
      left_type  _left;
      right_type _right;
    };

  public:
    /// An IndexMap allows to keep track of the changing indices during the coarsening
    /// procedure. It enables mapping between the indices of the coarse and fine mesh in
    /// both directions using BiMap as datastructures. It stores pairs in the form of
    /// (fine_id,coarse_id).
    /// The mapping is kept twice to enable updating the map without affecting the existing
    /// map. This for example allows to build the new map between a new coarse mesh and the
    /// original fine map while keeping the map between the current coarse mesh and the 
    /// original fine map. The second map is accessed or modified by the functions with a 
    /// trailing ...New in the function name.
    class IndexMap
    {
    public:
      /// Instantiates a new IndexMap with the specified number of entries
      ///
      /// *Arguments*
      ///   num_entries (uint)
      ///     The desired number of entries (usually the number of cells or vertices of the
      ///     fine mesh)
      explicit IndexMap(uint num_entries = 0)
      : _fine2coarse( map_type(num_entries) )
      , _fine2coarse_new( map_type(num_entries) )
      { init(num_entries); }

      /// Initializes the IndexMap. If necessary the datastructures are resized. The map is
      /// filled initially with identity pairs of the form 
      /// (0,0), (1,1), (2,2), ..., (num_entries-1,num_entries-1)
      ///
      /// *Arguments*
      ///   num_entries (uint)
      ///     The desired number of entries
      inline void init(uint num_entries) {
        _fine2coarse.resize(num_entries);
        _fine2coarse_new.resize(num_entries);

        for ( uint i(0) ; i < num_entries ; ++i )
        {
          set(i,i);
          setNew(i,i);
        }
      }

      /// Gives the number of entries in the map from fine indices to coarse indices
      inline uint sizeFineFromCoarse() const 
      { return _fine2coarse.right().size(); }

      /// Gives the number of entries in the map from coarse indices to fine indices
      inline uint sizeCoarseFromFine() const 
      { return _fine2coarse.left().size(); }

      /// Gives the number of entries in the second map from fine indices to coarse indices
      inline uint sizeFineFromCoarseNew() const 
      { return _fine2coarse_new.right().size(); }

      /// Gives the number of entries in the second map from coarse indices to fine indices
      inline uint sizeCoarseFromFineNew() const 
      { return _fine2coarse_new.left().size(); }

      /// Gives the corresponding fine index to a given coarse index
      ///
      /// *Arguments*
      ///   coarse_id (int)
      ///     Index in the coarse mesh
      ///
      /// *Returns*
      ///   int
      ///     Index in the fine mesh, i.e. left hand entry in the pair (fine_id, coarse_id).
      inline int getFineFromCoarse(int coarse_id) const 
      { return _fine2coarse.right().at(coarse_id); }

      /// Gives the corresponding coarse index to a given fine index
      ///
      /// *Arguments*
      ///   fine_id (int)
      ///     Index in the fine mesh
      ///
      /// *Returns*
      ///   int
      ///     Index in the coarse mesh, i.e. right hand entry in the pair (fine_id, coarse_id).
      inline int getCoarseFromFine(int fine_id) const 
      { return _fine2coarse.left().at(fine_id); }

      /// Gives the corresponding fine index to a given coarse index in the second map
      ///
      /// *Arguments*
      ///   coarse_id (int)
      ///     Index in the coarse mesh
      ///
      /// *Returns*
      ///   int
      ///     Index in the fine mesh, i.e. left hand entry in the pair (fine_id, coarse_id).
      inline int getNewFineFromCoarse(int coarse_id) const 
      { return _fine2coarse_new.right().at(coarse_id); }
      
      /// Gives the corresponding coarse index to a given fine index in the second map
      ///
      /// *Arguments*
      ///   fine_id (int)
      ///     Index in the fine mesh
      ///
      /// *Returns*
      ///   int
      ///     Index in the coarse mesh, i.e. right hand entry in the pair (fine_id, coarse_id).
      inline int getNewCoarseFromFine(int fine_id) const 
      { return _fine2coarse_new.left().at(fine_id); }

      /// Gives the corresponding coarse index in the second map to a given coarse index in
      /// the first map.
      ///
      /// *Arguments*
      ///   coarse_id (int)
      ///     Index in the coarse mesh
      ///
      /// *Returns*
      ///   int
      ///     Index in the second coarse mesh
      inline int getNewCoarseFromCoarse(int coarse_id) const 
      { return getNewCoarseFromFine(getFineFromCoarse(coarse_id)); }

      /// Gives the corresponding coarse index in the first map to a given coarse index in
      /// the second map.
      ///
      /// *Arguments*
      ///   coarse_id (int)
      ///     Index in the second coarse mesh
      ///
      /// *Returns*
      ///   int
      ///     Index in the coarse mesh
      inline int getCoarseFromNewCoarse(int coarse_id) const 
      { return getCoarseFromFine(getNewFineFromCoarse(coarse_id)); }

      /// Sets the values in the pair (fine_id,coarse_id)
      ///
      /// *Arguments*
      ///   fine_id (int)
      ///     The index in the fine mesh
      ///   coarse_id (int)
      ///     The corresponding index in the coarse mesh
      inline void set(int fine_id, int coarse_id) 
      { _fine2coarse.set(fine_id,coarse_id); }

      /// Sets the values in the pair (fine_id,coarse_id) in the second map
      ///
      /// *Arguments*
      ///   fine_id (int)
      ///     The index in the fine mesh
      ///   coarse_id (int)
      ///     The corresponding index in the coarse mesh
      inline void setNew(int fine_id, int coarse_id) 
      { _fine2coarse_new.set(fine_id,coarse_id); }

      /// Changes the size of the mapping that gives fine indices to given coarse indices. This
      /// allows to save some memory when the size of the coarse mesh decreases.
      ///
      /// *Arguments*
      ///   num_entries (uint)
      ///     The new number of entries in the map.
      inline void setFineFromCoarseSize(uint num_entries) 
      { _fine2coarse.right().resize(num_entries); }

      /// Changes the size of the second mapping that gives fine indices to given coarse indices. 
      /// This allows to save some memory when the size of the coarse mesh decreases.
      ///
      /// *Arguments*
      ///   num_entries (uint)
      ///     The new number of entries in the map.
      inline void setNewFineFromCoarseSize(uint num_entries) 
      { _fine2coarse_new.right().resize(num_entries); }

      /// Commits the mapping of the second map to the first (main) map. After that both maps
      /// will carry the pairs of the second map.
      inline void commit() { _fine2coarse = _fine2coarse_new; }

      /// Reverts the mapping of the second map to the first (main) map. After that both maps
      /// will carry the pairs of the first (main) map.
      inline void revert() { _fine2coarse_new = _fine2coarse; }

    private:
      typedef BiMap<int,int> map_type;
      typedef map_type::value_type pair_type;

      map_type _fine2coarse;
      map_type _fine2coarse_new;
    };

  private:
#ifdef ____USE_D_MESH____
    /// Dynamic mesh
    DMesh* _dmesh;
#endif // ____USE_D_MESH____

    /// Indicator for vertices on domain boundaries
    Array<bool> _bnd_vertices;

    /// Indicator for cells on domain boundaries
    Array<bool> _bnd_cells;

    /// Indicator for vertices on process boundaries
    Array<bool> _int_bnd_vertices;

    /// Indicator for cells on process boundaries
    Array<bool> _int_bnd_cells;

    /// Indicator for independent set of vertices
    Array<bool> _forbidden_vertices;

#ifndef ____USE_D_MESH____
    /// IndexMaps for relation between coarse and fine mesh indices
    IndexMap _vertex_map;
    IndexMap _cell_map;
#endif

#ifdef ____USE_D_MESH____
    /// List of cells to coarsen
    List<DCell*> _cells_to_coarsen;
#else // ____USE_D_MESH____
    /// List of cells to coarsen
    List<uint> _cells_to_coarsen;
#endif // ____USE_D_MESH____
    /// List of cells that need neighboring cells from other processes
    List<uint> _cells_to_request;

    /// List of vertices that need neighboring cells from other processes
    List<uint> _vertices_to_request;

    /// Number of cells that have been migrated away in the last iteration
    uint _migrated_cells;

    /// performs part of the initialization which is also used by migrate().
    /// MeshFunction is templated, because migrate() needs it for uint while
    /// init() provides bool.
    template<typename T>
    void initCommon(MeshFunction<T>& cell_marker);

    /// Extract process boundary information and store them in 
    /// _int_bnd_vertices and _bnd_cells.
    void findInteriorBoundaries(Mesh& mesh);

    /// Extract global boundary information and store them in _bnd_vertices
    /// and _bnd_cells.
    void findDomainBoundaries(Mesh& mesh);

    /// Extract a independent set of the vertices by a greedy algorithm and
    /// store it in _forbidden_vertices
    void findIndependentSet(Mesh& mesh, bool coarsen_boundary);

    /// Check if the given vertex has neighbours that are part of the
    /// independent set and returns false in that case, otherwise true.
    bool isIndependentVertex(Vertex& v);

    /// Find all cells that are marked for coarsening and put them into a list
    template<typename T>
    void findCellsToCoarsen(MeshFunction<T>& cell_marker);
  };
}

#endif 
