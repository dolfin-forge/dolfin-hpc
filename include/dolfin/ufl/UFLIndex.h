// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  
// Last changed: 

#ifndef __UFL_INDEX_H_
#define __UFL_INDEX_H_

//#include <string>
//#include <vector>

#include <dolfin/ufl/UFLObject.h>

#include <dolfin/common/types.h>

namespace ufl
{

  /**
   *  DOCUMENTATION:
   *
   *  @class  Class
   *
   *  @brief  Provides an interface complying with UFL IndexSum.
   */

  class IndexBase : public Object
  {

    public:

      ///
      IndexBase(dolfin::uint const& count);

      ///
      ~IndexBase();

      //--- INTERFACE -------------------------------------------------------------
      
      /// 
      dolfin::uint const count() const;


      //--- INTERFACE inherited from UFLClass -------------------------------------
      
      /// __repr__
      repr_t const repr() const;

      /// __str__
      std::string const str() const;

      ///
      void display() const;

      ///
      IndexBase const* create(repr_t const & repr) const;

    private:

      dolfin::uint const count_;

      mutable repr_t repr_;
      mutable std::string str_;

  };

  class Index: public IndexBase
  {

    public:

      ///
      Index(dolfin::uint const& count);

      ///
      ~Index();

      //--- INTERFACE -------------------------------------------------------------

      //--- INTERFACE inherited from UFLClass -------------------------------------
      
      /// __repr__
      repr_t const repr() const;

      /// __str__
      std::string const str() const;

      ///
      void display() const;

      ///
      Index const* create(repr_t const & repr) const;

    private:

      mutable repr_t repr_;
      mutable std::string str_;

  };

  class FixedIndex: public IndexBase
  {

    public:

      ///
      FixedIndex(dolfin::uint const& value);

      ///
      ~FixedIndex();

      //--- INTERFACE -------------------------------------------------------------

      //--- INTERFACE inherited from UFLClass -------------------------------------
      
      /// __repr__
      repr_t const repr() const;

      /// __str__
      std::string const str() const;

      ///
      void display() const;

      ///
      FixedIndex const* create(repr_t const & repr) const;

    private:

      mutable repr_t repr_;
      mutable std::string str_;

  };

  class MultiIndex: public IndexBase
  {

    public:

      ///
      MultiIndex(dolfin::uint const& value);

      ///
      MultiIndex(IndexBase const& index);

      ///
      ~MultiIndex();

      //--- INTERFACE -------------------------------------------------------------

      //--- INTERFACE inherited from UFLClass -------------------------------------
      
      /// __repr__
      repr_t const repr() const;

      /// __str__
      std::string const str() const;

      ///
      void display() const;

      ///
      MultiIndex const* create(repr_t const & repr) const;

    private:

      mutable repr_t repr_;
      mutable std::string str_;

  };
} /* namespace ufl */
#endif /* __UFL_INDEX_H_ */
