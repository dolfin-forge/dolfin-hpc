// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  
// Last changed: 

#ifndef __UFL_INDEX_H_
#define __UFL_INDEX_H_

//#include <string>
//#include <vector>

#include <dolfin/ufl/UFLClass.h>

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

  class IndexBase : public Class
  {

    public:

//      ///
//      IndexBase(dolfin::uint const& count);

      //--- INTERFACE -------------------------------------------------------------
      
      /// 
      dolfin::uint const & count() const;


      //--- INTERFACE inherited from UFLClass -------------------------------------
      
      /// __repr__
      virtual repr_t const repr() const = 0;

      /// __str__
      virtual std::string const str() const = 0;

      ///
      virtual void display() const;

    protected:

      ///
      IndexBase(std::string const& name,
          dolfin::uint const& count);

      ///
      IndexBase(std::string const& name,
          IndexBase const& index);

      ///
      IndexBase(std::string const & name, repr_t const & repr);

      ///
      virtual ~IndexBase();


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
      Index(repr_t const & repr);

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
      FixedIndex(repr_t const & repr);

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
      MultiIndex(repr_t const & repr);

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

    private:

      mutable repr_t repr_;
      mutable std::string str_;

  };
} /* namespace ufl */
#endif /* __UFL_INDEX_H_ */
