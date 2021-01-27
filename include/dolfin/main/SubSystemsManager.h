// Copyright (C) 2008 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_SUB_SYSTEMS_MANAGER_H
#define __DOLFIN_SUB_SYSTEMS_MANAGER_H

#include <dolfin/common/types.h>
#include <dolfin/main/alarm.h>

namespace dolfin
{

/// This is a singleton class which manages the initialization and
/// finalization of various sub systems, such as MPI and PETSc.

class SubSystemsManager
{

  // Singleton instance
  static auto instance() -> SubSystemsManager &
  {
    static SubSystemsManager instance_;
    return instance_;
  }

public:
  // Subsystem state mask
  enum Type
  {
    mpi      = 1,
    petsc    = 2,
    petscmpi = 4,
    janpack  = 8,
    zoltan   = 16
  };

  //-------------------------------------------------------------------------
  static auto start( int    argc    = 0,
                     char * argv[]  = nullptr,
                     size_t n       = 0,
                     long   w_limit = 0 ) -> int
  {
    return SubSystemsManager::instance().init( argc, argv, n, w_limit );
  }

  //-------------------------------------------------------------------------
  static void status()
  {
    SubSystemsManager::instance().disp();
  }

  //-------------------------------------------------------------------------
  static inline auto active( SubSystemsManager::Type s ) -> bool
  {
    return SubSystemsManager::instance().iset( s );
  }

  //-------------------------------------------------------------------------
  struct MPI
  {
    static SubSystemsManager::Type const flag = mpi;

    /// Initialize MPI
    static auto init( int argc = 0, char * argv[] = nullptr, size_t n = 0 )
      -> bool;

    /// Finalize MPI
    static auto fini() -> bool;

    // Check if MPI has been initialized
    static auto initialized() -> bool;

    ///
    static int sema;
  };

  //-------------------------------------------------------------------------
  struct PETSc
  {
    static SubSystemsManager::Type const flag = petsc;

    /// Initialize PETSc with command-line arguments
    static auto init( int argc = 0, char * argv[] = nullptr ) -> bool;

    /// Finalize PETSc
    static auto fini() -> bool;

    ///
    static int sema;
  };

  //-------------------------------------------------------------------------
  struct PETScMPI
  {
    static SubSystemsManager::Type const flag = petscmpi;
  };

  //-------------------------------------------------------------------------
  struct JANPACK
  {
    static SubSystemsManager::Type const flag = janpack;
  };

  //-------------------------------------------------------------------------
  struct Zoltan
  {
    static SubSystemsManager::Type const flag = zoltan;

    /// Initialize PETSc with command-line arguments
    static auto init( int argc = 0, char * argv[] = nullptr ) -> bool;

    /// Finalize PETSc
    static auto fini() -> bool;

    ///
    static int sema;
  };

  //-------------------------------------------------------------------------
  void disp() const;

  //--- ALARM ---------------------------------------------------------------

  static auto timer() -> alarm &
  {
    return SubSystemsManager::instance().alarm_handler();
  }

private:
  auto init( int    argc    = 0,
             char * argv[]  = nullptr,
             size_t n       = 0,
             long   w_limit = 0 ) -> int;
  auto fini() -> int;

  // Constructor
  SubSystemsManager() = default;

  // Copy constructor
  SubSystemsManager( SubSystemsManager const & other );

  // Destructor
  ~SubSystemsManager();

  /// State control
  void init( SubSystemsManager::Type s )
  {
    state_ = ( state_ & s ) + ( state_ | s );
  }

  void fini( SubSystemsManager::Type s )
  {
    state_ &= ( 1 ^ s );
  }

  auto iset( SubSystemsManager::Type s ) const -> bool
  {
    return ( state_ & s ) == s;
  }

  auto alarm_handler() -> alarm &
  {
    return timer_;
  }

  // State variable
  int count_ { 0 };
  int state_ { 0 };

  // Alarm handler
  alarm timer_;
};

} /* namespace dolfin */

#endif /* __DOLFIN_SUB_SYSTEMS_MANAGER_H */
