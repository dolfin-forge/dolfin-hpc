// Copyright (C) 2008 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson 2009-2015.
// Modified by Aurelien Larcher 2017.
//
// First added:  2008-01-07
// Last changed: 2017-02-23

#ifndef __DOLFIN_SUB_SYSTEMS_MANAGER_H
#define __DOLFIN_SUB_SYSTEMS_MANAGER_H

#include <dolfin/common/types.h>

namespace dolfin
{

/// This is a singleton class which manages the initialization and
/// finalization of various sub systems, such as MPI and PETSc.

class SubSystemsManager
{

  // Singleton instance
  static SubSystemsManager sub_systems_manager;

public:

  // Subsystem state mask
  enum Type { mpi       = 1,
              petsc     = 2,
              petscmpi  = 4,
              janpack   = 8,
              zoltan    = 16,
              slepc     = 32 };

  //-------------------------------------------------------------------------
  static int start(int argc = 0, char* argv[] = NULL, uint n = 0)
  {
    return sub_systems_manager.init(argc, argv, n);
  }

  //-------------------------------------------------------------------------
  static void status()
  {
    sub_systems_manager.disp();
  }

  //-------------------------------------------------------------------------
  static inline bool active(SubSystemsManager::Type s)
  {
    return sub_systems_manager.iset(s);
  }

  //-------------------------------------------------------------------------
  struct MPI
  {
    static SubSystemsManager::Type const flag = mpi;

    /// Initialize MPI
    static bool init(int argc = 0, char* argv[] = NULL, uint n = 0);

    /// Finalize MPI
    static bool fini();

    // Check if MPI has been initialized
    static bool initialized();

    ///
    static int sema;
  };

  //-------------------------------------------------------------------------
  struct PETSc
  {
    static SubSystemsManager::Type const flag = petsc;

    /// Initialize PETSc with command-line arguments
    static bool init(int argc = 0, char* argv[] = NULL);

    /// Finalize PETSc
    static bool fini();

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
    static bool init(int argc = 0, char* argv[] = NULL);

    /// Finalize PETSc
    static bool fini();

    ///
    static int sema;
  };

  //-------------------------------------------------------------------------
  struct SLEPc
  {
    static SubSystemsManager::Type const flag = slepc;
  };

  //-------------------------------------------------------------------------
  void disp() const;

private:

  int init(int argc = 0, char* argv[] = NULL, uint n = 0);
  int fini();

  // Constructor
  SubSystemsManager();

  // Copy constructor
  SubSystemsManager(SubSystemsManager const& other);

  // Destructor
  ~SubSystemsManager();

  /// State control
  void init(SubSystemsManager::Type s) { state_ = (state_ & s) + (state_ | s); }

  void fini(SubSystemsManager::Type s) { state_ &= (1 ^ s); }

  bool iset(SubSystemsManager::Type s) const { return (state_ & s) == s; }

  // State variable
  int count_;
  int state_;

};

} /* namespace dolfin */

#endif /* __DOLFIN_SUB_SYSTEMS_MANAGER_H */
