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
  static SubSystemsManager& instance()
  {
    static SubSystemsManager instance_;
    return instance_;
  }

public:

  // Subsystem state mask
  enum Type { mpi       = 1,
              petsc     = 2,
              petscmpi  = 4,
              janpack   = 8,
              zoltan    = 16 };

  //-------------------------------------------------------------------------
  static int start(int argc = 0, char* argv[] = nullptr, uint n = 0, long w_limit = 0)
  {
    return SubSystemsManager::instance().init(argc, argv, n, w_limit);
  }

  //-------------------------------------------------------------------------
  static void status()
  {
    SubSystemsManager::instance().disp();
  }

  //-------------------------------------------------------------------------
  static inline bool active(SubSystemsManager::Type s)
  {
    return SubSystemsManager::instance().iset(s);
  }

  //-------------------------------------------------------------------------
  struct MPI
  {
    static SubSystemsManager::Type const flag = mpi;

    /// Initialize MPI
    static bool init(int argc = 0, char* argv[] = nullptr, uint n = 0);

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
    static bool init(int argc = 0, char* argv[] = nullptr);

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
    static bool init(int argc = 0, char* argv[] = nullptr);

    /// Finalize PETSc
    static bool fini();

    ///
    static int sema;
  };

  //-------------------------------------------------------------------------
  void disp() const;

  //--- ALARM ---------------------------------------------------------------

  static alarm& timer()
  {
    return SubSystemsManager::instance().alarm_handler();
  }

private:

  int init(int argc = 0, char* argv[] = nullptr, uint n = 0, long w_limit = 0);
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

  alarm& alarm_handler() { return timer_; }

  // State variable
  int count_;
  int state_;

  // Alarm handler
  alarm timer_;

};

} /* namespace dolfin */

#endif /* __DOLFIN_SUB_SYSTEMS_MANAGER_H */
