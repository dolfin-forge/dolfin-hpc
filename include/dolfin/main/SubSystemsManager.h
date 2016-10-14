// Copyright (C) 2008 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson 2009-2015.
//
// First added:  2008-01-07
// Last changed: 2015-01-30

#ifndef __DOLFIN_SUB_SYSTEMS_MANAGER_H
#define __DOLFIN_SUB_SYSTEMS_MANAGER_H

namespace dolfin
{

  /// This is a singleton class which manages the initialisation and
  /// finalisation of various sub systems, such as MPI and PETSc.

  class SubSystemsManager
  {
  public:

    /// Initialise MPI and return if initialization
    static bool initMPI(int argc = 0, char* argv[] = 0, uint n = 0);

    /// Initialize PETSc without command-line arguments
    static bool initPETSc();

    /// Initialize PETSc with command-line arguments
    static bool initPETSc(int argc, char* argv[], bool cmd_line_args = true);

    /// Initialize Zoltan without command-line arguments
    static bool initZoltan();

    /// Initialize Zoltan without command-line arguments
    static bool initZoltan(int argc, char* argv[]);

    /// Finalize MPI
    static void finalizeMPI();

    /// Finalize PETSc
    static void finalizePETSc();

    // Check if MPI has been initialised (returns true if MPI has been
    //   initialised, even if it is later finalised)
    static bool MPIinitialized();

  private:

    // Constructor
    SubSystemsManager();

    // Copy construtor
    SubSystemsManager(SubSystemsManager const& sub_sys_manager);

    // Destructor
    ~SubSystemsManager();

    // Singleton instance
    static SubSystemsManager sub_systems_manager;

    // Static state variables
    static uint mpi_init_sema_;

    // State variables
    bool petsc_initialized;
    bool petsc_controls_mpi;
    bool zoltan_initialized;

  };

}

#endif
