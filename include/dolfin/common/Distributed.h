// Copyright (C) 2017. Aurelien Larcher
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2017-05-10
//

#ifndef __DOLFIN_COMMON_DISTRIBUTED_H_
#define __DOLFIN_COMMON_DISTRIBUTED_H_

#include <dolfin/main/MPI.h>
#include <dolfin/log/log.h>

namespace dolfin
{

template<class T>
class Distributed
{

public:

  ///
  Distributed(MPI::Communicator& comm) :
    comm_(DOLFIN_COMM_NULL)
  {
#if HAVE_MPI
    /*
     * MPI 1.1:
     *  "A null handle argument is an erroneous IN argument in MPI calls"
     */
    if(comm != DOLFIN_COMM_NULL) MPI_Comm_dup(comm, &comm_);
#endif
  }

  ///
  Distributed(Distributed const& other) :
      comm_(DOLFIN_COMM_NULL)
  {
    *this = other;
  }

  /// Swap instances
  void swap(Distributed& other)
  {
    std::swap(comm_, other.comm_);
  }

  ///
  MPI::Communicator& comm()
  {
    return comm_;
  }

  ///
  uint comm_rank() const
  {
    int ret = 0;
#if HAVE_MPI
    if(comm_ != DOLFIN_COMM_NULL) MPI_Comm_rank(comm_, &ret);
#endif
    return uint(ret);
  }

  ///
  uint comm_size() const
  {
    int ret = 1;
#if HAVE_MPI
    if(comm_ != DOLFIN_COMM_NULL) MPI_Comm_size(comm_, &ret);
#endif
    return uint(ret);
  }

  ///
  inline bool distributed() const
  {
    return (this->comm_size() > 1);
  }

protected:

  ///
  virtual ~Distributed()
  {
#if HAVE_MPI
    if(comm_ != DOLFIN_COMM_NULL)  MPI_Comm_free(&comm_);
#endif
  }

  ///
  Distributed& operator=(Distributed const& other)
  {
    if (this != &other)
    {
#if HAVE_MPI
      if(comm_ != DOLFIN_COMM_NULL) MPI_Comm_free(&comm_);
      if(other.comm_ != DOLFIN_COMM_NULL) MPI_Comm_dup(other.comm_, &comm_);
#endif
    }
    return *this;
  }

private:

  MPI::Communicator comm_;

};

}

#endif /* __DOLFIN_COMMON_DISTRIBUTED_H_ */
