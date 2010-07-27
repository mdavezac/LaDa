#include "LaDaConfig.h"
#include "FCMangle.h"

#include <boost/python/def.hpp>

#include <physics/physics.h>
#include <crystal/structure.h>
#include <opt/mpi.h>

#include "escan.hpp"

//! \cond
extern "C"
{
  void FC_GLOBAL_(iaga_set_mpi, IAGA_SET_MPI)( MPI_Fint * );
  void FC_GLOBAL_(getvlarg, GETVLARG)();
  void FC_GLOBAL_(iaga_just_call_escan, IAGA_JUST_CALL_ESCAN)();
}
//! \endcond

void just_call_escan(boost::mpi::communicator const &_c)
{
  MPI_Comm __commC = (MPI_Comm) ( _c ) ;
  MPI_Fint __commF = MPI_Comm_c2f( __commC );
  FC_GLOBAL_(iaga_set_mpi, IAGA_SET_MPI)( &__commF );

  FC_GLOBAL_(iaga_just_call_escan, IAGA_just_CALL_ESCAN)();
}
// void just_call_escan2(MPI_Comm &_comm)
// {
//   MPI_Fint __commF = MPI_Comm_c2f( _comm );
//   FC_GLOBAL_(iaga_set_mpi, IAGA_SET_MPI)( &__commF );
//
//   FC_GLOBAL_(iaga_just_call_escan, IAGA_just_CALL_ESCAN)();
// }
void just_call_genpot(boost::mpi::communicator const &_c)
{
  MPI_Comm __commC = (MPI_Comm) ( _c ) ;
  MPI_Fint __commF = MPI_Comm_c2f( __commC );
  FC_GLOBAL_(iaga_set_mpi, IAGA_SET_MPI)( &__commF );

  FC_GLOBAL_(getvlarg, GETVLARG)();
}
// void just_call_genpot2(MPI_Comm &_comm)
// {
//   MPI_Fint __commF = MPI_Comm_c2f( _comm );
//   FC_GLOBAL_(iaga_set_mpi, IAGA_SET_MPI)( &__commF );
//
//   FC_GLOBAL_(getvlarg, GETVLARG)();
// }

namespace LaDa
{
  namespace python
  {
    namespace bp = boost::python;
    types::t_unsigned nb_valence_states( Crystal::TStructure<std::string> const &_str ) 
    {
      Crystal::TStructure<std::string>::t_Atoms::const_iterator i_atom = _str.atoms.begin();
      Crystal::TStructure<std::string>::t_Atoms::const_iterator i_atom_end = _str.atoms.end();
      types::t_unsigned bgstates = 0;
      for(; i_atom != i_atom_end; ++i_atom)
        bgstates += Physics::Atomic::Charge( i_atom->type );
      return bgstates;
    }
    void expose_escan()
    {
      bp::def("_call_escan", &just_call_escan, "Calls escan, accepts a boost.mpi.communicator.");
      bp::def("_call_genpot", &just_call_genpot, "Calls genpot, accepts a boost.mpi.communicator.");
      bp::def( "nb_valence_states", &nb_valence_states, bp::arg("structure"), 
               "Returns the number of valence states in a structure." );
    }
  } // namespace python
} // namespace LaDa
