#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/python/class.hpp>
#include <boost/python/object.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/def.hpp>
#include <boost/python/str.hpp>
#include <boost/python/enum.hpp>
#include <boost/python/data_members.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/with_custodian_and_ward.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/make_constructor.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>


#include <python/numpy_types.h>
#include <opt/tuple_serialize.h>

#include "escan.hpp"

//! \cond
extern "C"
{
  void FC_FUNC_(iaga_set_mpi, IAGA_SET_MPI)( MPI_Fint * );
  void FC_FUNC_(getvlarg, GETVLARG)();
  void FC_FUNC_(iaga_just_call_escan, IAGA_JUST_CALL_ESCAN)();
}
//! \endcond

void just_call_escan(boost::mpi::communicator const &_c)
{
  MPI_Comm __commC = (MPI_Comm) ( _c ) ;
  MPI_Fint __commF = MPI_Comm_c2f( __commC );
  FC_FUNC_(iaga_set_mpi, IAGA_SET_MPI)( &__commF );

  FC_FUNC_(iaga_just_call_escan, IAGA_just_CALL_ESCAN)();
}
// void just_call_escan2(MPI_Comm &_comm)
// {
//   MPI_Fint __commF = MPI_Comm_c2f( _comm );
//   FC_FUNC_(iaga_set_mpi, IAGA_SET_MPI)( &__commF );
//
//   FC_FUNC_(iaga_just_call_escan, IAGA_just_CALL_ESCAN)();
// }
void just_call_genpot(boost::mpi::communicator const &_c)
{
  MPI_Comm __commC = (MPI_Comm) ( _c ) ;
  MPI_Fint __commF = MPI_Comm_c2f( __commC );
  FC_FUNC_(iaga_set_mpi, IAGA_SET_MPI)( &__commF );

  FC_FUNC_(getvlarg, GETVLARG)();
}
// void just_call_genpot2(MPI_Comm &_comm)
// {
//   MPI_Fint __commF = MPI_Comm_c2f( _comm );
//   FC_FUNC_(iaga_set_mpi, IAGA_SET_MPI)( &__commF );
//
//   FC_FUNC_(getvlarg, GETVLARG)();
// }

namespace LaDa
{
  namespace bp = boost::python;
    void expose_escan()
    {
      bp::def("_call_escan", &just_call_escan);
//     bp::def("_call_escan", &just_call_escan2, "Private interface. @see lada.escan.call_escan.");
      bp::def("_call_genpot", &just_call_genpot);
//     bp::def("_call_genpot", &just_call_genpot2, "Private interface. @see lada.escan.call_escan.");

      bp::def( "nb_valence_states", &nb_valence_states<Crystal::TStructure<std::string> >,
               bp::arg("structure"), "Returns the number of valence states in a structure." );
    }

  } // namespace python
} // namespace LaDa
