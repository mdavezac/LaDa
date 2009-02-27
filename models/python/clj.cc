//
//  Version: $Id$
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/python/self.hpp> 
#include <boost/python/class.hpp> 
#include <boost/python/def.hpp> 
#include <boost/python/dict.hpp> 
#include <boost/python/suite/indexing/map_indexing_suite.hpp> 


namespace Python
{
  class Clj :: public Models :: Clj
  {
    public:
      typedef Clj :: LennardJones :: Bond t_Bond; 
      typedef Clj :: LennardJones :: t_Bonds t_Bonds; 
      typedef Clj :: Ewald :: t_Charges t_Charges; 
      typedef std::map< Key, types::t_int > t_ChargeMap; 

      const t_Bonds& get_bonds() const { return Clj :: LennardJones :: bonds_; }
      const t_Charges& get_charges() const { return Clj :: Ewald :: charges_; }
      void set_bonds( const t_Bonds& _bonds ) const
        { Clj :: LennardJones :: bonds_ = _bonds; }
      void set_charges( const t_Charges& _charges ) const
        { Clj :: Ewald :: charges_ = _charges; }

  };

  void expose_clj()
  {
    namespace bp = boost :: python;

    bp::class_< Clj :: t_Charges >( "Charges", "Dictionary of charges" )
      .def( bp :: map_indexing_suite< Clj :: t_Charges, true >() );
 
    bp::class_< Clj :: t_Bond >( "LJBond", "Holds bond parameters for Lennard-Jones." )
      .def_readwrite( "hardsphere", &Bond::hard_sphere, "+r^12 factor." )
      .def_readwrite( "vandderwalls", &Bond::van_der_walls, "-r^6 factor." );

    bp::class_< Clj :: t_Bonds >( "LJBonds", "Dictionary of bonds" )
      .def( bp :: map_indexing_suite< Clj :: t_Bonds, true >() );

    bp::class_< Clj >( "Clj", "Coulomb + LennardJones functional.\n" )
      .add_property( "charges", &Clj::get_charges, Clj::set_charges )
      .add_property( "bonds", &Clj::get_bonds, Clj::set_bonds )
      .def("__call__", &Clj::operator() )
      .def("gradient", &Clj::gradient );
  }
}
