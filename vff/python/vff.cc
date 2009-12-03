//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif


#include <boost/python.hpp>
#include <boost/python/tuple.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_same.hpp>

#include "../functional.h"
#include "../layered.h"
#include "../va.h"

#include "vff.hpp"

namespace LaDa
{
  namespace Python
  {
    namespace XML
    {
      template<class T_TYPE>
        void Vff_from_XML(T_TYPE &_type, const std::string &_filename )
        {
          TiXmlDocument doc( _filename ); 
          TiXmlHandle docHandle( &doc ); 
        
          __DOASSERT( not doc.LoadFile(), 
                         doc.ErrorDesc() << "\n"  
                      << "Could not load input file " << _filename  
                      << ".\nAborting.\n" ) 
          __DOASSERT( not docHandle.FirstChild("Job").Element(),
                      "Could not find <Job> tag in " << _filename << ".\n" )
         
          __DOASSERT( not _type.Load( *docHandle.FirstChild("Job").Element() ),
                         "Could not load Vff functional from " + _filename + ".\n" )
        }
    }

    struct center_iterator
    {
      typdef vff::AtomicCenter::t_Centers::const_iterator type;

      type first_;
      type end_;
      bool is_first;

      center_iterator   (type const &_f, type const &_end)
                      : first_(_f), end_(_end), is_first(true) {}
      center_iterator   (center_iterator const &_c)
                      : first_(_c.first_), end_(_c.end_), is_first(_c.is_first_) {}

      center_iterator iter() const { return *this; }
      Value next() const 
      {
        if( first_ == end_ ) 
        {
          Py
        }
        if( is_first ) 
      };
    };

    template<class T_VFF>
      struct WithIterators : public T_VFF
      {
        WithIterator(Crystal::Structure &_str) : T_VFF(_str) {}
        WithIterator(WithIterators const &_c) : T_VFF(_c) {}
        virtual ~WithIterator() {}
        //! iterator over first neighbor tree.
        center_iterator iter() const { return center_iterator( centers_.begin(), centers_.end() ); }

        protected:
          using T_VFF::centers_;
          using T_VFF::operator();
          using T_VFF::print_escan_input;
          using T_VFF::init;
#         ifdef _MPI
            using T_VFF::set_mpi;
#         endif
      };


    //! Assumes ownership of the Crystal::Structure object needed by vff.
    template< class T_VFF > class Vff 
    {
      public:
        //! Type of the functional.
        typedef LaDa::Vff::VABase< WithIterators<T_VFF> > t_Functional;
        //! Constructor.
        Vff() { functional_.reset( new t_Functional( structure ) ); }
        //! Copy Constructor.
        Vff( const Vff& _c ) : structure( _c.structure )
          { functional_.reset( new t_Functional( structure ) ); }

        //! Redo initialization of vff from scratch.
        void init() { functional_->init( true ); }

        //! fast evaluation with no reinitialization.
        types::t_real operator()() const { return functional_->evaluate() / 16.0217733; }
        //! Returns the stress.
        atat::rMatrix3d get_stress() const { return functional_->get_stress(); }

        //! Loads from an XML input file.
        bool Load( const TiXmlElement &_node ) { return functional_->Load( _node ); }

        //! Prints escan input.
        void print_escan_input( const std::string &_path ) const
          { functional_->print_escan_input( _path ); }
          
        //! The owned structure.
        Crystal::Structure structure;
#       ifdef _MPI
          //! Sets mpi pointer.
          void set_mpi( boost::mpi::communicator* _c )
            { functional_->set_mpi( _c ); }
#       endif

        boost::python::tuple get_bond( const std::string& _bond ) const;
        void set_bond( const std::string& _bond, const boost::python::tuple& _t );
        boost::python::tuple get_angle( const std::string& _bond ) const;
        void set_angle( const std::string& _bond, const boost::python::tuple& _t );



      protected:
        //! The functional.
        boost::shared_ptr< LaDa::Vff::VABase< T_VFF > > functional_;
    };

    template< class T_VFF >
      boost::python::tuple Vff<T_VFF>::get_bond( const std::string &_str ) const 
      {
        try
        {
          typedef boost::tuples::tuple
          < 
            const types::t_real&, const types::t_real&,
            const types::t_real&, const types::t_real&,
            const types::t_real&, const types::t_real& 
          > t_Tuple;
          const t_Tuple result( functional_->get_bond( _str ) );
          return boost::python::make_tuple
          (
            boost::tuples::get<0>( result ), boost::tuples::get<1>( result ), 
            boost::tuples::get<2>( result ), boost::tuples::get<3>( result ), 
            boost::tuples::get<4>( result ), boost::tuples::get<5>( result )
          );
        }
        catch(...)
        {
          PyErr_SetString(PyExc_RuntimeError,
                          ("could not find parameters of bond " + _str).c_str() );
          boost::python::throw_error_already_set();
        }
      };
    template< class T_VFF >
      boost::python::tuple Vff<T_VFF>::get_angle( const std::string &_str ) const 
      {
        try
        {
          typedef boost::tuples::tuple
          < 
            const types::t_real&, const types::t_real&,
            const types::t_real&, const types::t_real&,
            const types::t_real&, const types::t_real&,
            const types::t_real&
          > t_Tuple;
          const t_Tuple result( functional_->get_angle( _str ) );
          return boost::python::make_tuple
          (
            boost::tuples::get<0>( result ), boost::tuples::get<1>( result ), 
            boost::tuples::get<2>( result ), boost::tuples::get<3>( result ), 
            boost::tuples::get<4>( result ), boost::tuples::get<5>( result ),
            boost::tuples::get<6>( result )
          );
        }
        catch(...)
        {
          PyErr_SetString(PyExc_RuntimeError,
                          ("could not find parameters of angle " + _str).c_str() );
          boost::python::throw_error_already_set();
        }
      };
    template< class T_VFF >
      void Vff<T_VFF>::set_bond( const std::string &_str, const boost::python::tuple &_t ) 
      {
        namespace bp = boost::python;
        const size_t N(bp::len( _t )); 
        if( N == 0 or bp::len( _t ) > 6 )
        {
          PyErr_SetString( PyExc_RuntimeError,
                           "Too many or to few bond parameters.\n"  );
          bp::throw_error_already_set();
          return;
        }
        try
        {
          const types::t_real
            a = N ? (types::t_real) (bp::extract<types::t_real>( _t[0] ) ): 0e0,
            b = N > 1 ? (types::t_real) (bp::extract<types::t_real>( _t[1] ) ): 0e0,
            c = N > 2 ? (types::t_real) (bp::extract<types::t_real>( _t[2] ) ): 0e0,
            d = N > 3 ? (types::t_real) (bp::extract<types::t_real>( _t[3] ) ): 0e0,
            e = N > 4 ? (types::t_real) (bp::extract<types::t_real>( _t[4] ) ): 0e0,
            f = N > 5 ? (types::t_real) (bp::extract<types::t_real>( _t[5] ) ): 0e0;
          functional_->set_bond( _str, boost::tuples::make_tuple( a, b, c, d, e, f ) );
        }
        catch(...)
        {
          PyErr_SetString(PyExc_RuntimeError,
                          ("could not find parameters of bond " + _str).c_str() );
          bp::throw_error_already_set();
        }
      };
    template< class T_VFF >
      void Vff<T_VFF>::set_angle( const std::string &_str, const boost::python::tuple &_t )
      {
        namespace bp = boost::python;
        const size_t N(bp::len( _t )); 
        if( N == 0 or bp::len( _t ) > 7 )
        {
          PyErr_SetString( PyExc_RuntimeError,
                           "Too many or to few bond parameters.\n"  );
          bp::throw_error_already_set();
          return;
        }
        try
        {
          const types::t_real 
            a = N ? (types::t_real) (bp::extract<types::t_real>( _t[0] ) ): 0e0,
            b = N > 1 ? (types::t_real) (bp::extract<types::t_real>( _t[1] ) ): 0e0,
            c = N > 2 ? (types::t_real) (bp::extract<types::t_real>( _t[2] ) ): 0e0,
            d = N > 3 ? (types::t_real) (bp::extract<types::t_real>( _t[3] ) ): 0e0,
            e = N > 4 ? (types::t_real) (bp::extract<types::t_real>( _t[4] ) ): 0e0,
            f = N > 5 ? (types::t_real) (bp::extract<types::t_real>( _t[5] ) ): 0e0,
            g = N > 6 ? (types::t_real) (bp::extract<types::t_real>( _t[6] ) ): 0e0;
          functional_->set_angle( _str, boost::tuples::make_tuple( a, b, c, d, e, f, g ) );
        }
        catch(...)
        {
          PyErr_SetString(PyExc_RuntimeError,
                          ("could not find parameters of bond " + _str).c_str() );
          bp::throw_error_already_set();
        }
      };

    class LVff : public LaDa::Vff::Layered
    {
      public:
        //! Constructor.
        LVff( LaDa::Crystal::Structure& _str ) : LaDa::Vff::Layered( _str ) {}
        //! Exposes direction setters.
        atat::rVector3d get_direction() const { return direction; } 
        //! Returns direction.
        void set_direction( const atat::rVector3d& _direction)
        {
          is_fixed_by_input = true;
          direction = _direction;
          create_template_strain();
        } 
    };

    class LayeredVff : public Vff< LVff >
    {
      public:
        //! Exposes direction setters.
        atat::rVector3d get_direction() const
          { return functional_->Vff().get_direction(); } 
        //! Returns direction.
        void set_direction( const atat::rVector3d& _direction)
          { functional_->Vff().set_direction( _direction ); } 
    };


#   ifdef EXPOSEVFF 
#     error Macro EXPOSEVFF already exists.
#   endif 
#   ifdef _MPI
#     ifdef SETMPI
#       error Macro SETMPI already exists.
#     endif 
#     define SETMPI(b) .def("set_mpi", &b::set_mpi, "Sets the boost.mpi communicator.")
#   else
#     define SETMPI(b)
#   endif
#   define EXPOSEVFF( a, b, c ) \
      bp::class_< b >( a, c ) \
        .def( bp::init< b >() ) \
        .def_readwrite( "structure",    &b::structure, \
                        "The structure to minimize. Input and Output to the functional." ) \
        .def( "fromXML",  &XML::Vff_from_XML<b>, bp::arg("file"),\
              "Loads the vff parameters from an XML file." ) \
        .def( "evaluate",  &b::operator(), \
              "Minimizes the current structure and returns the energy in eV." ) \
        .def( "init",  &b::init, "Initializes the functional for the current structure." ) \
        .def( "print_escan_input",  &b::print_escan_input, bp::arg("file"), \
              "Outputs the current structure in a format suitable for pescan." ) \
        .def( "set_bond",  &b::set_bond, ( bp::arg("bond"), bp::arg("params") ), \
              "Sets the parameters of bond from a tuple where:\n" \
              "  _ the first component is the bond length.\n" \
              "  _ the second to sixth components are the gammas at order 2-6.\n" \
              "The tuple can be shortened to include only the first few parameters.\n" ) \
        .def( "get_bond",  &b::get_bond, bp::arg("angle"), \
              "Returns the parameters of bond in form of a tuple.\n" \
              "  _ the first component is the bond length.\n" \
              "  _ the second to sixth components are the gammas at order 2-6.\n" ) \
        .def( "set_angle",  &b::set_angle, ( bp::arg("angle"), bp::arg("params") ), \
              "Sets the parameters of angle from a tuple where:\n" \
              "  _ the first component is gamma.\n" \
              "  _ the second component is sigma.\n" \
              "  _ the third to seventh components are the betas at order 2-6.\n" \
              "The tuple can be shortened to include only the first few parameters.\n" ) \
        .def( "get_angle",  &b::get_angle, bp::arg("angle"), \
              "Returns the parameters of angle in form of a tuple.\n" \
              "  _ the first component is gamma.\n" \
              "  _ the second component is sigma.\n" \
              "  _ the third to seventh components are the betas at order 2-6.\n" ) \
        .add_property( "stress",  &b::get_stress, \
                       "Returns the stress. Meaningfull only "\
                       "following a call to Vff.evaluate()." ) \
        SETMPI(b)

    void expose_vff()
    {
      typedef Vff< LaDa::Vff::Functional > t_Vff;
      namespace bp = boost::python;
      EXPOSEVFF
      ( 
        "Vff", 
        t_Vff, 
        "A Valence Force Field Functional.\n"
        "Prior to use, the parameters must be loaded from an XML file, "
        "The structure must be itself initialized, and LaDa.Vff.init() must be called."
        "Order does count :)."
      );
    }

    void expose_layeredvff()
    {
      typedef LayeredVff t_Vff;
      namespace bp = boost::python;
      EXPOSEVFF
      ( 
        "LayeredVff", 
        t_Vff,
        "A Valence Force Field Functional with epitaxial constraints.\n"
        "Prior to use, the parameters must be loaded from an XML file, "
        "The structure must be itself initialized, and LaDa.Vff.init() must be called."
        "Order does count :)."
      ).add_property( "direction",  &t_Vff::get_direction, &t_Vff::set_direction,
                      "Defines the direction of growth." );
    }

#   undef EXPOSEVFF
#   undef SETMPI

  }
} // namespace LaDa
