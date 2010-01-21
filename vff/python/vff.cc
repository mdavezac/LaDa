//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif


#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/make_constructor.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_same.hpp>

#include "../functional.h"
#include "../layered.h"
#include "../va.h"

#include "vff.hpp"

namespace LaDa
{
  namespace bp = boost::python;
  namespace Python
  {
    namespace XML
    {
      template<class T_TYPE>
        void Vff_from_XML(T_TYPE &_type, const std::string &_filename )
        {
          TiXmlDocument doc( _filename ); 
          TiXmlHandle docHandle( &doc ); 
        
          if( not doc.LoadFile() )
          {
            PyErr_SetString
            (
              PyExc_RuntimeError,
              (   "Could not parse xml file" + _filename + ".\n" 
                + doc.ErrorDesc() + "\n" ).c_str()
            );
            bp::throw_error_already_set();
          }
          else if( not docHandle.FirstChild("Job").Element() )
          {
            PyErr_SetString
            (
              PyExc_RuntimeError,
              ("Could not find <Job> in xml file" + _filename + ".\n").c_str()
            );
            bp::throw_error_already_set();
          }
          else
          {
            try
            {
              if(not _type.Load( *docHandle.FirstChild("Job").Element() ) )
              {
                PyErr_SetString
                (
                  PyExc_RuntimeError,
                  ("Could not load VFF functional from xml file" + _filename + ".\n").c_str()
                );
                bp::throw_error_already_set();
              }
            }
            catch(std::exception &_e)
            {
              PyErr_SetString
              (
                PyExc_RuntimeError,
                (
                     "Encountered error while loading  VFF functional from xml file"
                   + _filename + ".\n"
                   + _e.what() + "\n"
                ).c_str()
              );
              bp::throw_error_already_set();
            }
          }
        }
    }

//     struct center_iterator
//     {
//       typedef Vff::AtomicCenter::t_Centers::const_iterator type;
//
//       type first_;
//       type end_;
//       bool is_first;
//
//       center_iterator   (type const &_f, type const &_end)
//                       : first_(_f), end_(_end), is_first(true) {}
//       center_iterator   (center_iterator const &_c)
//                       : first_(_c.first_), end_(_c.end_), is_first(_c.is_first_) {}
//
//       center_iterator iter() const { return *this; }
//       Value next() const 
//       {
//         if( first_ == end_ ) 
//         {
//           Py
//         }
//         if( is_first ) 
//       };
//     };
//
//     template<class T_VFF>
//       struct WithIterators : public T_VFF
//       {
//         WithIterator(Crystal::Structure &_str) : T_VFF(_str) {}
//         WithIterator(WithIterators const &_c) : T_VFF(_c) {}
//         virtual ~WithIterator() {}
//         //! iterator over first neighbor tree.
//         center_iterator iter() const { return center_iterator( centers_.begin(), centers_.end() ); }
//
//         protected:
//           using T_VFF::centers_;
//           using T_VFF::operator();
//           using T_VFF::print_escan_input;
//           using T_VFF::init;
// #         ifdef _MPI
//             using T_VFF::set_mpi;
// #         endif
//       };

    template<class T> T* cload(std::string const &_filename) 
    {
      T* result(new T);
      if( not result )
      {
        PyErr_SetString(PyExc_RuntimeError, "Memory error.\n");
        bp::throw_error_already_set();
        return NULL;
      }
      XML::Vff_from_XML(*result, _filename); 
      if( PyErr_Occurred() != NULL ) 
      {
        delete result;
        return NULL;
      }
      return result;
    }
#   ifdef _MPI
      template<class T> T* mpi_create(boost::mpi::communicator *_c)
      {
        T* result(new T);
        if( not result )
        {
          PyErr_SetString(PyExc_RuntimeError, "Memory error.\n");
          bp::throw_error_already_set();
          return NULL;
        }
        result->set_mpi(_c);
        return result;
      }
      template<class T> T* mpi_create2( std::string const &_filename,
                                        boost::mpi::communicator *_c )
                                       
      {
        T* result(new T);
        if( not result )
        {
          PyErr_SetString(PyExc_RuntimeError, "Memory error.\n");
          bp::throw_error_already_set();
          return NULL;
        }
        result->set_mpi(_c);
        XML::Vff_from_XML(*result, _filename); 
        if( PyErr_Occurred() != NULL ) 
        {
          delete result;
          return NULL;
        }
        return result;
      }
#   else 
      template<class T> T* mpi_create(bp::object const &)
      {
        T* result(new T);
        if( not result )
        {
          PyErr_SetString(PyExc_RuntimeError, "Memory error.\n");
          bp::throw_error_already_set();
          return NULL;
        }
        return result;
      }
      template<class T> T* mpi_create2(bp::object const &, std::string const &_filename)
        { return cload(_filename); }
#   endif

    //! Assumes ownership of the Crystal::Structure object needed by vff.
    template< class T_VFF > class Vff 
    {
      public:
        //! Type of the functional.
        typedef LaDa::Vff::VABase< T_VFF > t_Functional;
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
        math::rMatrix3d get_stress() const { return functional_->get_stress(); }

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

        bp::tuple get_bond( const std::string& _bond ) const;
        void set_bond( const std::string& _bond, const bp::tuple& _t );
        bp::tuple get_angle( const std::string& _bond ) const;
        void set_angle( const std::string& _bond, const bp::tuple& _t );



      protected:
        //! The functional.
        boost::shared_ptr< LaDa::Vff::VABase< T_VFF > > functional_;
    };

    template< class T_VFF >
      bp::tuple Vff<T_VFF>::get_bond( const std::string &_str ) const 
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
          return bp::make_tuple
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
          bp::throw_error_already_set();
        }
      };
    template< class T_VFF >
      bp::tuple Vff<T_VFF>::get_angle( const std::string &_str ) const 
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
          return bp::make_tuple
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
          bp::throw_error_already_set();
        }
      };
    template< class T_VFF >
      void Vff<T_VFF>::set_bond( const std::string &_str, const bp::tuple &_t ) 
      {
        namespace bp = bp;
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
      void Vff<T_VFF>::set_angle( const std::string &_str, const bp::tuple &_t )
      {
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
        math::rVector3d get_direction() const { return direction; } 
        //! Returns direction.
        void set_direction( const math::rVector3d& _direction)
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
        math::rVector3d get_direction() const
          { return functional_->Vff().get_direction(); } 
        //! Returns direction.
        void set_direction( const math::rVector3d& _direction)
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
      bp::class_< b >\
      ( \
        a,\
        (std::string(c) + "\n\nThis object can be created with: \n"\
            "  - No argument.\n"\
            "  - A single string argument representing the path to an XML input file.\n" \
            "  - A single boost.mpi communicator.\n"\
            "  - A string argument (see above), followed by a boost.mpi communicator.\n\n" \
            "If compiled without mpi, including the communicator will have no effect.\n").c_str()\
       ) \
        .def( bp::init< b const &>() ) \
        .def( "__init__", bp::make_constructor(&cload<b>) )\
        .def( "__init__", bp::make_constructor(&mpi_create<b>) )\
        .def( "__init__", bp::make_constructor(&mpi_create2<b>) )\
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
              "  - the first component is the bond length.\n" \
              "  - the second to sixth components are the gammas at order 2-6.\n" \
              "The tuple can be shortened to include only the first few parameters.\n" ) \
        .def( "get_bond",  &b::get_bond, bp::arg("angle"), \
              "Returns the parameters of bond in form of a tuple.\n" \
              "  - the first component is the bond length.\n" \
              "  - the second to sixth components are the gammas at order 2-6.\n" ) \
        .def( "set_angle",  &b::set_angle, ( bp::arg("angle"), bp::arg("params") ), \
              "Sets the parameters of angle from a tuple where:\n" \
              "  - the first component is gamma.\n" \
              "  - the second component is sigma.\n" \
              "  - the third to seventh components are the betas at order 2-6.\n" \
              "The tuple can be shortened to include only the first few parameters.\n" ) \
        .def( "get_angle",  &b::get_angle, bp::arg("angle"), \
              "Returns the parameters of angle in form of a tuple.\n" \
              "  - the first component is gamma.\n" \
              "  - the second component is sigma.\n" \
              "  - the third to seventh components are the betas at order 2-6.\n" ) \
        .add_property( "stress",  &b::get_stress, \
                       "Returns the stress. Meaningfull only "\
                       "following a call to Vff.evaluate()." ) \
        SETMPI(b)

    void expose_vff()
    {
      typedef Vff< LaDa::Vff::Functional > t_Vff;
      EXPOSEVFF
      ( 
        "Vff", 
        t_Vff, 
        "A Valence Force Field Functional.\n\n"
        "Prior to use, the parameters must be loaded from an XML file, "
        "The structure must be itself initialized, and LaDa.Vff.init() must be called."
        "Order does count :).\n\n"
        ">>> vff = lada.vff.Vff(\"input.xml\", boost.mpi.world)\n"
        ">>> vff.structure.fromXML(\"input.xml\")\n"
        ">>> vff.init()\n"
        ">>> vff.evaluate()\n\n"
        "For deep copy, one may use the default constructor: >>> vff2 = lada.vff.Vff(vff1)\n"

      );
    }

    void expose_layeredvff()
    {
      typedef LayeredVff t_Vff;
      EXPOSEVFF
      ( 
        "LayeredVff", 
        t_Vff,
        "A Valence Force Field Functional with epitaxial constraints.\n\n"
        "Prior to use, the parameters must be loaded from an XML file, "
        "The structure must be itself initialized, and LaDa.Vff.init() must be called."
        "Order does count :).\n\n"
        ">>> vff = lada.vff.Vff()\n"
        ">>> vff.set_mpi(boost.mpi.world) # if compiled with mpi only!\n"
        ">>> vff.fromXML(\"input.xml\")\n"
        ">>> vff.structure.fromXML(\"input.xml\")\n"
        ">>> # if next line commented out, uses vff.structure[:,0] by default.\n" 
        ">>> vff.direction = numpy.array([0,0,1], dtype=\"float64\")"
        ">>> vff.init()\n"
        ">>> vff.evaluate()\n\n"
        "For deep copy, one may use the default constructor: C{vff2 = lada.vff.Vff(vff1)}\n"
      ).add_property
       (
         "direction",
         bp::make_function(&t_Vff::get_direction, bp::return_value_policy<bp::return_by_value>()),
         &t_Vff::set_direction, "Growth/Epitaxial direction.\n\n3x1 float64 numpy array.\n" 
       ); 
    }

#   undef EXPOSEVFF
#   undef SETMPI

  }
} // namespace LaDa
