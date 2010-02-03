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
              (   "Could not parse xml file " + _filename + ".\n" 
                + doc.ErrorDesc() + "\n" ).c_str()
            );
            bp::throw_error_already_set();
          }
          else if( not docHandle.FirstChild("Job").Element() )
          {
            PyErr_SetString
            (
              PyExc_RuntimeError,
              ("Could not find <Job> in xml file " + _filename + ".\n").c_str()
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
                  ("Could not load VFF functional from xml file " + _filename + ".\n").c_str()
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
                     "Encountered error while loading  VFF functional from xml file "
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
        { return cload<T>(_filename); }
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
        void init(bool _redo, bool _verbose) { functional_->init( _redo, _verbose ); }

        //! fast evaluation with no reinitialization.
        boost::shared_ptr< Crystal::TStructure<std::string> >
          operator()(Crystal::TStructure<std::string> const& _str, bool _doinit) 
          { 
            Crystal::convert_string_to_real_structure(_str, structure);     
            functional_->init(_doinit, false);
            types::t_real const energy = functional_->evaluate() / 16.0217733;

            boost::shared_ptr< Crystal::TStructure<std::string> >
              result( new Crystal::TStructure<std::string> );
            Crystal::convert_real_to_string_structure(structure, *result);     
            result->energy = energy;
            return result;
          }
        //! Returns the stress.
        math::rMatrix3d get_stress() const { return functional_->get_stress(); }

        //! Loads from an XML input file.
        bool Load( const TiXmlElement &_node ) { return functional_->Load( _node ); }

        //! Prints escan input.
        void print_escan_input( const std::string &_path,
                                Crystal::TStructure<std::string> const &_str) 
        {
          Crystal::convert_string_to_real_structure(_str, structure);     
          functional_->print_escan_input( _path ); 
        }
          
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


    template<class T> bp::class_<T> expose_vff_(std::string const &_name, std::string const &_doc)
    {
      return bp::class_< T >
      ( 
        _name.c_str(),
        (
          _doc
          + "\n\nThis object can be created with: \n"
            "  - No argument.\n"
            "  - A single string argument representing the path to an XML input file.\n" 
            "  - A single boost.mpi communicator.\n"
            "  - A string argument (see above), followed by a boost.mpi communicator.\n\n" 
            "If compiled without mpi, including the communicator will have no effect.\n"
        ).c_str()
      ).def( bp::init< T const &>() ) 
       .def( "__init__", bp::make_constructor(&cload<T>) )
       .def( "__init__", bp::make_constructor(&mpi_create<T>) )
       .def( "__init__", bp::make_constructor(&mpi_create2<T>) )
       .def( "fromXML",  &XML::Vff_from_XML<T>, bp::arg("file"),
             "Loads the vff parameters from an XML file." ) 
       .def
       ( 
         "__call__",  
         &T::operator(), 
         (bp::arg("structure"), bp::arg("doinit") = true),
         "Minimizes structure.\n\n"
         "@param structure: structure to evaluate.\n"
         "@type structure: L{crystal.Structure}.\n"
         "@param doinit: If true, reconstructs first-neighbor tree. Default: true.\n"
         "@return: the relaxed structure. The energy is in "
         "structure.L{energy<crystal.Structure.energy>}.\n"
       )
       .def( "_init",  &T::init, (bp::arg("redo_tree") = true, bp::arg("verbose")=false), 
             "Initializes the functional for the current structure." ) 
       .def( "print_escan_input",  &T::print_escan_input, (bp::arg("file"), bp::arg("structure")), 
             "Outputs the current structure in a format suitable for pescan." ) 
       .def( "set_bond",  &T::set_bond, ( bp::arg("bond"), bp::arg("params") ), 
             "Sets the parameters of bond from a tuple where:\n" 
             "  - the first component is the bond length.\n" 
             "  - the second to sixth components are the gammas at order 2-6.\n" 
             "The tuple can be shortened to include only the first few parameters.\n" ) 
       .def( "get_bond",  &T::get_bond, bp::arg("angle"), 
             "Returns the parameters of bond in form of a tuple.\n" 
             "  - the first component is the bond length.\n" 
             "  - the second to sixth components are the gammas at order 2-6.\n" ) 
       .def( "set_angle",  &T::set_angle, ( bp::arg("angle"), bp::arg("params") ), 
             "Sets the parameters of angle from a tuple where:\n" 
             "  - the first component is gamma.\n" 
             "  - the second component is sigma.\n" 
             "  - the third to seventh components are the betas at order 2-6.\n" 
             "The tuple can be shortened to include only the first few parameters.\n" ) 
       .def( "get_angle",  &T::get_angle, bp::arg("angle"), 
             "Returns the parameters of angle in form of a tuple.\n" 
             "  - the first component is gamma.\n" 
             "  - the second component is sigma.\n" 
             "  - the third to seventh components are the betas at order 2-6.\n" ) 
#      ifdef _MPI
         .def("set_mpi", &T::set_mpi, "Sets the boost.mpi communicator.")
#      endif
       .add_property( "stress",  &T::get_stress, 
                      "Returns the stress. Meaningfull only "
                      "following a call to Vff.evaluate()." );
    }

    void expose_vff()
    {
      expose_vff_< Vff<LaDa::Vff::Functional> > 
      ( 
        "Vff", 
        "A Valence Force Field Functional.\n\n"
        "Usage:\n"
        ">>> vff = lada.vff.Vff(\"input.xml\", boost.mpi.world)\n"
        ">>> structure = crystal.Structure(\"input.xml\")\n"
        ">>> relaxed_structure = vff(structure)\n\n"
      );
    }

    void expose_layeredvff()
    {
      expose_vff_<LayeredVff>
      ( 
        "LayeredVff", 
        "A Valence Force Field Functional with epitaxial constraints.\n\n"
        "Usage:\n"
        ">>> vff = lada.vff.LayeredVff(\"input.xml\", boost.mpi.world)\n"
        ">>> vff.direction = numpy.array([0,0,1], dtype=\"float64\")"
        ">>> structure = crystal.Structure(\"input.xml\")\n"
        ">>> relaxed_structure = vff(structure)\n\n"
      ).add_property
       (
         "direction",
         bp::make_function
         (
           &LayeredVff::get_direction, 
           bp::return_value_policy<bp::return_by_value>()
         ),
         &LayeredVff::set_direction, "Growth/Epitaxial direction.\n\n3x1 float64 numpy array.\n" 
       ); 
    }

#   undef EXPOSEVFF
#   undef SETMPI

  }
} // namespace LaDa
