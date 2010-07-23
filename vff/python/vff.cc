#include "LaDaConfig.h"


#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/make_constructor.hpp>
#include <boost/python/data_members.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/python/register_ptr_to_python.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#ifdef LADA_MPI
# include <boost/mpi/collectives.hpp>
#endif


#include "vff.hpp"

namespace LaDa
{
  namespace bp = boost::python;
  namespace python
  {
    template<class T> 
      boost::shared_ptr<T> create()
      {
        typename T::second_type second;
        return boost::shared_ptr<T>(new T( typename T::first_type(second), second )); 
      }
    template<class T> 
#     ifndef _MPI
        boost::shared_ptr<T> create_mpi(bp::object const& _object)
        { 
          try
          {
            std::string const filename = bp::extract<std::string>(_object);
            return create_input(filename);
          }
          catch(...){}
          return create<T>(); 
        }
#     else
        boost::shared_ptr<T> create_mpi(boost::mpi::communicator *_c )
        {
          boost::shared_ptr<T> result = create<T>();
          result->first.set_mpi(_c);
          return result;
        }
#     endif
    template<class T> 
      boost::shared_ptr<T> create_docnhand(TiXmlDocument &_doc, TiXmlHandle &_docHandle)
      {
        boost::shared_ptr<T> result = create<T>();
        typename T::first_type &_type = result->first;
      
        if( not _docHandle.FirstChild("Job").Element() )
        {
          PyErr_SetString(PyExc_RuntimeError, "Could not find <Job> in xml file.\n");
          bp::throw_error_already_set();
        }
        else
        {
          try
          {
            if(not _type.Load( *_docHandle.FirstChild("Job").Element() ) )
            {
              PyErr_SetString(PyExc_RuntimeError, "Could not load VFF functional.");
              bp::throw_error_already_set();
            }
          }
          catch(std::exception &_e)
          {
            PyErr_SetString(PyExc_RuntimeError, "Could not create VFF functional.");
            bp::throw_error_already_set();
          }
        }
        return result;
      }
    template<class T> 
      boost::shared_ptr<T> create_input(std::string const &_filename ) 
      {
        TiXmlDocument doc( _filename ); 
        TiXmlHandle docHandle( &doc ); 
        if( not doc.LoadFile() )
        {
          std::ostringstream sstr;
          sstr << "Could not parse xml file.\n" << doc.ErrorDesc() << "\n";
          PyErr_SetString(PyExc_RuntimeError, sstr.str().c_str());
          bp::throw_error_already_set();
        }
        return create_docnhand<T>(doc, docHandle);
      }
#     ifndef _MPI
      template<class T>
        boost::shared_ptr<T> create_inputmpi(std::string const& _f, bp::object const&) 
          { return create_input(_f); }
#     else
      template<class T>
        boost::shared_ptr<T> create_inputmpi(std::string const& _f, boost::mpi::communicator *_c )
        {
          namespace bfs = boost::filesystem;
          namespace bm = boost::mpi;
          std::string error = "", filestring;
          if( _c->rank() == 0 )
          {
            if( not bfs::exists(_f) ) error = _f + " does not exist.";
            else try
            {
              std::ifstream file( _f.c_str(), std::ios_base::in );
              std::istream_iterator< std::string > i_file( file );
              std::istream_iterator< std::string > i_file_end;
              for(; i_file != i_file_end; ++i_file )
              {
                filestring.append( *i_file );
                filestring.append(" ");  
              }
              file.close();
            }
            catch (std::exception const &_e)
            {
              std::ostringstream sstr;
              sstr << "Encountered error while creating Vff C++ functional.\n"
                   << _e.what();
              error =  sstr.str();
            }
          }
          bm::broadcast(*_c, error, 0);
          if( error.size() != 0 )
          {
            PyErr_SetString(PyExc_RuntimeError, error.c_str());
            return boost::shared_ptr<T>();
          }
          bm::broadcast(*_c, filestring, 0);
          TiXmlDocument doc;
          doc.Parse(filestring.c_str()); 
          TiXmlHandle handle(&doc); 
          boost::shared_ptr<T> result = create_docnhand<T>(doc, handle); 
          result->first.set_mpi(_c);
          return result;
        }
#     endif

    template<class T>  
        boost::shared_ptr<T> create_copy(T const &_f) { return boost::shared_ptr<T>( new T(_f) ); }
    template<class T> 
      void init( T &_self, bool _redo, bool _verbose ) { _self.first.init(_redo, _verbose); }
    template<class T> 
      bp::tuple __call__( T &_self, Crystal::TStructure<std::string> const &_str,
                          bool _doinit, bool relax )
      { 
        Crystal::convert_string_to_real_structure(_str, _self.second);     
        init(_self, _doinit, false);
        types::t_real const energy = _self.first.evaluate(relax) / 16.0217733;

        boost::shared_ptr< Crystal::TStructure<std::string> >
          result( new Crystal::TStructure<std::string> );
        Crystal::convert_real_to_string_structure(_self.second, *result);     
        result->energy = energy;
        return bp::make_tuple(result, math::rMatrix3d(_self.first.get_stress()));
      }
    template<class T>
      bp::tuple get_bond(T const &_self, const std::string &_str ) 
      {
        try
        {
          typedef boost::tuples::tuple
          < 
            const types::t_real&, const types::t_real&,
            const types::t_real&, const types::t_real&,
            const types::t_real&, const types::t_real& 
          > t_Tuple;
          const t_Tuple result( _self.first.get_bond( _str ) );
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
        return bp::tuple();
      };
    template<class T>
      bp::tuple get_angle(T const &_self, const std::string &_str ) 
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
          const t_Tuple result( _self.first.get_angle( _str ) );
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
        return bp::tuple();
      };
    template<class T>
      void set_bond(T &_self, const std::string &_str, const bp::tuple &_t ) 
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
          _self.first.set_bond( _str, boost::tuples::make_tuple( a, b, c, d, e, f ) );
        }
        catch(...)
        {
          PyErr_SetString(PyExc_RuntimeError,
                          ("could not find parameters of bond " + _str).c_str() );
          bp::throw_error_already_set();
        }
      };
    template<class T>
      void set_angle(T &_self, const std::string &_str, const bp::tuple &_t )
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
          _self.first.set_angle( _str, boost::tuples::make_tuple( a, b, c, d, e, f, g ) );
        }
        catch(...)
        {
          PyErr_SetString(PyExc_RuntimeError,
                          ("could not find parameters of bond " + _str).c_str() );
          bp::throw_error_already_set();
        }
      };

#   ifdef _MPI
      template<class T> boost::mpi::communicator const &
         get_mpi(T const &_self) { return _self.first.comm(); }
      template<class T> void set_mpi(T &_self, boost::mpi::communicator * _c)
        { _self.first.set_mpi(_c); }
#   endif

    math::rVector3d get_direction(t_LayeredVff const &_self)
      { return _self.first.Vff().get_direction(); } 
    void set_direction(t_LayeredVff &_self, math::rVector3d const &_dir)
      { _self.first.Vff().set_direction(_dir); } 

    template<class T> void set_minimizer(T &_self, typename T::first_type::t_Minimizer const &_minimizer)
      { _self.first.set_minimizer(_minimizer); }

    template<class T> bp::class_<T> expose_vff_(std::string const &_name, std::string const &_doc)
    {
      return bp::class_< T >
      ( 
        _name.c_str(),
        (
          _doc
          + "\n\nThis object can be created with: \n"
            "  - No argument.\n"
           //"  - A single string argument representing the path to an XML input file.\n" 
            "  - A single boost.mpi communicator.\n"
            "  - Another functional, eg deepcopy.\n" 
            "  - A string argument (see above), followed by a boost.mpi communicator.\n\n" 
            "If compiled without mpi, including the communicator will have no effect.\n"
        ).c_str(),
        bp::no_init
      ).def( "__init__", bp::make_constructor(&create<T>) )
       .def( "__init__", bp::make_constructor(&create_mpi<T>) )
       .def( "__init__", bp::make_constructor(&create_inputmpi<T>) )
       .def( "__init__", bp::make_constructor(&create_copy<T>) )
       .def
       ( 
         "__call__",  
         &__call__<T>,
         (bp::arg("structure"), bp::arg("doinit") = true, bp::arg("relax")=true),
         "Minimizes structure.\n\n"
         "@param structure: structure to evaluate.\n"
         "@type structure: L{crystal.Structure}.\n"
         "@param doinit: If true, reconstructs first-neighbor tree. Default: true.\n"
         "@return: A 2-tuple (structure, stress) where the first is the relaxed "
         "structure, and the second a matrix with the stress. The energy is in "
         "structure.L{energy<crystal.Structure.energy>}.\n"
       )
       .def("set_minimizer", &set_minimizer<T>)
#      ifdef _MPI
         .add_property
         (
           "comm", 
           bp::make_function(&get_mpi<T>, bp::return_internal_reference<>()),
           &set_mpi<T>, 
           "Sets the boost.mpi communicator."
         )
#      endif
       .def( "_init",  &init<T>, (bp::arg("redo_tree") = true, bp::arg("verbose")=false), 
             "Initializes the functional for the current structure." ) 
       .def( "print_escan_input",  &print_escan_input<typename T::first_type>,
             (bp::arg("file"), bp::arg("structure")), 
             "Outputs the current structure in a format suitable for pescan." ) 
       .def( "set_bond",  &set_bond<T>, ( bp::arg("bond"), bp::arg("params") ), 
             "Sets the parameters of bond from a tuple where:\n" 
             "  - the first component is the bond length.\n" 
             "  - the second to sixth components are the gammas at order 2-6.\n" 
             "The tuple can be shortened to include only the first few parameters.\n" ) 
       .def( "get_bond",  &get_bond<T>, bp::arg("angle"), 
             "Returns the parameters of bond in form of a tuple.\n" 
             "  - the first component is the bond length.\n" 
             "  - the second to sixth components are the gammas at order 2-6.\n" ) 
       .def( "set_angle",  &set_angle<T>, ( bp::arg("angle"), bp::arg("params") ), 
             "Sets the parameters of angle from a tuple where:\n" 
             "  - the first component is gamma.\n" 
             "  - the second component is sigma.\n" 
             "  - the third to seventh components are the betas at order 2-6.\n" 
             "The tuple can be shortened to include only the first few parameters.\n" ) 
       .def( "get_angle",  &get_angle<T>, bp::arg("angle"), 
             "Returns the parameters of angle in form of a tuple.\n" 
             "  - the first component is gamma.\n" 
             "  - the second component is sigma.\n" 
             "  - the third to seventh components are the betas at order 2-6.\n" );
    }

    void expose_vff()
    {
      expose_vff_<t_Vff>
      ( 
        "Vff", 
        "A Valence Force Field Functional.\n\n"
        "Usage:\n"
        ">>> vff = lada.vff.Vff(\"input.xml\", boost.mpi.world)\n"
        ">>> structure = crystal.Structure(\"input.xml\")\n"
        ">>> relaxed_structure = vff(structure)\n\n"
      );
      bp::register_ptr_to_python< boost::shared_ptr<t_Vff> >();
    }

    void expose_layeredvff()
    {
      expose_vff_<t_LayeredVff>
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
           &get_direction, 
           bp::return_value_policy<bp::return_by_value>()
         ),
         &set_direction, "Growth/Epitaxial direction.\n\n3x1 float64 numpy array.\n" 
       ); 
      bp::register_ptr_to_python< boost::shared_ptr<t_LayeredVff> >();
    }

  }
} // namespace LaDa
