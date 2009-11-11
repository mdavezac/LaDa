//
//  Version: $Id$
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <mpi/mpi_object.h>

#include "../functional_builder.h"
#include "../constituent_strain.h"
#include "../harmonic.h"


namespace LaDa
{
  namespace Python
  {
    namespace detailsCE
    {
      template< class T_HARMONIC >
        class Functional : public CE::Builder<T_HARMONIC> __MPICODE( __COMMA__ public MPI_COMMDEC ) 
        {
          public:
            using CE::Builder<T_HARMONIC>::add_equivalent_clusters;
            using CE::Builder<T_HARMONIC>::clusters;
            using CE::Builder<T_HARMONIC>::lattice;
            void read_clusters( std::string const& _path );
            void clear_clusters() { clusters->clear(); }
            std::vector<CE::Cluster> const & get_clusters() const { return *clusters; }
        };
      
      template< class T_HARMONIC >
        void Functional<T_HARMONIC> :: read_clusters( std::string const& _path )
        {
          namespace fs = boost::filesystem;  

          std::ifstream file( _path.c_str(), std::ifstream::in );
          std::string line;
        
          CE::Cluster cluster;
          while( CE::read_cluster( *lattice, file, cluster ) )
            { clusters->push_back( cluster ); }
          add_equivalent_clusters();
        }
        
      template<class T_HARMONIC>
        void load_builder( Functional<T_HARMONIC>& _functional, const std::string &_filename )
        {
          namespace bp = boost::python;
          TiXmlDocument doc( _filename ); 
          TiXmlHandle docHandle( &doc ); 
        
          if( not doc.LoadFile() ) 
          {
            PyErr_SetString
            (
              PyExc_IOError, 
              ("Could not open/parse " + _filename + ": " + doc.ErrorDesc() + "\n" ).c_str() 
            );
            bp::throw_error_already_set();
            return;
          }
          const TiXmlElement* parent = docHandle.FirstChild("Job").Element();
          if( not parent )
          {
            PyErr_SetString
            (
              PyExc_IOError, 
              ("Could not find <Job> </Job> tags in " + _filename + "\n").c_str() 
            );
            bp::throw_error_already_set();
            return;
          }

          try
          {
            if( not _functional.Load( *parent ) )
            {
              PyErr_SetString
              ( 
                PyExc_IOError, 
                ("Could not load ce functional from " + _filename + "\n").c_str() 
              );
              bp::throw_error_already_set();
              return;
            }
            _functional.add_equivalent_clusters();
          }
          catch( std::exception &_e )
          {
            PyErr_SetString
            (
              PyExc_IOError, 
              (   "Could not create ce functional from input-file " + _filename + ".\n" 
                + std::string( _e.what() ) + "\n" ).c_str() 
            );
            bp::throw_error_already_set();
          }
          catch( ... )
          {
            PyErr_SetString
            (
              PyExc_IOError, 
              ("Could not create ce functional from input-file " + _filename + ".\n").c_str() 
            );
            bp::throw_error_already_set();
          }
        }

     
      template< class T_HARMONIC >
        std::pair
        <
          typename Functional<T_HARMONIC>::t_Chemical*,
          typename Functional<T_HARMONIC>::t_CS*
        > create( const Functional<T_HARMONIC> &_functional, const Crystal::Structure &_str )
        {
          namespace bp = boost::python;
          typedef std::pair 
                  <
                    typename Functional<T_HARMONIC>::t_Chemical*,
                    typename Functional<T_HARMONIC>::t_CS*
                  > t_Pair;
     
          if( not _str.atoms.size() ) 
          {
            PyErr_SetString( PyExc_RuntimeError, "Structure is empty.\n" );
            bp::throw_error_already_set();
            return t_Pair(NULL, NULL);
          }
          try
          {
            if ( not _str.k_vecs.size() )
            {
              std::auto_ptr<Crystal::Structure> str( new Crystal::Structure( _str ) );
              str->find_k_vectors();
              return create( _functional, *str );
            }
          }
          catch(...)
          { 
            bp::throw_error_already_set();
            return t_Pair(NULL, NULL);
          }
     
          try { return _functional.generate_functional( _str ); }
          catch( std::exception &_e )
          {
            PyErr_SetString( PyExc_IOError, "Could not evaluate CE functional.\n" ); 
            bp::throw_error_already_set();
            return t_Pair(NULL, NULL);
          }
        }
     
      template< class T_HARMONIC >
        types::t_real call_which( Functional<T_HARMONIC> &_functional,
                                  const Crystal::Structure &_str, size_t _which )
        {
          namespace bp = boost::python;
          typedef typename Functional<T_HARMONIC>::t_Chemical :: t_Container t_Container;
          typedef std::pair 
                  <
                    typename Functional<T_HARMONIC>::t_Chemical*,
                    typename Functional<T_HARMONIC>::t_CS*
                  > t_Pair;
          t_Pair pair = create( _functional, _str );
          if( pair.first == NULL ) 
          {
            bp::throw_error_already_set();
            return 0;
          }
          try
          {
            t_Container container;
            foreach( const Crystal::Structure::t_Atom &atom, _str.atoms )
              container.push_back( atom.type );
            types::t_real result(0);
            if( _which & 1 ) 
            {
              __MPICODE( pair.first->set_mpi( &_functional.comm() )  );
              pair.first->set_variables( &container );
              result += pair.first->evaluate();
            }
            if( _which & 2 )
            {
              __MPICODE( pair.second->set_mpi( &_functional.comm() )  );
              pair.second->set_variables( &container );
              result += pair.second->evaluate();
            }
            delete pair.first;
            delete pair.second;
            return result;
          }
          catch(std::exception &_e)
          {
            delete pair.first;
            delete pair.second;
            PyErr_SetString
            ( 
              PyExc_IOError, 
              ("Encountered error while computing CE: " + std::string(_e.what()) + "\n" ).c_str()
            );
            bp::throw_error_already_set();
            return 0;
          }
          catch(...)
          {
            delete pair.first;
            delete pair.second;
            PyErr_SetString
            ( 
              PyExc_IOError, 
              "Encountered error while computing CE.\n"
            );
            bp::throw_error_already_set();
            return 0;
          }
        }

    template< class T_HARMONIC, int WHICH >
      types::t_real call( Functional<T_HARMONIC> &_functional,
                          const Crystal::Structure &_str )
      { return call_which( _functional, _str, WHICH ); }
    template< class T_HARMONIC, int WHICH >
      types::t_real call_str( Functional<T_HARMONIC> &_functional,
                              const Crystal::TStructure<std::string> &_str )
      {
        Crystal::Structure str;
        Crystal::convert_string_to_real_structure(_str, str);
        return call<T_HARMONIC, WHICH>( _functional, str ); 
      }

#     ifndef _MPI
        template<class T_HARMONIC>
          void set_mpi( Functional<T_HARMONIC> const &, boost::python::object const & ) {}
#     endif
     
      template<class T_HARMONIC>
        void expose_ce_functional( const std::string &_name, const std::string &_docstring )
        {
          namespace bp = boost::python;
          typedef Functional<T_HARMONIC> t_Builder;
          bp::class_< t_Builder >( _name.c_str(), _docstring.c_str() )
            .def( bp::init<const t_Builder&>() )
#           ifdef _MPI
              .def( "set_mpi", &Functional<T_HARMONIC>::set_mpi )
#           else
              .def( "set_mpi", &set_mpi<T_HARMONIC> )
#           endif
            .def( "__call__", &call<T_HARMONIC, 3>, "Computes CE + CS." )
            .def( "cs", &call<T_HARMONIC, 2>, "Computes  CS.")
            .def( "chemical", &call<T_HARMONIC, 1>, "Computes CE." )
            .def( "__call__", &call_str<T_HARMONIC, 3>, "Computes CE + CS." )
            .def( "cs", &call_str<T_HARMONIC, 2>, "Computes  CS.")
            .def( "chemical", &call_str<T_HARMONIC, 1>, "Computes CE." )
            .def( "load", &load_builder<T_HARMONIC>, "Loads functional from XML.")
            .def( "read_clusters", &Functional<T_HARMONIC>::read_clusters )
            .def( "clear_clusters", &Functional<T_HARMONIC>::read_clusters )
            .add_property
            ( 
               "clusters", 
               bp::make_function( &Functional<T_HARMONIC>::get_clusters, 
                                  bp::return_internal_reference<>() )
            );
        }
    } // namespace detailsCE
  } // namespace Python
} // namespace LaDa
