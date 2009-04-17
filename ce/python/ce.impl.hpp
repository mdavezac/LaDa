//
//  Version: $Id$
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "../functional_builder.h"
#include "../constituent_strain.h"
#include "../harmonic.h"


namespace LaDa
{
  namespace Python
  {
    namespace detailsCE
    {
      template<class T_HARMONIC>
        void load_builder( CE::Builder<T_HARMONIC>& _functional, const std::string &_filename )
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
          typename CE::Builder<T_HARMONIC>::t_Chemical*,
          typename CE::Builder<T_HARMONIC>::t_CS*
        > create( const CE::Builder<T_HARMONIC> &_functional, const Crystal::Structure &_str )
        {
          namespace bp = boost::python;
          typedef std::pair 
                  <
                    typename CE::Builder<T_HARMONIC>::t_Chemical*,
                    typename CE::Builder<T_HARMONIC>::t_CS*
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
        types::t_real call_which( const CE::Builder<T_HARMONIC> &_functional,
                                  const Crystal::Structure &_str, size_t _which )
        {
          namespace bp = boost::python;
          typedef typename CE::Builder<T_HARMONIC>::t_Chemical :: t_Container t_Container;
          typedef std::pair 
                  <
                    typename CE::Builder<T_HARMONIC>::t_Chemical*,
                    typename CE::Builder<T_HARMONIC>::t_CS*
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
              pair.first->set_variables( &container );
              result += pair.first->evaluate();
            }
            if( _which & 2 )
            {
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

    template< class T_HARMONIC >
      types::t_real call_all( const CE::Builder<T_HARMONIC> &_functional,
                              const Crystal::Structure &_str )
      { return call_which( _functional, _str, 3 ); }
    template< class T_HARMONIC >
      types::t_real call_chem( const CE::Builder<T_HARMONIC> &_functional,
                              const Crystal::Structure &_str )
      { return call_which( _functional, _str, 1 ); }
    template< class T_HARMONIC >
      types::t_real call_cs( const CE::Builder<T_HARMONIC> &_functional,
                             const Crystal::Structure &_str )
      { return call_which( _functional, _str, 2 ); }

     
      template<class T_HARMONIC>
        void expose_ce_functional( const std::string &_name, const std::string &_docstring )
        {
          namespace bp = boost::python;
          typedef CE::Builder<T_HARMONIC> t_Builder;
          bp::class_< t_Builder >( _name.c_str(), _docstring.c_str() )
            .def( bp::init<const t_Builder&>() )
            .def( "__call__", &call_all<T_HARMONIC>, "Computes CE + CS." )
            .def( "chemical", &call_chem<T_HARMONIC>, "Computes CE." )
            .def( "cs", &call_cs<T_HARMONIC>, "Computes  CS.")
            .def( "load", &load_builder<T_HARMONIC>, "Loads functional from XML.");
        }
    } // namespace detailsCE
  } // namespace Python
} // namespace LaDa
