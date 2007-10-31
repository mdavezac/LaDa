//
//  Version: $Id$
//
#ifndef _CE_H_
#define _CE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <algorithm>
#include <functional>
#ifndef __PGI
  #include<ext/algorithm>
  using __gnu_cxx::copy_n;
#endif

#include <tinyxml/tinyxml.h>
#include <eo/eoOp.h>

#include "lamarck/structure.h"
#include "lamarck/functional_builder.h"
#include "opt/opt_function_base.h"
#include "opt/types.h"

#include "evaluator.h"
#include "individual.h"
#include "single_site.h"

#ifdef _MPI
#include "mpi/mpi_object.h"
#endif

namespace CE
{
  struct Object: public SingleSite::Object
  {
#ifdef _MPI
    friend bool mpi::BroadCast::serialize<CE::Object>( CE::Object & );
#endif
  };

  typedef Individual::Types< SingleSite::Object, 
                             SingleSite::Concentration, 
                             SingleSite::Fourier        > :: Scalar t_Individual;

  class Evaluator : public SingleSite::Evaluator< t_Individual >,
                    public VA_CE :: Functional_Builder
  {
    public:
      typedef CE::t_Individual t_Individual;
      typedef Traits::GA< Evaluator > t_GATraits;
    protected:
      typedef Evaluator t_This;
      typedef SingleSite::Evaluator< t_Individual > t_Base;
      typedef VA_CE::Functional_Builder::t_VA_Functional t_Functional;
      
    public:
      using t_Base::Load;
      using t_Base::Save;

    protected:
      using t_Base :: current_individual;
      using t_Base :: current_object;

    protected:
      t_Functional functional;

    public:
      Evaluator() : functional() {}; 
      ~Evaluator() 
      {
        if ( functional.get_functional1() ) 
          delete functional.get_functional1();
        if ( functional.get_functional2() ) 
          delete functional.get_functional2();
      };

      bool Save( const t_Individual &_indiv, TiXmlElement &_node, bool _type ) const;
      bool Load( t_Individual &_indiv, const TiXmlElement &_node, bool _type );
      bool Load( const TiXmlElement &_node );
      void* Load_Niche( const TiXmlElement &_node )
        { return (void *) SingleSite::new_Niche_from_xml<t_GATraits, 1>( _node ); }

      void init( t_Individual &_indiv )
        { t_Base :: init( _indiv ); functional.set_variables( &_indiv.Object().bitstring ); }

      void evaluate()
        { current_individual->quantities() = functional.evaluate(); }
      void evaluate_gradient( t_QuantityGradients& _grad )
      {
        // note that t_Function sets gradient to 0
        types::t_unsigned N = current_object->Container().size();
        types::t_real *gradient = new types::t_real[N];
        types::t_real *keep = gradient;
        if ( not gradient ) throw std::runtime_error( "Could not allocate memory" );
        functional.evaluate_gradient( gradient );
        copy_n( gradient, N, _grad.begin() );
        delete[] keep;
      }
      void evaluate_with_gradient( t_QuantityGradients& _grad )
      {
        types::t_unsigned N = current_object->Container().size();
        types::t_real *gradient = new types::t_real[N];
        types::t_real *keep = gradient;
        if ( not gradient ) throw std::runtime_error( "Could not allocate memory" );
        current_individual->quantities() = functional.evaluate_with_gradient( gradient );
        copy_n( gradient, N, _grad.begin() );
        delete[] keep;
      }
      void evaluate_one_gradient( t_QuantityGradients& _grad, types::t_unsigned _pos ) 
      {
        _grad[_pos] = functional.evaluate_one_gradient( _pos );
      }
  };
} // namespace CE

#ifdef _MPI
namespace mpi
{
  template<>
  inline bool mpi::BroadCast::serialize<CE::Object>( CE::Object & _object )
  {
    return serialize<SingleSite::Object>( _object );
  }
}
#endif


#endif // _CE_OBJECT_H_
