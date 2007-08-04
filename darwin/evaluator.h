#ifndef _DARWIN_EVALUATOR_H_
#define _DARWIN_EVALUATOR_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>

#include <eo/eoGenOp.h>

#include "opt/types.h"
#include "opt/opt_function_base.h"
#include "loadsave.h"

namespace darwin 
{
  template< class T_INDIVIDUAL >
  class Evaluator
  {
    public:
      typedef T_INDIVIDUAL t_Individual;
      typedef typename t_Individual::t_Object t_Object;
      // DON'T FORGET THE FOLLOWING
      // add the following type in your derived class
      // typedef function::Function<?,?> t_Functional;
      // where ? is probably types::t_real or std::complex<...>
      // and ? some container, probably an std::vector

    public:
      Evaluator() {};
      virtual ~Evaluator() {}

      // Should load t_Individual and funtional related stuff
      // except attributes in <GA > tag
      bool Load ( std::string const &_f )
      {
        TiXmlDocument doc( _f.c_str() ); 
        TiXmlHandle docHandle( &doc ); 
        if  ( !doc.LoadFile() )
        { 
          std::cout << doc.ErrorDesc() << std::endl; 
          throw "Could not load input file in CE::Evaluator ";
        } 

        TiXmlElement *child = docHandle.FirstChild("Job").Element();
        if (not child)
          return false;
        return Load(*child);
      }
      virtual bool Load ( const TiXmlElement &_node ) = 0;
      // Load and Save individuals
      // _type can be either darwin::LOADSAVE_SHORT or
      // darwin::LOADSAVE_LONG. Results are save as latter, and darwin
      // internal stuff as the former. You need both only if the
      // ga object and the user-expected result object are different (say
      // bitstring versus a decorated lattice structure )
      virtual bool Load ( t_Individual &_indiv, const TiXmlElement &_node, bool _type ) = 0;
      virtual bool Save ( const t_Individual &_indiv, TiXmlElement &_node, bool _type ) const = 0;
      // attributes from <GA > tag in input.xml are passed to this
      // function from Darwin::Load(...)
      virtual void LoadAttribute ( const TiXmlAttribute &_att ) = 0;
      // returns a pointer to an eoOp object
      // pointer is owned by darwin::Darwin::eostates !!
      // don't deallocate yourself
      virtual eoOp<t_Individual>* LoadGaOp(const TiXmlElement &_el ) = 0;
      // returns a pointer to a eoF<bool> object
      // pointer is owned by darwin::Darwin::eostates !!
      // don't deallocate yourself
      virtual eoF<bool>* LoadContinue(const TiXmlElement &_el ) { return NULL; }
      // returns a pointer to a eoMonOp<const t_Individual> object
      // pointer is owned by darwin::Darwin::eostates !!
      // don't deallocate yourself
      virtual eoMonOp<const t_Individual>* LoadTaboo(const TiXmlElement &_el ) { return NULL; }
      // Initializes object before call to functional 
      // i.e. transforms t_Individual format to
      // function::Function<...>::variables format if necessary
      virtual bool initialize( t_Individual &_indiv ) = 0;
      // Called before objective function is evaluated
      // must return a void pointer to functional
      virtual void init( t_Individual &_indiv ) = 0;
      virtual void evaluate( t_Individual &_indiv );
  };


}

#endif // _EVALUATOR_H_
