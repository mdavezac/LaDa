#ifndef _DARWIN_H
#define _DARWIN_H

//-----------------------------------------------------------------------------

#include <apply.h>
#include <eoAlgo.h>
#include <eoPopAlgo.h>
#include <eoPopEvalFunc.h>
#include <eoContinue.h>
#include <eoSelect.h>
#include <eoTransform.h>
#include <eoBreed.h>
#include <eoMergeReduce.h>
#include <eoReplacement.h>



template<class Object>
class Darwin: public eoAlgo<Object>
{
  public:

     Darwin   ( eoContinue<Object>& _continuator, eoPopEvalFunc<Object>& _eval, eoBreed<Object>& _breed,
                eoReplacement<Object>& _replace,  eoPopAlgo<Object>  & _extra_popalgo )
            : continuator(_continuator),
              popEval(_eval),
              breed(_breed),
              replace(_replace),
              extra_popalgo(_extra_popalgo )
            {}

    /// Apply a few generation of evolution to the population.
    virtual void operator()(eoPop<Object>& _pop)
    {
      eoPop<Object> offspring, empty_pop; 
      popEval(empty_pop, _pop); // A first eval of pop.
      do
      {
        try
        {
           unsigned pSize = _pop.size();
           offspring.clear(); // new offspring
    
           breed(_pop, offspring);
    
           popEval(_pop, offspring); // eval of parents + offspring if necessary
    
           replace(_pop, offspring); // after replace, the new pop. is in _pop

           extra_popalgo(_pop); // minimizes best for instance
    
           if (pSize > _pop.size())
               throw std::runtime_error("Population shrinking!");
           else if (pSize < _pop.size())
               throw std::runtime_error("Population growing!");

        }
        catch (std::exception& e)
        {
              std::string s = e.what();
              s.append( " in eoEasyEA");
              throw std::runtime_error( s );
        }
      } while ( continuator( _pop ) );
    }

  protected :
  
    eoContinue<Object>&          continuator;
    eoPopEvalFunc<Object>&       popEval;
    eoBreed<Object>&             breed;
    eoReplacement<Object>&       replace;
    eoPopAlgo<Object>&           extra_popalgo;
  
};

//-----------------------------------------------------------------------------

#endif

