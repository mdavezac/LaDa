//
//  Version: $Id$
//

#ifndef _LADA_GA_OPERATOR_EOADAPTER_H_
#define _LADA_GA_OPERATOR_EOADAPTER_H_

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif


#include <eo/eoGenOp.h>
#include <eo/utils/eoState.h>

namespace LaDa
{
  namespace GA
  {
    //! Holds general stuff for GA operators.
    namespace Operator
    {
      template< class T_FACTORY >
        eoGenOp< typename T_FACTORY :: t_Individual >*
          create_eo_operator( const TiXmlElement&, T_FACTORY&, eoState& );

      template< class T_INDIVIDUAL, class T_POPULATOR = eoPopulator<T_INDIVIDUAL> >
        class EOAdapter : public eoGenOp<T_INDIVIDUAL >
        {
          template<class T_FACTORY >
          friend eoGenOp<typename T_FACTORY::t_Individual>*
            create_eo_operator( const TiXmlElement&_node, T_FACTORY&, eoState &_eostate );
          public:
            //! Type of the individual.
            typedef T_INDIVIDUAL t_Individual;
            //! Type of the populator.
            typedef T_POPULATOR t_Populator;
            //! Constructor.
            EOAdapter() {}
            //! Copy Constructor.
            EOAdapter( const EOAdapter& _c ) : callback_( _c.callback_) {}
            //! Number of produced offspring
            virtual types::t_unsigned max_production() { return 1; }

            //! \brief Tries to create an non-tabooed object by applying _op
            //! \details After max tries, creates a random untaboo object
            virtual void apply( t_Populator &_indiv )
              { callback_( _indiv ); }
            
            //! returns GA::TabooOp
            virtual std::string className () const { return "GA::Operator::Adapter"; }

            //! Connects the callback.
            template< class T_FUNCTOR > void connect( T_FUNCTOR& _functor )
              { callback_ = _functor; } 
            
          private:
            //! The callback object.
            boost::function<void( t_Populator& ) > callback_;
        };

      //! Creates an eoGenOp functor from a factory.
      template< class T_FACTORY >
        eoGenOp< typename T_FACTORY :: t_Individual >*
          create_eo_operator( const TiXmlElement& _node, T_FACTORY& _factory, eoState &_eostate )
          {
            typedef EOAdapter< typename T_FACTORY::t_Individual,
                               typename T_FACTORY::t_Populator > t_Adapter;
            t_Adapter *result = new t_Adapter;
            _eostate.storeFunctor( result );
            _factory( result->callback_, _node ); 
            return result;  
          }
    }
  }
}

#endif
