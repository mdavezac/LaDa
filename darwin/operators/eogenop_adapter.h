namespace LaDa
{
  namespace GA
  {
    //! Holds general stuff for GA operators.
    namespace Operator
    {
      template< class T_INDIVIDUAL >
        class EOAdapter : public eoGenOp<T_INDIVIDUAL >
        {
          public:
            //! Type of the individual.
            typedef T_INDIVIDUAL t_Individual;
            //! Constructor.
            EOAdapter() {}
            //! Copy Constructor.
            EOAdapter( const EOAdapter& _c ) : callback_( _c.callback_) {}
            //! Number of produced offspring
            virtual types::t_unsigned max_production() { return 1; }

            //! \brief Tries to create an non-tabooed object by applying _op
            //! \details After max tries, creates a random untaboo object
            virtual void apply( eoPopulator<t_Individual> &_indiv )
              { callback_( _indiv ); }
            
            //! returns GA::TabooOp
            virtual std::string className () const { return "GA::Operator::Adapter"; }

            //! Connects the callback.
            template< class T_FUNCTOR > void connect( T_FUNCTOR& _functor )
              { callback_ = _functor; } 
            
          private:
            //! The callback object.
            boost::function<bool(t_Individual&, const t_Individual&) > callback_;
        };
    }
  }
}
