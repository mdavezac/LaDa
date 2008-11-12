//
//  Version: $Id$
//
#ifndef _DARWIN_TABOO_H_
#define _DARWIN_TABOO_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<boost/function.hpp>
#include <list>
#include <algorithm>

#include <eo/eoFunctor.h>
#include <eo/eoOp.h>
#include <eo/eoGenOp.h>

#include <opt/types.h>
#include <mpi/mpi_object.h>

#include "functors.h"
#include "gatraits.h"

namespace LaDa
{
  namespace GA
  {

    //! \brief %Taboo base class
    //! \details %Taboos are functors which simply returns true if the individual
    //!          in argument should not be allowed to proceed in the population.
    //!          This class merely declares the relevant virtual interface.
    template<class T_INDIVIDUAL>
    class Taboo_Base : public eoUF<const T_INDIVIDUAL&, bool>
    {
      protected:
        //! The type of the individual
        typedef T_INDIVIDUAL t_Individual;

      public:
        //! Constructor
        Taboo_Base() {}
        //! Copy constructor
        Taboo_Base( const Taboo_Base<t_Individual> &_taboo ) {}
        //! Destructor
        virtual ~Taboo_Base() {};

        //! The number of produced individuals... Not sure why it's one...
        types::t_unsigned max_production(void) { return 1; }

        //! \cond
        virtual bool is_problematic() const {return false;}
        virtual void set_problematic(bool _p = false) {return; }
        //! \endcond

        //! prints out all tabooed individuals or whatever
        virtual void print_out( std::ostream &str ) const {return;}
        //! Loads a taboo from XML
        virtual bool Load( const TiXmlElement &_node ) { return true; }
    };

    //! A taboo class defined around a list of tabooed individuals
    template< class T_INDIVIDUAL, class T_CONTAINER = std::list<T_INDIVIDUAL> >
    class Taboo : public Taboo_Base<T_INDIVIDUAL>
    {
      public:
        //! The type of the individual
        typedef T_INDIVIDUAL t_Individual;
        //! The type of the container of tabooed individual
        typedef T_CONTAINER t_Container;
        //! The type of the individual
        typedef t_Individual value_type;

      protected:
        //! \cond
        mutable bool problematic;
        //! \endcond

        //! Whether the container is owned by this instance
        bool owns_pointer;
        //! The container of tabooed individuals
        t_Container *taboo_list;

      public:
        //! Constructs taboo around \a _list
        Taboo   ( t_Container *_list )
              : problematic(false), owns_pointer( false ),
                taboo_list(_list) {};
        //! Constructor. Creates a new Taboo::taboo_list object.
        Taboo() : problematic(false), owns_pointer(true)
          { taboo_list = new t_Container; }
        //! Copy Constructor
        Taboo   ( Taboo<t_Individual, t_Container> & _taboo )
              : owns_pointer( false ),
                taboo_list(_taboo.taboo_list) {};
        //! Destructor. If owned, deletes Taboo::taboo_list;
        virtual ~Taboo();

        //! returns true if _indiv is in taboo_list
        virtual bool operator()( const t_Individual& _indiv );

        //! \brief Adds an individual to the tabooed list. 
        //! \details If \a add_fast is false, then \a _indiv is added only if it
        //!          is not already there.
        void add( const t_Individual &_indiv, bool add_fast = true );
        //! \brief Adds an individual to the tabooed list after checking it is
        //!        not already there.
        void push_back( const t_Individual &_indiv )
          { add( _indiv, false ); }

        //! Prints all tabooed individuals in the list
        virtual void print_out( std::ostream &str ) const;

        //! \cond
        virtual bool is_problematic() const
          { return problematic; }
        virtual void set_problematic( bool _p = false ) const
          { problematic = _p; }
        //! \endcond

        //! Appends a complete container of individuals to the tabooed list.
        template<class tt_Container>
        void append( const tt_Container &_pop );
        //! returns a constant iterator to the beginning of the list
        typename t_Container :: const_iterator begin() const
          { return taboo_list->begin(); } 
        //! returns an iterator to the beginning of the list
        typename t_Container :: iterator begin() 
          { return taboo_list->begin(); } 
        //! returns a contant iterator to the end of the list
        typename t_Container :: const_iterator end() const
          { return taboo_list->end(); } 
        //! returns an iterator to the end of the list
        typename t_Container :: iterator end() 
          { return taboo_list->end(); } 
        //! Returns the number of tabooed individuals
        types::t_unsigned size() const { return taboo_list->size(); }
        //! Clears the list of tabooed individuals
        void clear() { taboo_list->clear(); }
    };

    //! Creates a taboo list from the offspring population
    template<class T_GATRAITS>
    class OffspringTaboo : public Taboo<typename T_GATRAITS::t_Individual,
                                        typename T_GATRAITS :: t_Population >
    {
      public:
        //! All relevant %GA traits
        typedef T_GATRAITS t_GATraits;
        //! The individual type
        typedef typename t_GATraits::t_Individual value_type;
      protected:
        //! The individual type
        typedef typename t_GATraits::t_Individual t_Individual;
        //! The population type
        typedef typename t_GATraits::t_Population t_Population;
        using Taboo<t_Individual, t_Population> :: taboo_list;

      public:
        //! Constructor 
        OffspringTaboo   ( t_Population *_list )
                       : Taboo<t_Individual, t_Population>( _list ) {}
  //     OffspringTaboo () : Taboo<t_Individual, t_Population>() {}
        //! Destructor
        virtual ~OffspringTaboo() {};
         
        //! returns true if _indiv is in taboo_list 
        virtual bool operator()( const t_Individual& _indiv );
    };

    //! \brief Creates a list of historically previously individuals. 
    //! \details Can be used directly as a taboo. Can also be used to recover the
    //!          quantities of previously assessed individuals.
    template<class T_INDIVIDUAL, class T_CONTAINER = std::list<T_INDIVIDUAL> >
    class History : public Taboo<T_INDIVIDUAL, T_CONTAINER>
    {
      public:
        //! Type of the individuals
        typedef T_INDIVIDUAL t_Individual;
        //! Type of the contaienr of previously assessed individuals
        typedef T_CONTAINER t_Container;
        //! Type of the individuals
        typedef t_Individual value_type;
      private:
        //! All %types relevant to an individual
        typedef typename t_Individual::t_IndivTraits t_IndivTraits;
      protected:
        using Taboo<t_Individual, t_Container> :: taboo_list;
        using Taboo<t_Individual, t_Container> :: owns_pointer;
        using Taboo<t_Individual, t_Container> :: problematic;

      public:
        //! Constructor
        History() : Taboo<t_Individual, t_Container>() {}
        //! Destructor
        virtual ~History() {};

        //! If \a _indiv already exists, copy the quantities and fitness()
        virtual bool clone(t_Individual &_indiv);
    };

    //! Container class for multiple taboos.
    template<class T_INDIVIDUAL>
    class Taboos : public Taboo_Base<T_INDIVIDUAL>
    {
      public:
        //! Type of the individuals
        typedef T_INDIVIDUAL t_Individual;

      protected:
        //! typeded to taboo base
        typedef Taboo_Base<t_Individual> t_Type;
        //! Type of the container of taboos 
        typedef std::list< t_Type* > t_Container;

      protected: 
        //! Container with the taboos
        t_Container taboos;

      public:
        //! Constructor
        Taboos() {};
        //! Copy Constructor
        Taboos( const Taboos<t_Individual> &_taboo ) : taboos( _taboo.taboos ) {};
        //! Destructor
        virtual ~Taboos(){};

        //! Number of taboos in the container
        types::t_unsigned size() const { return taboos.size(); }
        //! Returns a pointer to the first taboo
        t_Type* front() { return taboos.front(); }

        //! Adds a taboo 
        void add( Taboo_Base<t_Individual> * _taboo );
        //! Removes all taboo
        void clear() { taboos.clear(); } 

        //! Returns true as soon as one taboo operator returns true,
        virtual bool operator()( const t_Individual &_indiv );

        //! \cond
        virtual bool is_problematic() const;
        void set_problematic( bool _p = false );
        //! \endcond

        //! Forwars print out request to each taboo
        virtual void print_out( std::ostream &str ) const;
    };
    
    //! Taboo for a list of populations
    template<class T_GATRAITS>
    class IslandsTaboos : public Taboo_Base<typename T_GATRAITS::t_Individual>
    {
      public:
        //! All relevant %GA traits
        typedef T_GATRAITS t_GATraits;
      private:
        //! Type of this class
        typedef IslandsTaboos<t_GATraits>  t_This;
        //! Type of the individual
        typedef typename t_GATraits::t_Individual t_Individual;
        //! Type of the population
        typedef typename t_GATraits :: t_Population  t_Container;
        //! Type of the population container
        typedef typename t_GATraits :: t_Islands     t_Islands;

      protected: 
        //! \cond
        bool problematic;
        //! \endcond

        //! Reference to a list of tabooed populations
        t_Islands &populations;

      public:
        //! Constructor
        IslandsTaboos   ( t_Islands &_islands )
                      : problematic(false), 
                        populations( _islands ) {};
        //! Copy Constructor
        IslandsTaboos   ( const t_This &_taboos )
                      : problematic(_taboos.is_problematic), 
                        populations( _taboos.populations ) {};
        //! Destrcutor
        virtual ~IslandsTaboos(){};

        //! Returns true as soon as \a _indiv is found in one of the populations
        virtual bool operator()( const t_Individual &_indiv );
        //! Does nothing
        virtual void print_out( std::ostream &str ) const {}

        //! \cond
        virtual bool is_problematic() const
          { return problematic; }
        virtual void set_problematic( bool _p = false ) 
          { problematic = _p; }
        // \endcond
    };

    //! \brief Genetic operator which chechs wether a new offspring is tabooed
    //! \details When a newly generated offspring is found to be taboo, the
    //!          %function loops to generating a new one. This goes on for a
    //!          while. If no new individual can be found, a warning is issued and
    //!          the fucntor tries to create a non-taboo individual with
    //!          TabooOp::utterrandom. This also goes for while. After which, the
    //!          functor throws an error.
    template<class T_INDIVIDUAL>
    class TabooOp : public eoGenOp<T_INDIVIDUAL>
    {
      public:
        //! Type of the individual
        typedef T_INDIVIDUAL t_Individual; 

      protected:
        //! Reference to the taboos
        Taboo_Base<t_Individual> &taboo;
        //! Reference to a random initializer
        eoMonOp<t_Individual> &utterrandom;
        //! Maximum number of tries before going random
        types::t_unsigned max;
        //! The genetic operators for generating offsprings
        eoGenOp<t_Individual> &op;

      public:
        //! Constructor
        TabooOp   ( eoGenOp<t_Individual> &_op, 
                    Taboo_Base<t_Individual> &_taboo,
                    types::t_unsigned _max,
                    eoMonOp<t_Individual> &_ur )
                : taboo(_taboo), utterrandom(_ur), max(_max), op(_op) {}

        //! Number of produced offspring
        virtual types::t_unsigned max_production()
          { return op.max_production(); }

        //! \brief Tries to create an non-tabooed object by applying _op
        //! \details After max tries, creates a random untaboo object
        virtual void apply( eoPopulator<t_Individual> &_indiv );

        //! returns GA::TabooOp
        virtual std::string className () const { return "GA::TabooOp"; }

    };

    //! \brief Wraps a taboo over a %GA::Evaluator member %function.
    //! \details This is useful for implementing application specific taboos, say
    //!          based on phenotype.
    template< class T_EVALUATOR >
    class TabooFunction : public Taboo_Base< typename T_EVALUATOR::t_Individual >
    {
      public:
        //! Type of the evaluator
        typedef T_EVALUATOR t_Evaluator;
        //! Pointer to the member %function
        typedef bool ( t_Evaluator::*t_Function )(typename t_Evaluator::t_Individual &);
      protected:
        //! Type of the individual
        typedef typename t_Evaluator::t_Individual t_Individual;

      protected:
        //! Reference to the evaluator
        t_Evaluator &evaluator;
        //! Function over which to wrap a taboo
        t_Function member_func;
        //! String characterizing this taboo
        std::string class_name;

      public:
        //! Constructor
        explicit
          TabooFunction   ( t_Evaluator &_eval, t_Function _func, const std::string &_cn )
                        : evaluator(_eval), member_func(_func), class_name(_cn) {};
        //! Destructor
        ~TabooFunction() {}

        //! Returns the string characterizing this taboo
        std::string className() const { return class_name; }

        //! Returns the application specific taboo.
        bool operator()( const t_Individual& _indiv )
          { return ( (evaluator.*member_func) )( _indiv); }
    };

    template< class T_INDIVIDUAL, class T_POPULATOR >
      void taboo_op( T_POPULATOR &_populator, 
                     const boost::function< Taboo_Base<T_INDIVIDUAL>*() >&,
                     const boost::function<void(T_POPULATOR& )>& _inner,
                     const boost::function<void(T_POPULATOR& )>& _random );

    namespace Factory
    {
      //! Creates a taboo operator.
      template< class T_FACTORY >
        void taboo_op( T_FACTORY &_factory,
                       boost::function<void( typename T_FACTORY::t_Populator& )>&,
                       const TiXmlElement &_node,
                       const boost::function
                             < Taboo_Base<typename T_FACTORY::t_Individual>*() >&,
                       const std::string& _random );
    }

  } // namespace GA
} // namespace LaDa

#include "taboos.impl.h"

#endif
