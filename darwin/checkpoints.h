//
//  Version: $Id$
//
#ifndef _DARWIN_CHECKPOINT_H_
#define _DARWIN_CHECKPOINT_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <sstream>
#include <iomanip>
#include <algorithm>

#include <eo/eoPop.h>
#include <eo/utils/eoPopStat.h>
#include <eo/utils/eoStat.h>
#include <eo/utils/eoMonitor.h>
#include <eo/utils/eoUpdater.h>
#include <eo/eoContinue.h>
#include <eo/utils/eoHowMany.h>

#include <opt/types.h>
#include <print/xmg.h>
#include <print/stdout.h>
#include <mpi/mpi_object.h>

#include "taboos.h"
#include "operators.h"
#include "gencount.h"
#include "store.h"

namespace GA
{
  //! \brief Prints data to Print::xmg. 
  //! \details This functor is called at the end of each generation. It prints
  //! out new results to Pring::xmg if there are any, as well as the generation
  //! number and the number of calls to the functional (using
  //! T_EVALUATION::nb_eval, see Evaluation::Base::nb_eval). 
  template< class T_STORE, class T_EVALUATION >
  class PrintGA : public eoUpdater
  {
    public:
      typedef T_STORE t_Store; //!< Type of storage class
      typedef T_EVALUATION t_Evaluation;  //!< Type of evaluation class

    protected:
      const t_Store &store; //!< Reference to a storage class
      const t_Evaluation &evaluation; //!< Reference to an evalution class 
      const GenCount &age; //!< Reference to a generational counter
      //! \brief true if should print at each call. 
      //! \details You may not want to print if there are no new results.
      bool do_print_each_call; 

    public:
      //! Constructor and Initializor
      PrintGA   ( const t_Store &_store, const t_Evaluation &_eval, 
                  const GenCount &_age, bool _each )
              : store(_store), evaluation(_eval), age(_age), do_print_each_call(_each) {}
      //! Copy Constructor
      PrintGA   ( const PrintGA &_c )
              : store(_c.results), evaluation(_c.evaluation), age(_c.age),
                do_print_each_call(_c.do_print_each_call)  {}

      //! \brief This class is a functor
      virtual void operator()();

      //! Eo required 
      void printOn( std::ostream &__os ) const {};
      //! Eo required 
      void readFrom( std::istream &__os ) const {};
      //! Prints results to Print::xmg prior to termination
      void lastCall();

      //! EO required 
      virtual std::string className(void) const { return "GA::PrintGA"; }
  };

  //! \brief Prints out the offsprings to Print::out.
  //! \details The printout is sorted according to the \b scalar fitness
  template< class T_GATRAITS>
  class PrintOffspring : public eoSortedStatBase<typename T_GATRAITS::t_Individual> 
  {
    public:
      typedef T_GATRAITS t_GATraits; //!< All %GA %types
    protected:
      typedef typename t_GATraits::t_Individual t_Individual; //!< Type of the indiviudals
      //! \brief Type of the \e sorted population. 
      //! \details The population is sorted using the mechanism offered by EO.
      //!          As a result, the \e sorted population is actually a vector
      //!          of pointers to individuals, not the usual eoPop.
      typedef typename std::vector<const t_Individual*> t_Population;

    protected:
      GenCount &age; //!<  Reference to a generational counter

    public:
      //! Constructor and Initializor
      PrintOffspring   ( GenCount &_age )
                   : age(_age) {}
      //! Copy Constructor
      PrintOffspring   ( const PrintOffspring &_update )
                   : age(_update.age) {}

      //! \brief This class is a functor
      virtual void operator()( const t_Population &_pop );

      //! Eo required 
      void printOn( std::ostream &__os ) const {};
      //! Eo required 
      void readFrom( std::istream &__os ) const {};
      //! Eo required 
      void lastCall( const eoPop<t_Individual> &_pop) {}

      //! Eo required 
      virtual std::string className(void) const { return "GA::PrintOffspring"; }
  };
  
  //! \brief Prints out current population to Print::out
  //! \details The printout is sorted according to the \b scalar fitness
  template< class T_GATRAITS>
  class PrintPop : public eoSortedStatBase<typename T_GATRAITS::t_Individual> 
  {
    public:
      typedef T_GATRAITS t_GATraits; //!< All %GA %types
    protected:
      typedef typename t_GATraits::t_Individual t_Individual; //!< Type of the indiviudals
      //! \brief Type of the \e sorted population. 
      //! \details The population is sorted using the mechanism offered by EO.
      //!          As a result, the \e sorted population is actually a vector
      //!          of pointers to individuals, not the usual eoPop.
      typedef typename std::vector<const t_Individual*> t_Population;


    public:
      //! Constructor 
      PrintPop() {};
      //! Destructor 
      ~PrintPop() {}

      //! \brief This class is a functor
      //! \details Watch it! in \e sorted statistics t_Population is not quite
      //!          the container we are used to...
      virtual void operator()( const t_Population &_pop );

      //! Eo required 
      void printOn( std::ostream &__os ) const {};
      //! Eo required 
      void readFrom( std::istream &__os ) const {};
      //! Eo required 
      void lastCall( const eoPop<t_Individual> &_pop) {}

      //! Eo required 
      virtual std::string className(void) const { return "GA::PrintPop"; }
  };

  //! \brief In case of taboo diverges, starts a period of high mutations.
  //! \details The idea is that if no new individual can be found, than the %GA
  //! is stuck in a niche. A periode of high mutation could possibly dislodge it.
  //! \warning Everything in this class is a dirty hack. Its best simply not to use it.
  template< class T_GATRAITS>
  class NuclearWinter : public eoStatBase<typename T_GATRAITS::t_Individual>
  {
    public:
      typedef T_GATRAITS t_GATraits; //!< All %GA triats
    protected:
      typedef typename t_GATraits::t_Individual t_Individual; //!< Type of an individual
      typedef typename t_GATraits::t_Population t_Population; //!< Type of the population

    protected:
      Taboo_Base<t_Individual> &taboo; //!< Reference to the taboos to check
      eoGenOp<t_Individual> &normal_ops; //!< Reference to the "normal" operators
      SequentialOp<t_Individual> nuclear_ops; //!< The high-mutation operators
      eoGenOp<t_Individual> **breeding_ops; //!< The breeding operators. Watch out! Hack!
      types::t_unsigned nuclear_winter_length; //!< Number of generations of high mutation period
      types::t_unsigned nuclear_winter_age;//!< Starting date of the high mutation period
      bool is_gone_nuclear; //!< True if within a high mutation period
      eoHowMany nuclear_howmany; //!< Number of offsprings to create during high mutation period
      eoHowMany normal_howmany; //!< Number of offsprings to create during normal period
      eoHowMany **breeding_howmany; //!< Dirty Hack!

    public:
      //! Constructor and Initialisor
      NuclearWinter   ( Taboo_Base<t_Individual> &_taboo, 
                        eoGenOp<t_Individual> &_nops,
                        eoGenOp<t_Individual> &_nwops,
                        types::t_real &_normal_howmany )
                    : taboo(_taboo),
                      normal_ops(_nops),
                      breeding_ops(NULL),
                      nuclear_winter_length(2),
                      is_gone_nuclear(false),
                      nuclear_howmany(1.0),
                      normal_howmany(_normal_howmany),
                      breeding_howmany(NULL)
      {
        nuclear_ops.add(_nops, 1.0);
        nuclear_ops.add(_nwops, 1.0);
      }

      //! Destructor
      virtual ~NuclearWinter(){};

      
      //! Dirty Hack!
      void set_howmany( eoHowMany **_howmany);

      //! \brief This class is a dirty hack of a functor
      virtual void operator()( const t_Population &_pop );

      //! Dirty Hack!
      void set_op_address( eoGenOp<t_Individual> ** _ops ) { breeding_ops = _ops; } 

      //! Dirty Hack!
      eoGenOp<t_Individual>* get_op_address() const;

      //! EO required 
      void printOn( std::ostream &__os ) const {};
      //! EO required 
      void readFrom( std::istream &__os ) const {};
      //! EO required 
      virtual void lastCall( const t_Population &_pop) {};
      //! EO required 
      virtual std::string className(void) const { return "LaDa::NuclearWinter"; }

  };

  //! \brief adds "old" individuals to a taboo list
  //! \details Individuals are added to a taboo list here, but certainly not removed
  //!          from the population.
  template< class T_GATRAITS>
  class UpdateAgeTaboo : public eoStatBase<typename T_GATRAITS::t_Individual>
  {
    public:
      typedef T_GATRAITS t_GATraits; //!< All %GA traits.
    protected:
      typedef typename t_GATraits::t_Individual t_Individual; //!< Type of the indiviudals
      typedef typename t_GATraits::t_Population t_Population; //!< Type of the population

    protected:
      Taboo< t_Individual, std::list<t_Individual> > & taboo; //!< Reference to a Taboo list
      GenCount &age; //!< Generational counter
      types::t_unsigned max_age; //!< Maximum allowed age
      types::t_unsigned check_every; //!< Checks for the old every n generations
      bool do_print_out; //!< Prints out newly tabooed individual

    public:
      //! Constructor and Initialisor
      UpdateAgeTaboo  ( Taboo< t_Individual, std::list<t_Individual> > & _taboo,
                        GenCount &_age, types::t_unsigned _max_age, bool _do_print_out = false )
                     : taboo(_taboo), age(_age), max_age(_max_age),
                       do_print_out( _do_print_out )
      {
        check_every = max_age / 10;
        if ( check_every == 0 )
          check_every = 1;
      };
      //! Copy Constructor
      UpdateAgeTaboo  ( const UpdateAgeTaboo<t_Individual> & _update )
                     : taboo(_update.taboo), age(_update.age), 
                       max_age(_update.max_age), check_every( _update.check_every ),
                       do_print_out(_update.do_print_out) {};

      //! \brief This class is a functor
      //! \details Check for old individuals in \a _pop and adds them to the
      //! taboo list. If required, prints out the tabooed individual to
      //! Print::xmg as a comment.
      virtual void operator()( const t_Population &_pop );

      //! EO required
      void printOn( std::ostream &__os ) const {};
      //! EO required
      void readFrom( std::istream &__os ) const {};
      //! EO required
      virtual void lastCall( const eoPop<t_Individual> &_pop) {}

      //! EO required
      virtual std::string className(void) const { return "Darwin::UpdateAgeTaboo"; }
  };

  //! \brief Terminates %GA when condition is found to be true
  //! \details The condition should be a functor of type T_BINOP which takes
  //! two parameters of type T_VALUE, and returns true if the condition is
  //! fulfilled. Upon true, the %GA terminates. Here is an example taken from
  //! darwin.impl.h (at revision 313)
  //! \code
  //  eoContinue<t_Individual> *terminator = NULL;
  //  terminator = new Terminator< types::t_unsigned, std::less<types::t_unsigned>, t_GATraits >
  //                             ( evaluation->nb_eval, (types::t_unsigned) abs(max),
  //                               std::less<types::t_unsigned>(), "nb_eval < term" );
  //! \endcode
  //! \param T_VALUE type of the binary functor parameters
  //! \param T_BINOP type of the binary functor with which Terminator::ref and
  //!                Terminator::term are compared 
  //! \param T_GATRAITS All %GA traits
  template< class T_VALUE, class T_BINOP, class T_GATRAITS>
  class Terminator : public eoContinue<typename T_GATRAITS::t_Individual>
  {
    public:
      //! type of the binary functor parameters.
      typedef T_VALUE t_Value; 
      //! \brief type of the binary functor with which Terminator::ref and
      //!  Terminator::term are compared.
      typedef T_BINOP t_BinOp;         
      typedef T_GATRAITS t_GATraits; //!< All %GA traits.
    protected:
      typedef typename t_GATraits::t_Individual t_Individual; //!< Type of the indiviudals
      typedef typename t_GATraits::t_Population t_Population; //!< Type of the population

    protected:
      t_Value &ref;  //!< \brief Reference to the variable to be watched
      t_Value term;  //!< \brief Comparison parameter
      t_BinOp binop; //!< \brief Comparison binary functor 
      std::string type; //!< Terminator type

    public:
      //! Constructor and Initialisor
      //! \param _ref is the variable to watch
      //! \param _term a parameter to fulfilling the condition
      //! \param _op an instance of the condition functor
      //! \param _type a string describing the condition, for output purposes
      Terminator   ( t_Value &_ref, t_Value _term, t_BinOp _op, 
                     std::string _type)
                 : ref(_ref), term(_term), binop(_op),
                   type(_type) {}
      //! Copy Constructor
      Terminator   ( const Terminator & _copy)
                 : ref(_copy.ref), term(_copy.term), binop(_copy.op), 
                   type(_copy.type) {}
      //! Destructor
      virtual ~Terminator() {}

      //! \brief This class is a functor
      //! \details returns true if either:
      //!     - Terminator::term is false
      //!     - Terminator::term is true <STRONG >and</STRONG> Condition
      //!       Terminator::binop return true
      //!     .
      virtual bool operator() (const t_Population &_pop );
     
      //! EO required
      virtual std::string className(void) const { return "Darwin::Terminator"; }

      //! Prints out current status of the condition
      void lastCall() const;
  };

  //! \brief Collects all continuators and such, and works for range of population.
  //! \details This continuator expects on input an iterator range of whole
  //! populations, eg islands. It then applies the collected continuators,
  //! stats, updaters to each population. 
  template< class T_GATRAITS>
  class IslandsContinuator : public eoContinue<typename T_GATRAITS::t_Individual> 
  {
    public:
      typedef T_GATRAITS t_GATraits; //!< All %GA traits
    protected:
      typedef typename t_GATraits::t_Individual t_Individual; //!< Type of the indiviudals
      typedef typename t_GATraits::t_Population t_Population; //!< Type of the population
      typedef typename t_GATraits :: t_Islands t_Islands; //!< Type of island collection
      //! Type to an iterator to an island collection
      typedef typename t_Islands :: iterator iterator; 
      //! Constant type to an iterator to an island collection
      typedef typename t_Islands :: const_iterator const_iterator;
    protected:
      //! list of continuators. see EO library
      std::list < eoContinue<t_Individual>* >       continuators;
      //! list of sorted statistics. see EO library
      std::list < eoSortedStatBase<t_Individual>* > sorted;
      //! list of statistics. see EO library
      std::list < eoStatBase<t_Individual>* >       stats;
      //! list of monitors. see EO library
      std::list < eoMonitor* >                      monitors;
      //! list of updaters. see EO library
      std::list < eoUpdater* >                      updaters;
      //! A generational counter
      GenCount generation_counter;
      //! \brief Maximum number of generations before quitting
      //! \details May quit before if a something in IslandsContinuator ::
      //! continuators says quit. Actually, its generally more usefull to use a
      //! Terminator acting upon the number of evaluations than this variable. 
      types::t_unsigned  max_generations;
      //! %GA stops on finding this file
      std::string stop_filename;

    public:
      //! Constructor and Initialisor
      IslandsContinuator   ( types::t_unsigned _max, std::string _f = "stop" ) 
                         : generation_counter(0),
                           max_generations(_max), stop_filename(_f)
      {
        if ( stop_filename.empty() ) stop_filename = "stop";
        Print::xmg << Print::Xmg::comment << "Will stop on finding file "
                   << stop_filename << Print::endl;
      }


      //! adds a continuator
      void add(eoContinue<t_Individual>& _cont) { continuators.push_back(&_cont); }
      //! adds a sorted statistic
      void add(eoSortedStatBase<t_Individual>& _stat) { sorted.push_back(&_stat); }
      //! adds an unsorted statistic
      void add(eoStatBase<t_Individual>& _stat) { stats.push_back(&_stat); }
      //! adds a monitor
      void add(eoMonitor& _mon) { monitors.push_back(&_mon); }
      //! adds an updater
      void add(eoUpdater& _upd) { updaters.push_back(&_upd); }

      //! \brief Don't use
      //! \details EO declares this function virtual, but we don't need it.
      //!          use IslandsContinuator::apply() instead
      bool operator()(const t_Population& _pop) { return true; }

      //! \brief Applies sorted and unsorted statistics.
      void apply_stats(const t_Population& _pop) const;

      //! \brief Applies monitors and updaters
      void apply_monitors_updaters();

      //! \brief Applies continuators
      bool apply_continuators( const t_Population & _pop ) const;

      //! \brief When terminating, calls lastCall member function of statistics,
      //! monitors and updaters.
      void lastCall( const t_Population& _pop );

      //! \brief Applies all the apply_something member functions  to the iterator range.
      //! \details The populations in [\a _i_begin, \a _i_end) should all be valid.
      bool apply ( iterator &_i_begin, iterator &_i_end );

      //! EO required
      virtual std::string className(void) const { return "LaDa::IslandsContinuator"; }
      //! Returns a referenec to the generation counter
      GenCount& get_generation_counter()  { return generation_counter; }
      //! Returns the current number of GA generations
      types::t_unsigned age() const { return generation_counter(); }
  };

  //! \brief Calls SaveEvery::func every SaveEvery::every generations
  //! \details Used to save results every so many generations. 
  //! You can use this functor as follows (from darwin.iml.h at revision 313)
  //! \code
  //    SaveEvery<t_Darwin> *save = new SaveEvery<t_Darwin>( *this, &Darwin::Save, std::abs(n) );
  //! \endcode
  template < class T_CLASS >
  class SaveEvery: public eoUpdater
  {
    public:
      typedef T_CLASS t_Class; //!< Class of member function to call
    protected:
      typedef bool ( t_Class::*t_Function )(); //!< Member function type
    protected:
      t_Class &object; //!< Object of which to call member function
      t_Function func; //!< Address of member function to call
      types::t_unsigned every; //!< Perfom function every "n" calls
      types::t_unsigned age; //!< Call counter

    public:
      //! \brief Constructor and Initializer
      //! \param _object for which member function is called
      //! \param _func pointer to the member function to call 
      //! \param _n Initializes SaveEvery::every
      SaveEvery   ( t_Class &_object, t_Function _func, types::t_unsigned _n ) 
                : object(_object), func(_func), every( _n ), age(0) 
      {
        if ( not every ) every = 1;
      }
      //! Copy Constructor
      SaveEvery   ( const SaveEvery<t_Class> &_copy )
                : object(_copy.object), func(_copy.func), every( _copy.n ), age(_copy.age) {}

      //! \brief This class is a functor
      void operator()();
  };

#ifdef _MPI
  //! \brief Synchronizes scalars across processor at the end of each generation
  //! \details Basically does an mpi::Base::all_sum_all on Synchronize::object.
  //! No, nothing gets counted twice, so don't worry.
  template< class T_TYPE >
  class Synchronize : public eoUpdater
  {
    public:
      typedef T_TYPE t_Type; //!< The type of scalar to sync
    protected:
      t_Type &object; //!< A reference to the scalar to sync
      t_Type current_value; //!< A place holder
    public:
      //! Constructor and Initializer
      explicit
        Synchronize( t_Type &_object ) : object(_object), current_value(_object) {}
      //! Destructor
      ~Synchronize(){}
      //! This is a functor
      void operator()();
  };
#endif


  //! \brief Applies a functor to all stored individuals.
  //! \details Generally (depending on the overloading of
  //!          Store::Manip::apply_all), this should mean applying to whatever
  //!          container exists in Apply2Stored::store. This functor is not
  //!          permitted to change anything in the stored invididuals.
  template< class T_GATRAITS >
  class Apply2Stored: public eoUpdater
  {
    public:
      //! All relevant GA traits
      typedef T_GATRAITS t_GATraits;

    protected:
      //! Type of the individual
      typedef typename t_GATraits :: t_Individual t_Individual;
      //! Type of the abstract base storage class
      typedef Store::Base<t_GATraits> t_Store;
      //! Type of the functor.
      typedef eoMonOp<const t_Individual> t_Functor;
 
    protected:
      const t_Store &store; //!< Reference to the storage class
      t_Functor *functor; //!< Pointer to the functor
 
    public:
      //! \brief Constructor and (partial) Initializer.
      //! \details Apply2Stored::set_functor() still needs to be called prior
      //!          to use.
      Apply2Stored ( const t_Store &_store ) : store(_store ) {}
      //! Sets the functor to call for each stored individual.
      void set_functor( t_Functor *_functor ) { functor = _functor; }
 
      //! Functor. Reroutes calls to Store::Manip::apply_all().
      void operator()() { store.apply_all( *functor ); }
  };
 
  //! \brief Applies a functor to best stored individual.
  //! \details Generally (depending on the overloading of
  //!          Store::Manip::apply_best), this should mean applying to whatever
  //!          optimum exists in Apply2Best::store. This functor is not
  //!          permitted to change anything in the stored invididuals.
  template< class T_GATRAITS >
  class Apply2Best: public eoUpdater
  {
    public:
      //! All relevant GA traits
      typedef T_GATRAITS t_GATraits;

    protected:
      //! Type of the individual
      typedef typename t_GATraits :: t_Individual t_Individual;
      //! Type of the abstract base storage class
      typedef Store::Base<t_GATraits> t_Store;
      //! Type of the functor.
      typedef eoMonOp<const t_Individual> t_Functor;
 
    protected:
      const t_Store &store;
      t_Functor *functor; //!< Pointer to the functor
 
    public:
      //! \brief Constructor and (partial) Initializer.
      //! \details Apply2Best::set_functor() still needs to be called prior 
      //!          to use.
      Apply2Best   ( const t_Store &_store )
                 : store(_store ) {}
      //! Sets the functor to call for best stored individual.
      void set_functor( t_Functor *_functor ) { functor = _functor; }
 
      //! Functor. Reroutes calls to Store::Manip::apply_best().
      void operator()() { store.apply_best( functor ); }
  };
} // namespace GA

#include "checkpoints.impl.h"

#endif //  _CHECKPOINT_H_
