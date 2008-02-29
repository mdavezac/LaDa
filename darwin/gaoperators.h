//
//  Version: $Id$
//

#ifndef _GAOPERATORS_H_
#define _GAOPERATORS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <algorithm>
#include <functional>
#include <string>
#include <sstream>

#include <eo/eoOp.h>

#include <tinyxml/tinyxml.h>

#include <lamarck/structure.h>
#include <print/stdout.h>
#include <opt/types.h>
#include <mpi/mpi_object.h>

#include "functors.h"
#include "gatraits.h"


/** \ingroup Genetic 
 * @{*/

//! \brief Contains all Genetic Algorithm that is not in another namespace ;)
//! \todo make GA a more homogeneous namespace. Possibly move all %GA related
//! namespaces inside GA
namespace GA
{
  //! \brief Applies reciprocal-space %crossover. See <A
  //! HREF="http://dx.doi.org/10.1088/0953-8984/19/40/402201"> J. Phys.: Cond.
  //! Matt. <STRONG>19</STRONG>, 402201 (2007) </A>.
  //! \details Krossover refers to a crossover operation in reciprocal-space
  //! where the values of the structure factors of two parents of identical shape
  //! are interchanged to create a new individual. This process should be more
  //! seamless than the real-space crossover practiced in GA::Crossover.
  //! Krossover can be done over all <STRONG>k</STRONG>-vectors indifferently,
  //! or the parents can exchange a contiguous range of vectors with increasing
  //! wavelengths. Whether range or not does not seem to affect the results.
  //! The concentration after krossover is set (or not) using a
  //! \a T_INDIVIDUAL::t_IndivTraits::t_Concentration functor.
  //! \param T_INDIVIDUAL is a \e physical individual
  //! \todo The back Fourier transform is done in (!) the t_Concentration
  //! object. Which means both this functor and the (at present) non-template
  //! t_Concentration classes need to know independantly about going real space
  //! to reciprocal space and back. Not quite logical. Probably only the
  //! template functor should have to know about this type of operation.
  template<class T_INDIVIDUAL> 
  class Krossover : public eoGenOp<T_INDIVIDUAL>
  {
    public:
      //! Type of a \e physical individual
      typedef T_INDIVIDUAL t_Individual; 
    protected:
      //! Contains all types pertaining to a \e physicla individual
      typedef typename t_Individual :: t_IndivTraits t_IndivTraits; 
      typedef typename t_IndivTraits :: t_Concentration t_Concentration; //!< Concentration type
      typedef typename t_IndivTraits :: t_FourierRtoK t_FourierRtoK; //!< Direct Fourier Functor
      typedef typename t_IndivTraits :: t_FourierKtoR t_FourierKtoR; //!< Back Fourier Functor
      typedef typename t_IndivTraits :: t_Object t_Object; //!< Object type
      typedef eoGenOp<t_Individual> t_Base; //!< Base class

    protected:
      t_Concentration &concentration; //!< reference to a t_Concentration object
      Ising_CE::Structure &structure; //!< reference to a structure
      types::t_real rate; //!< crossover rate (0.5 by default)
      bool do_range; //!< Whether to a range or indifferent krossover

    public:
      using t_Base::operator();

    public:
      //! Constructor and Initializer
      Krossover   ( t_Concentration &_c, Ising_CE::Structure &_str )
                : concentration(_c), structure(_str), rate(0.5), do_range(false) {}
      //! Copy Constructor
      Krossover   ( const Krossover &_k )
                : concentration(_k.concentration), structure(_k.structure),
                  rate(_k.rate), do_range(_k.do_range) {}
      //! Destructor
      ~Krossover() {}

      //! Loads Krossover::rate and Krossover::do_range parameters from XML input
      bool Load( const TiXmlElement &_node );

      //! Eo required bullshit
      virtual std::string className() const { return "TwoSites::Krossover"; }
      //! Eo required bullshit
      unsigned max_production(void) { return 1; } 

      //! \brief Does the actual krossover
      //! \param  _indiv first parent on input, offspring on output
      //! \param _parent second parent
      bool operator()( t_Individual &_indiv, const t_Individual &_parent );

      //! \brief wrapper aroung Krossover::operator()() which uses eoPopulator
      //! to obtain two parents
      //! \note It seems that using a single line, as in 
      //! \code 
      // if ( operator()( *_pop, _pop.select() ) ) (*_pop).invalidate(); 
      //! \endcode  
      //! leads to a bug where the wrong individual is copied from the selector.
      //! See Revision 313 and earlier. 
      //! Hence, the current implementation uses a more explicit approach:
      //! \code
      //  t_Individual &offspring = *_pop;
      //  const t_Individual &parent = _pop.select();
      //  if ( operator()( offspring, parent ) ) (*_pop).invalidate(); 
      //!\endcode 
      void apply(eoPopulator<t_Individual>& _pop);

      //! \brief prints out parameters
      std::string print_out() const;
  };


  //! \brief %Mutation in reciprocal space
  //! \details Changes the intensity a random number of
  //! <STRONG>k</STRONG>-vectors to random values. More specifically, each
  //! <STRONG>k</STRONG>-vectors is susceptible to be mutated with a
  //! mutation_rate of KMutation::rate. If a <STRONG>k</STRONG>-vectors is
  //! mutated, then its intensity is changed to a random complex number in ([-M,M],[-M,M]) where
  //! M is the square root of the maximum absolute intensity in reciprocal-space
  //! before any mutation. 
  //! The concentration after mutation is set (or not) using a
  //! \a T_INDIVIDUAL::t_IndivTraits::t_Concentration functor.
  //! \param T_INDIVIDUAL is a \e physical individual.
  //! \todo The back Fourier transform is done in (!) the t_Concentration
  //! object. Which means both this functor and the (at present) non-template
  //! t_Concentration classes need to know independantly about going real space
  //! to reciprocal space and back. Not quite logical. Probably only the
  //! template functor should have to know about this type of operation.
  template<class T_INDIVIDUAL>
  class KMutation : public eoGenOp<T_INDIVIDUAL>
  {
    public:
      //! Type of a \e physical individual
      typedef T_INDIVIDUAL t_Individual; 
    protected:
      //! Contains all types pertaining to a \e physicla individual
      typedef typename t_Individual :: t_IndivTraits t_IndivTraits; 
      typedef typename t_IndivTraits :: t_Concentration t_Concentration; //!< Concentration type
      typedef typename t_IndivTraits :: t_FourierRtoK t_FourierRtoK; //!< Direct Fourier Functor
      typedef typename t_IndivTraits :: t_FourierKtoR t_FourierKtoR; //!< Back Fourier Functor
      typedef typename t_IndivTraits :: t_Object t_Object; //!< Object type
      typedef eoGenOp<t_Individual> t_Base; //!< Base class

    protected:
      t_Concentration &concentration; //!< reference to a t_Concentration object
      Ising_CE::Structure &structure; //!< reference to a structure
      types::t_real rate; //!< crossover rate (0.5 by default)

    public:
      using t_Base::operator();

    public:
      //! Constructor and Initializer
      KMutation   ( t_Concentration &_c, Ising_CE::Structure &_str )
                : concentration(_c), structure(_str), rate(0.5) {}
      //! Copy Constructor
      KMutation   ( const KMutation &_k )
                : concentration(_k.concentration), structure(_k.structure),
                  rate(_k.rate) {}
      //! Destructor
      ~KMutation() {}

      //! Loads KMutation::rate parameter from XML input
      bool Load( const TiXmlElement &_node );

      //! Eo required bullshit
      virtual std::string className() const { return "TwoSites::KMutation"; }
      //! Eo required bullshit
      unsigned max_production(void) { return 1; } 

      //! \brief Does the actual KMutation
      //! \param _indiv individual to be mutated
      bool operator()( t_Individual &_indiv );
      //! \brief wrapper aroung KMutation::operator()() which uses eoPopulator
      //! to obtain an individual
      void apply(eoPopulator<t_Individual>& _pop)
        { if ( operator()( *_pop ) ) (*_pop).invalidate(); }
      //! \brief prints out parameters
      std::string print_out() const;
  };


  //! \brief creates a random individidual from random values in reciprocal-space
  //! \details The concentration is set (or not) using a
  //! \a T_INDIVIDUAL::t_IndivTraits::t_Concentration functor.
  //! \param T_INDIVIDUAL is a \e physical individual.
  //! \todo The back Fourier transform is done in (!) the t_Concentration
  //! object. Which means both this functor and the (at present) non-template
  //! t_Concentration classes need to know independantly about going real space
  //! to reciprocal space and back. Not quite logical. Probably only the
  //! template functor should have to know about this type of operation.
  template<class T_INDIVIDUAL> 
  class KRandom : public eoGenOp<T_INDIVIDUAL>
  {
    public:
      //! Type of a \e physical individual
      typedef T_INDIVIDUAL t_Individual; 
    protected:
      //! Contains all types pertaining to a \e physicla individual
      typedef typename t_Individual :: t_IndivTraits t_IndivTraits; 
      typedef typename t_IndivTraits :: t_Concentration t_Concentration; //!< Concentration type
      typedef typename t_IndivTraits :: t_Object t_Object; //!< Object type
      typedef eoGenOp<t_Individual> t_Base; //!< Base class

    protected:
      t_Concentration &concentration; //!< reference to a t_Concentration object
      Ising_CE::Structure &structure; //!< reference to a structure

    public:
      using t_Base::operator();

    public:
      //! Constructor and Initializer
      KRandom   ( t_Concentration &_c, Ising_CE::Structure &_str )
              : concentration(_c), structure(_str) {}
      //! Acts a functor. For convenience
      KRandom   ( t_Concentration &_c, Ising_CE::Structure &_str, t_Individual &_indiv )
              : concentration(_c), structure(_str) { operator()(_indiv); }
      //! Copy Constructor
      KRandom   ( const KRandom &_k )
              : concentration(_k.crossover), structure(_k.structure) {}
      //! Destructor
      ~KRandom() {}

      //! Doesn't load anything from anywhere and returns true
      bool Load( const TiXmlElement &_node ) { return true; }

      //! Eo required bullshit
      virtual std::string className() const { return "Darwin::Random"; }
      //! Eo required bullshit
      unsigned max_production(void) { return 1; } 

      //! \brief Creates a random bitstring individual
      //! \param _indiv individual to be initialized
      bool operator()( t_Individual &_indiv );

      //! \brief wrapper aroung KRandom::operator()() which uses eoPopulator
      //! to obtain two parents
      void apply(eoPopulator<t_Individual>& _pop)
        { if ( operator()( *_pop ) ) (*_pop).invalidate(); }

      //! \brief prints out parameters
      std::string print_out() const { return "Darwin::KRandom"; }
  };


  //! \brief Applies a standard bitstring crossover to a bitstring (of form \f$b_i=\pm1\f$)
  //! \details The concentration after krossover is set (or not) using a
  //! \a T_INDIVIDUAL::t_IndivTraits::t_Concentration functor.
  //! \param T_INDIVIDUAL is a \e physical individual.
  //! \related BitString::Krossover
  template<class T_INDIVIDUAL> 
  class Crossover : public eoGenOp<T_INDIVIDUAL>
  {
    public:
      //! Type of a \e physical individual
      typedef T_INDIVIDUAL t_Individual; 
    protected:
      //! Contains all types pertaining to a \e physicla individual
      typedef typename t_Individual :: t_IndivTraits t_IndivTraits; 
      typedef typename t_IndivTraits :: t_Concentration t_Concentration; //!< Concentration type
      typedef typename t_IndivTraits :: t_Object t_Object; //!< Object type
      typedef eoGenOp<t_Individual> t_Base; //!< Base class

    protected:
      t_Concentration &concentration; //!< reference to a t_Concentration object
      BitString::Crossover<t_Object> op; //!< the bitstring crossover functor

    public:
      using t_Base::operator();

    public:
      //! Constructor and Initializer
      Crossover   ( t_Concentration &_c )
             : concentration(_c), op() {}
      //! Copy Constructor
      Crossover   ( const Crossover &_k )
             : concentration(_k.crossover), op(_k.op) {}
      //! Destructor
      ~Crossover() {}

      //! Loads parameters from XML input
      bool Load( const TiXmlElement &_node ) { return op.Load( _node ); }

      //! Eo required bullshit
      virtual std::string className() const { return op.className(); }
      //! Eo required bullshit
      unsigned max_production(void) { return 1; } 

      //! \brief Does the actual crossover
      //! \param _indiv first parent on input, offspring on output
      //! \param _parent second parent
      bool operator()(t_Individual& _indiv, const t_Individual _parent );

      //! \brief wrapper aroung Crossover::operator()() which uses eoPopulator
      //! to obtain two parents
      //! \note It seems that using a single line, as in 
      //! \code 
      // if ( operator()( *_pop, _pop.select() ) ) (*_pop).invalidate(); 
      //! \endcode  
      //! leads to a bug where the wrong individual is copied from the selector.
      //! See Revision 313 and earlier. 
      //! Hence, the current implementation uses a more explicit approach:
      //! \code
      //  t_Individual &offspring = *_pop;
      //  const t_Individual &parent = _pop.select();
      //  if ( operator()( offspring, parent ) ) (*_pop).invalidate(); 
      //!\endcode 
      void apply(eoPopulator<t_Individual>& _pop);

      //! \brief prints out parameters
      std::string print_out() const;
  };


  //! \brief Applies a standard bitstring mutation to a bitstring (of form \f$b_i=\pm1\f$)
  //! \details The concentration after mutation is set (or not) using a
  //! \a T_INDIVIDUAL::t_IndivTraits::t_Concentration functor.
  //! \param T_INDIVIDUAL is a \e physical individual.
  //! \related BitString::Mutation
  template<class T_INDIVIDUAL> 
  class Mutation : public eoGenOp<T_INDIVIDUAL>
  {
    public:
      //! Type of a \e physical individual
      typedef T_INDIVIDUAL t_Individual; 
    protected:
      //! Contains all types pertaining to a \e physicla individual
      typedef typename t_Individual :: t_IndivTraits t_IndivTraits; 
      typedef typename t_IndivTraits :: t_Concentration t_Concentration; //!< Concentration type
      typedef typename t_IndivTraits :: t_Object t_Object; //!< Object type
      typedef eoGenOp<t_Individual> t_Base; //!< Base class

    protected:
      t_Concentration &concentration; //!< reference to a t_Concentration object
      BitString::Mutation<t_Object> op; //!< the bitstring crossover functor

    public:
      using t_Base::operator();

    public:
      //! Constructor and Initializer
      Mutation   ( t_Concentration &_c )
             : concentration(_c), op() {}
      //! Copy Constructor
      Mutation   ( const Mutation &_k )
             : concentration(_k.crossover), op(_k.op) {}
      //! Destructor
      ~Mutation() {}

      //! Loads parameters from XML input
      bool Load( const TiXmlElement &_node ) { return op.Load( _node ); }

      //! Eo required bullshit
      virtual std::string className() const { return op.className(); }
      //! Eo required bullshit
      unsigned max_production(void) { return 1; } 

      //! \brief Does the actual crossover
      //! \param _indiv individual on which the mutation is performed
      bool operator()( t_Individual &_indiv );

      //! \brief wrapper aroung Mutation::operator()() which uses eoPopulator
      //! to obtain an individual
      void apply(eoPopulator<t_Individual>& _pop)
        { if ( operator()( *_pop ) ) (*_pop).invalidate(); }

      //! \brief prints out parameters
      std::string print_out() const;
  };
  

  //! \brief Creates a random bitstring of \f$b_i=\pm1\f$ 
  //! \details The concentration after random generation is set (or not) using a
  //! \a T_INDIVIDUAL::t_IndivTraits::t_Concentration functor. This functor needs to know
  //! about the structure in order to create a bitstring of correct length.
  //! \param T_INDIVIDUAL is a \e physical individual.
  //! \todo Remove Random::structure? should bitstring length become an
  //! adjustable parameter?
  template<class T_INDIVIDUAL> 
  class Random : public eoGenOp<T_INDIVIDUAL>
  {
    public:
      //! Type of a \e physical individual
      typedef T_INDIVIDUAL t_Individual; 
    protected:
      //! Contains all types pertaining to a \e physicla individual
      typedef typename t_Individual :: t_IndivTraits t_IndivTraits; 
      typedef typename t_IndivTraits :: t_Concentration t_Concentration; //!< Concentration type
      typedef typename t_IndivTraits :: t_Object t_Object; //!< Object type
      typedef eoGenOp<t_Individual> t_Base; //!< Base class

    protected:
      t_Concentration &concentration; //!< reference to a t_Concentration object
      Ising_CE::Structure &structure; //!< reference to a structure

    public:
      using t_Base::operator();

    public:
      //! Constructor and Initializer
      Random   ( t_Concentration &_c, Ising_CE::Structure &_str )
             : concentration(_c), structure(_str) {}
      //! Acts a functor. For convenience
      Random   ( t_Concentration &_c, Ising_CE::Structure &_str, t_Individual &_indiv )
             : concentration(_c), structure(_str) { operator()(_indiv); }
      //! Copy Constructor
      Random   ( const Random &_k )
             : concentration(_k.crossover), structure(_k.structure) {}
      //! Destructor
      ~Random() {}

      //! Doesn't load anything from anywhere and returns true
      bool Load( const TiXmlElement &_node ) { return true; }

      //! Eo required bullshit
      virtual std::string className() const { return "Darwin::Random"; }
      //! Eo required bullshit
      unsigned max_production(void) { return 1; } 

      //! \brief Creates a random bitstring individual
      //! \param _indiv individual to be initialized
      bool operator()( t_Individual &_indiv );

      //! \brief wrapper aroung KRandom::operator()() which uses eoPopulator
      //! to obtain two parents
      void apply(eoPopulator<t_Individual>& _pop)
        { if ( operator()( *_pop ) ) (*_pop).invalidate(); }

      //! \brief prints out parameters
      std::string print_out() const { return "Darwin::Random"; }
  };
  

  //! \brief Creates GA operators for \e physical individuals from XML input
  //! \details Can create at present anyone of the following operators:
  //! - Krossover
  //! - KMutation
  //! - KRandom
  //! - Crossover
  //! - Mutation
  //! - Random
  //! . 
  //! \param _el XML node from which to load
  //! \param _structure needed by some functors 
  //! \param _concentration needed by some functors 
  //! \return a pointer to a eoGenOp object, or NULL if could not load from input
  //! \warning Checking the returned pointer, as well as storing and destroying it, is
  //! the responsability of the caller, not the callee.
  template<class T_INDIVIDUAL>
    eoGenOp<T_INDIVIDUAL>*
      LoadGaOp(const TiXmlElement &_el, Ising_CE::Structure &_structure, 
               typename T_INDIVIDUAL :: t_IndivTraits :: t_Concentration &_concentration );


  //! \brief Create a functor which returns false if the individual is outside
  //! a given concentration range
  //! \details This functor can be added to a GA::Taboos object. It can be used
  //! as any other taboo.
  //! \param T_INDIVIDUAL is a physical individual
  template< class T_INDIVIDUAL >
  class xTaboo : public Taboo_Base<T_INDIVIDUAL>
  {
    public:
      //! Type of a \e physical individual
      typedef T_INDIVIDUAL t_Individual; 
    protected:
      //! Contains all types pertaining to a \e physical individual
      typedef typename t_Individual :: t_IndivTraits t_IndivTraits; 
      typedef typename t_IndivTraits :: t_Concentration t_Concentration; //!< Concentration type
      typedef typename t_IndivTraits :: t_Object t_Object; //!< Object type

    protected:
      t_Concentration &concentration; //!< reference to a t_Concentration object
      types::t_real morethan; //!< lower limit of allowed range
      types::t_real lessthan; //!< upper limit of allowed range

    public:
      //! Constructor and Initializer 
       xTaboo   ( t_Concentration &_c) 
              : concentration(_c), morethan(-1.0), lessthan(-1.0) {}
      //! Copy Constructor
       xTaboo   ( const xTaboo &_c)
              : concentration(_c.concentration), morethan(_c.morethan),
                lessthan(_c.lessthan) {}
      //! Destructor
      ~xTaboo() {}

      //! Loads parameters from XML input
      bool Load( const TiXmlElement &_el );

      //! Eo required bullshit
      std::string className() const { return "Darwin::xTaboo"; }

      //! \brief returns true if \a _indiv is within allowed concentration range
      bool operator()( const t_Individual& _indiv );
  };


  //! \brief Creates a functor which writes to file an individual as an
  //!        xyz animation.
  //! \details It is expected that
  //!          operator<<( Ising_CE::Structure&, const t_Object )
  //!          is appropriately overloaded.
  template< class T_INDIVIDUAL >
  class XYZAnim : public eoMonOp<const T_INDIVIDUAL>
  {
    public:
      //! Type of a \e physical individual
      typedef T_INDIVIDUAL t_Individual; 
    protected:
      //! Contains all types pertaining to a \e physical individual
      typedef typename t_Individual :: t_IndivTraits t_IndivTraits; 
       //! Object type.
      typedef typename t_IndivTraits :: t_Object t_Object;
      //! The type of the eoMonOp base
      typedef eoMonOp<const t_Individual> t_OpBase;

    protected:
      //! \brief Pointer to a file for printing.
      //! \details Is a pointer for easy copy of the functor.
      Print::StdOut *file;
      //! Whether XYZAnim::file should be deleted by destructor.
      bool owns_pointer;
      //! Reference to a structure. 
      Ising_CE::Structure &structure;
      //! The number of calls made to XYZAnim::operator()().
      types::t_unsigned n;

    public:
      //! Constructor.
      XYZAnim   (Ising_CE::Structure &_str)
              : t_OpBase(), file(NULL), owns_pointer(false),
                structure(_str), n(0) {}
      //! Copy Constructor.
      XYZAnim   ( const XYZAnim &_c )
              : t_OpBase(_c), file( _c.file ), owns_pointer( false),
                structure(_c.structure), n(_c.n) {}
      //! This class contains a virtual function from its eoMonOp base.
      virtual ~XYZAnim() { if( owns_pointer and file ) delete file; }

      //! Prints \e _indiv to the file
      bool operator()( const t_Individual & _indiv);

      //! Loads this object from XML.
      bool Load( const TiXmlElement &_node );

    protected:
      //! \brief Opens the file for printing.
      //! \details Checks if XYZAnim::file is assigned and creates it if not.
      //!          Then relays to Print::StdOut::is_open(). 
      bool open();

      //! \brief Sets the filename and begins the file
      //! \details Print::StdOut::init(const std::string& ) is hidden from the
      //!          outside, so this function is necessary.
      void init( const std::string &_f );

  };


}

#include "gaoperators.impl.h"
/*@}*/
#endif
