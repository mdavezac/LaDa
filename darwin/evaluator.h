//
//  Version: $Id$
//
#ifndef _DARWIN_EVALUATOR_H_
#define _DARWIN_EVALUATOR_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <list>
#include <string>
#include <algorithm>

#include <eo/eoGenOp.h>

#include <opt/types.h>
#include <opt/function_base.h>
#include <mpi/mpi_object.h>

#include "loadsave.h"

#include "gatraits.h"

namespace LaDa
{
  namespace GA 
  {
    //! \brief base class for future functionals. See \ref secnew. 
    //! \details Implements the expected behavior of an evaluator class, eg an
    //! application-specific class capable of returning application-specific mating
    //! operators, application-specific evaluations... Note that only one of the
    //! behaviors implemented below is virtual. Indeed, it is expected that they
    //! will be overridden explicitely.
    //!
    //! Evaluator contains all behaviors required by %GA for an evaluator class,
    //! except for one, a traits class describing all %GA %types.
    //! The implementation of an actual could start with defining an individual.
    //! One easy way is to use the convenience %types afforded by Individual::Types,
    //! \code
    // typedef Individual::Types< SingleSite::Object, 
    //                            SingleSite::Concentration, 
    //                            SingleSite::Fourier        > :: t_Scalar t_Individual;
    //!  \endcode 
    //! You can use whatever else you want, as long as the required behaviors are there.
    //! The next step is to derive a class from GA::Evaluator as follows,
    //! \code
    // class Evaluator : public GA::Evaluator< t_Individual >
    // {
    //   public:
    //     typedef t_Individual t_Individual;
    //! \endcode
    //! where the instanciation of GA::Evaluator is  explicitely given in the derived class. 
    //! Just as important, the full traits of %GA must be defined, say  using Traits::GA,
    //! \code
    //     typedef Traits::GA< Evaluator > t_GATraits;
    //  \endcode 
    //! Then, override all the member %functions of GA::Evaluator that you want
    //! overridden. Finally,
    //! \code
    //  };
    //! \endcode 
    //! close your class. A job well done!
    //! For an example, see CE::Evaluator::t_GATraits and Traits::GA
    //! \param T_INDIVIDUAL is the type of individual to be used in this GA. This
    //! type will be passed on to rest of the %GA through Evaluator :: t_Individual
    //! and the derived class' t_GATraits :: t_Individual.
    template< class T_INDIVIDUAL>
    class Evaluator __MPICODE( : public MPI_COMMDEC )
    {
      public:
        typedef T_INDIVIDUAL   t_Individual; //!< The type of individual
      protected:
        //! %Traits of the individual
        typedef typename t_Individual :: t_IndivTraits        t_IndivTraits;
        //! Genetic Object type
        typedef typename t_IndivTraits :: t_Object            t_Object;
        //! %Traits of the quantity
        typedef typename t_IndivTraits :: t_QuantityTraits    t_QuantityTraits;
        //! Type of the quantity
        typedef typename t_QuantityTraits :: t_Quantity       t_Quantity;
        //! Type of the scalar quantity derived from t_Quantity
        typedef typename t_QuantityTraits :: t_ScalarQuantity t_ScalarQuantity;
        //! %Traits for Lamarckian behaviors
        typedef typename t_IndivTraits :: t_VA_Traits         t_VATraits;
        //! Type of the Larmarckian gradients
        typedef typename t_VATraits :: t_QuantityGradients    t_QuantityGradients;
        //! Type of the Larmarckian variables
        typedef typename t_VATraits :: t_Type                 t_VA_Type;
 
      protected:
        //! \brief Individual currently being assessed. 
        //! \details Evaluator::current_individual is set in Evaluator::init (and
        //! nowhere else, please!). It is usefull to keep track of the individual
        //! when doing Lamarckian stuff (and hence mutliple calls by, say,
        //! MinimizerGenOp to evaluation members of your derived class),  so that
        //! initialization is done once only, at the beginning of the local search.
        t_Individual *current_individual; 
        //! \brief see Evaluator::current_individual
        //! \details Provided for convenience.
        t_Object *current_object;
        //! Suffix string for files/directories
        __MPICODE( std::string suffix; ) 
 
      public:
        //! Constructor
        Evaluator() : __MPICODE( MPI_COMMCOPY( *::LaDa::mpi::main ) __COMMA__ )
                      current_individual(NULL),
                      current_object(NULL)
                      __MPICODE( __COMMA__ suffix("") ) {};
        //! Copy Constructor
        Evaluator   ( const Evaluator &_c )
                  : __MPICODE( MPI_COMMCOPY( _c ) __COMMA__ )
                    current_individual( _c.current_individual ),
                    current_object( _c.current_object )
                    __MPICODE( __COMMA__ suffix( _c.suffix ) ) {};
        //! Destructor
        virtual ~Evaluator() {}
 
      public:
        //! \brief opens XML file \a _f and calls Evaluator::Load(const TiXmlElement &_node )
        //! \details Should load t_Individual and funtional related stuff, except
        //! for the attributes of the \<GA\> tag.
        bool Load ( std::string const &_f );
        //! \brief Loads from XML input
        //! \details Is made virtual so that the correct member %function is
        //! called in Evaluator::Load( std::string const &)
        virtual bool Load ( const TiXmlElement &_node ) { return true; }
        //! \brief  Loads an  individual
        //! \details Results are expected in GA::LOADSAVE_LONG format, and GA
        //! internal stuff in GA::LOADSAVE_SHORT. You need both only if the
        //! ga object and the user-expected result object are different (say
        //! bitstring versus a decorated lattice structure )
        //! \param _indiv is the individual to load. More explicitely, load the
        //! indiviudal's t_Object instance!
        //! \param _node The XML node to load from. Should be an \<Individual\> Tag...
        //! \param _type can be either GA::LOADSAVE_SHORT or
        //!              GA::LOADSAVE_LONG.
        bool Load ( t_Individual &_indiv, const TiXmlElement &_node, bool _type ) {return true;};
        //! \brief  Saves an individual
        //! \details Results are expected in GA::LOADSAVE_LONG format, and GA
        //! internal stuff in GA::LOADSAVE_SHORT. You need both only if the
        //! ga object and the user-expected result object are different (say
        //! bitstring versus a decorated lattice structure )
        //! \param _indiv is the individual to save. More explicitely, save the
        //! indiviudal's t_Object instance!
        //! \param _node The XML node to save to. Should be an \<Individual\> Tag...
        //! \param _type can be either GA::LOADSAVE_SHORT or
        //!              GA::LOADSAVE_LONG.
        bool Save ( const t_Individual &_indiv, TiXmlElement &_node, bool _type ) const {return true;};
        //! \brief attributes from \<GA\> tag in input.xml are passed to this
        //! %function from Darwin::Load(...)
        void LoadAttribute ( const TiXmlAttribute &_att ) {};
        //! \brief Loads application-specific Mating operators
        //! \return  a pointer to an eoOp object. This pointer is owned by
        //! GA::Darwin::eostates, so don't deallocate yourself.
        eoGenOp<t_Individual>* LoadGaOp(const TiXmlElement &_el ) { return NULL; };
        //! \brief Loads application-specific continuators
        //! \details For a description  of continuators, see the EO library.
        //! \return  a pointer to an eoF<bool> object. This pointer is owned by
        //! GA::Darwin::eostates, so don't deallocate yourself.
        eoF<bool>* LoadContinue(const TiXmlElement &_el ) { return NULL; }
        //! \brief Loads an application-specific Scaling
        //! \details Scalings depend on t_GATraits which is undefined at this
        //!          level. Hence this function returns a void pointer.
        //! \returns a Scaling object which will be owned by the callee.
        void* Load_Scaling( const TiXmlElement &_node ) { return NULL; }
        //! \brief Loads a niche with an application-specific Distance
        //! \details Niches depend on the type of Distance which is undefined at this
        //!          level. Hence this function returns a void pointer.
        //! \return a pointer to a niche object, eg Sharing::Triangular< T_DISTANCE >. 
        //!         This pointer is owned by the callee.
        void* Load_Niche( const TiXmlElement &_node ) { return NULL; }
        //! \brief Returns a functor specifying the type of printout todo for
        //!        best individual.
        //! \details This functor will be called at the end of each generation
        //!          and applied to the best stored individual.
        //! \return a pointer to an eoOp<const t_Individual>.
        //!         This pointer is owned by the callee.
        eoMonOp<const t_Individual>* LoadPrintBest( const TiXmlElement &_node ) { return NULL; }
        //! \brief Random initialization of the t_object instance of \a _indiv
        //! \return true if \a _indiv should be invalidated (eg, \a _indiv has changed)
        bool initialize( t_Individual &_indiv ) {return false; }; 
        //! \brief Called before objective %function is evaluated. 
        //! \details See implementation below.
        void init( t_Individual &_indiv );
        //! \brief Evaluate Evaluator::current_indiv and stores the results in
        //! its quantities.
        void evaluate() {};
        //! \brief Evaluates the gradient of Evaluator::current_individual
        //! \details Only needed for Lamarckian %GA
        void evaluate_gradient( t_QuantityGradients& _grad ) { Traits::zero_out( _grad ); }
        //! \brief Evaluates Evaluator::current_individual and its gradients
        //! \details Only needed for Lamarckian %GA
        void evaluate_with_gradient( t_QuantityGradients& _grad );
        //! \brief Evaluates the gradient of Evaluator::current_individual in direction \a _pos
        //! \details Only needed for Lamarckian %GA
        void evaluate_one_gradient( t_QuantityGradients& _grad, types::t_unsigned _pos) 
          { Traits :: zero_out( _grad[_pos] ); }
        //! \brief Submits individuals to history, etc, prior to starting %GA
        //! \details initializes the endopoints of a convex-hull, for instance.
        //! Presubmitted individuals are not put into the population.
        //! \note GA::Evaluator does not know about Traits::GA (or affiliated).
        //! Hence it does not know about Traits::GA::t_Population. So we input this
        //! population an eoPop<t_Individual> directly.
        void presubmit( std::list<t_Individual>& ) { return; }
        //! \brief Prints parameters.
        //! \details Any and all, but is meant mostly for GA parameters.
        std::string print() const{  return ""; }
 
#       ifdef _MPI
          //! Sets communicator and suffix for mpi stuff.
          void set_mpi( boost::mpi::communicator *_comm, const std::string &_suffix )
            { MPI_COMMDEC::set_mpi( _comm ); suffix = _suffix; }
          //! Allows derived classes to have access to ::mpi::AddCommunicator members. 
          boost::mpi::communicator &comm() { return MPI_COMMDEC::comm(); } 
          //! Allows derived classes to have access to ::mpi::AddCommunicator members. 
          const boost::mpi::communicator &comm() const { return MPI_COMMDEC::comm(); } 
#       endif
    };
 
    template< class T_INDIVIDUAL>
    inline bool Evaluator<T_INDIVIDUAL> :: Load ( std::string const &_f )
    {
      TiXmlDocument doc( _f.c_str() ); 
      TiXmlHandle docHandle( &doc ); 
      if  ( !doc.LoadFile() )
      { 
        std::cerr << doc.ErrorDesc() << std::endl; 
        throw "Could not load input file in CE::Evaluator ";
      } 
 
      TiXmlElement *child = docHandle.FirstChild("Job").Element();
      if (not child)
        return false;
      return Load(*child);
    }
 
    template< class T_INDIVIDUAL>
    inline void Evaluator<T_INDIVIDUAL> :: init( t_Individual &_indiv )
    {
      current_individual = &_indiv;
      current_object     = &_indiv.Object();
    };
 
    template< class T_INDIVIDUAL>
    inline void Evaluator<T_INDIVIDUAL> :: evaluate_with_gradient( t_QuantityGradients& _grad )
    {
      evaluate_gradient( _grad );
      evaluate();
    }
 
  }
} // namespace LaDa

#endif // _EVALUATOR_H_
