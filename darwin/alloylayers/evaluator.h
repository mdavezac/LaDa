//
//  Version: $Id$
//
#ifndef _DARWIN_ALLOY_LAYERS_EVALUATOR_H_
#define _DARWIN_ALLOY_LAYERS_EVALUATOR_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/lexical_cast.hpp>
#include <string>
#include <ostream>

#include <eo/eoOp.h>

#include <tinyxml/tinyxml.h>

#include <crystal/structure.h>
#include <opt/types.h>
#include <mpi/mpi_object.h>
#include <vff/layered.h>

#include "../evaluator.h"
#include "../individual.h"
#include "../bitstring.h"
#include "../vff.h"
#include "../bandgap_stubs.h"
#include "../electric_dipole.h"

#include "policies.h"

namespace LaDa
{
  namespace GA
  {
    /** \ingroup Genetic
     * @{ */
    //! Optimization for superlattices of random alloy layers.
    namespace AlloyLayers
    {
      //! \brief Partially overrides GA::Evaluator default to implement behaviors
      //!        relevant to random-alloy layers. 
      //! \tparam T_INDIVIDUAL the type of the GA individual.
      //! \tparam T_TRANSLATE a policy to transforms a GA object into a structure.
      //!                     It should contain a void translate( const t_Object&,
      //!                     Crystal::Structure& ) member, a void translate(
      //!                     const Crystal::Structure&, t_Object& ) member, 
      //!                     a void translate( const std::string&, t_Object&)
      //!                     member, an a void translate( const t_Object&,
      //!                     std::string ).  These member must be
      //!                     public.  However, inheritance is further public
      //!                     interfaces. They allow transformation from objects
      //!                     to structure (for use in functional evaluations )
      //!                     and back, and from objects to strings (for
      //!                     serialization) and back. The string need not be
      //!                     human-readable, it is meant  to save the
      //!                     individual's genome in a compressed format. Results
      //!                     are serialized in a human readable manner by using
      //!                     the first two member functions.
      //! \tparam T_ASSIGN a policy which assigns functional values into GA
      //!                  quantities. E.g. sets Individual::quantities to
      //!                  whatever is needed. This extra layer of indirection
      //!                  should allow for us to input into the quantities
      //!                  different functional values as needed, such as
      //!                  bandgap, bandedges, strain...  In other words, it
      //!                  allows static and dynamic setting of multi-valued
      //!                  quantities. This policy must contain a public void
      //!                  assign( const t_Object&, t_Quantity ) member.
      //!                  Inheritance is public to allow further public
      //!                  interfaces.
      template< class T_INDIVIDUAL,
                template<class> class T_TRANSLATE,
                template<class,class> class T_ASSIGN >
      class EvaluatorBase : public GA::Evaluator< T_INDIVIDUAL >,
                            public T_TRANSLATE
                                     < typename T_INDIVIDUAL :: t_IndivTraits :: t_Object >,
                            public T_ASSIGN
                                     < typename T_INDIVIDUAL :: t_IndivTraits :: t_Object,
                                       typename T_INDIVIDUAL :: t_IndivTraits
                                                             :: t_QuantityTraits :: t_Quantity >
      {
        // Allow the policies to call private members of the derived class (eg this class).
        friend class T_TRANSLATE< typename T_INDIVIDUAL :: t_IndivTraits :: t_Object >;
        friend class T_ASSIGN
                     <
                       typename T_INDIVIDUAL :: t_IndivTraits :: t_Object,
                       typename T_INDIVIDUAL :: t_IndivTraits
                                             :: t_QuantityTraits :: t_Quantity
                     >;
        public:
          //! Type of the individual.
          typedef T_INDIVIDUAL t_Individual;
          //! Type of the translation policy
          typedef T_TRANSLATE< typename t_Individual :: t_IndivTraits :: t_Object >
                  t_Translate;
          //! Type of the Assignement policy
          typedef T_ASSIGN
                  < 
                    typename t_Individual :: t_IndivTraits :: t_Object,
                    typename t_Individual :: t_IndivTraits :: t_QuantityTraits :: t_Quantity
                  > t_Assign;
          //! All pertinent %GA traits
          typedef Traits::GA< EvaluatorBase > t_GATraits;

        protected:
          //! Traits of the individual.
          typedef typename t_Individual ::t_IndivTraits t_IndivTraits;
          //! Type of the %GA object.
          typedef typename t_IndivTraits :: t_Object t_Object;
          //! Evaluator base.
          typedef GA::Evaluator<t_Individual> t_Base;
          //! Type of this class.
          typedef EvaluatorBase<t_Individual, T_TRANSLATE, T_ASSIGN> t_This;

        protected:
          using t_Base :: current_individual;
          using t_Base :: current_object;
          using t_Translate :: translate;
          using t_Assign :: assign;

        public:
          using t_Base::Load;

        public:
          //! Default Constructor
          EvaluatorBase() {}
          //! Copy Constructor
          EvaluatorBase   ( const EvaluatorBase &_c )
                    : t_Base(_c), t_Translate(_c), t_Assign(_c),
                      structure( _c.structure), direction(_c.direction) {}
          //! Destructor
          virtual ~EvaluatorBase() {}

          //! \brief Loads the lattice, the epitaxial parameters from XML, and
          //!        constructs the structure.
          bool Load( const TiXmlElement &_node );
          //! Loads an individual from XML.
          bool Load( t_Individual &_indiv, const TiXmlElement &_node, bool _type );
          //! Saves an individual to XML.
          bool Save( const t_Individual &_indiv, 
                     TiXmlElement &_node, bool _type ) const;;

          //! Sets the current_individual pointer.
          void init( t_Individual &_indiv );

          //! Prints epitaxial parameters and structure
          std::string print() const;

          //! Creates random individuals using GA::Random.
          bool initialize( t_Individual &_indiv );

          //! Returns direction.
          const atat::rVector3d& get_direction() const { return direction; }

          //! Returns a constant reference to a structure.
          const Crystal::Structure& get_structure() const { return structure; }

        protected:
          //! The structure (cell-shape) for which decoration search is done
          mutable Crystal :: Structure structure;

          //! The growth direction.
          atat::rVector3d direction;
          //! A (static) pointer to a lattice object.
          static boost::shared_ptr< Crystal::Lattice > lattice;
      };

      //! An evaluator class for pescan/vff in epitaxial alloy configuration space.
      template< class T_INDIVIDUAL,
                template<class> class T_TRANSLATE,
                template<class,class> class T_ASSIGN >
        class Evaluator : public EvaluatorBase< T_INDIVIDUAL, T_TRANSLATE, T_ASSIGN >
        {
            //! Type of the base class.
            typedef EvaluatorBase< T_INDIVIDUAL, T_TRANSLATE, T_ASSIGN > t_Base;
          public:
            //! Type of the individual.
            typedef typename t_Base :: t_Individual t_Individual;
            //! Type of the translation policy.
            typedef typename t_Base :: t_Translate t_Translate;
            //! Type of the assignement policy.
            typedef typename t_Base :: t_Assign t_Assign;
            //! All pertinent %GA traits
            typedef Traits::GA< Evaluator > t_GATraits;
       
            //! Constructor
            Evaluator() : t_Base(), bandgap(structure),
                          edipole( structure ), do_dipole_(false) {}
            //! Copy Constructor
            Evaluator   ( const Evaluator &_c )
                      : t_Base(_c), bandgap(_c.bandgap),
                        edipole( structure ), do_dipole_(false) {}
            //! Destructor
            virtual ~Evaluator() {};
       
            //! Loads structure, lattice, bandgap, vff from XML
            bool Load( const TiXmlElement &_node );
       
            //! Initializes the bandgap.
            void init( t_Individual& _indiv ) { t_Base::init( _indiv ); bandgap.init(); }
       
            //! \brief Evaluates the band gap after strain minimization, and
            //!        optionally the electic dipole.
            void evaluate();

       
  #         ifdef _MPI
              //! forwards comm and suffix to bandgap.
              void set_mpi( boost::mpi::communicator *_comm, const std::string &_str )
                { bandgap.set_mpi( _comm, _str ); edipole.set_mpi( _comm, _str ); }
  #         endif
            using t_Base::Load;

            //! Returns true if dipole transition will be computed.
            bool do_dipole() const { return do_dipole; }
            //! Sets pescan to compute dipoles on \a _do true.
            void do_dipole( bool _do );
       
            //! Type of the band gap cum vff all-in-one functional.
            typedef BandGap::Darwin<Vff::Layered> t_BandGap;
            //! Type of the electric dipole.
            typedef OscStrength::Darwin t_ElectricDipole;
       
            //! BandGap/Vff functional
            t_BandGap bandgap; 
            //! Electric dipole functional.
            t_ElectricDipole edipole;
          protected:
            //! Wether to perform oscilator strength evaluations.
            bool do_dipole_;
            
            using t_Base :: structure;
        };



    } // namespace Layered
    /** @} */


 }
} // namespace LaDa
 

#include "evaluator.impl.h"

#endif // _LAYERED_H_
