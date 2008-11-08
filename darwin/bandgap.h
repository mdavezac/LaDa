//
//  Version: $Id$
//
#ifndef _BANDGAP_H_
#define _BANDGAP_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>

#include <eo/eoOp.h>

#include <tinyxml/tinyxml.h>

#include <crystal/structure.h>
#include <opt/function_base.h>
#include <opt/gsl_minimizers.h>
#include <opt/types.h>
#include <mpi/mpi_object.h>

#include "two_sites.h"
#include "evaluator.h"
#include "concentration.h"
#include "functors.h"
#include "individual.h"
#include "bandgap_stubs.h"


namespace LaDa
{
  /** \ingroup Genetic
   * @{ */
  //! \brief Implements single-cell-shape decoration search for band gaps.
  //! \details The strain is relaxed using the Valence Force Field developped in
  //!          namespace Vff. The band-gap is obtained from pescan \e via the
  //!          Pescan::Interface class. 
  //! 
  //!          The band gap can be obtained from a complete all-electron
  //!          diagonalization or, more efficiently, through a partial
  //!          diagonalization using the folded spectrum method (see <A
  //!          HREF="http://dx.doi.org/10.1103/PhysRevB.51.17398"> L-W Wang and A.
  //!          Zunger PRB \b 51, 17398 (1995) </A>). The latter needs two
  //!          reference energies to do find the band-gap. If these are \b not
  //!          given on input, BandGap::Darwin will first perform a full
  //!          diagonalization and peg the reference energies to the HOMO and the
  //!          LUMO. From there, only Folded-Spectrum methods are performed
  //!          unless one of two situations arise: 
  //!            - A metallic "band gap" is found 
  //!            - The user has requested all-electron calculations every \a N
  //!              generations
  //!            .
  //! \see For more details on the functionals, see Pescan::Interface, 
  //                Vff::Functional and Pescan::BandGap.
  namespace BandGap
  {

    //! \brief BitString Object with BandGap capacity
    //! \details Other than the bitstring itself and the ::GA::Keepers::BandGap
    //           variables, this object also stores the x and y concentrations of
    //           a quaternary. It overloads dumping an object to a stream.
    //! \see BandGap::operator<<( std::ostream &, const Object& )
    struct Object : public TwoSites::Object,
                    public ::GA::Keepers::BandGap,
                    public ::GA::Keepers::ConcTwo
    {
      friend class boost::serialization::access;
      //! The type of the BitString container
      typedef TwoSites::Object :: t_Container t_Container;

      //! Constructor
      Object() : TwoSites::Object(), ::GA::Keepers::BandGap(),
                 ::GA::Keepers::ConcTwo() {}
      //! Copy Constructor
      Object   (const Object &_c)
             : TwoSites::Object(_c), ::GA::Keepers::BandGap(_c),
               ::GA::Keepers::ConcTwo(_c) {};
      //! Loads from \a _node.
      bool Load( const TiXmlElement &_node )
        { return     ::GA::Keepers::BandGap::Load( _node ) 
                 and ::GA::Keepers::ConcTwo::Load( _node ); }
      //! Saves to \a _node.
      bool Save( TiXmlElement &_node ) const
        { return     ::GA::Keepers::BandGap::Save( _node ) 
                 and ::GA::Keepers::ConcTwo::Save( _node ); }
      //! Destructor
      virtual ~Object() {};
      private:
        //! Serializes a scalar individual.
        template<class Archive> void serialize(Archive & _ar, const unsigned int _version);
    };

    //! \brief Dumps a BandGap::Object to a stream.
    //! \details  This routine is used for printing results only, and never to
    //!           serialize in XML format.
    inline std::ostream& operator<<(std::ostream &_stream, const Object &_o)
    { 
      if( _o.Container().size() <= 30 )
        _stream << (const SingleSite::Object& ) _o << " ";
      _stream << "x=" << _o.x << " y="  << _o.y 
              << (const ::GA::Keepers::BandGap&) _o;
      return _stream; 
    } 

    //! \brief Type of the \e physical BandGap individual
    //! \details In addition to BandGap object, the individual uses
    //!          TwoSites::Concentration and TwoSites Fourier. By default, it
    //!          declares a \e scalar fitness.
    typedef Individual::Types< BandGap::Object, 
                               TwoSites::Concentration, 
                               TwoSites::Fourier > :: t_Scalar t_Individual;

    //! \brief %Evaluator class for band-gap decoration search
    //! \details A BandGap::Darwin and a Vff::Darwin<Vff::Functional> objects are
    //!          declared which allow for strain and bandgap computation.
    //!          The \e physical ga operators of gaoperators.h are used to mate
    //!          individuals. \e A \e priori, I can't see what we would need to
    //!          compute a convex-hull for, so GA::Evaluator::presubmit() is
    //!          overriden to simply clear the population in its argument.
    class Evaluator : public TwoSites::Evaluator< BandGap::t_Individual >
    {
      public:
        //! \brief The type of the \e physical individual used in this decoration
        //! search.
        typedef BandGap::t_Individual t_Individual;
        //! All pertinent %GA traits
        typedef Traits::GA< Evaluator > t_GATraits;
      protected:
        //! \cond
        typedef Evaluator t_This;
        typedef TwoSites::Evaluator< t_Individual > t_Base;
        typedef Crystal::Structure::t_kAtoms t_kvecs;
        typedef Crystal::Structure::t_Atoms t_rvecs;
        //! \endcond

        //! Type of the band gap cum vff all-in-one functional.
        typedef BandGap::Darwin<Vff::Functional> t_BandGap;

      public:
        using t_Base :: Load;
        using t_Base :: Save;
      protected:
        using t_Base :: current_individual;
        using t_Base :: current_object;

      protected:
        //! BandGap/Vff functional
        t_BandGap bandgap; 

      public:
        //! Constructor
        Evaluator() : t_Base(), bandgap(structure) {}
        //! Copy Constructor
        Evaluator   ( const Evaluator &_c )
                  : t_Base(_c), bandgap(_c.bandgap) {}
        //! Destructor
        ~Evaluator() {};

        //! Saves an individual to XML
        bool Save ( const t_Individual &_indiv, TiXmlElement &_node, bool _type ) const;
        //! Loads an individual from XML
        bool Load ( t_Individual &_indiv, const TiXmlElement &_node, bool _type );
        //! Loads structure, lattice, bandgap, vff from XML
        bool Load( const TiXmlElement &_node );
        //! Allows bandgap all-electron recomputation.
        eoF<bool>* LoadContinue(const TiXmlElement &_el )
          { return new GA::mem_zerop_t<t_BandGap>( bandgap, &t_BandGap::Continue,
                                                        "BandGap::Continue" );     }
       
        //! \brief Intializes before calls to evaluation member routines
        //! \details The bandgap does not need explicit initialization, since it
        //!          will act upon the structure as minimized by vff. More
        //!          explicitely, its "initialization" is carried out in the body
        //!          of Darwin::evaluate().
        void init( t_Individual &_indiv )
          { t_Base :: init( _indiv ); bandgap.init(); }

        //! Evaluates the band gap after strain minimization
        void evaluate();

        //! \brief Pescan is costly and I'm not sure why we would want to do a
        //!        convex-hull, so presubmit does nothing in this implementation.
        //! Presubmitted individuals are not put into the population.
        //! \see GA::Evaluator::presubmit(), TwoSites::Evaluator::presubmit()
        void presubmit( std::list<t_Individual> &_pop ) { _pop.clear(); }
  #ifdef _MPI
        //! forwards comm and suffix to bandgap.
        void set_mpi( boost::mpi::communicator *_comm, const std::string &_str )
          { t_Base::set_mpi( _comm, _str ); bandgap.set_mpi( _comm, _str ); }
  #endif
    };
    
    template<class Archive>
      void Object :: serialize(Archive & _ar, const unsigned int _version)
      {
        _ar & boost::serialization::base_object< TwoSites::Object >( *this ); 
        _ar & boost::serialization::base_object< ::GA::Keepers::BandGap >( *this ); 
        _ar & boost::serialization::base_object< ::GA::Keepers::ConcTwo >( *this ); 
      }
    
  } // namespace BandGap
} // namespace LaDa

/** @} */


#endif // _BANDGAP_H_
