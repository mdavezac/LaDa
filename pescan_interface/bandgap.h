//
//  Version: $Id$
//
#ifndef _PESCAN_BANDGAP_H_
#define _PESCAN_BANDGAP_H_


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <mpi/mpi_object.h>
#include <crystal/structure.h>
#include <print/stdout.h>

#include "interface.h"

namespace Pescan 
{

  //! Keeps track of the HOMO and LUMO
  struct Bands
  { 
    types::t_real cbm;  //!< Conduction Band Minimum
    types::t_real vbm;  //!< Valence Band minimum
    Bands() : cbm(0), vbm(0) {}; //!< Constructor
    //! Constructor and Initializer
    Bands( const types :: t_real _cbm, types::t_real _vbm ) : cbm(_cbm), vbm(_vbm) {};
    //! Copy Constructor.
    Bands( const Bands &_bands ) : cbm(_bands.cbm), vbm(_bands.vbm) {};
    //! Returns the gap.
    types::t_real gap() const { return cbm - vbm; }
  };

  //! Computes the band-gap using the pescan interface
  class BandGap : public Interface
  {
    protected:
      //! For folded spectrum, whcih band is being computed.
      enum t_computation
      { 
        VBM,  //!< Valence Band Maximum is being computed.
        CBM   //!< Conduction Band Minimum is being computed.
      };

    public:
      //! Stores results
      Bands bands;
      //! Reference energy for folded spectrum method
      Bands Eref;

    protected:
      //! For folded spectrum, which computation (VBM or CBM) is being performed.
      t_computation computation;
      //! Wether to apply correction scheme when metallic band-gap is found
      bool do_correct;
      //! Below this, a band-gap is assumed "metallic"
      const types::t_real metallicity;
      //! Amounts by which to increase/decrease references.
      const types::t_real inc_dec;

    public:
      //! Constructor
      BandGap  () 
             : bands(0,0), Eref(0,0), computation( CBM ), do_correct(true),
               metallicity(0.001), inc_dec(0.1) {}
      //! Copy Constructor
      BandGap   ( const BandGap & _c ) 
              : Interface( _c ), bands( _c.bands ), Eref( _c.Eref ),
                computation( CBM ), do_correct( _c.do_correct ),
                metallicity(0.001), inc_dec(0.1) {} 
      //! Destructor
      ~BandGap() {};


      //! \brief Loads all parameters from XML.
      //! \details Adds band-gap references
      bool Load( const TiXmlElement &_node );

      //! Launches a calculation for structure \a _str
      types::t_real operator()( const Crystal::Structure &_str ); 

      // Forwards mpi::AddCommunicator members.
      MPI_FORWARD_MEMBERS( Interface )

    protected:
      //! Folded spectrum computation
      types::t_real folded_spectrum(const Crystal::Structure& );
      //! All-electron spectrum computation
      types::t_real all_electron( const Crystal::Structure &_str );
      //! Returns the closest eigenvalue to \a _ref
      types::t_real find_closest_eig( types::t_real _ref );
      //! \brief Makes sure that a non-metallic bandgap has been obtained.
      //! \details If a metallic band-gap is found, then the eigenvalue
      //!          farthest from the fermi energy is moved away from it (its
      //!          value is decreased or increased depending on which side of
      //!          the fermi energy it is ) by BandGap::move_gap. After this, a
      //!          folded_spectrum is launched again, which itself will call
      //!          BandGap::correct. Once a true band-gap has been found, all
      //!          escan parameters are reset (including references).
      void correct( const std::string &_dir );
      //!  Read results + Throw error if eigenvalues cannot be found. 
      void read_and_throw();
  };

  //! \cond
  inline types::t_real BandGap :: operator()( const Crystal::Structure &_str )
  {
    set_scale( _str );

    return escan.method == ALL_ELECTRON ? all_electron( _str ):
                                          folded_spectrum( _str );
  }

  inline types::t_real BandGap :: find_closest_eig( types::t_real _ref )
  {
    if( eigenvalues.empty() ) return 0;
    std::vector<types::t_real> :: const_iterator i_eig = eigenvalues.begin();
    std::vector<types::t_real> :: const_iterator i_eig_end = eigenvalues.end();
    std::vector<types::t_real> :: const_iterator i_eig_result = i_eig;
    types::t_real mini = std::abs(*i_eig-_ref); 
    for(++i_eig; i_eig != i_eig_end; ++i_eig)
      if( std::abs( *i_eig - _ref ) < mini )
      {
        mini = std::abs( *i_eig - _ref );
        i_eig_result = i_eig;
      }
    return *i_eig_result;
  }
  //! \endcond

} // namespace Pescan

#endif // _PESCAN_BANDGAP_H_
