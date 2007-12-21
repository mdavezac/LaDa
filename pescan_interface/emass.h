//
//  Version: $Id$
//
#ifndef _PESCAN_EMASS_H_
#define _PESCAN_EMASS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef _MPI
#include <mpi/mpi_object.h>
#endif

#include <lamarck/structure.h>
#include <atat/vectmac.h>
#include <physics/physics.h>

#include "interface.h"

namespace Pescan 
{

  //! electronic effective mass for ternary and quaternary super-lattices in (001) direction
  class eMassSL : public Interface
  {
    public:
      //! Inverse Effective Mass tensors
      atat::rMatrix3d inverse;
      //! Effective mass tensor;
      atat::rMatrix3d tensor;

    protected:
      //! Amplitude of the numerical derivative
      static const types::t_real amplitude = 0.002;
      //! \f$\frac{\sqrt(2)}{2}\f$
      static const types::t_real sqrt2_over2 = 0.707106781186547524400844362104849039284835937;
      //! \f$-\frac{\sqrt(2)}{2}\f$
      static const types::t_real msqrt2_over2 = -0.707106781186547524400844362104849039284835937;
      //! Kpoint Gamma
      static const atat::rVector3d Gamma;
      //! Kpoint \f$(100)\f$
      static const atat::rVector3d Lx;
      //! Kpoint \f$(010)\f$
      static const atat::rVector3d Ly;
      //! Kpoint \f$(001)\f$
      static const atat::rVector3d Lz;
      //! Kpoint \f$(110)\f$ normalized
      static const atat::rVector3d Hxy;
      //! Kpoint \f$(1\bar{1}0)\f$
      static const atat::rVector3d Hxmy;
      //! Kpoint \f$(011)\f$
      static const atat::rVector3d Hyz;
      //! Kpoint \f$(01\bar{1})\f$
      static const atat::rVector3d Hymz;
      //! Kpoint \f$(101)\f$
      static const atat::rVector3d Hxz;
      //! Kpoint \f$(\bar{1}01)\f$
      static const atat::rVector3d Hmxz;

    protected:
      //! Eigenvalues at gamma (degenerate)
      types::t_real eig_Gamma;
      //! Eigenvalues at \f$\Gamma + \epsilon(100)\f$
      types::t_real eig_Lx[2];
      //! Eigenvalues at \f$\Gamma + \epsilon(010)\f$
      types::t_real eig_Ly[2];
      //! Eigenvalues at \f$\Gamma + \epsilon(001)\f$
      types::t_real eig_Lz[2];
      //! Eigenvalues at \f$\Gamma + \epsilon(110)\f$
      types::t_real eig_Hxy[2];
      //! Eigenvalues at \f$\Gamma + \epsilon(1\bar{1}0)\f$
      types::t_real eig_Hxmy[2];
      //! Eigenvalues at \f$\Gamma + \epsilon(011)\f$
      types::t_real eig_Hyz[2];
      //! Eigenvalues at \f$\Gamma + \epsilon(01\bar{1})\f$
      types::t_real eig_Hymz[2];
      //! Eigenvalues at \f$\Gamma + \epsilon(101)\f$
      types::t_real eig_Hxz[2];
      //! Eigenvalues at \f$\Gamma + \epsilon(\bar{1}01)\f$
      types::t_real eig_Hmxz[2];
      

    public:
      //! Constructor
      eMassSL() { tensor.zero(); escan.kpoint = Gamma; }
      //! Copy Constructor
      eMassSL( const eMassSL & _c ) : tensor( _c.tensor ) {}
      //! Destructor
      ~eMassSL() {};


      //! Launches a calculation for structure \a _str
      bool operator()( const Ising_CE::Structure &_str ); 
      //! Loads Functional from XML
      bool Load( const TiXmlElement &_node );


    protected:
      //! Folded Spectrum gamma computation
      types::t_real gamma_folded_spectrum();
      //! All Electron gamma computation
      types::t_real gamma_all_electron( const Ising_CE :: Structure& );
      //! Kpoints other than Gamm
      void other_kpoints( const atat::rVector3d &_v, types :: t_real *_eig );
      //! Regroups all escan computations
      void escan_calls( const Ising_CE::Structure &_str );
  };

  inline void eMassSL :: other_kpoints( const atat::rVector3d &_v, 
                                        types :: t_real *_eig )
  {
    escan.kpoint = _v;
    launch_pescan();
    read_result();
    *_eig     = eigenvalues.front();
    *(_eig+1) = eigenvalues.back();
  }

  inline types::t_real eMassSL::gamma_folded_spectrum()
  {
    escan.kpoint = Gamma;
    escan.nbstates = 1;
    launch_pescan();
    read_result();

    return eigenvalues.back();
  }
  

} // namespace Pescan

#endif // _PESCAN_EMASS_H_
