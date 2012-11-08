#ifndef LADA_VFF_DATA_H_
#define LADA_VFF_DATA_H_

#include "LaDaConfig.h"

#include <misc/types.h>
#include <opt/debug.h>
#include "exceptions.h"

namespace LaDa
{
  namespace vff
  {
    //! Returns name of a bond.
    inline std::string bond_type(std::string const &_A, std::string const &_B)
    {
      LADA_BASSERT(not _A.empty(), exceptions::bond_input()
                                      << exceptions::string("Empty string for atomic specie."));
      LADA_BASSERT(not _B.empty(), exceptions::bond_input()
                                      << exceptions::string("Empty string for atomic specie."));
      return _A[0] < _B[0] ? _A + _B: _B + _A;
    }
    
    //! Returns name of an angle.
    inline std::string angle_type(std::string const &_origin, std::string const &_A, std::string const &_B)
    {
      LADA_BASSERT( not _origin.empty(), 
                    exceptions::angle_input() << exceptions::string("Empty string for atomic specie."));
      LADA_BASSERT(not _A.empty(), exceptions::angle_input()
                                      << exceptions::string("Empty string for atomic specie."));
      LADA_BASSERT(not _B.empty(), exceptions::angle_input()
                                      << exceptions::string("Empty string for atomic specie."));
      return _A[0] < _B[0] ? _A + _origin + _B: _B + _origin +_A;
    }
    
    size_t const max_vff_expansion = 5;
    
    //! \brief Holds bond information (length and stretching parameters).
    //! \note This is a POD structure.
    struct BondData
    {
      //! Bond-length at equilibrium (Angstrom).
      types::t_real length;
      //! Bond-stretching parameters.
      types::t_real alphas[max_vff_expansion];
    };
    
    //! \brief Holds bond-angle and angle parameters.
    //! \note This is a POD structure.
    struct AngleData
    {
      //! sigma parameter.
      types::t_real sigma;
      //! \brief Bond-angle deformation parameter. 
      //! \details Taylor expansion term corresponding to first derivative versus
      //!          bond-length and angle. 
      types::t_real gamma;
      //! Angle deformation parameter.
      types::t_real betas[max_vff_expansion];
    };

  } // namespace vff 
} // namespace LaDa

#endif
