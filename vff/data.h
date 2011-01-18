#ifndef LADA_VFF_DATA_H_
#define LADA_VFF_DATA_H_

#include "LaDaConfig.h"

#include <opt/types.h>
#include <opt/debug.h>
#include "exceptions.h"

namespace LaDa
{
  namespace vff
  {
    //! Returns name of a bond.
    inline std::string bond_type(std::string const &_A, std::string const &_B)
    {
      LADA_BASSERT(not _A.empty(), input() << error_string("Empty string for atomic specie."));
      LADA_BASSERT(not _B.empty(), input() << error_string("Empty string for atomic specie."));
      return _A[0] < _B[0] ? _A + _B: _B + _A;
    }
    
    //! Returns name of an angle.
    inline std::string angle_type(std::string const &_A, std::string const &_B, std::string const &_C)
    {
      LADA_BASSERT(not _A.empty(), input() << error_string("Empty string for atomic specie."));
      LADA_BASSERT(not _B.empty(), input() << error_string("Empty string for atomic specie."));
      LADA_BASSERT(not _C.empty(), input() << error_string("Empty string for atomic specie."));
      if(_A[0] < _B[0])
      {
        if(_A[0] < _C[0]) return _B[0] < _C[0] ? _A + _B + _C: _A + _C + _B;
        return _C + _A + _B;
      }
      else if(_B[0] < _C[0]) return _A[0] < _C[0] ? _B + _A + _C: _B + _C + _A;
      return _C + _B + _A;
    }
    
    size_t const max_vff_expansion = 6;
    
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
