#ifndef _ISING_CE_STRUCTURE_H_
#define _ISING_CE_STRUCTURE_H_

#include "LaDaConfig.h"

#include <vector>
#include <ostream>
#include <boost/lambda/lambda.hpp>
#include <boost/serialization/serialization.hpp>


#include <opt/types.h>
#include <math/eigen.h>
#ifdef LADA_WITH_LNS
# include <load_n_save/lns.h>
#endif

#include "atom.h"
#include "add_atom.h"


namespace LaDa 
{
  namespace crystal
  {
    template<class TYPE> struct StructureData : 
       public details::call_add_atom< StructureData<TYPE> >
    {
      friend class boost::serialization::access;
#     ifdef LADA_WITH_LNS
        friend class load_n_save::access;
#     endif
      //! Namespace for frozen enum.
      struct frozen
      {
        enum type
        {
          NONE =   0, //!< Freeze no coordinate of the unit-cell
          XX   =   1, //!< Freeze (0,0)
          XY   =   2, //!< Freeze (0,1) 
          XZ   =   4, //!< Freeze (0,2)
          YX   =   8, //!< Freeze (0,1) 
          YY   =  16, //!< Freeze (1,1) 
          YZ   =  32, //!< Freeze (1,2) 
          ZX   =  64, //!< Freeze (2,2)
          ZY   = 128, //!< Freeze (2,2)
          ZZ   = 256, //!< Freeze (2,2)
          ALL  = 511, //!< Freeze all coordinates 
          A0   =  73, //!< Freeze all coordinates 
          A1   = 146, //!< Freeze all coordinates 
          A2   = 292, //!< Freeze all coordinates 
        };
      };

      //! The type of the collection of atoms. 
      typedef std::vector< Atom<TYPE> > t_Atoms;
      //! The type of the collection of atoms. 
      typedef Atom<TYPE> value_type;
      //! Type of the species
      typedef TYPE t_Type;

      //! The unit-cell of the structure in cartesian coordinate.
      math::rMatrix3d cell;
      //! Just an old variable with the number of the structure in those NREL PI files.
      std::string name;
      //! The energy of the structure, whatever (and if) it is computed with.
      types::t_real energy;
      //! Weight of structure in "some" set.
      types::t_real weight;
      //! The scale in which cartesian units are given.
      types::t_real scale;
      //! The frozen coordinates of the unit-cell.
      types::t_unsigned freeze;
      //! The atomic position in cartesian unit and their occupation.
      t_Atoms atoms;

      //! Constructor
      StructureData() : name(""), energy(0), weight(1), freeze(frozen::NONE) {};
      //! Copy Constructor
      StructureData   (const StructureData &_str)
                    : cell(_str.cell), atoms(_str.atoms), name(_str.name), 
                      energy(_str.energy), weight(_str.weight),
                      freeze(_str.freeze), scale( _str.scale ) {}
      //! Serializes a structure.
      template<class ARCHIVE> void serialize(ARCHIVE & _ar, const unsigned int _version);
#     ifdef LADA_WITH_LNS
        //! To load and save to xml-like input.
        template<class T_ARCHIVE>
          bool lns_access(T_ARCHIVE &_ar, load_n_save::version_type const _version);
#     endif
      //! Initializer for cell.
      details::SetCell< boost::mpl::int_<1> > set_cell(types::t_real _x, types::t_real _y, types::t_real _z)
        { return details::SetCell< boost::mpl::int_<0> >(cell)(_x, _y, _z); }
      //! Initializer for cell.
      details::SetCell< boost::mpl::int_<1> > set_cell(math::rVector3d _pos)
        { return details::SetCell< boost::mpl::int_<0> >(cell)(_pos); }
    };

#     ifdef LADA_WITH_LNS
        template<class TYPE> template<class T_ARCHIVE>
          bool StructureData<TYPE>::lns_access(T_ARCHIVE &_ar, load_n_save::version_type const _version) 
          {
            namespace lns = LaDa :: load_n_save;
            namespace bf = boost::fusion;
            std::map<std::string, LaDa::types::t_unsigned> freeze_map;
            freeze_map["none"] = frozen::NONE;
            freeze_map["a0"]   = frozen::A0;
            freeze_map["a1"]   = frozen::A1;
            freeze_map["a2"]   = frozen::A2;
            freeze_map["all"]  = frozen::ALL;
#           ifdef LADA_TIE
#              error LADA_TIE already defined.
#           endif
#           define LADA_TIE(i) bf::tie(cell(i,0), cell(i,1), cell(i,2))
#           ifdef LADA_TOE
#              error LADA_TOE already defined.
#           endif
#           define LADA_TOE(i) bf::tie(cell(0,i), cell(1,i), cell(2,i))
            lns::xpr::Section const seccell = lns::section("Cell") 
              << (
                      (    lns::option("r0", lns::tag=lns::required, lns::action=LADA_TIE(0))
                        && lns::option("r1", lns::tag=lns::required, lns::action=LADA_TIE(1))
                        && lns::option("r2", lns::tag=lns::required, lns::action=LADA_TIE(2)) )
                   || (    lns::option("a0", lns::tag=lns::required, lns::action=LADA_TOE(0))
                        && lns::option("a1", lns::tag=lns::required, lns::action=LADA_TOE(1))
                        && lns::option("a2", lns::tag=lns::required, lns::action=LADA_TOE(2))  )
                 );
#           undef LADA_TIE
#           undef LADA_TOE
            lns::xpr::Section const section =
              lns::section("Structure")  
                << ( seccell ); //  && lns::push_back(atoms) );
              // << lns::option("name", lns::action=name, lns::default_="")
              // << lns::option("energy", lns::action=energy, lns::default_=0)
              // << lns::option("weight", lns::action=weight, lns::default_=0)
              // << lns::option("freeze", lns::action=lns::enum_(freeze, freeze_map),
              //                lns::default_=FREEZE_NONE)
              // << lns::option("scale", lns::action=scale, lns::default_=1e0);
            return _ar & section;
          }
#     endif

    template< class T_TYPE >
      std::ostream& operator<<(std::ostream &_stream, StructureData<T_TYPE> const &_str)
      {
        namespace bl = boost::lambda;
        _stream << "\n Structure, scale: " << _str.scale << ", Volume: "
                << _str.cell.determinant()
                << ", Cell\n"
                << std::fixed << std::setprecision(5)
                << "   " << std::setw(9) << _str.cell(0,0)
                << "   " << std::setw(9) << _str.cell(0,1)
                << "   " << std::setw(9) << _str.cell(0,2) << "\n"
                << "   " << std::setw(9) << _str.cell(1,0)
                << "   " << std::setw(9) << _str.cell(1,1)
                << "   " << std::setw(9) << _str.cell(1,2) << "\n"
                << "   " << std::setw(9) << _str.cell(2,0)
                << "   " << std::setw(9) << _str.cell(2,1)
                << "   " << std::setw(9) << _str.cell(2,2) << "\n"
                << "\n   Atoms:\n";
                
        
        _stream << " Structure atoms\n";
        std::for_each( _str.atoms.begin(), _str.atoms.end(),
                       bl::var(_stream) << bl::_1 << "\n" );
        return _stream;
      }

    

  } // namespace Crystal

} // namespace LaDa

#endif
