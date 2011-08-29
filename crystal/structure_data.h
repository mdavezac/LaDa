#ifndef LADA_CRYSTAL_STRUCTUREDATA_H
#define LADA_CRYSTAL_STRUCTUREDATA_H

#include "LaDaConfig.h"

#include <vector>
#include <ostream>
#include <boost/lambda/lambda.hpp>
#include <boost/serialization/serialization.hpp>


#include <math/set_cell.h>
#ifdef LADA_WITH_LNS
# include <load_n_save/xpr/push_back.h>
# include <load_n_save/action/fusion.h>
# include <load_n_save/action/vector.h>
#endif

#include "add_atom.h"
#include "atom.h"
#include "traits.h"

namespace LaDa 
{
  namespace crystal
  {
    //! Namespace for frozen structure enum.
    struct frozenstr
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
    template<class T_TYPE> struct StructureData 
    {
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW
      friend class boost::serialization::access;
#     ifdef LADA_WITH_LNS
        friend class load_n_save::access;
#     endif
      
      //! \typedef Type of the species
      typedef typename traits::StructureData<T_TYPE>::t_Type t_Type;
      //! \typedef The type of the collection of atoms. 
      typedef typename traits::StructureData<T_TYPE>::t_Atom t_Atom;
      //! \typedef The type of the collection of atoms. 
      typedef typename traits::StructureData<T_TYPE>::t_Atoms t_Atoms;

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
      StructureData() : name(""), energy(0), weight(1), scale(1), freeze(frozenstr::NONE) {};
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
      math::details::SetCell< boost::mpl::int_<1> >
        set_cell(types::t_real _x, types::t_real _y, types::t_real _z)
          { return math::details::SetCell< boost::mpl::int_<0> >(cell)(_x, _y, _z); }
      //! Initializer for cell.
      math::details::SetCell< boost::mpl::int_<1> >
        set_cell(math::rVector3d _pos)
          { return math::details::SetCell< boost::mpl::int_<0> >(cell)(_pos); }
      //! Returns nth atom.
      typename t_Atoms::reference operator[](size_t _n) { return atoms[_n]; }
      //! Returns nth atom.
      typename t_Atoms:: const_reference operator[](size_t _n) const { return atoms[_n]; }
      //! Access to cell parameters
      types::t_real operator()(size_t i, size_t j) const { return cell(i,j); }
      //! Access to cell parameters
      types::t_real& operator()(size_t i, size_t j) { return cell(i,j); }
    };

#   ifdef LADA_WITH_LNS
      template<class T_TYPE> template<class T_ARCHIVE>
        bool StructureData<T_TYPE>::lns_access(T_ARCHIVE &_ar, load_n_save::version_type const _version) 
        {
          namespace lns = LaDa :: load_n_save;
          namespace bf = boost::fusion;
          std::map<std::string, LaDa::types::t_unsigned> freeze_map;
          freeze_map["none"] = frozenstr::NONE;
          freeze_map["a0"]   = frozenstr::A0;
          freeze_map["a1"]   = frozenstr::A1;
          freeze_map["a2"]   = frozenstr::A2;
          freeze_map["all"]  = frozenstr::ALL;
#         ifdef LADA_TIE
#            error LADA_TIE already defined.
#         endif
#         define LADA_TIE(i) bf::tie(cell(i,0), cell(i,1), cell(i,2))
#         ifdef LADA_TOE
#            error LADA_TOE already defined.
#         endif
#         define LADA_TOE(i) bf::tie(cell(0,i), cell(1,i), cell(2,i))
          lns::xpr::Section const seccell = lns::section("Cell") 
            << (
                    (    lns::option("r0", lns::action=LADA_TIE(0))
                      && lns::option("r1", lns::action=LADA_TIE(1))
                      && lns::option("r2", lns::action=LADA_TIE(2)) )
                 || (    lns::option("a0", lns::action=LADA_TOE(0))
                      && lns::option("a1", lns::action=LADA_TOE(1))
                      && lns::option("a2", lns::action=LADA_TOE(2))  )
               );
#         undef LADA_TIE
#         undef LADA_TOE
          lns::xpr::Section const section =
            lns::section("Structure")  
              << lns::option("name", lns::action=name, lns::default_="")
              << lns::option("energy", lns::action=energy, lns::default_=0)
              << lns::option("weight", lns::action=weight, lns::default_=0)
              << lns::option("freeze", lns::action=lns::enum_(freeze, freeze_map),
                             lns::default_=frozenstr::NONE)
              << lns::option("scale", lns::action=scale, lns::default_=1e0)
              << ( ( seccell )  && lns::push_back(atoms) );
          return _ar & section;
        }
#   endif

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
