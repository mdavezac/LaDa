#ifndef LADA_CRYSTAL_CAST_H
#define LADA_CRYSTAL_CAST_H

#include "LaDaConfig.h"


#include "structure.h"


namespace LaDa 
{
  namespace crystal
  {
    //! \brief Casts a sequence atom to a string atom.
    //! \param[in] _in: Atom to cast.
    //! \param[inout] _in: resulting cast atom.
    //! \param[_in] _sep: Separator between atomic type in resulting atom.
    void cast( AtomData< std::vector<std::string> > const &_in, 
               AtomData< std::string> &_out, std::string const &_sep = ", ");
    //! \brief Casts a sequence atom to a string atom.
    //! \param[in] _in: Atom to cast.
    //! \param[_in] _sep: Separator between atomic type in resulting atom.
    Atom<std::string> cast(Atom< std::vector<std::string> > const &_in, std::string const &_sep = ", ");
    //! \brief Casts a string atom to a sequence atom.
    //! \param[in] _in: Atom to cast.
    //! \param[inout] _in: resulting cast atom.
    //! \param[_in] _sep: Separator between atomic type in resulting atom.
    void cast( AtomData<std::string> const &_in, 
               AtomData< std::vector<std::string> > &_out, std::string const &_sep = ",");
    //! \brief Casts a string atom to a sequence.
    //! \param[in] _in: Atom to cast.
    //! \param[_in] _sep: Different species will be split using this separator.
    Atom< std::vector<std::string> > cast(Atom<std::string> const &_in, std::string const &_sep = ",");
    //! \brief Casts a string structure to a sequence structure.
    //! \param[in] _in: Structure to cast.
    //! \param[inout] _in: resulting cast structure.
    //! \param[_in] _sep: Separator between atomic type in resulting structure.
    void cast( StructureData< std::vector<std::string> > const &_in, 
               StructureData<std::string> &_out, std::string const &_sep = ",");
    //! \brief Casts a sequence structure to a string structure.
    //! \param[in] _in: Structure to cast.
    //! \param[_in] _sep: Separator between atomic type in resulting structure.
    Structure<std::string> cast( Structure< std::vector<std::string> > const &_in,
                                 std::string const &_sep = ", " );
    //! \brief Casts a sequence structure to a structure structure.
    //! \param[in] _in: Structure to cast.
    //! \param[inout] _in: resulting cast structure.
    //! \param[_in] _sep: Different species will be split using this separator.
    void cast( StructureData<std::string> const &_in, 
               StructureData< std::vector<std::string> > &_out, std::string const &_sep = ",");
    //! \brief Casts a string structure to a sequence structure.
    //! \param[in] _in: Structure to cast.
    //! \param[_in] _sep: Different species will be split using this separator.
    Structure< std::vector<std::string> > cast(Structure<std::string> const &_in, std::string const &_sep = ", ");
  } // namespace crystal
} // namespace LaDa

#endif
