#include "LaDaConfig.h"
#include "FCMangle.h"

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/trim.hpp>

#include <physics/physics.h>

#include "ewald.h"

extern "C" void FC_GLOBAL( ewaldf, EWALDF )
                (
                  const int *const,    // verbosity
                  double * const,      // Energy
                  double * const,      // forces (reduced)
                  double * const,      // forces (cartesian)
                  const double *const, // stress
                  const int *const,    // number of atoms
                  const double *const, // reduced atomic coordinates.
                  const double *const, // atomic charges
                  const double *const, // real space cutoff
                  const double *const, // cell vectors
                  const int *const,    // dimension of arrays.
                  int * const          // ERROR
                );
namespace LaDa
{
  namespace Models
  {
    bool fill_poscharge( crystal::Structure const &_input, PyObject* _cmap,
                         math::rMatrix3d const &_invcell, 
                         std::vector<double> &_positions, 
                         std::vector<double> &_charges )
    {
      python::Object charge_str = PyString_FromString("charge");
      if(not charge_str) return NULL;
      crystal::Structure::const_iterator i_atom = _input.begin();
      crystal::Structure::const_iterator const i_atom_end = _input.end();
      _positions.reserve(_input.size()*3)
      _charges.reserve(_input.size())
      for(; i_atom != i_atom_end; ++i_atom)
      {
        PyObject *item = PyDict_GetItem(i_atom->dict(), charge_str.borrowed());
        if(item == NULL and cmap != NULL and i_atom->pytype() != NULL)
          item = PyDict_GetItem(_cmap, i_atom->pytype());
        if(item == NULL)
        {
          LADA_PYERROR(ValueError, "Could not find charge associated with specie.");
          return false;
        }
        types::t_real const charge = math::convert_toreal(item, "elementary_charge");
        if(std::abs(charge) < 1e-12 and PyErr_Occurred()) return false;
        _charges.push_back(charge);
        math::rVector3d const frac(_invcell * i_atom->pos());
        _positions.push_back(pos[0]);
        _positions.push_back(pos[1]);
        _positions.push_back(pos[2]);
      }
      return true;
    }

    PyObject ewald(PyObject *_module, PyObject* _args, PyObject* _kwargs)  
    {
      static char *kwlist[] = { const_cast<char*>("crystal"),
                                const_cast<char*>("cutoff"),
                                const_cast<char*>("map"), NULL };
      PyObject *_in = NULL;
      PyObject *_cutoff = NULL;
      PyObject *cmap = NULL;
      if(not PyArg_ParseTupleAndKeywords( _args, _kwargs, "O|OO:ewals",
                                          kwlist, &_in, &_cutoff, &cmap) )
        return NULL;
      if(not PyStructure_Check(_in))
      {
        LADA_PYERROR(TypeError, "First argument to ewal should be a structure.");
        return NULL;
      }
      crystal::Structure input = crystal::Structure::acquire(_in);

      double const cutoff( convert(_cutoff, "Ry", 5e0) );
      if(cutoff < 0e0) return NULL;

      if(cmap != NULL and cmap != Py_None and (not PyDict_Check(cmap)))
      {
        LADA_PYERROR(TypeError, "cmap should be None or a dictionary of charges.");
        return NULL;
      }

      std::vector<double> charges, positions, forces, cforces;
      double cell[ 9 ], stress[ 6 ];
      const int verbosity(0);
      double energy(0);
      int error;
      const int n( input.size() );

      if(not fill_poscharge(input, cmap, _input.cell().inverse(), positions, charges))
        return NULL;


      types::t_real scale = convert_toreal(input.scale(), "a0")
      for( size_t i(0); i < 3; ++i )
        for( size_t j(0); j < 3; ++j )
          cell[ i + j*3 ] = _in.cell(i,j) * _in.scale / Physics::a0("A");
      forces.resize(3*input.size());
      cforces.resize(3*input.size());
      std::fill(forces.begin(), forces.end(), 0e0);
      std::fill(cforces.begin(), cforces.end(), 0e0);
      std::fill(stress, stress+6, 0e0);

      FC_GLOBAL( ewaldf, EWALDF )
      (
        &verbosity,   // verbosity
        &energy,      // Energy
        forces,       // forces (reduced)
        cforces,      // forces (cartesian)
        stress,       // stress
        &n,           // number of atoms
        positions,    // reduced atomic coordinates.
        Qs,           // atomic charges
        &gcut,        // g-space cutoff in Ry.
        cell,         // cell vectors
        &n,           // dimension of arrays.
        &error
      );
      LADA_ASSERT(error != 1, "Could not find optimal alpha for ewald summation.")

      // Copy (reduced) forces.
      types::t_real ffac = math::Quantity(1e0, "Ry").get("eV")    
                           / math::Quantity(1e0, "a0").get("angstrom");
      crystal::Structure :: iterator i_atom = input.begin();
      crystal::Structure :: iterator const i_atom_end = input.end();
      std::vector<double>::const_iterator i_force = cforces.begin();
      int const nptype = math::numpy::type<types::t_real>::value;
      python::Object np_forces = math::numpy::create_array<nptype>(input.size(), 3);
      PyObject * const inp_forces = np_forces.borrowed()
      if(not np_forces) return NULL;
      npy_intp dims[1] = {3}
      for(size_t i(0); i_atom != i_atom_end; ++i_atom, ++i_force, ++i )
      {
        *((types::t_real*)PyArray_GETPTR2(inp_forces, i, 0)) = *i_force;
        *((types::t_real*)PyArray_GETPTR2(inp_forces, i, 1)) = *(++i_force);
        *((types::t_real*)PyArray_GETPTR2(inp_forces, i, 2)) = *(++i_force);
        python::Object npforce
          = PyArray_NewFromData(1, dims, nptype, PyArray_GETPTR2(inp_forces, i, 0));
        if(not npforce) return NULL;
        ((PyArrayObject*)npforce.borrowed())->base = inp_forces;
        Py_INCREF(inp_forces);
        if( PyDict_SetItemString(i_atom->dict(), "forces", npforce.borrowed()) != 0)
          return NULL;
      } // loop over atoms
      // copy stress.
      python::Object npstress = math::numpy::create_array<nptype>(3, 3);
      if(not npstress) return NULL;
#     ifdef LADA_MACRO
#       error LADA_MACRO already defined
#     endif
#     define LADA_MACRO(I, J, T)                                               \
        *((types::t_real*)PyArray_GETPTR2(stress.borrowed(), I, J))            \
            =  stress[t] * sfac;
      LADA_MACRO(0,0,0)
      LADA_MACRO(1,1,1)
      LADA_MACRO(2,2,2)
      LADA_MACRO(0,1,3)
      LADA_MACRO(1,0,3)
      LADA_MACRO(1,2,4)
      LADA_MACRO(2,1,4)
      LADA_MACRO(0,2,5)
      LADA_MACRO(2,0,5)
#     undef LADA_MACRO
      
      return energy * Physics::Rydberg("eV");
      
    }


    Ewald :: t_Return Ewald :: operator()( const t_Arg& _in, t_Arg &_out ) const
    {
      LADA_DO_NASSERT( _in.atoms.size() != _out.atoms.size(), "Incoherent structure size.\n" )

      const size_t natoms( _in.atoms.size() );
      double Qs[ natoms ];
      double positions[ natoms * 3 ], forces[ natoms * 3 ], cforces[ natoms * 3 ];
      double cell[ 9 ], stress[ 6 ];
      const int verbosity(0);
      double energy(0);
      int error;
      const int n( natoms );
      const double gcut( cutoff_ / Physics::Rydberg("eV") );
      math::rMatrix3d const inv(!_in.cell);

      typedef t_Arg :: t_Atoms :: const_iterator t_cit;
      t_cit i_atom = _in.atoms.begin();
      for( size_t i(0); i < natoms; ++i, ++i_atom )
      {
        LADA_NASSERT( charges.find( i_atom->type ) == charges.end(),
                  "Atomic charge does not exist.\n" )
        Qs[i] = charges.find(i_atom->type)->second;
        math::rVector3d const frac(inv * i_atom->pos);
        positions[ i*3 ]     = frac[0];
        positions[ i*3 + 1 ] = frac[1];
        positions[ i*3 + 2 ] = frac[2];
      } // loop over atoms

      types::t_real const scale = input.scale() / math::Quantity(1e0, "a0").get("angstrom");
      for( size_t i(0); i < 3; ++i )
        for( size_t j(0); j < 3; ++j )
          cell[ i + j*3 ] = _in.cell(i,j) * scale;

      FC_GLOBAL( ewaldf, EWALDF )
      (
        &verbosity,   // verbosity
        &energy,      // Energy
        forces,       // forces (reduced)
        cforces,      // forces (cartesian)
        stress,       // stress
        &n,           // number of atoms
        positions,    // reduced atomic coordinates.
        Qs,           // atomic charges
        &gcut,        // g-space cutoff in Ry.
        cell,         // cell vectors
        &n,           // dimension of arrays.
        &error
      );
      LADA_ASSERT(error != 1, "Could not find optimal alpha for ewald summation.")

      // Copy (reduced) forces.
      crystal::iterator i_force = output.begin();
      crystal::const_iterator i_force = output.end();
      for( size_t i(0); i < natoms; ++i, ++i_force )
      {
        i_force->pos[0] += cforces[ i*3 ]     * ffac;
        i_force->pos[1] += cforces[ i*3 + 1 ] * ffac;
        i_force->pos[2] += cforces[ i*3 + 2 ] * ffac;
      } // loop over atoms
      // copy stress.
      _out.scale = 1e0;
      _out.cell(0,0) += stress[0] * Physics::Rydberg("eV");
      _out.cell(1,1) += stress[1] * Physics::Rydberg("eV");
      _out.cell(2,2) += stress[2] * Physics::Rydberg("eV");
      _out.cell(0,1) += stress[3] * Physics::Rydberg("eV");
      _out.cell(1,0) += stress[3] * Physics::Rydberg("eV");
      _out.cell(1,2) += stress[4] * Physics::Rydberg("eV");
      _out.cell(2,1) += stress[4] * Physics::Rydberg("eV");
      _out.cell(0,2) += stress[5] * Physics::Rydberg("eV");
      _out.cell(2,0) += stress[5] * Physics::Rydberg("eV");
      
      return energy * Physics::Rydberg("eV");
    }

    bool Ewald :: Load( const TiXmlElement& _node )
    {
      const TiXmlElement* const parent = opt::find_node( _node, "Functional", "Ewald" );
      if( not parent ) return false;
      const TiXmlElement* child = parent->FirstChildElement( "Atom" );
      for(; child; child = child->NextSiblingElement("Atom") )
      {
        if( child->Attribute("type") and child->Attribute("Charge") ) continue;
        const std::string type = boost::algorithm::trim_copy(std::string(child->Attribute("type")));
        const types::t_real charge
          = boost::lexical_cast<types::t_real>( child->Attribute("charge") );
        LADA_DO_NASSERT( charges.find( type ) != charges.end(),
                    "Duplicate entry in Ewald functional for " + type + "\n" )
        charges[type] = charge;
      }

      return true;
    }

  } // namespace Models
} // namespace LaDa
