#if LADA_MATH_MODULE != 1
  namespace LaDa
  {
    namespace error
    {
      //! Root of math errors.
      struct math : virtual root {};
      //! Thrown two arrays have different sizes.
      struct array_of_different_sizes: virtual math {};
      //! \brief Thrown when an structure is not the expected supercell of a lattice.
      //! \details This error may occur in methods which expects structures to be
      //!          supercells of a lattice, without allowance for relaxation. 
      struct unideal_lattice: virtual math, virtual root {};
      //! Thrown when an atomic position is unexpectedly off-lattice.
      struct off_lattice_position: virtual math, virtual unideal_lattice {};
      //! Thrown when a structure is not a supercell of an ideal lattice.
      struct not_a_supercell: virtual math, virtual unideal_lattice {};
      //! Thrown when a matrix is singular.
      struct singular_matrix: virtual math, virtual input {};
      //! Thrown when a matrix has a negative determinant.
      struct negative_volume: virtual math, virtual input {};
    }
  }
# endif
  
  
  
  
  
  
  
  
  
  
  
  
