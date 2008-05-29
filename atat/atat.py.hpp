
   class_< Matrix3d<_TYPE_> >( _PYTHONNAME_(Matrix3d) )
     .def( self + Matrix3d<_TYPE_>() )
     .def( self * Vector3d<_TYPE_>() );

   class_< Vector3d<_TYPE_> >( _PYTHONNAME_(Vector3d) )
     .def( self + Vector3d<_TYPE_>() )
     .def( Matrix3d<_TYPE_>() * self );

   def("to_attatvector", to_atatvector<_TYPE_>, return_value_policy<manage_new_object>() );
   def("to_attatmatrix", to_atatmatrix<_TYPE_>, return_value_policy<manage_new_object>() );

#undef _TYPE_
#undef _PYTHONNAME_
