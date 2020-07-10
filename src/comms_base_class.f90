Module comms_base_class_module

  Use constants, Only : wp

  Implicit None

  Type, Public, Abstract :: comms_base_class
     Integer :: communicator
   Contains
     Procedure( reduce_real    ), Deferred, Public :: reduce
     Procedure( reduce_real    ), Deferred, Public :: max
     Procedure( set_comm       ), Deferred, Public :: set_comm
     Procedure( get_comm       ), Deferred, Public :: get_comm
     Procedure( inquiry_1d     ), Deferred, Public :: get_rank
     Procedure( inquiry_1d     ), Deferred, Public :: get_size
     Procedure( inquiry_3d     ), Deferred, Public :: get_proc_coords
     Procedure( inquiry_3d     ), Deferred, Public :: get_proc_grid
  End type comms_base_class

  Abstract Interface

     Subroutine set_comm( c, comm, error )
       Import :: comms_base_class
       Class( comms_base_class ), Intent(   Out ) :: c
       Integer,                   Intent( In    ) :: comm
       Integer,                   Intent(   Out ) :: error
     End Subroutine set_comm

     Subroutine get_comm( c, comm )
       Import :: comms_base_class
       Class( comms_base_class ), Intent( In    ) :: c
       Integer,                   Intent(   Out ) :: comm
     End Subroutine get_comm

     Subroutine inquiry_1d( c, data )
       Import :: comms_base_class
       Class( comms_base_class ), Intent( In    ) :: c
       Integer,                   Intent(   Out ) :: data
     End Subroutine inquiry_1d

     Subroutine inquiry_3d( c, data )
       Import :: comms_base_class
       Class( comms_base_class ), Intent( In    ) :: c
       Integer, Dimension( : ),   Intent(   Out ) :: data
     End Subroutine inquiry_3d

     Subroutine reduce_int( c, data )
       Import :: comms_base_class
       Class( comms_base_class ), Intent( In    ) :: c
       Integer,                   Intent( InOut ) :: data
     End Subroutine reduce_int

     Subroutine reduce_int_1d( c, data )
       Import :: wp
       Import :: comms_base_class
       Class( comms_base_class ),  Intent( In    ) :: c
       Integer,    Dimension( : ), Intent( InOut ) :: data
     End Subroutine reduce_int_1d

     Subroutine reduce_real( c, data )
       Import :: wp
       Import :: comms_base_class
       Class( comms_base_class ), Intent( In    ) :: c
       Real( wp ),                Intent( InOut ) :: data
     End Subroutine reduce_real

     Subroutine reduce_real_1d( c, data )
       Import :: wp
       Import :: comms_base_class
       Class( comms_base_class ),  Intent( In    ) :: c
       Real( wp ), Dimension( : ), Intent( InOut ) :: data
     End Subroutine reduce_real_1d

  End Interface

  Private

Contains

End Module comms_base_class_module
