Module comms_serial_module

  Use constants, Only : wp
  Use comms_base_class_module, Only : comms_base_class

  Implicit None

  Type, Public, Extends( comms_base_class ) :: comms_serial
   Contains
     Procedure, Public  :: reduce          => reduce_real
     Procedure, Public  :: max             => reduce_real
     Procedure, Public  :: set_comm        => set_comm
     Procedure, Public  :: get_comm        => get_comm
     Procedure, Public  :: get_rank        => get_rank
     Procedure, Public  :: get_size        => get_size
     Procedure, Public  :: get_proc_coords => get_proc_coords
     Procedure, Public  :: get_proc_grid   => get_proc_grid
  End type comms_serial

Contains

  Subroutine set_comm( c, comm, error )

    Implicit None

    Class( comms_serial ), Intent(   Out ) :: c
    Integer,               Intent( In    ) :: comm
    Integer,               Intent(   Out ) :: error

    c%communicator = comm

    error = 0

  End Subroutine set_comm

  Subroutine get_comm( c, comm )

    Implicit None

    Class( comms_serial ), Intent( In    ) :: c
    Integer,               Intent(   Out ) :: comm

    comm = c%communicator

  End Subroutine get_comm

  Subroutine get_rank( c, data )

    Class( comms_serial ), Intent( In    ) :: c
    Integer,               Intent(   Out ) :: data

    data = 0

  End Subroutine get_rank

  Subroutine get_size( c, data )

    Class( comms_serial ), Intent( In    ) :: c
    Integer,               Intent(   Out ) :: data

    data = 1

  End Subroutine get_size

  Subroutine get_proc_coords( c, data )

    Class( comms_serial ),   Intent( In    ) :: c
    Integer, Dimension( : ), Intent(   Out ) :: data

    data = 0

  End Subroutine get_proc_coords

  Subroutine get_proc_grid( c, data )

    Class( comms_serial ),   Intent( In    ) :: c
    Integer, Dimension( : ), Intent(   Out ) :: data

    data = 1

  End Subroutine get_proc_grid

  Subroutine reduce_real( c, data )

    Class( comms_serial ), Intent( In    ) :: c
    Real( wp ),            Intent( InOut ) :: data

    data = data

  End Subroutine reduce_real

End Module comms_serial_module
