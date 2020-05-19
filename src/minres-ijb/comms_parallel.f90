Module comms_parallel_module

  Use comms_base_class_module, Only : comms_base_class

  Implicit None

  Type, Public, Extends( comms_base_class ) :: comms_parallel
   Contains
     Procedure, Public  :: reduce          => reduce_real
     Procedure, Public  :: set_comm        => set_comm
     Procedure, Public  :: get_comm        => get_comm
     Procedure, Public  :: get_rank        => get_rank
     Procedure, Public  :: get_size        => get_size
     Procedure, Public  :: get_proc_coords => get_proc_coords
     Procedure, Public  :: get_proc_grid   => get_proc_grid
  End type comms_parallel

Contains

  Subroutine set_comm( c, comm, error )

    Use mpi_f08, Only : mpi_comm, mpi_topo_test, mpi_cart

    Implicit None

    Class( comms_parallel ), Intent(   Out ) :: c
    Integer                , Intent( In    ) :: comm
    Integer                , Intent(   Out ) :: error

    Type( mpi_comm ) :: comm_f08

    Integer :: comm_type
    
    error = 0

    comm_f08%mpi_val = comm
    Call mpi_topo_test( comm_f08, comm_type )

    If( comm_type == mpi_cart ) Then
       c%communicator = comm
    Else
       error = 1
    End If

  End Subroutine set_comm

  Subroutine get_comm( c, comm )

    Implicit None

    Class( comms_parallel ), Intent( In    ) :: c
    Integer              , Intent(   Out ) :: comm

    comm = c%communicator
    
  End Subroutine get_comm

  Subroutine get_rank( c, data )

    Use mpi_f08, Only : mpi_comm, mpi_comm_rank

    Class( comms_parallel ), Intent( In    ) :: c
    Integer              , Intent(   Out ) :: data

    Type( mpi_comm ) :: comm_f08

    comm_f08%mpi_val = c%communicator
    Call mpi_comm_rank( comm_f08, data )
    
  End Subroutine get_rank

  Subroutine get_size( c, data )

    Use mpi_f08, Only : mpi_comm, mpi_comm_size

    Class( comms_parallel ), Intent( In    ) :: c
    Integer              , Intent(   Out ) :: data

    Type( mpi_comm ) :: comm_f08

    comm_f08%mpi_val = c%communicator
    Call mpi_comm_rank( comm_f08, data )

  End Subroutine get_size

  Subroutine get_proc_coords( c, data )

    Use mpi_f08, Only : mpi_comm, mpi_cart_get

    Class( comms_parallel ), Intent( In    ) :: c
    Integer, Dimension( : ), Intent(   Out ) :: data

    Type( mpi_comm ) :: comm_f08

    Integer, Dimension( 1:3 ) :: dims, coords

    Logical, Dimension( 1:3 ) :: periods

    comm_f08%mpi_val = c%communicator
    Call mpi_cart_get( comm_f08, Size( dims ), dims, periods, coords )

    data = coords

  End Subroutine get_proc_coords

  Subroutine get_proc_grid( c, data )

    Use mpi_f08, Only : mpi_comm, mpi_cart_get

    Class( comms_parallel ), Intent( In    ) :: c
    Integer, Dimension( : ), Intent(   Out ) :: data

    Type( mpi_comm ) :: comm_f08

    Integer, Dimension( 1:3 ) :: dims, coords

    Logical, Dimension( 1:3 ) :: periods

    comm_f08%mpi_val = c%communicator
    Call mpi_cart_get( comm_f08, Size( dims ), dims, periods, coords )

    data = coords

  End Subroutine get_proc_grid

  Subroutine reduce_real( c, data )

    Use mpi_f08, Only : mpi_comm, mpi_sizeof, mpi_type_match_size, mpi_typeclass_real, mpi_datatype, &
         mpi_allreduce, mpi_in_place, mpi_sum

    Use, Intrinsic :: iso_fortran_env, Only :  wp => real64
    
    Class( comms_parallel ), Intent( In    ) :: c
    Real( wp )             , Intent( InOut ) :: data

    Type( mpi_comm     ) :: comm_f08
    Type( mpi_datatype ) :: datatype

    Integer :: data_size

    Call mpi_sizeof( data, data_size )
    Call mpi_type_match_size( mpi_typeclass_real, data_size, datatype )

    comm_f08%mpi_val = c%communicator
    Call mpi_allreduce( mpi_in_place, data, 1, datatype, mpi_sum, comm_f08 )
    
  End Subroutine reduce_real


  
End Module comms_parallel_module
