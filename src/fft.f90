Module fft_module

  Use, Intrinsic :: iso_fortran_env, Only :  wp => real64

  Implicit None

  Public :: fft_fft3d

  Private

  Integer, Parameter :: long = Selected_int_kind( 12 )

  Interface
     Subroutine dfftw_plan_dft_3d( plan, k1, k2, k3, a, b, dir, method )
       Use, Intrinsic :: iso_fortran_env, Only :  wp => real64
       Implicit None
       Integer, Parameter :: long = Selected_int_kind( 12 )
       Integer( long )                             , Intent(   Out ) :: plan
       Integer                                     , Intent( In    ) :: k1
       Integer                                     , Intent( In    ) :: k2
       Integer                                     , Intent( In    ) :: k3
       Complex( wp ), Dimension( 1:k1, 1:k2, 1:k3 ), Intent( In    ) :: a
       Complex( wp ), Dimension( 1:k1, 1:k2, 1:k3 ), Intent( In    ) :: b
       Integer                                     , Intent( In    ) :: dir
       Integer                                     , Intent( In    ) :: method
     End Subroutine dfftw_plan_dft_3d
     Subroutine dfftw_destroy_plan( plan )
       Use, Intrinsic :: iso_fortran_env, Only :  wp => real64
       Implicit None
       Integer, Parameter :: long = Selected_int_kind( 12 )
       Integer( long ), Intent( InOut ) :: plan
     End Subroutine dfftw_destroy_plan
     Subroutine dfftw_execute( plan )
       Use, Intrinsic :: iso_fortran_env, Only :  wp => real64
       Implicit None
       Integer, Parameter :: long = Selected_int_kind( 12 )
       Integer( long ), Intent( InOut ) :: plan
     End Subroutine dfftw_execute
  End Interface

Contains

  Subroutine fft_fft3d( direction, data )

    Integer                             , Intent( In    ) :: direction
    Complex( wp ), Dimension( :, :, :  ), Intent( InOut ) :: data

    Integer( long ) :: plan

    Integer :: k1, k2, k3

    k1 = Size( data, Dim = 1 )
    k2 = Size( data, Dim = 2 )
    k3 = Size( data, Dim = 3 )

    Select Case( direction )

    Case( 1 )
       Call dfftw_plan_dft_3d( plan, k1, k2, k3, data, data,  1, 64 )
       Call dfftw_execute( plan )
       Call dfftw_destroy_plan( plan )


    Case( -1 )
       Call dfftw_plan_dft_3d( plan, k1, k2, k3, data, data, -1, 64 )
       Call dfftw_execute( plan )
       Call dfftw_destroy_plan( plan )

    Case Default
       Stop "FFT does not know if it is coming or going !!"

    End Select

  End Subroutine fft_fft3d

End Module fft_module
