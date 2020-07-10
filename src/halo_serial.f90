Module halo_serial_module

  Use constants, Only : wp
  Use halo_setter_base_module, Only : halo_setter_base_class

  Implicit None

  ! Mostly dummy type for the moment, eventually will contain things like MPI communicator
  Type, Public, Extends( halo_setter_base_class ) :: halo_serial_setter
     Private
   Contains
     Procedure, Public :: fill => halo_fill
  End Type halo_serial_setter

  Private

Contains

  Subroutine halo_fill( H, halo_width, hdlb, gin, hout, error )

    Class( halo_serial_setter ),                                 Intent( InOut ) :: H
    Integer,                                                     Intent( In    ) :: halo_width
    Integer,    Dimension( 1:3 ),                                Intent( In    ) :: hdlb
    Real( wp ), Dimension( 0:, 0:, 0: ),                         Intent( In    ) :: gin
    Real( wp ), Dimension( hdlb( 1 ):, hdlb( 2 ):, hdlb( 3 ): ), Intent(   Out ) :: hout
    Integer,                                                     Intent(   Out ) :: error

    Integer, Dimension( 1:3 ) :: gub, gn
    Integer, Dimension( 1:3 ) :: hub, hlb

    Integer, Dimension( 1:3 ) :: hi, gi

    Integer :: hi1, hi2, hi3

    Logical, Parameter :: debug = .True.

    error = 0

    Call H%inc_n_calls()

    ! Bounds on grid without halo - we know the lower bound is zero from declaration
    gub = Ubound( gin )

    ! Bounds on grid with halo
    ! Have to be careful as HOUT might be declared bigger than needed
    hlb =     - halo_width
    hub = gub + halo_width

    ! Size
    gn = gub + 1

    ! Fill the central section
    hout( 0:gub( 1 ), 0:gub( 2 ), 0:gub( 3 ) ) = gin

    ! Need to finish actually filling halos!!!
    ! First all the stuff below the bottom face of the non-haloed grid
    Do hi3 = hlb( 3 ), -1
       Do hi2 = hlb( 2 ), hub( 2 )
          Do hi1 = hlb( 1 ), hub( 1 )
             hi = [ hi1, hi2, hi3 ]
             gi = Modulo( hi, gn )
             hout( hi( 1 ), hi( 2 ), hi( 3 ) ) = gin( gi( 1 ), gi( 2 ), gi( 3 ) )
          End Do
       End Do
    End Do

    ! Now above the top face
    Do hi3 = gub( 3 ) + 1, hub( 3 )
       Do hi2 = hlb( 2 ), hub( 2 )
          Do hi1 = hlb( 1 ), hub( 1 )
             hi = [ hi1, hi2, hi3 ]
             gi = Modulo( hi, gn )
             hout( hi( 1 ), hi( 2 ), hi( 3 ) ) = gin( gi( 1 ), gi( 2 ), gi( 3 ) )
          End Do
       End Do
    End Do

    ! Now the front face
    Do hi3 = 0, gub( 3 )
       Do hi2 = hlb( 2 ), -1
          Do hi1 = hlb( 1 ), hub( 1 )
             hi = [ hi1, hi2, hi3 ]
             gi = Modulo( hi, gn )
             hout( hi( 1 ), hi( 2 ), hi( 3 ) ) = gin( gi( 1 ), gi( 2 ), gi( 3 ) )
          End Do
       End Do
    End Do

    ! Back face
    Do hi3 = 0, gub( 3 )
       Do hi2 = gub( 2 ) + 1, hub( 2 )
          Do hi1 = hlb( 1 ), hub( 1 )
             hi = [ hi1, hi2, hi3 ]
             gi = Modulo( hi, gn )
             hout( hi( 1 ), hi( 2 ), hi( 3 ) ) = gin( gi( 1 ), gi( 2 ), gi( 3 ) )
          End Do
       End Do
    End Do

    ! Left Face
    Do hi3 = 0, gub( 3 )
       Do hi2 = 0, gub( 2 )
          Do hi1 = hlb( 1 ), -1
             hi = [ hi1, hi2, hi3 ]
             gi = Modulo( hi, gn )
             hout( hi( 1 ), hi( 2 ), hi( 3 ) ) = gin( gi( 1 ), gi( 2 ), gi( 3 ) )
          End Do
       End Do
    End Do

    ! Right Face
    Do hi3 = 0, gub( 3 )
       Do hi2 = 0, gub( 2 )
          Do hi1 = gub( 1 ) + 1, hub( 1 )
             hi = [ hi1, hi2, hi3 ]
             gi = Modulo( hi, gn )
             hout( hi( 1 ), hi( 2 ), hi( 3 ) ) = gin( gi( 1 ), gi( 2 ), gi( 3 ) )
          End Do
       End Do
    End Do

    ! Sanity check on fiddly routine
    If( debug ) Then
       Do hi3 = hlb( 3 ), hub( 3 )
          Do hi2 = hlb( 2 ), hub( 2 )
             Do hi1 = hlb( 1 ), hub( 1 )
                hi = [ hi1, hi2, hi3 ]
                gi = Modulo( hi, gn )
                If( Abs( hout( hi( 1 ), hi( 2 ), hi( 3 ) ) - &
                     gin( gi( 1 ), gi( 2 ), gi( 3 ) ) ) > Tiny( hout ) ) Then
                   Write( *, '( a, 1x, 2( 3( i5, 1x ), g20.12, 5x ) )'  ) &
                        'halo filling problem',                           &
                        hi, hout( hi( 1 ), hi( 2 ), hi( 3 ) ),            &
                        gi,  gin( gi( 1 ), gi( 2 ), gi( 3 ) )
                End If
             End Do
          End Do
       End Do
    End If

  End Subroutine halo_fill

End Module halo_serial_module
