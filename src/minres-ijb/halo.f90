Module halo_module
  
  Use, Intrinsic :: iso_fortran_env, Only :  wp => real64

  Implicit None

  ! Mostly dummy type for the moment, eventually will contain things like MPI communicator
  Type, Public :: halo_setter
     Private
   Contains
     Procedure, Public :: allocate => halo_allocate
     Procedure, Public :: fill     => halo_fill
     Procedure, Public :: free     => halo_free
  End Type halo_setter

  Private
  
Contains

  Pure Subroutine halo_allocate( H, glb, gub, halo_width, hout )

    Class( halo_setter )            ,              Intent( In    ) :: H
    Integer   , Dimension( 1:3     ),              Intent( In    ) :: glb
    Integer   , Dimension( 1:3     ),              Intent( In    ) :: gub
    Integer   ,                                    Intent( In    ) :: halo_width
    Real( wp ), Dimension( :, :, : ), Allocatable, Intent(   Out ) :: hout

    Integer, Dimension( 1:3 ) :: hlb, hub

    hlb       = glb - halo_width
    hub       = gub + halo_width
    Allocate( hout( hlb( 1 ):hub( 1 ), hlb( 2 ):hub( 2 ), hlb( 3 ):hub( 3 ) ) )
  
  End Subroutine halo_allocate

  Subroutine halo_fill( H, halo_width, hdlb, gin, hout )

    Class( halo_setter ),                                        Intent( In    ) :: H
    Integer   ,                                                  Intent( In    ) :: halo_width
    Integer   , Dimension( 1:3 ),                                Intent( In    ) :: hdlb
    Real( wp ), Dimension( 0:, 0:, 0: )                        , Intent( In    ) :: gin
    Real( wp ), Dimension( hdlb( 1 ):, hdlb( 2 ):, hdlb( 3 ): ), Intent(   Out ) :: hout

    Integer, Dimension( 1:3 ) :: gub, gn
    Integer, Dimension( 1:3 ) :: hub, hlb

    Integer, Dimension( 1:3 ) :: hi, gi

    Integer :: hi1, hi2, hi3

    Logical, Parameter :: debug = .True.

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

  Pure Subroutine halo_free( H, hout )

    Class( halo_setter )            ,              Intent( In    ) :: H
    Real( wp ), Dimension( :, :, : ), Allocatable, Intent(   Out ) :: hout

    Deallocate( hout )

  End Subroutine halo_free

End Module halo_module
