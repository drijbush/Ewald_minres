Module sort_module

  Use, Intrinsic :: iso_fortran_env, Only :  wp => real64

  Implicit None

  Integer, Parameter, Public :: SORT_ASCEND   = 1
  Integer, Parameter, Public :: SORT_DESCEND  = -1

  Public :: sort
  Public :: sort_index
  Public :: sort_move_by_index

  Private

  Interface sort
     Module Procedure sort_real_1d
     Module Procedure sort_integer_1d
     Module Procedure sort_character_1d
  End Interface

  Interface sort_index
     Module procedure index_real_1d
     Module procedure index_integer_1d
     Module procedure index_character_1d
  End Interface

  Interface sort_move_by_index
     Module Procedure move_by_index_real_1d
     Module Procedure move_by_index_integer_1d
     Module Procedure move_by_index_character_1d
  End Interface

  Integer, Parameter :: insert_size = 10

Contains

  Pure Subroutine sort_real_1d( a, dir )

    Real( wp ), Dimension( 1: ), Intent( InOut ) :: a
    Integer                    , Intent( In    ) :: dir

    ! Stack big enough for all arrays indexed by 64 bit integers
    Integer, Dimension( 1:2, 1:64 ) :: qsort_stack

    Real( wp ) :: part_val
    
    Integer :: n
    Integer :: stack_size
    Integer :: part, bot, top
    Integer :: i, j

    n = Size( a )

    stack_size = 0

    bot = 1
    top = n

    Do

       If( top - bot < insert_size ) Then

          Call insert_sort( bot, top, a )

          If( stack_size == 0 ) Then
             Exit
          End If

          Call pop_stack( qsort_stack, stack_size, bot, top )

       Else

          part = ( top + bot ) / 2

          Call swap( part, bot + 1, a )
          If( a( bot ) > a( top ) ) Then
             Call swap( bot, top, a )
          End If
          If( a( bot + 1 ) > a( top ) ) Then
             Call swap( bot + 1, top, a )
          End If
          If( a( bot ) > a( bot + 1 ) ) Then
             Call swap( bot, bot + 1, a )
          End If

          i = bot + 1
          j = top

          part_val = a( bot + 1 )

          Do

             Do 
                i = i + 1
                If( a( i ) >= part_val ) Then
                   Exit
                End If
             End Do

             Do 
                j = j - 1
                If( a( j ) <= part_val ) Then
                   Exit
                End If
             End Do

             If( j < i ) Then
                Exit
             End If

             Call swap( i, j, a )

          End Do

          a( bot + 1 ) = a( j )
          a( j       ) = part_val

          If( top - i + 1 >= j - bot ) Then
             Call push_stack(   i,   top, stack_size, qsort_stack )
             top = j - 1
          Else
             Call push_stack( bot, j - 1, stack_size, qsort_stack )
             bot = i
          End If

       End If

    End Do

    If( dir < 0 ) Then
       Do i = 1, ( n + 1 ) / 2
          Call swap( i, n - i + 1, a )
       End Do
    End If

  Contains

    Pure Subroutine insert_sort( bot, top, a )

      Integer                    , Intent( In    ) :: bot
      Integer                    , Intent( In    ) :: top
      Real( wp ), Dimension( 1: ), Intent( InOut ) :: a

      Real( wp ) :: tmp

      Integer :: i, j

      Do j = bot + 1, top
         tmp = a( j )
         Do i = j - 1, bot, -1
            If( a( i ) < tmp ) Then
               Exit
            End If
            a( i + 1 ) = a( i )
         End Do
         a( i + 1 ) = tmp
      End Do

    End Subroutine insert_sort

    Pure Subroutine swap( i, j, a )
      
      Integer                    , Intent( In    ) :: i
      Integer                    , Intent( In    ) :: j
      Real( wp ), Dimension( 1: ), Intent( InOut ) :: a

      Real( wp ) :: tmp

      tmp    = a( i )
      a( i ) = a( j )
      a( j ) = tmp

    End Subroutine swap

  End Subroutine sort_real_1d

  Pure Subroutine sort_integer_1d( a, dir )

    Integer, Dimension( 1: ), Intent( InOut ) :: a
    Integer                 , Intent( In    ) :: dir

    ! Stack big enough for all arrays indexed by 64 bit integers
    Integer, Dimension( 1:2, 1:64 ) :: qsort_stack

    Integer :: part_val
    Integer :: n
    Integer :: stack_size
    Integer :: part, bot, top
    Integer :: i, j

    n = Size( a )

    stack_size = 0

    bot = 1
    top = n

    Do

       If( top - bot < insert_size ) Then

          Call insert_sort( bot, top, a )

          If( stack_size == 0 ) Then
             Exit
          End If

          Call pop_stack( qsort_stack, stack_size, bot, top )

       Else

          part = ( top + bot ) / 2

          Call swap( part, bot + 1, a )
          If( a( bot ) > a( top ) ) Then
             Call swap( bot, top, a )
          End If
          If( a( bot + 1 ) > a( top ) ) Then
             Call swap( bot + 1, top, a )
          End If
          If( a( bot ) > a( bot + 1 ) ) Then
             Call swap( bot, bot + 1, a )
          End If

          i = bot + 1
          j = top

          part_val = a( bot + 1 )

          Do

             Do 
                i = i + 1
                If( a( i ) >= part_val ) Then
                   Exit
                End If
             End Do

             Do 
                j = j - 1
                If( a( j ) <= part_val ) Then
                   Exit
                End If
             End Do

             If( j < i ) Then
                Exit
             End If

             Call swap( i, j, a )

          End Do

          a( bot + 1 ) = a( j )
          a( j       ) = part_val

          If( top - i + 1 >= j - bot ) Then
             Call push_stack(   i,   top, stack_size, qsort_stack )
             top = j - 1
          Else
             Call push_stack( bot, j - 1, stack_size, qsort_stack )
             bot = i
          End If

       End If

    End Do

    If( dir < 0 ) Then
       Do i = 1, ( n + 1 ) / 2
          Call swap( i, n - i + 1, a )
       End Do
    End If

  Contains

    Pure Subroutine insert_sort( bot, top, a )

      Integer                 , Intent( In    ) :: bot
      Integer                 , Intent( In    ) :: top
      Integer, Dimension( 1: ), Intent( InOut ) :: a

      Integer :: tmp

      Integer :: i, j

      Do j = bot + 1, top
         tmp = a( j )
         Do i = j - 1, bot, -1
            If( a( i ) < tmp ) Then
               Exit
            End If
            a( i + 1 ) = a( i )
         End Do
         a( i + 1 ) = tmp
      End Do

    End Subroutine insert_sort

    Pure Subroutine swap( i, j, a )
      
      Integer                 , Intent( In    ) :: i
      Integer                 , Intent( In    ) :: j
      Integer, Dimension( 1: ), Intent( InOut ) :: a

      Integer :: tmp

      tmp    = a( i )
      a( i ) = a( j )
      a( j ) = tmp

    End Subroutine swap

  End Subroutine sort_integer_1d

  Pure Subroutine sort_character_1d( a, dir )

    Character( Len = * ), Dimension( 1: ), Intent( InOut ) :: a
    Integer                              , Intent( In    ) :: dir

    ! Stack big enough for all arrays indexed by 64 bit integers
    Integer, Dimension( 1:2, 1:64 ) :: qsort_stack

    Character( Len = Len( a ) ) :: part_val
    
    Integer :: n
    Integer :: stack_size
    Integer :: part, bot, top
    Integer :: i, j

    n = Size( a )

    stack_size = 0

    bot = 1
    top = n

    Do

       If( top - bot < insert_size ) Then

          Call insert_sort( bot, top, a )

          If( stack_size == 0 ) Then
             Exit
          End If

          Call pop_stack( qsort_stack, stack_size, bot, top )

       Else

          part = ( top + bot ) / 2

          Call swap( part, bot + 1, a )
          If( a( bot ) > a( top ) ) Then
             Call swap( bot, top, a )
          End If
          If( a( bot + 1 ) > a( top ) ) Then
             Call swap( bot + 1, top, a )
          End If
          If( a( bot ) > a( bot + 1 ) ) Then
             Call swap( bot, bot + 1, a )
          End If

          i = bot + 1
          j = top

          part_val = a( bot + 1 )

          Do

             Do 
                i = i + 1
                If( a( i ) >= part_val ) Then
                   Exit
                End If
             End Do

             Do 
                j = j - 1
                If( a( j ) <= part_val ) Then
                   Exit
                End If
             End Do

             If( j < i ) Then
                Exit
             End If

             Call swap( i, j, a )

          End Do

          a( bot + 1 ) = a( j )
          a( j       ) = part_val

          If( top - i + 1 >= j - bot ) Then
             Call push_stack(   i,   top, stack_size, qsort_stack )
             top = j - 1
          Else
             Call push_stack( bot, j - 1, stack_size, qsort_stack )
             bot = i
          End If

       End If

    End Do

    If( dir < 0 ) Then
       Do i = 1, ( n + 1 ) / 2
          Call swap( i, n - i + 1, a )
       End Do
    End If

  Contains

    Pure Subroutine insert_sort( bot, top, a )

      Integer                              , Intent( In    ) :: bot
      Integer                              , Intent( In    ) :: top
      Character( Len = * ), Dimension( 1: ), Intent( InOut ) :: a

      Character( Len = Len( a ) ) :: tmp

      Integer :: i, j

      Do j = bot + 1, top
         tmp = a( j )
         Do i = j - 1, bot, -1
            If( a( i ) < tmp ) Then
               Exit
            End If
            a( i + 1 ) = a( i )
         End Do
         a( i + 1 ) = tmp
      End Do

    End Subroutine insert_sort

    Pure Subroutine swap( i, j, a )
      
      Integer                              , Intent( In    ) :: i
      Integer                              , Intent( In    ) :: j
      Character( Len = * ), Dimension( 1: ), Intent( InOut ) :: a

      Character( Len = Len( a ) ) :: tmp

      tmp    = a( i )
      a( i ) = a( j )
      a( j ) = tmp

    End Subroutine swap

  End Subroutine sort_character_1d

  Pure Subroutine index_real_1d( a, dir, ind )

    Real( wp ), Dimension( 1: ), Intent( In    ) :: a
    Integer                    , Intent( In    ) :: dir
    Integer   , Dimension( 1: ), Intent(   Out ) :: ind

    ! Stack big enough for all arrays indexed by 64 bit integers
    Integer, Dimension( 1:2, 1:64 ) :: qsort_stack

    Real( wp ) :: part_val
    
    Integer :: n
    Integer :: ind_tmp
    Integer :: stack_size
    Integer :: part, bot, top
    Integer :: i, j

    n = Size( a )

    ind = (/ ( i, i = 1, n ) /)

    stack_size = 0

    bot = 1
    top = n

    Do

       If( top - bot < insert_size ) Then

          Call insert_index( bot, top, a, ind )

          If( stack_size == 0 ) Then
             Exit
          End If

          Call pop_stack( qsort_stack, stack_size, bot, top )

       Else

          part = ( top + bot ) / 2

          Call swap( part, bot + 1, ind )
          If( a( ind( bot ) ) > a( ind( top ) ) ) Then
             Call swap( bot, top, ind )
          End If
          If( a( ind( bot + 1 ) ) > a( ind( top ) ) ) Then
             Call swap( bot + 1, top, ind )
          End If
          If( a( ind( bot ) ) > a( ind( bot + 1 ) ) ) Then
             Call swap( bot, bot + 1, ind )
          End If

          i = bot + 1
          j = top

          ind_tmp  = ind( bot + 1 )
          part_val = a( ind( bot + 1 ) )

          Do

             Do 
                i = i + 1
                If( a( ind( i ) ) >= part_val ) Then
                   Exit
                End If
             End Do

             Do 
                j = j - 1
                If( a( ind( j ) ) <= part_val ) Then
                   Exit
                End If
             End Do

             If( j < i ) Then
                Exit
             End If

             Call swap( i, j, ind )

          End Do

          ind( bot + 1 ) = ind( j )
          ind( j       ) = ind_tmp

          If( top - i + 1 >= j - bot ) Then
             Call push_stack(   i,   top, stack_size, qsort_stack )
             top = j - 1
          Else
             Call push_stack( bot, j - 1, stack_size, qsort_stack )
             bot = i
          End If

       End If

    End Do

    If( dir < 0 ) Then
       Do i = 1, ( n + 1 ) / 2
          Call swap( i, n - i + 1, ind )
       End Do
    End If

  Contains

    Pure Subroutine insert_index( bot, top, a, ind )

      Integer                    , Intent( In    ) :: bot
      Integer                    , Intent( In    ) :: top
      Real( wp ), Dimension( 1: ), Intent( In    ) :: a
      Integer   , Dimension( 1: ), Intent( InOut ) :: ind

      Real( wp ) :: tmp

      Integer :: ind_tmp
      Integer :: i, j

      Do j = bot + 1, top
         ind_tmp = ind( j )
         tmp     = a( ind_tmp )
         Do i = j - 1, bot, -1
            If( a( ind( i ) ) < tmp ) Then
               Exit
            End If
            ind( i + 1 ) = ind( i )
         End Do
         ind( i + 1 ) = ind_tmp
      End Do

    End Subroutine insert_index

    Pure Subroutine swap( i, j, a )
      
      Integer                 , Intent( In    ) :: i
      Integer                 , Intent( In    ) :: j
      Integer, Dimension( 1: ), Intent( InOut ) :: a

      Integer :: tmp

      tmp    = a( i )
      a( i ) = a( j )
      a( j ) = tmp

    End Subroutine swap

  End Subroutine index_real_1d

  Pure Subroutine index_integer_1d( a, dir, ind )

    Integer, Dimension( 1: ), Intent( In    ) :: a
    Integer                 , Intent( In    ) :: dir
    Integer, Dimension( 1: ), Intent(   Out ) :: ind

    ! Stack big enough for all arrays indexed by 64 bit integers
    Integer, Dimension( 1:2, 1:64 ) :: qsort_stack

    Integer :: part_val
    
    Integer :: n
    Integer :: ind_tmp
    Integer :: stack_size
    Integer :: part, bot, top
    Integer :: i, j

    n = Size( a )

    ind = (/ ( i, i = 1, n ) /)

    stack_size = 0

    bot = 1
    top = n

    Do

       If( top - bot < insert_size ) Then

          Call insert_index( bot, top, a, ind )

          If( stack_size == 0 ) Then
             Exit
          End If

          Call pop_stack( qsort_stack, stack_size, bot, top )

       Else

          part = ( top + bot ) / 2

          Call swap( part, bot + 1, ind )
          If( a( ind( bot ) ) > a( ind( top ) ) ) Then
             Call swap( bot, top, ind )
          End If
          If( a( ind( bot + 1 ) ) > a( ind( top ) ) ) Then
             Call swap( bot + 1, top, ind )
          End If
          If( a( ind( bot ) ) > a( ind( bot + 1 ) ) ) Then
             Call swap( bot, bot + 1, ind )
          End If

          i = bot + 1
          j = top

          ind_tmp  = ind( bot + 1 )
          part_val = a( ind( bot + 1 ) )

          Do

             Do 
                i = i + 1
                If( a( ind( i ) ) >= part_val ) Then
                   Exit
                End If
             End Do

             Do 
                j = j - 1
                If( a( ind( j ) ) <= part_val ) Then
                   Exit
                End If
             End Do

             If( j < i ) Then
                Exit
             End If

             Call swap( i, j, ind )

          End Do

          ind( bot + 1 ) = ind( j )
          ind( j       ) = ind_tmp

          If( top - i + 1 >= j - bot ) Then
             Call push_stack(   i,   top, stack_size, qsort_stack )
             top = j - 1
          Else
             Call push_stack( bot, j - 1, stack_size, qsort_stack )
             bot = i
          End If

       End If

    End Do

    If( dir < 0 ) Then
       Do i = 1, ( n + 1 ) / 2
          Call swap( i, n - i + 1, ind )
       End Do
    End If

  Contains

    Pure Subroutine insert_index( bot, top, a, ind )

      Integer                 , Intent( In    ) :: bot
      Integer                 , Intent( In    ) :: top
      Integer, Dimension( 1: ), Intent( In    ) :: a
      Integer, Dimension( 1: ), Intent( InOut ) :: ind

      Integer :: tmp
      Integer :: ind_tmp
      Integer :: i, j

      Do j = bot + 1, top
         ind_tmp = ind( j )
         tmp     = a( ind_tmp )
         Do i = j - 1, bot, -1
            If( a( ind( i ) ) < tmp ) Then
               Exit
            End If
            ind( i + 1 ) = ind( i )
         End Do
         ind( i + 1 ) = ind_tmp
      End Do

    End Subroutine insert_index

    Pure Subroutine swap( i, j, a )
      
      Integer                 , Intent( In    ) :: i
      Integer                 , Intent( In    ) :: j
      Integer, Dimension( 1: ), Intent( InOut ) :: a

      Integer :: tmp

      tmp    = a( i )
      a( i ) = a( j )
      a( j ) = tmp

    End Subroutine swap

  End Subroutine index_integer_1d

  Pure Subroutine index_character_1d( a, dir, ind )

    Character( Len = * ), Dimension( 1: ), Intent( In    ) :: a
    Integer                              , Intent( In    ) :: dir
    Integer             , Dimension( 1: ), Intent(   Out ) :: ind

    ! Stack big enough for all arrays indexed by 64 bit integers
    Integer, Dimension( 1:2, 1:64 ) :: qsort_stack

    Character( Len = Len( a ) ) :: part_val
    
    Integer :: n
    Integer :: ind_tmp
    Integer :: stack_size
    Integer :: part, bot, top
    Integer :: i, j

    n = Size( a )

    ind = (/ ( i, i = 1, n ) /)

    stack_size = 0

    bot = 1
    top = n

    Do

       If( top - bot < insert_size ) Then

          Call insert_index( bot, top, a, ind )

          If( stack_size == 0 ) Then
             Exit
          End If

          Call pop_stack( qsort_stack, stack_size, bot, top )

       Else

          part = ( top + bot ) / 2

          Call swap( part, bot + 1, ind )
          If( a( ind( bot ) ) > a( ind( top ) ) ) Then
             Call swap( bot, top, ind )
          End If
          If( a( ind( bot + 1 ) ) > a( ind( top ) ) ) Then
             Call swap( bot + 1, top, ind )
          End If
          If( a( ind( bot ) ) > a( ind( bot + 1 ) ) ) Then
             Call swap( bot, bot + 1, ind )
          End If

          i = bot + 1
          j = top

          ind_tmp  = ind( bot + 1 )
          part_val = a( ind( bot + 1 ) )

          Do

             Do 
                i = i + 1
                If( a( ind( i ) ) >= part_val ) Then
                   Exit
                End If
             End Do

             Do 
                j = j - 1
                If( a( ind( j ) ) <= part_val ) Then
                   Exit
                End If
             End Do

             If( j < i ) Then
                Exit
             End If

             Call swap( i, j, ind )

          End Do

          ind( bot + 1 ) = ind( j )
          ind( j       ) = ind_tmp

          If( top - i + 1 >= j - bot ) Then
             Call push_stack(   i,   top, stack_size, qsort_stack )
             top = j - 1
          Else
             Call push_stack( bot, j - 1, stack_size, qsort_stack )
             bot = i
          End If

       End If

    End Do

    If( dir < 0 ) Then
       Do i = 1, ( n + 1 ) / 2
          Call swap( i, n - i + 1, ind )
       End Do
    End If

  Contains

    Pure Subroutine insert_index( bot, top, a, ind )

      Integer                              , Intent( In    ) :: bot
      Integer                              , Intent( In    ) :: top
      Character( Len = * ), Dimension( 1: ), Intent( In    ) :: a
      Integer             , Dimension( 1: ), Intent( InOut ) :: ind

      Character( Len = Len( a ) ) :: tmp

      Integer :: ind_tmp
      Integer :: i, j

      Do j = bot + 1, top
         ind_tmp = ind( j )
         tmp     = a( ind_tmp )
         Do i = j - 1, bot, -1
            If( a( ind( i ) ) < tmp ) Then
               Exit
            End If
            ind( i + 1 ) = ind( i )
         End Do
         ind( i + 1 ) = ind_tmp
      End Do

    End Subroutine insert_index

    Pure Subroutine swap( i, j, a )
      
      Integer                 , Intent( In    ) :: i
      Integer                 , Intent( In    ) :: j
      Integer, Dimension( 1: ), Intent( InOut ) :: a

      Integer :: tmp

      tmp    = a( i )
      a( i ) = a( j )
      a( j ) = tmp

    End Subroutine swap

  End Subroutine index_character_1d

  Pure Subroutine move_by_index_real_1d( ind, a )

    Integer   , Dimension( : ), Intent( In    ) :: ind
    Real( wp ), Dimension( : ), Intent( InOut ) :: a

    a = a( (/ ind /) )

  End Subroutine move_by_index_real_1d

  Pure Subroutine move_by_index_integer_1d( ind, a )

    Integer, Dimension( : ), Intent( In    ) :: ind
    Integer, Dimension( : ), Intent( InOut ) :: a

    a = a( (/ ind /) )

  End Subroutine move_by_index_integer_1d

  Pure Subroutine move_by_index_character_1d( ind, a )

    Integer             , Dimension( : ), Intent( In    ) :: ind
    Character( Len = * ), Dimension( : ), Intent( InOut ) :: a

    a = a( (/ ind /) )

  End Subroutine move_by_index_character_1d

  Pure Subroutine pop_stack( stack, stack_size, bot, top )

    Integer, Dimension( :, : ), Intent( In    ) :: stack
    Integer                   , Intent( InOut ) :: stack_size
    Integer                   , Intent(   Out ) :: bot
    Integer                   , Intent(   Out ) :: top

    bot = stack( 1, stack_size )
    top = stack( 2, stack_size )
      
    stack_size = stack_size - 1

  End Subroutine pop_stack

  Pure Subroutine push_stack( bot, top, stack_size, stack )
    
    Integer                   , Intent( In    ) :: bot
    Integer                   , Intent( In    ) :: top
    Integer                   , Intent( InOut ) :: stack_size
    Integer, Dimension( :, : ), Intent( InOut ) :: stack
    
    stack_size = stack_size + 1
    
    stack( 1, stack_size ) = bot
    stack( 2, stack_size ) = top
    
  End Subroutine push_stack

End Module sort_module

