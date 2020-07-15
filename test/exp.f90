program test_exp

  Use iso_fortran_env, Only : wp => real64
  Implicit None

  Real( wp ), Dimension(3) :: r = 0.0_wp
  Real( wp ) :: g_r
  Real( wp ), Dimension(3) :: g_2dr, g_dr, f_rdr, f_drdr
  Real( wp ), Dimension(3,3), Parameter :: dr = Reshape([&
       0.01_wp, 0.00_wp, 0.00_wp, &
       0.00_wp, 0.00_wp, 0.00_wp, &
       0.00_wp, 0.00_wp, 0.00_wp], [3,3])
  Real( wp ), Parameter :: alpha = 0.4_wp
  Integer,    Parameter :: nSamp = 10
  Integer :: i, j, k

  g_r = g(r, alpha)
  do i = 1, 3
    g_dr(i) = g(dr(:, i), alpha)
    g_2dr(i) = g( Sqrt( 2.0_wp ) * dr(:, i), alpha)
    f_rdr(i) = f(r, dr(:, i), alpha)
    f_drdr(i) = f(dr(:, i), dr(:, i), alpha)
  end do

  print*, 'i   j   k     rx     ry     rz           apprx                exact                 diff'
  do k = 0, nSamp
    do j = 0, nSamp
      do i = 0, nSamp
        print('(3(i3.1, 1X),1X,3(f6.3,1X),1X,3(g21.15, 1x) )'), i, j, k, r, g_r, g(r, alpha), g_r - g(r, alpha)
        g_r = g_r * f_rdr(1) * g_dr(1)
        f_rdr(1) = f_rdr(1) * g_2dr(1)
        r = r + dr(:, 1)
      end do
      g_r = g_r * f_rdr(2) * g_dr(2)
      f_rdr(2) = f_rdr(2) * g_2dr(2)
      r = r + dr(:, 2)
    end do
    g_r = g_r * f_rdr(3) * g_dr(3)
    f_rdr(3) = f_rdr(3) * g_2dr(3)
    r = r + dr(:, 3)
  end do

contains

  function g(r, alpha)
    Real( wp ) :: g
    Real( wp ), Dimension(3), Intent( In    ) :: r
    Real( wp ), Intent( In    ) :: alpha

    g = exp(- (alpha**2) * dot_product(r, r))

  end function g

  function f(r, dr, alpha)
    Real( wp ) :: f
    Real( wp ), Dimension(3), Intent( In    ) :: r
    Real( wp ), Dimension(3), Intent( In    ) :: dr
    Real( wp ), Intent( In    ) :: alpha

    f = exp(-2.0_wp * alpha**2 * dot_product(r, dr))

  end function f

end program test_exp
