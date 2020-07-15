program test_exp

  Use iso_fortran_env, Only : wp => real64
  Implicit None

  Real( wp ), Dimension(3) :: r_0 = 0.0_wp, r
  Real( wp ) :: g_r
  Real( wp ), Dimension(3) :: g_r0, g_2dr, g_dr, f_rdr, f_rdr0, f_drdr
  Real( wp ), Dimension(3,3), Parameter :: dr = Reshape([&
       0.01_wp, 0.00_wp, 0.00_wp, &
       0.00_wp, 0.01_wp, 0.00_wp, &
       0.00_wp, 0.00_wp, 0.01_wp], [3,3])
  Real( wp ), Parameter :: alpha = 0.4_wp
  Integer,    Parameter :: nSamp = 10
  Integer :: i, j, k

  g_r0 = g(r, alpha)
  g_r = g_r0(1)
  do i = 1, 3
    print*, dr(:, i)
    g_dr(i) = g(dr(:, i), alpha)
    g_2dr(i) = g( Sqrt( 2.0_wp ) * dr(:, i), alpha)
    f_rdr0(i) = f(r, dr(:, i), alpha)
    f_drdr(i) = f(dr(:, i), dr(:, i), alpha)
  end do

  r = r_0
  f_rdr = f_rdr0

  i=0;j=0;k=0

  print*, 'i   j   k     rx     ry     rz           apprx                exact                 diff'
  do k = 0, nSamp
    do j = 0, nSamp
      do i = 0, nSamp
        r = r_0 + i*dr(:, 1) + j*dr(:, 2) + k*dr(:, 3)
        print('(3(i3.1, 1X),1X,3(f6.3,1X),1X,3(g21.15, 1x) )'), i, j, k, r, g_r, g(r, alpha), g_r - g(r, alpha)
        g_r = g_r * f_rdr(1) * g_dr(1)
        f_rdr(1) = f_rdr(1) * g_2dr(1)
      end do
      f_rdr(1) = f_rdr0(1)
      g_r = g_r0(2)
      g_r = g_r * f_rdr(2) * g_dr(2)
      g_r0(2) = g_r
      f_rdr(2) = f_rdr(2) * g_2dr(2)
    end do
    f_rdr(2) = f_rdr0(2)
    g_r = g_r0(1)
    g_r = g_r * f_rdr(3) * g_dr(3)
    g_r0(1) = g_r
    g_r0(2) = g_r
    f_rdr(3) = f_rdr(3) * g_2dr(3)
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
