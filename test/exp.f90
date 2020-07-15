program test_exp

  Use iso_fortran_env, Only : wp => real64
  Implicit None

  Real( wp ), Dimension(3) :: r = 0.0_wp
  Real( wp ), Dimension(3) :: g_2dr, g_dr, g_r, f_rdr, f_drdr
  Real( wp ), Dimension(3), Parameter :: dr = [0.02_wp, 0.01_wp, 0.01_wp]
  Real( wp ), Parameter :: alpha = 0.4_wp
  Integer,    Parameter :: nSamp = 100
  Integer :: i

  g_dr = g(dr, alpha)
  g_2dr = g( Sqrt( 2.0_wp ) * dr, alpha)
  g_r = g(r, alpha)
  f_rdr = f(r, dr, alpha)
  f_drdr = f(dr, dr, alpha)

  do i = 0, nSamp
    print('(1x, a,1x,i0.1,3(1x,f8.3))'), "i, r ", i, r
    print*, "apprx", g_r
    print*, "exact", exact(r, alpha)
    print*, "diff ", g_r - exact(r, alpha)
    g_r = g_r * f_rdr * g_dr
    f_rdr = f_rdr * g_2dr
    r = r + dr
  end do

contains

  function g(r, alpha)
    Real( wp ), Dimension(3) :: g
    Real( wp ), Dimension(3), Intent( In    ) :: r
    Real( wp ), Intent( In    ) :: alpha

    g = exp(- (alpha**2) * dot_product(r, r))

  end function g

  function f(r, dr, alpha)
    Real( wp ), Dimension(3) :: f
    Real( wp ), Dimension(3), Intent( In    ) :: r
    Real( wp ), Dimension(3), Intent( In    ) :: dr
    Real( wp ), Intent( In    ) :: alpha

    f = exp(- (alpha**2) * dot_product(r, dr))

  end function f

  function exact(r, alpha)
    Real( wp ), Dimension(3) :: exact
    Real( wp ), Dimension(3), Intent( In    ) :: r
    Real( wp ), Intent( In    ) :: alpha

    exact = exp(- (alpha**2) * dot_product(r,r))
  end function exact


end program test_exp
