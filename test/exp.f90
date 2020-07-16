program test_exp

  Use iso_fortran_env, Only : wp => real64
  Implicit None

  Real( wp ), Dimension(3) :: r_0 = 0.0_wp, r
  Real( wp ) :: g_r
  Real( wp ), Dimension(3) :: g_r0, g_dr, f_rdr, f_rdr0, f_rdr00, f_drdr
  Real( wp ), Dimension(3,3) :: g_drdr
  Real( wp ), Dimension(3,3), Parameter :: dr = Reshape([&
       0.01_wp, 0.01_wp, 0.03_wp, &
       0.02_wp, 0.04_wp, 0.07_wp, &
       0.03_wp, -0.05_wp, 0.02_wp], [3,3])
  Real( wp ), Parameter :: alpha = 0.4_wp
  Integer,    Parameter :: nSamp = 16
  Integer :: i, j, k

  r = r_0 - sum(nSamp*dr, dim=2)

  g_r0 = g(r, alpha)
  g_r = g_r0(1)

  do i = 1, 3
    g_dr(i) = g(dr(:, i), alpha)
    f_rdr0(i) = f(r, dr(:, i), alpha)
    f_drdr(i) = f(dr(:, i), dr(:, i), alpha)
  end do

  do j = 1, 3
    do i = 1, 3
      g_drdr(i, j) = f(dr(:,i), dr(:,j), alpha)
    end do
  end do

  f_rdr00 = f_rdr0
  f_rdr = f_rdr0

  do k = -nSamp, nSamp
     do j = -nSamp, nSamp
       do i = -nSamp, nSamp
         ! r used only for printing
         r = r_0 + i*dr(:, 1) + j*dr(:, 2) + k*dr(:, 3)
         print('(3(i3.1, 1X),1X,3(f6.3,1X),1X,3(g21.15, 1x) )'), &
              i, j, k, r, g_r, g(r, alpha), g_r - g(r, alpha)
        ! print*, "INFO ",i, j, k, g_r - g(r, alpha)
        ! print*, "F_RDR", f_rdr
        ! print*, "EXACT", [(f(r, dr(:, i), alpha), i=1,3)]
        ! print*, "DIFF ", f_rdr - [(f(r, dr(:, i), alpha), i=1,3)]

        g_r = g_r * f_rdr(1) * g_dr(1)
        f_rdr = f_rdr * g_drdr(:, 1)
      end do

      g_r = g_r0(2) * f_rdr0(2) * g_dr(2)
      g_r0(2) = g_r

      f_rdr = f_rdr0 * g_drdr(:, 2)
      f_rdr0 = f_rdr
    end do

    g_r = g_r0(1) * f_rdr00(3) * g_dr(3)
    g_r0 = g_r

    f_rdr = f_rdr00 * g_drdr(:, 3)
    f_rdr0 = f_rdr
    f_rdr00 = f_rdr
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
