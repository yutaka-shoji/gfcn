module cmils
  use ISO_C_binding
  use omp_lib
  implicit none
  private
  public calc_responsefcn, set_params
  ! integrand parameter
  double precision :: RA, RB, RD
  double precision :: alpha, k
  double precision :: Fo
  integer :: n_i
  !$omp threadprivate(RA, RB, RD, alpha, k, Fo, n_i)
contains
  !----------------------------------------------------------------------
  subroutine set_params(RA_in, RB_in, RD_in, alpha_in, k_in) &
      bind(C, name='set_params')
    !----------------------------------------------------------------------
    double precision, intent(in) :: RA_in, RB_in, RD_in
    double precision, intent(in) :: alpha_in, k_in
    !----------------------------------------------------------------------
    RA = RA_in
    RB = RB_in
    RD = RD_in
    alpha = alpha_in
    k = k_in
  end subroutine set_params
  !----------------------------------------------------------------------
  subroutine calc_responsefcn(nt, Fo_arr, fcn) &
      bind(C, name='calc_responsefcn')
    !----------------------------------------------------------------------
    integer :: i, j, jj
    !----------------------------------------------------------------------
    integer, intent(in) :: nt
    double precision, dimension(nt), intent(in) :: Fo_arr
    double precision, dimension(nt), intent(out) :: fcn
    !----------------------------------------------------------------------
    integer, parameter :: n_int = 6, nn_int = n_int/2+1
    integer, parameter, dimension(nn_int) :: n_arr = (/ (i*2,i=0,nn_int-1) /)
    double precision, dimension(nt) :: I1, I2
    double precision, dimension(nt, nn_int) :: tmpI2
    double precision :: abserr
    double precision, parameter :: ep = 1.0d-3, V = 6.0d2
    ! for dqag-------------------------------------------------------------
    double precision :: s
    double precision, parameter :: epsabs = 1.0d-6, epsrel = 0.0d0
    integer, parameter :: key = 6
    integer, parameter :: limit = 400, lenw = limit*4
    integer :: neval, ier, last
    integer, dimension(limit) :: iwork
    double precision, dimension(limit*4) :: work
    !----------------------------------------------------------------------
    I1 = ep*ep * Fo_arr / (2.0d0 * k)

    !$omp parallel do private(j, jj, s, abserr, neval, ier, last, iwork, work) copyin(RA, RB, RD, alpha, k)
    do j = 1, nt
      Fo = Fo_arr(j)
      do jj = 1, nn_int
        n_i = n_arr(jj)
        call dqag(integrand,ep,V,epsabs,epsrel,key,s,abserr, &
          neval,ier,limit,lenw,last,iwork,work)
        tmpI2(j,jj) = s
      end do
    end do
    !$omp end parallel do

    tmpI2(:,2:) = 2 * tmpI2(:,2:)
    I2 = sum(tmpI2, dim=2)
    fcn = I1 + I2
  end subroutine calc_responsefcn
  !----------------------------------------------------------------------
  subroutine calc_responsefcn_dblU(nt, Fo_arr, fcn) &
      bind(C, name='calc_responsefcn_dblU')
    !----------------------------------------------------------------------
    integer :: i, j, jj
    !----------------------------------------------------------------------
    integer, intent(in) :: nt
    double precision, dimension(nt), intent(in) :: Fo_arr
    double precision, dimension(nt), intent(out) :: fcn
    !----------------------------------------------------------------------
    integer, parameter :: n_int = 6, nn_int = n_int/2+1
    integer, parameter, dimension(nn_int) :: n_arr = (/ (i*4,i=0,nn_int-1) /)
    double precision, dimension(nt) :: I1, I2
    double precision, dimension(nt, nn_int) :: tmpI2
    double precision :: abserr
    double precision, parameter :: ep = 1.0d-3, V = 6.0d2
    ! for dqag-------------------------------------------------------------
    double precision :: s
    double precision, parameter :: epsabs = 1.0d-6, epsrel = 0.0d0
    integer, parameter :: key = 6
    integer, parameter :: limit = 400, lenw = limit*4
    integer :: neval, ier, last
    integer, dimension(limit) :: iwork
    double precision, dimension(limit*4) :: work
    !----------------------------------------------------------------------
    I1 = ep*ep * Fo_arr / (2.0d0 * k)

    !$omp parallel do private(j, jj, s, abserr, neval, ier, last, iwork, work) copyin(RA, RB, RD, alpha, k)
    do j = 1, nt
      Fo = Fo_arr(j)
      do jj = 1, nn_int
        n_i = n_arr(jj)
        call dqag(integrand,ep,V,epsabs,epsrel,key,s,abserr, &
          neval,ier,limit,lenw,last,iwork,work)
        tmpI2(j,jj) = s
      end do
    end do
    !$omp end parallel do

    tmpI2(:,2:) = 2 * tmpI2(:,2:)
    I2 = sum(tmpI2, dim=2)
    fcn = I1 + I2
  end subroutine calc_responsefcn_dblU
  !----------------------------------------------------------------------
  function integrand(x)
    double precision :: integrand
    double precision, intent(in) :: x

    double precision :: phi, psi, f, g
    double precision :: jn_x, jnp_x, jn_ax, jnp_ax
    double precision :: yn_x, ynp_x, yn_ax, ynp_ax
    double precision :: ax
    !-------------------------------------------------------------------------
    ax = alpha * x

    jn_x  = bessel_jn(n_i, x)
    jnp_x = n_i/x * jn_x - bessel_jn(n_i+1, x)
    jn_ax  = bessel_jn(n_i, ax)
    jnp_ax = n_i/ax * jn_ax - bessel_jn(n_i+1, ax)
    yn_x  = bessel_yn(n_i, x)
    ynp_x = n_i/x * yn_x - bessel_yn(n_i+1, x)
    yn_ax  = bessel_yn(n_i, ax)
    ynp_ax = n_i/ax * yn_ax - bessel_yn(n_i+1, ax)

    phi = alpha * k * jn_x * jnp_ax - jnp_x * jn_ax
    psi = alpha * k * jn_x * ynp_ax - jnp_x * yn_ax
    f   = alpha * k * yn_x * jnp_ax - ynp_x * jn_ax
    g   = alpha * k * yn_x * ynp_ax - ynp_x * yn_ax

    integrand = ( 1.0d0 - exp(-x*x * Fo) ) &
      * ( bessel_jn(n_i, RA*x) + bessel_jn(n_i, RB*x) ) * 0.5d0 &
      * bessel_jn(n_i, RD*x) &
      * ( phi * g - psi * f ) / ( x * ( phi*phi + psi*psi ) )
  end function integrand
  !----------------------------------------------------------------------
end module cmils
