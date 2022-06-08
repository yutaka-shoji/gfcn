module ics
  use ISO_C_binding
  use omp_lib
  implicit none
  private
  public calc_responsefcn, set_params
  ! integrand parameter
  double precision :: R
  double precision :: Fo
  !$omp threadprivate(R, Fo)
contains
  !----------------------------------------------------------------------
  subroutine set_params(R_in) &
      bind(C, name='set_params')
    !----------------------------------------------------------------------
    double precision, intent(in) :: R_in
    !----------------------------------------------------------------------
    R = R_in
  end subroutine set_params
  !----------------------------------------------------------------------
  subroutine calc_responsefcn(nt, Fo_arr, fcn) &
      bind(C, name='calc_responsefcn')
    !----------------------------------------------------------------------
    integer :: j
    !----------------------------------------------------------------------
    integer, intent(in) :: nt
    double precision, dimension(nt), intent(in) :: Fo_arr
    double precision, dimension(nt), intent(out) :: fcn
    !----------------------------------------------------------------------
    double precision :: abserr
    double precision :: ep = 0.0d0, V = 1.0d3
    ! for dqag-------------------------------------------------------------
    double precision :: s
    double precision, parameter :: epsabs = 1.0d-6, epsrel = 0.0d0
    integer, parameter :: key = 6
    integer, parameter :: limit = 400, lenw = limit*4
    integer :: neval, ier, last
    integer, dimension(limit) :: iwork
    double precision, dimension(limit*4) :: work
    !----------------------------------------------------------------------
    !$omp parallel do private(j, s, abserr, neval, ier, last, iwork, work) copyin(R)
    do j = 1, nt
      Fo = Fo_arr(j)
      call dqag(integrand,0,V,epsabs,epsrel,key,s,abserr, &
        neval,ier,limit,lenw,last,iwork,work)
      fcn(j) = s
    end do
    !$omp end parallel do
  end subroutine calc_responsefcn
  !----------------------------------------------------------------------
  function integrand(x)
    double precision :: integrand
    double precision, intent(in) :: x
    !-------------------------------------------------------------------------
    integrand = ( exp( - x*x * Fo ) - 1 ) &
      * ( bessel_j0( x * R ) * bessel_y1( x ) - bessel_y0( x * R ) * bessel_j1( x ) )&
      / ( x*x * ( bessel_j1( x ) ** 2 + bessel_y1( x ) ** 2 ) )
  end function integrand
  !----------------------------------------------------------------------
end module ics
