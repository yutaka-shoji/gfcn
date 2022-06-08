module ils
  use ISO_C_binding
  implicit none
  private
  public calc_responsefcn
contains
  !----------------------------------------------------------------------
  subroutine calc_responsefcn(n, Fo, fcn) bind(C, name='calc_responsefcn')
    integer, intent(in) :: n
    double precision, dimension(n), intent(in) :: Fo
    double precision, dimension(n), intent(out) :: fcn
    integer :: j
    interface
      function de1(x)
        double precision :: de1
        double precision :: x
      end function de1
    end interface
    !-------------------------------------------------------------------------
    do j = 1, n
      fcn(j) = de1(1/(4*Fo(j)))
    end do
  end subroutine calc_responsefcn
  !----------------------------------------------------------------------
end module ils
