module output

  implicit none

  private

  public :: writing

contains

  subroutine writing(x, u_exact, u)

    use, non_intrinsic :: kinds, only: wp

    implicit none
    
    real(wp), intent(in) :: x(:), u_exact(:), u(:)

    integer :: i, N, file_unit

    N = size(x)

    open(newunit = file_unit, file = 'result.dat', status = 'replace', action = 'write')
    do i = 1, N, 1
      write(unit = file_unit, fmt = *) x(i), u_exact(i), u(i)
    end do
    close(unit = file_unit)

  end subroutine writing

end module output
