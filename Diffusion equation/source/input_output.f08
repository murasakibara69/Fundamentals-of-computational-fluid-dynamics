module input_output

  implicit none

  private

  public :: reading
  public :: writing

contains

  subroutine reading(N, scheme_number, L, A0, AL, lambda, rho, c, alpha, qL, Te, T0, r, time)

    use, non_intrinsic :: kinds, only: wp
    
    implicit none

    integer,  intent(out) :: N, scheme_number
    real(wp), intent(out) :: L, A0, AL, lambda, rho, c, alpha, qL, Te, T0, r, time

    namelist /INPUT_PARAMETERS/ N, scheme_number, L, A0, AL, lambda, rho, c, alpha, qL, Te, T0, r, time

    integer :: input

    open(newunit = input, file = 'conditions.nml', status = 'old', action = 'read')
    read(unit = input, nml = INPUT_PARAMETERS)
    close(unit = input)
    write(unit = *, nml = INPUT_PARAMETERS)
    
  end subroutine reading

  subroutine writing(x, u)

    use, non_intrinsic :: kinds, only: wp

    implicit none
    
    real(wp), intent(in) :: x(:), u(:)

    integer :: i, N, output

    N = size(x)
    
    open(newunit = output, file = "result.dat", action = "write", status = "replace")
    do i = 1, N, 1
      write(unit = output, fmt = *) x(i), u(i)
    end do
    close(unit = output)

  end subroutine writing

end module input_output
