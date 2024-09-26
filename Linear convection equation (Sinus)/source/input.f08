module input

  implicit none
  
  private

  public :: reading

contains

  subroutine reading(scheme_number, m, N, c, nu, time)

    use, non_intrinsic :: kinds, only: wp

    implicit none
    
    integer,  intent(out) :: scheme_number, N, m
    real(wp), intent(out) :: c, nu, time

    integer :: file_unit

    namelist /INPUT_PARAMETERS/ scheme_number, m, N, c, nu, time

    open(newunit = file_unit, file = 'conditions.nml', status = 'old', action = 'read')
    read(unit = file_unit, nml = INPUT_PARAMETERS)
    close(unit = file_unit)

  end subroutine reading

end module input
