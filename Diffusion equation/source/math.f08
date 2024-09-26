module math

  use, non_intrinsic :: kinds, only: wp

  implicit none
  
  private

  real(wp), parameter :: pi = acos(- 1.0_wp)

  public :: pi
  public :: tridiagonal_matrix_algorithm

contains

  pure function tridiagonal_matrix_algorithm(a, b, c, d) result(x)

    use, non_intrinsic :: kinds, only: wp

    implicit none
    
    real(wp), intent(in) :: a(:), b(:), c(:), d(:)

    real(wp) :: x(1:size(b))

    integer  :: i, N
    real(wp) :: new_c(1:size(b)-1), new_d(1:size(b))

    N = size(b)

    new_c(1) = c(1) / b(1)
    do i = 2, N - 1, 1
      new_c(i) = c(i) / ( b(i) - a(i) * new_c(i-1) )
    end do

    new_d(1) = d(1) / b(1)
    do i = 2, N, 1
      new_d(i) = ( d(i) - a(i) * new_d(i-1) ) / ( b(i) - a(i) * new_c(i-1) )
    end do

    x(N) = new_d(N)
    do i = N - 1, 1, - 1
      x(i) = new_d(i) - new_c(i) * x(i+1)
    end do

  end function tridiagonal_matrix_algorithm

end module math
