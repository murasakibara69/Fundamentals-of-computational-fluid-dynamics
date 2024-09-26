module math

  use, non_intrinsic :: kinds, only: wp

  implicit none

  private

  real(wp), parameter :: pi = acos(- 1.0_wp)

  public :: pi
  public :: periodic_run_through_method
  
contains
  
  pure function periodic_run_through_method(a, b, c, d) result(x)

    use, non_intrinsic :: kinds, only: wp

    real(wp), intent(in) :: a(:), b(:), c(:), d(:)

    real(wp) :: x(1:size(b))

    integer  :: i, n
    real(wp) :: new_c(1:size(b)), u(1:size(b)), alpha, beta, omega, temp

    n = size(b)

    x = d

    alpha =   a(1)
    beta  =   c(n)
    omega = - b(1)

    new_c(1) = alpha / (b(1) - omega)
    u(1)     = omega / (b(1) - omega)
    x(1)     = x(1)  / (b(1) - omega)

    do i = 2, n - 1, 1
      temp = 1.0_wp / ( b(i) - a(i) * new_c(i - 1) )
      new_c(i) = temp * c(i)
      u(i) = temp * ( 0.0_wp - a(i) * u(i - 1) )
      x(i) = temp * ( x(i)   - a(i) * x(i - 1) )
    end do

    temp = 1.0_wp / ( b(N) - alpha * beta / omega - beta * new_c(N - 1) )
    u(n) = temp * ( alpha - a(n) * u(n - 1) )
    x(n) = temp * ( x(n)  - a(n) * x(n - 1) )

    do i = n - 1, 1, - 1
      u(i) = u(i) - new_c(i) * u(i + 1)
      x(i) = x(i) - new_c(i) * x(i + 1)
    end do

    temp = (x(1) + x(n) * beta / omega) / (1.0_wp + u(1) + u(n) * beta / omega)

    x = x - temp * u

  end function periodic_run_through_method

end module math
