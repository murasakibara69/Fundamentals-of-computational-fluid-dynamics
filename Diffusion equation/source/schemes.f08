module schemes

  implicit none

  private

  abstract interface
    
    pure function scheme(u, A, rhs, r, dt, B1, C1, D1, AN, BN, DN) result(u_next)

      use, non_intrinsic :: kinds, only: wp

      implicit none
      
      real(wp), intent(in) :: u(:), A(:), rhs(:), r, dt, B1, C1, D1, AN, BN, DN

      real(wp) :: u_next(1:size(u))

    end function scheme

  end interface

  public :: scheme
  public :: scheme_setter

contains

  pure function scheme_setter(scheme_number) result(p_scheme)

    implicit none
    
    integer, intent(in) :: scheme_number

    procedure(scheme), pointer :: p_scheme

    select case (scheme_number)
      case (1)
        p_scheme => central_differencing_scheme
      case (10)
        p_scheme => combined_method
      case default
        p_scheme => null()
    end select

  end function scheme_setter

  pure function central_differencing_scheme(u, A, rhs, r, dt, B1, C1, D1, AN, BN, DN) result(u_next)

    use, non_intrinsic :: kinds, only: wp

    implicit none
    
    real(wp), intent(in) :: u(:), A(:), rhs(:), r, dt, B1, C1, D1, AN, BN, DN

    real(wp) :: u_next(1:size(u))

    integer  :: i, N
    real(wp) :: A_plus, A_minus

    N = size(u)

    do i = 2, N - 1, 1
      A_plus    = 0.5_wp * (A(i+1) + A(i))
      A_minus   = 0.5_wp * (A(i-1) + A(i))
      u_next(i) = ( rhs(i) * dt + A(i) * u(i) + r * (A_plus * (u(i+1) - u(i)) - A_minus * (u(i) - u(i-1))) ) / A(i)
    end do
    u_next(1) = (D1 - C1 * u_next(2)) / B1
    u_next(N) = (DN - AN * u_next(N-1)) / BN

  end function central_differencing_scheme

  pure function combined_method(u, A, rhs, r, dt, B1, C1, D1, AN, BN, DN) result(u_next)

    use, non_intrinsic :: kinds, only: wp
    use, non_intrinsic :: math,  only: tridiagonal_matrix_algorithm

    real(wp), intent(in) :: u(:), A(:), rhs(:), r, dt, B1, C1, D1, AN, BN, DN

    real(wp) :: u_next(1:size(u))

    real(wp), parameter :: teta = 1.0_wp / 3.0_wp

    integer  :: i, N
    real(wp) :: aa(1:size(u)), bb(1:size(u)), cc(1:size(u)), dd(1:size(u)), A_plus, A_minus

    N = size(u)

    aa(1) = 0.0_wp
    bb(1) = B1
    cc(1) = C1
    dd(1) = D1
    do i = 2, N - 1, 1
      A_plus = 0.5_wp * (A(i+1) + A(i))
      A_minus = 0.5_wp * (A(i-1) + A(i))
      aa(i) = - teta * r * A_minus
      bb(i) = A(i) + teta * r * (A_plus + A_minus)
      cc(i) = - teta * r * A_plus
      dd(i) = rhs(i) * dt + A(i) * u(i) + (1.0_wp - teta) * r * (A_plus * (u(i+1) - u(i)) - A_minus * (u(i) - u(i-1)))
    end do
    aa(N) = AN
    bb(N) = BN
    cc(N) = 0.0_wp
    dd(N) = DN

    u_next = tridiagonal_matrix_algorithm(aa, bb, cc, dd)

  end function combined_method

end module schemes
