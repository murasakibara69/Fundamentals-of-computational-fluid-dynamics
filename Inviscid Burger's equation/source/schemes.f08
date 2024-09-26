module schemes

  implicit none

  private

  abstract interface

    pure function scheme(u, nu) result(u_next)

      use, non_intrinsic :: kinds, only: wp

      real(wp), intent(in) :: u(:), nu

      real(wp) :: u_next(1:size(u))

    end function scheme

  end interface

  public :: scheme
  public :: scheme_setter

contains

  pure function scheme_setter(scheme_number) result(pointer_to_scheme)

    implicit none
    
    integer, intent(in) :: scheme_number

    procedure(scheme), pointer :: pointer_to_scheme

    select case (scheme_number)
      case (1)
        pointer_to_scheme => upwind_scheme
      case (6)
        pointer_to_scheme => Crank_Nicolson_scheme
      case default
        pointer_to_scheme => null()
    end select

  end function scheme_setter

  pure function upwind_scheme(u, nu) result(u_next)

    use, non_intrinsic :: kinds, only: wp

    implicit none
    
    real(wp), intent(in) :: u(:), nu

    real(wp) :: u_next(1:size(u))

    integer  :: i, N

    N = size(u)

    if (u(1) > 0.0_wp) then
      u_next(1) = u(1) - 0.5_wp * nu * ( u(1) ** 2.0_wp - u(N) ** 2.0_wp )
    else
      u_next(1) = u(1) - 0.5_wp * nu * ( u(2) ** 2.0_wp - u(1) ** 2.0_wp )
    end if
    do i = 2, N - 1, 1
      if (u(i) > 0.0_wp) then
        u_next(i) = u(i) - 0.5_wp * nu * ( u(i) ** 2.0_wp - u(i - 1) ** 2.0_wp )
      else
        u_next(i) = u(i) - 0.5_wp * nu * ( u(i + 1) ** 2.0_wp - u(i) ** 2.0_wp )
      end if
    end do
    if (u(N) > 0.0_wp) then
      u_next(N) = u(N) - 0.5_wp * nu * ( u(N) ** 2.0_wp - u(N - 1) ** 2.0_wp )
    else
      u_next(N) = u(N) - 0.5_wp * nu * ( u(1) ** 2.0_wp - u(N) ** 2.0_wp )
    end if

  end function upwind_scheme

  pure function Crank_Nicolson_scheme(u, nu) result(u_next)

    use, non_intrinsic :: kinds, only: wp
    use, non_intrinsic :: math,  only: periodic_run_through_method

    implicit none
    
    real(wp), intent(in) :: u(:), nu

    real(wp) :: u_next(1:size(u))

    integer  :: i, N
    real(wp) :: a(1:size(u)), b(1:size(u)), c(1:size(u)), d(1:size(u))
    ! real(wp) :: A(1:size(u), 1:size(u)), B(1:size(u))

    N = size(u)

    a(1) = 0.25_wp * nu * u(1)
    b(1) = 1.0_wp
    c(1) = 0.25_wp * nu * u(2)
    do i = 2, N - 1, 1
      a(i) = - 0.25_wp * nu * u(i - 1)
      b(i) =   1.0_wp
      c(i) =   0.25_wp * nu * u(i + 1)
    end do
    a(N) = - 0.25_wp * nu * u(N - 1)
    b(N) =   1.0_wp
    c(N) = - 0.25_wp * nu * u(N)

    d = u

    ! B = u
    ! A = 0.0_wp
    !
    ! A(1, 1) =   1.0_wp
    ! A(1, 2) =   0.25_wp * nu * u(2)
    ! A(1, N) = - 0.25_wp * nu * u(N)
    ! do i = 2, N - 1, 1
    !   A(i, i - 1) = - 0.25_wp * nu * u(i - 1)
    !   A(i, i    ) =   1.0_wp
    !   A(i, i + 1) =   0.25_wp * nu * u(i + 1)
    ! end do
    ! A(N, 1    ) =   0.25_wp * nu * u(1)
    ! A(N, N - 1) = - 0.25_wp * nu * u(N - 1)
    ! A(N, N    ) =   1.0_wp
    !
    ! call Gaussian_elimination(A, B, u_next)

    u_next = periodic_run_through_method(a, b, c, d)

  end function Crank_Nicolson_scheme

end module schemes
