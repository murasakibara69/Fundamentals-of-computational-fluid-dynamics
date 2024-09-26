module schemes

  implicit none

  private

  abstract interface
    
    pure function scheme(u, nu) result(u_next)

      use, non_intrinsic :: kinds, only: wp

      implicit none
      
      real(wp), intent(in) :: u(:), nu

      real(wp) :: u_next(1:size(u))

    end function scheme

  end interface

  public :: scheme
  public :: upwind_scheme__positive_c
  public :: upwind_scheme__negative_c
  public :: Crank_Nicolson_scheme
  
contains

  pure function upwind_scheme__positive_c(u, nu) result(u_next)

    use, non_intrinsic :: kinds, only: wp

    implicit none

    real(wp), intent(in) :: u(:), nu

    real(wp) :: u_next(1:size(u))

    integer  :: i, N

    N = size(u)

    u_next(1) = u(1) - nu * ( u(1) - u(N) )
    do i = 2, N - 1, 1
      u_next(i) = u(i) - nu * ( u(i) - u(i - 1) )
    end do
    u_next(N) = u(N) - nu * ( u(N) - u(N - 1) )
    
  end function upwind_scheme__positive_c

  pure function upwind_scheme__negative_c(u, nu) result(u_next)

    use, non_intrinsic :: kinds, only: wp

    implicit none

    real(wp), intent(in) :: u(:), nu

    real(wp) :: u_next(1:size(u))

    integer  :: i, N

    N = size(u)

    u_next(1) = u(1) - nu * ( u(2) - u(1) )
    do i = 2, N - 1, 1
      u_next(i) = u(i) - nu * ( u(i + 1) - u(i) )
    end do
    u_next(N) = u(N) - nu * ( u(1) - u(N) )
    
  end function upwind_scheme__negative_c

  pure function Crank_Nicolson_scheme(u, nu) result(u_next)

    use, non_intrinsic :: kinds, only: wp
    use, non_intrinsic :: math,  only: periodic_run_through_method

    implicit none

    real(wp), intent(in) :: u(:), nu

    real(wp) :: u_next(1:size(u))

    integer  :: i, N
    real(wp) :: a(1:size(u)), b(1:size(u)), c(1:size(u)), d(1:size(u))

    N = size(u)

    a(1) = 0.25_wp * nu
    b(1) = 1.0_wp
    c(1) = 0.25_wp * nu
    do i = 2, N - 1, 1
      a(i) = - 0.25_wp * nu
      b(i) =   1.0_wp
      c(i) =   0.25_wp * nu
    end do
    a(N) = - 0.25_wp * nu
    b(N) =   1.0_wp
    c(N) = - 0.25_wp * nu

    d(1) = u(1) - 0.25_wp * nu * ( u(2) - u(N) )
    do i = 2, N - 1, 1
      d(i) = u(i) - 0.25_wp * nu * ( u(i + 1) - u(i - 1) )
    end do
    d(N) = u(N) - 0.25_wp * nu * ( u(1) - u(N - 1) )

    u_next = periodic_run_through_method(a, b, c, d)

  end function Crank_Nicolson_scheme

end module schemes
