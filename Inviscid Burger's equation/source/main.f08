program main

  use, non_intrinsic :: kinds,   only: wp
  use, non_intrinsic :: input,   only: reading
  use, non_intrinsic :: schemes, only: scheme_interface => scheme, &
                                       scheme_setter
  use, non_intrinsic :: output,  only: writing

  implicit none

  real(wp), parameter :: L = 1.0_wp

  integer  :: scheme_number, N
  real(wp) :: nu, time

  procedure(scheme_interface), pointer :: scheme

  integer               :: i
  real(wp)              :: t, delta_t, delta_x
  real(wp), allocatable :: x(:), u(:), u_next(:), u_analytical(:)

  real(wp), target  :: u1, u2, u3
  real(wp), pointer :: uL, uR
  real(wp)          :: S

  call reading(scheme_number, N, nu, time)

  scheme => scheme_setter(scheme_number)

  delta_x = L / real(N - 1, wp)

  allocate(x(1:N), source = [(delta_x * (i - 1), i = 1, N)])

  u1 = 0.6_wp
  u2 = 0.2_wp
  u3 = 0.6_wp
  allocate(u(1:N))
  where (x < 0.4_wp)
    u = u1
  elsewhere (x < 0.6_wp)
    u = u2
  elsewhere
    u = u3
  end where

  allocate(u_next(1:N))


  t = 0.0_wp
  do while (t < time)
    delta_t = nu * delta_x / maxval(abs(u))
    u_next = scheme(u, delta_t / delta_x)
    u = u_next
    t = t + delta_t
  end do

  deallocate(u_next)

  allocate(u_analytical(1:N))
  do i = 1, N, 1
    if (x(i) < 0.6_wp) then
      uL => u1
      uR => u2
      S = 0.5_wp * (uL + uR)
      if (x(i) - 0.4_wp - S * t < 0.0_wp) then
        u_analytical(i) = uL
      else
        u_analytical(i) = uR
      end if
    else
      uL => u2
      uR => u3
      if ((x(i) - 0.6_wp) / t <= uL) then
        u_analytical(i) = uL
      elseif ((x(i) - 0.6_wp) / t < uR) then
        u_analytical(i) = (x(i) - 0.6_wp) / t
      else
        u_analytical(i) = uR
      end if
    end if
  end do

  call writing(x, u_analytical, u)

  deallocate(x)
  deallocate(u)

end program main
