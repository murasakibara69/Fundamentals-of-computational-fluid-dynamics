program main

  use, non_intrinsic :: kinds,   only: wp
  use, non_intrinsic :: math,    only: pi
  use, non_intrinsic :: input,   only: reading
  use, non_intrinsic :: schemes, only: scheme_interface => scheme, &
                                       scheme_setter
  use, non_intrinsic :: output,  only: writing

  implicit none

  real(wp), parameter :: L = 1.0_wp

  integer  :: task_number, scheme_number, N, m
  real(wp) :: nu, time

  procedure(scheme_interface), pointer :: scheme

  integer               :: i
  real(wp)              :: t, delta_t, delta_x, k
  real(wp), allocatable :: x(:), u_initial(:), u(:), u_next(:)

  call reading(scheme_number, N, m, nu, time)

  scheme => scheme_setter(scheme_number)

  delta_x = L / real(N - 1, wp)
  k       = m * pi

  allocate(x(1:N),         source = [(delta_x * (i - 1), i = 1, N)])
  allocate(u_initial(1:N), source = sin(k * x))
  allocate(u(1:N),         source = u_initial)
  allocate(u_next(1:N))

  t = 0.0_wp
  do while (t < time)
    delta_t = nu * delta_x / maxval(abs(u))
    u_next = scheme(u, delta_t / delta_x)
    u = u_next
    t = t + delta_t
  end do

  deallocate(u_next)

  call writing(x, u_initial, u)

  deallocate(x)
  deallocate(u_initial)
  deallocate(u)

end program main
