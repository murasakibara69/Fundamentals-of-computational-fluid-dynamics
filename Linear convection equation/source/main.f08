program main

  use, non_intrinsic :: kinds,   only: wp
  use, non_intrinsic :: input,   only: reading
  use, non_intrinsic :: signals, only: signal_interface => signal, &
                                       signal_setter
  use, non_intrinsic :: schemes, only: scheme_interface => scheme, &
                                       upwind_scheme__positive_c,  &
                                       upwind_scheme__negative_c,  &
                                       Crank_Nicolson_scheme
  use, non_intrinsic :: output,  only: writing

  implicit none

  real(wp), parameter :: L = 1.0_wp

  integer  :: task_number, scheme_number, N 
  real(wp) :: c, nu, time

  procedure(signal_interface), pointer :: signal
  procedure(scheme_interface), pointer :: scheme

  integer               :: i
  real(wp)              :: t, delta_t, delta_x 
  real(wp), allocatable :: x(:), u(:), u_next(:)

  call reading(task_number, scheme_number, N, c, nu, time)

  delta_x = L / real(N - 1, wp)
  delta_t = nu * delta_x / c

  signal => signal_setter(task_number)

  select case (scheme_number)
    case (1)
      if (c > 0) then
        scheme => upwind_scheme__positive_c
      else
        scheme => upwind_scheme__negative_c
      end if
    case (6)
      scheme => Crank_Nicolson_scheme
    case default
      stop 'Invalid scheme_number!'
  end select

  allocate(x(1:N), source = [(delta_x * (i - 1), i = 1, N)])
  allocate(u(1:N), source = signal(x))
  allocate(u_next(1:N))

  t = 0.0_wp
  do while (t < time)
    u_next = scheme(u, nu)
    u = u_next
    t = t + delta_t
  end do

  deallocate(u_next)

  call writing(x, signal(x - c * t), u)

  deallocate(u)
  deallocate(x)

end program main
