program main

  use, non_intrinsic :: kinds,   only: wp
  use, non_intrinsic :: math,    only: pi
  use, non_intrinsic :: input,   only: reading
  use, non_intrinsic :: schemes, only: scheme_interface => scheme, &
                                       upwind_scheme__positive_c,  &
                                       upwind_scheme__negative_c,  &
                                       Crank_Nicolson_scheme
  use, non_intrinsic :: output,  only: writing

  implicit none

  complex(wp), parameter :: imaginary_unit = (0.0_wp, 1.0_wp)
  real(wp),    parameter :: L = 1.0_wp

  integer  :: scheme_number, m, N
  real(wp) :: c, nu, time

  procedure(scheme_interface), pointer :: scheme

  integer               :: i, iteration
  real(wp)              :: t, delta_t, delta_x, k
  real(wp), allocatable :: x(:), u(:), u_next(:)

  complex(wp) :: G
  real(wp)    :: beta, phi, phi_e, error_phase, new_amplitude

  call reading(scheme_number, m, N, c, nu, time)

  delta_x = L / real(N - 1, wp)
  delta_t = nu * delta_x / c

  k = m * pi
  beta = k * delta_x

  select case (scheme_number)
    case (1)
      if (c > 0) then
        scheme => upwind_scheme__positive_c
        G = 1.0_wp - nu * ( 1.0_wp - cos(beta) ) - imaginary_unit * nu * sin(beta)
      else
        scheme => upwind_scheme__negative_c
        G = 1.0_wp + nu * ( 1.0_wp - cos(beta) ) - imaginary_unit * nu * sin(beta)
      end if
    case (6)
      scheme => Crank_Nicolson_scheme
      G = ( 2.0_wp - imaginary_unit * nu * beta ) / ( 2.0_wp + imaginary_unit * nu * beta )
    case default
      stop 'Invalid scheme_number!'
  end select

  allocate(x(1:N), source = [(delta_x * (i - 1), i = 1, N)])
  allocate(u(1:N), source = sin(k * x))
  allocate(u_next(1:N))

  t = 0.0_wp
  iteration = 0
  do while (t < time)
    u_next = scheme(u, nu)
    u = u_next
    t = t + delta_t
    iteration = iteration + 1
  end do

  deallocate(u_next)

  phi_e = - beta * nu
  phi = atan( G%im / G%re )
  error_phase = iteration * (phi_e - phi)
  new_amplitude = abs(G) ** iteration

  call writing(x, sin(k * (x - c * t)), new_amplitude * sin(k * (x - c * t) - error_phase), u)

  deallocate(u)
  deallocate(x)

end program main
