program main

  use, non_intrinsic :: kinds,        only: wp
  use, non_intrinsic :: math,         only: pi
  use, non_intrinsic :: input_output, only: reading, &
                                            writing
  use, non_intrinsic :: schemes,      only: scheme, &
                                            scheme_setter

  implicit none

  procedure(scheme), pointer :: p_scheme

  integer  :: N, scheme_number
  real(wp) :: L, A0, AL, lambda, rho, c, alpha, qL, Te, T0, r, time

  integer               :: i, io
  real(wp), allocatable :: x(:), area(:), perimeter(:), u(:), u_next(:)
  real(wp)              :: a, t, dt, dx, B1, C1, D1, AN, BN, DN

  call reading(N, scheme_number, L, A0, AL, lambda, rho, c, alpha, qL, Te, T0, r, time)

  p_scheme => scheme_setter(scheme_number)

  a  = lambda / (rho * c)
  dx = L / real(N - 1, wp)
  dt = r * dx ** 2.0_wp / a

  allocate(x(1:N),         source = [(dx * (i - 1), i = 1, N)])
  allocate(area(1:N),      source = A0 + (AL - A0) * x / L)
  allocate(perimeter(1:N), source = sqrt(4.0_wp * pi * area))
  allocate(u(1:N),         source = T0)
  allocate(u_next(1:N))

  B1 = 1.0_wp
  C1 = - 1.0_wp
  D1 = qL * dx / lambda

  AN = - 1.0_wp
  BN = 1.0_wp + alpha * dx / lambda
  DN = alpha * Te * dx / lambda

  u(1) = ( D1 - C1 * u(2)   ) / B1
  u(N) = ( DN - AN * u(N-1) ) / BN
  t = 0.0_wp
  open(newunit = io, file = 'Monitoring point.dat', status = 'replace', action = 'write')
  do while (t < time)
    write(unit = io, fmt = *) t, u(201)
    u_next = p_scheme(u, area, perimeter * alpha * (Te - u), r, dt, B1, C1, D1, AN, BN, DN)
    u = u_next
    t = t + dt
  end do
  close(unit = io)

  deallocate(area)
  deallocate(perimeter)
  deallocate(u_next)

  call writing(x, u)

  deallocate(x)
  deallocate(u)

end program main
