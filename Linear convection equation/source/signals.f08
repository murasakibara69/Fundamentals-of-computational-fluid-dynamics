module signals

  implicit none

  private

  abstract interface

    pure function signal(x) result(y)

      use, non_intrinsic :: kinds, only: wp

      implicit none

      real(wp), intent(in) :: x(:)

      real(wp) :: y(1:size(x))

    end function signal

  end interface

  public :: signal
  public :: signal_setter
  
contains

  pure function signal_setter(task_number) result(pointer_to_signal)

    implicit none

    integer, intent(in) :: task_number

    procedure(signal), pointer :: pointer_to_signal

    select case (task_number)
      case (15)
        pointer_to_signal => another_signal
      case default
        pointer_to_signal => null()
    end select
    
  end function signal_setter

  pure function another_signal(x) result(y)

    use, non_intrinsic :: kinds, only: wp

    implicit none

    real(wp), intent(in) :: x(:)

    real(wp) :: y(1:size(x))

    real(wp), parameter :: L = 1.0_wp

    real(wp) :: x_in_range(1:size(x))

    where (x > L)
      x_in_range = x - int(x / L)
    elsewhere (x < 0.0_wp)
      x_in_range = L + x - int(x / L)
    elsewhere
      x_in_range = x
    end where

    where (x_in_range < 0.4_wp)
      y = 0.6_wp
    elsewhere (x_in_range >= 0.4_wp .and. x_in_range < 0.6_wp)
      y = 0.2_wp
    elsewhere
      y = 0.6_wp
    end where

  end function another_signal

end module signals
