!===============================================================================
! fortran-ballistics-sim
! 2D projectile motion with quadratic air drag.
!
!   Model     : F_d = 1/2 * rho * Cd * A * v^2, opposite to the velocity
!   Method    : fourth-order Runge-Kutta (RK4) with a fixed time step
!   Precision : real64 (iso_fortran_env), Fortran 2008 standard
!
! Build (from the project root):
!   gfortran -Wall -Wextra -Wimplicit-interface -fcheck=all -O2 \
!            src/ballistics.f90 -o ballistics
!
! Usage:
!   ./ballistics            main simulation, writes results/trajectory.csv
!   ./ballistics --verify   physics verification tests (V1 to V3)
!   ./ballistics --study    RK4 convergence study and optimum elevation angle
!   ./ballistics --help     show help
!===============================================================================
module ballistics_core
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private

   integer, parameter, public :: dp = real64
   real(dp), parameter, public :: PI = acos(-1.0_dp)

   ! Component indices of the state vector s = [x, y, vx, vy]
   integer, parameter, public :: IX = 1, IY = 2, IVX = 3, IVY = 4

   !==================== DEFAULT PARAMETERS (edit here) =======================
   ! The values below are used by the main simulation. Example: a ball with a
   ! diameter of 7.4 cm and a mass of 145 g thrown at 50 m/s and 45 degrees.
   type, public :: sim_params
      real(dp) :: v0        = 50.0_dp      ! initial speed               [m/s]
      real(dp) :: theta_deg = 45.0_dp      ! elevation angle             [deg]
      real(dp) :: y0        = 0.0_dp       ! initial height              [m]
      real(dp) :: mass      = 0.145_dp     ! projectile mass             [kg]
      real(dp) :: diameter  = 0.074_dp     ! projectile diameter         [m]
      real(dp) :: cd        = 0.47_dp      ! drag coefficient (sphere)   [-]
      real(dp) :: rho       = 1.225_dp     ! air density (ISA)           [kg/m^3]
      real(dp) :: g         = 9.80665_dp   ! gravitational acceleration  [m/s^2]
      real(dp) :: dt        = 1.0e-3_dp    ! integration time step       [s]
      real(dp) :: t_max     = 600.0_dp     ! safety time limit           [s]
   end type sim_params
   !===========================================================================

   ! Summary of a single simulation run
   type, public :: sim_summary
      real(dp) :: range     = 0.0_dp       ! horizontal distance at y = 0 [m]
      real(dp) :: h_max     = 0.0_dp       ! maximum height              [m]
      real(dp) :: x_apex    = 0.0_dp       ! x position at the apex      [m]
      real(dp) :: t_apex    = 0.0_dp       ! time to reach the apex      [s]
      real(dp) :: t_flight  = 0.0_dp       ! total time of flight        [s]
      real(dp) :: v_impact  = 0.0_dp       ! speed at ground impact      [m/s]
      integer  :: n_points  = 0            ! number of trajectory rows
      logical  :: landed    = .false.      ! .true. once y = 0 is reached
   end type sim_summary

   public :: validate_params, simulate, rk4_step
   public :: vacuum_range, vacuum_height, vacuum_time

contains

   !----------------------------------------------------------------------------
   ! Check the parameters and stop the program if any is physically invalid.
   !----------------------------------------------------------------------------
   subroutine validate_params(p)
      type(sim_params), intent(in) :: p

      if (p%v0 <= 0.0_dp)       error stop 'Invalid parameter: v0 must be > 0'
      if (p%mass <= 0.0_dp)     error stop 'Invalid parameter: mass must be > 0'
      if (p%diameter <= 0.0_dp) error stop 'Invalid parameter: diameter must be > 0'
      if (p%cd < 0.0_dp)        error stop 'Invalid parameter: cd must be >= 0'
      if (p%rho < 0.0_dp)       error stop 'Invalid parameter: rho must be >= 0'
      if (p%g <= 0.0_dp)        error stop 'Invalid parameter: g must be > 0'
      if (p%dt <= 0.0_dp)       error stop 'Invalid parameter: dt must be > 0'
      if (p%t_max <= p%dt)      error stop 'Invalid parameter: t_max must be > dt'
      if (p%y0 < 0.0_dp)        error stop 'Invalid parameter: y0 must be >= 0'
      if (p%theta_deg <= -90.0_dp .or. p%theta_deg > 90.0_dp) &
         error stop 'Invalid parameter: theta_deg must be in (-90, 90]'
      if (p%y0 <= 0.0_dp .and. p%theta_deg <= 0.0_dp) &
         error stop 'Invalid parameter: from the ground (y0 = 0) theta_deg must be > 0'
   end subroutine validate_params

   !----------------------------------------------------------------------------
   ! Right-hand side of the first-order ODE system ds/dt = f(s)
   !   dx/dt  = vx
   !   dy/dt  = vy
   !   dvx/dt = -k |v| vx
   !   dvy/dt = -g - k |v| vy,     with k = rho * Cd * A / (2 m)
   !----------------------------------------------------------------------------
   pure function deriv(s, k, g) result(ds)
      real(dp), intent(in) :: s(4), k, g
      real(dp) :: ds(4)
      real(dp) :: speed

      speed   = hypot(s(IVX), s(IVY))
      ds(IX)  = s(IVX)
      ds(IY)  = s(IVY)
      ds(IVX) = -k * speed * s(IVX)
      ds(IVY) = -g - k * speed * s(IVY)
   end function deriv

   !----------------------------------------------------------------------------
   ! One fourth-order Runge-Kutta step of size h from state s.
   !----------------------------------------------------------------------------
   pure function rk4_step(s, h, k, g) result(s_new)
      real(dp), intent(in) :: s(4), h, k, g
      real(dp) :: s_new(4)
      real(dp) :: k1(4), k2(4), k3(4), k4(4)

      k1 = deriv(s, k, g)
      k2 = deriv(s + 0.5_dp * h * k1, k, g)
      k3 = deriv(s + 0.5_dp * h * k2, k, g)
      k4 = deriv(s + h * k3, k, g)
      s_new = s + (h / 6.0_dp) * (k1 + 2.0_dp * k2 + 2.0_dp * k3 + k4)
   end function rk4_step

   !----------------------------------------------------------------------------
   ! Event location inside one step.
   ! Find h in (0, dt] such that component 'comp' of rk4_step(s, h) = 0,
   ! given that the component is > 0 at h = 0 and <= 0 at h = dt.
   ! Bisection stops when the interval cannot be halved again (machine precision).
   !----------------------------------------------------------------------------
   pure function locate_zero(s, dt, comp, k, g) result(h)
      real(dp), intent(in) :: s(4), dt, k, g
      integer,  intent(in) :: comp
      real(dp) :: h
      real(dp) :: lo, hi, mid, s_mid(4)
      integer  :: iter

      lo = 0.0_dp
      hi = dt
      do iter = 1, 200
         mid = 0.5_dp * (lo + hi)
         if (mid <= lo .or. mid >= hi) exit
         s_mid = rk4_step(s, mid, k, g)
         if (s_mid(comp) > 0.0_dp) then
            lo = mid
         else
            hi = mid
         end if
      end do
      h = hi
   end function locate_zero

   !----------------------------------------------------------------------------
   ! Integrate the trajectory from t = 0 until the projectile hits the ground
   ! (y <= 0). If csv_unit is present, every point (t, x, y, vx, vy) is written
   ! to that unit. The last point is the impact point from event location.
   !----------------------------------------------------------------------------
   subroutine simulate(p, res, csv_unit)
      type(sim_params),  intent(in)  :: p
      type(sim_summary), intent(out) :: res
      integer, optional, intent(in)  :: csv_unit

      real(dp) :: s(4), s_new(4), s_evt(4)
      real(dp) :: k, area, theta, t, h
      integer  :: n, n_max

      area  = 0.25_dp * PI * p%diameter**2
      k     = 0.5_dp * p%rho * p%cd * area / p%mass
      theta = p%theta_deg * PI / 180.0_dp
      n_max = ceiling(p%t_max / p%dt)

      s = [0.0_dp, p%y0, p%v0 * cos(theta), p%v0 * sin(theta)]
      t = 0.0_dp

      res%h_max  = p%y0
      res%x_apex = 0.0_dp
      res%t_apex = 0.0_dp
      call emit(t, s)

      do n = 1, n_max
         s_new = rk4_step(s, p%dt, k, p%g)

         ! Apex: vy changes sign from positive to non-positive
         if (s(IVY) > 0.0_dp .and. s_new(IVY) <= 0.0_dp) then
            h     = locate_zero(s, p%dt, IVY, k, p%g)
            s_evt = rk4_step(s, h, k, p%g)
            if (s_evt(IY) > res%h_max) then
               res%h_max  = s_evt(IY)
               res%x_apex = s_evt(IX)
               res%t_apex = t + h
            end if
         end if

         ! Impact: y crosses zero within this step
         if (s_new(IY) <= 0.0_dp) then
            h          = locate_zero(s, p%dt, IY, k, p%g)
            s_evt      = rk4_step(s, h, k, p%g)
            t          = t + h
            res%landed = .true.
            call emit(t, s_evt)
            res%range    = s_evt(IX)
            res%t_flight = t
            res%v_impact = hypot(s_evt(IVX), s_evt(IVY))
            return
         end if

         s = s_new
         t = real(n, dp) * p%dt          ! avoids error build-up from t = t + dt
         call emit(t, s)
      end do

      ! t_max reached before hitting the ground
      res%landed   = .false.
      res%range    = s(IX)
      res%t_flight = t
      res%v_impact = hypot(s(IVX), s(IVY))

   contains

      subroutine emit(t_now, state)
         real(dp), intent(in) :: t_now, state(4)

         res%n_points = res%n_points + 1
         if (present(csv_unit)) then
            write(csv_unit, '(a)') fmt_real(t_now) // ',' //          &
                                   fmt_real(state(IX)) // ',' //      &
                                   fmt_real(state(IY)) // ',' //      &
                                   fmt_real(state(IVX)) // ',' //     &
                                   fmt_real(state(IVY))
         end if
      end subroutine emit

   end subroutine simulate

   !----------------------------------------------------------------------------
   ! Format a real64 value with 17 significant digits and no leading spaces.
   !----------------------------------------------------------------------------
   function fmt_real(x) result(str)
      real(dp), intent(in) :: x
      character(len=:), allocatable :: str
      character(len=32) :: buf

      write(buf, '(es24.16e3)') x
      str = trim(adjustl(buf))
   end function fmt_real

   !----------------------------------------------------------------------------
   ! Vacuum analytic solution (Cd = 0) for a launch from the ground (y0 = 0).
   !----------------------------------------------------------------------------
   pure function vacuum_range(p) result(r)
      type(sim_params), intent(in) :: p
      real(dp) :: r
      r = p%v0**2 * sin(2.0_dp * p%theta_deg * PI / 180.0_dp) / p%g
   end function vacuum_range

   pure function vacuum_height(p) result(hm)
      type(sim_params), intent(in) :: p
      real(dp) :: hm
      hm = (p%v0 * sin(p%theta_deg * PI / 180.0_dp))**2 / (2.0_dp * p%g)
   end function vacuum_height

   pure function vacuum_time(p) result(tf)
      type(sim_params), intent(in) :: p
      real(dp) :: tf
      tf = 2.0_dp * p%v0 * sin(p%theta_deg * PI / 180.0_dp) / p%g
   end function vacuum_time

end module ballistics_core


!===============================================================================
program ballistics
   use, intrinsic :: iso_fortran_env, only: error_unit
   use ballistics_core
   implicit none

   character(len=*), parameter :: OUT_FILE = 'results/trajectory.csv'
   character(len=64) :: arg

   if (command_argument_count() == 0) then
      call run_simulation()
   else
      call get_command_argument(1, arg)
      select case (trim(arg))
      case ('--verify', '-v')
         call run_verification()
      case ('--study', '-s')
         call run_study()
      case ('--help', '-h')
         call print_usage()
      case default
         write(error_unit, '(a)') 'Unknown argument: ' // trim(arg)
         call print_usage()
         error stop 2
      end select
   end if

contains

   !----------------------------------------------------------------------------
   subroutine print_usage()
      write(*, '(a)') 'Usage: ballistics [--verify | --study | --help]'
      write(*, '(a)') '  (no argument)    main simulation, writes ' // OUT_FILE
      write(*, '(a)') '  --verify, -v     verification tests V1 to V3'
      write(*, '(a)') '  --study,  -s     RK4 convergence study and optimum elevation angle'
      write(*, '(a)') '  --help,   -h     show this help'
   end subroutine print_usage

   !----------------------------------------------------------------------------
   ! Main simulation: write the trajectory CSV and print a summary.
   !----------------------------------------------------------------------------
   subroutine run_simulation()
      type(sim_params)  :: p
      type(sim_summary) :: res
      integer :: u, ios
      character(len=256) :: msg

      call validate_params(p)

      open(newunit=u, file=OUT_FILE, status='replace', action='write', &
           iostat=ios, iomsg=msg)
      if (ios /= 0) then
         write(error_unit, '(a)') 'Cannot open ' // OUT_FILE // ': ' // trim(msg)
         write(error_unit, '(a)') 'Run the program from the project root. The results/ folder must exist.'
         error stop 1
      end if

      write(u, '(a)') 't,x,y,vx,vy'
      call simulate(p, res, u)
      close(u)

      call print_summary(p, res)
   end subroutine run_simulation

   !----------------------------------------------------------------------------
   subroutine print_summary(p, res)
      type(sim_params),  intent(in) :: p
      type(sim_summary), intent(in) :: res
      character(len=*), parameter :: LINE = repeat('=', 58)
      character(len=*), parameter :: SEP  = repeat('-', 58)
      character(len=*), parameter :: FR   = '(2x, a, f14.4, 1x, a)'

      write(*, '(a)') LINE
      write(*, '(a)') '  2D PROJECTILE SIMULATION  |  RK4, real64'
      write(*, '(a)') LINE
      write(*, FR) 'Initial speed       v0    =', p%v0,        'm/s'
      write(*, FR) 'Elevation angle     theta =', p%theta_deg, 'deg'
      write(*, FR) 'Initial height      y0    =', p%y0,        'm'
      write(*, FR) 'Mass                m     =', p%mass,      'kg'
      write(*, FR) 'Diameter            d     =', p%diameter,  'm'
      write(*, FR) 'Drag coefficient    Cd    =', p%cd,        '-'
      write(*, FR) 'Air density         rho   =', p%rho,       'kg/m^3'
      write(*, '(2x, a, es14.4, 1x, a)') 'Time step           dt    =', p%dt, 's'
      write(*, '(a)') SEP
      write(*, FR) 'Range               R     =', res%range,    'm'
      write(*, FR) 'Maximum height      H     =', res%h_max,    'm'
      write(*, FR) 'Time of flight      T     =', res%t_flight, 's'
      write(*, FR) 'Time to apex        t_H   =', res%t_apex,   's'
      write(*, FR) 'Impact speed        v_imp =', res%v_impact, 'm/s'
      write(*, '(2x, a, i14)') 'Number of data points     =', res%n_points
      if (p%y0 <= 0.0_dp) then
         write(*, '(a)') SEP
         write(*, '(a)') '  Vacuum comparison (no air drag):'
         write(*, FR) 'R vacuum                  =', vacuum_range(p),  'm'
         write(*, FR) 'H vacuum                  =', vacuum_height(p), 'm'
         write(*, FR) 'Range reduction by drag   =', &
                      100.0_dp * (1.0_dp - res%range / vacuum_range(p)), '%'
      end if
      write(*, '(a)') SEP
      if (res%landed) then
         write(*, '(a)') '  Trajectory data written to ' // OUT_FILE
      else
         write(*, '(a)') '  WARNING: t_max reached before the projectile hit the ground.'
      end if
      write(*, '(a)') LINE
   end subroutine print_summary

   !----------------------------------------------------------------------------
   ! Physics verification tests.
   !  V1: Cd = 0 compared with the vacuum analytic formulas
   !      R = v0^2 sin(2 theta) / g,  H = v0^2 sin^2(theta) / (2 g),
   !      T = 2 v0 sin(theta) / g     for several angles.
   !  V2: with Cd > 0, R and H must be smaller than the vacuum values.
   !  V3: RK4 convergence, R at dt = 1e-2 s and 1e-3 s must be very close.
   ! The program stops with exit code 1 if any test fails.
   !----------------------------------------------------------------------------
   subroutine run_verification()
      real(dp), parameter :: ANGLES(5) = [15.0_dp, 30.0_dp, 45.0_dp, 60.0_dp, 75.0_dp]
      real(dp), parameter :: TOL_VAC  = 1.0e-10_dp   ! max relative error, V1
      real(dp), parameter :: TOL_CONV = 1.0e-6_dp    ! max relative error, V3
      type(sim_params)  :: p
      type(sim_summary) :: res, res_coarse
      real(dp) :: r_ref, h_ref, t_ref, e_r, e_h, e_t, e_conv
      integer  :: i, n_fail
      logical  :: ok

      n_fail = 0

      write(*, '(a)') repeat('=', 78)
      write(*, '(a)') '  TEST V1: Cd = 0 vs vacuum analytic solution  (v0 = 50 m/s, dt = 1e-3 s)'
      write(*, '(a)') repeat('=', 78)
      write(*, '(a)') '  theta     R_sim [m]      R_ref [m]      err_R      err_H      err_T  status'
      write(*, '(a)') repeat('-', 78)

      p%cd = 0.0_dp
      p%y0 = 0.0_dp
      do i = 1, size(ANGLES)
         p%theta_deg = ANGLES(i)
         call validate_params(p)
         call simulate(p, res)

         r_ref = vacuum_range(p)
         h_ref = vacuum_height(p)
         t_ref = vacuum_time(p)
         e_r = abs(res%range    - r_ref) / r_ref
         e_h = abs(res%h_max    - h_ref) / h_ref
         e_t = abs(res%t_flight - t_ref) / t_ref

         ok = res%landed .and. max(e_r, e_h, e_t) <= TOL_VAC
         if (.not. ok) n_fail = n_fail + 1
         write(*, '(2x, f5.1, 2f15.9, 3es11.2, 2x, a)') &
            p%theta_deg, res%range, r_ref, e_r, e_h, e_t, status_text(ok)
      end do
      write(*, '(a, es9.1)') '  Relative error tolerance: ', TOL_VAC

      ! V2: air drag must reduce both range and height
      write(*, '(a)') repeat('-', 78)
      p = sim_params()
      call simulate(p, res)
      ok = res%landed .and. res%range < vacuum_range(p) .and. res%h_max < vacuum_height(p)
      if (.not. ok) n_fail = n_fail + 1
      write(*, '(a)') '  TEST V2: Cd = 0.47 gives R and H below the vacuum values'
      write(*, '(4x, a, f10.4, a, f10.4, a)') 'R = ', res%range, ' m  <  ', vacuum_range(p), ' m'
      write(*, '(4x, a, f10.4, a, f10.4, a, 2x, a)') 'H = ', res%h_max, ' m  <  ', &
         vacuum_height(p), ' m', status_text(ok)

      ! V3: convergence with respect to the step size (Cd > 0)
      write(*, '(a)') repeat('-', 78)
      p%dt = 1.0e-2_dp
      call simulate(p, res_coarse)
      e_conv = abs(res_coarse%range - res%range) / res%range
      ok = res_coarse%landed .and. e_conv <= TOL_CONV
      if (.not. ok) n_fail = n_fail + 1
      write(*, '(a)') '  TEST V3: RK4 convergence, R(dt = 1e-2 s) vs R(dt = 1e-3 s)'
      write(*, '(4x, a, es10.2, a, es8.1, a, 2x, a)') 'relative difference = ', e_conv, &
         '  (limit', TOL_CONV, ')', status_text(ok)

      write(*, '(a)') repeat('=', 78)
      if (n_fail > 0) then
         write(*, '(a, i0, a)') '  RESULT: ', n_fail, ' test(s) FAILED'
         error stop 1
      end if
      write(*, '(a)') '  RESULT: all tests PASSED'
   end subroutine run_verification

   !----------------------------------------------------------------------------
   ! Numerical study with the default parameters.
   !  1) Convergence: relative error of R, H, T with respect to a reference
   !     solution at dt = 1e-4 s for dt = 0.2 / 2^j (j = 0..4), with empirical
   !     order p = log2( e(dt) / e(dt/2) ).
   !  2) Elevation angle that maximizes the range (sweep from 25 to 45 degrees,
   !     resolution 0.05 degrees).
   !----------------------------------------------------------------------------
   subroutine run_study()
      integer,  parameter :: N_DT = 5
      real(dp), parameter :: DT_REF = 1.0e-4_dp
      type(sim_params)  :: p
      type(sim_summary) :: ref, res
      real(dp) :: dt, k, e(3), e_prev(3), best_r, best_th
      integer  :: j

      k = 0.5_dp * p%rho * p%cd * (0.25_dp * PI * p%diameter**2) / p%mass

      write(*, '(a)') repeat('=', 74)
      write(*, '(a)') '  RK4 CONVERGENCE STUDY  (default parameters, reference dt = 1e-4 s)'
      write(*, '(a)') repeat('=', 74)
      write(*, '(a)') '      dt [s]    steps      err_R      err_H      err_T   p_R   p_H   p_T'
      write(*, '(a)') repeat('-', 74)

      p%dt = DT_REF
      call simulate(p, ref)

      e_prev = 0.0_dp
      do j = 0, N_DT - 1
         dt   = 0.2_dp / 2.0_dp**j
         p%dt = dt
         call simulate(p, res)
         e = [abs(res%range    - ref%range)    / ref%range,  &
              abs(res%h_max    - ref%h_max)    / ref%h_max,  &
              abs(res%t_flight - ref%t_flight) / ref%t_flight]
         if (j == 0) then
            write(*, '(2x, f10.4, i9, 3es11.2)') dt, res%n_points - 1, e
         else
            write(*, '(2x, f10.4, i9, 3es11.2, 3f6.2)') dt, res%n_points - 1, e, &
               log(e_prev / e) / log(2.0_dp)
         end if
         e_prev = e
      end do

      write(*, '(a)') repeat('-', 74)
      write(*, '(a, es12.4, a)') '  Drag constant     k   = ', k, ' 1/m'
      write(*, '(a, f12.4, a)')  '  Terminal speed    v_t = ', sqrt(p%g / k), ' m/s'
      write(*, '(a, f12.4)')     '  Ratio F_d / (m g) at t = 0    = ', k * p%v0**2 / p%g

      p%dt    = 1.0e-3_dp
      best_r  = 0.0_dp
      best_th = p%theta_deg
      do j = 0, 400
         p%theta_deg = 25.0_dp + 0.05_dp * real(j, dp)
         call simulate(p, res)
         if (res%range > best_r) then
            best_r  = res%range
            best_th = p%theta_deg
         end if
      end do
      write(*, '(a, f6.2, a, f10.4, a)') '  Optimum angle     theta* = ', best_th, &
         ' deg  (R = ', best_r, ' m)'
      write(*, '(a)') repeat('=', 74)
   end subroutine run_study

   !----------------------------------------------------------------------------
   pure function status_text(ok) result(txt)
      logical, intent(in) :: ok
      character(len=:), allocatable :: txt
      if (ok) then
         txt = 'PASS'
      else
         txt = 'FAIL'
      end if
   end function status_text

end program ballistics
