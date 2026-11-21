!===============================================================================
! fortran-ballistics-sim
! Simulasi gerak peluru 2D dengan hambatan udara kuadratik.
!
!   Model   : F_d = 1/2 * rho * Cd * A * v^2, berlawanan arah kecepatan
!   Metode  : Runge-Kutta orde 4 (RK4) dengan langkah waktu tetap
!   Presisi : real64 (iso_fortran_env), standar Fortran 2008
!
! Kompilasi (dari akar proyek):
!   gfortran -Wall -Wextra -Wimplicit-interface -fcheck=all -O2 \
!            src/ballistics.f90 -o ballistics
!
! Pemakaian:
!   ./ballistics            simulasi utama, tulis results/trajectory.csv
!   ./ballistics --verify   uji verifikasi fisika (Cd = 0 vs solusi vakum)
!   ./ballistics --study    studi konvergensi RK4 dan sudut elevasi optimum
!   ./ballistics --help     tampilkan bantuan
!===============================================================================
module ballistics_core
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private

   integer, parameter, public :: dp = real64
   real(dp), parameter, public :: PI = acos(-1.0_dp)

   ! Indeks komponen vektor keadaan s = [x, y, vx, vy]
   integer, parameter, public :: IX = 1, IY = 2, IVX = 3, IVY = 4

   !=================== PARAMETER AWAL (silakan ubah) =========================
   ! Nilai bawaan di bawah dipakai oleh simulasi utama. Contoh: bola berdiameter
   ! 7,4 cm dan bermassa 145 g yang dilempar 50 m/s pada sudut 45 derajat.
   type, public :: sim_params
      real(dp) :: v0        = 50.0_dp      ! kecepatan awal              [m/s]
      real(dp) :: theta_deg = 45.0_dp      ! sudut elevasi               [derajat]
      real(dp) :: y0        = 0.0_dp       ! ketinggian awal             [m]
      real(dp) :: mass      = 0.145_dp     ! massa peluru                [kg]
      real(dp) :: diameter  = 0.074_dp     ! diameter peluru             [m]
      real(dp) :: cd        = 0.47_dp      ! koefisien hambat (bola)     [-]
      real(dp) :: rho       = 1.225_dp     ! massa jenis udara (ISA)     [kg/m^3]
      real(dp) :: g         = 9.80665_dp   ! percepatan gravitasi        [m/s^2]
      real(dp) :: dt        = 1.0e-3_dp    ! langkah waktu integrasi     [s]
      real(dp) :: t_max     = 600.0_dp     ! batas waktu pengaman        [s]
   end type sim_params
   !===========================================================================

   ! Ringkasan hasil satu kali simulasi
   type, public :: sim_summary
      real(dp) :: range     = 0.0_dp       ! jarak horizontal saat y = 0 [m]
      real(dp) :: h_max     = 0.0_dp       ! tinggi maksimum             [m]
      real(dp) :: x_apex    = 0.0_dp       ! posisi x di titik puncak    [m]
      real(dp) :: t_apex    = 0.0_dp       ! waktu mencapai puncak       [s]
      real(dp) :: t_flight  = 0.0_dp       ! waktu tempuh total          [s]
      real(dp) :: v_impact  = 0.0_dp       ! laju saat menyentuh tanah   [m/s]
      integer  :: n_points  = 0            ! jumlah baris data lintasan
      logical  :: landed    = .false.      ! .true. bila y = 0 tercapai
   end type sim_summary

   public :: validate_params, simulate, rk4_step
   public :: vacuum_range, vacuum_height, vacuum_time

contains

   !----------------------------------------------------------------------------
   ! Periksa parameter; hentikan program bila tidak masuk akal secara fisik.
   !----------------------------------------------------------------------------
   subroutine validate_params(p)
      type(sim_params), intent(in) :: p

      if (p%v0 <= 0.0_dp)       error stop 'Parameter tidak valid: v0 harus > 0'
      if (p%mass <= 0.0_dp)     error stop 'Parameter tidak valid: mass harus > 0'
      if (p%diameter <= 0.0_dp) error stop 'Parameter tidak valid: diameter harus > 0'
      if (p%cd < 0.0_dp)        error stop 'Parameter tidak valid: cd harus >= 0'
      if (p%rho < 0.0_dp)       error stop 'Parameter tidak valid: rho harus >= 0'
      if (p%g <= 0.0_dp)        error stop 'Parameter tidak valid: g harus > 0'
      if (p%dt <= 0.0_dp)       error stop 'Parameter tidak valid: dt harus > 0'
      if (p%t_max <= p%dt)      error stop 'Parameter tidak valid: t_max harus > dt'
      if (p%y0 < 0.0_dp)        error stop 'Parameter tidak valid: y0 harus >= 0'
      if (p%theta_deg <= -90.0_dp .or. p%theta_deg > 90.0_dp) &
         error stop 'Parameter tidak valid: theta_deg harus di (-90, 90]'
      if (p%y0 <= 0.0_dp .and. p%theta_deg <= 0.0_dp) &
         error stop 'Parameter tidak valid: dari tanah (y0 = 0) theta_deg harus > 0'
   end subroutine validate_params

   !----------------------------------------------------------------------------
   ! Ruas kanan sistem PDB orde satu: ds/dt = f(s)
   !   dx/dt  = vx
   !   dy/dt  = vy
   !   dvx/dt = -k |v| vx
   !   dvy/dt = -g - k |v| vy,     dengan k = rho * Cd * A / (2 m)
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
   ! Satu langkah Runge-Kutta orde 4 sebesar h dari keadaan s.
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
   ! Lokasi kejadian (event location) di dalam satu langkah.
   ! Cari h di (0, dt] sehingga komponen 'comp' dari rk4_step(s, h) = 0,
   ! dengan syarat komponen tersebut > 0 di h = 0 dan <= 0 di h = dt.
   ! Bisection berhenti saat interval tidak bisa dibelah lagi (presisi mesin).
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
   ! Integrasi lintasan dari t = 0 hingga peluru menyentuh tanah (y <= 0).
   ! Bila csv_unit diberikan, setiap titik (t, x, y, vx, vy) ditulis ke unit
   ! tersebut; titik terakhir adalah titik tumbukan hasil lokasi kejadian.
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

         ! Titik puncak: vy berganti tanda dari positif ke non-positif
         if (s(IVY) > 0.0_dp .and. s_new(IVY) <= 0.0_dp) then
            h     = locate_zero(s, p%dt, IVY, k, p%g)
            s_evt = rk4_step(s, h, k, p%g)
            if (s_evt(IY) > res%h_max) then
               res%h_max  = s_evt(IY)
               res%x_apex = s_evt(IX)
               res%t_apex = t + h
            end if
         end if

         ! Tumbukan: y melewati nol di dalam langkah ini
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
         t = real(n, dp) * p%dt          ! hindari akumulasi galat t = t + dt
         call emit(t, s)
      end do

      ! Batas t_max tercapai tanpa menyentuh tanah
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
   ! Format bilangan real64 dengan 17 digit signifikan tanpa spasi di depan.
   !----------------------------------------------------------------------------
   function fmt_real(x) result(str)
      real(dp), intent(in) :: x
      character(len=:), allocatable :: str
      character(len=32) :: buf

      write(buf, '(es24.16e3)') x
      str = trim(adjustl(buf))
   end function fmt_real

   !----------------------------------------------------------------------------
   ! Solusi analitik vakum (Cd = 0) untuk peluncuran dari tanah (y0 = 0).
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
         write(error_unit, '(a)') 'Argumen tidak dikenal: ' // trim(arg)
         call print_usage()
         error stop 2
      end select
   end if

contains

   !----------------------------------------------------------------------------
   subroutine print_usage()
      write(*, '(a)') 'Pemakaian: ballistics [--verify | --study | --help]'
      write(*, '(a)') '  (tanpa argumen)  simulasi utama, tulis ' // OUT_FILE
      write(*, '(a)') '  --verify, -v     uji verifikasi Cd = 0 vs solusi analitik vakum'
      write(*, '(a)') '  --study,  -s     studi konvergensi RK4 dan sudut elevasi optimum'
      write(*, '(a)') '  --help,   -h     tampilkan bantuan ini'
   end subroutine print_usage

   !----------------------------------------------------------------------------
   ! Simulasi utama: tulis CSV lintasan dan cetak ringkasan.
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
         write(error_unit, '(a)') 'Gagal membuka ' // OUT_FILE // ': ' // trim(msg)
         write(error_unit, '(a)') 'Jalankan program dari akar proyek; folder results/ harus ada.'
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
      write(*, '(a)') '  SIMULASI GERAK PELURU 2D  |  RK4, real64'
      write(*, '(a)') LINE
      write(*, FR) 'Kecepatan awal      v0    =', p%v0,        'm/s'
      write(*, FR) 'Sudut elevasi       theta =', p%theta_deg, 'derajat'
      write(*, FR) 'Ketinggian awal     y0    =', p%y0,        'm'
      write(*, FR) 'Massa               m     =', p%mass,      'kg'
      write(*, FR) 'Diameter            d     =', p%diameter,  'm'
      write(*, FR) 'Koefisien hambat    Cd    =', p%cd,        '-'
      write(*, FR) 'Massa jenis udara   rho   =', p%rho,       'kg/m^3'
      write(*, '(2x, a, es14.4, 1x, a)') 'Langkah waktu       dt    =', p%dt, 's'
      write(*, '(a)') SEP
      write(*, FR) 'Jarak terjauh       R     =', res%range,    'm'
      write(*, FR) 'Tinggi maksimum     H     =', res%h_max,    'm'
      write(*, FR) 'Waktu tempuh        T     =', res%t_flight, 's'
      write(*, FR) 'Waktu ke puncak     t_H   =', res%t_apex,   's'
      write(*, FR) 'Laju tumbukan       v_imp =', res%v_impact, 'm/s'
      write(*, '(2x, a, i14)') 'Jumlah titik data         =', res%n_points
      if (p%y0 <= 0.0_dp) then
         write(*, '(a)') SEP
         write(*, '(a)') '  Pembanding tanpa hambatan udara (vakum):'
         write(*, FR) 'R vakum                   =', vacuum_range(p),  'm'
         write(*, FR) 'H vakum                   =', vacuum_height(p), 'm'
         write(*, FR) 'Reduksi jarak oleh hambat =', &
                      100.0_dp * (1.0_dp - res%range / vacuum_range(p)), '%'
      end if
      write(*, '(a)') SEP
      if (res%landed) then
         write(*, '(a)') '  Data lintasan disimpan ke ' // OUT_FILE
      else
         write(*, '(a)') '  PERINGATAN: t_max tercapai sebelum peluru menyentuh tanah.'
      end if
      write(*, '(a)') LINE
   end subroutine print_summary

   !----------------------------------------------------------------------------
   ! Uji verifikasi fisika.
   !  Uji 1: Cd = 0 dibandingkan dengan rumus analitik vakum
   !         R = v0^2 sin(2 theta) / g,  H = v0^2 sin^2(theta) / (2 g),
   !         T = 2 v0 sin(theta) / g     untuk beberapa sudut.
   !  Uji 2: dengan Cd > 0, R dan H harus lebih kecil dari nilai vakum.
   !  Uji 3: konvergensi RK4, R pada dt = 1e-2 s dan 1e-3 s harus sangat dekat.
   ! Program berhenti dengan kode keluar 1 bila ada uji yang gagal.
   !----------------------------------------------------------------------------
   subroutine run_verification()
      real(dp), parameter :: ANGLES(5) = [15.0_dp, 30.0_dp, 45.0_dp, 60.0_dp, 75.0_dp]
      real(dp), parameter :: TOL_VAC  = 1.0e-10_dp   ! galat relatif maks. uji 1
      real(dp), parameter :: TOL_CONV = 1.0e-6_dp    ! galat relatif maks. uji 3
      type(sim_params)  :: p
      type(sim_summary) :: res, res_coarse
      real(dp) :: r_ref, h_ref, t_ref, e_r, e_h, e_t, e_conv
      integer  :: i, n_fail
      logical  :: ok

      n_fail = 0

      write(*, '(a)') repeat('=', 78)
      write(*, '(a)') '  UJI 1: Cd = 0 vs solusi analitik vakum  (v0 = 50 m/s, dt = 1e-3 s)'
      write(*, '(a)') repeat('=', 78)
      write(*, '(a)') '  theta     R_sim [m]      R_ref [m]    galat_R    galat_H    galat_T  status'
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
      write(*, '(a, es9.1)') '  Toleransi galat relatif: ', TOL_VAC

      ! Uji 2: hambatan udara harus mengurangi jarak dan tinggi
      write(*, '(a)') repeat('-', 78)
      p = sim_params()
      call simulate(p, res)
      ok = res%landed .and. res%range < vacuum_range(p) .and. res%h_max < vacuum_height(p)
      if (.not. ok) n_fail = n_fail + 1
      write(*, '(a)') '  UJI 2: Cd = 0.47 menghasilkan R dan H lebih kecil dari vakum'
      write(*, '(4x, a, f10.4, a, f10.4, a)') 'R = ', res%range, ' m  <  ', vacuum_range(p), ' m'
      write(*, '(4x, a, f10.4, a, f10.4, a, 2x, a)') 'H = ', res%h_max, ' m  <  ', &
         vacuum_height(p), ' m', status_text(ok)

      ! Uji 3: konvergensi terhadap ukuran langkah (Cd > 0)
      write(*, '(a)') repeat('-', 78)
      p%dt = 1.0e-2_dp
      call simulate(p, res_coarse)
      e_conv = abs(res_coarse%range - res%range) / res%range
      ok = res_coarse%landed .and. e_conv <= TOL_CONV
      if (.not. ok) n_fail = n_fail + 1
      write(*, '(a)') '  UJI 3: konvergensi RK4, R(dt = 1e-2 s) vs R(dt = 1e-3 s)'
      write(*, '(4x, a, es10.2, a, es8.1, a, 2x, a)') 'selisih relatif = ', e_conv, &
         '  (batas', TOL_CONV, ')', status_text(ok)

      write(*, '(a)') repeat('=', 78)
      if (n_fail > 0) then
         write(*, '(a, i0, a)') '  HASIL: ', n_fail, ' uji GAGAL'
         error stop 1
      end if
      write(*, '(a)') '  HASIL: semua uji LULUS'
   end subroutine run_verification

   !----------------------------------------------------------------------------
   ! Studi numerik dengan parameter bawaan.
   !  1) Konvergensi: galat relatif R, H, T terhadap solusi acuan dt = 1e-4 s
   !     untuk dt = 0,2 / 2^j (j = 0..4), dengan orde empiris
   !     p = log2( e(dt) / e(dt/2) ).
   !  2) Sudut elevasi yang memaksimalkan jarak (sapuan 25..45 derajat,
   !     resolusi 0,05 derajat).
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
      write(*, '(a)') '  STUDI KONVERGENSI RK4  (parameter bawaan, acuan dt = 1e-4 s)'
      write(*, '(a)') repeat('=', 74)
      write(*, '(a)') '      dt [s]  langkah    galat_R    galat_H    galat_T   p_R   p_H   p_T'
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
      write(*, '(a, es12.4, a)') '  Konstanta hambat  k   = ', k, ' 1/m'
      write(*, '(a, f12.4, a)')  '  Laju terminal     v_t = ', sqrt(p%g / k), ' m/s'
      write(*, '(a, f12.4)')     '  Rasio F_d / (m g) pada t = 0  = ', k * p%v0**2 / p%g

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
      write(*, '(a, f6.2, a, f10.4, a)') '  Sudut optimum     theta* = ', best_th, &
         ' derajat  (R = ', best_r, ' m)'
      write(*, '(a)') repeat('=', 74)
   end subroutine run_study

   !----------------------------------------------------------------------------
   pure function status_text(ok) result(txt)
      logical, intent(in) :: ok
      character(len=:), allocatable :: txt
      if (ok) then
         txt = 'LULUS'
      else
         txt = 'GAGAL'
      end if
   end function status_text

end program ballistics
