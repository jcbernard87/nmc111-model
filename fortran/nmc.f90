!> NMC111 half-cell model: Li foil | separator | porous NMC111 cathode, with a choice of
!> particle model: uniform particles, or porous agglomerates of crystals.
!>
!> A self-contained program: edit the input file, run, and it writes Time_Voltage.txt.
!> The equations are described in docs/model.md; parameters in docs/parameters.md.
!> The linear block system of each time step is solved with bandsolver (Newman's BAND).
!>
!>   usage:  nmc_f [input.nml]        (default input file: nmc.nml)
!>
!> mode = 'faithful' reproduces the original research code, including the defects
!> listed in docs/deviations.md; mode = 'corrected' applies the fixes.
!>
!> SPDX-License-Identifier: BSD-3-Clause
program nmc
    use, intrinsic :: iso_fortran_env, only: dp => real64, sp => real32, error_unit
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite, ieee_value, ieee_quiet_nan, ieee_positive_inf
    use bandsolver_kernel, only: band_solve, BAND_OK, PIVOT_LEGACY, SINGULAR_EXACT
    use bandsolver_factor, only: band_factorization, band_factor, band_factor_solve
    implicit none

    integer, parameter :: NV = 4                  ! unknowns per node
    integer, parameter :: IC = 1, IP1 = 2, IP2 = 3, ICS = 4

    ! ---------------- parameters (defaults = original research code) ----------------
    character(len=16) :: particle_model = 'uniform'    ! 'uniform' or 'agglomerate'
    real(dp) :: L_cath_um = 24.0_dp, L_sep = 25.0e-4_dp   ! cathode thickness in um
    real(dp) :: L_cath                                    ! [cm], derived in setup
    integer  :: nj = 43, sep_node = 22
    real(dp) :: eps = 0.5_dp, eps_AM = 0.4_dp, eps_sep = 0.39_dp, tau_sep = 4.8_dp, bruggeman = -0.5_dp
    real(dp) :: D = 2.0e-6_dp, t_plus = 0.25_dp, c_bulk = 1.0e-3_dp, z_plus = 1.0_dp, z_minus = -1.0_dp
    real(dp) :: sigma = -1.0_dp, M = 96.46_dp, rho = 4.6_dp, Q_th = 0.150_dp, R_p = 200.0e-7_dp
    real(dp) :: k_rxn = -1.0_dp                    ! < 0: use the default for the chosen mode
    real(dp) :: alpha_a = 0.5_dp, alpha_c = 0.5_dp, k_Li = 1.0e-6_dp, c_Li_ref = 1.0e-3_dp
    real(dp) :: R = 8.314_dp, T = 298.0_dp, F = 96485.0_dp
    real(dp) :: C_rate = 1.0_dp, phi1_init = 4.1_dp, phi2_init = 0.0_dp, cs_init = 1.0e-5_dp
    real(dp) :: t_max = 36000.0_dp
    integer  :: n_steps = 36000
    real(dp) :: V_min = 2.5_dp, V_max = 4.2_dp        ! corrected-mode cutoffs (D-3)
    real(dp) :: fd_step = 1.0e-6_dp
    real(dp) :: newton_tol = 1.0e-10_dp                ! corrected-mode Newton tolerance (D-7)
    integer  :: newton_max_iter = 25
    character(len=16)  :: mode = 'faithful'
    character(len=2048) :: steps = ''                  ! corrected-mode protocol (docs/protocol.md)
    integer  :: cycles = 1
    real(dp) :: write_interval = 18.0_dp               ! [s], corrected mode
    character(len=256) :: file = 'Time_Voltage.txt'

    namelist /model/ particle_model
    namelist /cell/ L_cath_um, L_sep, nj, sep_node, eps, eps_AM, eps_sep, tau_sep, bruggeman
    real(dp) :: kappa_bg = 1.0e-8_dp                   ! corrected: background (solvent) ionic conductivity [S/cm]
    namelist /electrolyte/ D, t_plus, c_bulk, z_plus, z_minus, kappa_bg
    namelist /active/ sigma, M, rho, Q_th, R_p, k_rxn, alpha_a, alpha_c, k_Li, c_Li_ref
    namelist /constants/ R, T, F
    namelist /operation/ C_rate, phi1_init, phi2_init, cs_init, t_max, n_steps, V_min, V_max
    namelist /numerics/ fd_step, newton_tol, newton_max_iter, mode
    namelist /protocol/ steps, cycles
    namelist /output/ file, write_interval

    ! ---------------- derived quantities and mode switches ----------------
    logical  :: faithful
    real(dp) :: lit36, x_max, spec_a, tortuosity, i_spec, i_app, eps_sep_face, phi1_sign
    logical  :: full_current
    real(dp) :: dcat_s, dan_s, ucat_s, uan_s, dcat_c, dan_c, ucat_c, uan_c
    real(dp) :: mass_area, i_1C, vf_AM   ! vf_AM: active volume fraction
    real(dp), parameter :: THETA_REG = 1.0e-6_dp       ! D-13 regularization threshold

    ! ---------------- agglomerate model (particle_model = 'agglomerate') ----------------
    integer  :: nja = 33                          ! nodes along an agglomerate radius
    real(dp) :: time_mod = 20.0_dp                ! number of steps = 3600*C_rate*time_mod
    real(dp) :: R_agg = 1.0e-4_dp*5.0_dp          ! agglomerate radius [cm]
    real(dp) :: R_xtal = 200.0e-7_dp              ! crystal radius inside agglomerates [cm]
    real(dp) :: eps_agg = 0.2_dp                  ! agglomerate porosity
    real(dp) :: D_agg = -1.0_dp                   ! faithful only: salt diffusivity in agglomerates [cm2/s]
    real(dp) :: tortuosity_e = -1.0_dp            ! faithful only: cathode tortuosity (D_e = eps*D/tortuosity_e)
    real(dp) :: mass_loading = 0.020_dp           ! [g/cm2]
    real(dp) :: percent_active = 0.95_dp
    real(dp) :: mol_vol = 0.0476881609_dp         ! [mol/cm3], a literal in the original
    real(dp) :: c0_init = 1.0e-3_dp
    real(dp) :: tau_agg = -1.0_dp                 ! corrected: agglomerate-pore tortuosity (< 0: eps_agg**-0.5)
    real(dp) :: sigma_agg = -1.0_dp               ! corrected: agglomerate solid conductivity (< 0: sigma)
    real(dp) :: dt_s = 1.0_dp                     ! corrected: time step [s]
    namelist /agglomerate/ nja, time_mod, R_agg, R_xtal, eps_agg, D_agg, tortuosity_e, mass_loading, &
        percent_active, mol_vol, c0_init, tau_agg, sigma_agg, dt_s
    ! corrected agglomerate model (aggc_*)
    logical :: aggc = .false.
    integer :: nl                                  ! cathode volumes, one agglomerate each
    real(dp), allocatable :: cag(:,:,:), cag_start(:,:,:), Aq(:,:,:), Bq(:,:,:), Dq(:,:,:), Gq(:,:)
    real(dp) :: v_agg, s_agg, a_x, x_max_a, csmax_a, sig_a
    real(dp), parameter :: U_REF_A = 3.8637058886774844_dp
    real(dp), parameter :: U_REF_U = 3.8685682447595453_dp
    real(dp), parameter :: AK_U(0:10) = [-0.2018059457910574_dp, 0.1123408808723528_dp, -0.0483699097647364_dp, &
        0.0231624989428732_dp, -0.0377897311905149_dp, -0.3307806975105846_dp, 0.2392976745148739_dp, &
        0.7787126945566982_dp, -0.2599451275866008_dp, -0.5898456896544948_dp, 0.0520147453263591_dp]
    real(dp) :: dplus, dminus                      ! ion diffusivities (corrected mode)
    real(dp), allocatable :: a_area(:), a_g(:), a_gs(:)
    real(dp), parameter :: AK_A(0:11) = [-0.255139064974728_dp, 0.0691287746986728_dp, -0.1178158454270744_dp, &
        -0.0444434841626702_dp, 0.243569591966704_dp, 0.0775338354167729_dp, -1.0934643144519782_dp, &
        -0.8893166395840808_dp, 1.7690915896916977_dp, 1.8213923583001588_dp, -1.2074949744867922_dp, &
        -1.3952076583801158_dp]
    ! working storage of the agglomerate model
    real(dp), allocatable :: cel(:,:), dcel(:,:), ca(:,:,:), dcag(:,:)
    real(dp), allocatable :: Ael(:,:,:), Bel(:,:,:), Del(:,:,:), Gel(:,:)
    real(dp), allocatable :: Aag(:,:,:), Bag(:,:,:), Dag(:,:,:), Gag(:,:), iel(:)
    real(dp), allocatable :: dxa(:), xa(:), aWa(:), aEa(:), bWa(:), bEa(:), AWs(:), AEs(:), dV(:)
    real(dp) :: PIr, pi4, pi43, pi43b, vol_agg, area_agg, v_AM, spec_a_e, spec_a_agg, i_final, i_specific
    real(dp) :: tmx, dt_nom, write_every, i_now, ramp, to_electrons, thr_dep, thr_dep2, thr_x, lit36a
    real(dp) :: h_sep, h_cat, h_c, t_an, diff_e, dcat, dan, u0, ucat, uan, uagg, dt_a, tt
    integer(8) :: nsteps
    character(len=1) :: stt
    logical :: first
    integer :: na, last_write, un

    ! ---------------- protocol (corrected mode) ----------------
    integer, parameter :: K_CC = 1, K_CV = 2, K_REST = 3
    integer :: nstep
    integer, allocatable :: skind(:)
    real(dp), allocatable :: sC(:), sV(:), sT(:), sVmin(:), sVmax(:), sImin(:)   ! sT, sImin < 0: unset

    ! ---------------- mesh and state ----------------
    integer :: s                                  ! interface node index (1-based = sep_node)
    real(dp), allocatable :: dx(:), aW(:), aE(:), bW(:), bE(:)
    real(dp), allocatable :: c(:,:), dc(:,:), A(:,:,:), B(:,:,:), Dm(:,:,:), G(:,:)

    ! ---------------- time loop ----------------
    integer  :: ounit, nsolve
    real(dp) :: time, dt, mAhg
    character(len=1) :: state
    character(len=64) :: exit_reason
    character(len=256) :: input_file

    ! ---------------- read input ----------------
    input_file = 'nmc.nml'
    if (command_argument_count() >= 1) call get_command_argument(1, input_file)
    call read_input(trim(input_file))
    select case (trim(particle_model))
    case ('uniform')
    case ('agglomerate')
        if (trim(mode) == 'faithful') then
            call agg_main()
            stop
        end if
        aggc = .true.
        eps_AM = 0.0_dp            ! the electrode scale has no active material of its own
        if (k_rxn < 0) k_rxn = 2.5e-6_dp       ! before setup(), which would use the uniform default
    case default
        write(error_unit,'(A)') 'particle_model must be ''uniform'' or ''agglomerate'''
        error stop 2
    end select
    call setup()

    open(newunit=ounit, file=trim(file), status='replace', action='write')

    allocate(c(NV,nj), dc(NV,nj), A(NV,NV,nj), B(NV,NV,nj), Dm(NV,NV,nj), G(NV,nj))
    c(IC,:) = c_bulk
    c(IP1,:) = phi1_init
    c(IP2,:) = phi2_init
    c(ICS,:) = cs_init
    dc = 0.0_dp
    nsolve = 0
    time = 0.0_dp
    mAhg = 0.0_dp
    dt = t_max/real(n_steps, dp)
    if (aggc) call aggc_setup()
    if (.not. faithful) then                   ! corrected mode: u = ln(c/c_bulk) and the particles' log-odds
        c(IC,:) = 0.0_dp
        if (.not. aggc) c(ICS,:) = log((cs_init/cs_max())/(1.0_dp - cs_init/cs_max()))
    end if
    exit_reason = 'max_steps'

    if (faithful) then
        call run_faithful()
    else
        call parse_protocol()
        call run_protocol()
    end if

    close(ounit)
    write(*,'(A,G0,A,I0,A)') trim(mode)//' run, C-rate ', C_rate, ': exit '//trim(exit_reason)//' after ', &
        nsolve, ' steps; wrote '//trim(file)

contains

    ! =============================== faithful run ===============================
    subroutine run_faithful()
        !! The original program's constant-current discharge, step for step.
        integer :: it, status
        integer :: last_write                 ! an integer in the original (D-4)
        real(dp) :: write_every
        last_write = 0
        write_every = t_max/n_steps/200
        state = 'D'
        if (C_rate < 0) state = 'C'
        do it = 1, n_steps
            if (it == 1) then
                call write_row(.true.)
            else if ((time - last_write)/3600 >= write_every) then
                call write_row(.false.)
                last_write = int(time - dt)
            else if (it >= n_steps) then
                call write_row(.false.)
            else if (c(IP1,nj) >= 99.0_dp .and. state == 'C') then
                call write_row(.false.)
                exit_reason = 'end_of_charge'
                exit
            else if (ieee_is_nan(dc(IC,1))) then
                call write_row(.false.)
                exit_reason = 'nan'
                exit
            else if (time >= 99.0_dp*3600.0_dp) then
                call write_row(.false.)
                exit_reason = 'max_time'
                exit
            end if
            if (state == 'D') then
                mAhg = mAhg + 1000.0_dp*i_spec*dt/3600.0_dp
            else if (state == 'C') then
                mAhg = mAhg - 1000.0_dp*i_spec*dt/3600.0_dp
            end if
            call assemble(dt)
            call band_solve(NV, nj, A, B, Dm, G, dc, status, pivot=PIVOT_LEGACY)
            if (status /= BAND_OK) dc = ieee_value(1.0_dp, ieee_quiet_nan)
            c = c + dc
            nsolve = it
            if (state == 'R') then
                dt = dt*1.0001_dp
            else
                dt = t_max/real(n_steps, dp)
            end if
            time = time + dt
        end do
    end subroutine run_faithful

    ! =============================== protocol (corrected mode) ===============================
    subroutine parse_protocol()
        !! Parse `steps` (docs/protocol.md) and expand it `cycles` times.
        character(len=len(steps)) :: txt, part
        character(len=64) :: word, key
        integer :: n1, k, pos, semi, cyc, w0, w1, eqp, ios
        real(dp) :: val
        real(dp), allocatable :: tC(:), tV(:), tT(:), tVmin(:), tVmax(:), tImin(:)
        integer, allocatable :: tk(:)
        logical :: hasC, hasV

        txt = adjustl(steps)
        if (len_trim(txt) == 0) then
            allocate(skind(1), sC(1), sV(1), sT(1), sVmin(1), sVmax(1), sImin(1))
            nstep = 1
            skind = K_CC; sC = C_rate; sV = 0.0_dp; sT = -1.0_dp; sVmin = V_min; sVmax = V_max; sImin = -1.0_dp
            return
        end if
        n1 = count_steps(txt)
        allocate(tk(n1), tC(n1), tV(n1), tT(n1), tVmin(n1), tVmax(n1), tImin(n1))
        pos = 1
        k = 0
        do while (pos <= len_trim(txt))
            semi = index(txt(pos:), ';')
            if (semi == 0) then
                part = txt(pos:)
                pos = len_trim(txt) + 1
            else
                part = txt(pos:pos+semi-2)
                pos = pos + semi
            end if
            part = adjustl(part)
            if (len_trim(part) == 0) cycle
            k = k + 1
            tC(k) = 0.0_dp; tV(k) = 0.0_dp; tT(k) = -1.0_dp; tVmin(k) = V_min; tVmax(k) = V_max; tImin(k) = -1.0_dp
            hasC = .false.; hasV = .false.
            w0 = 1
            call next_word(part, w0, w1, word)
            select case (lower(word))
            case ('cc');   tk(k) = K_CC
            case ('cv');   tk(k) = K_CV
            case ('rest'); tk(k) = K_REST
            case default;  call proto_error(k, 'unknown step type '//trim(word))
            end select
            w0 = w1
            do
                call next_word(part, w0, w1, word)
                if (len_trim(word) == 0) exit
                w0 = w1
                eqp = index(word, '=')
                if (eqp == 0) call proto_error(k, 'expected key=value, got '//trim(word))
                key = lower(word(:eqp-1))
                word = word(eqp+1:)
                call fix_exponent(word)
                read(word, *, iostat=ios) val
                if (ios /= 0) call proto_error(k, 'bad number '//trim(word))
                select case (tk(k))
                case (K_CC)
                    select case (key)
                    case ('c');    tC(k) = val; hasC = .true.
                    case ('t');    tT(k) = val
                    case ('vmin'); tVmin(k) = val
                    case ('vmax'); tVmax(k) = val
                    case default;  call proto_error(k, 'cc does not take '//trim(key))
                    end select
                case (K_CV)
                    select case (key)
                    case ('v');    tV(k) = val; hasV = .true.
                    case ('t');    tT(k) = val
                    case ('imin'); tImin(k) = val
                    case default;  call proto_error(k, 'cv does not take '//trim(key))
                    end select
                case (K_REST)
                    if (key /= 't') call proto_error(k, 'rest does not take '//trim(key))
                    tT(k) = val
                end select
                if (key == 't' .and. val <= 0) call proto_error(k, 't must be positive')
            end do
            if (tk(k) == K_CC .and. .not. hasC) call proto_error(k, 'cc needs C=')
            if (tk(k) == K_CV .and. .not. hasV) call proto_error(k, 'cv needs V=')
            if (tk(k) == K_CV .and. tT(k) < 0 .and. tImin(k) < 0) call proto_error(k, 'cv needs t= or Imin= to end')
            if (tk(k) == K_REST .and. tT(k) < 0) call proto_error(k, 'rest needs t=')
        end do
        cyc = max(1, cycles)
        nstep = k*cyc
        allocate(skind(nstep), sC(nstep), sV(nstep), sT(nstep), sVmin(nstep), sVmax(nstep), sImin(nstep))
        skind = [(tk(1:k), n1 = 1, cyc)]
        sC = [(tC(1:k), n1 = 1, cyc)]
        sV = [(tV(1:k), n1 = 1, cyc)]
        sT = [(tT(1:k), n1 = 1, cyc)]
        sVmin = [(tVmin(1:k), n1 = 1, cyc)]
        sVmax = [(tVmax(1:k), n1 = 1, cyc)]
        sImin = [(tImin(1:k), n1 = 1, cyc)]
    end subroutine parse_protocol

    integer function count_steps(txt)
        character(len=*), intent(in) :: txt
        integer :: i
        count_steps = 1
        do i = 1, len_trim(txt)
            if (txt(i:i) == ';') count_steps = count_steps + 1
        end do
    end function count_steps

    subroutine next_word(str, i0, i1, word)
        !! The next blank-separated word of str at or after position i0; i1 is the position after it.
        character(len=*), intent(in) :: str
        integer, intent(in) :: i0
        integer, intent(out) :: i1
        character(len=*), intent(out) :: word
        integer :: i, j
        word = ''
        i = i0
        do while (i <= len_trim(str))
            if (str(i:i) /= ' ') exit
            i = i + 1
        end do
        j = i
        do while (j <= len_trim(str))
            if (str(j:j) == ' ') exit
            j = j + 1
        end do
        if (j > i) word = str(i:j-1)
        i1 = j
    end subroutine next_word

    function lower(str) result(out)
        character(len=*), intent(in) :: str
        character(len=len(str)) :: out
        integer :: i
        out = str
        do i = 1, len(str)
            if (str(i:i) >= 'A' .and. str(i:i) <= 'Z') out(i:i) = achar(iachar(str(i:i)) + 32)
        end do
    end function lower

    subroutine fix_exponent(word)
        !! Accept Fortran d exponents written in lower case as well.
        character(len=*), intent(inout) :: word
        integer :: i
        do i = 1, len_trim(word)
            if (word(i:i) == 'd' .or. word(i:i) == 'D') word(i:i) = 'e'
        end do
    end subroutine fix_exponent

    subroutine proto_error(k, msg)
        integer, intent(in) :: k
        character(len=*), intent(in) :: msg
        write(error_unit,'(A,I0,A)') 'protocol step ', k, ': '//msg
        error stop 2
    end subroutine proto_error

    subroutine run_protocol()
        !! Run the protocol steps in order (docs/protocol.md).
        integer :: k, nsteps_done
        real(dp) :: I, h, h_done, t_step, last_write
        logical :: stopped, ok
        character(len=16) :: why, reason

        I = 0.0_dp
        if (skind(1) == K_CC) I = sC(1)*i_1C
        i_app = I
        call write_row_c(.true., 1)
        last_write = time
        nsteps_done = 0
        reason = ''
        do k = 1, nstep
            t_step = 0.0_dp
            if (skind(k) == K_CC) then
                I = sC(k)*i_1C
            else if (skind(k) == K_REST) then
                I = 0.0_dp
            end if
            do
                h = dt
                if (sT(k) >= 0) h = min(dt, sT(k) - t_step)
                if (skind(k) == K_CV) then
                    call cv_step(h, sV(k), I, ok)
                    h_done = h
                    stopped = sImin(k) >= 0 .and. abs(I) <= sImin(k)*i_1C
                    why = 'current_limit'
                else
                    i_app = I
                    call advance(h, h_done, stopped, ok, skind(k) == K_CC, sVmin(k), sVmax(k))
                    why = 'cutoff_high'
                    if (stopped) then
                        if (cell_voltage() <= sVmin(k)) why = 'cutoff_low'
                    end if
                end if
                i_app = I
                if (.not. ok) then
                    call write_row_c(.false., k)
                    exit_reason = limit_reason()
                    nsolve = nsteps_done
                    return
                end if
                mAhg = mAhg + 1000.0_dp*(I/mass_area)*h_done/3600.0_dp
                time = time + h_done
                t_step = t_step + h_done
                nsteps_done = nsteps_done + 1
                if (any(ieee_is_nan(c)) .or. aggc_nan()) then
                    call write_row_c(.false., k)
                    exit_reason = 'nan'
                    nsolve = nsteps_done
                    return
                end if
                if (stopped .or. (sT(k) >= 0 .and. t_step >= sT(k)*(1.0_dp - 1.0e-12_dp))) then
                    call write_row_c(.false., k)
                    last_write = time
                    if (stopped) then
                        reason = why
                    else
                        reason = 'duration'
                    end if
                    exit
                end if
                if (time - last_write >= write_interval) then
                    call write_row_c(.false., k)
                    last_write = time
                end if
                if (time >= 99.0_dp*3600.0_dp) then
                    call write_row_c(.false., k)
                    exit_reason = 'max_time'
                    nsolve = nsteps_done
                    return
                end if
            end do
        end do
        if (nstep == 1) then
            exit_reason = reason
        else
            exit_reason = 'end_of_protocol'
        end if
        nsolve = nsteps_done
    end subroutine run_protocol

    subroutine cv_step(h, V_set, I, ok)
        !! One constant-voltage time step: find I with V(I) = V_set (see simulate.cv_step in Python).
        real(dp), intent(in) :: h, V_set
        real(dp), intent(inout) :: I
        logical, intent(out) :: ok
        real(dp), parameter :: tol = 1.0e-9_dp
        real(dp) :: c_start(NV,nj), fI, grow, a, b, fa, fb
        logical :: have_a, have_b, good
        integer :: it, side
        c_start = c
        if (aggc) cag_start = cag
        ok = .false.
        call cv_feval(h, V_set, c_start, I, fI, good)
        if (good .and. abs(fI) <= tol) then
            ok = .true.
            return
        end if
        grow = max(abs(I), 1.0e-2_dp*i_1C)
        have_a = .false.; have_b = .false.
        a = 0; b = 0; fa = 0; fb = 0
        do it = 1, 60
            if (fI > 0) then
                a = I; fa = fI; have_a = .true.
                if (have_b) exit
                if (I < 0) then
                    I = 0.0_dp
                else
                    I = I + grow
                end if
            else
                b = I; fb = fI; have_b = .true.
                if (have_a) exit
                if (I > 0) then
                    I = 0.0_dp
                else
                    I = I - grow
                end if
            end if
            grow = grow*2.0_dp
            call cv_feval(h, V_set, c_start, I, fI, good)
            if (good .and. abs(fI) <= tol) then
                ok = .true.
                return
            end if
        end do
        if (.not. (have_a .and. have_b)) then
            c = c_start
            if (aggc) cag = cag_start
            return
        end if
        side = 0
        do it = 1, 200
            if (ieee_is_finite(fa) .and. ieee_is_finite(fb)) then
                I = (a*fb - b*fa)/(fb - fa)
                if (.not. (a < I .and. I < b)) I = 0.5_dp*(a + b)
            else
                I = 0.5_dp*(a + b)
            end if
            call cv_feval(h, V_set, c_start, I, fI, good)
            if (abs(fI) <= tol .or. (b - a) <= 1.0e-14_dp*i_1C) then
                ok = good
                if (.not. good) then
                    c = c_start
                    if (aggc) cag = cag_start
                end if
                return
            end if
            if (fI > 0) then
                a = I; fa = fI
                if (side == 1 .and. ieee_is_finite(fb)) fb = fb*0.5_dp
                side = 1
            else
                b = I; fb = fI
                if (side == -1 .and. ieee_is_finite(fa)) fa = fa*0.5_dp
                side = -1
            end if
        end do
        c = c_start
        if (aggc) cag = cag_start
    end subroutine cv_step

    subroutine cv_feval(h, V_set, c_start, Itry, fval, success)
        !! f(I) = V - V_set after one Newton step from c_start; +/-inf if the step fails.
        real(dp), intent(in) :: h, V_set, c_start(NV,nj), Itry
        real(dp), intent(out) :: fval
        logical, intent(out) :: success
        c = c_start
        if (aggc) cag = cag_start
        i_app = Itry
        call newton_step(h, success)
        if (success) then
            fval = cell_voltage() - V_set
        else
            fval = ieee_value(1.0_dp, ieee_positive_inf)
            if (Itry >= 0) fval = -fval
        end if
    end subroutine cv_feval

    ! =============================== input ===============================
    subroutine read_input(path)
        character(len=*), intent(in) :: path
        integer :: u, ios
        logical :: exists
        inquire(file=path, exist=exists)
        if (.not. exists) then
            write(error_unit,'(A)') 'input file not found: '//path
            error stop 2
        end if
        open(newunit=u, file=path, status='old', action='read')
        rewind(u); read(u, nml=model, iostat=ios);       call check(ios, 'model')
        if (trim(particle_model) == 'agglomerate') call agg_defaults()
        rewind(u); read(u, nml=agglomerate, iostat=ios); call check(ios, 'agglomerate')
        rewind(u); read(u, nml=cell, iostat=ios);        call check(ios, 'cell')
        rewind(u); read(u, nml=electrolyte, iostat=ios); call check(ios, 'electrolyte')
        rewind(u); read(u, nml=active, iostat=ios);      call check(ios, 'active')
        rewind(u); read(u, nml=constants, iostat=ios);   call check(ios, 'constants')
        rewind(u); read(u, nml=operation, iostat=ios);   call check(ios, 'operation')
        rewind(u); read(u, nml=numerics, iostat=ios);    call check(ios, 'numerics')
        rewind(u); read(u, nml=protocol, iostat=ios);    call check(ios, 'protocol')
        rewind(u); read(u, nml=output, iostat=ios);      call check(ios, 'output')
        close(u)
    end subroutine read_input

    subroutine check(ios, group)
        integer, intent(in) :: ios
        character(len=*), intent(in) :: group
        if (ios > 0) then
            write(error_unit,'(A)') 'error reading namelist group &'//group
            error stop 2
        end if
        ! ios < 0: group absent, keep defaults
    end subroutine check

    pure real(dp) function r32(x)
        !! Round to single precision (reproduces the original's single-precision literals).
        real(dp), intent(in) :: x
        r32 = real(real(x, sp), dp)
    end function r32

    ! =============================== setup ===============================
    subroutine setup()
        real(dp) :: h_sep, h_cat, t_an, d_cat, d_an, u_cat, u_an
        integer :: j

        select case (trim(mode))
        case ('faithful')
            faithful = .true.
        case ('corrected')
            faithful = .false.
        case default
            write(error_unit,'(A)') 'mode must be ''faithful'' or ''corrected'''
            error stop 2
        end select

        if (faithful) then
            R = r32(R); c_bulk = r32(c_bulk); Q_th = r32(Q_th); M = r32(M); rho = r32(rho)
            phi1_init = r32(phi1_init); eps_sep = r32(eps_sep); eps_AM = r32(eps_AM); C_rate = r32(C_rate)
            tau_sep = r32(tau_sep)
            ! the original's fitted values are not distributed: a faithful run must supply them
            if (k_rxn < 0 .or. sigma < 0) then
                write(error_unit,'(A)') 'faithful uniform-particle runs need k_rxn and sigma in the input'
                error stop 2
            end if
            lit36 = r32(3.6_dp)
            x_max = r32(0.55_dp)
            eps_sep_face = eps          ! D-2
            phi1_sign = -1.0_dp         ! D-1
            full_current = .false.      ! D-11
        else
            if (k_rxn < 0) k_rxn = 2.5e-6_dp           ! generic defaults (docs/parameters.md)
            if (sigma < 0) sigma = 0.1_dp
            lit36 = 3.6_dp
            x_max = 0.55_dp
            eps_sep_face = eps_sep
            phi1_sign = 1.0_dp
            full_current = .true.
        end if

        L_cath = L_cath_um*1.0e-4_dp     ! the uniform original writes 24 * 1.0d-4
        vf_AM = eps_AM
        spec_a = 3*vf_AM/R_p
        tortuosity = eps**bruggeman
        i_spec = Q_th*C_rate
        i_app = i_spec*L_cath*vf_AM*rho
        mass_area = L_cath*vf_AM*rho
        i_1C = Q_th*mass_area

        ! ion diffusivities and mobilities, then effective values per region
        t_an = 1.0_dp - t_plus
        d_cat = D*(1.0_dp + (t_an/t_plus))/(2.0_dp*t_an/t_plus)
        d_an = d_cat*t_an/t_plus
        u_cat = d_cat/(R*T)
        u_an = d_an/(R*T)
        dplus = d_cat; dminus = d_an
        dcat_s = d_cat/tau_sep;    dan_s = d_an/tau_sep;    ucat_s = u_cat/tau_sep;    uan_s = u_an/tau_sep
        dcat_c = d_cat/tortuosity; dan_c = d_an/tortuosity; ucat_c = u_cat/tortuosity; uan_c = u_an/tortuosity

        s = sep_node
        allocate(dx(nj), aW(nj), aE(nj), bW(nj), bE(nj))
        h_sep = L_sep/real(sep_node - 2, dp)
        h_cat = L_cath/real(nj - sep_node - 1, dp)
        dx = 0.0_dp
        dx(2:s-1) = h_sep
        dx(s+1:nj-1) = h_cat
        aW = 0.0_dp; bW = 0.0_dp; aE = 0.0_dp; bE = 0.0_dp
        do j = 2, nj
            aW(j) = dx(j-1)/(dx(j-1) + dx(j))
            bW(j) = 2.0_dp/(dx(j-1) + dx(j))
        end do
        do j = 1, nj - 1
            aE(j) = dx(j)/(dx(j+1) + dx(j))
            bE(j) = 2.0_dp/(dx(j) + dx(j+1))
        end do
    end subroutine setup

    ! =============================== kinetics ===============================
    real(dp) function cs_max()
        cs_max = x_max*(rho/M)
    end function cs_max

    real(dp) function ocp(cc, cs)
        !! U = U_ref + (RT/F) ln[(c/c_bulk)(1-theta)/theta] + Redlich-Kister sum (docs/model.md section 2).
        real(dp), intent(in) :: cc, cs
        real(dp) :: th
        th = (cs/(rho/M))/x_max
        ocp = u_ref() + R*T/F*log(cc/c_bulk*(1.0_dp - th)/th) + rk_sum(th)
    end function ocp

    real(dp) function u_ref()
        if (faithful) then
            u_ref = r32(3.8685682447595453_dp)
        else
            u_ref = 3.8685682447595453_dp
        end if
    end function u_ref

    real(dp) function rk_sum(th)
        !! Redlich-Kister sum; accumulated in single precision in faithful mode (the original's Vint, D-4).
        real(dp), intent(in) :: th
        real(dp), parameter :: AK(0:10) = [-0.2018059457910574_dp, 0.1123408808723528_dp, -0.0483699097647364_dp, &
            0.0231624989428732_dp, -0.0377897311905149_dp, -0.3307806975105846_dp, 0.2392976745148739_dp, &
            0.7787126945566982_dp, -0.2599451275866008_dp, -0.5898456896544948_dp, 0.0520147453263591_dp]
        real(sp) :: vint
        real(dp) :: vd, ak_k
        integer :: kk
        vint = 0.0_sp
        vd = 0.0_dp
        do kk = 0, 10
            ak_k = AK(kk)
            if (faithful) then
                ak_k = r32(ak_k)
                vint = real(vint + ak_k*((2*th - 1)**(kk + 1) - (2*th*kk*(1 - th))/(2*th - 1)**(1 - kk)), sp)
            else
                vd = vd + ak_k*((2*th - 1)**(kk + 1) - (2*th*kk*(1 - th))/(2*th - 1)**(1 - kk))
            end if
        end do
        if (faithful) then
            rk_sum = real(vint, dp)
        else
            rk_sum = vd
        end if
    end function rk_sum

    subroutine ocp_slopes(cc, cs, du_dc, du_dcs)
        !! dU/dc and dU/dcs (corrected mode).
        real(dp), intent(in) :: cc, cs
        real(dp), intent(out) :: du_dc, du_dcs
        real(dp), parameter :: AK(0:10) = [-0.2018059457910574_dp, 0.1123408808723528_dp, -0.0483699097647364_dp, &
            0.0231624989428732_dp, -0.0377897311905149_dp, -0.3307806975105846_dp, 0.2392976745148739_dp, &
            0.7787126945566982_dp, -0.2599451275866008_dp, -0.5898456896544948_dp, 0.0520147453263591_dp]
        real(dp) :: th, x, drk, term, rtf
        integer :: k
        th = (cs/(rho/M))/x_max
        x = 2*th - 1
        rtf = R*T/F
        drk = 0.0_dp
        do k = 0, 10
            term = 2.0_dp*(2*k + 1)*x**k
            if (k >= 2) term = term - 4.0_dp*k*(k - 1)*th*(1 - th)*x**(k - 2)
            drk = drk + AK(k)*term
        end do
        du_dc = rtf/cc
        du_dcs = (rtf*(-1.0_dp/(1.0_dp - th) - 1.0_dp/th) + drk)/((rho/M)*x_max)
    end subroutine ocp_slopes

    real(dp) function rate(cc, cs, p1, p2)
        !! Butler-Volmer current density per interfacial area [A/cm2], anodic positive.
        real(dp), intent(in) :: cc, cs, p1, p2
        real(dp) :: eta, i0
        eta = p1 - p2 - ocp(cc, cs)
        i0 = F*k_rxn*(cc**alpha_a)*((cs_max() - cs)**alpha_a)*(cs**alpha_c)
        if (faithful) i0 = r32(i0)     ! D-4: implicitly single precision in the original
        rate = i0*(exp(alpha_a*F*eta/(R*T)) - exp(-(alpha_c*F*eta/(R*T))))
    end function rate

    subroutine rate_derivs(cc, cs, p1, p2, i, di)
        !! Rate and finite-difference derivatives w.r.t. (c, phi1, phi2, cs) (D-6).
        real(dp), intent(in) :: cc, cs, p1, p2
        real(dp), intent(out) :: i, di(NV)
        real(dp) :: h
        if (.not. faithful) then
            call rate_derivs_exact(cc, cs, p1, p2, i, di)
            return
        end if
        h = fd_step
        i = rate(cc, cs, p1, p2)
        if (cc <= h) then
            di(IC) = (rate(cc + h, cs, p1, p2) - i)/h
        else
            di(IC) = (rate(cc + h, cs, p1, p2) - rate(cc - h, cs, p1, p2))/(2.0_dp*h)
        end if
        if (cs <= h) then
            di(ICS) = (rate(cc, cs + h, p1, p2) - i)/h
        else
            di(ICS) = (rate(cc, cs + h, p1, p2) - rate(cc, cs - h, p1, p2))/(2.0_dp*h)
        end if
        di(IP1) = (rate(cc, cs, p1 + h, p2) - rate(cc, cs, p1 - h, p2))/(2.0_dp*h)
        di(IP2) = (rate(cc, cs, p1, p2 + h) - rate(cc, cs, p1, p2 - h))/(2.0_dp*h)
    end subroutine rate_derivs

    subroutine power_reg(x, alpha, delta, g, dg)
        !! x**alpha, replaced below delta by a C1 quadratic with g(0) = 0 and a finite slope (D-13).
        real(dp), intent(in) :: x, alpha, delta
        real(dp), intent(out) :: g, dg
        real(dp) :: u
        if (x < delta) then
            u = x/delta
            g = delta**alpha*((2.0_dp - alpha)*u + (alpha - 1.0_dp)*u*u)
            dg = delta**(alpha - 1.0_dp)*((2.0_dp - alpha) + 2.0_dp*(alpha - 1.0_dp)*u)
        else
            g = x**alpha
            dg = alpha*x**(alpha - 1.0_dp)
        end if
    end subroutine power_reg

    subroutine rate_derivs_exact(cc, cs, p1, p2, i, di)
        !! Rate and exact derivatives w.r.t. (c, phi1, phi2, cs); corrected mode (fixes D-6).
        real(dp), intent(in) :: cc, cs, p1, p2
        real(dp), intent(out) :: i, di(NV)
        real(dp) :: rt, aa, bb, eta, i0, ea, ec, di_deta, gv, dgv, gs, dgs, pre, di0, du_dc, du_dcs
        rt = R*T
        aa = alpha_a*F/rt
        bb = alpha_c*F/rt
        eta = p1 - p2 - ocp(cc, cs)
        call ocp_slopes(cc, cs, du_dc, du_dcs)
        call power_reg(cs_max() - cs, alpha_a, THETA_REG*cs_max(), gv, dgv)
        call power_reg(cs, alpha_c, THETA_REG*cs_max(), gs, dgs)
        pre = F*k_rxn*(cc**alpha_a)
        i0 = pre*gv*gs
        di0 = pre*(gs*(-dgv) + gv*dgs)
        ea = exp(aa*eta)
        ec = exp(-bb*eta)
        i = i0*(ea - ec)
        di_deta = i0*(aa*ea + bb*ec)
        di(IC) = alpha_a*i/cc - di_deta*du_dc
        di(ICS) = di0*(ea - ec) - di_deta*du_dcs
        di(IP1) = di_deta
        di(IP2) = -di_deta
    end subroutine rate_derivs_exact

    ! =============================== corrected-mode time step ===============================
    subroutine equilibrate()
        !! Scale every equation by the largest entry of its row in B (the solution is unchanged).
        integer :: j, r
        real(dp) :: sc
        do j = 1, nj
            do r = 1, NV
                sc = maxval(abs(B(r,:,j)))
                if (sc == 0.0_dp) sc = 1.0_dp
                A(r,:,j) = A(r,:,j)/sc
                B(r,:,j) = B(r,:,j)/sc
                Dm(r,:,j) = Dm(r,:,j)/sc
                G(r,j) = G(r,j)/sc
            end do
        end do
    end subroutine equilibrate

    subroutine newton_step(h, ok)
        !! Corrected mode: one backward-Euler step of length h, solved with Newton's method (log variables).
        real(dp), intent(in) :: h
        logical, intent(out) :: ok
        if (aggc) then
            call aggc_newton(h, ok)
        else
            call lu_newton(h, ok)
        end if
    end subroutine newton_step

    logical function converged(upd, prev)
        !! Newton convergence: the scaled update is below newton_tol, or it has stagnated at the
        !! round-off floor (within 1e3*newton_tol and down by less than half since the last iteration).
        real(dp), intent(in) :: upd, prev
        converged = upd <= newton_tol .or. (upd <= 1.0e3_dp*newton_tol .and. upd >= 0.5_dp*prev)
    end function converged

    subroutine advance(dt, t_done, stopped, ok, check, vlo, vhi)
        !! Advance by dt at the current i_app, halving the sub-step on Newton failure. With `check`,
        !! a sub-step that crosses a voltage cutoff by more than 0.1 mV is halved, so the step ends
        !! within 0.1 mV of the cutoff (see simulate.advance in Python).
        real(dp), intent(in) :: dt, vlo, vhi
        logical, intent(in) :: check
        real(dp), intent(out) :: t_done
        logical, intent(out) :: stopped, ok
        real(dp), parameter :: min_dt = 1.0e-10_dp, event_dv = 1.0e-4_dp, event_min_dt = 1.0e-12_dp
        real(dp) :: hh, vv, mg, c_save(NV,nj), c_begin(NV,nj)
        real(dp), allocatable :: cag_save(:,:,:), cag_begin(:,:,:)
        integer, parameter :: max_failures = 200     ! Newton failures allowed within one step
        integer :: failures
        logical :: good
        if (aggc) then
            allocate(cag_save(NV,na,nl), cag_begin(NV,na,nl))
            cag_begin = cag
        else
            allocate(cag_save(NV,0,0), cag_begin(NV,0,0))
        end if
        c_begin = c
        t_done = 0.0_dp
        hh = dt
        failures = 0
        stopped = .false.
        ok = .true.
        do while (t_done < dt)
            hh = min(hh, dt - t_done)
            c_save = c
            if (aggc) cag_save = cag
            call newton_step(hh, good)
            if (.not. good) then
                failures = failures + 1
                if (hh/2 < min_dt .or. failures >= max_failures) then
                    ok = .false.
                    c = c_begin                ! give up: report the state at the start of the step
                    if (aggc) cag = cag_begin
                    return
                end if
                hh = hh/2
                cycle
            end if
            if (check) then
                vv = cell_voltage()
                mg = min(vv - vlo, vhi - vv)
                if (mg < 0.0_dp) then
                    if (mg < -event_dv .and. hh/2 >= event_min_dt) then
                        c = c_save
                        if (aggc) cag = cag_save
                        hh = hh/2
                        cycle
                    end if
                    t_done = t_done + hh
                    stopped = .true.
                    return
                end if
            end if
            t_done = t_done + hh
            hh = 2.0_dp*hh             ! grow back after a success (up to dt, by the min above)
        end do
    end subroutine advance

    ! =============================== assembly ===============================
    subroutine face_coeffs(e, dcat, ucat, dan, uan, cface, gphi2, dd, ff)
        !! Cation-flux (row 1) and ionic-current (row 3) coefficients of one face.
        real(dp), intent(in) :: e, dcat, ucat, dan, uan, cface, gphi2
        real(dp), intent(inout) :: dd(NV,NV), ff(NV,NV)
        real(dp) :: k
        k = z_plus**2*ucat + z_minus**2*uan
        dd(IC,IC)   = -(e*dcat)
        ff(IC,IC)   = -(e*z_plus*ucat*F*gphi2)
        dd(IC,IP2)  = -(e*z_plus*ucat*F*cface)
        dd(IP2,IC)  = -(e*F*(z_plus*dcat + z_minus*dan))
        ff(IP2,IC)  = -(e*F**2*k*gphi2)
        dd(IP2,IP2) = -(e*F**2*k*cface)
    end subroutine face_coeffs

    subroutine assemble(dt)
        real(dp), intent(in) :: dt
        real(dp) :: dW(NV,NV), dE(NV,NV), fW(NV,NV), fE(NV,NV), rj(NV,NV), gg(NV)
        real(dp) :: cW(NV), cE(NV), gW(NV), gE(NV), i, di(NV)
        integer :: j, k

        A = 0.0_dp; B = 0.0_dp; Dm = 0.0_dp; G = 0.0_dp
        do j = 1, nj
            dW = 0.0_dp; dE = 0.0_dp; fW = 0.0_dp; fE = 0.0_dp; rj = 0.0_dp; gg = 0.0_dp
            cW = 0.0_dp; cE = 0.0_dp; gW = 0.0_dp; gE = 0.0_dp
            if (j > 1) then
                cW = aW(j)*c(:,j) + (1.0_dp - aW(j))*c(:,j-1)
                gW = bW(j)*(c(:,j) - c(:,j-1))
            end if
            if (j < nj) then
                cE = aE(j)*c(:,j+1) + (1.0_dp - aE(j))*c(:,j)
                gE = bE(j)*(c(:,j+1) - c(:,j))
            end if
            call rate_derivs(c(IC,j), c(ICS,j), c(IP1,j), c(IP2,j), i, di)

            if (j == 1) then
                ! ---- Li-foil face ----
                dE(IC,IC)  = -(eps_sep_face*dcat_s)
                fE(IC,IC)  = -(eps_sep_face*z_plus*ucat_s*F*gE(IP2))
                dE(IC,IP2) = -(eps_sep_face*z_plus*ucat_s*F*cE(IC))
                gg(IC) = -i_app/F + (dE(IC,IC)*gE(IC) + fE(IC,IC)*cE(IC))
                rj(ICS,ICS) = 0.0_dp - 1.0_dp*vf_AM/dt
                dE(IP1,IP1) = -(1.0_dp - eps_sep_face)*sigma
                gg(IP1) = phi1_sign*dE(IP1,IP1)*gE(IP1)
                rj(IP2,IP2) = 1.0_dp
                gg(IP2) = 0.0_dp - c(IP2,j)

            else if (j < s) then
                ! ---- separator interior ----
                call face_coeffs(eps_sep, dcat_s, ucat_s, dan_s, uan_s, cW(IC), gW(IP2), dW, fW)
                call face_coeffs(eps_sep, dcat_s, ucat_s, dan_s, uan_s, cE(IC), gE(IP2), dE, fE)
                rj(IC,IC) = -(eps_sep/dt*dx(j))
                gg(IC) = 0.0_dp - (fW(IC,IC)*cW(IC) + dW(IC,IC)*gW(IC)) + (fE(IC,IC)*cE(IC) + dE(IC,IC)*gE(IC))
                rj(ICS,ICS) = -(1.0_dp*(1.0_dp - eps_sep)/dt)
                dW(IP1,IP1) = -(1.0_dp - eps_sep)*sigma
                dE(IP1,IP1) = -(1.0_dp - eps_sep)*sigma
                gg(IP1) = 0.0_dp - (fW(IP1,IP1)*cW(IP1) + dW(IP1,IP1)*gW(IP1)) &
                                 + (fE(IP1,IP1)*cE(IP1) + dE(IP1,IP1)*gE(IP1))
                gg(IP2) = 0.0_dp - current(fW, dW, cW, gW) + current(fE, dE, cE, gE)

            else if (j == s) then
                ! ---- separator/cathode interface ----
                call face_coeffs(eps_sep_face, dcat_s, ucat_s, dan_s, uan_s, cW(IC), gW(IP2), dW, fW)
                call face_coeffs(eps, dcat_c, ucat_c, dan_c, uan_c, cE(IC), gE(IP2), dE, fE)
                gg(IC) = 0.0_dp - (fW(IC,IC)*cW(IC) + dW(IC,IC)*gW(IC)) + (fE(IC,IC)*cE(IC) + dE(IC,IC)*gE(IC))
                dW(IP1,IP1) = -(1.0_dp - eps_sep_face)*sigma
                dE(IP1,IP1) = -(1.0_dp - eps)*sigma
                gg(IP1) = 0.0_dp - (fW(IP1,IP1)*cW(IP1) + dW(IP1,IP1)*gW(IP1)) &
                                 + (fE(IP1,IP1)*cE(IP1) + dE(IP1,IP1)*gE(IP1))
                gg(IP2) = 0.0_dp - current(fW, dW, cW, gW) + current(fE, dE, cE, gE)
                call solid_row(i, di, dt, rj, gg)

            else if (j < nj) then
                ! ---- cathode interior ----
                call face_coeffs(eps, dcat_c, ucat_c, dan_c, uan_c, cW(IC), gW(IP2), dW, fW)
                call face_coeffs(eps, dcat_c, ucat_c, dan_c, uan_c, cE(IC), gE(IP2), dE, fE)
                do k = 1, NV
                    rj(IC,k) = (spec_a*di(k)/F)*dx(j)
                    rj(IP1,k) = -((spec_a*di(k))*dx(j))
                    rj(IP2,k) = (spec_a*di(k))*dx(j)
                end do
                rj(IC,IC) = (spec_a*di(IC)/F)*dx(j) - (eps/dt)*dx(j)
                gg(IC) = -((spec_a*i/F)*dx(j)) - (fW(IC,IC)*cW(IC) + dW(IC,IC)*gW(IC)) &
                                               + (fE(IC,IC)*cE(IC) + dE(IC,IC)*gE(IC))
                dW(IP1,IP1) = -(1.0_dp - eps)*sigma
                dE(IP1,IP1) = -(1.0_dp - eps)*sigma
                gg(IP1) = (spec_a*i)*dx(j) - (fW(IP1,IP1)*cW(IP1) + dW(IP1,IP1)*gW(IP1)) &
                                           + (fE(IP1,IP1)*cE(IP1) + dE(IP1,IP1)*gE(IP1))
                gg(IP2) = -((spec_a*i)*dx(j)) - current(fW, dW, cW, gW) + current(fE, dE, cE, gE)
                call solid_row(i, di, dt, rj, gg)

            else
                ! ---- current collector ----
                call face_coeffs(eps, dcat_c, ucat_c, dan_c, uan_c, cW(IC), gW(IP2), dW, fW)
                gg(IC) = 0.0_dp - dW(IC,IC)*gW(IC) - fW(IC,IC)*cW(IC)
                dW(IP1,IP1) = -(1.0_dp - eps)*sigma
                gg(IP1) = i_app - dW(IP1,IP1)*gW(IP1)
                gg(IP2) = 0.0_dp - dW(IP2,IC)*gW(IC) - fW(IP2,IC)*cW(IC)
                call solid_row(i, di, dt, rj, gg)
            end if

            ! control-volume blocks: A dc(j-1) + B dc(j) + D dc(j+1) = G
            if (j == 1) then
                B(:,:,j) = rj - (1.0_dp - aE(j))*fE + bE(j)*dE
                Dm(:,:,j) = -(aE(j)*fE) - bE(j)*dE
            else if (j == nj) then
                A(:,:,j) = (1.0_dp - aW(j))*fW - bW(j)*dW
                B(:,:,j) = rj + bW(j)*dW + aW(j)*fW
            else
                A(:,:,j) = (1.0_dp - aW(j))*fW - bW(j)*dW
                B(:,:,j) = rj + bW(j)*dW + aW(j)*fW - (1.0_dp - aE(j))*fE + bE(j)*dE
                Dm(:,:,j) = -(aE(j)*fE) - bE(j)*dE
            end if
            G(:,j) = gg
        end do
    end subroutine assemble

    real(dp) function current(ff, dd, cf, gf)
        !! Ionic-current face value used in the residual (D-11: faithful mode omits diffusion).
        real(dp), intent(in) :: ff(NV,NV), dd(NV,NV), cf(NV), gf(NV)
        current = ff(IP2,IP2)*cf(IP2) + dd(IP2,IP2)*gf(IP2)
        if (full_current) current = current + dd(IP2,IC)*gf(IC)
    end function current

    subroutine solid_row(i, di, dt, rj, gg)
        !! vf_AM dcs/dt = -a i_n / F  (uniform particles)
        real(dp), intent(in) :: i, di(NV), dt
        real(dp), intent(inout) :: rj(NV,NV), gg(NV)
        integer :: k
        do k = 1, NV
            rj(ICS,k) = -(spec_a*di(k)/F)
        end do
        rj(ICS,ICS) = -(spec_a*di(ICS)/F) - 1.0_dp*vf_AM/dt
        gg(ICS) = +(spec_a*i/F)
    end subroutine solid_row

    ! =============================== output ===============================
    real(dp) function li_eta()
        !! Overpotential of the lithium counter electrode (output only).
        real(dp) :: i0_li, alpha
        i0_li = F*k_Li*(c_foil()**0.5_dp)*(c_Li_ref**0.5_dp)
        alpha = 0.5_dp
        if (.not. faithful) then           ! symmetric Butler-Volmer (D-12); i_app is the present current
            li_eta = -(R*T/(alpha*F))*asinh(i_app/(2.0_dp*i0_li))
        else if (state == 'C') then
            li_eta = 0.5_dp*log(i_app/i0_li)/(alpha*F/(R*T))
        else if (state == 'D') then
            li_eta = -(0.5_dp*log(i_app/i0_li))/(alpha*F/(R*T))
        else
            li_eta = 0.0_dp
        end if
    end function li_eta

    real(dp) function li_nernst()
        !! Corrected mode: Nernst potential of the lithium foil, (RT/F) ln(c/c_Li_ref) [V].
        li_nernst = R*T/F*log(c_foil()/c_Li_ref)
    end function li_nernst

    real(dp) function c_foil()
        !! Electrolyte concentration at the foil (corrected mode stores u = ln(c/c_bulk)).
        if (faithful) then
            c_foil = c(IC,1)
        else
            c_foil = c_bulk*exp(c(IC,1))
        end if
    end function c_foil

    real(dp) function cell_voltage()
        !! Corrected mode: voltage against the lithium foil (0 V). The solver fixes the gauge with
        !! phi2 = 0 at the foil face; the equations depend only on potential differences, so the
        !! foil-referenced potentials are the solved ones minus U_Li + eta_Li (D-16).
        if (faithful) then
            cell_voltage = c(IP1,nj) + li_eta()
        else
            cell_voltage = c(IP1,nj) + li_eta() - li_nernst()
        end if
    end function cell_voltage

    character(len=32) function limit_reason()
        !! The physical limit the state has reached, reported as the exit reason when a step cannot be
        !! solved: electrolyte below 1e-3*c_bulk anywhere, or particles within 1e-3 of full or empty.
        real(dp) :: cmin, thmin, thmax
        cmin = c_bulk*exp(minval(c(IC,:)))
        if (aggc) then
            cmin = min(cmin, c_bulk*exp(minval(cag(IC,:,:))))
            thmin = sigm(minval(cag(ICS,:,:)))
            thmax = sigm(maxval(cag(ICS,:,:)))
        else
            thmin = sigm(minval(c(ICS,s:nj)))
            thmax = sigm(maxval(c(ICS,s:nj)))
        end if
        if (cmin < 1.0e-3_dp*c_bulk) then
            limit_reason = 'electrolyte_depleted'
        else if (thmax > 1.0_dp - 1.0e-3_dp) then
            limit_reason = 'particles_full'
        else if (thmin < 1.0e-3_dp) then
            limit_reason = 'particles_empty'
        else
            limit_reason = 'solver_fail'
        end if
    end function limit_reason

    logical function aggc_nan()
        aggc_nan = .false.
        if (aggc) aggc_nan = any(ieee_is_nan(cag))
    end function aggc_nan

    subroutine write_row_c(header, step)
        !! Corrected-mode output row: the original columns plus the current and the step index.
        logical, intent(in) :: header
        integer, intent(in) :: step
        real(dp) :: i0_li, eta
        character(len=1) :: st
        if (header .and. aggc) then
            write(ounit,'(A5,1X,2(A12,1X),20(A15,1X))') 'State', 'Time', 'Voltage', 'Equivalence', 'Anode_Eta', &
                'anode_exchange_c', 'Edge_c0', 'Current', 'Step', 'Li_Nernst', 'x_front', 'c_collector'
            write(ounit,'(A5,1X,2(A12,1X),20(A15,1X))') 'CDR', 'hours', 'Volts', 'electron_equivs', 'mV', &
                'mA/cm2', 'mol/cm3', 'mA/cm2', '#', 'mV', 'LixNMC', 'mol/cm3'
        else if (header) then
            write(ounit,'(A5,1X,2(A12,1X),20(A15,1X))') 'State', 'Time', 'Voltage', 'Equivalence', 'Anode_Eta', &
                'anode_exchange_c', 'Edge_c0', 'Current', 'Step', 'Li_Nernst'
            write(ounit,'(A5,1X,2(A12,1X),20(A15,1X))') 'CDR', 'hours', 'Volts', 'electron_equivs', 'mV', &
                'mA/cm2', 'mol/cm3', 'mA/cm2', '#', 'mV'
        end if
        st = 'R'
        if (i_app > 0) st = 'D'
        if (i_app < 0) st = 'C'
        i0_li = F*k_Li*(c_foil()**0.5_dp)*(c_Li_ref**0.5_dp)
        eta = li_eta()
        if (aggc) then
            write(ounit,'(A5,1X,2(F12.5,1X),5(ES15.5,1X),I15,3(1X,ES15.5))') st, time/3600.0_dp, cell_voltage(), &
                mAhg*M*3.6_dp/F, eta*1.0e3_dp, i0_li*1.0e3_dp, c_foil(), i_app*1.0e3_dp, step, li_nernst()*1.0e3_dp, &
                csmax_a*sigm(cag(ICS,na,1))/mol_vol, c_bulk*exp(c(IC,nj))
        else
            write(ounit,'(A5,1X,2(F12.5,1X),5(ES15.5,1X),I15,1X,ES15.5)') st, time/3600.0_dp, cell_voltage(), &
                mAhg*M*3.6_dp/F, eta*1.0e3_dp, i0_li*1.0e3_dp, c_foil(), i_app*1.0e3_dp, step, li_nernst()*1.0e3_dp
        end if
    end subroutine write_row_c

    subroutine write_row(header)
        logical, intent(in) :: header
        real(dp) :: equiv, i0_li, eta
        if (header) then
            write(ounit,'(A5,1X,2(A12,1X),20(A15,1X))') 'State', 'Time', 'Voltage', 'Equivalence', 'Anode_Eta', &
                'anode_exchange_c', 'Edge_c0'
            write(ounit,'(A5,1X,2(A12,1X),20(A15,1X))') 'CDR', 'hours', 'Volts', 'electron_equivs', 'mV', &
                'mA/cm2', 'mol/cm3'
        end if
        equiv = mAhg*M*lit36/F
        i0_li = F*k_Li*(c(IC,1)**0.5_dp)*(c_Li_ref**0.5_dp)
        eta = li_eta()
        write(ounit,'(A5,1X,2(F12.5,1X),20(ES15.5,1X))') state, time/real(3600, dp), c(IP1,nj) + eta, equiv, &
            eta*1.0e3_dp, i0_li*1.0e3_dp, c(IC,1)
    end subroutine write_row


    ! ================================================================================
    ! =============================== agglomerate model ===============================
    ! ================================================================================
    ! Porous electrode whose particles are porous spherical agglomerates of uniform crystals
    ! (docs/model.md section 4). Faithful mode reproduces NMC111_agg.f95: single-precision
    ! constants, one linearized solve per step, and the two scales solved one after the other.

    subroutine agg_defaults()
        !! Agglomerate-model defaults: the original template's geometry and materials, with generic
        !! loading, thickness, conductivity and rate constant (docs/parameters.md).
        nj = 75; sep_node = 22; L_sep = 25.0e-4_dp; L_cath_um = 100.0_dp
        eps = 0.4_dp; eps_sep = 0.39_dp; tau_sep = 4.0_dp
        D = 2.89e-6_dp; t_plus = 0.375_dp; c_bulk = 1.0e-3_dp
        sigma = 0.1_dp; M = 96.46_dp; rho = 4.7_dp; Q_th = 0.155_dp
        phi1_init = 4.3_dp; phi2_init = 0.0_dp; cs_init = 1.0e-5_dp
        C_rate = 1.0_dp; t_max = 72000.0_dp
        V_min = 3.0_dp; V_max = 4.4_dp
    end subroutine agg_defaults

    subroutine agg_main()
        !! Faithful agglomerate run (the corrected model is aggc_*).
        integer :: it, j, l

        faithful = .true.
        ! ---- single-precision literals of the original (D-4) ----
        R = r32(R); c_bulk = r32(c_bulk); c0_init = r32(c0_init); eps_sep = r32(eps_sep)
        M = r32(M); rho = r32(rho); mol_vol = r32(mol_vol); Q_th = r32(Q_th); eps_agg = r32(eps_agg)
        mass_loading = r32(mass_loading); percent_active = r32(percent_active); phi1_init = r32(phi1_init)
        eps = r32(eps); sigma = r32(sigma); C_rate = r32(C_rate); time_mod = r32(time_mod)
        ! the original's fitted values are not distributed: a faithful run must supply them
        if (D_agg < 0 .or. k_rxn < 0 .or. tortuosity_e < 0) then
            write(error_unit,'(A)') 'faithful agglomerate runs need D_agg, k_rxn and tortuosity_e in the input'
            error stop 2
        end if
        PIr = r32(3.141592654_dp)
        lit36a = r32(3.6_dp)
        L_cath = real(real(L_cath_um, sp)/10000.0_sp, dp)      ! THICKNESS/10000.0
        tmx = real(real(t_max, sp), dp)                         ! tmax is an implicitly REAL parameter
        nsteps = int(3.6e3_dp*C_rate*time_mod, 8)
        dt_nom = real(real(tmx, sp)/real(nsteps, sp), dp)
        write_every = real(real(tmx, sp)/real(nsteps, sp), dp)
        thr_dep = r32(0.0001_dp); thr_dep2 = r32(0.00001_dp); thr_x = r32(0.55_dp)

        v_AM = percent_active*mass_loading/(rho*L_cath)
        spec_a_e = 3.0_dp*v_AM/R_agg
        spec_a_agg = 3.0_dp*(1.0_dp - eps_agg)/R_xtal
        i_specific = Q_th*C_rate
        i_final = i_specific*L_cath*v_AM*rho
        t_an = 1.0_dp - t_plus
        diff_e = eps*D/tortuosity_e
        dcat = diff_e*(1.0_dp + (t_an/t_plus))/(2.0_dp*t_an/t_plus)
        dan = dcat*t_an/t_plus
        u0 = D/(R*T)
        ucat = dcat/(R*T)
        uan = dan/(R*T)
        uagg = D_agg/(R*T)

        ! ---- electrode grid ----
        s = sep_node
        allocate(dx(nj), aW(nj), aE(nj), bW(nj), bE(nj))
        h_sep = L_sep/real(sep_node - 2, sp)
        h_cat = L_cath/real(nj - sep_node - 1, sp)
        dx = 0.0_dp; dx(2:s-1) = h_sep; dx(s+1:nj-1) = h_cat
        call agg_faces(nj, dx, aW, aE, bW, bE)
        ! ---- agglomerate grid (spherical control volumes) ----
        na = nja
        allocate(dxa(na), xa(na), aWa(na), aEa(na), bWa(na), bEa(na), AWs(na), AEs(na), dV(na))
        h_c = R_agg/real(na - 2, sp)
        xa = 0.0_dp
        do j = 2, na - 1
            xa(j) = h_c*real(j - 1, sp) - h_c/2.0_dp
        end do
        xa(na) = R_agg
        dxa = 0.0_dp; dxa(2:na-1) = h_c
        call agg_faces(na, dxa, aWa, aEa, bWa, bEa)
        pi4 = real(4.0_sp*real(PIr, sp), dp)
        pi43 = real((4.0_sp*real(PIr, sp))/3.0_sp, dp)
        pi43b = real((4.0_sp/3.0_sp)*real(PIr, sp), dp)
        do j = 1, na
            AWs(j) = pi4*(xa(j) - dxa(j)/2.0_dp)**2.0_dp
            AEs(j) = pi4*(xa(j) + dxa(j)/2.0_dp)**2.0_dp
            dV(j) = pi43*((xa(j) + dxa(j)/2.0_dp)**3.0_dp - (xa(j) - dxa(j)/2.0_dp)**3.0_dp)
        end do
        vol_agg = pi43b*R_agg**3
        area_agg = pi4*R_agg**2

        allocate(cel(NV,nj), dcel(NV,nj), ca(NV,na,nj), dcag(NV,na), iel(nj))
        allocate(Ael(NV,NV,nj), Bel(NV,NV,nj), Del(NV,NV,nj), Gel(NV,nj))
        allocate(Aag(NV,NV,na), Bag(NV,NV,na), Dag(NV,NV,na), Gag(NV,na))
        cel(IC,:) = c_bulk; cel(IP1,:) = phi1_init; cel(IP2,:) = 0.0_dp; cel(ICS,:) = cs_init
        ca(IC,:,:) = c0_init; ca(ICS,:,:) = cs_init; ca(IP1,:,:) = phi1_init; ca(IP2,:,:) = 0.0_dp
        dcel = 0.0_dp

        open(newunit=un, file=trim(file), status='replace', action='write')
        tt = 0.0_dp; mAhg = 0.0_dp; last_write = 0     ! last_write_time is uninitialized in the original (D-19)
        dt_a = dt_nom; stt = 'D'; ramp = 1.0_dp; i_now = 0.0_dp
        to_electrons = 1.0_dp/rho*M
        exit_reason = 'max_steps'; nsolve = 0; first = .true.

        do it = 1, int(nsteps)
            if (it == 1) then
                call agg_write(.true.)
            else if ((tt - last_write) >= write_every) then
                call agg_write(.false.)
                last_write = int(tt - dt_a)
            else if (it >= nsteps) then
                call agg_write(.false.)
            end if
            if (cel(IP1,nj) >= 99.0_dp .and. stt == 'C') then
                call agg_write(.false.); exit_reason = 'end_of_charge'; exit
            else if (ca(ICS,na,s)*to_electrons >= thr_x) then
                call agg_write(.false.); exit_reason = 'x_limit'; exit
            else if (ieee_is_nan(dcel(IC,1))) then
                call agg_write(.false.); exit_reason = 'nan'; exit
            else if (tt >= 99.0_dp*3600.0_dp) then
                call agg_write(.false.); exit_reason = 'max_time'; exit
            end if

            if (stt == 'R') then
                dt_a = dt_a*1.0001_dp
            else if (cel(IC,nj) <= thr_dep) then
                dt_a = real(real(tmx, sp)/(real(nsteps, sp)*10.0_sp), dp)
            else if (cel(IC,nj) <= thr_dep2) then
                dt_a = real(real(tmx, sp)/(real(nsteps, sp)*100.0_sp), dp)
            else
                dt_a = dt_nom
            end if
            tt = tt + dt_a

            if (ramp == 1.0_dp) then
                i_now = i_final/50.0_dp
                ramp = 2.0_dp
            else if (ramp == 2.0_dp .and. abs(i_now*1.5_dp) < abs(i_final)) then
                i_now = i_now*1.5_dp
            else
                i_now = i_final
                ramp = 0.0_dp
            end if

            call agg_electrode(dt_a, i_now)
            mAhg = mAhg + 1000.0_dp*i_specific*dt_a/3600.0_dp
            call agg_solve(nj, Ael, Bel, Del, Gel, dcel)
            cel = cel + dcel
            do l = s, nj
                call agg_particle(l, dt_a, dcel(ICS,l)/dt_a)
                call agg_solve(na, Aag, Bag, Dag, Gag, dcag)
                ca(:,:,l) = ca(:,:,l) + dcag
            end do
            nsolve = it
        end do
        close(un)
        write(*,'(A,G0,A,I0,A)') 'faithful agglomerate run, C-rate ', C_rate, ': exit '//trim(exit_reason)// &
            ' after ', nsolve, ' steps; wrote '//trim(file)

    end subroutine agg_main


    subroutine agg_write(header)
        logical, intent(in) :: header
        real(dp) :: u, p1, p2, c0, cs, eta_contact
        if (header) then
            write(un,'(A5,1X,2(A12,1X),20(A15,1X))') 'State', 'Time', 'Voltage', 'mAhg', 'Equivalence', &
                'Solid_Conc', 'current_density', 'iloc', 'Solution_Pot', 'c0', 'cs_edge', 'cs', 'eta_contact', &
                'U', 'eta_rxn', 'i0', 'current', 'ramp_current'
            write(un,'(A5,1X,2(A12,1X),20(A15,1X))') 'CDR', 'hours', 'Volts', 'mAh/g', 'LixNMC', 'LixNMC', &
                'A/cm2', 'A/cm2', 'Volts', 'mol/cm3', 'mol/cm3', 'mol/cm3', 'Volts', 'Volts', 'Volts', 'A/cm2', &
                'A/cm2', ','
        end if
        p1 = cel(IP1,nj); p2 = cel(IP2,nj); c0 = cel(IC,nj); cs = ca(ICS,na,nj)
        u = agg_ocp(c0, cs)
        eta_contact = i_now*0.0_dp
        write(un,'(A5,1X,2(F12.5,1X),20(ES15.5,1X))') stt, tt/real(3600, sp), p1 - eta_contact, mAhg, &
            mAhg*M*lit36a/F, cs*to_electrons, i_final, agg_rate(c0, cs, p1, p2), p2, c0, cs, cel(ICS,nj), &
            eta_contact, u, p1 - p2 - u, agg_ex(c0, cs), i_now, ramp
    end subroutine agg_write

    subroutine agg_electrode(h, I)
        !! Electrode-scale rows (docs/model.md sections 1 and 4); the kinetics use the
        !! agglomerate-surface crystal concentration.
        real(dp), intent(in) :: h, I
        real(dp) :: dW(NV,NV), dE(NV,NV), fW(NV,NV), fE(NV,NV), rj(NV,NV), gg(NV)
        real(dp) :: cW(NV), cE(NV), gW(NV), gE(NV), di(NV), Ds, us, ksep, ksep_b, kcat, a
        integer :: j, k
        Ds = D/tau_sep; us = u0/tau_sep
        ksep = 1.0_dp**2*us + (-1.0_dp)**2*us
        ksep_b = 1.0_dp**2*u0 + (-1.0_dp)**2*u0
        kcat = 1.0_dp**2*ucat + (-1.0_dp)**2*uan
        a = spec_a_e
        Ael = 0.0_dp; Bel = 0.0_dp; Del = 0.0_dp; Gel = 0.0_dp
        do j = 1, nj
            dW = 0.0_dp; dE = 0.0_dp; fW = 0.0_dp; fE = 0.0_dp; rj = 0.0_dp; gg = 0.0_dp
            cW = 0.0_dp; cE = 0.0_dp; gW = 0.0_dp; gE = 0.0_dp
            if (j > 1) then
                cW = aW(j)*cel(:,j) + (1.0_dp - aW(j))*cel(:,j-1)
                gW = bW(j)*(cel(:,j) - cel(:,j-1))
            end if
            if (j < nj) then
                cE = aE(j)*cel(:,j+1) + (1.0_dp - aE(j))*cel(:,j)
                gE = bE(j)*(cel(:,j+1) - cel(:,j))
            end if
            call agg_rate_fd(cel(IC,j), ca(ICS,na,j), cel(IP1,j), cel(IP2,j), iel(j), di)
            if (j == 1) then
                dE(IC,IC) = -(eps_sep*Ds)
                fE(IC,IC) = -(eps_sep*1.0_dp*us*F*gE(IP2))
                dE(IC,IP2) = -(eps_sep*1.0_dp*us*F*cE(IC))
                gg(IC) = -I/F + (dE(IC,IC)*gE(IC) + fE(IC,IC)*cE(IC))
                rj(ICS,ICS) = 0.0_dp - 1.0_dp*(1.0_dp - eps_sep)/h
                dE(IP1,IP1) = -(1.0_dp - eps_sep)*sigma
                gg(IP1) = 0.0_dp - dE(IP1,IP1)*gE(IP1)                       ! D-1
                rj(IP2,IP2) = 1.0_dp
                gg(IP2) = 0.0_dp - cel(IP2,j)
            else if (j < s) then
                dW(IC,IC) = -(eps_sep*Ds); dE(IC,IC) = -(eps_sep*Ds)
                fW(IC,IC) = -(eps_sep*1.0_dp*us*F*gW(IP2)); fE(IC,IC) = -(eps_sep*1.0_dp*us*F*gE(IP2))
                dW(IC,IP2) = -(eps_sep*1.0_dp*us*F*cW(IC)); dE(IC,IP2) = -(eps_sep*1.0_dp*us*F*cE(IC))
                rj(IC,IC) = -(eps_sep/h*dx(j))
                gg(IC) = 0.0_dp - (fW(IC,IC)*cW(IC) + dW(IC,IC)*gW(IC)) + (fE(IC,IC)*cE(IC) + dE(IC,IC)*gE(IC))
                rj(ICS,ICS) = -(1.0_dp*(1.0_dp - eps_sep)/h)
                dW(IP1,IP1) = -(1.0_dp - eps_sep)*sigma; dE(IP1,IP1) = -(1.0_dp - eps_sep)*sigma
                gg(IP1) = 0.0_dp - (fW(IP1,IP1)*cW(IP1) + dW(IP1,IP1)*gW(IP1)) &
                                 + (fE(IP1,IP1)*cE(IP1) + dE(IP1,IP1)*gE(IP1))
                dW(IP2,IC) = -(eps_sep*F*(1.0_dp*Ds + (-1.0_dp)*Ds)); dE(IP2,IC) = dW(IP2,IC)
                fW(IP2,IC) = -((eps_sep/tau_sep)*F**2*ksep_b*gW(IP2))
                fE(IP2,IC) = -((eps_sep/tau_sep)*F**2*ksep_b*gE(IP2))
                dW(IP2,IP2) = -(eps_sep*F**2*ksep*cW(IC)); dE(IP2,IP2) = -(eps_sep*F**2*ksep*cE(IC))
                gg(IP2) = 0.0_dp - (fW(IP2,IP2)*cW(IP2) + dW(IP2,IP2)*gW(IP2)) &
                                 + (fE(IP2,IP2)*cE(IP2) + dE(IP2,IP2)*gE(IP2))
            else if (j == s) then
                dW(IC,IC) = -(eps_sep*Ds); dE(IC,IC) = -dcat
                fW(IC,IC) = -(eps_sep*1.0_dp*us*F*gW(IP2)); fE(IC,IC) = -(1.0_dp*ucat*F*gE(IP2))
                dW(IC,IP2) = -(eps_sep*1.0_dp*us*F*cW(IC)); dE(IC,IP2) = -(1.0_dp*ucat*F*cE(IC))
                gg(IC) = 0.0_dp - (fW(IC,IC)*cW(IC) + dW(IC,IC)*gW(IC)) + (fE(IC,IC)*cE(IC) + dE(IC,IC)*gE(IC))
                dW(IP1,IP1) = -(1.0_dp - eps)*sigma; dE(IP1,IP1) = -(1.0_dp - eps)*sigma
                gg(IP1) = 0.0_dp - (fW(IP1,IP1)*cW(IP1) + dW(IP1,IP1)*gW(IP1)) &
                                 + (fE(IP1,IP1)*cE(IP1) + dE(IP1,IP1)*gE(IP1))
                dW(IP2,IC) = -(eps_sep*F*(1.0_dp*Ds + (-1.0_dp)*Ds)); dE(IP2,IC) = -(F*(1.0_dp*dcat + (-1.0_dp)*dan))
                fW(IP2,IC) = -(eps_sep*F**2*ksep*gW(IP2)); fE(IP2,IC) = -(F**2*kcat*gE(IP2))
                dW(IP2,IP2) = -(eps_sep*F**2*ksep*cW(IC)); dE(IP2,IP2) = -(F**2*kcat*cE(IC))
                gg(IP2) = 0.0_dp - (fW(IP2,IP2)*cW(IP2) + dW(IP2,IP2)*gW(IP2)) &
                                 + (fE(IP2,IP2)*cE(IP2) + dE(IP2,IP2)*gE(IP2))
                call agg_solid(a, di, iel(j), h, rj, gg)
            else if (j < nj) then
                dW(IC,IC) = -dcat; dE(IC,IC) = -dcat
                fW(IC,IC) = -(1.0_dp*ucat*F*gW(IP2)); fE(IC,IC) = -(1.0_dp*ucat*F*gE(IP2))
                dW(IC,IP2) = -(1.0_dp*ucat*F*cW(IC)); dE(IC,IP2) = -(1.0_dp*ucat*F*cE(IC))
                do k = 1, NV
                    rj(IC,k) = (a*di(k)/F)*dx(j)
                    rj(IP1,k) = -((a*di(k))*dx(j))
                    rj(IP2,k) = (a*di(k))*dx(j)
                end do
                rj(IC,IC) = (a*di(IC)/F)*dx(j) - (eps/h)*dx(j)
                gg(IC) = -((a*iel(j)/F)*dx(j)) - (fW(IC,IC)*cW(IC) + dW(IC,IC)*gW(IC)) &
                                             + (fE(IC,IC)*cE(IC) + dE(IC,IC)*gE(IC))
                dW(IP1,IP1) = -(1.0_dp - eps)*sigma; dE(IP1,IP1) = -(1.0_dp - eps)*sigma
                gg(IP1) = (a*iel(j))*dx(j) - (fW(IP1,IP1)*cW(IP1) + dW(IP1,IP1)*gW(IP1)) &
                                          + (fE(IP1,IP1)*cE(IP1) + dE(IP1,IP1)*gE(IP1))
                dW(IP2,IC) = -(F*(1.0_dp*dcat + (-1.0_dp)*dan)); dE(IP2,IC) = dW(IP2,IC)
                fW(IP2,IC) = -(F**2*kcat*gW(IP2)); fE(IP2,IC) = -(F**2*kcat*gE(IP2))
                dW(IP2,IP2) = -(F**2*kcat*cW(IC)); dE(IP2,IP2) = -(F**2*kcat*cE(IC))
                gg(IP2) = -((a*iel(j))*dx(j)) - (fW(IP2,IP2)*cW(IP2) + dW(IP2,IP2)*gW(IP2)) &
                                             + (fE(IP2,IP2)*cE(IP2) + dE(IP2,IP2)*gE(IP2))   ! D-11
                call agg_solid(a, di, iel(j), h, rj, gg)
            else
                dW(IC,IC) = -dcat
                fW(IC,IC) = -(1.0_dp*ucat*F*gW(IP2))
                dW(IC,IP2) = -(1.0_dp*ucat*F*cW(IC))
                gg(IC) = 0.0_dp - dW(IC,IC)*gW(IC) - fW(IC,IC)*cW(IC)
                dW(IP1,IP1) = -(1.0_dp - eps)*sigma
                gg(IP1) = I - dW(IP1,IP1)*gW(IP1)
                dW(IP2,IC) = -(F*(1.0_dp*dcat + (-1.0_dp)*dan))
                fW(IP2,IC) = -(F**2*kcat*gW(IP2))
                dW(IP2,IP2) = -(F**2*kcat*cW(IC))
                gg(IP2) = 0.0_dp - dW(IP2,IC)*gW(IC) - fW(IP2,IC)*cW(IC)
                call agg_solid(a, di, iel(j), h, rj, gg)
            end if
            call agg_blocks(j, nj, rj, dW, dE, fW, fE, aW(j), aE(j), bW(j), bE(j), Ael(:,:,j), Bel(:,:,j), Del(:,:,j))
            Gel(:,j) = gg
        end do
    end subroutine agg_electrode

    subroutine agg_solid(a, di, iv, h, rj, gg)
        !! Solid balance of the agglomerate model's electrode scale: v_AM dcs/dt = -a i_n/F.
        real(dp), intent(in) :: a, di(NV), iv, h
        real(dp), intent(inout) :: rj(NV,NV), gg(NV)
        integer :: kk
        do kk = 1, NV
            rj(ICS,kk) = -(a*di(kk)/F)
        end do
        rj(ICS,ICS) = -(a*di(ICS)/F) - 1.0_dp*v_AM/h
        gg(ICS) = +(a*iv/F)
    end subroutine agg_solid

    subroutine agg_particle(l, h, dcs_dt)
        !! Rows of the agglomerate at electrode node l (spherical control volumes).
        integer, intent(in) :: l
        real(dp), intent(in) :: h, dcs_dt
        real(dp) :: dW(NV,NV), dE(NV,NV), fW(NV,NV), fE(NV,NV), rj(NV,NV), gg(NV)
        real(dp) :: cW(NV), cE(NV), gW(NV), gE(NV), di(NV), i, a, kagg, i_agg, c_spec_agg
        integer :: j, k
        a = spec_a_agg
        kagg = 1.0_dp**2*uagg + (-1.0_dp)**2*uagg
        c_spec_agg = dcs_dt*F/rho
        i_agg = c_spec_agg*vol_agg*(1 - eps_agg)*rho/area_agg
        Aag = 0.0_dp; Bag = 0.0_dp; Dag = 0.0_dp; Gag = 0.0_dp
        do j = 1, na
            dW = 0.0_dp; dE = 0.0_dp; fW = 0.0_dp; fE = 0.0_dp; rj = 0.0_dp; gg = 0.0_dp
            cW = 0.0_dp; cE = 0.0_dp; gW = 0.0_dp; gE = 0.0_dp
            if (j > 1) then
                cW = aWa(j)*ca(:,j,l) + (1.0_dp - aWa(j))*ca(:,j-1,l)
                gW = bWa(j)*(ca(:,j,l) - ca(:,j-1,l))
            end if
            if (j < na) then
                cE = aEa(j)*ca(:,j+1,l) + (1.0_dp - aEa(j))*ca(:,j,l)
                gE = bEa(j)*(ca(:,j+1,l) - ca(:,j,l))
            end if
            call agg_rate_fd(ca(IC,j,l), ca(ICS,j,l), ca(IP1,j,l), ca(IP2,j,l), i, di)
            do k = 1, NV
                rj(ICS,k) = -(a*di(k)/F)
            end do
            gg(ICS) = +(a*i/F)
            if (j == 1) then
                dE(IC,IC) = 1.0_dp
                gg(IC) = 0.0_dp - dE(IC,IC)*gE(IC)                           ! D-1 (center)
                dE(IP1,IP1) = -(1.0_dp - eps_agg)*sigma
                gg(IP1) = 0.0_dp - dE(IP1,IP1)*gE(IP1)
                dE(IP2,IP2) = 1.0_dp
                gg(IP2) = 0.0_dp - dE(IP2,IP2)*gE(IP2)
                rj(ICS,ICS) = -(a*di(ICS)/F) - 1.0_dp*(1.0_dp - eps_agg)/h
            else if (j < na) then
                dW(IC,IC) = AWs(j)*(-(eps_agg*D_agg)); dE(IC,IC) = AEs(j)*(-(eps_agg*D_agg))
                do k = 1, NV
                    rj(IC,k) = (a*di(k)/F)*dV(j)
                    rj(IP1,k) = -((a*di(k))*dV(j))
                    rj(IP2,k) = (a*di(k))*dV(j)
                end do
                rj(IC,IC) = (a*di(IC)/F)*dV(j) - (eps_agg/h)*dV(j)
                gg(IC) = -((a*i/F)*dV(j)) - (fW(IC,IC)*cW(IC) + dW(IC,IC)*gW(IC)) + (fE(IC,IC)*cE(IC) + dE(IC,IC)*gE(IC))
                rj(ICS,ICS) = -(a*di(ICS)/F) - (1.0_dp - eps_agg)/h
                dW(IP1,IP1) = AWs(j)*(-(1.0_dp - eps_agg)*sigma); dE(IP1,IP1) = AEs(j)*(-(1.0_dp - eps_agg)*sigma)
                gg(IP1) = (a*i)*dV(j) - (fW(IP1,IP1)*cW(IP1) + dW(IP1,IP1)*gW(IP1)) &
                                      + (fE(IP1,IP1)*cE(IP1) + dE(IP1,IP1)*gE(IP1))
                dW(IP2,IC) = AWs(j)*(-(eps_agg*F*(1.0_dp*D_agg + (-1.0_dp)*D_agg)))
                dE(IP2,IC) = AEs(j)*(-(eps_agg*F*(1.0_dp*D_agg + (-1.0_dp)*D_agg)))
                fW(IP2,IC) = AWs(j)*(-(eps_agg*F**2*kagg*gW(IP2))); fE(IP2,IC) = AEs(j)*(-(eps_agg*F**2*kagg*gE(IP2)))
                dW(IP2,IP2) = AWs(j)*(-(eps_agg*F**2*kagg*cW(IC))); dE(IP2,IP2) = AEs(j)*(-(eps_agg*F**2*kagg*cE(IC)))
                fW(IP2,IP2) = AWs(j)*0.0_dp; fE(IP2,IP2) = AEs(j)*0.0_dp
                gg(IP2) = -((a*i)*dV(j)) - (fW(IP2,IP2)*cW(IP2) + dW(IP2,IP2)*gW(IP2)) &
                                         + (fE(IP2,IP2)*cE(IP2) + dE(IP2,IP2)*gE(IP2))
            else
                rj(IC,IC) = 1.0_dp
                gg(IC) = c_bulk - ca(IC,j,l)                                 ! D-17
                rj(ICS,ICS) = -(a*di(ICS)/F) - (1.0_dp - eps_agg)/h
                dW(IP1,IP1) = -(1.0_dp - eps_agg)*sigma
                gg(IP1) = i_agg - dW(IP1,IP1)*gW(IP1)
                rj(IP2,IP2) = 1.0_dp
                gg(IP2) = cel(IP2,l) - ca(IP2,j,l)
            end if
            call agg_blocks(j, na, rj, dW, dE, fW, fE, aWa(j), aEa(j), bWa(j), bEa(j), Aag(:,:,j), Bag(:,:,j), Dag(:,:,j))
            Gag(:,j) = gg
        end do
    end subroutine agg_particle

    real(dp) function agg_ocp(cc, cs)
        !! U = U_ref + (RT/F) ln[(c/c_bulk)(1-theta)/theta] + 12-term Redlich-Kister sum, theta = x/x_max.
        real(dp), intent(in) :: cc, cs
        real(dp), parameter :: AK(0:11) = [-0.255139064974728_dp, 0.0691287746986728_dp, -0.1178158454270744_dp, &
            -0.0444434841626702_dp, 0.243569591966704_dp, 0.0775338354167729_dp, -1.0934643144519782_dp, &
            -0.8893166395840808_dp, 1.7690915896916977_dp, 1.8213923583001588_dp, -1.2074949744867922_dp, &
            -1.3952076583801158_dp]
        real(dp) :: th, xish_max
        real(sp) :: vint
        integer :: kk
        xish_max = M*Q_th*3600/F
        th = (cs/mol_vol)/xish_max
        vint = 0.0_sp
        do kk = 0, 11
            vint = real(vint + r32(AK(kk))*((2*th - 1)**(kk + 1) - (2*th*kk*(1 - th))/(2*th - 1)**(1 - kk)), sp)
        end do
        agg_ocp = r32(3.8637058886774844_dp) + R*T/F*log(cc/c_bulk*(1.0_dp - th)/th) + vint
    end function agg_ocp

    real(dp) function agg_ex(cc, cs)
        real(dp), intent(in) :: cc, cs
        real(dp) :: cimax
        cimax = mol_vol*M*Q_th*3600/F
        agg_ex = F*k_rxn*(cc**alpha_a)*((cimax - cs)**alpha_a)*(cs**alpha_c)
    end function agg_ex

    real(dp) function agg_rate(cc, cs, p1, p2)
        real(dp), intent(in) :: cc, cs, p1, p2
        real(dp) :: eta, ex
        eta = p1 - p2 - agg_ocp(cc, cs)
        ex = r32(agg_ex(cc, cs))                                          ! ex_curr is REAL (D-4)
        agg_rate = ex*(exp(alpha_a*F*eta/(R*T)) - exp(-(alpha_c*F*eta/(R*T))))
    end function agg_rate

    subroutine agg_rate_fd(cc, cs, p1, p2, i, di)
        real(dp), intent(in) :: cc, cs, p1, p2
        real(dp), intent(out) :: i, di(NV)
        real(dp) :: hh
        hh = fd_step
        i = agg_rate(cc, cs, p1, p2)
        if (cc <= hh) then
            di(IC) = (agg_rate(cc + hh, cs, p1, p2) - i)/hh
        else
            di(IC) = (agg_rate(cc + hh, cs, p1, p2) - agg_rate(cc - hh, cs, p1, p2))/(2.0_dp*hh)
        end if
        if (cs <= hh) then
            di(ICS) = (agg_rate(cc, cs + hh, p1, p2) - i)/hh
        else
            di(ICS) = (agg_rate(cc, cs + hh, p1, p2) - agg_rate(cc, cs - hh, p1, p2))/(2.0_dp*hh)
        end if
        di(IP1) = (agg_rate(cc, cs, p1 + hh, p2) - agg_rate(cc, cs, p1 - hh, p2))/(2.0_dp*hh)
        di(IP2) = (agg_rate(cc, cs, p1, p2 + hh) - agg_rate(cc, cs, p1, p2 - hh))/(2.0_dp*hh)
    end subroutine agg_rate_fd


    subroutine agg_faces(n, dxx, aWx, aEx, bWx, bEx)
        integer, intent(in) :: n
        real(dp), intent(in) :: dxx(n)
        real(dp), intent(out) :: aWx(n), aEx(n), bWx(n), bEx(n)
        integer :: j
        aWx = 0.0_dp; bWx = 0.0_dp; aEx = 0.0_dp; bEx = 0.0_dp
        do j = 2, n
            aWx(j) = dxx(j-1)/(dxx(j-1) + dxx(j))
            bWx(j) = 2.0_dp/(dxx(j-1) + dxx(j))
        end do
        do j = 1, n - 1
            aEx(j) = dxx(j)/(dxx(j+1) + dxx(j))
            bEx(j) = 2.0_dp/(dxx(j) + dxx(j+1))
        end do
    end subroutine agg_faces

    subroutine agg_blocks(j, n, rj, dW, dE, fW, fE, aWj, aEj, bWj, bEj, Aj, Bj, Dj)
        !! Control-volume coefficients of node j -> BAND blocks.
        integer, intent(in) :: j, n
        real(dp), intent(in) :: rj(NV,NV), dW(NV,NV), dE(NV,NV), fW(NV,NV), fE(NV,NV), aWj, aEj, bWj, bEj
        real(dp), intent(out) :: Aj(NV,NV), Bj(NV,NV), Dj(NV,NV)
        Aj = 0.0_dp; Dj = 0.0_dp
        if (j == 1) then
            Bj = rj - (1.0_dp - aEj)*fE + bEj*dE
            Dj = -(aEj*fE) - bEj*dE
        else if (j == n) then
            Aj = (1.0_dp - aWj)*fW - bWj*dW
            Bj = rj + bWj*dW + aWj*fW
        else
            Aj = (1.0_dp - aWj)*fW - bWj*dW
            Bj = rj + bWj*dW + aWj*fW - (1.0_dp - aEj)*fE + bEj*dE
            Dj = -(aEj*fE) - bEj*dE
        end if
    end subroutine agg_blocks

    subroutine agg_solve(n, Ax, Bx, Dx, Gx, dcx)
        !! The archival MATINV: legacy pivot, and a block is singular only when no nonzero pivot
        !! is left (faithful mode). The original then continued with an undefined result; here the
        !! update is NaN, which ends the run as the original's NaN did.
        integer, intent(in) :: n
        real(dp), intent(in) :: Ax(NV,NV,n), Bx(NV,NV,n), Dx(NV,NV,n), Gx(NV,n)
        real(dp), intent(out) :: dcx(NV,n)
        integer :: st
        call band_solve(NV, n, Ax, Bx, Dx, Gx, dcx, st, pivot=PIVOT_LEGACY, singular=SINGULAR_EXACT)
        if (st /= BAND_OK) dcx = ieee_value(1.0_dp, ieee_quiet_nan)
    end subroutine agg_solve

    ! ================================================================================
    ! ============== corrected mode: log variables and Scharfetter-Gummel ==============
    ! ================================================================================
    ! Both models' corrected mode (docs/model.md sections 6-8). The unknowns per node are
    ! (u, phi1, phi2, s): u = ln(c/c_bulk) in column IC and the particles' log-odds
    ! s = ln(theta/(1-theta)) in column ICS, so c > 0 and 0 < theta < 1 by construction. Ion fluxes
    ! use exponential fitting; kappa_bg keeps phi2 defined where the salt is exhausted. Residuals
    ! are R(x) = 0 with blocks of dR/dx; G = -R. The agglomerate model condenses its agglomerates
    ! onto the electrode's diagonal blocks at every Newton iteration (a fully coupled step).

    pure real(dp) function bern(x)
        !! B(x) = x/(e^x - 1), with its series near 0.
        real(dp), intent(in) :: x
        if (abs(x) < 1.0e-3_dp) then
            bern = 1.0_dp - x/2.0_dp + x*x/12.0_dp
        else if (x > 700.0_dp) then
            bern = x*exp(-x)
        else
            bern = x/(exp(x) - 1.0_dp)
        end if
    end function bern

    pure real(dp) function bern_p(x)
        !! B'(x) = B(x) (1 - B(-x)) / x, with its series near 0.
        real(dp), intent(in) :: x
        if (abs(x) < 1.0e-3_dp) then
            bern_p = -0.5_dp + x/6.0_dp - x**3/180.0_dp
        else
            bern_p = bern(x)*(1.0_dp - bern(-x))/x
        end if
    end function bern_p

    pure real(dp) function softplus(x)
        real(dp), intent(in) :: x
        softplus = max(x, 0.0_dp) + log(1.0_dp + exp(-abs(x)))
    end function softplus

    pure real(dp) function sigm(x)
        real(dp), intent(in) :: x
        sigm = 0.5_dp*(1.0_dp + tanh(0.5_dp*x))
    end function sigm

    subroutine lrate(x, agg, i, di)
        !! Butler-Volmer in the log variables: i_n and d i_n / d(u, phi1, phi2, s).
        !! agg selects the agglomerate model's OCP fit and c_s,max, else the uniform model's.
        real(dp), intent(in) :: x(NV)
        logical, intent(in) :: agg
        real(dp), intent(out) :: i, di(NV)
        real(dp) :: rtf, ff, th, xx, rk, drk, term, uref, csm, dU_ds, ln_i0, eta, ea, ec, di_deta
        integer :: k, nk
        rtf = R*T/F
        ff = 1.0_dp/rtf
        th = sigm(x(ICS))
        xx = 2*th - 1
        rk = 0.0_dp; drk = 0.0_dp
        if (agg) then
            nk = 11; uref = U_REF_A; csm = csmax_a
        else
            nk = 10; uref = U_REF_U; csm = cs_max()
        end if
        do k = 0, nk
            if (agg) then
                rk = rk + AK_A(k)*(xx**(k + 1) - (2*th*k*(1 - th))/xx**(1 - k))
            else
                rk = rk + AK_U(k)*(xx**(k + 1) - (2*th*k*(1 - th))/xx**(1 - k))
            end if
            term = 2.0_dp*(2*k + 1)*xx**k
            if (k >= 2) term = term - 4.0_dp*k*(k - 1)*th*(1 - th)*xx**(k - 2)
            if (agg) then
                drk = drk + AK_A(k)*term
            else
                drk = drk + AK_U(k)*term
            end if
        end do
        dU_ds = -rtf + drk*th*(1 - th)
        ln_i0 = log(F*k_rxn*c_bulk**alpha_a*csm**(alpha_a + alpha_c)) + alpha_a*x(IC) &
                - alpha_a*softplus(x(ICS)) - alpha_c*softplus(-x(ICS))
        eta = x(IP1) - x(IP2) - (uref + rtf*(x(IC) - x(ICS)) + rk)
        ea = exp(ln_i0 + alpha_a*ff*eta)
        ec = exp(ln_i0 - alpha_c*ff*eta)
        i = ea - ec
        di_deta = alpha_a*ff*ea + alpha_c*ff*ec
        di(IC) = i*alpha_a - di_deta*rtf
        di(IP1) = di_deta
        di(IP2) = -di_deta
        di(ICS) = i*(-alpha_a*th + alpha_c*(1 - th)) - di_deta*dU_ds
    end subroutine lrate

    subroutine sg_face(xa, xb, g, gs, Fv, dFa, dFb)
        !! Face fluxes (N+, i1, i2) from state xa to state xb and their derivatives (Scharfetter-Gummel).
        !! g = eps/(tau h), gs = (1-eps) sigma/h.
        real(dp), intent(in) :: xa(NV), xb(NV), g, gs
        real(dp), intent(out) :: Fv(3), dFa(3,NV), dFb(3,NV)
        real(dp) :: ff, ca, cb, d, Bp, Bm, dBp, dBm, Np, Nm, dNp_dd, dNm_dd, kb
        ff = F/(R*T)
        ca = c_bulk*exp(xa(IC)); cb = c_bulk*exp(xb(IC))
        d = ff*(xb(IP2) - xa(IP2))
        Bp = bern(d); Bm = bern(-d); dBp = bern_p(d); dBm = bern_p(-d)
        Np = g*dplus*(Bp*ca - Bm*cb)
        Nm = g*dminus*(Bm*ca - Bp*cb)
        dNp_dd = g*dplus*(dBp*ca + dBm*cb)
        dNm_dd = g*dminus*(-dBm*ca - dBp*cb)
        kb = g*kappa_bg
        Fv(1) = Np
        Fv(2) = -gs*(xb(IP1) - xa(IP1))
        Fv(3) = F*(Np - Nm) - kb*(xb(IP2) - xa(IP2))
        dFa = 0.0_dp; dFb = 0.0_dp
        dFa(1,IC) = g*dplus*Bp*ca
        dFb(1,IC) = -g*dplus*Bm*cb
        dFa(1,IP2) = -ff*dNp_dd
        dFb(1,IP2) = ff*dNp_dd
        dFa(2,IP1) = gs
        dFb(2,IP1) = -gs
        dFa(3,IC) = F*(dFa(1,IC) - g*dminus*Bm*ca)
        dFb(3,IC) = F*(dFb(1,IC) + g*dminus*Bp*cb)
        dFa(3,IP2) = -F*ff*(dNp_dd - dNm_dd) + kb
        dFb(3,IP2) = F*ff*(dNp_dd - dNm_dd) - kb
    end subroutine sg_face

    subroutine lel_assemble(x, xold, dt, Ia, particles)
        !! Electrode residual and blocks (A, B, Dm, G = -R) in the log variables, without agglomerate
        !! sources. particles: the uniform model's particles at nodes s..nj (else the S column is fixed).
        real(dp), intent(in) :: x(NV,nj), xold(NV,nj), dt, Ia
        logical, intent(in) :: particles
        real(dp) :: Fv(3,nj-1), dFa(3,NV,nj-1), dFb(3,NV,nj-1), Rr(NV,nj), h, gf, gsf, e, cj, cjo
        real(dp) :: i, di(NV), th, tho
        integer :: j, k, fl, row
        A = 0.0_dp; B = 0.0_dp; Dm = 0.0_dp; Rr = 0.0_dp
        do k = 1, nj - 1
            h = (dx(k) + dx(k+1))/2.0_dp
            if (k < s) then
                gf = eps_sep/tau_sep/h; gsf = (1.0_dp - eps_sep)*sigma/h
            else
                gf = eps/tortuosity/h; gsf = (1.0_dp - eps)*sigma/h
            end if
            call sg_face(x(:,k), x(:,k+1), gf, gsf, Fv(:,k), dFa(:,:,k), dFb(:,:,k))
        end do
        ! foil face: N+ = I/F, zero electronic current, phi2 = 0 (the gauge)
        Rr(IC,1) = Fv(1,1) - Ia/F
        B(IC,:,1) = dFa(1,:,1); Dm(IC,:,1) = dFb(1,:,1)
        Rr(IP1,1) = x(IP1,2) - x(IP1,1)
        B(IP1,IP1,1) = -1.0_dp; Dm(IP1,IP1,1) = 1.0_dp
        Rr(IP2,1) = x(IP2,1)
        B(IP2,IP2,1) = 1.0_dp
        ! separator, interface and cathode
        do j = 2, nj - 1
            do fl = 1, 3
                row = fl           ! flux fl -> row IC, IP1, IP2
                Rr(row,j) = Fv(fl,j) - Fv(fl,j-1)
                B(row,:,j) = B(row,:,j) + dFa(fl,:,j) - dFb(fl,:,j-1)
                Dm(row,:,j) = Dm(row,:,j) + dFb(fl,:,j)
                A(row,:,j) = A(row,:,j) - dFa(fl,:,j-1)
            end do
            e = eps
            if (j < s) e = eps_sep
            cj = c_bulk*exp(x(IC,j)); cjo = c_bulk*exp(xold(IC,j))
            Rr(IC,j) = Rr(IC,j) + e*dx(j)*(cj - cjo)/dt
            B(IC,IC,j) = B(IC,IC,j) + e*dx(j)*cj/dt
        end do
        ! collector: no salt flux, no ionic current, electronic current I
        Rr(IC,nj) = -Fv(1,nj-1); A(IC,:,nj) = -dFa(1,:,nj-1); B(IC,:,nj) = -dFb(1,:,nj-1)
        Rr(IP1,nj) = Fv(2,nj-1) - Ia; A(IP1,:,nj) = dFa(2,:,nj-1); B(IP1,:,nj) = dFb(2,:,nj-1)
        Rr(IP2,nj) = Fv(3,nj-1); A(IP2,:,nj) = dFa(3,:,nj-1); B(IP2,:,nj) = dFb(3,:,nj-1)
        ! the S column: fixed, or the uniform model's particles
        do j = 1, nj
            Rr(ICS,j) = x(ICS,j) - xold(ICS,j)
            B(ICS,:,j) = 0.0_dp
            B(ICS,ICS,j) = 1.0_dp
        end do
        if (particles) then
            do j = s, nj
                call lrate(x(:,j), .false., i, di)
                Rr(IC,j) = Rr(IC,j) - spec_a*i*dx(j)/F
                B(IC,:,j) = B(IC,:,j) - spec_a*di*dx(j)/F
                Rr(IP1,j) = Rr(IP1,j) + spec_a*i*dx(j)
                B(IP1,:,j) = B(IP1,:,j) + spec_a*di*dx(j)
                Rr(IP2,j) = Rr(IP2,j) - spec_a*i*dx(j)
                B(IP2,:,j) = B(IP2,:,j) - spec_a*di*dx(j)
                th = sigm(x(ICS,j)); tho = sigm(xold(ICS,j))
                Rr(ICS,j) = vf_AM*cs_max()*(th - tho)/dt + spec_a*i/F
                B(ICS,:,j) = spec_a*di/F
                B(ICS,ICS,j) = B(ICS,ICS,j) + vf_AM*cs_max()*th*(1 - th)/dt
            end do
        end if
        G = -Rr
    end subroutine lel_assemble

    real(dp) function lbounded(d)
        !! Newton step length <= 1 limiting |du| <= 1, |dphi| <= 0.1 V and |ds| <= 2 per iteration.
        real(dp), intent(in) :: d(:,:)
        real(dp) :: caps(NV), mx
        integer :: k
        caps = [1.0_dp, 0.1_dp, 0.1_dp, 2.0_dp]
        lbounded = 1.0_dp
        do k = 1, NV
            mx = maxval(abs(d(k,:)))
            if (mx > caps(k)) lbounded = min(lbounded, caps(k)/mx)
        end do
    end function lbounded

    real(dp) function phys_update(x, d, frozen_s)
        !! The Newton update in the physical variables: max of e^u |du|, |dphi| and theta(1-theta) |ds|.
        real(dp), intent(in) :: x(:,:), d(:,:)
        logical, intent(in) :: frozen_s
        integer :: j
        real(dp) :: th
        phys_update = 0.0_dp
        do j = 1, size(x, 2)
            phys_update = max(phys_update, exp(x(IC,j))*abs(d(IC,j)), abs(d(IP1,j)), abs(d(IP2,j)))
            if (.not. frozen_s) then
                th = sigm(x(ICS,j))
                phys_update = max(phys_update, th*(1 - th)*abs(d(ICS,j)))
            end if
        end do
    end function phys_update

    subroutine lu_newton(h, ok)
        !! Uniform model, corrected mode: one backward-Euler step of length h (log variables).
        real(dp), intent(in) :: h
        logical, intent(out) :: ok
        real(dp) :: c_old(NV,nj), lam, upd, prev, raw
        integer :: k, st
        c_old = c
        ok = .false.
        prev = huge(1.0_dp)
        do k = 1, newton_max_iter
            call lel_assemble(c, c_old, h, i_app, .true.)
            call equilibrate()
            call band_solve(NV, nj, A, B, Dm, G, dc, st)
            if (st /= BAND_OK) exit
            lam = lbounded(dc)
            raw = maxval(abs(dc))
            c = c + lam*dc
            upd = phys_update(c, dc, .false.)
            if (.not. ieee_is_finite(raw) .or. raw > 1.0e3_dp) exit
            if (lam == 1.0_dp .and. converged(upd, prev)) then
                ok = .true.
                return
            end if
            prev = upd
        end do
        c = c_old
    end subroutine lu_newton

    ! ------------------------------------------------------------------ agglomerate model
    subroutine aggc_setup()
        !! After setup(): the corrected agglomerate model's grid, transport and state.
        real(dp) :: hr, tau_a, r_W, r_E, rf
        integer :: j
        mass_area = percent_active*mass_loading          ! = L_cath*v_AM*rho
        i_1C = Q_th*mass_area
        v_AM = percent_active*mass_loading/(rho*L_cath)
        v_agg = v_AM/(1.0_dp - eps_agg)
        s_agg = 3.0_dp*v_agg/R_agg
        a_x = 3.0_dp*(1.0_dp - eps_agg)/R_xtal
        if (eps + v_agg > 1.0_dp + 1.0e-12_dp) then
            write(error_unit,'(A)') 'porosity + agglomerate volume fraction exceed 1'
            error stop 2
        end if
        x_max_a = M*Q_th*3600.0_dp/F
        csmax_a = mol_vol*x_max_a
        tau_a = tau_agg
        if (tau_a < 0) tau_a = eps_agg**(-0.5_dp)
        sig_a = sigma_agg
        if (sig_a < 0) sig_a = sigma
        ! radial grid (zero-volume center and surface nodes)
        na = nja
        nl = nj - s - 1
        allocate(dxa(na), xa(na), dV(na), a_area(na-1), a_g(na-1), a_gs(na-1))
        hr = R_agg/real(na - 2, dp)
        xa = 0.0_dp
        do j = 2, na - 1
            xa(j) = hr*real(j - 1, dp) - hr/2.0_dp
        end do
        xa(na) = R_agg
        dxa = 0.0_dp; dxa(2:na-1) = hr
        do j = 1, na
            r_W = xa(j) - dxa(j)/2.0_dp; r_E = xa(j) + dxa(j)/2.0_dp
            if (j == 1) then
                r_W = 0.0_dp; r_E = 0.0_dp
            end if
            if (j == na) then
                r_W = R_agg; r_E = R_agg
            end if
            dV(j) = 4.0_dp*acos(-1.0_dp)/3.0_dp*(r_E**3 - r_W**3)
        end do
        do j = 1, na - 1                          ! face j joins nodes j and j+1
            rf = xa(j) + dxa(j)/2.0_dp
            if (j == 1) rf = 0.0_dp
            a_area(j) = 4.0_dp*acos(-1.0_dp)*rf**2
            a_g(j) = eps_agg/tau_a/((dxa(j) + dxa(j+1))/2.0_dp)
            a_gs(j) = (1.0_dp - eps_agg)*sig_a/((dxa(j) + dxa(j+1))/2.0_dp)
        end do
        allocate(cag(NV,na,nl), cag_start(NV,na,nl), Aq(NV,NV,na*nl), Bq(NV,NV,na*nl), Dq(NV,NV,na*nl), Gq(NV,na*nl))
        cag(IC,:,:) = log(c0_init/c_bulk); cag(IP1,:,:) = phi1_init; cag(IP2,:,:) = 0.0_dp
        cag(ICS,:,:) = log((cs_init/csmax_a)/(1.0_dp - cs_init/csmax_a))
        dt = dt_s
    end subroutine aggc_setup

    subroutine aggc_assemble(h, cold)
        !! Stacked blocks Aq, Bq, Dq and G = -R (Gq) of every agglomerate; surface values from c.
        real(dp), intent(in) :: h, cold(NV,na,nl)
        real(dp) :: Fv(3,na-1), dFa(3,NV,na-1), dFb(3,NV,na-1), Rr(NV), i, di(NV), V, cj, cjo, th, tho
        integer :: l, j, k, q, col, fl
        do l = 1, nl
            do k = 1, na - 1
                call sg_face(cag(:,k,l), cag(:,k+1,l), a_g(k), a_gs(k), Fv(:,k), dFa(:,:,k), dFb(:,:,k))
            end do
            do j = 1, na
                q = (l - 1)*na + j
                Aq(:,:,q) = 0.0_dp; Bq(:,:,q) = 0.0_dp; Dq(:,:,q) = 0.0_dp; Rr = 0.0_dp
                call lrate(cag(:,j,l), .true., i, di)
                if (j == 1) then
                    ! center: zero gradient
                    do col = IC, IP2
                        Rr(col) = cag(col,2,l) - cag(col,1,l)
                        Bq(col,col,q) = -1.0_dp; Dq(col,col,q) = 1.0_dp
                    end do
                else if (j == na) then
                    ! surface: the electrode's u, phi1, phi2
                    do col = IC, IP2
                        Rr(col) = cag(col,na,l) - c(col, s + l)
                        Bq(col,col,q) = 1.0_dp
                    end do
                else
                    V = dV(j)
                    do fl = 1, 3
                        Rr(fl) = a_area(j)*Fv(fl,j) - a_area(j-1)*Fv(fl,j-1)
                        Bq(fl,:,q) = a_area(j)*dFa(fl,:,j) - a_area(j-1)*dFb(fl,:,j-1)
                        Dq(fl,:,q) = a_area(j)*dFb(fl,:,j)
                        Aq(fl,:,q) = -a_area(j-1)*dFa(fl,:,j-1)
                    end do
                    cj = c_bulk*exp(cag(IC,j,l)); cjo = c_bulk*exp(cold(IC,j,l))
                    Rr(IC) = Rr(IC) + eps_agg*V*(cj - cjo)/h - a_x*i*V/F
                    Bq(IC,IC,q) = Bq(IC,IC,q) + eps_agg*V*cj/h
                    Bq(IC,:,q) = Bq(IC,:,q) - a_x*di*V/F
                    Rr(IP1) = Rr(IP1) + a_x*i*V
                    Bq(IP1,:,q) = Bq(IP1,:,q) + a_x*di*V
                    Rr(IP2) = Rr(IP2) - a_x*i*V
                    Bq(IP2,:,q) = Bq(IP2,:,q) - a_x*di*V
                end if
                ! crystals, every node: (1 - eps_agg) dcs/dt = -a i/F
                th = sigm(cag(ICS,j,l)); tho = sigm(cold(ICS,j,l))
                Rr(ICS) = (1.0_dp - eps_agg)*csmax_a*(th - tho)/h + a_x*i/F
                Bq(ICS,:,q) = a_x*di/F
                Bq(ICS,ICS,q) = Bq(ICS,ICS,q) + (1.0_dp - eps_agg)*csmax_a*th*(1 - th)/h
                Gq(:,q) = -Rr
            end do
        end do
    end subroutine aggc_assemble

    subroutine aggc_surface(l, qv, Jq)
        !! Flux into agglomerate l through r = R: q = (N+, i1, i2), and dq/dx at nodes na-1 and na.
        integer, intent(in) :: l
        real(dp), intent(out) :: qv(3), Jq(3, 2*NV)
        real(dp) :: Fv(3), dFa(3,NV), dFb(3,NV)
        call sg_face(cag(:,na-1,l), cag(:,na,l), a_g(na-1), a_gs(na-1), Fv, dFa, dFb)
        qv = -Fv
        Jq(:,1:NV) = -dFa
        Jq(:,NV+1:2*NV) = -dFb
    end subroutine aggc_surface

    subroutine aggc_newton(h, ok)
        !! One backward-Euler step of length h for both scales, by the condensed Newton iteration.
        real(dp), intent(in) :: h
        logical, intent(out) :: ok
        real(dp) :: c_old(NV,nj), cag_old(NV,na,nl), lam, sc, w, upd, prev, raw
        real(dp) :: x0(NV,na*nl), Z(NV,na*nl,3), E(NV,na*nl), qv(3), Jq(3,2*NV), dca(NV,na,nl), Zl(2*NV,3), x0l(2*NV)
        real(dp) :: rsurf(3,nl), dca_flat(NV,na*nl)
        type(band_factorization) :: fac
        integer :: it, st, k, q, r, l, j, col, nq
        nq = na*nl
        c_old = c; cag_old = cag
        ok = .false.
        prev = huge(1.0_dp)
        do it = 1, newton_max_iter
            ! ---- agglomerates: factor once, solve for the update and the three surface responses ----
            call aggc_assemble(h, cag_old)
            do q = 1, nq
                do r = 1, NV
                    sc = maxval(abs(Bq(r,:,q)))
                    if (sc == 0.0_dp) sc = 1.0_dp
                    Aq(r,:,q) = Aq(r,:,q)/sc; Bq(r,:,q) = Bq(r,:,q)/sc; Dq(r,:,q) = Dq(r,:,q)/sc
                    Gq(r,q) = Gq(r,q)/sc
                    if (mod(q, na) == 0 .and. r <= IP2) rsurf(r, q/na) = 1.0_dp/sc
                end do
            end do
            call band_factor(NV, nq, Aq, Bq, Dq, fac)
            if (fac%status /= BAND_OK) exit
            call band_factor_solve(fac, Gq, x0, st)
            if (st /= BAND_OK) exit
            do col = 1, 3
                E = 0.0_dp
                do l = 1, nl
                    E(col, l*na) = rsurf(col, l)
                end do
                call band_factor_solve(fac, E, Z(:,:,col), st)
                if (st /= BAND_OK) exit
            end do
            if (st /= BAND_OK) exit
            ! ---- electrode with the condensed agglomerate sources: R_e + w q ----
            call lel_assemble(c, c_old, h, i_app, .false.)
            do l = 1, nl
                j = s + l
                w = s_agg*dx(j)
                call aggc_surface(l, qv, Jq)
                x0l(1:NV) = x0(:, l*na - 1); x0l(NV+1:2*NV) = x0(:, l*na)
                do col = 1, 3
                    Zl(1:NV, col) = Z(:, l*na - 1, col); Zl(NV+1:2*NV, col) = Z(:, l*na, col)
                end do
                G(1:3, j) = G(1:3, j) - w*qv - w*matmul(Jq, x0l)
                B(1:3, 1:3, j) = B(1:3, 1:3, j) + w*matmul(Jq, Zl)
            end do
            call equilibrate()
            call band_solve(NV, nj, A, B, Dm, G, dc, st)
            if (st /= BAND_OK) exit
            do l = 1, nl
                j = s + l
                do k = 1, na
                    q = (l - 1)*na + k
                    dca(:,k,l) = x0(:,q) + Z(:,q,1)*dc(IC,j) + Z(:,q,2)*dc(IP1,j) + Z(:,q,3)*dc(IP2,j)
                end do
            end do
            dca_flat = reshape(dca, [NV, nq])
            lam = min(lbounded(dc), lbounded(dca_flat))
            raw = max(maxval(abs(dc)), maxval(abs(dca)))
            c = c + lam*dc
            cag = cag + lam*dca
            upd = max(phys_update(c, dc, .true.), phys_update(reshape(cag, [NV, nq]), dca_flat, .false.))
            if (.not. ieee_is_finite(raw) .or. raw > 1.0e3_dp) exit
            if (lam == 1.0_dp .and. converged(upd, prev)) then
                ok = .true.
                return
            end if
            prev = upd
        end do
        c = c_old
        cag = cag_old
    end subroutine aggc_newton

end program nmc
