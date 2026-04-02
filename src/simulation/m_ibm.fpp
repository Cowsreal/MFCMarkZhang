!>
!! @file
!! @brief Contains module m_ibm

#:include 'macros.fpp'

!> @brief Ghost-node immersed boundary method: locates ghost/image points, computes interpolation coefficients, and corrects the flow state
module m_ibm

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_variables_conversion !< State variables type conversion procedures

    use m_helper

    use m_helper_basic         !< Functions to compare floating point numbers

    use m_constants

    use m_compute_levelset

    use m_ib_patches

    use m_viscous

    use m_model

    use m_boundary_common

    implicit none

    private :: s_compute_image_points, &
               s_compute_interpolation_coeffs, &
               s_interpolate_image_point, &
               s_find_ghost_points, &
               s_find_num_ghost_points
    ; public :: s_initialize_ibm_module, &
 s_ibm_setup, &
 s_ibm_correct_state, &
 s_finalize_ibm_module, &
 s_interpolate_sigma_igr

    type(integer_field), public :: ib_markers
    $:GPU_DECLARE(create='[ib_markers]')

    type(ghost_point), dimension(:), allocatable :: ghost_points
    type(ghost_point), dimension(:), allocatable :: inner_points
    type(ghost_point), dimension(:), allocatable :: domain_points
    $:GPU_DECLARE(create='[ghost_points,inner_points,domain_points]')

    integer :: num_gps !< Number of ghost points
    integer :: num_inner_gps !< Number of ghost points
    integer :: gp_total
    integer, allocatable :: gp_ib_patch_id(:)
    integer, allocatable :: gp_loc(:, :)
    integer, allocatable :: gp_M(:)
    integer, allocatable :: gp_offset(:)
    integer, allocatable :: gp_stencil(:, :)
    real(wp), allocatable :: gp_A_neum(:)
    real(wp), allocatable :: gp_A_dirich(:)
#if defined(MFC_OpenACC)
    $:GPU_DECLARE(create='[num_gps,num_inner_gps]')
    $:GPU_DECLARE(create='[gp_total,gp_ib_patch_id,gp_loc,gp_M,gp_offset,gp_stencil,gp_A_neum,gp_A_dirich]')
#elif defined(MFC_OpenMP)
    $:GPU_DECLARE(create='[num_gps,num_inner_gps]')
#endif
    logical :: moving_immersed_boundary_flag

contains

    !>  Allocates memory for the variables in the IBM module
    impure subroutine s_initialize_ibm_module()

        if (p > 0) then
            @:ALLOCATE(ib_markers%sf(-buff_size:m+buff_size, &
                -buff_size:n+buff_size, -buff_size:p+buff_size))
        else
            @:ALLOCATE(ib_markers%sf(-buff_size:m+buff_size, &
                -buff_size:n+buff_size, 0:0))
        end if

        @:ALLOCATE(models(num_ibs))

        @:ACC_SETUP_SFs(ib_markers)

        $:GPU_ENTER_DATA(copyin='[num_gps,num_inner_gps]')

    end subroutine s_initialize_ibm_module

    !> Initializes the values of various IBM variables, such as ghost points and
    !! image points.
    impure subroutine s_ibm_setup()

        integer :: i, j, k
        integer :: max_num_gps, max_num_inner_gps

        call nvtxStartRange("SETUP-IBM-MODULE")

        ! do all set up for moving immersed boundaries
        moving_immersed_boundary_flag = .false.
        do i = 1, num_ibs
            if (patch_ib(i)%moving_ibm /= 0) then
                call s_compute_moment_of_inertia(i, patch_ib(i)%angular_vel)
                moving_immersed_boundary_flag = .true.
            end if
            call s_update_ib_rotation_matrix(i)
        end do
        $:GPU_UPDATE(device='[patch_ib(1:num_ibs)]')

        ! GPU routines require updated cell centers
        $:GPU_UPDATE(device='[num_ibs, x_cc, y_cc, dx, dy, x_domain, y_domain]')
        if (p /= 0) then
            $:GPU_UPDATE(device='[z_cc, dz, z_domain]')
        end if

        ! allocate STL models
        call s_instantiate_STL_models()

        ! recompute the new ib_patch locations and broadcast them.
        ib_markers%sf = 0._wp
        $:GPU_UPDATE(device='[ib_markers%sf]')
        call s_apply_ib_patches(ib_markers)
        $:GPU_UPDATE(host='[ib_markers%sf]')
        do i = 1, num_ibs
            if (patch_ib(i)%moving_ibm /= 0) call s_compute_centroid_offset(i) ! offsets are computed after IB markers are generated
            $:GPU_UPDATE(device='[patch_ib(i)]')
        end do

        ! find the number of ghost points and set them to be the maximum total across ranks
        call s_find_num_ghost_points(num_gps, num_inner_gps)
        call s_mpi_allreduce_integer_sum(num_gps, max_num_gps)
        call s_mpi_allreduce_integer_sum(num_inner_gps, max_num_inner_gps)
        max_num_gps = min(max_num_gps, (m + 1)*(n + 1)*(p + 1)/2)
        max_num_inner_gps = min(max_num_inner_gps, (m + 1)*(n + 1)*(p + 1)/2)

        ! set the size of the ghost point arrays to be the amount of points total, plus a factor of 2 buffer
        $:GPU_UPDATE(device='[num_gps, num_inner_gps]')
        if (moving_immersed_boundary_flag) then
            @:ALLOCATE(ghost_points(1:int((max_num_gps + max_num_inner_gps) * 2.0)))
            @:ALLOCATE(inner_points(1:int((max_num_gps + max_num_inner_gps) * 2.0)))
        else
            @:ALLOCATE(ghost_points(1:int(max_num_gps + max_num_inner_gps)))
            @:ALLOCATE(inner_points(1:int(max_num_gps + max_num_inner_gps)))
        end if

        $:GPU_ENTER_DATA(copyin='[ghost_points,inner_points]')
        call s_find_ghost_points(ghost_points, inner_points)
        call s_apply_levelset(ghost_points, num_gps)

        $:GPU_UPDATE(host='[ghost_points]')
        if(num_gps > 0) then

            if (patch_ib(1)%hybrid) then
                do i = 1, num_gps
                    if (ghost_points(i)%first_layer) then
                        call s_radial_search(ghost_points(i))
                    else
                        call s_compute_image_points(ghost_points(i))
                        call s_compute_interpolation_coeffs(ghost_points(i))
                    end if
                end do
            else
                do i = 1, num_gps
                    if(patch_ib(1)%high_order /= 0) then
                        call s_radial_search(ghost_points(i))
                    else
                        call s_compute_image_points(ghost_points(i))
                        call s_compute_interpolation_coeffs(ghost_points(i))
                    end if
                end do
            end if

            call s_build_gp_csr()
        end if

        call nvtxEndRange

    end subroutine s_ibm_setup

    subroutine s_build_gp_csr()
        integer :: i, q, offs, total,j , k
        logical :: bilinear

        @:ALLOCATE(gp_loc(1:3, num_gps))
        @:ALLOCATE(gp_ib_patch_id(num_gps))
        @:ALLOCATE(gp_M(1:num_gps))
        @:ALLOCATE(gp_offset(1:num_gps + 1))

        ! Define offsets arrays
        gp_offset(1) = 0
        do i = 1, num_gps
            bilinear = (patch_ib(1)%high_order == 0) .or. &
                        (patch_ib(1)%hybrid .and. .not. ghost_points(i)%first_layer)
            if(bilinear) then
                gp_M(i) = 4
            else
                gp_M(i) = ghost_points(i)%M
            end if
            gp_offset(i + 1) = gp_offset(i) + gp_M(i)
        end do

        gp_total = gp_offset(num_gps + 1)

        ! Allocate stencil related arrays
        @:ALLOCATE(gp_stencil(1:3, 1:gp_total))
        @:ALLOCATE(gp_A_neum(1:gp_total))
        @:ALLOCATE(gp_A_dirich(1:gp_total))

        ! Now copy all the data
        do i = 1, num_gps
            gp_loc(1, i) = ghost_points(i)%loc(1)
            gp_loc(2, i) = ghost_points(i)%loc(2)
            gp_loc(3, i) = ghost_points(i)%loc(3)
            gp_ib_patch_id(i) = ghost_points(i)%ib_patch_id

            offs = gp_offset(i)
            bilinear = (patch_ib(1)%high_order == 0) .or. &
                        (patch_ib(1)%hybrid .and. .not. ghost_points(i)%first_layer)

            ! Copy each stencil point
            if(bilinear) then
                do q = 1, gp_M(i)
                    j = (q - 1) / 2
                    k = mod(q - 1, 2)
                    gp_stencil(1, offs + q) = ghost_points(i)%ip_grid(1) + j
                    gp_stencil(2, offs + q) = ghost_points(i)%ip_grid(2) + k
                    !gp_stencil(3, offs + q) = ghost_points(i)%ip_grid(3)
                    gp_stencil(3, offs + q) = 0
                    gp_A_neum(offs + q) = ghost_points(i)%interp_coeffs(j + 1, k + 1, 1)
                    gp_A_dirich(offs + q) = ghost_points(i)%interp_coeffs(j + 1, k + 1, 1)
                end do
            else
                do q = 1, gp_M(i)
                    gp_stencil(1, offs + q) = ghost_points(i)%stencil(q, 1)
                    gp_stencil(2, offs + q) = ghost_points(i)%stencil(q, 2)
                    gp_stencil(3, offs + q) = ghost_points(i)%stencil(q, 3)
                    gp_A_neum(offs + q) = ghost_points(i)%A_temp(1, q)
                    gp_A_dirich(offs + q) = ghost_points(i)%A_temp(2, q)
                end do
            end if
        end do

        $:GPU_UPDATE(device='[gp_total,gp_loc,gp_ib_patch_id,gp_M,gp_offset, &
            gp_stencil,gp_A_neum,gp_A_dirich]')

    end subroutine s_build_gp_csr

    subroutine s_smooth_ib_boundaries(bc_type, q_cons_vf)
        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf
        integer :: i, j, k, l, smooth_radius,d,num_smooth_points, local_idx
        real(wp) :: r,alpha,steep,shift
        real(wp) :: r_bound, search_radius, dynPresOld, dynPresNew
        type(ghost_point) :: gp
        type(integer_field), dimension(1:num_dims, 1:2), intent(in) :: bc_type
        integer :: count

        search_radius = patch_ib(1)%radius * (1._wp + patch_ib(1)%smooth + 0.5_wp)
        r_bound = patch_ib(1)%radius * patch_ib(1)%smooth
        count = 0
        $:GPU_PARALLEL_LOOP(private='[j,k,l,gp]', copyin='[r_bound,search_radius]', copy='[count]', reduction='[[count]]',reductionOp='[+]', collapse=3)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    if (abs(x_cc(j) - patch_ib(1)%x_centroid) < search_radius &
                        .and. abs(y_cc(k) - patch_ib(1)%y_centroid) < search_radius &
                        .and. (p == 0 .or. abs(z_cc(l) - patch_ib(1)%z_centroid) < search_radius)) then
                        gp%ib_patch_id = 1
                        gp%loc = [j, k, l]
                        gp%x_periodicity = 0
                        gp%y_periodicity = 0
                        gp%z_periodicity = 0
                        call s_compute_levelset(gp)
                        if (gp%levelset < r_bound .and. gp%levelset >= 0._wp) then
                            count = count + 1
                        end if
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        @:ALLOCATE(domain_points(1:count))
        $:GPU_ENTER_DATA(copyin='[domain_points]')

        num_smooth_points = count
        count = 0
        $:GPU_PARALLEL_LOOP(private='[j,k,l,local_idx,gp]', copyin='[search_radius,r_bound,count]', collapse=3)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    if (abs(x_cc(j) - patch_ib(1)%x_centroid) < search_radius &
                        .and. abs(y_cc(k) - patch_ib(1)%y_centroid) < search_radius &
                        .and. (p == 0 .or. abs(z_cc(l) - patch_ib(1)%z_centroid) < search_radius)) then
                        gp%ib_patch_id = 1
                        gp%loc = [j, k, l]
                        gp%x_periodicity = 0
                        gp%y_periodicity = 0
                        gp%z_periodicity = 0
                        call s_compute_levelset(gp)
                        if (gp%levelset < r_bound .and. gp%levelset >= 0._wp) then
                            $:GPU_ATOMIC(atomic='capture')
                            count = count + 1
                            local_idx = count
                            $:END_GPU_ATOMIC_CAPTURE()
                            domain_points(local_idx)%loc = [j, k, l]
                            domain_points(local_idx)%levelset = gp%levelset
                        end if
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        steep = patch_ib(1)%steep
        shift = patch_ib(1)%shift
        call s_populate_variables_buffers(bc_type, q_cons_vf)
        $:GPU_PARALLEL_LOOP(private='[dynPresOld,dynPresNew,alpha,i,j,k,l,d,gp]', &
            copyin='[r_bound,num_smooth_points]', firstprivate='[steep,shift]')
        do i = 1, num_smooth_points
            gp = domain_points(i)
            j = domain_points(i)%loc(1)
            k = domain_points(i)%loc(2)
            l = domain_points(i)%loc(3)
            alpha = domain_points(i)%levelset / r_bound
            alpha = 0.5_wp * (1.0_wp + tanh(steep * (alpha - shift)) / tanh(steep / 2._wp))
            dynPresOld = 0._wp
            dynPresNew = 0._wp
            do d = 1, num_dims
                dynPresOld = dynPresOld + q_cons_vf(momxb + d - 1)%sf(j, k, l) ** 2
                dynPresNew = dynPresNew + (alpha * q_cons_vf(momxb + d - 1)%sf(j, k, l)) ** 2
                q_cons_vf(momxb + d - 1)%sf(j, k, l) = alpha * q_cons_vf(momxb + d - 1)%sf(j, k, l)
            end do
            dynPresOld = dynPresOld / (2._wp * q_cons_vf(contxb)%sf(j, k, l))
            dynPresNew = dynPresNew / (2._wp * q_cons_vf(contxb)%sf(j, k, l))
            q_cons_vf(E_idx)%sf(j, k, l) = q_cons_vf(E_idx)%sf(j, k, l) - dynPresOld + dynPresNew
            !q_cons_vf(E_idx)%sf(j, k, l) = domain_points(i)%levelset
        end do
        $:END_GPU_PARALLEL_LOOP()
        ! do i = 1, num_gps
        !     gp = ghost_points(i)
        !     j = gp%loc(1)
        !     k = gp%loc(2)
        !     l = gp%loc(3)
        !     q_cons_vf(E_idx)%sf(j, k, l) = gp%levelset
        ! end do

    @:DEALLOCATE(domain_points)
    end subroutine s_smooth_ib_boundaries

    !>  Interpolates sigma from the m_igr module at all ghost points
        !!  @param jac: Sigma, Entropic pressure
    subroutine s_interpolate_sigma_igr(jac)

        real(wp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:), intent(inout) :: jac
        integer :: j, k, l, r, s, t, i, q, offs, stencil_idx
        integer :: j1, j2, k1, k2, l1, l2
        real(wp) :: coeff, jac_IP, neum
        type(ghost_point) :: gp

        ! At all ghost points, use its image point to interpolate sigma
        if (num_gps > 0) then
            $:GPU_PARALLEL_LOOP(private='[q,offs,stencil_idx,neum,i, j, k, l, &
                j1, j2, k1, k2, l1, l2, r, s, t, coeff, jac_IP]')
            do i = 1, num_gps
                !jac_IP = 0._wp
                !gp = ghost_points(i)
                r = gp_loc(1, i)
                s = gp_loc(2, i)
                t = gp_loc(3, i)


                !j1 = gp%ip_grid(1); j2 = j1 + 1
                !k1 = gp%ip_grid(2); k2 = k1 + 1
                !l1 = gp%ip_grid(3); l2 = l1 + 1

                !if (p == 0) then
                !    l1 = 0
                !    l2 = 0
                !end if

                if(patch_ib(1)%high_order /= 0) then
                    jac(r, s, t) = 0._wp
                    ! $:GPU_LOOP(parallelism='[seq]')
                    ! do q = 2, ghost_points(i)%M
                    !     j = ghost_points(i)%stencil(q, 1)
                    !     k = ghost_points(i)%stencil(q, 2)
                    !     l = ghost_points(i)%stencil(q, 3)
                    !     jac(r, s, t) = jac(r, s, t) &
                    !         + ghost_points(i)%A_temp(1, q) * jac(j, k, l)
                    ! end do
                    offs = gp_offset(i)
                    $:GPU_LOOP(parallelism='[seq]')
                    do q = 1, gp_M(i)
                        stencil_idx = offs + q
                        j = gp_stencil(1, stencil_idx)
                        k = gp_stencil(2, stencil_idx)
                        l = gp_stencil(3, stencil_idx)
                        jac(r, s, t) = jac(r, s, t) &
                            + gp_A_neum(stencil_idx) * jac(j, k, l)
                    end do
                else
                    ! Interpolate sigma using coeffs at the points around the corresponding image point
                    ! $:GPU_LOOP(parallelism='[seq]')
                    ! do l = l1, l2
                    !     $:GPU_LOOP(parallelism='[seq]')
                    !     do k = k1, k2
                    !         $:GPU_LOOP(parallelism='[seq]')
                    !         do j = j1, j2
                    !             coeff = gp%interp_coeffs(j - j1 + 1, k - k1 + 1, l - l1 + 1)
                    !             jac_IP = jac_IP + coeff*jac(j, k, l)
                    !         end do
                    !     end do
                    ! end do
                    ! jac(r, s, t) = jac_IP
                end if
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

    end subroutine s_interpolate_sigma_igr

    !>  Subroutine that updates the conservative variables at the ghost points
        !!  @param pb_in Internal bubble pressure
        !!  @param mv_in Mass of vapor in bubble
    subroutine s_ibm_correct_state(q_cons_vf, q_prim_vf, pb_in, mv_in)

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: q_cons_vf !< Primitive Variables

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: q_prim_vf !< Primitive Variables

        real(stp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:, 1:), optional, intent(INOUT) :: pb_in, mv_in

        integer :: i, j, k, l, q, r!< Iterator variables
        integer :: u, v
        integer :: patch_id !< Patch ID of ghost point
        real(wp) :: rho, gamma, pi_inf, dyn_pres !< Mixture variables
        real(wp), dimension(2) :: Re_K
        real(wp) :: G_K
        real(wp) :: qv_K
        integer :: stencil_idx, offs

        real(wp) :: pres_IP, E_IP
        real(wp), dimension(3) :: vel_IP
        real(wp), dimension(3) :: vel_wall, rel_norm_vel
        real(wp) :: c_IP
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3) :: Gs
            real(wp), dimension(3) :: alpha_rho_IP, alpha_IP
            real(wp), dimension(3) :: r_IP, v_IP, pb_IP, mv_IP
            real(wp), dimension(18) :: nmom_IP
            real(wp), dimension(12) :: presb_IP, massv_IP
        #:else
            real(wp), dimension(num_fluids) :: Gs
            real(wp), dimension(num_fluids) :: alpha_rho_IP, alpha_IP
            real(wp), dimension(nb) :: r_IP, v_IP, pb_IP, mv_IP
            real(wp), dimension(nb*nmom) :: nmom_IP
            real(wp), dimension(nb*nnode) :: presb_IP, massv_IP
        #:endif
        !! Primitive variables at the image point associated with a ghost point,
        !! interpolated from surrounding fluid cells.

        real(wp), dimension(3) :: norm !< Normal vector from GP to IP
        real(wp), dimension(3) :: physical_loc !< Physical loc of GP
        real(wp), dimension(3) :: vel_g !< Velocity of GP
        real(wp), dimension(3) :: radial_vector !< vector from centroid to ghost point
        real(wp), dimension(3) :: rotation_velocity !< speed of the ghost point due to rotation

        real(wp) :: nbub
        type(ghost_point) :: gp, innerp
        real(wp) :: dirich,neum,rho_inv


        if(.not. igr .or. dummy) then
            ! set the Moving IBM interior Pressure Values
            $:GPU_PARALLEL_LOOP(private='[i,j,k,patch_id,rho]', copyin='[E_idx,momxb]', collapse=3)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        patch_id = ib_markers%sf(j, k, l)
                        if (patch_id /= 0) then
                            q_prim_vf(E_idx)%sf(j, k, l) = 1._wp
                            if (patch_ib(patch_id)%moving_ibm > 0) then
                                rho = 0._wp
                                do i = 1, num_fluids
                                    rho = rho + q_prim_vf(contxb + i - 1)%sf(j, k, l)
                                end do

                                ! Sets the momentum
                                do i = 1, num_dims
                                    q_cons_vf(momxb + i - 1)%sf(j, k, l) = patch_ib(patch_id)%vel(i)*rho
                                    q_prim_vf(momxb + i - 1)%sf(j, k, l) = patch_ib(patch_id)%vel(i)
                                end do
                            end if
                        end if
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        if (num_gps > 0) then
            $:GPU_PARALLEL_LOOP(private='[rho_inv,neum,dirich,u,v,i,physical_loc,dyn_pres,alpha_rho_IP,alpha_IP,pres_IP,vel_IP,vel_g, &
                vel_wall,rel_norm_vel,r_IP,E_IP,v_IP,pb_IP,mv_IP,nmom_IP,presb_IP,massv_IP,rho,gamma,pi_inf,Re_K,G_K,Gs, &
                innerp,norm, radial_vector, rotation_velocity, j,k,l,q,qv_K,c_IP,nbub,patch_id,offs,stencil_idx]')
            do i = 1, num_gps
                j = gp_loc(1, i)
                k = gp_loc(2, i)
                l = gp_loc(3, i)
                patch_id = gp_ib_patch_id(i)
                vel_IP(1) = 0._wp
                vel_IP(2) = 0._wp
                rho = 0._wp
                offs = gp_offset(i)

                if(igr) then
                    E_IP = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do q = 1, gp_M(i)
                        stencil_idx = offs + q
                        u = gp_stencil(1, stencil_idx)
                        v = gp_stencil(2, stencil_idx)
                        neum = gp_A_neum(stencil_idx)
                        dirich = gp_A_dirich(stencil_idx)
                        !rho_inv = 1.0_wp / q_cons_vf(contxb)%sf(u, v, 0)
                        vel_IP(1) = vel_IP(1) + dirich * q_cons_vf(momxb)%sf(u, v, 0) / q_cons_vf(contxb)%sf(u, v, 0)
                        vel_IP(2) = vel_IP(2) + dirich * q_cons_vf(momxb + 1)%sf(u, v, 0) / q_cons_vf(contxb)%sf(u, v, 0)

                        rho = rho + neum * q_cons_vf(contxb)%sf(u, v, 0)
                        if ()
                        E_IP = E_IP + neum * q_cons_vf(E_idx)%sf(u, v, 0)

                    end do
                    q_cons_vf(momxb)%sf(j, k, 0) = vel_IP(1) * rho
                    q_cons_vf(momxb + 1)%sf(j, k, 0) = vel_IP(2) * rho
                    q_cons_vf(contxb)%sf(j, k, 0) = rho
                    q_cons_vf(E_idx)%sf(j, k, 0) = E_IP

                    ! TEMPORARY DEBUGGING SECTION
                    if(patch_ib(1)%print_cond) then
                        q_cons_vf(E_idx)%sf(j, k, 0) = ghost_points(i)%cond

                        if (mod(i, int(num_gps / 30)) == 0) then
                            do q = 1, ghost_points(i)%M
                                u = ghost_points(i)%stencil(q, 1)
                                v = ghost_points(i)%stencil(q, 2)
                                q_cons_vf(contxb)%sf(u, v, 0) = i
                            end do
                        end if
                    end if
                else
                    pres_IP = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do q = 1, gp_M(i)
                        stencil_idx = offs + q
                        u = gp_stencil(1, stencil_idx)
                        v = gp_stencil(2, stencil_idx)
                        neum = gp_A_neum(stencil_idx)
                        dirich = gp_A_dirich(stencil_idx)
                        vel_IP(1) = vel_IP(1) + dirich * q_prim_vf(momxb)%sf(u, v, 0)
                        vel_IP(2) = vel_IP(2) + dirich * q_prim_vf(momxb + 1)%sf(u, v, 0)
                        rho = rho + neum * q_prim_vf(contxb)%sf(u, v, 0)
                        pres_IP = pres_IP + neum * q_prim_vf(E_idx)%sf(u, v, 0)
                    end do
                    q_prim_vf(contxb)%sf(j, k, 0) = rho
                    q_prim_vf(momxb)%sf(j, k, 0) = vel_IP(1)
                    q_prim_vf(momxb + 1)%sf(j, k, 0) = vel_IP(2)
                    q_prim_vf(E_idx)%sf(j, k, 0) = pres_IP
                    q_cons_vf(contxb)%sf(j, k, 0) = q_prim_vf(contxb)%sf(j, k, 0)
                    q_cons_vf(momxb)%sf(j, k, 0) = q_prim_vf(contxb)%sf(j, k, 0) * q_prim_vf(momxb)%sf(j, k, 0)
                    q_cons_vf(momxb + 1)%sf(j, k, 0) = q_prim_vf(contxb)%sf(j, k, 0) * q_prim_vf(momxb + 1)%sf(j, k, 0)
                    q_cons_vf(E_idx)%sf(j, k, 0) = q_prim_vf(E_idx)%sf(j, k, 0) * 2.5_wp + 0.5_wp * q_prim_vf(contxb)%sf(j, k, 0) * (q_prim_vf(momxb)%sf(j, k, 0) ** 2 + q_prim_vf(momxb + 1)%sf(j, k, 0) ** 2)
                end if
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        !Correct the state of the inner points in IBs
        ! if (num_inner_gps > 0) then
        !     $:GPU_PARALLEL_LOOP(private='[i,j,k,l,q,innerp]')
        !     do i = 1, num_inner_gps
        !
        !         innerp = inner_points(i)
        !         j = innerp%loc(1)
        !         k = innerp%loc(2)
        !         l = innerp%loc(3)
        !
        !         $:GPU_LOOP(parallelism='[seq]')
        !         do q = momxb, momxe
        !             q_cons_vf(q)%sf(j, k, l) = 0._wp
        !         end do
        !     end do
        !     $:END_GPU_PARALLEL_LOOP()
        ! end if

    end subroutine s_ibm_correct_state

    !>  Function that computes the image points for each ghost point
        !!  @param ghost_points_in Ghost Points
    subroutine s_compute_image_points(gp)
        type(ghost_point), intent(INOUT) :: gp

        real(wp) :: dist
        real(wp), dimension(3) :: norm
        real(wp), dimension(3) :: physical_loc
        real(wp) :: temp_loc
        real(wp), pointer, dimension(:) :: s_cc => null()
        integer :: bound

        integer :: q, dim !< Iterator variables
        integer :: i, j, k, l !< Location indexes
        integer :: patch_id !< IB Patch ID
        integer :: dir
        integer :: index
        logical :: bounds_error

        bounds_error = .false.

        i = gp%loc(1)
        j = gp%loc(2)
        k = gp%loc(3)

        ! Calculate physical location of ghost point
        if (p > 0) then
            physical_loc = [x_cc(i), y_cc(j), z_cc(k)]
        else
            physical_loc = [x_cc(i), y_cc(j), 0._wp]
        end if

        ! Calculate and store the precise location of the image point
        patch_id = gp%ib_patch_id
        dist = abs(real(gp%levelset, kind=wp))
        norm(:) = gp%levelset_norm
        gp%ip_loc(:) = physical_loc(:) + 2*dist*norm(:)

        ! Find the closest grid point to the image point
        do dim = 1, num_dims

            ! s_cc points to the dim array we need
            if (dim == 1) then
                s_cc => x_cc
                bound = m + buff_size - 1
            elseif (dim == 2) then
                s_cc => y_cc
                bound = n + buff_size - 1
            else
                s_cc => z_cc
                bound = p + buff_size - 1
            end if

            if (f_approx_equal(norm(dim), 0._wp)) then
                ! if the ghost point is almost equal to a cell location, we set it equal and continue
                gp%ip_grid(dim) = gp%loc(dim)
            else
                if (norm(dim) > 0) then
                    dir = 1
                else
                    dir = -1
                end if

                index = gp%loc(dim)
                temp_loc = gp%ip_loc(dim)
                do while ((temp_loc < s_cc(index) &
                            .or. temp_loc > s_cc(index + 1)) .and. (.not. bounds_error))
                    index = index + dir
                    if (index < -buff_size .or. index > bound) then
#if !defined(MFC_OpenACC) && !defined(MFC_OpenMP)
                        print *, "A required image point is not located in this computational domain."
                        print *, "Ghost Point is located at :"
                        if (p == 0) then
                            print *, [x_cc(i), y_cc(j)]
                        else
                            print *, [x_cc(i), y_cc(j), z_cc(k)]
                        end if
                        print *, "We are searching in dimension ", dim, " for image point at ", gp%ip_loc(:)
                        print *, "Domain size: ", [x_cc(-buff_size), y_cc(-buff_size), z_cc(-buff_size)]
                        print *, "x: ", x_cc(-buff_size), " to: ", x_cc(m + buff_size - 1)
                        print *, "y: ", y_cc(-buff_size), " to: ", y_cc(n + buff_size - 1)
                        if (p /= 0) print *, "z: ", z_cc(-buff_size), " to: ", z_cc(p + buff_size - 1)
                        print *, "Image point is located approximately ", (gp%loc(dim) - gp%ip_loc(dim))/(s_cc(1) - s_cc(0)), " grid cells away"
                        print *, "Levelset ", dist, " and Norm: ", norm(:)
                        print *, "A short term fix may include increasing buff_size further in m_helper_basic (currently set to a minimum of 10)"
#endif
                        bounds_error = .true.
                    end if
                end do

                gp%ip_grid(dim) = index
                if (gp%DB(dim) == -1) then
                    gp%ip_grid(dim) = gp%loc(dim) + 1
                else if (gp%DB(dim) == 1) then
                    gp%ip_grid(dim) = gp%loc(dim) - 1
                end if
            end if
        end do

        if (bounds_error) error stop "Ghost Point and Image Point on Different Processors. Exiting"

    end subroutine s_compute_image_points

    !> Subroutine that finds the number of ghost points, used for allocating
    !! memory.
    subroutine s_find_num_ghost_points(num_gps_out, num_inner_gps_out)

        integer, intent(out) :: num_gps_out
        integer, intent(out) :: num_inner_gps_out

        integer :: i, j, k, ii, jj, kk, gp_layers_z !< Iterator variables
        integer :: num_gps_local, num_inner_gps_local !< local copies of the gp count to support GPU compute
        logical :: is_gp

        num_gps_local = 0
        num_inner_gps_local = 0
        gp_layers_z = gp_layers
        if (p == 0) gp_layers_z = 0

        $:GPU_PARALLEL_LOOP(private='[i,j,k,ii,jj,kk,is_gp]', copy='[num_gps_local,num_inner_gps_local]', firstprivate='[gp_layers,gp_layers_z]', collapse=3)
        do i = 0, m
            do j = 0, n
                do k = 0, p
                    if (ib_markers%sf(i, j, k) /= 0) then
                        is_gp = .false.
                        marker_search: do ii = i - gp_layers, i + gp_layers
                            do jj = j - gp_layers, j + gp_layers
                                do kk = k - gp_layers_z, k + gp_layers_z
                                    if (ib_markers%sf(ii, jj, kk) == 0) then
                                        ! if any neighbors are not in the IB, it is a ghost point
                                        is_gp = .true.
                                        exit marker_search
                                    end if
                                end do
                            end do
                        end do marker_search

                        if (is_gp) then
                            $:GPU_ATOMIC(atomic='update')
                            num_gps_local = num_gps_local + 1
                        else
                            $:GPU_ATOMIC(atomic='update')
                            num_inner_gps_local = num_inner_gps_local + 1
                        end if
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        num_gps_out = num_gps_local
        num_inner_gps_out = num_inner_gps_local

    end subroutine s_find_num_ghost_points

    !> Function that finds the ghost points
    subroutine s_find_ghost_points(ghost_points_in, inner_points_in)

        type(ghost_point), dimension(num_gps), intent(INOUT) :: ghost_points_in
        type(ghost_point), dimension(num_inner_gps), intent(INOUT) :: inner_points_in
        integer :: i, j, k, ii, jj, kk, gp_layers_z !< Iterator variables
        integer :: xp, yp, zp !< periodicities
        integer :: count, count_i, local_idx
        integer :: patch_id, encoded_patch_id
        logical :: is_gp

        count = 0
        count_i = 0
        gp_layers_z = gp_layers
        if (p == 0) gp_layers_z = 0

        $:GPU_PARALLEL_LOOP(private='[i,j,k,ii,jj,kk,is_gp,local_idx,patch_id,encoded_patch_id,xp,yp,zp]', copyin='[count,count_i, x_domain, y_domain, z_domain]', firstprivate='[gp_layers,gp_layers_z]', collapse=3)
        do i = 0, m
            do j = 0, n
                do k = 0, p
                    if (ib_markers%sf(i, j, k) /= 0) then
                        is_gp = .false.
                        marker_search: do ii = i - gp_layers, i + gp_layers
                            do jj = j - gp_layers, j + gp_layers
                                do kk = k - gp_layers_z, k + gp_layers_z
                                    if (ib_markers%sf(ii, jj, kk) == 0) then
                                        ! if any neighbors are not in the IB, it is a ghost point
                                        is_gp = .true.
                                        exit marker_search
                                    end if
                                end do
                            end do
                        end do marker_search

                        if (is_gp) then
                            $:GPU_ATOMIC(atomic='capture')
                            count = count + 1
                            local_idx = count
                            $:END_GPU_ATOMIC_CAPTURE()

                            ghost_points_in(local_idx)%loc = [i, j, k]
                            encoded_patch_id = ib_markers%sf(i, j, k)
                            call s_decode_patch_periodicity(encoded_patch_id, patch_id, xp, yp, zp)
                            ghost_points_in(local_idx)%ib_patch_id = patch_id
                            ib_markers%sf(i, j, k) = patch_id
                            ghost_points_in(local_idx)%x_periodicity = xp
                            ghost_points_in(local_idx)%y_periodicity = yp
                            ghost_points_in(local_idx)%z_periodicity = zp
                            ghost_points_in(local_idx)%slip = patch_ib(patch_id)%slip

                            ghost_points_in(local_idx)%first_layer = .false.
                            first_layer: do ii = i - 1, i + 1
                                do jj = j - 1, j + 1
                                    do kk = k - 1, k + 1
                                        if (ib_markers%sf(ii, jj, kk) == 0) then
                                            ! if any neighbors are not in the IB, it is a ghost point
                                            ghost_points_in(local_idx)%first_layer = .true.
                                            exit first_layer
                                        end if
                                    end do
                                end do
                            end do first_layer

                            if ((x_cc(i) - dx(i)) < x_domain%beg) then
                                ghost_points_in(local_idx)%DB(1) = -1
                            else if ((x_cc(i) + dx(i)) > x_domain%end) then
                                ghost_points_in(local_idx)%DB(1) = 1
                            else
                                ghost_points_in(local_idx)%DB(1) = 0
                            end if

                            if ((y_cc(j) - dy(j)) < y_domain%beg) then
                                ghost_points_in(local_idx)%DB(2) = -1
                            else if ((y_cc(j) + dy(j)) > y_domain%end) then
                                ghost_points_in(local_idx)%DB(2) = 1
                            else
                                ghost_points_in(local_idx)%DB(2) = 0
                            end if

                            if (p /= 0) then
                                if ((z_cc(k) - dz(k)) < z_domain%beg) then
                                    ghost_points_in(local_idx)%DB(3) = -1
                                else if ((z_cc(k) + dz(k)) > z_domain%end) then
                                    ghost_points_in(local_idx)%DB(3) = 1
                                else
                                    ghost_points_in(local_idx)%DB(3) = 0
                                end if
                            end if

                        else
                            $:GPU_ATOMIC(atomic='capture')
                            count_i = count_i + 1
                            local_idx = count_i
                            $:END_GPU_ATOMIC_CAPTURE()

                            inner_points_in(local_idx)%loc = [i, j, k]
                            encoded_patch_id = ib_markers%sf(i, j, k)
                            call s_decode_patch_periodicity(encoded_patch_id, patch_id, xp, yp, zp)
                            inner_points_in(local_idx)%ib_patch_id = patch_id
                            ib_markers%sf(i, j, k) = patch_id
                            inner_points_in(local_idx)%slip = patch_ib(patch_id)%slip

                        end if
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_find_ghost_points

    subroutine s_radial_search(gp)
        type(ghost_point), intent(inout) :: gp
        integer :: p
        integer :: count, maxCount, i, j, k, r_step, i0, j0, q, l, col
        real(wp) :: R, r_loc, r_prev, dr
        real(wp) :: x_p, y_p, z_p
        real(wp), dimension(3) :: center
        integer :: M, order, LN
        real(wp), dimension(200) :: weights
        real(wp), allocatable, dimension(:) :: gp_poly
        real(wp), dimension(:, :), allocatable :: WV, WV_PINV, A
        real(wp) :: cond, maxCond
        real(wp) :: finalRad, rad_ratio

        maxCond = 0._wp
        order = patch_ib(1)%high_order

        select case (order)
            case(1)
                LN = 3
            case(2)
                LN = 6
            case(3)
                LN = 10
            case(4)
                LN = 15
        end select
        rad_ratio = 0.0_wp

        allocate(gp_poly(LN))
        ! center = [x_cc(gp%loc(1)), y_cc(gp%loc(2)), 0._wp]
        center = [gp%bound_loc(1), gp%bound_loc(2), 0._wp]
        !maxCount = int(10 * patch_ib(1)%smooth)
        maxCount = int(LN * patch_ib(1)%pts)
        i0 = gp%loc(1)
        j0 = gp%loc(2)
        !k0 = gp%loc(3)
        dr = dx(i0) ** 2 + dy(j0) ** 2
        if (num_dims == 3) then
        !    dr = dr + dz(k0) ** 2
        end if
        dr = sqrt(dr)
        !dr = 0.3_wp * dr
        dr = patch_ib(1)%dr_fac * dr

        ! entry 1 is ghost cell
        gp%stencil(1, :) = [i0, j0, 0]
        count = 1

        r_prev = -1e-6_wp
        R = dr
        r_step = 10
        finalRad = -99.0_wp

        do while (count < maxCount)
            do i = i0 - r_step, i0 + r_step
                do j = j0 - r_step, j0 + r_step
                    r_loc = sqrt((x_cc(i) - center(1))**2 + (y_cc(j) - center(2))**2)
                    if (r_loc <= R .and. r_loc > r_prev) then
                        if (ib_markers%sf(i, j, 0) == 0) then
                            count = count + 1
                            gp%stencil(count, :) = [i, j, 0]
                            finalRad = max(finalRad, r_loc)
                        end if
                    end if
                end do
            end do
            r_prev = R
            R = R + dr
        end do
        !print *, "FOUND STENCIL"
        finalRad = max(finalRad, sqrt((x_cc(i0) - center(1)) ** 2 + (y_cc(j0) - center(2)) ** 2))

        R = finalRad * 1.05_wp
        rad_ratio = max(rad_ratio, finalRad / R)

        M = count
        allocate(WV(M, LN))
        allocate(WV_PINV(LN, M))

        ! compute weights
        do q = 1, M
            r_loc = (x_cc(gp%stencil(q, 1)) - center(1)) ** 2 &
            + (y_cc(gp%stencil(q, 2)) - center(2)) ** 2
            if (num_dims == 3) then
                r_loc = r_loc + (z_cc(gp%stencil(q, 3)) - center(3)) ** 2
                end if
                r_loc = sqrt(r_loc)
                weights(q) = 0.5_wp * (1 + cos((pi * r_loc / R)))
            end do
            weights(1) = 1.0e6_wp
            ! Construct WV
            do q = 1, M
                x_p = x_cc(gp%stencil(q, 1)) - center(1) / R
                y_p = y_cc(gp%stencil(q, 2)) - center(2) / R
                if (num_dims == 3) then
                    z_p = z_cc(gp%stencil(q, 3)) - center(3)
                end if
                col = 1
                do l = 0, order
                    do i = l, 0, -1
                        j = l - i
                        if (q == 1) then
                            if (col == 2) then
                                WV(q, col) = weights(q) * gp%levelset_norm(1) / R
                            else if (col == 3) then
                                WV(q, col) = weights(q) * gp%levelset_norm(2) / R
                            else 
                                WV(q, col) = 0._wp
                            end if
                            gp_poly(col) = x_p ** i * y_p ** j
                        else
                            WV(q, col) = weights(q) * x_p ** i * y_p ** j
                        end if
                        col = col + 1
                    end do
                end do
            end do

            call find_p_inv(M, LN, WV, WV_PINV, cond)

            maxCond = max(maxCond, cond)
            do q = 1, M
                gp%A_temp(1, q) = weights(q) * dot_product(gp_poly, WV_PINV(:, q))
            end do
            gp%M = M
            gp%cond = cond

            do q = 1, M
                x_p = (x_cc(gp%stencil(q, 1)) - center(1)) / R
                y_p = (y_cc(gp%stencil(q, 2)) - center(2)) / R
                if (num_dims == 3) then
                    z_p = z_cc(gp%stencil(q, 3)) - center(3)
                end if
                col = 1
                do l = 0, order
                    do i = l, 0, -1
                        j = l - i
                        if (q == 1) then
                            if(col == 1) then
                                WV(q, col) = weights(q)
                            else
                                WV(q, col) = 0.0_wp
                            end if
                        else
                            WV(q, col) = weights(q) * x_p ** i * y_p ** j
                        end if
                        col = col + 1
                    end do
                end do
            end do

            call find_p_inv(M, LN, WV, WV_PINV, cond)

            do q = 1, M
                gp%A_temp(2, q) = weights(q) * dot_product(gp_poly, WV_PINV(:, q))
            end do
            gp%A_temp(1, 1) = 0._wp
            gp%A_temp(2, 1) = 0._wp

        deallocate(WV)
        deallocate(WV_PINV)

    end subroutine s_radial_search

    subroutine find_p_inv(M, N, A, PINV,cond)
        integer, intent(in) :: M, N
        real(wp), dimension(M, N), intent(inout) :: A
        real(wp), dimension(N, M), intent(out) :: PINV
        integer :: K, LDU, LDVT, LDA, LWORK, INFO, j
        real(wp) :: TOL
        real(wp), dimension(:), allocatable :: S, WORK
        real(wp), dimension(:, :), allocatable :: U, VT
        real(wp), dimension(1) :: DUMMY_WORK
        real(wp), intent(out) :: cond

        LDA = M
        K = min(M, N)
        LDU = M
        LDVT = K

        allocate(S(K))
        allocate(U(LDU, K))
        allocate(VT(LDVT, N))

        LWORK = -1
        call dgesvd('S', 'S', M, N, A, LDA, S, U, LDU, VT, LDVT, DUMMY_WORK, LWORK, INFO)

        LWORK = int(DUMMY_WORK(1))
        allocate(WORK(LWORK))

        call dgesvd('S', 'S', M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO)
        if (INFO /= 0) then
            print *, "SVD DID NOT CONVERGE", INFO
        end if
        
        tol = max(M, N) * epsilon(1._wp) * S(1)

        do j = 1, K
            if (S(j) > tol) then
                U(1:M, j) = U(1:M, j) / S(j)
            else
                U(1:M, j) = 0._wp
            end if
        end do

        call dgemm('T', 'T', N, M, K, 1._wp, VT, LDVT, U, LDU, 0._wp, PINV, N)

        cond = S(1) / S(K)
        deallocate(S, U, VT, WORK)
    end subroutine find_p_inv

    !>  Function that computes the interpolation coefficients of image points
    subroutine s_compute_interpolation_coeffs(gp)

        type(ghost_point), intent(INOUT) :: gp

        real(wp), dimension(2, 2, 2) :: dist
        real(wp), dimension(2, 2, 2) :: alpha
        real(wp), dimension(2, 2, 2) :: interp_coeffs
        real(wp) :: buf
        real(wp), dimension(2, 2, 2) :: eta
        integer :: q, i, j, k, ii, jj, kk !< Grid indexes and iterators
        integer :: patch_id
        logical is_cell_center

        ! Get the interpolation points
        i = gp%ip_grid(1)
        j = gp%ip_grid(2)
        if (p /= 0) then
            k = gp%ip_grid(3)
        else
            k = 0; 
        end if

        ! get the distance to a cell in each direction
        dist = 0._wp
        buf = 1._wp
        do ii = 0, 1
            do jj = 0, 1
                if (p == 0) then
                    dist(1 + ii, 1 + jj, 1) = sqrt( &
                                                (x_cc(i + ii) - gp%ip_loc(1))**2 + &
                                                (y_cc(j + jj) - gp%ip_loc(2))**2)
                else
                    do kk = 0, 1
                        dist(1 + ii, 1 + jj, 1 + kk) = sqrt( &
                                                        (x_cc(i + ii) - gp%ip_loc(1))**2 + &
                                                        (y_cc(j + jj) - gp%ip_loc(2))**2 + &
                                                        (z_cc(k + kk) - gp%ip_loc(3))**2)
                    end do
                end if
            end do
        end do

        ! check if we are arbitrarily close to a cell center
        interp_coeffs = 0._wp
        is_cell_center = .false.
        check_is_cell_center: do ii = 0, 1
            do jj = 0, 1
                if (dist(ii + 1, jj + 1, 1) <= 1.e-16_wp) then
                    interp_coeffs(ii + 1, jj + 1, 1) = 1._wp
                    is_cell_center = .true.
                    exit check_is_cell_center
                else
                    if (p /= 0) then
                        if (dist(ii + 1, jj + 1, 2) <= 1.e-16_wp) then
                            interp_coeffs(ii + 1, jj + 1, 2) = 1._wp
                            is_cell_center = .true.
                            exit check_is_cell_center
                        end if
                    end if
                end if
            end do
        end do check_is_cell_center

        if (.not. is_cell_center) then
            ! if we are not arbitrarily close, interpolate
            alpha = 1._wp
            patch_id = gp%ib_patch_id
            if (ib_markers%sf(i, j, k) /= 0) alpha(1, 1, 1) = 0._wp
            if (ib_markers%sf(i + 1, j, k) /= 0) alpha(2, 1, 1) = 0._wp
            if (ib_markers%sf(i, j + 1, k) /= 0) alpha(1, 2, 1) = 0._wp
            if (ib_markers%sf(i + 1, j + 1, k) /= 0) alpha(2, 2, 1) = 0._wp
            if (p == 0) then
                eta(:, :, 1) = 1._wp/dist(:, :, 1)**2
                buf = sum(alpha(:, :, 1)*eta(:, :, 1))
                if (buf > 0._wp) then
                    interp_coeffs(:, :, 1) = alpha(:, :, 1)*eta(:, :, 1)/buf
                else
                    buf = sum(eta(:, :, 1))
                    interp_coeffs(:, :, 1) = eta(:, :, 1)/buf
                end if
            else

                if (ib_markers%sf(i, j, k + 1) /= 0) alpha(1, 1, 2) = 0._wp
                if (ib_markers%sf(i + 1, j, k + 1) /= 0) alpha(2, 1, 2) = 0._wp
                if (ib_markers%sf(i, j + 1, k + 1) /= 0) alpha(1, 2, 2) = 0._wp
                if (ib_markers%sf(i + 1, j + 1, k + 1) /= 0) alpha(2, 2, 2) = 0._wp
                eta = 1._wp/dist**2
                buf = sum(alpha*eta)

                if (buf > 0._wp) then
                    interp_coeffs = alpha*eta/buf
                else
                    buf = sum(eta)
                    interp_coeffs = eta/buf
                end if
            end if

        end if

        gp%interp_coeffs = interp_coeffs

    end subroutine s_compute_interpolation_coeffs

    !> Function that uses the interpolation coefficients and the current state
    !! at the cell centers in order to estimate the state at the image point
    !! @param gp Ghost point data structure
    !> @brief Interpolates primitive variables from the fluid domain to a ghost point's image point using bilinear or trilinear interpolation.
    !! @param alpha_rho_IP Partial density at image point
    !! @param alpha_IP Volume fraction at image point
    !! @param pres_IP Pressure at image point
    !! @param vel_IP Velocity at image point
    !! @param c_IP Speed of sound at image point
    !! @param r_IP Bubble radius at image point
    !! @param v_IP Bubble radial velocity at image point
    !! @param pb_IP Bubble pressure at image point
    !! @param mv_IP Bubble vapor mass at image point
    !! @param nmom_IP Bubble moment at image point
    !! @param pb_in Internal bubble pressure array
    !! @param mv_in Mass of vapor in bubble array
    !! @param presb_IP Bubble node pressure at image point
    !! @param massv_IP Bubble node vapor mass at image point
    subroutine s_interpolate_image_point(q_vf, gp, alpha_rho_IP, alpha_IP, &
                                         pres_IP, vel_IP, c_IP, r_IP, v_IP, pb_IP, &
                                         mv_IP, nmom_IP, pb_in, mv_in, presb_IP, massv_IP, E_IP)
        type(scalar_field), &
            dimension(sys_size), &
            intent(IN) :: q_vf !< State Variables

        real(stp), optional, dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:, 1:), intent(IN) :: pb_in, mv_in

        type(ghost_point), intent(IN) :: gp
        real(wp), intent(INOUT) :: pres_IP
        real(wp), optional, intent(INOUT) :: E_IP !< IGR uses conservative variables
        real(wp), dimension(3), intent(INOUT) :: vel_IP
        real(wp), intent(INOUT) :: c_IP
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3), intent(INOUT) :: alpha_IP, alpha_rho_IP
        #:else
            real(wp), dimension(num_fluids), intent(INOUT) :: alpha_IP, alpha_rho_IP
        #:endif
        real(wp), optional, dimension(:), intent(INOUT) :: r_IP, v_IP, pb_IP, mv_IP
        real(wp), optional, dimension(:), intent(INOUT) :: nmom_IP
        real(wp), optional, dimension(:), intent(INOUT) :: presb_IP, massv_IP

        integer :: i, j, k, l, q !< Iterator variables
        integer :: i1, i2, j1, j2, k1, k2 !< Iterator variables
        real(wp) :: coeff, alpha_sum, rho

        i1 = gp%ip_grid(1); i2 = i1 + 1
        j1 = gp%ip_grid(2); j2 = j1 + 1
        k1 = gp%ip_grid(3); k2 = k1 + 1

        if (p == 0) then
            k1 = 0
            k2 = 0
        end if

        alpha_rho_IP = 0._wp
        alpha_IP = 0._wp
        pres_IP = 0._wp
        vel_IP = 0._wp

        if (igr) E_IP = 0._wp

        if (surface_tension) c_IP = 0._wp

        if (bubbles_euler) then
            r_IP = 0._wp
            v_IP = 0._wp
            if (.not. polytropic) then
                mv_IP = 0._wp
                pb_IP = 0._wp
            end if
        end if

        if (qbmm) then
            nmom_IP = 0._wp
            if (.not. polytropic) then
                presb_IP = 0._wp
                massv_IP = 0._wp
            end if
        end if

        $:GPU_LOOP(parallelism='[seq]')
        do i = i1, i2
            $:GPU_LOOP(parallelism='[seq]')
            do j = j1, j2
                $:GPU_LOOP(parallelism='[seq]')
                do k = k1, k2

                    coeff = gp%interp_coeffs(i - i1 + 1, j - j1 + 1, k - k1 + 1)

                    if (igr) then
                        E_IP = E_IP + coeff* &
                            q_vf(E_idx)%sf(i, j, k)
                        rho = 0._wp
                        if (num_fluids == 1) then
                            alpha_rho_IP(1) = alpha_rho_IP(1) + coeff*q_vf(contxb)%sf(i, j, k)
                            alpha_IP(1) = alpha_IP(1) + coeff
                            rho = q_vf(contxb)%sf(i, j, k)
                        else
                            alpha_sum = 0._wp
                            $:GPU_LOOP(parallelism='[seq]')
                            do l = 1, num_fluids - 1
                                alpha_rho_IP(l) = alpha_rho_IP(l) + coeff*q_vf(contxb + l - 1)%sf(i, j, k)
                                alpha_IP(l) = alpha_IP(l) + coeff*q_vf(E_idx + l)%sf(i, j, k)
                                alpha_sum = alpha_sum + q_vf(E_idx + l)%sf(i, j, k)
                                rho = rho + q_vf(l)%sf(i, j, k)
                            end do
                            alpha_rho_IP(num_fluids) = alpha_rho_IP(num_fluids) + coeff*q_vf(num_fluids)%sf(i, j, k)
                            alpha_IP(num_fluids) = alpha_IP(num_fluids) + coeff*(1._wp - alpha_sum)
                            rho = rho + q_vf(num_fluids)%sf(i, j, k)
                        end if
                        $:GPU_LOOP(parallelism='[seq]')
                        do q = momxb, momxe
                            vel_IP(q + 1 - momxb) = vel_IP(q + 1 - momxb) + coeff* &
                                                    q_vf(q)%sf(i, j, k) / rho
                        end do

                    else
                        pres_IP = pres_IP + coeff* &
                                q_vf(E_idx)%sf(i, j, k)

                        $:GPU_LOOP(parallelism='[seq]')
                        do q = momxb, momxe
                            vel_IP(q + 1 - momxb) = vel_IP(q + 1 - momxb) + coeff* &
                                                    q_vf(q)%sf(i, j, k)
                        end do
                        
                        $:GPU_LOOP(parallelism='[seq]')
                        do l = contxb, contxe
                            alpha_rho_IP(l) = alpha_rho_IP(l) + coeff* &
                                            q_vf(l)%sf(i, j, k)
                            alpha_IP(l) = alpha_IP(l) + coeff* &
                                        q_vf(advxb + l - 1)%sf(i, j, k)
                        end do
                    end if

                    if (surface_tension) then
                        c_IP = c_IP + coeff*q_vf(c_idx)%sf(i, j, k)
                    end if

                    if (bubbles_euler .and. .not. qbmm) then
                        $:GPU_LOOP(parallelism='[seq]')
                        do l = 1, nb
                            if (polytropic) then
                                r_IP(l) = r_IP(l) + coeff*q_vf(bubxb + (l - 1)*2)%sf(i, j, k)
                                v_IP(l) = v_IP(l) + coeff*q_vf(bubxb + 1 + (l - 1)*2)%sf(i, j, k)
                            else
                                r_IP(l) = r_IP(l) + coeff*q_vf(bubxb + (l - 1)*4)%sf(i, j, k)
                                v_IP(l) = v_IP(l) + coeff*q_vf(bubxb + 1 + (l - 1)*4)%sf(i, j, k)
                                pb_IP(l) = pb_IP(l) + coeff*q_vf(bubxb + 2 + (l - 1)*4)%sf(i, j, k)
                                mv_IP(l) = mv_IP(l) + coeff*q_vf(bubxb + 3 + (l - 1)*4)%sf(i, j, k)
                            end if
                        end do
                    end if

                    if (qbmm) then
                        do l = 1, nb*nmom
                            nmom_IP(l) = nmom_IP(l) + coeff*q_vf(bubxb - 1 + l)%sf(i, j, k)
                        end do
                        if (.not. polytropic) then
                            do q = 1, nb
                                do l = 1, nnode
                                    presb_IP((q - 1)*nnode + l) = presb_IP((q - 1)*nnode + l) + &
                                                                  coeff*real(pb_in(i, j, k, l, q), kind=wp)
                                    massv_IP((q - 1)*nnode + l) = massv_IP((q - 1)*nnode + l) + &
                                                                  coeff*real(mv_in(i, j, k, l, q), kind=wp)
                                end do
                            end do
                        end if
                    end if
                end do
            end do
        end do

    end subroutine s_interpolate_image_point

    !> Resets the current indexes of immersed boundaries and replaces them after updating
    !> the position of each moving immersed boundary
    impure subroutine s_update_mib(num_ibs)

        integer, intent(in) :: num_ibs

        integer :: i, j, k, ierr, z_gp_layers

        call nvtxStartRange("UPDATE-MIBM")

        ! Clears the existing immersed boundary indices
        z_gp_layers = 0; if (p /= 0) z_gp_layers = gp_layers + 1
        $:GPU_PARALLEL_LOOP(private='[i,j,k]')
        do i = -gp_layers - 1, m + gp_layers + 1; do j = -gp_layers - 1, n + gp_layers + 1; do k = -z_gp_layers, p + z_gp_layers
                    ib_markers%sf(i, j, k) = 0._wp
                end do; end do; end do
        $:END_GPU_PARALLEL_LOOP()

        ! recalulcate the rotation matrix based upon the new angles
        do i = 1, num_ibs
            if (patch_ib(i)%moving_ibm /= 0) then
                call s_update_ib_rotation_matrix(i)
            end if
        end do

        $:GPU_UPDATE(device='[patch_ib]')

        ! recompute the new ib_patch locations and broadcast them.
        call nvtxStartRange("APPLY-IB-PATCHES")
        call s_apply_ib_patches(ib_markers)
        call nvtxEndRange

        call nvtxStartRange("COMPUTE-GHOST-POINTS")
        ! recalculate the ghost point locations and coefficients
        call s_find_num_ghost_points(num_gps, num_inner_gps)
        call s_find_ghost_points(ghost_points, inner_points)
        call nvtxEndRange

        call nvtxStartRange("COMPUTE-IMAGE-POINTS")
        call s_apply_levelset(ghost_points, num_gps)
            do i = 1, num_gps
                call s_compute_image_points(ghost_points(i))
                call s_compute_interpolation_coeffs(ghost_points(i))
            end do
        call nvtxEndRange

        call nvtxEndRange

    end subroutine s_update_mib

    subroutine s_compute_ib_forces(q_prim_vf, fluid_pp)

        ! real(wp), dimension(idwbuff(1)%beg:idwbuff(1)%end, &
        !             idwbuff(2)%beg:idwbuff(2)%end, &
        !             idwbuff(3)%beg:idwbuff(3)%end), intent(in) :: pressure
        type(scalar_field), dimension(1:sys_size), intent(in) :: q_prim_vf
        type(physical_parameters), dimension(1:num_fluids), intent(in) :: fluid_pp

        integer :: gp_id, i, j, k, l, q, ib_idx, fluid_idx
        real(wp), dimension(num_ibs, 3) :: forces, torques
        real(wp), dimension(num_ibs, 3) :: forces_viscous
        real(wp), dimension(1:3, 1:3) :: viscous_stress_div, viscous_stress_div_1, viscous_stress_div_2, viscous_cross_1, viscous_cross_2 ! viscous stress tensor with temp vectors to hold divergence calculations
        real(wp), dimension(1:3) :: local_force_contribution, radial_vector, local_torque_contribution, vel
        real(wp), dimension(1:3) :: local_force_contribution_viscous
        real(wp) :: cell_volume, dx, dy, dz, dynamic_viscosity
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3) :: dynamic_viscosities
        #:else
            real(wp), dimension(num_fluids) :: dynamic_viscosities
        #:endif
        real(wp), dimension(1:3, 2) :: pres, dynPres
        character(len=50) :: file_loc
        call nvtxStartRange("COMPUTE-IB-FORCES")

        forces = 0._wp
        forces_viscous = 0._wp
        torques = 0._wp

        if (viscous) then
            do fluid_idx = 1, num_fluids
                if (fluid_pp(fluid_idx)%Re(1) /= 0._wp) then
                    dynamic_viscosities(fluid_idx) = 1._wp/fluid_pp(fluid_idx)%Re(1)
                else
                    dynamic_viscosities(fluid_idx) = 0._wp
                end if
            end do
        end if

        $:GPU_PARALLEL_LOOP(private='[l,pres,dynPres,ib_idx,fluid_idx,radial_vector,local_force_contribution,local_force_contribution_viscous,cell_volume,local_torque_contribution, dynamic_viscosity, viscous_stress_div, viscous_stress_div_1, viscous_stress_div_2, viscous_cross_1, viscous_cross_2, dx, dy, dz]', copy='[forces,torques,forces_viscous]', copyin='[ib_markers,patch_ib,dynamic_viscosities,fluid_pp(1)]', collapse=3)
        do i = 0, m
            do j = 0, n
                do k = 0, p
                    ib_idx = ib_markers%sf(i, j, k)
                    if (ib_idx /= 0) then
                        ! get the vector pointing to the grid cell from the IB centroid
                        if (num_dims == 3) then
                            radial_vector = [x_cc(i), y_cc(j), z_cc(k)] - [patch_ib(ib_idx)%x_centroid, patch_ib(ib_idx)%y_centroid, patch_ib(ib_idx)%z_centroid]
                        else
                            radial_vector = [x_cc(i), y_cc(j), 0._wp] - [patch_ib(ib_idx)%x_centroid, patch_ib(ib_idx)%y_centroid, 0._wp]
                        end if
                        dx = x_cc(i + 1) - x_cc(i)
                        dy = y_cc(j + 1) - y_cc(j)

                        local_force_contribution(:) = 0._wp
                        local_force_contribution_viscous(:) = 0._wp

                        if(igr .or. dummy) then
                            ! Q_PRIM_VF IS CONSERVATIVE WHEN IGR
                            pres(:, :) = 0._wp
             
                            !pres(1, 1) = s_compute_pressure_igr(q_prim_vf, i - 1, j, k, pres(1, 1))
                            pres(1, 1) = (q_prim_vf(E_idx)%sf(i-1,j,k) - 0.5_wp / q_prim_vf(contxb)%sf(i-1,j,k) &
                                            * (q_prim_vf(momxb)%sf(i-1,j,k) * q_prim_vf(momxb)%sf(i-1,j,k) + &
                                            q_prim_vf(momxb + 1)%sf(i-1,j,k) * q_prim_vf(momxb + 1)%sf(i-1,j,k))) / fluid_pp(1)%gamma

                            pres(1, 2) = (q_prim_vf(E_idx)%sf(i+1,j,k) - 0.5_wp / q_prim_vf(contxb)%sf(i+1,j,k) &
                                            * (q_prim_vf(momxb)%sf(i+1,j,k) * q_prim_vf(momxb)%sf(i+1,j,k) + &
                                            q_prim_vf(momxb + 1)%sf(i+1,j,k) * q_prim_vf(momxb + 1)%sf(i+1,j,k))) / fluid_pp(1)%gamma

                            pres(2, 1) = (q_prim_vf(E_idx)%sf(i,j-1,k) - 0.5_wp / q_prim_vf(contxb)%sf(i,j-1,k) &
                                            * (q_prim_vf(momxb)%sf(i,j-1,k) * q_prim_vf(momxb)%sf(i,j-1,k) + &
                                            q_prim_vf(momxb + 1)%sf(i,j-1,k) * q_prim_vf(momxb + 1)%sf(i,j-1,k))) / fluid_pp(1)%gamma
                            pres(2, 2) = (q_prim_vf(E_idx)%sf(i,j+1,k) - 0.5_wp / q_prim_vf(contxb)%sf(i,j+1,k) &
                                            * (q_prim_vf(momxb)%sf(i,j+1,k) * q_prim_vf(momxb)%sf(i,j+1,k) + &
                                            q_prim_vf(momxb + 1)%sf(i,j+1,k) * q_prim_vf(momxb + 1)%sf(i,j+1,k))) / fluid_pp(1)%gamma
                        end if

                        do fluid_idx = 0, num_fluids - 1
                            ! Get the pressure contribution to force via a finite difference to compute the 2D components of the gradient of the pressure and cell volume
                            if(igr .or. dummy) then
                                local_force_contribution(1) = local_force_contribution(1) - (pres(1,2) - pres(1,1))/(2._wp*dx) ! force is the negative pressure gradient
                                local_force_contribution(2) = local_force_contribution(2) - (pres(2,2) - pres(2,1))/(2._wp*dy)
                            else
                                local_force_contribution(1) = local_force_contribution(1) - (q_prim_vf(E_idx + fluid_idx)%sf(i + 1, j, k) - q_prim_vf(E_idx + fluid_idx)%sf(i - 1, j, k))/(2._wp*dx) ! force is the negative pressure gradient
                                local_force_contribution(2) = local_force_contribution(2) - (q_prim_vf(E_idx + fluid_idx)%sf(i, j + 1, k) - q_prim_vf(E_idx + fluid_idx)%sf(i, j - 1, k))/(2._wp*dy)
                            end if
                            cell_volume = abs(dx*dy)
                            ! add the 3D component of the pressure gradient, if we are working in 3 dimensions
                            if (num_dims == 3) then
                                dz = z_cc(k + 1) - z_cc(k)
                                local_force_contribution(3) = local_force_contribution(3) - (q_prim_vf(E_idx + fluid_idx)%sf(i, j, k + 1) - q_prim_vf(E_idx + fluid_idx)%sf(i, j, k - 1))/(2._wp*dz)
                                cell_volume = abs(cell_volume*dz)
                            end if
                        end do

                        ! Update the force values atomically to prevent race conditions
!                        call s_cross_product(radial_vector, local_force_contribution, local_torque_contribution)

                        ! get the viscous stress and add its contribution if that is considered
                        ! TODO :: This is really bad code
                        if (viscous) then
                            ! compute the volume-weighted local dynamic viscosity
                            dynamic_viscosity = 0._wp
                            do fluid_idx = 1, num_fluids
                                ! local dynamic viscosity is the dynamic viscosity of the fluid times alpha of the fluid
                                if(igr .or. dummy) then
                                    dynamic_viscosity = dynamic_viscosity + 1.0_wp * dynamic_viscosities(fluid_idx)
                                else
                                    dynamic_viscosity = dynamic_viscosity + (q_prim_vf(fluid_idx + advxb - 1)%sf(i, j, k)*dynamic_viscosities(fluid_idx))
                                end if
                            end do

                            ! get the linear force component first
                            call s_compute_viscous_stress_tensor(viscous_stress_div_1, q_prim_vf, dynamic_viscosity, i - 1, j, k)
                            call s_compute_viscous_stress_tensor(viscous_stress_div_2, q_prim_vf, dynamic_viscosity, i + 1, j, k)
                            viscous_stress_div = (viscous_stress_div_2 - viscous_stress_div_1)/(2._wp*dx) ! get the x derivative of the viscous stress tensor
                            local_force_contribution_viscous(1:3) = local_force_contribution_viscous(1:3) + viscous_stress_div(1, 1:3) ! add te x components of the derivative to the force
                            do l = 1, 3
                                ! take the cross products for the torque component
                                !call s_cross_product(radial_vector, viscous_stress_div_1(l, 1:3), viscous_cross_1(l, 1:3))
                                !call s_cross_product(radial_vector, viscous_stress_div_2(l, 1:3), viscous_cross_2(l, 1:3))
                            end do

                            !viscous_stress_div = (viscous_cross_2 - viscous_cross_1)/(2._wp*dx) ! get the x derivative of the cross product
                            !local_torque_contribution(1:3) = local_torque_contribution(1:3) + viscous_stress_div(1, 1:3) ! apply the cross product derivative to the torque

                            call s_compute_viscous_stress_tensor(viscous_stress_div_1, q_prim_vf, dynamic_viscosity, i, j - 1, k)
                            call s_compute_viscous_stress_tensor(viscous_stress_div_2, q_prim_vf, dynamic_viscosity, i, j + 1, k)
                            viscous_stress_div = (viscous_stress_div_2 - viscous_stress_div_1)/(2._wp*dy)
                            local_force_contribution_viscous(1:3) = local_force_contribution_viscous(1:3) + viscous_stress_div(2, 1:3)
                            do l = 1, 3
                                !    call s_cross_product(radial_vector, viscous_stress_div_1(l, 1:3), viscous_cross_1(l, 1:3))
                                !    call s_cross_product(radial_vector, viscous_stress_div_2(l, 1:3), viscous_cross_2(l, 1:3))
                            end do

                            !viscous_stress_div = (viscous_cross_2 - viscous_cross_1)/(2._wp*dy)
                            !local_torque_contribution(1:3) = local_torque_contribution(1:3) + viscous_stress_div(2, 1:3)

                            if (num_dims == 3) then
                                call s_compute_viscous_stress_tensor(viscous_stress_div_1, q_prim_vf, dynamic_viscosity, i, j, k - 1)
                                call s_compute_viscous_stress_tensor(viscous_stress_div_2, q_prim_vf, dynamic_viscosity, i, j, k + 1)
                                viscous_stress_div = (viscous_stress_div_2 - viscous_stress_div_1)/(2._wp*dz)
                                local_force_contribution_viscous(1:3) = local_force_contribution_viscous(1:3) + viscous_stress_div(3, 1:3)
                                do l = 1, 3
                                    !        call s_cross_product(radial_vector, viscous_stress_div_1(l, 1:3), viscous_cross_1(l, 1:3))
                                    !        call s_cross_product(radial_vector, viscous_stress_div_2(l, 1:3), viscous_cross_2(l, 1:3))
                                end do
                                !    viscous_stress_div = (viscous_cross_2 - viscous_cross_1)/(2._wp*dz)
                                !    local_torque_contribution(1:3) = local_torque_contribution(1:3) + viscous_stress_div(3, 1:3)
                            end if
                        end if

                        do l = 1, 3
                            $:GPU_ATOMIC(atomic='update')
                            forces(ib_idx, l) = forces(ib_idx, l) + (local_force_contribution(l)*cell_volume)
                            $:GPU_ATOMIC(atomic='update')
                            forces_viscous(ib_idx, l) = forces_viscous(ib_idx, l) + (local_force_contribution_viscous(l)*cell_volume)
                            !$:GPU_ATOMIC(atomic='update')
!                           torques(ib_idx, l) = torques(ib_idx, l) + local_torque_contribution(l)*cell_volume
                        end do
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! reduce the forces across all MPI ranks
        call s_mpi_allreduce_vectors_sum(forces, forces, num_ibs, 3)
        call s_mpi_allreduce_vectors_sum(forces_viscous, forces_viscous, num_ibs, 3)
!        call s_mpi_allreduce_vectors_sum(torques, torques, num_ibs, 3)

        ! consider body forces after reducing to avoid double counting
        do i = 1, num_ibs
            if (bf_x) then
                forces(i, 1) = forces(i, 1) + accel_bf(1)*patch_ib(i)%mass
            end if
            if (bf_y) then
                forces(i, 2) = forces(i, 2) + accel_bf(2)*patch_ib(i)%mass
            end if
            if (bf_z) then
                forces(i, 3) = forces(i, 3) + accel_bf(3)*patch_ib(i)%mass
            end if
        end do

        ! apply the summed forces
        do i = 1, num_ibs
!            patch_ib(i)%force(:) = forces(i, :)
!            patch_ib(i)%torque(:) = matmul(patch_ib(i)%rotation_matrix_inverse, torques(i, :)) ! torques must be converted to the local coordinates of the IB
        end do

        call nvtxEndRange
        if (proc_rank == 0) then
            file_loc = trim(case_dir)//'/forces.csv'
            open (unit=92, file=trim(file_loc), status='UNKNOWN', position='APPEND', action='WRITE')
            write (92, '(ES16.8, A1, ES16.8, A1, ES16.8, A1, ES16.8, A1, ES16.8)') mytime, ",", forces(1, 1), ",", forces_viscous(1, 1), ",", forces(1, 2), ",", forces_viscous(1, 2)
            close (92)
        end if

    end subroutine s_compute_ib_forces

   !> @brief Computes pressure and viscous forces and torques on immersed bodies via a volume integration method.
    ! subroutine s_compute_ib_forces(q_prim_vf, fluid_pp)
    !
    !     ! real(wp), dimension(idwbuff(1)%beg:idwbuff(1)%end, &
    !     !             idwbuff(2)%beg:idwbuff(2)%end, &
    !     !             idwbuff(3)%beg:idwbuff(3)%end), intent(in) :: pressure
    !     type(scalar_field), dimension(1:sys_size), intent(in) :: q_prim_vf
    !     type(physical_parameters), dimension(1:num_fluids), intent(in) :: fluid_pp
    !
    !     integer :: gp_id, i, j, k, l, q, ib_idx, fluid_idx
    !     real(wp), dimension(num_ibs, 3) :: forces, torques
    !     real(wp), dimension(1:3, 1:3) :: viscous_stress_div, viscous_stress_div_1, viscous_stress_div_2, viscous_cross_1, viscous_cross_2 ! viscous stress tensor with temp vectors to hold divergence calculations
    !     real(wp), dimension(1:3) :: local_force_contribution, radial_vector, local_torque_contribution, vel
    !     real(wp), dimension(2) :: dynamic_viscosity
    !     real(wp) :: cell_volume, dx, dy, dz
    !     #:if not MFC_CASE_OPTIMIZATION and USING_AMD
    !         real(wp), dimension(3) :: dynamic_viscosities
    !     #:else
    !         real(wp), dimension(num_fluids) :: dynamic_viscosities
    !     #:endif
    !     forces = 0._wp
    !     torques = 0._wp
    !
    !     if (viscous) then
    !         do fluid_idx = 1, num_fluids
    !             if (fluid_pp(fluid_idx)%Re(1) /= 0._wp) then
    !                 dynamic_viscosities(fluid_idx) = 1._wp/fluid_pp(fluid_idx)%Re(1)
    !             else
    !                 dynamic_viscosities(fluid_idx) = 0._wp
    !             end if
    !         end do
    !     end if
    !
    !     $:GPU_PARALLEL_LOOP(private='[l,ib_idx,fluid_idx, radial_vector,local_force_contribution,cell_volume,local_torque_contribution, dynamic_viscosity, viscous_stress_div, viscous_stress_div_1, viscous_stress_div_2, viscous_cross_1, viscous_cross_2, dx, dy, dz]', copy='[forces,torques]', copyin='[ib_markers,patch_ib,dynamic_viscosities]', collapse=3)
    !     do i = 0, m
    !         do j = 0, n
    !             do k = 0, p
    !                 ib_idx = ib_markers%sf(i, j, k)
    !                 if (ib_idx /= 0) then
    !                     ! get the vector pointing to the grid cell from the IB centroid
    !                     if (num_dims == 3) then
    !                         radial_vector = [x_cc(i), y_cc(j), z_cc(k)] - [patch_ib(ib_idx)%x_centroid, patch_ib(ib_idx)%y_centroid, patch_ib(ib_idx)%z_centroid]
    !                     else
    !                         radial_vector = [x_cc(i), y_cc(j), 0._wp] - [patch_ib(ib_idx)%x_centroid, patch_ib(ib_idx)%y_centroid, 0._wp]
    !                     end if
    !                     dx = x_cc(i + 1) - x_cc(i)
    !                     dy = y_cc(j + 1) - y_cc(j)
    !
    !                     local_force_contribution(:) = 0._wp
    !                     local_torque_contribution(:) = 0._wp
    !                     do fluid_idx = 0, num_fluids - 1
    !                         ! Get the pressure contribution to force via a finite difference to compute the 2D components of the gradient of the pressure and cell volume
    !                         local_force_contribution(1) = local_force_contribution(1) - (q_prim_vf(E_idx + fluid_idx)%sf(i + 1, j, k) - q_prim_vf(E_idx + fluid_idx)%sf(i - 1, j, k))/(2._wp*dx) ! force is the negative pressure gradient
    !                         local_force_contribution(2) = local_force_contribution(2) - (q_prim_vf(E_idx + fluid_idx)%sf(i, j + 1, k) - q_prim_vf(E_idx + fluid_idx)%sf(i, j - 1, k))/(2._wp*dy)
    !                         cell_volume = abs(dx*dy)
    !                         ! add the 3D component of the pressure gradient, if we are working in 3 dimensions
    !                         if (num_dims == 3) then
    !                             dz = z_cc(k + 1) - z_cc(k)
    !                             local_force_contribution(3) = local_force_contribution(3) - (q_prim_vf(E_idx + fluid_idx)%sf(i, j, k + 1) - q_prim_vf(E_idx + fluid_idx)%sf(i, j, k - 1))/(2._wp*dz)
    !                             cell_volume = abs(cell_volume*dz)
    !                         end if
    !                     end do
    !
    !                     call s_cross_product(radial_vector, local_force_contribution, local_torque_contribution)
    !
    !                     ! get the viscous stress and add its contribution if that is considered
    !                     ! TODO :: This is really bad code
    !                     if (viscous) then
    !                         ! compute the volume-weighted local dynamic viscosity
    !                         dynamic_viscosity = 0._wp
    !                         do fluid_idx = 1, num_fluids
    !                             ! local dynamic viscosity is the dynamic viscosity of the fluid times alpha of the fluid
    !                             dynamic_viscosity(1) = dynamic_viscosity(1) + (q_prim_vf(fluid_idx + advxb - 1)%sf(i - 1, j, k)*dynamic_viscosities(fluid_idx))
    !                             dynamic_viscosity(2) = dynamic_viscosity(2) + (q_prim_vf(fluid_idx + advxb - 1)%sf(i + 1, j, k)*dynamic_viscosities(fluid_idx))
    !                         end do
    !
    !                         ! get the linear force component first
    !                         call s_compute_viscous_stress_tensor(viscous_stress_div_1, q_prim_vf, dynamic_viscosity(1), i - 1, j, k)
    !                         call s_compute_viscous_stress_tensor(viscous_stress_div_2, q_prim_vf, dynamic_viscosity(2), i + 1, j, k)
    !                         viscous_stress_div = (viscous_stress_div_2 - viscous_stress_div_1)/(2._wp*dx) ! get the x derivative of the viscous stress tensor
    !                         local_force_contribution(1:3) = local_force_contribution(1:3) + viscous_stress_div(1, 1:3) ! add te x components of the derivative to the force
    !                         do l = 1, 3
    !                             ! take the cross products for the torque component
    !                             call s_cross_product(radial_vector, viscous_stress_div_1(l, 1:3), viscous_cross_1(l, 1:3))
    !                             call s_cross_product(radial_vector, viscous_stress_div_2(l, 1:3), viscous_cross_2(l, 1:3))
    !                         end do
    !
    !                         viscous_stress_div = (viscous_cross_2 - viscous_cross_1)/(2._wp*dx) ! get the x derivative of the cross product
    !                         local_torque_contribution(1:3) = local_torque_contribution(1:3) + viscous_stress_div(1, 1:3) ! apply the cross product derivative to the torque
    !
    !                         dynamic_viscosity = 0._wp
    !                         do fluid_idx = 1, num_fluids
    !                             dynamic_viscosity(1) = dynamic_viscosity(1) + (q_prim_vf(fluid_idx + advxb - 1)%sf(i, j - 1, k)*dynamic_viscosities(fluid_idx))
    !                             dynamic_viscosity(2) = dynamic_viscosity(2) + (q_prim_vf(fluid_idx + advxb - 1)%sf(i, j + 1, k)*dynamic_viscosities(fluid_idx))
    !                         end do
    !
    !                         call s_compute_viscous_stress_tensor(viscous_stress_div_1, q_prim_vf, dynamic_viscosity(1), i, j - 1, k)
    !                         call s_compute_viscous_stress_tensor(viscous_stress_div_2, q_prim_vf, dynamic_viscosity(2), i, j + 1, k)
    !                         viscous_stress_div = (viscous_stress_div_2 - viscous_stress_div_1)/(2._wp*dy)
    !                         local_force_contribution(1:3) = local_force_contribution(1:3) + viscous_stress_div(2, 1:3)
    !                         do l = 1, 3
    !                             call s_cross_product(radial_vector, viscous_stress_div_1(l, 1:3), viscous_cross_1(l, 1:3))
    !                             call s_cross_product(radial_vector, viscous_stress_div_2(l, 1:3), viscous_cross_2(l, 1:3))
    !                         end do
    !
    !                         viscous_stress_div = (viscous_cross_2 - viscous_cross_1)/(2._wp*dy)
    !                         local_torque_contribution(1:3) = local_torque_contribution(1:3) + viscous_stress_div(2, 1:3)
    !
    !                         if (num_dims == 3) then
    !                             dynamic_viscosity = 0._wp
    !                             do fluid_idx = 1, num_fluids
    !                                 dynamic_viscosity(1) = dynamic_viscosity(1) + (q_prim_vf(fluid_idx + advxb - 1)%sf(i, j, k - 1)*dynamic_viscosities(fluid_idx))
    !                                 dynamic_viscosity(2) = dynamic_viscosity(2) + (q_prim_vf(fluid_idx + advxb - 1)%sf(i, j, k + 1)*dynamic_viscosities(fluid_idx))
    !                             end do
    !
    !                             call s_compute_viscous_stress_tensor(viscous_stress_div_1, q_prim_vf, dynamic_viscosity(1), i, j, k - 1)
    !                             call s_compute_viscous_stress_tensor(viscous_stress_div_2, q_prim_vf, dynamic_viscosity(2), i, j, k + 1)
    !                             viscous_stress_div = (viscous_stress_div_2 - viscous_stress_div_1)/(2._wp*dz)
    !                             local_force_contribution(1:3) = local_force_contribution(1:3) + viscous_stress_div(3, 1:3)
    !                             do l = 1, 3
    !                                 call s_cross_product(radial_vector, viscous_stress_div_1(l, 1:3), viscous_cross_1(l, 1:3))
    !                                 call s_cross_product(radial_vector, viscous_stress_div_2(l, 1:3), viscous_cross_2(l, 1:3))
    !                             end do
    !                             viscous_stress_div = (viscous_cross_2 - viscous_cross_1)/(2._wp*dz)
    !                             local_torque_contribution(1:3) = local_torque_contribution(1:3) + viscous_stress_div(3, 1:3)
    !                         end if
    !                     end if
    !
    !                     ! Update the force values atomically to prevent race conditions
    !                     do l = 1, 3
    !                         $:GPU_ATOMIC(atomic='update')
    !                         forces(ib_idx, l) = forces(ib_idx, l) + local_force_contribution(l)*cell_volume
    !                         $:GPU_ATOMIC(atomic='update')
    !                         torques(ib_idx, l) = torques(ib_idx, l) + local_torque_contribution(l)*cell_volume
    !                     end do
    !                 end if
    !             end do
    !         end do
    !     end do
    !     $:END_GPU_PARALLEL_LOOP()
    !
    !     ! reduce the forces across all MPI ranks
    !     call s_mpi_allreduce_vectors_sum(forces, forces, num_ibs, 3)
    !     call s_mpi_allreduce_vectors_sum(torques, torques, num_ibs, 3)
    !
    !     ! consider body forces after reducing to avoid double counting
    !     do i = 1, num_ibs
    !         if (bf_x) then
    !             forces(i, 1) = forces(i, 1) + accel_bf(1)*patch_ib(i)%mass
    !         end if
    !         if (bf_y) then
    !             forces(i, 2) = forces(i, 2) + accel_bf(2)*patch_ib(i)%mass
    !         end if
    !         if (bf_z) then
    !             forces(i, 3) = forces(i, 3) + accel_bf(3)*patch_ib(i)%mass
    !         end if
    !     end do
    !
    !     ! apply the summed forces
    !     do i = 1, num_ibs
    !         patch_ib(i)%force(:) = forces(i, :)
    !         patch_ib(i)%torque(:) = matmul(patch_ib(i)%rotation_matrix_inverse, torques(i, :)) ! torques must be converted to the local coordinates of the IB
    !     end do
    !
    ! end subroutine s_compute_ib_forces

    !> Subroutine to deallocate memory reserved for the IBM module
    impure subroutine s_finalize_ibm_module()

        @:DEALLOCATE(ib_markers%sf)
        if (allocated(airfoil_grid_u)) then
            @:DEALLOCATE(airfoil_grid_u)
            @:DEALLOCATE(airfoil_grid_l)
        end if

        @:DEALLOCATE(gp_loc, gp_ib_patch_id, gp_offset, gp_M)
        @:DEALLOCATE(gp_stencil, gp_A_dirich, gp_A_neum)

    end subroutine s_finalize_ibm_module

    !> Computes the center of mass for IB patch types where we are unable to determine their center of mass analytically.
    !> These patches include things like NACA airfoils and STL models
    subroutine s_compute_centroid_offset(ib_marker)

        integer, intent(in) :: ib_marker

        integer :: i, j, k, num_cells, num_cells_local
        real(wp), dimension(1:3) :: center_of_mass, center_of_mass_local

        ! Offset only needs to be computes for specific geometries
        if (patch_ib(ib_marker)%geometry == 4 .or. &
            patch_ib(ib_marker)%geometry == 5 .or. &
            patch_ib(ib_marker)%geometry == 11 .or. &
            patch_ib(ib_marker)%geometry == 12) then

            center_of_mass_local = [0._wp, 0._wp, 0._wp]
            num_cells_local = 0

            ! get the summed mass distribution and number of cells to divide by
            do i = 0, m
                do j = 0, n
                    do k = 0, p
                        if (ib_markers%sf(i, j, k) == ib_marker) then
                            num_cells_local = num_cells_local + 1
                            center_of_mass_local = center_of_mass_local + [x_cc(i), y_cc(j), 0._wp]
                            if (num_dims == 3) center_of_mass_local(3) = center_of_mass_local(3) + z_cc(k)
                        end if
                    end do
                end do
            end do

            ! reduce the mass contribution over all MPI ranks and compute COM
            call s_mpi_allreduce_integer_sum(num_cells_local, num_cells)
            if (num_cells /= 0) then
                call s_mpi_allreduce_sum(center_of_mass_local(1), center_of_mass(1))
                call s_mpi_allreduce_sum(center_of_mass_local(2), center_of_mass(2))
                call s_mpi_allreduce_sum(center_of_mass_local(3), center_of_mass(3))
                center_of_mass = center_of_mass/real(num_cells, wp)
            else
                patch_ib(ib_marker)%centroid_offset = [0._wp, 0._wp, 0._wp]
                return
            end if

            ! assign the centroid offset as a vector pointing from the true COM to the "centroid" in the input file and replace the current centroid
            patch_ib(ib_marker)%centroid_offset = [patch_ib(ib_marker)%x_centroid, patch_ib(ib_marker)%y_centroid, patch_ib(ib_marker)%z_centroid] &
                                                  - center_of_mass
            patch_ib(ib_marker)%x_centroid = center_of_mass(1)
            patch_ib(ib_marker)%y_centroid = center_of_mass(2)
            patch_ib(ib_marker)%z_centroid = center_of_mass(3)

            ! rotate the centroid offset back into the local coords of the IB
            patch_ib(ib_marker)%centroid_offset = matmul(patch_ib(ib_marker)%rotation_matrix_inverse, patch_ib(ib_marker)%centroid_offset)
        else
            patch_ib(ib_marker)%centroid_offset(:) = [0._wp, 0._wp, 0._wp]
        end if

    end subroutine s_compute_centroid_offset

    !>  Computes the moment of inertia for an immersed boundary
        !!  @param ib_marker Immersed boundary marker index
    subroutine s_compute_moment_of_inertia(ib_marker, axis)

        real(wp), dimension(3), intent(in) :: axis !< the axis about which we compute the moment. Only required in 3D.
        integer, intent(in) :: ib_marker

        real(wp) :: moment, distance_to_axis, cell_volume
        real(wp), dimension(3) :: position, closest_point_along_axis, vector_to_axis, normal_axis
        integer :: i, j, k, count

        if (p == 0) then
            normal_axis = [0, 0, 1]
        else if (sqrt(sum(axis**2)) == 0) then
            ! if the object is not actually rotating at this time, return a dummy value and exit
            patch_ib(ib_marker)%moment = 1._wp
            return
        else
            normal_axis = axis/sqrt(sum(axis))
        end if

        ! if the IB is in 2D or a 3D sphere, we can compute this exactly
        if (patch_ib(ib_marker)%geometry == 2) then ! circle
            patch_ib(ib_marker)%moment = 0.5_wp*patch_ib(ib_marker)%mass*(patch_ib(ib_marker)%radius)**2
        elseif (patch_ib(ib_marker)%geometry == 3) then ! rectangle
            patch_ib(ib_marker)%moment = patch_ib(ib_marker)%mass*(patch_ib(ib_marker)%length_x**2 + patch_ib(ib_marker)%length_y**2)/6._wp
        elseif (patch_ib(ib_marker)%geometry == 6) then ! ellipse
            patch_ib(ib_marker)%moment = 0.0625_wp*patch_ib(ib_marker)%mass*(patch_ib(ib_marker)%length_x**2 + patch_ib(ib_marker)%length_y**2)
        elseif (patch_ib(ib_marker)%geometry == 8) then ! sphere
            patch_ib(ib_marker)%moment = 0.4*patch_ib(ib_marker)%mass*(patch_ib(ib_marker)%radius)**2

        else ! we do not have an analytic moment of inertia calculation and need to approximate it directly via a sum
            count = 0
            moment = 0._wp
            cell_volume = (x_cc(1) - x_cc(0))*(y_cc(1) - y_cc(0)) ! computed without grid stretching. Update in the loop to perform with stretching
            if (p /= 0) then
                cell_volume = cell_volume*(z_cc(1) - z_cc(0))
            end if

            $:GPU_PARALLEL_LOOP(private='[position,closest_point_along_axis,vector_to_axis,distance_to_axis]', copy='[moment,count]', copyin='[ib_marker,cell_volume,normal_axis]', collapse=3)
            do i = 0, m
                do j = 0, n
                    do k = 0, p
                        if (ib_markers%sf(i, j, k) == ib_marker) then
                            $:GPU_ATOMIC(atomic='update')
                            count = count + 1 ! increment the count of total cells in the boundary

                            ! get the position in local coordinates so that the axis passes through 0, 0, 0
                            if (p == 0) then
                                position = [x_cc(i), y_cc(j), 0._wp] - [patch_ib(ib_marker)%x_centroid, patch_ib(ib_marker)%y_centroid, 0._wp]
                            else
                                position = [x_cc(i), y_cc(j), z_cc(k)] - [patch_ib(ib_marker)%x_centroid, patch_ib(ib_marker)%y_centroid, patch_ib(ib_marker)%z_centroid]
                            end if

                            ! project the position along the axis to find the closest distance to the rotation axis
                            closest_point_along_axis = normal_axis*dot_product(normal_axis, position)
                            vector_to_axis = position - closest_point_along_axis
                            distance_to_axis = dot_product(vector_to_axis, vector_to_axis) ! saves the distance to the axis squared

                            ! compute the position component of the moment
                            $:GPU_ATOMIC(atomic='update')
                            moment = moment + distance_to_axis
                        end if
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

            ! write the final moment assuming the points are all uniform density
            patch_ib(ib_marker)%moment = moment*patch_ib(ib_marker)%mass/(count*cell_volume)
            $:GPU_UPDATE(device='[patch_ib(ib_marker)%moment]')
        end if

    end subroutine s_compute_moment_of_inertia

    !> @brief Checks for periodic boundary conditions in all directions, and if so, moves patch location if it left the domain
    subroutine s_wrap_periodic_ibs()

        integer :: patch_id

        do patch_id = 1, num_ibs
            ! check domain wraps in x, y,
            #:for X in [('x'), ('y')]
                ! check for periodicity
                if (bc_${X}$%beg == BC_PERIODIC) then
                    ! check if the boundary has left the domain, and then correct
                    if (patch_ib(patch_id)%${X}$_centroid < ${X}$_domain%beg) then
                        ! if the boundary exited "left", wrap it back around to the "right"
                        patch_ib(patch_id)%${X}$_centroid = patch_ib(patch_id)%${X}$_centroid + (${X}$_domain%end - ${X}$_domain%beg)
                    else if (patch_ib(patch_id)%${X}$_centroid > ${X}$_domain%end) then
                        ! if the boundary exited "right", wrap it back around to the "left"
                        patch_ib(patch_id)%${X}$_centroid = patch_ib(patch_id)%${X}$_centroid - (${X}$_domain%end - ${X}$_domain%beg)
                    end if
                end if
            #:endfor

            if (p /= 0) then
                ! check for periodicity
                if (bc_z%beg == BC_PERIODIC) then
                    ! check if the boundary has left the domain, and then correct
                    if (patch_ib(patch_id)%z_centroid < z_domain%beg) then
                        ! if the boundary exited "left", wrap it back around to the "right"
                        patch_ib(patch_id)%z_centroid = patch_ib(patch_id)%z_centroid + (z_domain%end - z_domain%beg)
                    else if (patch_ib(patch_id)%z_centroid > z_domain%end) then
                        ! if the boundary exited "right", wrap it back around to the "left"
                        patch_ib(patch_id)%z_centroid = patch_ib(patch_id)%z_centroid - (z_domain%end - z_domain%beg)
                    end if
                end if
            end if
        end do

    end subroutine s_wrap_periodic_ibs

    !> @brief Computes the cross product c = a x b of two 3D vectors.
    subroutine s_cross_product(a, b, c)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: a(3), b(3)
        real(wp), intent(out) :: c(3)

        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)
    end subroutine s_cross_product

end module m_ibm
