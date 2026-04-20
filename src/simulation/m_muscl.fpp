!>
!! @file
!! @brief Contains module m_muscl

#:include 'macros.fpp'

!> @brief MUSCL reconstruction with interface sharpening for contact-preserving advection
module m_muscl

    use m_derived_types
    use m_global_parameters
    use m_variables_conversion
#ifdef MFC_OpenACC
    use openacc
#endif

    use m_mpi_proxy
    use m_helper

    private; public :: s_initialize_muscl_module, s_muscl, s_finalize_muscl_module, s_pack_muscl_input_arr, s_interface_compression

    real(wp), allocatable, dimension(:,:,:,:) :: v_rs_muscl
    $:GPU_DECLARE(create='[v_rs_muscl]')

    type(int_bounds_info) :: is1_muscl, is2_muscl, is3_muscl
    $:GPU_DECLARE(create='[is1_muscl, is2_muscl, is3_muscl]')

    integer :: v_size
    $:GPU_DECLARE(create='[v_size]')

    !> @name The cell-average variables that will be MUSCL-reconstructed. Formerly, they are stored in v_vf. However, they are
    !! transferred to v_rs_wsL and v_rs_wsR as to be reshaped (RS) and/or characteristically decomposed. The reshaping allows the
    !! muscl procedure to be independent of the coordinate direction of the reconstruction. Lastly, notice that the left (L) and
    !! right (R) results of the characteristic decomposition are stored in custom-constructed muscl- stencils (WS) that are annexed
    !! to each position of a given scalar field.
    !> @{
    real(wp), allocatable, dimension(:,:,:,:) :: v_rs_ws_x_muscl, v_rs_ws_y_muscl, v_rs_ws_z_muscl
    !> @}
    $:GPU_DECLARE(create='[v_rs_ws_x_muscl, v_rs_ws_y_muscl, v_rs_ws_z_muscl]')

contains

    !> Allocate and initialize MUSCL reconstruction working arrays
    subroutine s_initialize_muscl_module()

        @:ALLOCATE(v_rs_muscl(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end, 1:sys_size))

    end subroutine s_initialize_muscl_module

    subroutine s_pack_muscl_input_arr(v_vf)

        type(scalar_field), dimension(1:), intent(in) :: v_vf
        integer                                       :: i, j, k, l

        $:GPU_PARALLEL_LOOP(collapse=4)
        do i = 1, sys_size
            do l = idwbuff(3)%beg, idwbuff(3)%end
                do k = idwbuff(2)%beg, idwbuff(2)%end
                    do j = idwbuff(1)%beg, idwbuff(1)%end
                        v_rs_muscl(j, k, l, i) = v_vf(i)%sf(j, k, l)
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_pack_muscl_input_arr

    !> Perform MUSCL reconstruction of left and right cell-boundary values from cell-averaged variables
    subroutine s_muscl(vL_rs_vf_x, vL_rs_vf_y, vL_rs_vf_z, vR_rs_vf_x, vR_rs_vf_y, vR_rs_vf_z, muscl_dir)

        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(inout) :: vL_rs_vf_x, vL_rs_vf_y, &
             & vL_rs_vf_z, vR_rs_vf_x, vR_rs_vf_y, vR_rs_vf_z
        integer, intent(in)               :: muscl_dir
        integer                           :: j, k, l, i
        real(wp)                          :: slopeL, slopeR, slope
        real(wp)                          :: v0

        v_size = ubound(vL_rs_vf_x, 1)
        $:GPU_UPDATE(device='[v_size]')

        if (muscl_order == 1 .or. dummy) then
            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, v_size
                do l = idwbuff(3)%beg, idwbuff(3)%end
                    do k = idwbuff(2)%beg, idwbuff(2)%end
                        do j = idwbuff(1)%beg, idwbuff(1)%end
                            vL_rs_x(j, k, l, i) = v_rs_muscl(j, k, l, i)
                            vR_rs_x(j, k, l, i) = v_rs_muscl(j, k, l, i)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        if (muscl_order == 2 .or. dummy) then
            ! MUSCL Reconstruction
            #:for MUSCL_DIR, XYZ, STENCIL_VAR, COORDS, X_BND_OFFS, Y_BND_OFFS, Z_BND_OFFS in &
                [(1, 'x', 'j', '{STENCIL_IDX}, k, l', 'muscl_polyn', '0', '0'), &
                 (2, 'y', 'k', 'j, {STENCIL_IDX}, l', '0', 'muscl_polyn', '0'), &
                 (3, 'z', 'l', 'j, k, {STENCIL_IDX}', '0', '0', 'muscl_polyn')]
                #:set SV = STENCIL_VAR
                #:set SF = lambda offs: COORDS.format(STENCIL_IDX = SV + offs)
                if (muscl_dir == ${MUSCL_DIR}$) then
                    $:GPU_PARALLEL_LOOP(collapse=4,private='[i, j, k, l, slopeL, slopeR, slope, v0]',copyin='[v_size]')
                    do l = idwbuff(3)%beg + ${Z_BND_OFFS}$, idwbuff(3)%end - ${Z_BND_OFFS}$
                        do k = idwbuff(2)%beg + ${Y_BND_OFFS}$, idwbuff(2)%end - ${Y_BND_OFFS}$
                            do j = idwbuff(1)%beg + ${X_BND_OFFS}$, idwbuff(1)%end - ${X_BND_OFFS}$
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, v_size
                                    v0 = v_rs_muscl(${SF('')}$, i)
                                    slopeL = v_rs_muscl(${SF(' + 1')}$, i) - v0
                                    slopeR = v0 - v_rs_muscl(${SF(' - 1')}$, i)
                                    slope = 0._wp

                                    if (muscl_lim == 1) then  ! minmod
                                        if (slopeL*slopeR > 1e-9_wp) then
                                            slope = min(abs(slopeL), abs(slopeR))
                                        end if
                                        if (slopeL < 0._wp) slope = -slope
                                    else if (muscl_lim == 2) then  ! MC
                                        if (slopeL*slopeR > 1e-9_wp) then
                                            slope = min(2._wp*abs(slopeL), 2._wp*abs(slopeR))
                                            slope = min(slope, 5e-1_wp*(abs(slopeL) + abs(slopeR)))
                                        end if
                                        if (slopeL < 0._wp) slope = -slope
                                    else if (muscl_lim == 3) then  ! Van Albada
                                        if (abs(slopeL) > 1e-6_wp .and. abs(slopeR) > 1e-6_wp .and. abs(slopeL + slopeR) &
                                            & > 1e-6_wp .and. slopeL*slopeR > 1e-6_wp) then
                                            slope = ((slopeL + slopeR)*slopeL*slopeR)/(slopeL**2._wp + slopeR**2._wp)
                                        end if
                                    else if (muscl_lim == 4) then  ! Van Leer
                                        if (abs(slopeL + slopeR) > 1.e-6_wp .and. slopeL*slopeR > 1.e-6_wp) then
                                            slope = 2._wp*slopeL*slopeR/(slopeL + slopeR)
                                        end if
                                    else if (muscl_lim == 5) then  ! SUPERBEE
                                        if (slopeL*slopeR > 1e-6_wp) then
                                            slope = -1._wp*min(-min(2._wp*abs(slopeL), abs(slopeR)), -min(abs(slopeL), &
                                                               & 2._wp*abs(slopeR)))
                                        end if
                                    end if

                                    ! reconstruct from left side
                                    vL_rs_vf_${XYZ}$ (j, k, l, i) = v0 - (5.e-1_wp*slope)

                                    ! reconstruct from the right side
                                    vR_rs_vf_${XYZ}$ (j, k, l, i) = v0 + (5.e-1_wp*slope)
                                end do
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end if
            #:endfor
        end if

        !> :TODO Apply the same NO RESHAPE changes to MUSCL before uncommenting this
        !if (int_comp) then
        !    call s_interface_compression(vL_rs_vf_x, vL_rs_vf_y, vL_rs_vf_z, vR_rs_vf_x, vR_rs_vf_y, vR_rs_vf_z, muscl_dir, &
        !                                 & is1_muscl_d, is2_muscl_d, is3_muscl_d)
        !end if

    end subroutine s_muscl

    !> Apply THINC interface-compression to sharpen volume-fraction reconstructions at material interfaces
    subroutine s_interface_compression(vL_rs_vf_x, vL_rs_vf_y, vL_rs_vf_z, vR_rs_vf_x, vR_rs_vf_y, vR_rs_vf_z, muscl_dir, &

        & is1_muscl_d, is2_muscl_d, is3_muscl_d)

        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(inout) :: vL_rs_vf_x, vL_rs_vf_y, &
             & vL_rs_vf_z, vR_rs_vf_x, vR_rs_vf_y, vR_rs_vf_z
        integer, intent(in)               :: muscl_dir
        type(int_bounds_info), intent(in) :: is1_muscl_d, is2_muscl_d, is3_muscl_d
        integer                           :: j, k, l
        real(wp)                          :: aCL, aCR, aC, aTHINC, qmin, qmax, A, B, C, sign, moncon

        #:for MUSCL_DIR, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]
            if (muscl_dir == ${MUSCL_DIR}$) then
                $:GPU_PARALLEL_LOOP(collapse=3,private='[j, k, l, aCL, aC, aCR, aTHINC, moncon, sign, qmin, qmax]')
                do l = is3_muscl%beg, is3_muscl%end
                    do k = is2_muscl%beg, is2_muscl%end
                        do j = is1_muscl%beg, is1_muscl%end
                            aCL = v_rs_ws_${XYZ}$_muscl(j - 1, k, l, eqn_idx%adv%beg)
                            aC = v_rs_ws_${XYZ}$_muscl(j, k, l, eqn_idx%adv%beg)
                            aCR = v_rs_ws_${XYZ}$_muscl(j + 1, k, l, eqn_idx%adv%beg)

                            moncon = (aCR - aC)*(aC - aCL)

                            if (aC >= ic_eps .and. aC <= 1._wp - ic_eps .and. moncon > moncon_cutoff) then  ! Interface cell

                                if (aCR - aCL > 0._wp) then
                                    sign = 1._wp
                                else
                                    sign = -1._wp
                                end if

                                qmin = min(aCR, aCL)
                                qmax = max(aCR, aCL) - qmin

                                C = (aC - qmin + sgm_eps)/(qmax + sgm_eps)
                                B = exp(sign*ic_beta*(2._wp*C - 1._wp))
                                A = (B/cosh(ic_beta) - 1._wp)/tanh(ic_beta)

                                ! Left reconstruction
                                aTHINC = qmin + 5e-1_wp*qmax*(1._wp + sign*A)
                                if (aTHINC < ic_eps) aTHINC = ic_eps
                                if (aTHINC > 1 - ic_eps) aTHINC = 1 - ic_eps
                                vL_rs_vf_${XYZ}$ (j, k, l, eqn_idx%cont%beg) = vL_rs_vf_${XYZ}$ (j, k, l, &
                                                  & eqn_idx%cont%beg)/vL_rs_vf_${XYZ}$ (j, k, l, eqn_idx%adv%beg)*aTHINC
                                vL_rs_vf_${XYZ}$ (j, k, l, eqn_idx%cont%end) = vL_rs_vf_${XYZ}$ (j, k, l, &
                                                  & eqn_idx%cont%end)/(1._wp - vL_rs_vf_${XYZ}$ (j, k, l, &
                                                  & eqn_idx%adv%beg))*(1._wp - aTHINC)
                                vL_rs_vf_${XYZ}$ (j, k, l, eqn_idx%adv%beg) = aTHINC
                                vL_rs_vf_${XYZ}$ (j, k, l, eqn_idx%adv%end) = 1 - aTHINC

                                ! Right reconstruction
                                aTHINC = qmin + 5e-1_wp*qmax*(1._wp + sign*(tanh(ic_beta) + A)/(1._wp + A*tanh(ic_beta)))
                                if (aTHINC < ic_eps) aTHINC = ic_eps
                                if (aTHINC > 1 - ic_eps) aTHINC = 1 - ic_eps
                                vR_rs_vf_${XYZ}$ (j, k, l, eqn_idx%cont%beg) = vL_rs_vf_${XYZ}$ (j, k, l, &
                                                  & eqn_idx%cont%beg)/vL_rs_vf_${XYZ}$ (j, k, l, eqn_idx%adv%beg)*aTHINC
                                vR_rs_vf_${XYZ}$ (j, k, l, eqn_idx%cont%end) = vL_rs_vf_${XYZ}$ (j, k, l, &
                                                  & eqn_idx%cont%end)/(1._wp - vL_rs_vf_${XYZ}$ (j, k, l, &
                                                  & eqn_idx%adv%beg))*(1._wp - aTHINC)
                                vR_rs_vf_${XYZ}$ (j, k, l, eqn_idx%adv%beg) = aTHINC
                                vR_rs_vf_${XYZ}$ (j, k, l, eqn_idx%adv%end) = 1 - aTHINC
                            end if
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endfor

    end subroutine s_interface_compression

    !> Finalize the MUSCL module
    subroutine s_finalize_muscl_module()

        @:DEALLOCATE(v_rs_muscl)
 
    end subroutine s_finalize_muscl_module

end module m_muscl
