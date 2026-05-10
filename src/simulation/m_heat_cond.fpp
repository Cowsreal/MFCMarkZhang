#:include 'case.fpp'
#:include 'macros.fpp'
#:include 'inline_riemann.fpp'

module m_heat_conduction

    use m_derived_types
    use m_global_parameters
    use m_ibm

    implicit none
    private; public :: s_compute_heat_conduction_flux
contains

    subroutine s_compute_heat_conduction_flux(q_prim_vf, flux_src_vf, norm_dir, isx, isy, isz)

        type(int_bounds_info), intent(in):: isx, isy, isz
        type(scalar_field), dimension(sys_size), intent(inout) :: flux_src_vf
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        integer, intent(in) :: norm_dir

        integer :: j, k, l
        real(wp), dimension(0:5) :: T
        real(wp) :: dT_dx, dxp, dyp

        #:for NORM_DIR, JR, KR, LR, DX in &
            [(1, 'j+1', 'k',   'l',   'x_cc(j+1) - x_cc(j)'), &
            (2, 'j',   'k+1', 'l',   'y_cc(k+1) - y_cc(k)'), &
            (3, 'j',   'k',   'l+1', 'z_cc(l+1) - z_cc(l)')]
            #:if not MFC_CASE_OPTIMIZATION or num_dims >= NORM_DIR
                if (norm_dir == ${NORM_DIR}$) then
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[dT_dx,j,k,l]')
                    do l = isz%beg, isz%end
                        do k = isy%beg, isy%end
                            do j = isx%beg, isx%end
                            if(ib_markers%sf(j, k, l) == 0 .and. ib_markers%sf(${JR}$,${KR}$,${LR}$) == 0) then

                                    dT_dx = gammas(1) * &
                                        (q_prim_vf(eqn_idx%E)%sf(${JR}$, ${KR}$, ${LR}$) &
                                            / (q_prim_vf(eqn_idx%cont%beg)%sf(${JR}$, ${KR}$, ${LR}$) * cvs(1)) - &
                                        q_prim_vf(eqn_idx%E)%sf(j, k, l) &
                                            / (q_prim_vf(eqn_idx%cont%beg)%sf(j, k, l) * cvs(1))) &
                                        / (${DX}$)
                                    flux_src_vf(eqn_idx%E)%sf(j, k, l) = &
                                        flux_src_vf(eqn_idx%E)%sf(j, k, l) - kaps(1) * dT_dx
                            end if
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end if
            #:endif
        #:endfor
        ! if (norm_dir == 1) then
        !     $:GPU_PARALLEL_LOOP(collapse=3, private='[dxp,dyp,T,j,k,l]')
        !     do l = isz%beg, isz%end
        !         do k = isy%beg, isy%end
        !             do j = isx%beg, isx%end
        !                 dxp = x_cc(j + 1) - x_cc(j)
        !                 dyp = dy(k)
        !                 T(0) = q_prim_vf(eqn_idx%E)%sf(j, k, l) &
        !                     / (q_prim_vf(eqn_idx%cont%beg)%sf(j, k, l) * cvs(1))
        !                 T(1) = q_prim_vf(eqn_idx%E)%sf(j + 1, k, l) &
        !                     / (q_prim_vf(eqn_idx%cont%beg)%sf(j + 1, k, l) * cvs(1))
        !                 T(2) = q_prim_vf(eqn_idx%E)%sf(j, k + 1, l) &
        !                     / (q_prim_vf(eqn_idx%cont%beg)%sf(j, k + 1, l) * cvs(1))
        !                 T(3) = q_prim_vf(eqn_idx%E)%sf(j + 1, k + 1, l) &
        !                     / (q_prim_vf(eqn_idx%cont%beg)%sf(j + 1, k + 1, l) * cvs(1))
        !                 T(4) = q_prim_vf(eqn_idx%E)%sf(j, k - 1, l) &
        !                     / (q_prim_vf(eqn_idx%cont%beg)%sf(j, k - 1, l) * cvs(1))
        !                 T(5) = q_prim_vf(eqn_idx%E)%sf(j + 1, k - 1, l) &
        !                     / (q_prim_vf(eqn_idx%cont%beg)%sf(j + 1, k - 1, l) * cvs(1))
        !
        !                 flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, l) - kaps(1) * gammas(1)&
        !                     ((T(0) + T(1)) / dxp - (0.25_wp * (T(0) + T(1) + T(2) + T(3)) - 0.25_wp * (T(0) + T(1) + T(4) + T(5)) /
        !                 dyp)
        !             end do
        !         end do
        !     end do
        !     $:END_GPU_PARALLEL_LOOP()
        ! end if
        !
        ! if (norm_dir == 2) then
        !     $:GPU_PARALLEL_LOOP(collapse=3, private='[dxp,dyp,T,j,k,l]')
        !     do l = isz%beg, isz%end
        !         do k = isy%beg, isy%end
        !             do j = isx%beg, isx%end
        !                 dxp = x_cc(j + 1) - x_cc(j)
        !                 dyp = dy(k)
        !                 T(0) = q_prim_vf(eqn_idx%E)%sf(j, k, l) &
        !                     / (q_prim_vf(eqn_idx%cont%beg)%sf(j, k, l) * cvs(1))
        !                 T(1) = q_prim_vf(eqn_idx%E)%sf(j, k + 1, l) &
        !                     / (q_prim_vf(eqn_idx%cont%beg)%sf(j, k + 1, l) * cvs(1))
        !                 T(2) = q_prim_vf(eqn_idx%E)%sf(j + 1, k, l) &
        !                     / (q_prim_vf(eqn_idx%cont%beg)%sf(j + 1, k, l) * cvs(1))
        !                 T(3) = q_prim_vf(eqn_idx%E)%sf(j + 1, k + 1, l) &
        !                     / (q_prim_vf(eqn_idx%cont%beg)%sf(j + 1, k + 1, l) * cvs(1))
        !                 T(4) = q_prim_vf(eqn_idx%E)%sf(j - 1, k, l) &
        !                     / (q_prim_vf(eqn_idx%cont%beg)%sf(j - 1, k, l) * cvs(1))
        !                 T(5) = q_prim_vf(eqn_idx%E)%sf(j - 1, k + 1, l) &
        !                     / (q_prim_vf(eqn_idx%cont%beg)%sf(j - 1, k + 1, l) * cvs(1))
        !
        !                 flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, l) - kaps(1) * gammas(1)&
        !                     ((T(0) + T(1)) / dyp - (0.25_wp * (T(0) + T(1) + T(2) + T(3)) - 0.25_wp * (T(0) + T(1) + T(4) + T(5)) /
        !                 dxp)
        !             end do
        !         end do
        !     end do
        !     $:END_GPU_PARALLEL_LOOP()
        ! end if


    end subroutine s_compute_heat_conduction_flux

end module m_heat_conduction
