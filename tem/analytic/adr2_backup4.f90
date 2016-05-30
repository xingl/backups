program res_itg

        use file_utils

        implicit none

        real (kind=8), parameter :: pi=3.141592653589793
        complex (kind=8), parameter :: zi=(0.0,1.0)

        integer (kind=8) :: nkz, nky, nkx
        real (kind=8) :: dkz, dky, dkx
        integer (kind=8) :: ikz, iky, ikx
        real (kind=8) :: ky0, kz0, kx0
        real (kind=8) :: kx1, ky1, kz1
        real (kind=8), dimension(:), allocatable :: ky_grid, kz_grid, kx_grid

        integer (kind=8) :: nfprime_i, ntprime_i, nfprime_e, ntprime_e
        integer (kind=8) :: nomd
        real (kind=8) :: dfprime_i, dtprime_i, dfprime_e, dtprime_e
        real (kind=8) :: domd
        integer (kind=8) :: ifprime_i, itprime_i, ifprime_e, itprime_e
        integer (kind=8) ::  iomd
        real (kind=8) :: tprime0_i, fprime0_i, tprime0_e, fprime0_e
        real (kind=8) :: omd0
        real (kind=8), dimension(:), allocatable :: fprime_i_grid, tprime_i_grid
        real (kind=8), dimension(:), allocatable :: fprime_e_grid, tprime_e_grid
        real (kind=8), dimension(:), allocatable :: omd_grid
        real (kind=8) :: fprime_i, fprime_e, tprime_i, tprime_e

        real (kind=8) :: vy_1, vy_2, vz_1, vz_2
        real (kind=8) :: eps, nu0, A, fprime, tprime, omd, ky, kz, kx, kv
        real (kind=8) :: vi, int_tol,sec_tol, fr, na_e, na_z, na_i, omd_kx, t_e
        real (kind=8) :: theta, Zeff, Z, mu_e, mu_z
        complex (kind=8) :: seed1, seed2
        real (kind=8) :: om1_re, om1_im
        integer (kind=8) :: sec_nsteps, sec_istep
        integer (kind=8) :: int_nsteps, int_istep

        integer(kind=8), dimension(:),allocatable :: ikx_ref, ikz_ref, iky_ref

        complex (kind=8) :: omega 
        complex (kind=8) :: root, root_ikx_ref, root_ikzx_ref, root_ikyzx_ref
        complex (kind=8) :: root_ikx_1, root_ikzx_1, root_ikyzx_1
        integer :: gamma_unit=101, Dmixing_unit=102
        integer :: out_unit=103
        real(kind=8) :: cut_buffer, seed_fr, lower_bnd
        
        call init_file_utils
        call read_input_file
        call init_grids
        call init_output_files

        do iomd = 1, nomd
           omd = omd_grid(iomd)
           print*, 'omd'
           print*, omd
           do ifprime_i = 1, nfprime_i
              fprime_i = fprime_i_grid(ifprime_i)
              print*,'omn'
              print*, fprime_i
              do itprime_i = 1, ntprime_i
                 tprime_i = tprime_i_grid(itprime_i)
                 print*,'omt_i'
                 print*, tprime_i
                 do ifprime_e = 1, nfprime_e
!                    fprime_e = fprime_e_grid(ifprime_e)
                    fprime_e = fprime_i
                    do itprime_e = 1, ntprime_e
                       tprime_e = tprime_e_grid(itprime_e)
                       print*,'omt_e'
                       print*,tprime_e
                       call ky_scan
                    end do
                 end do
              end do
           end do
        end do

        call close_output_file (gamma_unit)
        call close_output_file (Dmixing_unit)
        call close_output_file (out_unit)
        call finish_file_utils

contains

 subroutine ky_scan 

        implicit none

        complex (kind=8), dimension(:), allocatable :: omega_ky_scan
        real (kind=8) :: kz_gam,kz_Dmix,Dmixing_kz_max
        complex (kind=8) :: omega_kz_max

        allocate(omega_ky_scan(nky))
        omega_ky_scan = 0. -9999.*zi
        iky_ref = minloc(abs(ky_grid-ky1))
        do iky = 1, nky
           ky = ky_grid(iky)
           call kz_scan(kz_gam,omega_kz_max,kz_Dmix,Dmixing_kz_max)
           omega_ky_scan(iky) = omega_kz_max
           if (iky==iky_ref(1).and.aimag(omega_kz_max)>aimag(root_ikyzx_ref)) &
              root_ikyzx_ref = omega_kz_max
        end do

        do iky = iky_ref(1), nky
           ky = ky_grid(iky)
           call kz_scan(kz_gam,omega_kz_max,kz_Dmix,Dmixing_kz_max)
           if (aimag(omega_kz_max)>aimag(omega_ky_scan(iky))) &
             omega_ky_scan(iky) = omega_kz_max
        end do

        do iky = iky_ref(1),1,-1
           ky = ky_grid(iky)
           call kz_scan(kz_gam,omega_kz_max,kz_Dmix,Dmixing_kz_max)
           if (aimag(omega_kz_max)>aimag(omega_ky_scan(iky))) &
             omega_ky_scan(iky) = omega_kz_max
        end do

        root_ikyzx_ref = omega_ky_scan(iky_ref(1))
        root_ikyzx_1 = omega_ky_scan(1)
 end subroutine

 subroutine kz_scan (kz_maxgam,omega_maxgam,kz_maxDmix,Dmixing_max)

        implicit none

        real (kind=8),intent(out) :: kz_maxgam, kz_maxDmix, Dmixing_max
        complex (kind=8),intent(out) :: omega_maxgam
        integer (kind=8), dimension(:), allocatable :: ikz_maxgamma
        integer (kind=8), dimension(:), allocatable :: ikz_maxDmixing
        real (kind=8), dimension(:), allocatable :: Dmixing_kz_scan
        complex (kind=8), dimension(:), allocatable :: omega_kz_scan
        real (kind=8) :: kx_gam,kx_Dmix,Dmixing_kx_max
        complex (kind=8) :: omega_kx_max

        allocate(Dmixing_kz_scan(nkz),omega_kz_scan(nkz))
        allocate(ikz_maxgamma(1),ikz_maxDmixing(1))
        Dmixing_kz_scan = -9999.0
        omega_kz_scan = 0. -9999.*zi
        ikz_ref = minloc(abs(kz_grid-kz1))
        do ikz = 1, nkz
           kz = kz_grid(ikz)
           print*, "ky, kz"
           print*, ky, kz
           call kx_scan(kx_gam,omega_kx_max,kx_Dmix,Dmixing_kx_max)
           omega_kz_scan(ikz) = omega_kx_max
           Dmixing_kz_scan(ikz) = Dmixing_kx_max
           if (ikz==ikz_ref(1).and.aimag(omega_kx_max)>aimag(root_ikzx_ref)) &
              root_ikzx_ref = omega_kx_max
        end do

        do ikz = ikz_ref(1), nkz
           kz = kz_grid(ikz)
           print*, "ky, kz"
           print*, ky, kz
           call kx_scan(kx_gam,omega_kx_max,kx_Dmix,Dmixing_kx_max)
           if (aimag(omega_kx_max)>aimag(omega_kz_scan(ikz))) &
              omega_kz_scan(ikz) = omega_kx_max
           if (Dmixing_kx_max>Dmixing_kz_scan(ikz)) &
              Dmixing_kz_scan(ikz) = Dmixing_kx_max
          if (.not.aimag(omega_kz_scan(ikz))==-9999.0) &
            write (gamma_unit,&
            '(10e12.4)') tprime_i,fprime_i,tprime_e,fprime_e,omd, &
            kx_gam,ky,kz,omega_kz_scan(ikz)
!          if (.not.Dmixing_kz_scan(ikz)==-9999.0) &
!            write (Dmixing_unit, &
!            '(9e12.4)') tprime_i,fprime_i,tprime_e,fprime_e,omd, &
!            kx_Dmix,ky,kz,Dmixing_kz_scan(ikz)
        end do
        write (gamma_unit,*)
!        write (Dmixing_unit,*)

        do ikz = ikz_ref(1),1,-1
           kz = kz_grid(ikz)
           print*, "ky, kz"
           print*, ky, kz
           call kx_scan(kx_gam,omega_kx_max,kx_Dmix,Dmixing_kx_max)
           if (aimag(omega_kx_max)>aimag(omega_kz_scan(ikz))) &
              omega_kz_scan(ikz) = omega_kx_max
           if (Dmixing_kx_max>Dmixing_kz_scan(ikz)) &
              Dmixing_kz_scan(ikz) = Dmixing_kx_max
          if (.not.aimag(omega_kz_scan(ikz))==-9999.0) &
            write (gamma_unit,&
            '(10e12.4)') tprime_i,fprime_i,tprime_e,fprime_e,omd, &
            kx_gam,ky,kz,omega_kz_scan(ikz)
!          if (.not.Dmixing_kz_scan(ikz)==-9999.0) &
!            write (Dmixing_unit, &
!            '(9e12.4)') tprime_i,fprime_i,tprime_e,fprime_e,omd, &
!            kx_Dmix,ky,kz,Dmixing_kz_scan(ikz)
        end do
        write (gamma_unit,*)
!        write (Dmixing_unit,*)
        
        ikz_maxgamma = maxloc(aimag(omega_kz_scan))
        ikz_maxDmixing = maxloc(Dmixing_kz_scan)
        kz_maxgam = kz_grid(ikz_maxgamma(1))
        kz_maxDmix = kz_grid(ikz_maxDmixing(1))
        omega_maxgam = omega_kz_scan(ikz_maxgamma(1))
        Dmixing_max = Dmixing_kz_scan(ikz_maxDmixing(1))

        root_ikzx_ref = omega_kz_scan(ikz_ref(1))
        root_ikzx_1 = omega_kz_scan(1)

 end subroutine

 subroutine kx_scan (kx_maxgam,omega_maxgam,kx_maxDmix,Dmixing_max)

        implicit none

        real (kind=8), intent(out) :: kx_maxgam, kx_maxDmix, Dmixing_max
        complex (kind=8), intent(out) :: omega_maxgam
        integer (kind=8), dimension(:), allocatable :: ikx_maxgamma
        integer (kind=8), dimension(:), allocatable :: ikx_maxDmixing
        real (kind=8), dimension(:), allocatable :: Dmixing_kx_scan
        complex (kind=8), dimension(:), allocatable :: omega_kx_scan
        complex (kind=8) :: root_tmp,om_linr,om_linr_a,om_quadr,om_quadr_a
        real (kind=8) :: Dmixing_tmp


          allocate(Dmixing_kx_scan(nkx), omega_kx_scan(nkx))
          allocate(ikx_maxgamma(1),ikx_maxDmixing(1))
          Dmixing_kx_scan = -9999.0
          omega_kx_scan = 0. -9999.*zi
          ikx_ref = minloc(abs(kx_grid-kx1))
          do ikx = 1, nkx
             kx = kx_grid(ikx)
             kv = ((1-omd_kx)*ky + omd_kx*kx)*omd
             print*, "kx,kv"
             print*, kx,kv
             call init_roots
             print*,seed1,seed2
             call rootfinder(seed1,seed2,root_tmp)
             if (sec_istep==sec_nsteps.or.isnan(real(root_tmp))) then
              print*, 'No root is found.'
              root = seed1
              omega_kx_scan(ikx) = 0.-9999.*zi
             else if (abs(aimag(root_tmp))<lower_bnd) then
              print*, 'No root is found.'
              root = seed1
              omega_kx_scan(ikx) = 0. -9999.*zi
             else if (aimag(root_tmp)==aimag(root_tmp)-1..or.&
                      isnan(aimag(root_tmp))) then
              print*, 'No root is found.'
              root = seed1
              omega_kx_scan(ikx) = 0.-9999.*zi
             else
              print*, "Root is found:"
              print*, root_tmp    
              omega_kx_scan(ikx) = root_tmp
              root = root_tmp
              !write (out_unit, '(8e12.4)') tprime,fprime,&
              !               omd,kx,ky,kz,omega_kx_scan(ikx)
             end if
             if (ikx==ikx_ref(1).and.aimag(omega_kx_scan(ikx))>aimag(root_ikx_ref)) &
              root_ikx_ref = omega_kx_scan(ikx)
          end do
          do ikx = ikx_ref(1), 1, -1
             kx = kx_grid(ikx)
             kv = ((1-omd_kx)*ky + omd_kx*kx)*omd
             print*, "kx,kv"
             print*, kx,kv
             call init_roots
             print*,seed1,seed2
             call rootfinder(seed1,seed2,root_tmp)
             if (sec_istep==sec_nsteps.or.isnan(real(root_tmp))) then
              print*, 'No root is found.'
              root = seed1
              omega_kx_scan(ikx) = 0.-9999.*zi
             else if (abs(aimag(root_tmp))<lower_bnd) then
              print*, 'No root is found.'
              root = seed1
              omega_kx_scan(ikx) = 0. -9999.*zi
             else if (aimag(root_tmp)==aimag(root_tmp)-1..or.&
                      isnan(aimag(root_tmp))) then
              print*, 'No root is found.'
              root = seed1
              omega_kx_scan(ikx) = 0.-9999.*zi
             else
              print*, "Root is found:"
              print*, root_tmp    
              root = root_tmp
              if (aimag(root_tmp)>aimag(omega_kx_scan(ikx))) &
                omega_kx_scan(ikx) = root_tmp
              call quadratic_solution(om_quadr,om_quadr_a)
              call linear_solution(om_linr,om_linr_a)
              write (out_unit, '(14e12.4)') tprime_i,fprime_i,tprime_e,&
                             ky,omega_kx_scan(ikx),om_linr,om_linr_a,&
                             om_quadr,om_quadr_a
             end if
          end do
          write (out_unit,*) 
          do ikx = ikx_ref(1), nkx
             kx = kx_grid(ikx)
             kv = ((1-omd_kx)*ky + omd_kx*kx)*omd
             print*, "kx,kv"
             print*, kx,kv
             call init_roots
             print*,seed1,seed2
             call rootfinder(seed1,seed2,root_tmp)
             if (sec_istep==sec_nsteps.or.isnan(real(root_tmp))) then
              print*, 'No root is found.'
              root = seed1
              omega_kx_scan(ikx) = 0.-9999.*zi
             else if (abs(aimag(root_tmp))<lower_bnd) then
              print*, 'No root is found.'
              root = seed1
              omega_kx_scan(ikx) = 0.-9999.*zi
             else if (aimag(root_tmp)==aimag(root_tmp)-1..or.&
                      isnan(aimag(root_tmp))) then
              print*, 'No root is found.'
              root = seed1
              omega_kx_scan(ikx) = 0.-9999.*zi
             else
              print*, "Root is found:"
              print*, root_tmp    
              root = root_tmp
              if (aimag(root_tmp)>aimag(omega_kx_scan(ikx))) &
                omega_kx_scan(ikx) = root_tmp
              call quadratic_solution(om_quadr,om_quadr_a)
              call linear_solution(om_linr,om_linr_a)
              write (out_unit, '(14e12.4)') tprime_i,fprime_i,tprime_e,&
                             ky,omega_kx_scan(ikx),om_linr,om_linr_a,&
                             om_quadr,om_quadr_a
             end if
          end do
          write (out_unit, *)

          ikx_maxgamma = maxloc(aimag(omega_kx_scan))
          Dmixing_kx_scan = aimag(omega_kx_scan)/(kx_grid**2+ky**2)
          ikx_maxDmixing = maxloc(Dmixing_kx_scan)

          kx_maxgam = kx_grid(ikx_maxgamma(1))
          omega_maxgam = omega_kx_scan(ikx_maxgamma(1))
          kx_maxDmix = kx_grid(ikx_maxDmixing(1))
          Dmixing_max = Dmixing_kx_scan(ikx_maxDmixing(1))

          if (.not.aimag(omega_maxgam)==-9999.0) then  
              root_ikx_ref = omega_kx_scan(ikx_ref(1))
              root_ikx_1 = omega_kx_scan(1)
          end if
          deallocate(Dmixing_kx_scan,omega_kx_scan)
          deallocate(ikx_maxgamma,ikx_maxDmixing)

  end subroutine

subroutine quadratic_solution(omega_quadr, omega_quadr_a)

        implicit none
        
        complex(kind=8), intent(out) :: omega_quadr, omega_quadr_a
        complex(kind=8) :: omega1, omega2
        complex (kind=8) :: term1, term2, term3
        integer(kind=8) :: t,s
        
        s = 2   !fluid ion limit
        t = 1 
        call vy_integral_simpson(omega,vy_1,vy_2,term1,s,t)
        t = 2
        call vy_integral_simpson(omega,vy_1,vy_2,term2,s,t)
        t = 3 
        call vy_integral_simpson(omega,vy_1,vy_2,term3,s,t)

        term1 = -1./theta - Zeff + term1

        omega1 = (-term2+sqrt(term2**2-4.*term1*term3))/2./term1
        omega2 = (-term2-sqrt(term2**2-4.*term1*term3))/2./term1

        if (aimag(omega1)>aimag(omega2)) then
           omega_quadr = omega1
        else
           omega_quadr = omega2
        end if

        omega1 = -(fprime_i*ky-2.*omd*ky/a)/2. + &
                sqrt((fprime_i*ky-2.*omd*ky/a)**2 - &
                8*(fprime_i*ky+tprime_i*ky)*omd*ky/a)/2.
        omega1 = -(fprime_i*ky-2.*omd*ky/a)/2. - &
                sqrt((fprime_i*ky-2.*omd*ky/a)**2 - &
                8*(fprime_i*ky+tprime_i*ky)*omd*ky/a)/2.

        if (aimag(omega1)>aimag(omega2)) then
           omega_quadr_a = omega1
        else
           omega_quadr_a = omega2
        end if

end subroutine

subroutine linear_solution(omega_linr, omega_linr_a)

        implicit none
        
        complex(kind=8), intent(out) :: omega_linr, omega_linr_a
        complex(kind=8) :: omega1, omega2
        complex (kind=8) :: term1, term2, term3
        integer(kind=8) :: t,s
        
        s = 2   !fluid ion limit, omd=0.
        t = 1 
        call vy_integral_simpson(omega,vy_1,vy_2,term1,s,t)
        t = 2
        call vy_integral_simpson(omega,vy_1,vy_2,term2,s,t)

        term1 = 1./theta + Zeff - term1

        omega1 = term2/term1

        omega_linr = omega1

        omega_linr_a = -(fprime_i*ky*(1-ky**2+ky**4/2.)+tprime_i*ky*(-ky**2+ky**4))/(1+ky**2-ky**4/2.)

end subroutine

 subroutine init_roots
 
        implicit none

        complex(kind=8) :: omega_quadr

        !if (kz==0..and.ky<0.8.and.(na_e+na_i==1.)) then
        if (0.==1.) then
           !call quadratic_solution(omega_quadr)
           !seed1 = omega_quadr
        else if (iomd==1.and.iky==1.and.ikz==1.and.ikx==1) then
           seed1 = om1_re + om1_im*zi
        else if (iomd==1.and.iky==1.and.ikz==1.and.ikx==1.and.aimag(root_ikyzx_1)>aimag(root)) then
           seed1 = root_ikyzx_1
        else if (iomd==1.and.ikz==1.and.ikx==1.and.aimag(root_ikzx_1)>aimag(root)) then
           seed1 = root_ikzx_1
        else if (iomd==1.and.ikx==1.and.aimag(root_ikx_1)>aimag(root)) then
           seed1 = root_ikx_1
        else if (iomd==1.and.iky==iky_ref(1).and.ikz==ikz_ref(1).and.ikx==ikx_ref(1).and.aimag(root_ikyzx_ref)>aimag(root)) then
           seed1 = root_ikyzx_ref
        else if (iomd==1.and.ikz==ikz_ref(1).and.ikx==ikx_ref(1).and.aimag(root_ikzx_ref)>aimag(root)) then
           seed1 = root_ikzx_ref
        else if (iomd==1.and.ikx==ikx_ref(1).and.aimag(root_ikx_ref)>aimag(root)) then
           seed1 = root_ikx_ref
        else 
           seed1 = root
        end if

        if (aimag(seed1)>0.) then
           seed2 = seed1*(1.0-seed_fr)
        else
           seed2 = seed1*(1.0+seed_fr)
        end if

 end subroutine

 subroutine rootfinder(sd1,sd2,rt)

        implicit none

        complex (kind=8), intent(in) :: sd1,sd2
        complex (kind=8), intent(out) :: rt
        complex (kind=8) :: x_2,f_2,x_1,f_1
        complex (kind=8) :: x_tmp, f_tmp

        x_1 = sd1
        call dispersion_relation(x_1,f_1)
        x_2 = sd2
        call dispersion_relation(x_2,f_2)
        x_tmp = x_1
        sec_istep = 0
        do while (abs(f_1)>sec_tol .and. abs(f_2)>sec_tol .and. sec_istep<sec_nsteps) 
            x_tmp = x_1 - f_1 * ((x_1-x_2)/(f_1-f_2))
                call dispersion_relation(x_tmp,f_tmp)
                if (isnan(aimag(f_tmp)).or.isnan(real(f_tmp))) then
                   sec_istep = sec_nsteps
                   exit
                end if
                print*, 'root finder'
                print*, sec_istep, x_tmp, f_tmp
                f_2 = f_1
                f_1 = f_tmp
                x_2 = x_1
                x_1 = x_tmp
                sec_istep = sec_istep +1
        end do
        rt = x_tmp

  end subroutine

subroutine dispersion_relation(omega, rhs)
        
        implicit none

        complex(kind=8), intent(in) :: omega
        complex (kind=8), intent(out) :: rhs
        complex (kind=8) :: integral_e, integral_i, integral_z
        complex (kind=8) :: integral_te, integral_te_s_nu_a, integral_te_l_nu_a
        complex (kind=8) :: integral_te_l_nu, integral_te_s_nu
        integer(kind=8) :: s

        if (t_e==0.) then
            integral_te = 0.
        else
            s = 4
            call vy_integral(omega,integral_te,s)
            s = 5
            call vy_integral(omega,integral_te_l_nu,s)
            s = 6
            call vy_integral(omega,integral_te_s_nu,s)

            integral_te_s_nu_a = sqrt(2.*eps)*(1.-(-theta*fprime_i*ky)/omega)
            integral_te_l_nu_a = sqrt(2.*eps)*(-zi)*eps/3./nu0*&
                (8.*sqrt(2./pi)*(omega-(-theta*fprime_i*ky))-12.*sqrt(2./pi)&
                *(-theta*tprime_e*ky))
            
            write (Dmixing_unit, '(13e12.4)') nu0,omega, &
                integral_te,integral_te_s_nu,integral_te_s_nu_a,&
                integral_te_l_nu, integral_te_l_nu_a
        end if
        if (na_e==0.) then
            integral_e = 0.
        else
            s = 1
            call vy_integral(omega,integral_e,s)
            print*, 'integral_e'
            print*, integral_e
        end if
        if (na_i==0.) then
            integral_i = 0.
        else
            s = 2
            call vy_integral(omega,integral_i,s)
            print*, 'integral_i'
            print*, integral_i
        end if
        if (na_z==0.) then
            integral_z = 0.
        else
            s = 3
            call vy_integral(omega,integral_z,s)
        end if

        rhs = 1.0 + theta*Zeff - na_e*integral_e - &
              na_i*theta*(Z-Zeff)/(Z-1.0)*integral_i - &
              na_z*theta*(Zeff-1.0)*Z/(Z-1.0)*integral_z - &
              t_e*integral_te

end subroutine

subroutine vy_integral(omega,integral,s)

        implicit none

        integer(kind=8), intent(in) :: s
        complex(kind=8), intent(in) :: omega
        complex (kind=8), intent(out) :: integral
        complex (kind=8) :: integral1, integral2, integral3, integral4
        real (kind=8) :: discrim1, discrim2
        real (kind=8) :: vy_d, vy_l, vy_r, vy_b
        
        vy_b = cut_buffer
        integral1 = 0.
        integral2 = 0.
        integral3 = 0.
        integral4 = 0.

        if (.not.kv==0.) then
          if (s==2) then
           discrim1 = real(kz**2-4.*kv/A*(kv/2./A*vy_1**2-omega))
           discrim2 = real(kz**2-4.*kv/A*(kv/2./A*vy_2**2-omega))
           vy_d = sqrt(2.*A/kv*(real(omega)+kz**2*A/kv/4.))
          else if (s==1) then
           discrim1 = real((sqrt(theta/mu_e)*kz)**2-&
             (-theta)*4.*kv/A*(-theta*kv/2./A*vy_1**2-omega))
           discrim2 = real((sqrt(theta/mu_e)*kz)**2-&
             (-theta)*4.*kv/A*(-theta*kv/2./A*vy_2**2-omega))
           vy_d = sqrt(2.*A/kv/(-theta)*(real(omega)+&
             (sqrt(theta/mu_e)*kz)**2*A/kv/(-theta)/4.))
          else if (s==3) then
           discrim1 = real((sqrt(1./mu_z)*kz)**2-&
             (1./Z)*4.*kv/A*(-1./Z*kv/2./A*vy_1**2-omega))
           discrim2 = real((sqrt(1./mu_z)*kz)**2-&
             (1./Z)*4.*kv/A*(-1./Z*kv/2./A*vy_2**2-omega))
           vy_d = sqrt(2.*A/kv*Z*(real(omega)+&
             (sqrt(1./mu_z)*kz)**2*A/kv*Z/4.))
          else if (s==4.or.s==5.or.s==6) then
           discrim1 = 1.
           discrim2 = 1.
           vy_d = 0.
          end if
          if ((.not.aimag(omega)>0.).and.discrim1*discrim2<0.) then
           vy_l = vy_d - vy_b
           vy_r = vy_d + vy_b
           if (vy_l>vy_1 .and. vy_r<vy_2) then
              call vy_integral_simpson(omega,vy_1,vy_l,integral1,s,t=0)
              call vy_integral_midpoint(omega,vy_l,vy_d,integral2,s)
              call vy_integral_midpoint(omega,vy_d,vy_r,integral3,s)
              call vy_integral_simpson(omega,vy_r,vy_2,integral4,s,t=0)
              integral = integral1+integral2+integral3+integral4 
           else if (vy_l<vy_1 .and. vy_r<vy_2) then
              integral1 = 0.
              call vy_integral_midpoint(omega,vy_1,vy_d,integral2,s)
              call vy_integral_midpoint(omega,vy_d,vy_r,integral3,s)
              call vy_integral_simpson(omega,vy_r,vy_2,integral4,s,t=0)
              integral = integral1+integral2+integral3+integral4 
           else if (vy_l>vy_1 .and. vy_r>vy_2) then
              call vy_integral_simpson(omega,vy_1,vy_l,integral1,s,t=0)
              call vy_integral_midpoint(omega,vy_l,vy_d,integral2,s)
              call vy_integral_midpoint(omega,vy_d,vy_2,integral3,s)
              integral4 = 0.
              integral = integral1+integral2+integral3+integral4 
           end if
          else !s==4
           call vy_integral_simpson(omega,vy_1,vy_2,integral1,s,t=0) 
           integral2 = 0.
           integral3 = 0.
           integral4 = 0.
           integral = integral1+integral2+integral3+integral4 
           !print*,'vy_integral, omega, s'
           !print*, integral, omega, s
          end if
        end if

        !print*, 'integral'
        !print*, omega,integral,s

end subroutine

subroutine vy_integral_midpoint(omega,vy_a,vy_b,integral,s)

        implicit none

        integer(kind=8), intent(in) :: s
        integer(kind=8) :: t
        real(kind=8),intent(in) :: vy_a, vy_b
        complex(kind=8), intent(in) :: omega
        complex (kind=8), intent(out) :: integral
        real(kind=8) :: vy_m, dvy, vy_s, vy_f
        complex(kind=8) :: f_m

        integral = 0.
        dvy = 0.02*cut_buffer
        vy_s = vy_a
        vy_f = vy_s+dvy
        do while (vy_f<vy_b)
                vy_m = (vy_s+vy_f)/2.
                call vz_integral(vy_m,f_m,s,omega,t=0)
                integral = integral + f_m*dvy
                vy_s = vy_f
                vy_f = vy_f + dvy
                if (vy_f>vy_b) then
                    vy_f = vy_b
                    dvy = vy_f-vy_s
                end if
        end do

        !print*, 'vy_integral_midpoint'
        !print*, omega,integral,s,t

end subroutine

subroutine vy_integral_simpson(omega,vy_a,vy_b,integral,s,t)

        implicit none

        integer(kind=8), intent(in) :: s,t
        real(kind=8),intent(in) :: vy_a, vy_b
        complex(kind=8), intent(in) :: omega
        complex (kind=8), intent(out) :: integral
        real(kind=8) :: x_left, x_center, x_right, dx
        complex (kind=8) :: f_left, f_center, f_right
        complex (kind=8) :: i_trapezoid, i_simpson
        integer :: vy_istep

        vy_istep = 0.
        dx = 0.1
        x_left = vy_a
        x_right = x_left + dx
        call vz_integral(x_left,f_left,s,omega,t)
        integral = 0.0

        do while (x_left<vy_b.and.vy_istep<int_nsteps)
                x_center = 0.5*(x_left+x_right)
                call vz_integral(x_center,f_center,s,omega,t)
                call vz_integral(x_right,f_right,s,omega,t)
                i_trapezoid = 0.5*(f_left+2.0*f_center+f_right)*0.5*dx
                i_simpson = 1.0/6.0*dx*(f_left+4.0*f_center+f_right)
                do while (abs(i_trapezoid-i_simpson)>int_tol)
                        dx=fr*dx*(abs(i_trapezoid-i_simpson)/int_tol)**(-1.0/3.0)
                        if (dx>0.2) dx=0.2
                        x_right = x_left + dx
                        if (x_right>vy_b) then
                                x_right = vy_b
                                dx = x_right - x_left
                        end if
                        x_center = 0.5*(x_left+x_right)
                        call vz_integral(x_center,f_center,s,omega,t)
                        call vz_integral(x_right,f_right,s,omega,t)
                        i_trapezoid = 0.5*(f_left+2.0*f_center+f_right)*0.5*dx
                        i_simpson = 1.0/6.0*dx*(f_left+4.0*f_center+f_right)
                end do
                integral = integral + i_simpson
                x_left = x_right
                dx = fr*dx*(abs(i_trapezoid-i_simpson)/int_tol)**(-1.0/3.0)
                if (dx>0.2) dx=0.2
                x_right = x_left + dx
                if (x_right>vy_b) then
                        x_right = vy_b
                        dx = x_right - x_left
                end if
                f_left = f_right
                vy_istep = vy_istep+1
        end do
        if (vy_istep==int_nsteps) then
           int_istep = int_nsteps
           integral = 0.
        end if

        !print*, 'vy_integral_simpson'
        !print*, omega,integral

  end subroutine

  subroutine vz_integral(vy,integral,s,omega,t)

        implicit none

        integer(kind=8), intent(in) :: s,t
        real(kind=8), intent(in) :: vy
        complex(kind=8), intent(in) :: omega
        real(kind=8) :: vz_a, vz_b
        complex (kind=8), intent(out) :: integral
        real(kind=8) :: x_left, x_center, x_right, dx
        complex (kind=8) :: f_left, f_center, f_right
        complex (kind=8) :: i_trapezoid, i_simpson
        complex (kind=8) :: res_v1, res_v2
        complex (kind=8) :: shift_vz
        real(kind=8) :: rho,bj,dj,fj
        integer(kind=8) :: n
        integer(kind=8) :: vz_istep

        if (t==0.) then
           if (s==2) then
              call vz_integral_residue(vy,omega,s,res_v1,res_v2,&
                shift_vz)
           else if (s==1) then
              call vz_integral_residue(vy,omega,s,res_v1,res_v2,&
                shift_vz)
           else if (s==3) then
              call vz_integral_residue(vy,omega,s,res_v1,res_v2,&
                shift_vz)
           else if (s==4.or.s==5.or.s==6) then
              res_v1 = 0.
              res_v2 = 0.
              shift_vz = 0.
           end if
        else
           res_v1 = 0.
           res_v2 = 0.
           shift_vz = 0.
        end if

        if ((.not.s==4).and.(.not.s==5).and.(.not.s==6)) then
           vz_istep = 0
           vz_a = vz_1
           vz_b = vz_2
           dx = 0.1
           x_left = vz_a
           x_right = x_left + dx
           f_left = integrand(x_left,vy,omega,shift_vz,s,t)
           integral = 0.0

           do while (x_left<vz_b.and.vz_istep<int_nsteps)
                x_center = 0.5*(x_left+x_right)
                f_center = integrand(x_center,vy,omega,shift_vz,s,t)
                f_right = integrand(x_right,vy,omega,shift_vz,s,t)
                i_trapezoid = 0.5*(f_left+2.0*f_center+f_right)*0.5*dx
                i_simpson = 1.0/6.0*dx*(f_left+4.0*f_center+f_right)
                do while (abs(i_trapezoid-i_simpson)>int_tol)
                        dx=fr*dx*(abs(i_trapezoid-i_simpson)/int_tol)**(-1.0/3.0)
                        if (dx>0.2) dx=0.2
                        x_right = x_left + dx
                        if (x_right>vz_b) then
                                x_right = vz_b
                                dx = x_right - x_left
                        end if
                        x_center = 0.5*(x_left+x_right)
                        f_center = integrand(x_center,vy,omega,shift_vz,s,t)
                        f_right = integrand(x_right,vy,omega,shift_vz,s,t)
                        i_trapezoid = 0.5*(f_left+2.0*f_center+f_right)*0.5*dx
                        i_simpson = 1.0/6.0*dx*(f_left+4.0*f_center+f_right)
                end do
                integral = integral + i_simpson
                x_left = x_right
                dx = fr*dx*(abs(i_trapezoid-i_simpson)/int_tol)**(-1.0/3.0)
                if (dx>0.2) dx=0.2
                x_right = x_left + dx
                if (x_right>vz_b) then 
                        x_right = vz_b
                        dx = x_right - x_left
                end if
                f_left = f_right
                vz_istep = vz_istep+1
           end do
        else
           integral = integrand(vy,vy,omega,shift_vz,s,t)
        end if

        if (vz_istep==int_nsteps) then
           int_istep = int_nsteps
           integral = 0.
        else
           n=0
           if (s==2) then
              rho = sqrt(ky**2+kx**2)*vy
              call bjndd (n, rho, bj, dj, fj)
              if (isnan(bj)) bj = 1.
              integral = &
                 bj**2*sqrt(1.0/pi/2.0)*vy*exp(-vy**2/2.)*(integral+res_v1+res_v2)
           else if (s==1) then
              rho = sqrt(theta*mu_e)*sqrt(ky**2+kx**2)*vy
              call bjndd (n, rho, bj, dj, fj)
              if (isnan(bj)) bj = 1.
              integral = &
                 bj**2*sqrt(1.0/pi/2.0)*vy*exp(-vy**2/2.)*(integral+res_v1+res_v2)
           else if (s==3) then
              rho = sqrt(mu_z)/Z*sqrt(ky**2+kx**2)*vy
              call bjndd (n, rho, bj, dj, fj)
              if (isnan(bj)) bj = 1.
              integral = &
                 bj**2*sqrt(1.0/pi/2.0)*vy*exp(-vy**2/2.)*(integral+res_v1+res_v2)
           else if (s==4.or.s==5.or.s==6) then
              rho = sqrt(theta*mu_e)*sqrt(ky**2+kx**2)*vy
              call bjndd (n, rho, bj, dj, fj)
              if (isnan(bj)) bj = 1.
              integral = &
                 sqrt((2*eps)/(1+eps))*bj**2*sqrt(2.0/pi)*vy**2*exp(-vy**2/2.)*integral
           end if
        end if

  end subroutine

  subroutine vz_integral_residue(vy,omega,s,res_v1,res_v2,shift_vz)

        implicit none

        real(kind=8), intent(in) :: vy
        complex(kind=8), intent(in) :: omega
        integer(kind=8), intent(in) :: s
        complex (kind=8), intent(out) :: res_v1, res_v2, shift_vz
        real(kind=8) :: real_discrim
        real(kind=8) :: v0,coeff_a,coeff_b,coeff_d
        complex(kind=8) :: dv,v1,v2
        complex(kind=8) :: coeff_c, coeff_e, coeff_f, coeff_g, tmp_f, tmp_g
        complex(kind=8) :: tmp_res1, tmp_res2

        if (s==1) then
           fprime = fprime_e
           tprime = tprime_e
        else
           fprime = fprime_i
           tprime = tprime_i
        end if
        if (s==2) then
          if (kv==0.) then
           if (kz==0.) then
             v1 = 100.
             v2 = 100.
             tmp_res1 = 0.
             tmp_res2 = 0.
           else
             v1 = omega/kz
             v2 = 0.
             tmp_res1 = 2.*pi*zi*exp(-v1**2/2.)*&
                        (-1./kz*(omega-fprime*ky-&
                        (-3./2.+v1**2/2.+vy**2/2.)*&
                        tprime*ky))
             tmp_res2 = 0.
           end if
          else
           coeff_a = kv/A
           coeff_b = kz
           coeff_c = 0.5*kv/A*vy**2-omega
           coeff_d = 0.5*tprime*ky
           coeff_e = (-3./2.+vy**2/2.)*tprime*ky+fprime*ky-omega
           v0 = -coeff_b/2./coeff_a
           dv = sqrt(coeff_b**2-4.*coeff_a*coeff_c)/coeff_a
           real_discrim = real(coeff_b**2-4.*coeff_a*coeff_c)
           if (real_discrim<0..and.aimag(omega)<0.) then
              v1 = v0-0.5*dv
              v2 = v0+0.5*dv
           else
              v1 = v0+0.5*dv
              v2 = v0-0.5*dv
           end if
           tmp_res1 = 2*pi*zi*exp(-v1**2/2.)*(coeff_d*&
                      v1**2+coeff_e)/coeff_a/(v1-v2)
           tmp_res2 = 2*pi*zi*exp(-v2**2/2.)*(coeff_d*&
                      v2**2+coeff_e)/coeff_a/(v2-v1)
          end if
        else if (s==1) then
          if (kv==0.) then
           if (kz==0.) then
             v1 = 100.
             v2 = 100.
             tmp_res1 = 0.
             tmp_res2 = 0.
           else
             v1 = omega/kz/sqrt(theta/mu_e)
             v2 = 0.
             tmp_res1 = 2.*pi*zi*exp(-v1**2/2.)*&
                        (-1./kz/sqrt(theta/mu_e)*&
                        (omega+theta*fprime*ky+&
                        (-3./2.+v1**2/2.+vy**2/2.)*&
                        theta*tprime*ky))
             tmp_res2 = 0.
           end if
          else
           coeff_a = -theta*kv/A
           coeff_b = sqrt(theta/mu_e)*kz
           coeff_c = -theta*0.5*kv/A*vy**2-omega
           coeff_d = -theta*0.5*tprime*ky
           coeff_e = -theta*((-3./2.+vy**2/2.)*tprime*ky+fprime*ky)-omega
           v0 = -coeff_b/2./coeff_a
           dv = sqrt(coeff_b**2-4.*coeff_a*coeff_c)/coeff_a
           real_discrim = real(coeff_b**2-4.*coeff_a*coeff_c)
           if (real_discrim<0..and.aimag(omega)<0.) then
              v1 = v0-0.5*dv
              v2 = v0+0.5*dv
           else
              v1 = v0+0.5*dv
              v2 = v0-0.5*dv
           end if
           tmp_res1 = 2*pi*zi*exp(-v1**2/2.)*(coeff_d*&
                      v1**2+coeff_e)/coeff_a/(v1-v2)
           tmp_res2 = 2*pi*zi*exp(-v2**2/2.)*(coeff_d*&
                      v2**2+coeff_e)/coeff_a/(v2-v1)
          end if
        else if (s==3) then
          if (kv==0.) then
           if (kz==0.) then
             v1 = 100.
             v2 = 100.
             tmp_res1 = 0.
             tmp_res2 = 0.
           else
             v1 = omega/kz/sqrt(1./mu_z)
             v2 = 0.
             tmp_res1 = 2.*pi*zi*exp(-v1**2/2.)*&
                        (-1./kz/sqrt(1./mu_z)*&
                        (omega-1./Z*fprime*ky-&
                        (-3./2.+v1**2/2.+vy**2/2.)*&
                        1./Z*tprime*ky))
             tmp_res2 = 0.
           end if
          else
           coeff_a = 1./Z*kv/A
           coeff_b = sqrt(1./mu_z)*kz
           coeff_c = 1./Z*0.5*kv/A*vy**2-omega
           coeff_d = 1./Z*0.5*tprime*ky
           coeff_e = 1./Z*((-3./2.+vy**2/2.)*tprime*ky+fprime*ky)-omega
           v0 = -coeff_b/2./coeff_a
           dv = sqrt(coeff_b**2-4.*coeff_a*coeff_c)/coeff_a
           real_discrim = real(coeff_b**2-4.*coeff_a*coeff_c)
           if (real_discrim<0..and.aimag(omega)<0.) then
              v1 = v0-0.5*dv
              v2 = v0+0.5*dv
           else
              v1 = v0+0.5*dv
              v2 = v0-0.5*dv
           end if
           tmp_res1 = 2*pi*zi*exp(-v1**2/2.)*(coeff_d*&
                      v1**2+coeff_e)/coeff_a/(v1-v2)
           tmp_res2 = 2*pi*zi*exp(-v2**2/2.)*(coeff_d*&
                      v2**2+coeff_e)/coeff_a/(v2-v1)
          end if
        end if

        !print*, 'roots v1,v2'
        !print*,v1,v2
        if (aimag(v1) > 0.) then
           shift_vz = 0.
           res_v1 = 0.
           res_v2 = 0.
        else if (aimag(v1) > -vi) then
           shift_vz = -zi*2.*vi
           res_v1 = 0.
           res_v2 = -tmp_res2
        else
           shift_vz = 0.
           res_v1 = tmp_res1
           res_v2 = -tmp_res2
        end if
        !print*, 'residues'
        !print*,res_v1,res_v2

  end subroutine

  function integrand(vz_prime,vy,omega,shift_vz,s,t)

        real(kind=8) :: vz_prime,vy
        complex(kind=8) :: shift_vz,vz
        complex(kind=8) :: omega,integrand,int_tmp
        complex(kind=8) :: omega_star_n, omega_star_t, omega_d
        complex(kind=8) :: omega_parallel, coll_te
        integer(kind=8) :: s,t

        if (s==1.or.s==4.or.s==5.or.s==6) then
           fprime = fprime_e
           tprime = tprime_e
        else
           fprime = fprime_i
           tprime = tprime_i
        end if

        vz = vz_prime + shift_vz

        if (s==2) then
                omega_star_n = fprime*ky
                omega_star_t = 0.5*(-3.0+vy**2+&
                        vz**2)*tprime*ky
                omega_d = kv/A*(0.5*vy**2+vz**2)
                omega_parallel = kz*vz
        else if (s==1) then
                omega_star_n = -theta*fprime*ky
                omega_star_t = -theta*0.5*(-3.0+vy**2+&
                        vz**2)*tprime*ky
                omega_d = -theta*kv/A*(0.5*vy**2+vz**2)
                omega_parallel = sqrt(theta/mu_e)*kz*vz
        else if (s==3) then
                omega_star_n = 1.0/Z*fprime*ky
                omega_star_t = 1.0/Z*0.5*(-3.0+vy**2+&
                        vz**2)*tprime*ky
                omega_d = 1.0/Z*kv/A*(0.5*vy**2+vz**2)
                omega_parallel = sqrt(1.0/mu_z)*kz*vz
        else if (s==4.or.s==5.or.s==6) then
                omega_star_n = -theta*fprime*ky
                omega_star_t = -theta*(-3./2.+vy**2/2.)&
                        *tprime*ky
                omega_d = -theta*kv/A*vy**2/2.
                coll_te = zi*3.*nu0/eps/vy**3
                
        end if
        
        if (t==0.and.(.not.s==4).and.(.not.s==5).and.(.not.s==6)) then
                int_tmp = exp(-vz**2/2.)*&
                  (omega-omega_star_n-omega_star_t)/&
                  (omega-omega_parallel-omega_d)
        else if (t==0.and.s==4) then
                int_tmp = (omega-omega_star_n-omega_star_t)/&
                  (omega+coll_te)
        else if (t==0.and.s==5) then
                int_tmp = (omega-omega_star_n-omega_star_t)/&
                  coll_te
        else if (t==0.and.s==6) then
                int_tmp = (omega-omega_star_n-omega_star_t)/&
                  omega
        else if (t==1) then
                int_tmp = exp(-vz**2/2.)
        else if (t==2) then
                int_tmp = exp(-vz**2/2.)*&
                   (omega_d-omega_star_n-omega_star_t)
        else if (t==3) then
                int_tmp = exp(-vz**2/2.)*&
                   (-omega_d*(omega_star_n+omega_star_t))
        end if

        integrand = int_tmp

  end function

  subroutine read_input_file

        use file_utils, only: input_unit_exist, input_unit

        implicit none

        integer(kind=8) :: in_file
        logical :: exist

        namelist / parameters / sec_nsteps,int_nsteps,nkz,dkz,kz0,kz1, &
                                nky,dky,ky0,ky1,dfprime_i,nfprime_i, &
                                fprime0_i,dtprime_i,ntprime_i, &
                                tprime0_i,dfprime_e,nfprime_e, &
                                fprime0_e,dtprime_e,ntprime_e, &
                                tprime0_e, nomd,domd,omd0, &
                                vy_1,vy_2,vz_1,vz_2,om1_re,om1_im, &
                                eps,nu0,A,vi,int_tol,sec_tol,fr, &
                                theta,Zeff,Z,mu_e,mu_z,na_e,na_z,na_i, &
                                nkx, dkx, kx0, kx1, omd_kx, cut_buffer, &
                                seed_fr, int_nsteps, lower_bnd, t_e
        
        sec_nsteps = 10
        int_nsteps = 1000

        nkz = 1
        dkz = 0.1
        kz0 = 0.0

        nky = 1
        dky = 0.1
        ky0 = 0.0

        dfprime_i = 0.5
        nfprime_i = 1.0
        fprime0_i = 0.0

        ntprime_i = 1
        dtprime_i = 1.0
        tprime0_i = 0.0

        dfprime_e = 0.5
        nfprime_e = 1.0
        fprime0_e = 0.0

        ntprime_e = 1
        dtprime_e = 1.0
        tprime0_e = 0.0

        nomd = 1
        domd = 1.0
        omd0 = 0.0

        vy_1 = 0.0
        vy_2 = 0.0
        vz_1 = 8.0
        vz_2 = 8.0

        om1_re = 1.0
        om1_im = 0.1

        A = 3.0
        nu0 = 0.05
        eps = 0.2
        vi  = 1.0
        sec_tol = 1.0E-05
        int_tol = 1.0E-05
        fr = 0.5

        theta = 1.0
        Zeff = 1.65
        Z = 5
        mu_e = 5.45D-04
        mu_z = 10.8
        na_e = 0.0
        na_z = 1.0
        na_i = 1.0
        t_e = 1.0
        
        nkx = 1
        dkx = 0.5
        kx0 = 0.0
        kx1 = 10.0
        omd_kx = 0.0
        
        cut_buffer = 0.001
        seed_fr = 0.95
        lower_bnd = 1E-16
        
    in_file = input_unit_exist ("parameters", exist)
    if (exist) read (unit=input_unit("parameters"), nml=parameters)

  end subroutine read_input_file

subroutine init_grids

        implicit none

        allocate(ky_grid(nky), kz_grid(nkz), kx_grid(nkx))

        do ikz = 1, nkz
                kz_grid(ikz) = ikz*dkz+kz0
        end do
        do iky = 1, nky
                ky_grid(iky) = iky*dky+ky0
        end do
        do ikx = 1, nkx
                kx_grid(ikx) = ikx*dkx+kx0
        end do

        allocate(fprime_i_grid(nfprime_i), tprime_i_grid(ntprime_i))
        allocate(fprime_e_grid(nfprime_e), tprime_e_grid(ntprime_e))
        allocate(omd_grid(nomd))

        do ifprime_i = 1, nfprime_i
                fprime_i_grid(ifprime_i) = ifprime_i*dfprime_i+fprime0_i
        end do 
        do itprime_i = 1, ntprime_i
                tprime_i_grid(itprime_i) = itprime_i*dtprime_i+tprime0_i
        end do
        do ifprime_e = 1, nfprime_e
                fprime_e_grid(ifprime_e) = ifprime_e*dfprime_e+fprime0_e
        end do 
        do itprime_e = 1, ntprime_e
                tprime_e_grid(itprime_e) = itprime_e*dtprime_e+tprime0_e
        end do 
        do iomd = 1, nomd
                omd_grid(iomd) = iomd*domd+omd0
        end do

        allocate(ikx_ref(1), ikz_ref(1), iky_ref(1))


end subroutine

subroutine init_output_files

        call open_output_file(gamma_unit,'.datgam')
        write (gamma_unit,'(10a12)') "tprime_i","fprime_i","tprime_e","fprime_e",&
                "omd","kx","ky","kz","omega","gamma"
        !call open_output_file(Dmixing_unit,'.datDmixing')
        !write (Dmixing_unit,'(9a12)') "tprime_i","fprime_i","tprime_e","fprime_e",&
        !        "omd","kx","ky","kz","Dmixing"
        call open_output_file(Dmixing_unit,'.te')
        write (Dmixing_unit,'(13a12)') "nu0","omega","","full_int","",&
                "s_nu_int","","s_nu_int_a","","l_nu_int","","l_nu_int_a",""
        call open_output_file(out_unit,'.dat')
        !write (out_unit,'(10a12)') "tprime_i","fprime_i","tprime_e","fprime_e",&
        !        "omd","kx","ky","kz","omega","gamma"
        write (out_unit,'(14a12)') "tprime_i","fprime_i","tprime_e",&
                "ky","omega","gamma","omlin","gmlin","omlina","gmlina",&
                "omquad","gmquad","omquada","gmquada" 
end subroutine

end program
