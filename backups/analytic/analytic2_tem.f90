program res_itg

        use file_utils

        implicit none

        real (kind=8), parameter :: pi=3.141592653589793
        complex (kind=8), parameter :: zi=(0.0,1.0)

        integer (kind=8) :: nkz, nky, nkx
        real (kind=8) :: dkz, dky, dkx
        integer (kind=8) :: ikz, iky, ikx
        real (kind=8) :: ky_start, kz_start, kx_start
        real (kind=8), dimension(:), allocatable :: ky_grid, kz_grid, kx_grid

        integer (kind=8) :: nfprime, ntprime, nomd
        real (kind=8) :: dfprime, dtprime, domd
        integer (kind=8) :: ifprime, itprime, iomd
        real (kind=8) :: tprime_start, fprime_start, omd_start
        real (kind=8), dimension(:), allocatable :: fprime_grid, tprime_grid
        real (kind=8), dimension(:), allocatable :: omd_grid

        real (kind=8) :: vy_1, vy_2, vz_1, vz_2
        real (kind=8) :: nu0, A, fprime, tprime, omd, ky, kz, kx, kv
        real (kind=8) :: vi, int_tol,sec_tol, fr, na_e, na_z, na_i, omd_kx
        real (kind=8) :: theta, Zeff, Z, mu_e, mu_z
        complex (kind=8) :: seed1, seed2
        real (kind=8) :: om1_re, om1_im, om2_re, om2_im
        integer (kind=8) :: sec_nsteps, sec_istep
        integer (kind=8) :: int_nsteps, int_istep

        real (kind=8) :: kx_gam, kx_Dmix
        real (kind=8) :: Dmixing_kx_max
        complex (kind=8) :: omega_kx_max
        integer(kind=8), dimension(:),allocatable :: ikx_ref, ikz_ref

        complex (kind=8) :: omega 
        complex (kind=8) :: root, root_ky_kz, root_kz
        complex (kind=8), dimension(:),allocatable :: omega_kx_grid
        complex (kind=8), dimension(:),allocatable :: omega_kz_grid
        integer :: gamma_unit=101, Dmixing_unit=102
        integer :: out_unit=103
        real(kind=8) :: cut_buffer, seed_fr, lower_bnd
        real (kind=8) :: kz_maxgam
        
        call init_file_utils
        call read_input_file
        call init_grids

        call open_output_file(gamma_unit,'.datgam')
        write (gamma_unit,'(8a12)') "tprime","fprime","omd","kx",&
                                "ky","kz","omega","gamma"
        call open_output_file(Dmixing_unit,'.datDmixing')
        write (Dmixing_unit,'(7a12)') "tprime","fprime","omd","kx",&
                                "ky","kz","Dmixing"
        call open_output_file(out_unit,'.dat')
        write (out_unit,'(8a12)') "tprime","fprime","omd","kx",&
                                "ky","kz","omega","gamma"

        root = om1_re + om1_im*zi

        do iomd = 1, nomd
           omd = omd_grid(iomd)
           do ifprime = 1, nfprime
              fprime = fprime_grid(ifprime)
              do itprime = 1, ntprime
                 tprime = tprime_grid(itprime)
                 do iky = 1, nky
                    ky = ky_grid(iky)
                    omega_kz_grid = 0.-9999.0*zi

                    do ikz = ikz_ref(1), nkz
                       kz = kz_grid(ikz)
                       print*, "ky, kz"
                       print*, ky, kz
                       call kx_scan(kx_gam,omega_kx_max,kx_Dmix,Dmixing_kx_max)
                       omega_kz_grid(ikz) = omega_kx_max
                    end do
                    write (gamma_unit,*)
                    write (Dmixing_unit,*)

                    do ikz = ikz_ref(1),1,-1
                       kz = kz_grid(ikz)
                       print*, "ky, kz"
                       print*, ky, kz
                       call kx_scan(kx_gam,omega_kx_max,kx_Dmix,Dmixing_kx_max)
                       omega_kz_grid(ikz) = omega_kx_max
                    end do
                    ikz_ref = maxloc(aimag(omega_kz_grid))
                    kz_maxgam = kz_grid(ikz_ref(1))
                    root_kz = omega_kz_grid(ikz_ref(1))
                    print*, 'kz max max'
                    print*, kz_grid(ikz_ref(1))
                    write (gamma_unit,*)
                    write (Dmixing_unit,*)
                 end do
              end do
           end do
        end do

        call close_output_file (gamma_unit)
        call close_output_file (Dmixing_unit)
        call close_output_file (out_unit)
        call finish_file_utils

contains

 subroutine kx_scan (kx_maxgam,omega,kx_maxDmix,Dmixing)

        implicit none

        real (kind=8), intent(out) :: kx_maxgam, kx_maxDmix, Dmixing
        complex (kind=8), intent(out) :: omega
        integer (kind=8), dimension(:), allocatable :: ikx_maxgamma
        integer (kind=8), dimension(:), allocatable :: ikx_maxDmixing
        real (kind=8), dimension(:), allocatable :: Dmixing_kx_scan
        complex (kind=8) :: root_tmp
        real (kind=8) :: Dmixing_tmp


          allocate(Dmixing_kx_scan(nkx))
          allocate(ikx_maxgamma(1),ikx_maxDmixing(1))
          Dmixing_kx_scan = -9999.0
          omega_kx_grid = 0. -9999.*zi
          ikx_ref = minloc(abs(kx_grid-kx_start))
          !print*, ikx_ref
          do ikx = ikx_ref(1), 1, -1
             kx = kx_grid(ikx)
             kv = ((1-omd_kx)*ky + omd_kx*kx)*omd
             print*, "kx,kv"
             print*, kx,kv
             !call pick_seeds
             !print*,seed1,seed2
             !call rootfinder(seed1,seed2,root_tmp)
             call quadratic_solution(root_tmp)
             if (sec_istep==sec_nsteps.or.isnan(real(root_tmp))) then
              print*, 'No root is found.'
              root = seed1
              omega_kx_grid(ikx) = 0.-9999.*zi
             else if (abs(aimag(root_tmp))<lower_bnd) then
              print*, 'No root is found.'
              root = seed1
              omega_kx_grid(ikx) = 0. -9999.*zi
             else if (aimag(root_tmp)==aimag(root_tmp)-1..or.&
                      isnan(aimag(root_tmp))) then
              print*, 'No root is found.'
              root = seed1
              omega_kx_grid(ikx) = 0.-9999.*zi
             else
              print*, "Root is found:"
              print*, root_tmp    
              write (out_unit, '(8e12.4)') tprime,fprime,&
                                        omd,kx,ky,kz,root_tmp
              omega_kx_grid(ikx) = root_tmp
              root = root_tmp
             end if
          end do
          write (out_unit,*) 
          do ikx = ikx_ref(1), nkx
             kx = kx_grid(ikx)
             kv = ((1-omd_kx)*ky + omd_kx*kx)*omd
             print*, "kx,kv"
             print*, kx,kv
             !call pick_seeds
             !print*,seed1,seed2
             !call rootfinder(seed1,seed2,root_tmp)
             call quadratic_solution(root_tmp)
             if (sec_istep==sec_nsteps.or.isnan(real(root_tmp))) then
              print*, 'No root is found.'
              root = seed1
              omega_kx_grid(ikx) = 0.-9999.*zi
             !else if (abs(aimag(root_tmp))<lower_bnd) then
             ! print*, 'No root is found.'
             ! root = seed1
             ! omega_kx_grid(ikx) = 0.-9999.*zi
             else if (aimag(root_tmp)==aimag(root_tmp)-1..or.&
                      isnan(aimag(root_tmp))) then
              print*, 'No root is found.'
              root = seed1
              omega_kx_grid(ikx) = 0.-9999.*zi
             else
              print*, "Root is found:"
              print*, root_tmp    
              write (out_unit, '(8e12.4)') tprime,fprime,&
                                        omd,kx,ky,kz,root_tmp
              omega_kx_grid(ikx) = root_tmp
              root = root_tmp
             end if
          end do
          write (out_unit, *)

          ikx_maxgamma = maxloc(aimag(omega_kx_grid))
          Dmixing_kx_scan = aimag(omega_kx_grid)/(kx_grid**2+ky**2)
          ikx_maxDmixing = maxloc(Dmixing_kx_scan)

          kx_maxgam = kx_grid(ikx_maxgamma(1))
          omega = omega_kx_grid(ikx_maxgamma(1))
          kx_maxDmix = kx_grid(ikx_maxDmixing(1))
          Dmixing = Dmixing_kx_scan(ikx_maxDmixing(1))

          if (.not.aimag(omega)==-9999.0) &
              write (gamma_unit,&
              '(8e12.4)') tprime,fprime,omd, &
              kx_maxgam,ky,kz,omega
          if (.not.Dmixing==-9999.0) &
              write (Dmixing_unit, &
              '(7e12.4)') tprime,fprime,omd, &
              kx_maxDmix,ky,kz,Dmixing
          if (.not.aimag(omega)==-9999.0) then  
              ikx_ref = minloc(abs(kx_grid-kx_maxgam))
              root_ky_kz = omega
          end if
          !print*, kx_maxgam, omega, kx_maxDmix, Dmixing
          deallocate(Dmixing_kx_scan,ikx_maxgamma,ikx_maxDmixing)

  end subroutine

subroutine quadratic_solution(omega)

        implicit none

        complex(kind=8), intent(out) :: omega
        complex(kind=8) :: omega1, omega2
        complex (kind=8)::term1,term2,term3,term4,term5,term6
        complex (kind=8)::term7,term8
        complex (kind=8) :: Delta
        integer(kind=8) :: t,s
        
        s = 2
        t = 1 
        call vy_integral_simpson(omega,vy_1,vy_2,term1,s,t)
        print*,'term1'
        print*, term1
        t = 2
        call vy_integral_simpson(omega,vy_1,vy_2,term2,s,t)
        print*,'term2/omegad'
        print*, term2/(kv/A)
        !print*, term2
        t = 3 
        call vy_integral_simpson(omega,vy_1,vy_2,term3,s,t)
        print*,'term3/omega*n'
        print*, term3/fprime/ky
        !print*, term3
        t = 4
        call vy_integral_simpson(omega,vy_1,vy_2,term4,s,t)
        print*,'term4/omegad/omega*n'
        print*, term4/(kv/A)/fprime/ky
        !print*, term4
        t = 5 
        call vy_integral_simpson(omega,vy_1,vy_2,term5,s,t)
        print*,'term5/omegad/omega*t'
        print*, term5/(kv/A)/tprime/ky
        !print*, term5
        !t = 6 
        !call vy_integral_simpson(omega,vy_1,vy_2,term6,s,t)
        !print*,'term6/omegad**2'
        !print*, term6/(kv/A)**2
        !print*, term6
        !s = 1
        !t = 7
        !call vy_integral_simpson(omega,vy_1,vy_2,term7,s,t)
        !print*,'term7/omega*n/(nu0*A)'
        !print*, term7/fprime/ky/nu0/A
        !print*, term7
        !t = 8 
        !call vy_integral_simpson(omega,vy_1,vy_2,term8,s,t)
        !print*,'term8/omega*t/(nu0*A)'
        !print*, term8/tprime/ky/nu0/A
        !print*, term8

        !term1 = 1. + theta*Zeff - term1 - zi*(term7+term8)
        term1 = 1. + theta - term1 + zi*nu0

        !Delta = (-term2+term3)**2-4.*(term1+zi*(term7+term8))*(term4+term5-term6)
        Delta = (-term2+term3)**2-4.*(term1-zi*nu0)*(term4+term5)
        print*, 'Delta'
        print*, Delta

        !omega1 = (-(-term2+term3)+sqrt((-term2+term3)**2-4.*term1*(term4+term5-term6)))/2./term1
        omega1 = (-(-term2+term3)+sqrt((-term2+term3)**2-4.*term1*(term4+term5)))/2./term1
        print*, 'omega1'
        print*, omega1

        omega2 = (-(-term2+term3)-sqrt((-term2+term3)**2-4.*term1*(term4+term5)))/2./term1
        !omega2 = (-(-term2+term3)-sqrt((-term2+term3)**2-4.*term1*(term4+term5-term6)))/2./term1
        print*, 'omega2'
        print*, omega2

        if (aimag(omega1)>aimag(omega2)) then
           omega = omega1
        else
           omega = omega2
        end if

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

        !print*, 'vy_integral_simpson'
        !print*, omega,kx,ky,kz,kv
        vy_istep = 0.
        dx = 0.1
        x_left = vy_a
        x_right = x_left + dx
        call vz_integral(x_left,f_left,s,t,omega)
        integral = 0.0

        do while (x_left<vy_b.and.vy_istep<int_nsteps)
                x_center = 0.5*(x_left+x_right)
                call vz_integral(x_center,f_center,s,t,omega)
                call vz_integral(x_right,f_right,s,t,omega)
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
                        call vz_integral(x_center,f_center,s,t,omega)
                        call vz_integral(x_right,f_right,s,t,omega)
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
                !print*,vy_istep, x_left,integral
        end do
        if (vy_istep==int_nsteps) then
           int_istep = int_nsteps
           integral = 0.
        end if

  end subroutine

  subroutine vz_integral(vy,integral,s,t,omega)

        implicit none

        integer(kind=8), intent(in) :: s,t
        real(kind=8), intent(in) :: vy
        complex(kind=8), intent(in) :: omega
        real(kind=8) :: vz_a, vz_b
        complex (kind=8), intent(out) :: integral
        real(kind=8) :: x_left, x_center, x_right, dx
        complex (kind=8) :: f_left, f_center, f_right
        complex (kind=8) :: i_trapezoid, i_simpson
        complex (kind=8) :: shift_vz
        real(kind=8) :: rho,bj,dj,fj
        integer(kind=8) :: n
        integer(kind=8) :: vz_istep

        shift_vz = 0.
        vz_istep = 0
        vz_a = vz_1
        vz_b = vz_2
        dx = 0.1
        x_left = vz_a
        x_right = x_left + dx
        f_left = integrand(x_left,vy,omega,shift_vz,s,t)
        integral = 0.0

        if ((.not.t==7).and.(.not.t==8)) then
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
           else if (s==1) then
              rho = sqrt(theta*mu_e)*sqrt(ky**2+kx**2)*vy
              call bjndd (n, rho, bj, dj, fj)
              if (isnan(bj)) bj = 1.
           else if (s==3) then
              rho = sqrt(mu_z)/Z*sqrt(ky**2+kx**2)*vy
              call bjndd (n, rho, bj, dj, fj)
              if (isnan(bj)) bj = 1.
           end if
           if ((.not.t==7).and.(.not.t==8)) then
              integral = &
                 bj**2*sqrt(1.0/pi/2.0)*vy*exp(-0.5*vy**2)*integral
           else
              integral = &
                 sqrt(2.0/pi)*vy**2*exp(-0.5*vy**2)*integral*sqrt(2./A)
           end if
        end if

  end subroutine

  function integrand(vz_prime,vy,omega,shift_vz,s,t)

        real(kind=8) :: vz_prime,vy
        complex(kind=8) :: shift_vz,vz
        complex(kind=8) :: omega,integrand,int_tmp
        complex(kind=8) :: omega_star_n, omega_star_t, omega_d
        complex(kind=8) :: omega_parallel
        integer(kind=8) :: s,t

        vz = vz_prime + shift_vz
        if (s==2) then
                omega_star_n = fprime*ky
                omega_star_t = 0.5*(-3.0+vy**2+&
                        vz**2)*tprime*ky
                omega_d = kv/A*(0.5*vy**2+vz**2)
                omega_parallel = kz*vz
        else if (s==1) then
            if (.not.t==8) then
                omega_star_n = -theta*fprime*ky
                omega_star_t = -theta*0.5*(-3.0+vy**2+&
                        vz**2)*tprime*ky
                omega_d = -theta*kv/A*(0.5*vy**2+vz**2)
                omega_parallel = sqrt(theta/mu_e)*kz*vz
            else
                omega_star_n = -theta*fprime*ky
                omega_star_t = -theta*0.5*(-3.0+vy**2)&
                        *tprime*ky
                omega_d = 0.
                omega_parallel = 0.
            end if
        else if (s==3) then
                omega_star_n = 1.0/Z*fprime*ky
                omega_star_t = 1.0/Z*0.5*(-3.0+vy**2+&
                        vz**2)*tprime*ky
                omega_d = 1.0/Z*kv/A*(0.5*vy**2+vz**2)
                omega_parallel = sqrt(1.0/mu_z)*kz*vz
        end if
        
        if (t==1) then
                int_tmp = exp(-vz**2/2.)
        else if (t==2) then
                int_tmp = exp(-vz**2/2.)*&
                   omega_d
        else if (t==3) then
                int_tmp = exp(-vz**2/2.)*&
                   omega_star_n
        else if (t==4) then
                int_tmp = exp(-vz**2/2.)*&
                   omega_d*omega_star_n
        else if (t==5) then
                int_tmp = exp(-vz**2/2.)*&
                   omega_d*omega_star_t
        else if (t==6) then
                int_tmp = exp(-vz**2/2.)*&
                   omega_d*omega_d
        else if (t==7) then
                int_tmp = A*nu0*omega_star_n/vy**3
        else if (t==8) then
                int_tmp = A*nu0*omega_star_t/vy**3
        end if

        integrand = int_tmp

  end function

  subroutine read_input_file

        use file_utils, only: input_unit_exist, input_unit

        implicit none

        integer(kind=8) :: in_file
        logical :: exist

        namelist / parameters / sec_nsteps,nkz,dkz,kz_start, &
                                nky,dky,ky_start, dfprime,nfprime, &
                                fprime_start,dtprime,ntprime, &
                                tprime_start, nomd,domd,omd_start, &
                                vy_1,vy_2,vz_1,vz_2,om1_re,om1_im, &
                                om2_re,om2_im,A,vi,int_tol,sec_tol,fr, &
                                theta,Zeff,Z,mu_e,mu_z,na_e,na_z,na_i, &
                                nkx, dkx, kx_start, omd_kx, cut_buffer, &
                                seed_fr, int_nsteps, lower_bnd, nu0
        
        sec_nsteps = 10
        int_nsteps = 1000

        nkz = 1
        dkz = 0.1
        kz_start = 0.0

        nky = 1
        dky = 0.1
        ky_start = 0.0

        dfprime = 0.5
        nfprime = 1.0
        fprime_start = 0.0

        ntprime = 1
        dtprime = 1.0
        tprime_start = 0.0

        nomd = 1
        domd = 1.0
        omd_start = 0.0

        vy_1 = 0.0
        vy_2 = 0.0
        vz_1 = 8.0
        vz_2 = 8.0

        om1_re = 1.0
        om1_im = 0.1
        om2_re = 1.0
        om2_im = 0.1

        nu0 = 0.1
        A = 3.0
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
        
        nkx = 1
        dkx = 0.5
        kx_start = 0.0
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
                kz_grid(ikz) = ikz*dkz+kz_start
        end do
        do iky = 1, nky
                ky_grid(iky) = iky*dky+ky_start
        end do
        do ikx = 1, nkx
                kx_grid(ikx) = ikx*dkx+kx_start
        end do

        allocate(fprime_grid(nfprime), tprime_grid(ntprime), omd_grid(nomd))

        do ifprime = 1, nfprime
                fprime_grid(ifprime) = ifprime*dfprime+fprime_start
        end do 
        do itprime = 1, ntprime
                tprime_grid(itprime) = itprime*dtprime+tprime_start
        end do
        do iomd = 1, nomd
                omd_grid(iomd) = iomd*domd+omd_start
        end do

        allocate(omega_kx_grid(nkx))
        allocate(omega_kz_grid(nkz))
        allocate(ikx_ref(1), ikz_ref(1))
        ikz_ref = minloc(abs(kz_grid-kz_start))


end subroutine

end program
