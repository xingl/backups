program imp_tor_itg

        use file_utils

        implicit none

        real (kind=8), parameter :: pi=3.141592653589793
        complex (kind=8), parameter :: zi=(0.0,1.0)

        integer (kind=8) :: nkz, nky, nkx
        real (kind=8) :: dkz, dky, dkx
        integer (kind=8) :: ikz, iky, ikx
        real (kind=8) :: ky_start, kz_start, kx_start

        integer (kind=8) :: nfprime, ntprime, nomd
        real (kind=8) :: dfprime, dtprime, domd
        integer (kind=8) :: ifprime, itprime, iomd
        real (kind=8) :: tprime_start, fprime_start, omd_start

        real (kind=8) :: vy_1, vy_2, vz_1, vz_2

        complex (kind=8) :: omega 
        real (kind=8) :: A, fprime, tprime, omd, ky, kz, kx, kv
        real (kind=8) :: vi, tol, fr, na_e, na_z, na_i, omd_kx
        real (kind=8) :: theta, Zeff, Z, mu_e, mu_z
        integer (kind=8) :: nstep, steps
        integer :: gamma_unit=101, Dmixing_unit=102, out_unit=103

        real (kind=8), dimension(:), allocatable :: ky_grid, kz_grid, kx_grid
        real (kind=8), dimension(:), allocatable :: kx_lp, kx_lp_D
        real (kind=8), dimension(:), allocatable :: fprime_grid, tprime_grid
        real (kind=8), dimension(:), allocatable :: Dmixing_lp, omd_grid
        complex (kind=8), dimension(:), allocatable :: omega_lp
        integer (kind=8), dimension(:), allocatable :: ikx_lp_gamma,ikx_lp_Dmixing

        complex (kind=8) :: seed1, seed2
        real (kind=8) :: om1_re, om1_im, om2_re, om2_im
        integer :: nlp
        complex (kind=8) :: root, root_kz, root_ky_kz
        real(kind=8) :: kx_ref
        
        call init_file_utils
        call read_input_file
        call init_grids

        call open_output_file(gamma_unit,'.datgam')
        write (gamma_unit,'(8a12)') "tprime","fprime","omd","kx",&
                                "ky","kz","omega","gamma"
        call open_output_file(Dmixing_unit,'.datDmixing')
        write (Dmixing_unit,'(9a12)') "tprime","fprime","omd","kx",&
                                "ky","kz","omega","gamma","Dmixing"
        call open_output_file(out_unit,'.dat')
        write (out_unit,'(8a12)') "tprime","fprime","omd","kx",&
                                "ky","kz","omega","gamma"

        root = om1_re + om1_im*zi
        root_kz = om1_re + om1_im*zi
        root_ky_kz = om1_re + om1_im*zi

        do iomd = 1, nomd
           omd = omd_grid(iomd)
           do ifprime = 1, nfprime
              fprime = fprime_grid(ifprime)
              do itprime = 1, ntprime
                 tprime = tprime_grid(itprime)
                 do iky = 1, nky
                       ky = ky_grid(iky)
                    do ikz = 1, nkz
                       kz = kz_grid(ikz)
                       print*, "ky, kz"
                       print*, ky, kz

                       kx_lp = 0.0; omega_lp = 0.0; Dmixing_lp = 0.0
                       do nlp = 1,4
                          call kx_scan(nlp,kx_lp(nlp),omega_lp(nlp),kx_lp_D(nlp),Dmixing_lp(nlp))
                          write(out_unit,*)
                       end do
                       write(out_unit,*)
                       write(out_unit,*)

                       ikx_lp_gamma = maxloc(aimag(omega_lp))
                       ikx_lp_Dmixing = maxloc(Dmixing_lp)

                       if (.not.maxval(aimag(omega_lp))==-9999.0) write (gamma_unit,&
                           '(8e12.4)') tprime,fprime,omd, &
                           kx_lp(ikx_lp_gamma(1)),ky,kz,omega_lp(ikx_lp_gamma(1))
                       if (.not.maxval(Dmixing_lp)==-9999.0) write (Dmixing_unit, &
                           '(9e12.4)') tprime,fprime,omd, &
                           kx_lp_D(ikx_lp_Dmixing(1)),ky,kz,omega_lp(ikx_lp_Dmixing(1)),&
                           Dmixing_lp(ikx_lp_Dmixing(1))

                       if (.not.maxval(aimag(omega_lp))==-9999.0) then  
                           kx_ref = kx_lp(ikx_lp_gamma(1))
                           root_ky_kz = omega_lp(ikx_lp_gamma(1))
                           if (ikz==1) then
                              root_kz = omega_lp(ikx_lp_gamma(1))
                           end if
                       end if

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

 subroutine kx_scan (iloop,kx_maxgam,omega,kx_maxDmix,Dmixing)

        implicit none

        integer, intent(in) :: iloop
        real (kind=8), intent(out) :: kx_maxgam, kx_maxDmix, Dmixing
        complex (kind=8), intent(out) :: omega
        complex (kind=8), dimension(:), allocatable :: omega_kx_scan
        integer (kind=8), dimension(:), allocatable :: ikx_maxgamma,ikx_maxDmixing
        real (kind=8), dimension(:), allocatable :: Dmixing_kx_scan
        real (kind=8) :: kx_init


        allocate(omega_kx_scan(nkx),Dmixing_kx_scan(nkx),ikx_maxgamma(1),ikx_maxDmixing(1))
        omega_kx_scan = -9999.0*zi
        Dmixing_kx_scan = -9999.0

        if (iloop==3) then
            kx_init = 0.0
            do ikx = 1, nkx
                kx_grid(ikx) = (ikx-1)*dkx+kx_init
            end do
        else if (iloop==4) then
            kx_init = 0.0
            do ikx = 1, nkx
                kx_grid(ikx) = -(ikx-1)*dkx+kx_init
            end do
        else if (iloop==1) then
            kx_init = kx_ref
            if (iky==1.and.ikz==1) kx_init = kx_start
            do ikx = 1, nkx
                kx_grid(ikx) = (ikx-1)*dkx+kx_init
            end do
        else if (iloop==2) then
            kx_init = kx_ref
            if (iky==1.and.ikz==1) kx_init = kx_start
            do ikx = 1, nkx
                kx_grid(ikx) = -(ikx-1)*dkx+kx_init
            end do
       end if 

       
       do ikx = 1, nkx
          kx = kx_grid(ikx)
          kv = ((1-omd_kx)*ky + omd_kx*kx)*omd
          print*, "kx,kv"
          print*, kx,kv

          if (ifprime ==1 .and. itprime==1 .and. iomd==1 &
              .and. iky==1 .and. ikz==1 .and. ikx==1) then
            seed1 = om1_re + om1_im*zi
            seed2 = om2_re + om2_im*zi
          else if (ifprime ==1 .and. itprime==1 .and. iomd==1 &
              .and. ikz==1 .and. ikx==1) then
            seed1 = root_kz
            seed2 = root_kz*0.9
          else if (ifprime ==1 .and. itprime==1 .and. iomd==1 &
              .and. ikx==1) then
            seed1 = root_ky_kz
            seed2 = root_ky_kz*0.9
          else
            seed1 = root
            seed2 = root*0.9
          end if

          call rootfinder(seed1,seed2,root,steps)
          if (steps==nstep) then
              print*, 'No root is found.'
              root = seed1
              exit
          else if (aimag(root)<1.E-08) then
              print*, 'No root is found.'
              root = seed1
              exit
          else
              print*, "Root is found:"
              print*, root    
              write (out_unit, '(8e12.4)') tprime,fprime,&
                                        omd,kx,ky,kz,root
              omega_kx_scan(ikx) = root
          end if
       end do

       ikx_maxgamma = maxloc(aimag(omega_kx_scan))
       Dmixing_kx_scan = aimag(omega_kx_scan)/(kx_grid**2+ky**2)
       ikx_maxDmixing = maxloc(Dmixing_kx_scan)

       kx_maxgam = kx_grid(ikx_maxgamma(1))
       omega = omega_kx_scan(ikx_maxgamma(1))
       kx_maxDmix = kx_grid(ikx_maxDmixing(1))
       Dmixing = Dmixing_kx_scan(ikx_maxDmixing(1))

       !print*, "kx, max growth rate, kx max Dmixing"
       !print*, kx_maxgam, omega, kx_maxDmix,Dmixing    
       

       deallocate(omega_kx_scan,Dmixing_kx_scan,ikx_maxgamma,ikx_maxDmixing)

  end subroutine

 subroutine rootfinder(sd1,sd2,rt,istep)

        implicit none

        complex (kind=8), intent(in) :: sd1,sd2
        complex (kind=8), intent(out) :: rt
        integer (kind=8), intent(out) :: istep
        complex (kind=8) :: x_2,f_2,x_1,f_1
        complex (kind=8) :: x_tmp, f_tmp

        x_1 = sd1
        call dispersion_relation(x_1,f_1)
        x_2 = sd2
        call dispersion_relation(x_2,f_2)
        istep = 0
        do while (abs(f_1)>1.0E-04 .and. abs(f_2)>1.0E-04 .and. istep<nstep) 
                x_tmp = x_1 - f_1 * ((x_1-x_2)/(f_1-f_2))
                call dispersion_relation(x_tmp,f_tmp)
                f_2 = f_1
                f_1 = f_tmp
                x_2 = x_1
                x_1 = x_tmp
                istep = istep +1
        end do
        rt = x_tmp

  end subroutine

subroutine dispersion_relation(omega, rhs)
        
        implicit none

        complex(kind=8), intent(in) :: omega
        complex (kind=8), intent(out) :: rhs
        complex (kind=8) :: integral_e, integral_i, integral_z
        integer :: s
        !s = 1
        !call vy_integral(omega,integral_e,s)
        s = 2
        call vy_integral(omega,integral_i,s)
        !s = 3
        !call vy_integral(omega,integral_z,s)

        !rhs = 1.0 + theta*Zeff - na_e*integral_e - &
        !      na_i*theta*(Z-Zeff)/(z-1.0)*integral_i - &
        !      na_z*theta*(Zeff-1.0)/z/(z-1.0)*integral_z
        rhs = 1.0 + theta*Zeff - na_i*integral_i 


end subroutine

subroutine vy_integral(omega,integral,s)

        implicit none

        integer, intent(in) :: s
        complex(kind=8), intent(in) :: omega
        complex (kind=8), intent(out) :: integral
        complex (kind=8) :: integral1, integral2, integral3, integral4
        real (kind=8) :: discrim1, discrim2
        real (kind=8) :: vy_d, vy_l, vy_r, vy_b
        
        vy_b = 1.0E-03
        integral1 = 0.
        integral2 = 0.
        integral3 = 0.
        integral4 = 0.

        if (s==2.and.(.not.kv==0.)) then
           discrim1 = real(kz**2-4.*kv/A*(kv/2./A*vy_1**2-omega))
           discrim2 = real(kz**2-4.*kv/A*(kv/2./A*vy_2**2-omega))
          if (aimag(omega)<0..and.discrim1*discrim2<0.) then
           vy_d = sqrt(2.*A/kx*(real(omega)+kz**2*4.*A/kx))
           vy_l = vy_d - vy_b
           vy_r = vy_d + vy_b
           call vy_integral_simpson(omega,vy_1,vy_l,integral1,s)
           call vy_integral_midpoint(omega,vy_l,vy_d,integral2,s)
           call vy_integral_midpoint(omega,vy_d,vy_r,integral3,s)
           call vy_integral_simpson(omega,vy_r,vy_2,integral4,s)
          else
           call vy_integral_simpson(omega,vy_1,vy_2,integral1,s) 
          end if
        end if
        integral = integral1+integral2+integral3+integral4 

end subroutine

subroutine vy_integral_midpoint(omega,vy_a,vy_b,integral,s)

        implicit none

        integer, intent(in) :: s
        real(kind=8),intent(in) :: vy_a, vy_b
        complex(kind=8), intent(in) :: omega
        complex (kind=8), intent(out) :: integral
        real(kind=8) :: vy_m
        complex(kind=8) :: f_m

        vy_m = (vy_a+vy_b)/2.
        call vz_integral(vy_m,f_m,s,omega)
        integral = f_m*(vy_b-vy_a)

end subroutine

subroutine vy_integral_simpson(omega,vy_a,vy_b,integral,s)

        implicit none

        integer, intent(in) :: s
        real(kind=8),intent(in) :: vy_a, vy_b
        complex(kind=8), intent(in) :: omega
        complex (kind=8), intent(out) :: integral
        real(kind=8) :: x_left, x_center, x_right, dx
        complex (kind=8) :: f_left, f_center, f_right
        complex (kind=8) :: i_trapezoid, i_simpson

        dx = 0.1
        x_left = vy_a
        x_right = x_left + dx
        call vz_integral(x_left,f_left,s,omega)
        integral = 0.0

        do while (x_left<vy_b)
                x_center = 0.5*(x_left+x_right)
                call vz_integral(x_center,f_center,s,omega)
                call vz_integral(x_right,f_right,s,omega)
                i_trapezoid = 0.5*(f_left+2.0*f_center+f_right)*0.5*dx
                i_simpson = 1.0/6.0*dx*(f_left+4.0*f_center+f_right)
                do while (abs(i_trapezoid-i_simpson)>tol)
                        dx=fr*dx*(abs(i_trapezoid-i_simpson)/tol)**(-1.0/3.0)
                        if (dx>0.2) dx=0.2
                        x_right = x_left + dx
                        if (x_right>vy_b) then
                                x_right = vy_b
                                dx = x_right - x_left
                        end if
                        x_center = 0.5*(x_left+x_right)
                        call vz_integral(x_center,f_center,s,omega)
                        call vz_integral(x_right,f_right,s,omega)
                        i_trapezoid = 0.5*(f_left+2.0*f_center+f_right)*0.5*dx
                        i_simpson = 1.0/6.0*dx*(f_left+4.0*f_center+f_right)
                end do
                integral = integral + i_simpson
                x_left = x_right
                dx = fr*dx*(abs(i_trapezoid-i_simpson)/tol)**(-1.0/3.0)
                if (dx>0.2) dx=0.2
                x_right = x_left + dx
                if (x_right>vy_b) then
                        x_right = vy_b
                        dx = x_right - x_left
                end if
                f_left = f_right

        end do

  end subroutine

  subroutine vz_integral(vy,integral,s,omega)

        implicit none

        integer, intent(in) :: s
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
        integer :: n

        if (s==2) then
           call vz_integral_residue(vy,omega,s,res_v1,res_v2,&
                shift_vz)
        else if (s==1) then
           call vz_integral_residue(vy,omega,s,res_v1,res_v2,&
                shift_vz)
        else if (s==3) then
           call vz_integral_residue(vy,omega,s,res_v1,res_v2,&
                shift_vz)
        end if
        print*,'omega,vy,res_v1,res_v2'
        print*, omega,vy,res_v1,res_v2
        vz_a = vz_1
        vz_b = vz_2
        dx = 0.1
        x_left = vz_a
        x_right = x_left + dx
        f_left = integrand(x_left,vy,omega,shift_vz,s)
        integral = 0.0

        do while (x_left<vz_b)
                x_center = 0.5*(x_left+x_right)
                f_center = integrand(x_center,vy,omega,shift_vz,s)
                f_right = integrand(x_right,vy,omega,shift_vz,s)
                i_trapezoid = 0.5*(f_left+2.0*f_center+f_right)*0.5*dx
                i_simpson = 1.0/6.0*dx*(f_left+4.0*f_center+f_right)
                do while (abs(i_trapezoid-i_simpson)>tol)
                        dx=fr*dx*(abs(i_trapezoid-i_simpson)/tol)**(-1.0/3.0)
                        if (dx>0.2) dx=0.2
                        x_right = x_left + dx
                        if (x_right>vz_b) then
                                x_right = vz_b
                                dx = x_right - x_left
                        end if
                        x_center = 0.5*(x_left+x_right)
                        f_center = integrand(x_center,vy,omega,shift_vz,s)
                        f_right = integrand(x_right,vy,omega,shift_vz,s)
                        i_trapezoid = 0.5*(f_left+2.0*f_center+f_right)*0.5*dx
                        i_simpson = 1.0/6.0*dx*(f_left+4.0*f_center+f_right)
                end do
                integral = integral + i_simpson
                x_left = x_right
                dx = fr*dx*(abs(i_trapezoid-i_simpson)/tol)**(-1.0/3.0)
                if (dx>0.2) dx=0.2
                x_right = x_left + dx
                if (x_right>vz_b) then 
                        x_right = vz_b
                        dx = x_right - x_left
                end if
                f_left = f_right

        end do

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
        !res_v1 = 0.
        !res_v2 = 0.
        integral = &
           bj**2*sqrt(1.0/pi/2.0)*vy*exp(-0.5*vy**2)*(integral+res_v1+res_v2)

  end subroutine

  subroutine vz_integral_residue(vy,omega,s,res_v1,res_v2,shift_vz)

        implicit none

        real(kind=8), intent(in) :: vy
        complex(kind=8), intent(in) :: omega
        integer, intent(in) :: s
        complex (kind=8), intent(out) :: res_v1, res_v2, shift_vz
        real(kind=8) :: real_discrim
        real(kind=8) :: v0,coeff_a,coeff_b,coeff_d
        complex(kind=8) :: dv,v1,v2
        complex(kind=8) :: coeff_c, coeff_e, coeff_f, coeff_g, tmp_f, tmp_g
        complex(kind=8) :: tmp_res1, tmp_res2
        real(kind=8) :: vz_line

        vz_line = 0.1
        if (s==2) then
          if (kv==0.) then
           if (kz==0.) then
             v1 = 100.
             v2 = 100.
             tmp_res1 = 0.
             tmp_res2 = 0.
           else
             v1 = 0.
             v2 = omega/kz
             tmp_res1 = 0.
             tmp_res2 = 2.*pi*zi*exp(-v2**2/2.)*&
                        (-1./kz*(omega-fprime*ky-&
                        (-3./2.+v2**2/2.+vy**2/2.)*&
                        tprime*ky))
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
           tmp_res1 = 0.
           tmp_res2 = 2*pi*zi*exp(-v2**2/2.)*(coeff_d*&
                      v2**2+coeff_e)/coeff_a/(v2-v1)
          end if
        else if (s==1) then
          if (kv==0.) then
             v1 = 0.
             v2 = omega/kz/sqrt(theta/mu_e)
             tmp_res1 = 0.
             tmp_res2 = 2.*pi*zi*exp(-v2**2/2.)*&
                        (-1./kz/sqrt(theta/mu_e)*&
                        (omega+theta*fprime*ky+&
                        (-3./2.+v2**2/2.+vy**2/2.)*&
                        theta*tprime*ky))
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
           tmp_res1 = 0.
           tmp_res2 = 2*pi*zi*exp(-v2**2/2.)*(coeff_d*&
                      v2**2+coeff_e)/coeff_a/(v2-v1)
          end if
        else if (s==3) then
          if (kv==0.) then
             v1 = 0.
             v2 = omega/kz/sqrt(1./mu_z)
             tmp_res1 = 0.
             tmp_res2 = 2.*pi*zi*exp(-v2**2/2.)*&
                        (-1./kz/sqrt(1./mu_z)*&
                        (omega-1./Z*fprime*ky-&
                        (-3./2.+v2**2/2.+vy**2/2.)*&
                        1./Z*tprime*ky))
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
           tmp_res1 = 0.
           tmp_res2 = 2*pi*zi*exp(-v2**2/2.)*(coeff_d*&
                      v2**2+coeff_e)/coeff_a/(v2-v1)
          end if
        end if

        if (aimag(v1) > 0.) then
           shift_vz = 0.
           res_v1 = 0.
           res_v2 = 0.
        else if (aimag(v1) > -vz_line) then
           shift_vz = -zi*2.*vz_line
           res_v1 = 0.
           res_v2 = 0.
        else
           shift_vz = 0.
           res_v1 = 0.
           res_v2 = tmp_res2
        end if

  end subroutine

  function integrand(vz_prime,vy,omega,shift_vz,s)

        real(kind=8) :: vz_prime,vy
        complex(kind=8) :: shift_vz,vz
        complex(kind=8) :: omega,integrand,int_tmp
        complex(kind=8) :: omega_star_n, omega_star_t, omega_star_d
        complex(kind=8) :: omega_parallel
        integer :: s

        vz = vz_prime + shift_vz
        if (s==2) then
                omega_star_n = fprime*ky
                omega_star_t = 0.5*(-3.0+vy**2+&
                        vz**2)*tprime*ky
                omega_star_d = kv/A*(0.5*vy**2+vz**2)
                omega_parallel = kz*vz
        else if (s==1) then
                omega_star_n = -theta*fprime*ky
                omega_star_t = -theta*0.5*(-3.0+vy**2+&
                        vz**2)*tprime*ky
                omega_star_d = -theta*kv/A*(0.5*vy**2+vz**2)
                omega_parallel = sqrt(theta/mu_e)*kz*vz
        else if (s==3) then
                omega_star_n = 1.0/Z*fprime*ky
                omega_star_t = 1.0/Z*0.5*(-3.0+vy**2+&
                        vz**2)*tprime*ky
                omega_star_d = 1.0/Z*kv/A*(0.5*vy**2+vz**2)
                omega_parallel = sqrt(1.0/mu_z)*kz*vz
        end if
        
        int_tmp = exp(-vz**2/2.)*&
                  (omega-omega_star_n-omega_star_t)/&
                  (omega-omega_parallel-omega_star_d)

        integrand = int_tmp

  end function

  subroutine read_input_file

        use file_utils, only: input_unit_exist, input_unit

        implicit none

        integer :: in_file
        logical :: exist

        namelist / parameters / nstep,nkz,dkz,kz_start, &
                                nky,dky,ky_start, dfprime,nfprime, &
                                fprime_start,dtprime,ntprime, &
                                tprime_start, nomd,domd,omd_start, &
                                vy_1,vy_2,vz_1,vz_2,om1_re,om1_im, &
                                om2_re,om2_im,A,vi,tol,fr, &
                                theta,Zeff,Z,mu_e,mu_z,na_e,na_z,na_i, &
                                nkx, dkx, kx_start, omd_kx
        
        nstep = 10

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

        A = 3.0
        vi  = 1.0
        tol = 1.0E-05
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

        allocate(kx_lp(4), omega_lp(4), kx_lp_D(4), Dmixing_lp(4))
        
        allocate(ikx_lp_gamma(1),ikx_lp_Dmixing(1))

end subroutine

end program
