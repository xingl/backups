program vy_test

        use file_utils

        implicit none

        real (kind=8), parameter :: pi=3.141592653589793
        complex (kind=8), parameter :: zi=(0.0,1.0)
        integer (kind=8) :: nkz, nky, nkx
        real (kind=8) :: dkz, dky, dkx
        integer (kind=8) :: ikz, iky, ikx
        real (kind=8) :: ky_start, kz_start, kx_start
        real (kind=8) :: vy_1, vy_2, vz_1, vz_2
        complex (kind=8) :: omega,integral1,integral2,integral3,integral4
        real (kind=8) :: A, fprime, tprime, omd, ky, kz, kx, omd_kx, kv
        real (kind=8) :: vi, tol, fr, vz_line
        integer (kind=8) :: nstep, steps
        integer :: gamma_unit=101, Dmixing_unit=102, out_unit=103
        real (kind=8), dimension(:), allocatable :: ky_grid, kz_grid, kx_grid
        real (kind=8) :: om1_re, om1_im
        integer :: s

        call init_file_utils
        call read_input_file
        call init_grids

        call open_output_file(gamma_unit,'.dat')
        write (gamma_unit,'(7a12)') 'kx','ky','kz','omega','',&
                        'int2',''
!        write (gamma_unit,'(10a12)') 'kx','omega','','vy','res1','','res2','',&
!                        'int2',''
        call open_output_file(out_unit,'.res')
        write (out_unit,'(13a12)') 'omega','','vy','coeff_a','discrim','v1','',&
                        'v2','','res1','','res2',''

        do iky = 1, nky
           ky = ky_grid(iky)
           do ikz = 1, nkz
              kz = kz_grid(ikz)
              do ikx = 1, nkx
                 kx = kx_grid(ikx)
                 kv = (1-omd_kx)*ky + omd_kx*kx
                 omega = om1_re + zi*om1_im
                 do while (aimag(omega)>-om1_im)
!                     print*, omega
                     call vy_integral(omega,integral2,s=2)
!                     print*, integral2
                      
                     write(gamma_unit,'(7e12.4)') kx,ky,kz,omega,&
                                integral2
                     omega = omega - 0.1*zi
                 end do
              end do
           end do
        end do
        call close_output_file(gamma_unit)
        call close_output_file(out_unit)
        call finish_file_utils

contains

  subroutine vy_integral(omega,integral,s)

        implicit none

        integer, intent(in) :: s
        complex(kind=8), intent(in) :: omega
        complex (kind=8), intent(out) :: integral
        real(kind=8) :: x_left, x_center, x_right, dx
        complex (kind=8) :: f_left, f_center, f_right
        complex (kind=8) :: i_trapezoid, i_simpson
        real(kind=8) :: int_a,int_b
        integer :: vy_isteps

        int_a = vy_1
        int_b = vy_2
        dx = 0.1
        x_left = int_a
        x_right = x_left + dx
        call vz_integral(omega,x_left,f_left,s)
        integral = 0.0
        vy_isteps = 0

        do while (x_left<int_b)
                x_center = 0.5*(x_left+x_right)
                call vz_integral(omega,x_center,f_center,s)
                call vz_integral(omega,x_right,f_right,s)
                i_trapezoid = 0.5*(f_left+2.0*f_center+f_right)*0.5*dx
                i_simpson = 1.0/6.0*dx*(f_left+4.0*f_center+f_right)
                do while (abs(i_trapezoid-i_simpson)>tol)
                        dx=fr*dx*(abs(i_trapezoid-i_simpson)/tol)**(-1.0/3.0)
                        if (dx>0.2) dx=0.2
                        x_right = x_left + dx
                        if (x_right>int_b) then
                                x_right = int_b
                                dx = x_right - x_left
                        end if
                        x_center = 0.5*(x_left+x_right)
                        call vz_integral(omega,x_center,f_center,s)
                        call vz_integral(omega,x_right,f_right,s)
                        i_trapezoid = 0.5*(f_left+2.0*f_center+f_right)*0.5*dx
                        i_simpson = 1.0/6.0*dx*(f_left+4.0*f_center+f_right)
                end do
                integral = integral + i_simpson
                x_left = x_right
                vy_isteps = vy_isteps + 1
                !print*, 'vy,integral,omega,kx,ky,kz,vy_isteps'
                !print*, x_left, integral,omega,kx,ky,kz,vy_isteps
                if (dx>0.2) dx=0.2
                x_right = x_left + dx
                if (x_right>int_b) then
                        x_right = int_b
                        dx = x_right - x_left
                end if
                f_left = f_right

        end do

  end subroutine

  subroutine vz_integral(omega,vy,integral,s)

        implicit none

        integer, intent(in) :: s
        real(kind=8), intent(in) :: vy
        complex(kind=8), intent(in) :: omega
        complex (kind=8), intent(out) :: integral
        real(kind=8) :: x_left, x_center, x_right, dx
        complex (kind=8) :: f_left, f_center, f_right
        complex (kind=8) :: i_trapezoid, i_simpson
        real(kind=8) :: int_a,int_b
        real(kind=8) :: rho,bj,dj,fj
        integer :: n
        complex (kind=8) :: res_v1, res_v2, shift_vz

        if (s==2) then
           call vz_integral_residue(vy,omega,s,res_v1,res_v2,&
                shift_vz)
           !print*, vy, shift_vz
        end if

        int_a = vz_1
        int_b = vz_2
        dx = 0.1
        x_left = int_a
        x_right = x_left + dx
        call integrand(x_left,vy,omega,s,f_left,shift_vz)
        integral = 0.0

        do while (x_left<int_b)
                x_center = 0.5*(x_left+x_right)
                call integrand(x_center,vy,omega,s,f_left,shift_vz)
                call integrand(x_right,vy,omega,s,f_left,shift_vz)
                i_trapezoid = 0.5*(f_left+2.0*f_center+f_right)*0.5*dx
                i_simpson = 1.0/6.0*dx*(f_left+4.0*f_center+f_right)
                do while (abs(i_trapezoid-i_simpson)>tol)
                        dx=fr*dx*(abs(i_trapezoid-i_simpson)/tol)**(-1.0/3.0)
                        if (dx>0.2) dx=0.2
                        x_right = x_left + dx
                        if (x_right>int_b) then
                                x_right = int_b
                                dx = x_right - x_left
                        end if
                        x_center = 0.5*(x_left+x_right)
                        call integrand(x_center,vy,omega,s,f_left,shift_vz)
                        call integrand(x_right,vy,omega,s,f_left,shift_vz)
                        i_trapezoid = 0.5*(f_left+2.0*f_center+f_right)*0.5*dx
                        i_simpson = 1.0/6.0*dx*(f_left+4.0*f_center+f_right)
                end do
                integral = integral + i_simpson
                x_left = x_right
                dx = fr*dx*(abs(i_trapezoid-i_simpson)/tol)**(-1.0/3.0)
                if (dx>0.2) dx=0.2
                x_right = x_left + dx
                if (x_right>int_b) then 
                        x_right = int_b
                        dx = x_right - x_left
                end if
                f_left = f_right

        end do

        n=0
        if (s==2) then
           rho = sqrt(ky**2+kx**2)*vy
           call bjndd (n, rho, bj, dj, fj)
           if (isnan(bj)) bj = 1.
        end if
        integral = &
           bj**2*sqrt(1.0/pi/2.0)*vy*exp(-0.5*vy**2)*(integral+res_v1+res_v2)

!         write(gamma_unit,'(10e12.4)') kx,omega,vy,res_v1,res_v2,integral

  end subroutine

  subroutine vz_integral_residue(vy,omega,s,res_v1,res_v2,shift_vz)

        implicit none

        real(kind=8), intent(in) :: vy
        complex(kind=8), intent(in) :: omega
        integer, intent(in) :: s
        real(kind=8) :: real_discrim
        complex (kind=8), intent(out) :: res_v1, res_v2, shift_vz
        real(kind=8) :: v0,coeff_a,coeff_b,coeff_d
        complex(kind=8) :: dv,v1,v2
        complex(kind=8) :: coeff_c, coeff_e, coeff_f, coeff_g, tmp_f, tmp_g
        complex(kind=8) :: tmp_res1, tmp_res2

        if (s==2) then
          if (kv==0.) then
             v2 = omega/kz
             v1 = 0.
             tmp_res2 = 2*pi*zi*exp(-v2**2/2.)*(-1./kz*(omega-fprime*ky-(-3./2.+v2**2/2.+vy**2/2.)*&
                       fprime*ky))
             tmp_res1 = 0.
          else
           coeff_a = kv/A
           coeff_b = kz
           coeff_c = 0.5*kv/A*vy**2-omega
           coeff_d = 0.5*tprime*ky
           coeff_e = (-3./2.+vy**2/2.)*tprime*ky+fprime*ky-omega
           v0 = -coeff_b/2./coeff_a
           dv = sqrt(coeff_b**2-4.*coeff_a*coeff_c)/coeff_a
           real_discrim = real(coeff_b**2-4.*coeff_a*coeff_c)
           if (aimag(dv)<0.) then
              v1 = v0-dv/2.
              v2 = v0+dv/2.
           else
              v1 = v0+dv/2.
              v2 = v0-dv/2.
           end if

           tmp_res1 = 0.
           tmp_res2 = zi*exp(-v2**2/2.)*(coeff_d*v2**2+coeff_e)/coeff_a/(v2-v1)

          end if
        end if

        if (aimag(v2) > 0.) then
           shift_vz = 0.
           res_v1 = 0.
           res_v2 = 0.
        else if (aimag(v2) > -vz_line) then
           shift_vz = -zi*2.*vz_line
           res_v1 = 0.
           res_v2 = 0.
        else
           shift_vz = 0.
           res_v1 = 0.
           res_v2 = tmp_res2
        end if

        write(out_unit,'(13e12.4)') omega,vy,coeff_a,real_discrim,v1,v2,res_v1,res_v2
  end subroutine

    subroutine integrand(vz_prime,vy,omega,s,int_tmp,shift_vz)

        real(kind=8), intent(in) :: vz_prime,vy
        integer, intent(in) :: s
        complex(kind=8), intent(out) :: int_tmp
        complex(kind=8), intent(in) :: omega,shift_vz
        complex(kind=8) :: omega_star_n, omega_star_t, omega_star_d
        complex(kind=8) :: omega_parallel,vz
        real(kind=8) :: rho,bj,dj,fj
        integer :: n

        vz = vz_prime + shift_vz
        if (s==2) then
                omega_star_n = fprime*ky
                omega_star_t = 0.5*(-3.0+vy**2+&
                        vz**2)*tprime*ky
                omega_star_d = omd*kv/A*(0.5*vy**2+vz**2)
                omega_parallel = kz*vz
        end if


        int_tmp = exp(-vz**2/2.)*&
                  (omega-omega_star_n-omega_star_t)/&
                  (omega-omega_parallel-omega_star_d)


  end subroutine

  subroutine read_input_file

        use file_utils, only: input_unit_exist, input_unit

        implicit none

        integer :: in_file
        logical :: exist

        namelist / parameters / nstep,nkz,dkz,kz_start, &
                                nky,dky,ky_start, &
                                vy_1,vy_2,vz_1,vz_2,om1_re,om1_im, &
                                A,vi,tol,fr, &
                                nkx, dkx, kx_start, &
                                tprime, fprime, omd, omd_kx, vz_line
        
        nstep = 10

        nkz = 1
        dkz = 0.1
        kz_start = 0.0

        nky = 1
        dky = 0.1
        ky_start = 0.0

        vy_1 = 0.0
        vy_2 = 5.0
        vz_1 = 8.0
        vz_2 = 8.0

        om1_re = 1.0
        om1_im = 0.1

        A = 3.0
        vi  = 1.0
        tol = 1.0E-05
        fr = 0.5
        
        nkx = 1
        dkx = 0.5
        kx_start = 0.0

        tprime = 15.
        fprime = 7.5
        omd = 1.0
        omd_kx = 1.0

        vz_line = 0.1
        
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

end subroutine

end program
