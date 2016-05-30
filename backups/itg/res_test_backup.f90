program imp_tor_itg

        use file_utils

        implicit none

        real (kind=8), parameter :: pi=3.141592653589793
        complex (kind=8), parameter :: zi=(0.0,1.0)
        integer (kind=8) :: nkz, nky, nkx
        real (kind=8) :: dkz, dky, dkx
        integer (kind=8) :: ikz, iky, ikx
        real (kind=8) :: ky_start, kz_start, kx_start
        real (kind=8) :: vy_1, vz_1, vz_2
        complex (kind=8) :: omega,integral1,integral2,integral3 
        real (kind=8) :: A, fprime, tprime, omd, ky, kz, kx, omd_kx
        real (kind=8) :: vi, tol, fr
        integer (kind=8) :: nstep, steps
        integer :: gamma_unit=101, Dmixing_unit=102, out_unit=103
        real (kind=8), dimension(:), allocatable :: ky_grid, kz_grid, kx_grid
        real (kind=8) :: om1_re, om1_im
        integer :: s

        call init_file_utils
        call read_input_file
        call init_grids

        call open_output_file(gamma_unit,'.dat')
        write (gamma_unit,'(11a12)') 'kx','ky','kz','omega','','int1','','int2','','int3',''
        call open_output_file(out_unit,'.res')
        write (out_unit,'(14a12)') 'kx','ky','kz','omega','','vz','v1','',&
                                        'v2','','res1','','res2',''

        do iky = 1, nky
           ky = ky_grid(iky)
           do ikz = 1, nkz
              kz = kz_grid(ikz)
              do ikx = 1, nkx
                 kx = kx_grid(ikx)
                 omega = om1_re + 5.*zi
                 do while (aimag(omega)>-5.)
                     print*, omega
                     call vz_integral(omega,vy_1,integral1,s=1)
                     call vz_integral(omega,vy_1,integral2,s=2)
                     call vz_integral(omega,vy_1,integral3,s=3)
                     write(gamma_unit,'(11e12.4)') kx,ky,kz,omega,&
                                integral1,integral2,integral3
                     omega = omega - 0.02*zi
                 end do
              end do
           end do
        end do
        call close_output_file(gamma_unit)
        call close_output_file(out_unit)
        call finish_file_utils

contains

  subroutine vz_integral(omega,vy,integral,s)

        implicit none

        integer, intent(in) :: s
        real(kind=8), intent(in) :: vy
        complex(kind=8), intent(in) :: omega
        real(kind=8) :: a, b
        complex (kind=8), intent(out) :: integral
        real(kind=8) :: x_left, x_center, x_right, dx
        complex (kind=8) :: f_left, f_center, f_right
        complex (kind=8) :: i_trapezoid, i_simpson
        complex (kind=8) :: res_v1, res_v2

        a = vz_1
        b = vz_2
        dx = 0.1
        x_left = a
        x_right = x_left + dx
        call integrand(x_left,vy,omega,s,res_v1,res_v2,f_left)
        integral = 0.0

        do while (x_left<b)
                x_center = 0.5*(x_left+x_right)
                call integrand(x_center,vy,omega,s,res_v1,res_v2,f_center)
                call integrand(x_right,vy,omega,s,res_v1,res_v2,f_right)
                i_trapezoid = 0.5*(f_left+2.0*f_center+f_right)*0.5*dx
                i_simpson = 1.0/6.0*dx*(f_left+4.0*f_center+f_right)
                do while (abs(i_trapezoid-i_simpson)>tol)
                        dx=fr*dx*(abs(i_trapezoid-i_simpson)/tol)**(-1.0/3.0)
                        if (dx>0.2) dx=0.2
                        x_right = x_left + dx
                        if (x_right>b) then
                                x_right = b
                                dx = x_right - x_left
                        end if
                        x_center = 0.5*(x_left+x_right)
                        call integrand(x_center,vy,omega,s,res_v1,res_v2,f_center)
                        call integrand(x_right,vy,omega,s,res_v1,res_v2,f_right)
                        i_trapezoid = 0.5*(f_left+2.0*f_center+f_right)*0.5*dx
                        i_simpson = 1.0/6.0*dx*(f_left+4.0*f_center+f_right)
                end do
                integral = integral + i_simpson
                x_left = x_right
                dx = fr*dx*(abs(i_trapezoid-i_simpson)/tol)**(-1.0/3.0)
                if (dx>0.2) dx=0.2
                x_right = x_left + dx
                if (x_right>b) then 
                        x_right = b
                        dx = x_right - x_left
                end if
                f_left = f_right

        end do
        if (s==1) integral = integral + res_v1 - res_v2
        if (s==2) integral = integral + res_v1
        if (s==3) integral = integral - res_v2

  end subroutine

  subroutine integrand(vz,vy,omega,s,res1,res2,int_tmp)

        real(kind=8),intent(in) :: vz,vy
        complex(kind=8),intent(in) :: omega
        integer, intent(in) :: s
        complex(kind=8),intent(out) :: res1, res2, int_tmp
        complex(kind=8) :: omega_star_n, omega_star_t, omega_star_d
        complex(kind=8) :: omega_parallel
        real(kind=8) :: rho,bj,dj,fj
        integer :: n
        real(kind=8) :: kv,v0,alpha,coeff_a,coeff_b,coeff_d
        complex(kind=8) :: dv,v1,v2,vzprime
        complex(kind=8) :: coeff_c, coeff_e, coeff_f, coeff_g, tmp_f, tmp_g
        complex(kind=8) :: tmp_res1, tmp_res2

        kv = (1.0-omd_kx)*ky + omd_kx*kx
        coeff_a = kv/A
        coeff_b = kz
        coeff_c = 0.5*kv/A*vy**2-omega
        coeff_d = 0.5*tprime*ky
        coeff_e = (-3./2.+vy**2/2.)*tprime*ky+fprime*ky-omega 
        v0 = -coeff_b/2./coeff_a
        dv = sqrt(coeff_b**2-4.*coeff_a*coeff_c)/coeff_a
        alpha = atan(aimag(dv)/real(dv))
        if (coeff_a>0) then
                v1 = v0+0.5*dv
                v2 = v0-0.5*dv
        else
                v1 = v0-0.5*dv
                v2 = v0+0.5*dv
        end if
        coeff_f = (coeff_d/coeff_a*(coeff_b*v1+coeff_c)-coeff_e)/(v2-v1)
        coeff_g = (-coeff_d/coeff_a*(coeff_b*v2+coeff_c)+coeff_e)/(v2-v1)
        n = 0
        rho = sqrt(ky**2+kx**2)*vy
        call bjndd (n, rho, bj, dj, fj)
        tmp_res1 = 2.*pi*zi/coeff_a*coeff_f
        tmp_res2 = 2.*pi*zi/coeff_a*coeff_g

        res1 =  bj**2*sqrt(1.0/pi/2.0)*vy*exp(-0.5*vy**2-&
                0.5*v1**2)*tmp_res1
        res2 =  bj**2*sqrt(1.0/pi/2.0)*vy*exp(-0.5*vy**2-&
                0.5*v2**2)*tmp_res2

        if (s==1) then
                if (aimag(v1)>0.1) then
                        vzprime = vz
                        tmp_res1 = 0.
                        tmp_res2 = 0.
                else if (aimag(v1)>-0.1) then
                        vzprime = vz-zi*0.2
                        tmp_res1 = 0.
                        tmp_res2 = 2.*pi*zi/coeff_a*coeff_g
                else
                        vzprime = vz
                        tmp_res1 = 2.*pi*zi/coeff_a*coeff_f
                        tmp_res2 = 2.*pi*zi/coeff_a*coeff_g
                end if

                omega_star_n = fprime*ky
                omega_star_t = 0.5*(-3.0+vy**2+&
                        vzprime**2)*tprime*ky
                omega_star_d = (1.0-omd_kx)*(0.5*vy**2+&
                        vzprime**2)*omd*ky/A +&
                        omd_kx*(0.5*vy**2+vzprime**2)*&
                        omd*kx/A
                omega_parallel = kz*vzprime

                int_tmp = bj**2*sqrt(1.0/pi/2.0)*vy*exp(-0.5*vy**2-&
                          0.5*vzprime**2)*&
                          (omega-omega_star_n-omega_star_t)/&
                          (omega-omega_parallel-omega_star_d)

                res1 =  bj**2*sqrt(1.0/pi/2.0)*vy*exp(-0.5*vy**2-&
                        0.5*v1**2)*tmp_res1
                res2 =  bj**2*sqrt(1.0/pi/2.0)*vy*exp(-0.5*vy**2-&
                        0.5*v2**2)*tmp_res2
                write(out_unit,'(14e12.4)') kx,ky,kz,omega,vz,v1,v2,tmp_res1,tmp_res2

        else if (s==2) then
                vzprime = vz+zi*3.
                omega_star_n = fprime*ky
                omega_star_t = 0.5*(-3.0+vy**2+vzprime**2)*tprime*ky
                omega_star_d = (1.0-omd_kx)*(0.5*vy**2+&
                        vzprime**2)*omd*ky/A +&
                        omd_kx*(0.5*vy**2+vzprime**2)*&
                        omd*kx/A
                omega_parallel = kz*vzprime
                rho = sqrt(ky**2+kx**2)*vy
                call bjndd (n, rho, bj, dj, fj)
                int_tmp = bj**2*sqrt(1.0/pi/2.0)*vy*exp(-0.5*vy**2-&
                  0.5*vzprime**2)*&
                  (omega-omega_star_n-omega_star_t)/&
                  (omega-omega_parallel-omega_star_d)

        else if (s==3) then
                vzprime = vz-zi*3.
                omega_star_n = fprime*ky
                omega_star_t = 0.5*(-3.0+vy**2+vzprime**2)*tprime*ky
                omega_star_d = (1.0-omd_kx)*(0.5*vy**2+&
                        vzprime**2)*omd*ky/A +&
                        omd_kx*(0.5*vy**2+vzprime**2)*&
                        omd*kx/A
                omega_parallel = kz*vzprime
                rho = sqrt(ky**2+kx**2)*vy
                call bjndd (n, rho, bj, dj, fj)
                int_tmp = bj**2*sqrt(1.0/pi/2.0)*vy*exp(-0.5*vy**2-&
                  0.5*vzprime**2)*&
                  (omega-omega_star_n-omega_star_t)/&
                  (omega-omega_parallel-omega_star_d)

        end if
  end subroutine

  subroutine read_input_file

        use file_utils, only: input_unit_exist, input_unit

        implicit none

        integer :: in_file
        logical :: exist

        namelist / parameters / nstep,nkz,dkz,kz_start, &
                                nky,dky,ky_start, &
                                vy_1,vz_1,vz_2,om1_re,om1_im, &
                                A,vi,tol,fr, &
                                nkx, dkx, kx_start, &
                                tprime, fprime, omd, omd_kx
        
        nstep = 10

        nkz = 1
        dkz = 0.1
        kz_start = 0.0

        nky = 1
        dky = 0.1
        ky_start = 0.0

        vy_1 = 0.0
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
                kx_grid(ikx) = -ikx*dkx+kx_start
        end do

end subroutine

end program