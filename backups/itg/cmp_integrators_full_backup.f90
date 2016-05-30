program cmp_integrators_vy

        use file_utils

        implicit none

        real (kind=8), parameter :: pi=3.141592653589793
        complex (kind=8), parameter :: zi=(0.0,1.0)

        integer :: nvy, nvz
        integer :: iy, iz
        real (kind=8) :: lvy, lvz
        real (kind=8) :: ky, kz
        complex (kind=8) :: omega, integral1, integral2 
        real (kind=8) :: A, theta, fprime, tprime, omd
        real (kind=8) :: vi, tol, fr
        integer :: out_unit=101, bf_unit= 102, sa_unit=103

        real (kind=8), dimension(:), allocatable :: vy_grid, vz_grid

        real (kind=8) :: om1_re, om1_im
        complex (kind=8) :: seed1
        
        call init_file_utils
        call read_input_file

        allocate(vy_grid(nvy), vz_grid(nvz))
        do iy = 1, nvy
                vy_grid(iy) = iy*lvy/nvy
        end do
        do iz = 1, nvz
                vz_grid(iz) = -lvz/2.0 + (iz-1.0)*lvz/nvz
        end do

        seed1 = om1_re + om1_im*zi
        call open_output_file(bf_unit,'.bf')
        call open_output_file(sa_unit,'.sa')
        call open_output_file(out_unit,'.dat')
        write (out_unit,'(8a16)') 'omega','','int1','','int2','','diff',''
        call integrator(ky,kz,seed1,integral1)
        call simpson_adaptive(seed1,integral2,tol)
        write (out_unit, '(8e16.6)') seed1, integral1, integral2, &
                                        integral1-integral2
        call close_output_file (out_unit)
        call close_output_file (bf_unit)
        call close_output_file (sa_unit)
        call finish_file_utils

contains

  subroutine simpson_adaptive(omega,integral,tol)
        
        implicit none

        real(kind=8), intent(in) :: tol
        real(kind=8) :: a,b
        complex(kind=8), intent(in) :: omega
        complex (kind=8), intent(out) :: integral
        real(kind=8) :: x_left, x_center, x_right, dx
        complex (kind=8) :: f_left, f_center, f_right
        complex (kind=8) :: i_trapezoid, i_simpson, int_trap
        
        a = vy_grid(1)
        b = vy_grid(nvy)
        dx = 0.1
        x_left = a
        x_right = x_left + dx
        call vz_integral(x_left,f_left,tol,omega)
        integral = 0.0
        int_trap = 0.0

        do while (x_left<b)
                x_center = 0.5*(x_left+x_right)
                call vz_integral(x_center,f_center,tol,omega)
                call vz_integral(x_right,f_right,tol,omega)
                i_trapezoid = fr*(f_left+2.0*f_center+f_right)*0.5*dx
                i_simpson = 1.0/6.0*dx*(f_left+4.0*f_center+f_right)
                do while (abs(i_trapezoid-i_simpson)>tol)
                        dx=0.5*dx*(abs(i_trapezoid-i_simpson)/tol)**(-1.0/3.0)
                        if (dx>0.2) dx=0.2
                        x_right = x_left + dx
                        if (x_right>b) then
                                x_right = b
                                dx = x_right - x_left
                        end if
                        x_center = 0.5*(x_left+x_right)
                        call vz_integral(x_center,f_center,tol,omega)
                        call vz_integral(x_right,f_right,tol,omega)
                        i_trapezoid = 0.5*(f_left+2.0*f_center+f_right)*0.5*dx
                        i_simpson = 1.0/6.0*dx*(f_left+4.0*f_center+f_right)
                end do
                integral = integral + i_simpson
                int_trap = int_trap + i_trapezoid
                write(sa_unit,'(7e16.6)') x_right, integral, int_trap,&
                                        integral-int_trap
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

        integral = 1.0 + theta - theta*integral

  end subroutine

  subroutine vz_integral(vy,integral,tol,omega)
        
        implicit none

        real(kind=8), intent(in) :: tol, vy
        complex(kind=8), intent(in) :: omega
        real(kind=8) :: a, b
        complex (kind=8), intent(out) :: integral
        real(kind=8) :: x_left, x_center, x_right, dx
        complex (kind=8) :: f_left, f_center, f_right
        complex (kind=8) :: i_trapezoid, i_simpson
        
        a = vz_grid(1)
        b = vz_grid(nvz)
        dx = 0.1
        x_left = a
        x_right = x_left + dx
        f_left = integrand(x_left,vy,omega) 
        integral = 0.0

        do while (x_left<b)
                x_center = 0.5*(x_left+x_right)
                f_center = integrand(x_center,vy,omega)
                f_right = integrand(x_right,vy,omega)
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
                        f_center = integrand(x_center,vy,omega)
                        f_right = integrand(x_right,vy,omega)
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

  end subroutine

  function integrand(vz,vy,omega)

        real(kind=8) :: vz,vy
        complex(kind=8) :: omega,integrand,int_tmp
        complex(kind=8) :: omega_star_n, omega_star_t, omega_star_d
        complex(kind=8) :: omega_parallel
        real(kind=8) :: rho,bj,dj,fj
        integer :: n 

        omega_star_n = fprime*ky
        int_tmp = 0.0
        rho = ky*vy
        n = 0
        call bjndd (n, rho, bj, dj, fj)
        omega_star_t = 0.5*(-3.0+vy**2+&
                        (vz**2-zi*2.0*vz*vi-vi**2))*tprime*ky        
        omega_star_d = (0.5*vy**2+&
                        vz**2-zi*2.0*vz*vi-vi**2)*omd*ky/A
        omega_parallel = kz*(vz-zi*vi)

        int_tmp = int_tmp + bj**2*&
                        sqrt(1.0/pi/2.0)*vy*exp(-0.5*vy**2-&
                        0.5*(vz**2-zi*2.0*vz*vi-&
                        vi**2))*(omega-omega_star_n-omega_star_t)/&
                        (omega-omega_parallel+omega_star_d)
        integrand = int_tmp

  end function

  subroutine integrator ( ky, kz, omega, intgral)

        implicit none

        complex (kind=8), intent(in) :: omega
        real (kind=8), intent(in) :: ky, kz
        complex (kind=8), intent(out) :: intgral
        complex (kind=8) :: intgral_tmp
	real (kind=8) :: rho, bj, dj, fj
        integer :: n
        
        intgral_tmp = 0.0
        n = 0
        do iy = 1, nvy
                do iz = 1, nvz
                        rho = ky*vy_grid(iy)
                        call bjndd (n, rho, bj, dj, fj)
                        intgral_tmp = intgral_tmp + bj**2*&
                        sqrt(1.0/pi/2.0)*vy_grid(iy)*exp(-0.5*vy_grid(iy)**2-&
                        0.5*(vz_grid(iz)**2-zi*2.0*vz_grid(iz)*vi-&
                        vi**2))*(omega-fprime*ky-&
                        0.5*(-3.0+vy_grid(iy)**2+&
                        (vz_grid(iz)**2-zi*2.0*vz_grid(iz)*vi-vi**2))&
                        *tprime*ky)/&
                        (omega-kz*(vz_grid(iz)-zi*vi)+&
                        (0.5*vy_grid(iy)**2+&
                        vz_grid(iz)**2-zi*2.0*vz_grid(iz)*vi- &
                        vi**2)*omd*ky/A)*lvy*lvz/nvy/nvz
                end do  
                write(bf_unit,'(3e16.6)') vy_grid(iy), intgral_tmp
        end do
        
        intgral = 1.0+theta-theta*intgral_tmp

  end subroutine

  subroutine read_input_file

        use file_utils, only: input_unit_exist, input_unit

        implicit none

        integer :: in_file
        logical :: exist

        namelist / parameters / nvy,nvz,lvy,lvz,A,theta,fr, &
                                vi,om1_re,om1_im,kz,ky,fprime,tprime,omd,tol
        
        nvy = 1000
        nvz = 2000
        lvy = 4.0
        lvz = 8.0
        fprime = 0.5
        A = 3.0
        theta = 1.0
        vi  = 1.0
        om1_re = 1.0
        om1_im = 0.1
        kz = 0.1
        ky = 0.1
        tprime = 1.0
        omd = 1.0
        tol = 1.0D-05
        fr = 0.5

    in_file = input_unit_exist ("parameters", exist)
    if (exist) read (unit=input_unit("parameters"), nml=parameters)

  end subroutine read_input_file
end program
