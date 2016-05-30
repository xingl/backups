program itg_calculation

        use file_utils

        implicit none

        real (kind=8), parameter :: pi=3.141592653589793
        complex (kind=8), parameter :: zi=(0.0,1.0)

        integer :: nvy, nvz, nkz, nky 
        integer :: nfprime, ntprime, nomd
        integer :: iy, iz, ikz, iky
        integer :: ifprime, itprime, iomd
        real (kind=8) :: lvy, lvz, dkz, dky
        real (kind=8) :: dfprime, dtprime, domd
        complex (kind=8) :: omega 
        real (kind=8) :: A, theta, fprime, tprime, omd
        real (kind=8) :: vi, fprime_start, ky_start
        integer :: nstep
        integer :: istep
        integer :: out_unit=101, int_unit=102

        real (kind=8), dimension(:), allocatable :: vy_grid, vz_grid
        real (kind=8), dimension(:), allocatable :: ky_grid, kz_grid
        real (kind=8), dimension(:), allocatable :: fprime_grid, tprime_grid, omd_grid

        complex (kind=8) :: root
        real (kind=8) :: om1_re, om1_im, om2_re, om2_im
        complex (kind=8) :: seed1, seed2
        
        character(len=2) :: ky_direction, fprime_direction
        
        call init_file_utils
        call read_input_file

        allocate(vy_grid(nvy), vz_grid(nvz))
        do iy = 1, nvy
                vy_grid(iy) = iy*lvy/nvy
        end do
        do iz = 1, nvz
                vz_grid(iz) = -lvz/2.0 + (iz-1.0)*lvz/nvz
        end do

        allocate(ky_grid(nky), kz_grid(nkz))
        do ikz = 1, nkz
                kz_grid(ikz) = (nkz+1-ikz)*dkz
        end do
        if (ky_direction=='bw') then
                do iky = 1, nky
                        ky_grid(iky) = (nky+1-iky)*dky+ky_start
                end do
        else
                do iky = 1, nky
                        ky_grid(iky) = iky*dky+ky_start
                end do
                
        end if

        allocate(fprime_grid(nfprime), tprime_grid(ntprime), omd_grid(nomd))
        if (fprime_direction=='bw') then
                do ifprime = 1, nfprime
                        fprime_grid(ifprime)=(nfprime+1-ifprime)*dfprime+fprime_start
                end do 
        else
                do ifprime = 1, nfprime
                        fprime_grid(ifprime) = ifprime*dfprime+fprime_start
                end do 
        end if
                
        do itprime = 1, ntprime
                tprime_grid(itprime) = (ntprime+1-itprime)*dtprime
        end do
        do iomd = 1, nomd
                omd_grid(iomd) = (nomd+1-iomd)*domd
        end do

        call open_output_file(out_unit,'.dat')
        write (out_unit,'(7a12)') "tprime","fprime","omd","ky","kz","frequency","growth rate"
        do iomd = 1, nomd
           omd = omd_grid(iomd)
           do ifprime = 1, nfprime
              fprime = fprime_grid(ifprime)
              do itprime = 1, ntprime
                 tprime = tprime_grid(itprime)
                 print*, "tprime, fprime, omd"
                 print*, tprime, fprime, omd
                 do iky = 1, nky
                    do ikz = 1, nkz
                       print*, "ky, kz"
                       print*, ky_grid(iky), kz_grid(ikz)
                       if (ifprime ==1 .and. itprime==1 .and. iky==1 .and. ikz==1) then
                            print*, "Two initial omega values:"
                            seed1 = om1_re + om1_im*zi
                            seed2 = om2_re + om2_im*zi
                            print*, seed1, seed2
                            call rootfinder(tprime,fprime,omd,ky_grid(iky),kz_grid(ikz),seed1,seed2,root)
                            print*, "Root is found:"
                            print*, root    
                            write (out_unit, '(7e12.4)') tprime,fprime,omd/A,ky_grid(iky),kz_grid(ikz),root
                       else
                            print*, "Two initial omega values:"
                            seed1 = root 
                            seed2 = real(root)*0.9+aimag(root) 
                            print*, seed1, seed2
                            call rootfinder(tprime,fprime,omd,ky_grid(iky),kz_grid(ikz),seed1,seed2,root)
                            print*, "Root is found:"
                            print*, root
                            write (out_unit, '(7e12.4)') tprime,fprime,omd/A,ky_grid(iky),kz_grid(ikz),root
                       end if
                    end do
                 end do
              end do
           end do
        end do
        call close_output_file (out_unit)
        call finish_file_utils

contains

 subroutine rootfinder(tprime,fprime,omd,ky,kz,a,b,c)

        implicit none

        complex (kind=8), intent(in) :: a,b
        real (kind=8), intent(in) :: tprime,fprime,omd,ky,kz
        complex (kind=8), intent(out) :: c
        complex (kind=8), dimension(:), allocatable :: x, f
        complex (kind=8) :: x_tmp, f_tmp

        allocate(x(nstep),f(nstep))
        x = 0.0
        f = 0.0
        x(1) = a
        call integrator(tprime,fprime,omd, ky, kz, x(1),f(1))
        x(2) = b
        call integrator(tprime,fprime,omd, ky, kz, x(2),f(2))

        do istep = 3, nstep
                x(istep) = x(istep-1) - f(istep-1) * &
                           (x(istep-1)-x(istep-2))/(f(istep-1)-f(istep-2))
                x_tmp = x(istep)
                call integrator(tprime,fprime,omd, ky, kz, x_tmp,f_tmp)
                f(istep) = f_tmp
                print*, istep
                print*, x(istep), f(istep)
                if (abs(f(istep))<1.0E-04) then
                        c = x(istep)
                        exit
                end if
        end do
        deallocate(x,f)

  end subroutine

  subroutine integrator (tprime,fprime,omd, ky, kz, omega, intgral)

        implicit none

        complex (kind=8), intent(in) :: omega
        real (kind=8), intent(in) :: ky, kz, tprime,fprime,omd
        complex (kind=8), intent(out) :: intgral
        complex (kind=8) :: intgral_tmp
        integer :: n
	real (kind=8) :: rho, bj, dj, fj
        
        n = 0
        intgral_tmp = 0.0

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
        end do
        
        intgral = 1+theta-theta*intgral_tmp

  end subroutine

  subroutine read_input_file

        use file_utils, only: input_unit_exist, input_unit

        implicit none

        integer :: in_file
        logical :: exist

        namelist / parameters / nstep,nvy,nvz,lvy,lvz,A,theta,&
                                vi,om1_re,om1_im,om2_re,om2_im,nkz,dkz, &
                                nky,dky,nfprime,dfprime,ntprime,dtprime, &
                                nomd,domd,ky_direction,fprime_direction,&
                                ky_start, fprime_start 
        
        nstep = 100
        nvy = 1000
        nvz = 2000
        lvy = 4.0
        lvz = 8.0
        dfprime = 0.5
        nfprime = 1.0
        A = 3.0
        theta = 1.0
        vi  = 1.0
        om1_re = 1.0
        om1_im = 0.1
        om2_re = 1.0
        om2_im = 0.1
        nkz = 1
        dkz = 0.1
        nky = 1
        dky = 0.1
        ntprime = 1
        dtprime = 1.0
        nomd = 1
        domd = 1.0
        ky_direction = 'fw'
        fprime_direction = 'bw'
        ky_start = 0.0
        fprime_start = 0.0

    in_file = input_unit_exist ("parameters", exist)
    if (exist) read (unit=input_unit("parameters"), nml=parameters)

  end subroutine read_input_file
end program
