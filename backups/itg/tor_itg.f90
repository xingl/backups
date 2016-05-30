program itg_calculation

        use file_utils

        implicit none

        real (kind=8), parameter :: pi=3.141592653589793
        complex (kind=8), parameter :: zi=(0.0,1.0)

        integer :: nvy, nvz, nkz, nky, neta
        integer :: iy, iz, ikz, iky, ieta
        real (kind=8) :: lvy, lvz, dkz, dky, deta
        complex (kind=8) :: omega 
        real (kind=8) :: omt_e, a_r, theta, eta, omn, omd
        real (kind=8) :: vi
        integer :: nstep
        integer :: istep
        integer :: out_unit=101

        real (kind=8), dimension(:), allocatable :: vy_grid, vz_grid
        real (kind=8), dimension(:), allocatable :: ky_grid, kz_grid
        real (kind=8), dimension(:), allocatable :: eta_grid

        complex (kind=8) :: root
        real (kind=8) :: om1_re, om1_im, om2_re, om2_im
        complex (kind=8) :: seed1, seed2

        complex(kind=8) :: intgral

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

        do iky = 1, nky
                ky_grid(iky) = (nky+1-iky)*dky
        end do

        allocate(eta_grid(neta))

        do ieta = 1, neta
                eta_grid(ieta) = (neta+1-ieta)*deta
        end do

!        open (unit = out_unit, file = "slab_itg_rf08a.dat")
        call open_output_file(out_unit,'.dat')
        write (out_unit,'(7a12)') "eta","omt","omn","ky","kz","frequency","growth rate"
        do ieta = 1, neta
          omt_e = omn*eta_grid(ieta)
          print*, "eta, omt, omn"
          print*, eta_grid(ieta), omt_e, omn
          do iky = 1, nky
            do ikz = 1, nkz
                print*, "ky, kz"
                print*, ky_grid(iky), kz_grid(ikz)
                if (ieta==1 .and. iky==1 .and. ikz==1) then
                        print*, "Two initial omega values:"
                        seed1 = om1_re + om1_im*zi
                        seed2 = om2_re + om2_im*zi
                        print*, seed1
                        print*, seed2 
                        call rootfinder(omt_e,ky_grid(iky),kz_grid(ikz),seed1,seed2,root)
                        print*, "Root is found:"
                        print*, root    
                        write (out_unit, '(7e12.4)') eta_grid(ieta),omt_e,omn,ky_grid(iky), kz_grid(ikz),root
                else
                        print*, "Two initial omega values:"
                        seed1 = root 
                        seed2 = root*(0.9) 
                        print*, seed1
                        print*, seed2 
                        call rootfinder(omt_e,ky_grid(iky),kz_grid(ikz),seed1,seed2,root)
                        print*, "Root is found:"
                        print*, root
                        write (out_unit, '(7e12.4)') eta_grid(ieta),omt_e,omn,ky_grid(iky),kz_grid(ikz),root
                end if
            end do
          end do
        end do
        call close_output_file (out_unit)
        call finish_file_utils

contains

 subroutine rootfinder(omt,ky,kz,a,b,c)

        implicit none

        complex (kind=8), intent(in) :: a,b
        real (kind=8), intent(in) :: omt,ky,kz
        complex (kind=8), intent(out) :: c
        complex (kind=8), dimension(:), allocatable :: x, f
        complex (kind=8) :: x_tmp, f_tmp

        allocate(x(nstep),f(nstep))
        x = 0.0
        f = 0.0
        x(1) = a
        call integrator(omt, ky, kz, x(1),f(1))
        x(2) = b
        call integrator(omt, ky, kz, x(2),f(2))

        do istep = 3, nstep
                x(istep) = x(istep-1) - f(istep-1) * &
                           (x(istep-1)-x(istep-2))/(f(istep-1)-f(istep-2))
                x_tmp = x(istep)
                call integrator(omt, ky, kz, x_tmp,f_tmp)
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

  subroutine integrator (omt, ky, kz, omega, intgral)

        implicit none

        complex (kind=8), intent(in) :: omega
        real (kind=8), intent(in) :: ky, kz, omt
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
                        vi**2))*(omega-omn*ky-&
                        0.5*(-3.0+vy_grid(iy)**2+&
                        (vz_grid(iz)**2-zi*2.0*vz_grid(iz)*vi-vi**2))&
                        *omt*ky)/&
                        (omega-kz*(vz_grid(iz)-zi*vi)+&
                        (0.5*vy_grid(iy)**2+&
                        vz_grid(iz)**2-zi*2.0*vz_grid(iz)*vi- &
                        vi**2)*omd*ky/a_r)*lvy*lvz/nvy/nvz
                end do  
        end do
        
        intgral = 1+theta-theta*intgral_tmp

  end subroutine



  subroutine read_input_file

        use file_utils, only: input_unit_exist, input_unit

        implicit none

        integer :: in_file
        logical :: exist

        namelist / parameters / nstep,nvy,nvz,lvy,lvz,omn,a_r,theta,&
                                vi,om1_re,om1_im,om2_re,om2_im,nkz,dkz, &
                                nky,dky,neta,deta,omd
        
        nstep = 100
        nvy = 1000
        nvz = 2000
        lvy = 4.0
        lvz = 8.0
        omn = 0.5
        omd = 1.0
        a_r = 3.0
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
        neta = 4
        deta = 0.5

    in_file = input_unit_exist ("parameters", exist)
    if (exist) read (unit=input_unit("parameters"), nml=parameters)

  end subroutine read_input_file
end program
