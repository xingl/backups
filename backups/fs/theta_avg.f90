program theta_avg
        
        implicit none

        real (kind=8), parameter :: pi=3.141592653589793
        complex (kind=8), parameter :: zi=(0.0,1.0)
        integer (kind=8) :: ndelkz
        real (kind=8) :: ddelkz
        integer (kind=8) :: idelkz
        real (kind=8) :: delkz_start
        real (kind=8), dimension(:), allocatable :: delkz_grid
        real (kind=8), dimension(:), allocatable :: deltheta_grid
        real (kind=8) :: delta_kz,delta_theta
        real (kind=8) :: denom
        integer :: in_unit=201, out_kz_unit = 101, gamE_unit=202 
        integer :: out_ky_unit = 101
        integer (kind=8) :: n_length, n_skip, i, itheta
        real (kind=8) :: theta_i, tprime_i, fprime_i, kv_i, kx_i, ky_i, kz_i
        real (kind=8) :: freq_i, gamma_i
        real (kind=8) :: theta, theta_center, theta_left, theta_right
        integer(kind=8), dimension(:), allocatable :: itheta_center
        integer(kind=8), dimension(:), allocatable :: itheta_left
        integer(kind=8), dimension(:), allocatable :: itheta_right
        real (kind=8) :: gamE_w, gamE, gamE_avg
        real (kind=8) :: R, q, alpha
        real (kind=8), dimension(:),allocatable :: B_tor_grid,B_tot_grid
        real (kind=8), dimension(:),allocatable :: theta_grid, theta_grid_ext
        real (kind=8), dimension(:),allocatable :: gamE_grid,B_pol_grid

        alpha = 1.
        !R = 6.36
        R = 1.
        q = 1.07

        n_length=1708
        allocate(gamE_grid(n_length),B_pol_grid(n_length),B_tor_grid(n_length))
        allocate(B_tot_grid(n_length),theta_grid(n_length))

        open(unit=gamE_unit,file='gamE_iterbase.dat',action='read')
        do i = 1, n_length
           read(gamE_unit,'(5f12.4)') gamE_grid(i), B_pol_grid(i),&
                B_tor_grid(i), B_tot_grid(i), theta_grid(i)
           print*, i
        end do
        close(unit=gamE_unit)
        allocate(theta_grid_ext(3*n_length-2))
        do i = n_length, 2*n_length-1
           theta_grid_ext(i) = theta_grid(i-n_length+1)
        end do
        do i = 1, n_length-1
           theta_grid_ext(i) = -theta_grid(n_length-i+1)
        end do
        do i = 2*n_length, 3*n_length-2
           theta_grid_ext(i) = 2*pi-theta_grid(4*n_length-i-2)
        end do

        ndelkz = 20
        ddelkz = 0.1
        delkz_start = 0.

        allocate(delkz_grid(ndelkz),deltheta_grid(ndelkz))
        do idelkz = 1, ndelkz
             delkz_grid(idelkz) = idelkz*ddelkz+delkz_start
             deltheta_grid(idelkz) = B_pol_grid(1)/B_tot_grid(1)*&
                2./delkz_grid(idelkz)/R/q
        end do

        print*, 'deltheta_grid'
        print*, deltheta_grid
        
        allocate(itheta_center(1),itheta_left(1),itheta_right(1))

        open(unit=out_kz_unit,file='gamE_avg.dat',action='write')
        write(out_kz_unit,'(4a12)') 'theta','delkz','deltheta','gamE_avg'

        n_skip = 100
        do i = 1, n_length, n_skip
           theta_center = theta_grid(i)
           print*, 'index along fs'
           print*, i
           do idelkz = 1, ndelkz
              delta_theta = deltheta_grid(idelkz)
              print*, 'delta_theta'
              print*, delta_theta
              theta_left = theta_center - 3.*delta_theta
              theta_right = theta_center + 3.*delta_theta
              itheta_center = minloc(abs(theta_grid_ext-theta_center))
              itheta_left = minloc(abs(theta_grid_ext-theta_left))
              itheta_right = minloc(abs(theta_grid_ext-theta_right))
              print*, 'index of theta'
              print*, itheta_center, itheta_left, itheta_right
              print*, 'bounds of theta integral'
              print*, theta_grid_ext(itheta_center),&
                theta_grid_ext(itheta_left),&
                theta_grid_ext(itheta_right)
              print*, theta_center, theta_left, theta_right
              gamE_w = 0.
              denom = 0.
              do itheta = itheta_left(1), itheta_right(1)
                 if (itheta>(n_length-1).and.itheta<2*n_length) then
                    print*, 'this theta'
                    print*,theta_grid_ext(itheta)
                    theta = theta_grid(itheta-n_length+1)
                    print*, theta
                    gamE = gamE_grid(itheta-n_length+1)
                 else if (itheta<n_length) then
                    print*, 'this theta'
                    print*,theta_grid_ext(itheta)
                    theta = theta_grid(n_length-itheta+1)
                    print*, theta
                    gamE = gamE_grid(n_length-itheta+1)
                 else
                    print*, 'this theta'
                    print*,theta_grid_ext(itheta)
                    theta = theta_grid(3*n_length-itheta)
                    print*, theta
                    gamE = gamE_grid(3*n_length-itheta)
                 end if
                 gamE_w = gamE_w + gamE*exp(-(theta-theta_center)**2/&
                        delta_theta**2)
                 denom = denom + exp(-(theta-theta_center)**2/&
                        delta_theta**2)
                 print*, 'gamE_w, denom'
                 print*, gamE_w, denom
              end do
              gamE_w = 1./delta_theta/sqrt(pi)*gamE_w
              denom  = 1./delta_theta/sqrt(pi)*denom
              gamE_avg = gamE_w/denom
              write(out_kz_unit,'(4e12.4)') theta_center, &
                delkz_grid(idelkz),deltheta_grid(idelkz),gamE_avg
           end do
           write(out_kz_unit,*)
        end do
        close(unit=out_kz_unit)
end program
