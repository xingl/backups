program kz_avg
 
        implicit none

        real (kind=8), parameter :: pi=3.141592653589793
        complex (kind=8), parameter :: zi=(0.0,1.0)
        integer (kind=8) :: nkz,ndelkz
        real (kind=8) :: dkz,ddelkz
        integer (kind=8) :: ikz,idelkz
        real (kind=8) :: delkz_start
        real (kind=8), dimension(:), allocatable :: delkz_grid
        real (kind=8), dimension(:), allocatable :: deltheta_grid
        real (kind=8) :: kz,delta_kz,delta_theta
        real (kind=8) :: denom
        integer :: in_unit=201, out_unit = 101, gamE_unit=202 
        integer (kind=8) :: n_length, i, itheta
        real (kind=8) :: theta_i, tprime_i, fprime_i, kv_i, kx_i, ky_i, kz_i
        real (kind=8) :: freq_i, gamma_i
        real (kind=8) :: theta
        integer(kind=8), dimension(:), allocatable :: iky_this, ikz_this
        complex (kind=8), dimension(:,:), allocatable :: omega_array
        real (kind=8), dimension(:), allocatable :: ky_base,kz_base, kx_base
        real (kind=8) :: ky_start,kz_start
        integer (kind=8) :: iky, nky
        real (kind=8) :: ky,dky,kx
        real (kind=8) :: gam, gam_w, gam_avg, gamE_w,gamE, gamE_avg
        real (kind=8), dimension(:,:), allocatable :: Dmix_delkz, Dmix_gamE
        real (kind=8), dimension(:), allocatable :: Dmix_ky
        integer(kind=8), dimension(:), allocatable :: idelkz_max, itheta_ref
        real (kind=8) :: Dmix_tmp
        real (kind=8) :: R, q
        real (kind=8), dimension(:),allocatable :: B_tor_grid,B_tot_grid,theta_grid
        real (kind=8), dimension(:),allocatable :: gamE_grid,B_pol_grid


        R = 6.36
        q = 1.07

        nky = 25
        nkz = 25
        dkz = 0.1
        dky = 0.1
        kz_start = -0.1
        ky_start = -0.1

        allocate(ky_base(nky),kz_base(nkz),kx_base(nky))
        do iky = 1, nky
           ky_base(iky) = iky*dky+ky_start
        end do
        do ikz = 1,nkz
           kz_base(ikz) = ikz*dkz+kz_start
        end do
        kx_base = 0.
        
        allocate(omega_array(nky,nkz))
        omega_array = 0.

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
        

        ndelkz = 10
        ddelkz = 0.1
        delkz_start = 0.

        allocate(delkz_grid(ndelkz),deltheta_grid(ndelkz))
        do idelkz = 1, ndelkz
             delkz_grid(idelkz) = idelkz*ddelkz+delkz_start
             deltheta_grid(idelkz) = B_pol_grid(1)/B_tot_grid(1)*&
                2./delkz_grid(idelkz)/R/q
        end do

        allocate(iky_this(1),ikz_this(1))
        n_length=450
        open(unit=in_unit,file='fs03f.datgam',action='read')
        do i = 1,n_length
           read(in_unit,'(9f12.4)') theta_i,tprime_i,&
                fprime_i, kv_i,kx_i,ky_i, kz_i,freq_i,gamma_i

           iky_this = minloc(abs(ky_base-ky_i))
           ikz_this = minloc(abs(kz_base-kz_i))
           omega_array(iky_this(1),ikz_this(1)) = freq_i + zi*gamma_i
           kx_base(iky_this(1)) = kx_i
           print*, ky_base(iky_this(1)),kz_base(ikz_this(1))
           print*, omega_array(iky_this(1),ikz_this(1))

        end do
        close(unit=in_unit)

        open(unit=out_unit,file='fs03f.dm',action='write')
        write(out_unit,'(4a12)') 'ky','delkz','Dmix','Dmix_gam'

        allocate(Dmix_delkz(nky,ndelkz),Dmix_gamE(nky,ndelkz))
        allocate(Dmix_ky(nky),idelkz_max(1),itheta_ref(1))
        Dmix_delkz = 0.; Dmix_gamE = 0.; Dmix_ky = 0.
        
        do iky = 1, nky
           ky = ky_base(iky)
           kx = kx_base(iky)
           !print*, 'ky,kx'
           !print*, ky,kx
           do idelkz = 1, ndelkz
              delta_kz = delkz_grid(idelkz)
              !print*, 'delta_kz'
              !print*, delta_kz
              denom = 0.
              gam_w = 0.
              do ikz = 1, nkz
                 kz = kz_base(ikz)
                 gam = aimag(omega_array(iky,ikz))
                 if (gam==0.) cycle
                 gam_w = gam_w + gam*exp(-kz**2/delta_kz**2)*dkz
                 denom = denom + exp(-kz**2/delta_kz**2)*dkz
                 !print*, 'kz, gam, gam_w, denom'
                 !print*, kz, gam, gam_w, denom
              end do
              gam_w = 2./delta_kz/sqrt(pi)*gam_w
              denom = 2./delta_kz/sqrt(pi)*denom
              if (denom==0.) cycle
              gam_avg = gam_w/denom
              Dmix_delkz(iky,idelkz) = gam_avg/(ky**2+kx**2)


              delta_theta = deltheta_grid(idelkz)
              itheta_ref = minloc(abs(theta_grid+4.*delta_theta))
              gamE_w = 0.
              denom = 0.
              do itheta = 1, itheta_ref(1)
                 theta = theta_grid(itheta)
                 gamE = gamE_grid(itheta)
                 gamE_w = gamE_w + gamE*exp(-theta**2/&
                        delta_theta**2)*(theta_grid(itheta+1)-theta)
                 denom = denom + exp(-theta**2/&
                        delta_theta**2)*(theta_grid(itheta+1)-theta)
              end do
              gamE_w = 2./delta_theta/sqrt(pi)*gamE_w
              denom  = 2./delta_theta/sqrt(pi)*denom
              gamE_avg = gamE_w/denom

              Dmix_gamE(iky, idelkz) = (gam_avg-gamE_avg)/(ky**2+kx**2)
              print*,' ky,delkz,Dmix_delkz,Dmix_gamE'
              print*, ky,delta_kz,Dmix_delkz(iky,idelkz),Dmix_gamE(iky,idelkz) 
        
              write(out_unit,'(4e12.4)') ky,delta_kz,Dmix_delkz(iky,idelkz), &
                                Dmix_gamE(iky,idelkz)
           end do
           write(out_unit,*)
        end do

        close(unit=out_unit)

end program
