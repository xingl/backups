program kz_avg
 
!        use file_utils
        implicit none

        real (kind=8), parameter :: pi=3.141592653589793
        complex (kind=8), parameter :: zi=(0.0,1.0)
        integer (kind=8) :: nkz,ndelkz
        real (kind=8) :: dkz,ddelkz
        integer (kind=8) :: ikz,idelkz
        real (kind=8) :: delkz_start
        real (kind=8), dimension(:), allocatable :: kz_grid_test,delkz_grid
        real (kind=8) :: kz,delta_kz
        real (kind=8) :: denom
        integer :: in_unit=201, out_unit = 101 
        integer (kind=8) :: n_length, i
        real (kind=8) :: theta_i, tprime_i, fprime_i, kv_i, kx_i, ky_i, kz_i
        real (kind=8) :: freq_i, gamma_i
        integer(kind=8), dimension(:), allocatable :: iky_this, ikz_this
        complex (kind=8), dimension(:,:), allocatable :: omega_array
        real (kind=8), dimension(:), allocatable :: ky_base,kz_base, kx_base
        real (kind=8) :: ky_start,kz_start
        integer (kind=8) :: iky, nky
        real (kind=8) :: ky,dky,kx
        real (kind=8) :: gam, gam_w, gam_avg
        real (kind=8), dimension(:), allocatable :: Dmix_delkz, Dmix_ky
        integer(kind=8), dimension(:), allocatable :: idelkz_max
        real (kind=8) :: Dmix_tmp

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

        ndelkz = 10
        ddelkz = 0.1
        delkz_start = 0.

        allocate(delkz_grid(ndelkz))
        do idelkz = 1, ndelkz
             delkz_grid(idelkz) = idelkz*ddelkz+delkz_start
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

!        call open_output_file(out_unit,'.Dmix')
!        write(out_unit,'(3a12)') 'ky','delkz','Dmixing'

        allocate(Dmix_delkz(ndelkz),Dmix_ky(nky),idelkz_max(1))
        
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
                 gam_w = gam_w + gam*exp(-kz**2/delta_kz**2)*dkz
                 denom = denom + exp(-kz**2/delta_kz**2)*dkz
                 !print*, 'kz, gam, gam_w, denom'
                 !print*, kz, gam, gam_w, denom
              end do
              gam_w = 2./delta_kz/sqrt(pi)*gam_w
              denom = 2./delta_kz/sqrt(pi)*denom
              gam_avg = gam_w/denom
              Dmix_delkz(idelkz) = gam_avg/(ky**2+kx**2)
              !write(out_unit,'(3e12.4)') ky,delta_kz,Dmix_delkz(idelkz)
              !print*, Dmix_delkz
           end do
           idelkz_max = maxloc(Dmix_delkz)
           !print*, iky,idelkz_max(1)
           !print*,Dmix_delkz(idelkz_max(1))
           Dmix_tmp = Dmix_delkz(idelkz_max(1))
           !Dmix_ky(iky) = Dmix_tmp
        end do

!        call close_output_file(out_unit)

end program
