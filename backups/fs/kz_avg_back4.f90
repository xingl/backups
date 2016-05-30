program kz_avg
 
!        use file_utils
        implicit none

        real (kind=8), parameter :: pi=3.141592653589793
        complex (kind=8), parameter :: zi=(0.0,1.0)
        integer (kind=8) :: nkz,ndelkz
        real (kind=8) :: dkz,ddelkz
        integer (kind=8) :: ikz,idelkz
        real (kind=8) :: kz_start,delkz_start
        real (kind=8), dimension(:), allocatable :: kz_grid_test,delkz_grid
        real (kind=8) :: kz,delta_kz
        real (kind=8) :: denom
        integer :: in_unit=201, out_unit = 101 
        integer (kind=8) :: n_length, i
        real (kind=8), dimension(:), allocatable :: theta_fs_grid
        real (kind=8), dimension(:), allocatable :: tprime_grid, fprime_grid
        real (kind=8), dimension(:), allocatable :: kv_grid,kx_grid
        real (kind=8), dimension(:), allocatable :: ky_grid,kz_grid     
        real (kind=8), dimension(:), allocatable :: freq_grid, gamma_grid
        integer (kind=8) :: iky, nky
        real (kind=8) :: ky,kx
        real (kind=8) :: gam, gam_w, gam_avg
        real (kind=8), dimension(:), allocatable :: Dmix_delkz, Dmix_ky
        integer(kind=8), dimension(:), allocatable :: idelkz_max
        real (kind=8) :: Dmix_tmp

        n_length=929
        allocate(fprime_grid(n_length),tprime_grid(n_length),theta_fs_grid(n_length))
        allocate(kv_grid(n_length),kx_grid(n_length),ky_grid(n_length),kz_grid(n_length))
        allocate(freq_grid(n_length),gamma_grid(n_length))

        open(unit=in_unit,file='fs02f.datgam',action='read')
        do i = 1,n_length
           read(in_unit,'(9f12.4)') theta_fs_grid(i),tprime_grid(i),&
                fprime_grid(i), kv_grid(i),kx_grid(i),ky_grid(i), &
                kz_grid(i),freq_grid(i),gamma_grid(i)
        !   if(i==101) print*,theta_fs_grid(i),tprime_grid(i),&
        !        fprime_grid(i), kv_grid(i),kx_grid(i),ky_grid(i), &
        !        kz_grid(i),freq_grid(i),gamma_grid(i)
           print*, i
        end do
        close(unit=in_unit)
        
        nkz = 30

        ndelkz = 10
        ddelkz = 0.003
        delkz_start = 0.

       ! allocate(kz_grid_test(nkz))
        !dkz = 0.001
        !kz_start = -0.001
        !do ikz = 1, nkz
        !        kz_grid_test(ikz) = ikz*dkz+kz_start
        !end do

        allocate(delkz_grid(ndelkz))
        do idelkz = 1, ndelkz
             delkz_grid(idelkz) = idelkz*ddelkz+delkz_start
        end do

!        call open_output_file(out_unit,'.Dmix')
!        write(out_unit,'(3a12)') 'ky','delkz','Dmixing'

        allocate(Dmix_delkz(ndelkz),Dmix_ky(nky),idelkz_max(1))
        nky = 30
        nkz = 30
        dkz = 0.001
        do iky = 1, nky
           ky = ky_grid((iky-1)*(nkz+1)+1)
           kx = kx_grid((iky-1)*(nkz+1)+1)
           !print*, 'ky,kx'
           !print*, ky,kx
           do idelkz = 1, ndelkz
              delta_kz = delkz_grid(idelkz)
              !print*, 'delta_kz'
              !print*, delta_kz
              denom = 0.
              gam_w = 0.
              do ikz = 1, nkz
                 kz = kz_grid((iky-1)*(nkz+1)+ikz)   
                 gam = gamma_grid((iky-1)*(nkz+1)+ikz)   
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
           print*, iky,idelkz_max(1)
           print*,Dmix_delkz(idelkz_max(1))
           Dmix_tmp = Dmix_delkz(idelkz_max(1))
           !Dmix_ky(iky) = Dmix_tmp
        end do

!        call close_output_file(out_unit)

end program
