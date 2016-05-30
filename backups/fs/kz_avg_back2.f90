program kz_avg
 
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

        n_length=182700
        allocate(fprime_grid(n_length),tprime_grid(n_length),theta_fs_grid(n_length))
        allocate(kv_grid(n_length),kx_grid(n_length),ky_grid(n_length),kz_grid(n_length))
        allocate(freq_grid(n_length),gamma_grid(n_length))

        open(unit=in_unit,file='fs02f.dat',action='read')
        do i = 1,n_length
           read(in_unit,'(9f12.4)') theta_fs_grid(i),tprime_grid(i),&
                fprime_grid(i), kv_grid(i),kx_grid(i),ky_grid(i), &
                kz_grid(i),freq_grid(i),gamma_grid(i)
        !   if(i==101) print*,theta_fs_grid(i),tprime_grid(i),&
        !        fprime_grid(i), kv_grid(i),kx_grid(i),ky_grid(i), &
        !        kz_grid(i),freq_grid(i),gamma_grid(i)
        !   print*, i
        end do
        close(unit=in_unit)
        
        nkz = 30
        dkz = 0.001
        kz_start = -0.001
        ndelkz = 10
        ddelkz = 0.003
        delkz_start = 0.

        allocate(kz_grid_test(nkz),delkz_grid(ndelkz))

        do ikz = 1, nkz
                kz_grid_test(ikz) = ikz*dkz+kz_start
        end do
        do idelkz = 1, ndelkz
                delkz_grid(idelkz) = idelkz*ddelkz+delkz_start
        end do

        do idelkz = 1, ndelkz
           delta_kz = delkz_grid(idelkz)
           denom = 0.
           do ikz = 1, nkz
              kz = kz_grid_test(ikz)   
              denom = denom + exp(-kz**2/delta_kz**2)*dkz
           end do
           denom = 2./delta_kz/sqrt(pi)*denom
           print*, denom
        end do


end program
