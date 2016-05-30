program kz_avg
 
        implicit none

        real (kind=8), parameter :: pi=3.141592653589793
        complex (kind=8), parameter :: zi=(0.0,1.0)
        integer (kind=8) :: nkz,ndelkz
        real (kind=8) :: dkz,ddelkz
        integer (kind=8) :: ikz,idelkz
        real (kind=8) :: kz_start,delkz_start
        real (kind=8), dimension(:), allocatable :: kz_grid,delkz_grid
        real (kind=8) :: kz,delta_kz
        real (kind=8) :: denom
        
        nkz = 30
        dkz = 0.001
        kz_start = -0.001
        ndelkz = 10
        ddelkz = 0.003
        delkz_start = 0.

        allocate(kz_grid(nkz),delkz_grid(ndelkz))

        do ikz = 1, nkz
                kz_grid(ikz) = ikz*dkz+kz_start
        end do
        do idelkz = 1, ndelkz
                delkz_grid(idelkz) = idelkz*ddelkz+delkz_start
        end do

        do idelkz = 1, ndelkz
           delta_kz = delkz_grid(idelkz)
           denom = 0.
           do ikz = 1, nkz
              kz = kz_grid(ikz)   
              denom = denom + exp(-kz**2/delta_kz**2)*dkz
           end do
           denom = 2./delta_kz/sqrt(pi)*denom
           print*, denom
        end do

end program
