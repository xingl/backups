program test

        implicit none

        real :: a,b
        integer :: i
        real (kind=8), parameter :: pi=3.141592653589793
        
        a = pi/6.
        b = pi/4.
        print*, sin(a)
        print*, sin(b)

        
end program
