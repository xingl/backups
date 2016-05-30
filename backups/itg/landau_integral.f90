module landau_integral

        implicit none
        private

        public :: vz_integral

contains

  subroutine vz_integral(vy,integral,s,omega)

        implicit none

        integer, intent(in) :: s
        real(kind=8), intent(in) :: vy
        complex(kind=8), intent(in) :: omega
        real(kind=8) :: a, b
        complex (kind=8), intent(out) :: integral
        real(kind=8) :: x_left, x_center, x_right, dx
        complex (kind=8) :: f_left, f_center, f_right
        complex (kind=8) :: i_trapezoid, i_simpson

        a = vz_1
        b = vz_2
        dx = 0.1
        x_left = a
        x_right = x_left + dx
        f_left = integrand(x_left,vy,omega,s)
        integral = 0.0

        do while (x_left<b)
                x_center = 0.5*(x_left+x_right)
                f_center = integrand(x_center,vy,omega,s)
                f_right = integrand(x_right,vy,omega,s)
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
                        f_center = integrand(x_center,vy,omega,s)
                        f_right = integrand(x_right,vy,omega,s)
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
