program simpson_adaptive_integrator
        
        use file_utils

        implicit none

        complex (kind=8), parameter :: zi=(0.0,1.0)
        real (kind=8) :: a, b
        real (kind=8) :: x_left, x_center, x_right, dx
        complex (kind=8) :: f_left, f_center, f_right
        complex (kind=8) :: i_trapezoid, i_simpson, integral, i_exact
        complex (kind=8) :: integral_exact, ie_sum
        real (kind=8) :: tol, eps
        integer :: out_unit=101

        call init_file_utils
        call read_input_file

        call open_output_file(out_unit,'.dat')
        !write (out_unit, '(13a12)') 'x','i_simpson','','i_trapezoid','','i_exact',&
        !                                '','i_s-i_t','','i_s-i_e','','i_t-i_e',''
        write(out_unit,'(7a12)') 'x','int_s','','ie_sum','','diff',''
        x_left = a
        x_right = x_left + dx
        f_left = test_function(x_left, eps)
        integral = 0.0
        ie_sum = 0.0
        do while (x_left<b)
                x_center = 0.5*(x_left+x_right)
                f_center = test_function(x_center,eps) 
                f_right = test_function(x_right,eps)
                i_trapezoid = 0.5*(f_left+2.0*f_center+f_right)*0.5*dx
                i_simpson = 1.0/6.0*dx*(f_left+4.0*f_center+f_right)
                do while (abs(i_trapezoid-i_simpson)>tol)
                        dx=0.5*dx*(abs(i_trapezoid-i_simpson)/tol)**(-1.0/3.0)
                        if (dx>0.2) dx=0.2
                        x_right = x_left + dx
                        if (x_right>b) then
                                x_right = b
                                dx = x_right - x_left
                        end if
                        x_center = 0.5*(x_left+x_right)
                        f_center = test_function(x_center,eps) 
                        f_right = test_function(x_right,eps)
                        i_trapezoid = 0.5*(f_left+2.0*f_center+f_right)*0.5*dx
                        i_simpson = 1.0/6.0*dx*(f_left+4.0*f_center+f_right)
                end do
                print*, x_left
                i_exact = log(x_right+zi*eps) - log(x_left+zi*eps)
                ie_sum = ie_sum + i_exact
                !write(out_unit,'(13e12.4)') x_right,i_simpson,i_trapezoid,&
                !        i_exact, i_simpson-i_trapezoid, i_simpson-i_exact, &
                !        i_trapezoid-i_exact
                integral = integral + i_simpson
                write(out_unit,'(7e12.4)') x_right, integral, ie_sum, &
                                                integral-ie_sum
                x_left = x_right
                dx = 0.5*dx*(abs(i_trapezoid-i_simpson)/tol)**(-1.0/3.0)
                if (dx>0.2) dx=0.2
                x_right = x_left + dx
                if (x_right>b) then
                        x_right = b
                        dx = x_right - x_left
                end if
                f_left = f_right
        end do
        print*, integral
        integral_exact = log(b+zi*eps) - log(a+zi*eps)
        print*, integral_exact
        print*, ie_sum
        call close_output_file (out_unit)
        call finish_file_utils
                        
contains
   
function test_function(x,eps)
        real (kind=8) :: x, eps
        complex (kind=8) :: test_function
        test_function = 1.0/(x+zi*eps)
end function

subroutine read_input_file

        use file_utils, only: input_unit_exist, input_unit

        implicit none

        integer :: in_file
        logical :: exist

        namelist / parameters / a,b,eps,tol,dx

        a = -3.0
        b = 3.0
        eps = 1.0
        tol = 1.0D-05
        dx = 0.2

    in_file = input_unit_exist ("parameters", exist)
    if (exist) read (unit=input_unit("parameters"), nml=parameters)

  end subroutine read_input_file

end program
