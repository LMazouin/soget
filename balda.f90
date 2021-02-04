module balda

use global

implicit none

real(dp) :: bigu_tmp,t_tmp

private integrand, f
public v_balda, beta_parameter

contains

subroutine v_balda(n,t,bigu,beta,dbetadbigu,dbetadt, &
ec_balda,decdbigu_balda,decdt_balda,vhxc_balda)

real(dp),intent(in) :: beta,dbetadt,dbetadbigu
real(dp),intent(in) :: t, bigu
real(dp),intent(in) :: n
real(dp),intent(out) :: vhxc_balda,ec_balda,decdbigu_balda,decdt_balda


if (bigu .lt. 1.0e-6_dp) return

bigu_tmp = bigu
t_tmp = t

if (n .le. 1.0_dp) then
  ec_balda = - 2.0_dp*t_tmp*beta*sin(pi*n/beta)/pi + 4.0_dp*t_tmp*sin(pi*n/2.0_dp)/pi &
     - 0.25_dp*bigu_tmp*n**2
  decdbigu_balda = dbetadbigu*(-2.0_dp*t_tmp*sin(pi*n/beta)/pi + 2.0_dp*t_tmp*n*cos(pi*n/beta)/beta) - &
		0.25_dp*n**2
	decdt_balda = 2.0_dp*t_tmp*pi*n*dbetadt*cos(pi*n/beta)/(pi*beta) - & 
     (2.0_dp*t_tmp*dbetadt/pi + 2.0_dp*beta/pi)*sin(pi*n/beta) + 4.0_dp*sin(pi*n/2.0_dp)/pi
  vhxc_balda = -2.0_dp*t_tmp*cos(pi*n/beta) + 2.0_dp*t_tmp*cos(pi*n/2.0_dp)
else
  ec_balda = -2.0_dp*t_tmp*beta*sin((pi*(2.0_dp - n)/beta))/pi + &
     4.0_dp*t_tmp*sin((pi*(2.0_dp - n)/2.0_dp))/pi + bigu_tmp*(n - 1.0_dp) - 0.25_dp*bigu_tmp*n**2
  decdbigu_balda = dbetadbigu*(-2.0_dp*t_tmp*sin(pi*(2.0_dp - n)/beta)/pi + &
     2.0_dp*t_tmp*(2.0_dp - n)*cos(pi*(2.0_dp - n)/beta)/beta) + (n - 1.0_dp) - 0.25_dp*n**2
  decdt_balda = 2.0_dp*t_tmp*pi*(2.0_dp - n)*dbetadt*cos(pi*(2.0_dp - n)/beta)/(pi*beta) - &
     (2.0_dp*t_tmp*dbetadt/pi + 2.0_dp*beta/pi)*sin(pi*(2.0_dp - n)/beta) + &
     4.0_dp*sin(pi*(2.0_dp - n)/2.0_dp)/pi
  vhxc_balda =  2.0_dp*t_tmp*cos(pi*(2.0_dp - n)/beta) - 2.0_dp*t_tmp*cos(pi*(2.0_dp - n)/2.0_dp) + bigu_tmp
end if
!if (n .eq. 1.0_dp) vhxc_balda = 0.5_dp*bigu_tmp

return

end subroutine

subroutine beta_parameter(t, bigu, beta, dbetadt, dbetadbigu)

! This function calculates the beta parameters for the BALDA

real(dp) :: a
real(dp) :: b
real(dp) :: x
real(dp),parameter :: tol = 1.0e-9_dp
integer,parameter :: maxiter = 300
integer :: iter
integer :: ier

real(dp),intent(in) :: t, bigu
real(dp),intent(out) :: beta, dbetadt, dbetadbigu

a = 1.0_dp
b = 2.0_dp
bigu_tmp = bigu
t_tmp = t
call secant(f, a, b, x, tol, maxiter, iter, ier)
beta = x


a = 1.0_dp
b = 2.0_dp
bigu_tmp = bigu + 1.0e-6_dp
t_tmp = t
call secant(f, a, b, x, tol, maxiter, iter, ier)
dbetadbigu = (x - beta)/1.0e-6_dp
bigu_tmp = bigu

a = 1.0_dp
b = 2.0_dp
t_tmp = t + 1.0e-6_dp
call secant(f, a, b, x, tol, maxiter, iter, ier)
dbetadt = (x - beta)/1.0e-6_dp
t_tmp = t

return

end subroutine

function integrand(x) result(res)

real(dp),intent(in) :: x
real(dp) :: res
real(dp) :: num,den

num = bessel_j0(x)*bessel_j1(x)
den = x*(1.0_dp + exp(x*bigu_tmp/2.0_dp/t_tmp))

res = num/den

end function

function f(x) result(res)

real(dp),intent(in) :: x
real(dp) :: res
real(dp) :: a,b

integer,parameter :: limit=1000
integer,parameter :: lenw=limit*4
real(dp) :: abserr
real(dp),parameter :: epsabs=0.0_dp
real(dp),parameter :: epsrel=0.00001_dp
integer :: ier
integer :: iwork(limit)
integer,parameter :: inf=1
integer :: last
integer :: neval
real(dp) :: work(lenw)

call dqagi(integrand,0.0_dp,inf,epsabs,epsrel,a,abserr,neval, &
           ier,limit,lenw,last,iwork,work)

b = -sin(pi/x)*2.0_dp*t_tmp*x/pi

res = b + 4.0_dp*t_tmp*a

end function

end module
