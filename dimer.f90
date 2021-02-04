module dimer

!**************************************************!
! exact embedding potential after Carrascal et al. ! 
!**************************************************!

use global

private
public vhxc,dimer_param,gf_dimer,sfimpz,sfimpw

contains

subroutine vhxc(bigu,t,n,ts,ehx,ec,deltavhxc,decdbigu,decdt)

implicit none

real(dp),intent(in) :: bigu
real(dp),intent(in) :: n
real(dp),intent(in) :: t

! Parametrization of the correlation energy !
real(dp) :: rho,u,twot
real(dp) :: dtsdn
real(dp) :: g0,g1
real(dp) :: a12fun,a21fun,a11fun,a22fun
real(dp) :: a1fun,a2fun
real(dp) :: g0fun,g1fun,hfun,ffun
real(dp) :: dfdrho,dfdg
real(dp) :: dhfundg0,dhfundg1,dg0fundlbda,dfun
! Implementation of the potential, differentiation wrt n !
real(dp) :: dhfundrho
real(dp) :: dnfundrho,ddfundrho
real(dp) :: da1fundrho,da2fundrho
real(dp) :: dg0fundrho,dg1fundrho
real(dp) :: partone,parttwo,part1,part2
real(dp) :: partoned,parttwod
real(dp) :: dhfundgdrho
real(dp) :: da21drho,da11drho
real(dp) :: dpart1drho,dpart2drho
real(dp) :: dg0dlbdadrho
real(dp) :: dg1drho
real(dp) :: sd
!for the differentiation wrt U !
real(dp) :: da1fundu,da2fundu
real(dp) :: g,gnum,gdenom
real(dp) :: dgdu,dgnumdu,dgdenomdu
real(dp) :: dg0du,dg1du
real(dp) :: u1,u2,v1,v2,w2
real(dp) :: u1prime,u2prime,v1prime,v2prime,w2prime
real(dp) :: dhdgdu,dg0dlbdadu,dhdu,dfdu
real(dp) :: dhdbigu,dFdbigu
!for the differentiation wrt t !
real(dp) :: dtsdt
real(dp) :: da1dt,da2dt
real(dp) :: afunc,bfunc,cfunc
real(dp) :: dafuncdt,dbfuncdt,dcfuncdt
real(dp) :: num,denom
real(dp) :: dnumdt,ddenomdt
real(dp) :: dg0dt,dg1dt
real(dp) :: dg0dlbdadt,dhdgdt,dhdt,dfdt

real(dp),intent(out) :: ts,ehx
real(dp),intent(out) :: ec
real(dp),intent(out) :: deltavhxc
real(dp),intent(out) :: decdt
real(dp),intent(out) :: decdbigu


! Note that my implementation is different from Laurent's one.
! Indeed, I did it again with wolfram and I think this one is slightly better.
! Don't worry, the differences are very small.

twot = 2.0_dp*t
u = bigu/twot

!!! START WITH THE DIFFERENTIATION WITH RESPECT TO n !!!

dtsdn = -twot*(1.0_dp - n)/sqrt(n*(2.0_dp - n))
rho = abs(1.0_dp - n)

if (rho .lt. 1.0e-6_dp) rho = 1.0e-6_dp

a12fun = 0.5_dp*(1.0_dp - rho)
a21fun = 0.5_dp*sqrt((1.0_dp - rho)*rho*0.5_dp)
a11fun = (1.0_dp + 1.0_dp/rho)*a21fun
a22fun = 0.5_dp*a12fun
a1fun = a11fun + u*a12fun
a2fun = a21fun + u*a22fun
g0fun = sqrt( ( (1.0_dp - rho)*( 1.0_dp + rho*( 1.0_dp + (1.0_dp + rho)**3*u*a1fun ) ) ) / &
 ( 1.0_dp + (1.0_dp + rho)**3*u*a2fun ) )
g0 = g0fun
dfun = 1.0_dp + (1.0_dp + rho)**3*u*a2fun
da1fundrho = -0.5_dp*u + ( ( rho*(1.0_dp - 2.0_dp*rho) - 1.0_dp ) / (16.0*rho*a21fun) )
da2fundrho = -0.25_dp*u + ( (1.0_dp - 2.0_dp*rho)/(16.0*a21fun) )
dnfundrho = -1.0_dp+(1.0_dp-2.0_dp*rho)*(1.0_dp+((1.0_dp+rho)**3)*u*a1fun)+&
   rho*u*(1.0_dp-rho)*((1.0_dp+rho)**2)*(3.0_dp*a1fun+(1.0_dp+rho)*da1fundrho)
ddfundrho = u*(1.0_dp + rho)**2*( 3.0_dp*a2fun + (1.0_dp + rho)*da2fundrho )
dhfundg0 = -g0*(2.0_dp*rho**2 + g0**2*(1.0_dp - sqrt(1.0_dp - g0**2 - rho**2))) &
   /((g0**2+rho**2)**2) + (2.0_dp*g0*(1.0_dp - sqrt(1.0_dp - g0**2 - rho**2)) + &
   g0**3/sqrt(1.0_dp - g0**2 - rho**2))/(2.0_dp*(g0**2+rho**2))
dg0fundlbda = ( ( (1.0_dp - rho)*(1.0_dp + rho)**3*u**2 )/( 2.0_dp*g0fun*( 1.0_dp + &
   (1.0_dp + rho)**3*u*a2fun )**2 ) )*( ( 3.0_dp*rho/2.0_dp - 1.0_dp + &
   rho*(1.0_dp + rho)**3*u*a2fun )&
   *a12fun - rho*( 1.0_dp + (1.0_dp + rho)**3*u*a1fun )*a22fun )
g1fun = g0 + (u*dhfundg0 - 1.0_dp)*dg0fundlbda
g1 = g1fun
hfun = (1.0_dp/(2.0_dp*(g1**2+rho**2)))*(2.0_dp*rho**2+g1**2*(1.0_dp-sqrt(1.0_dp-g1**2-rho**2)))
dhfundg1 = -g1*(2.0_dp*rho**2 + g1**2*(1.0_dp - sqrt(1.0_dp - g1**2 - rho**2))) &
   /((g1**2+rho**2)**2) + (2.0_dp*g1*(1.0_dp - sqrt(1.0_dp - g1**2 - rho**2)) + &
   g1**3/sqrt(1.0_dp - g1**2 - rho**2))/(2.0_dp*(g1**2+rho**2))
dg0fundrho = ( 1.0_dp/( 2.0_dp*g0fun*dfun ) )*( dnfundrho - g0fun**2*ddfundrho )
partone = (-1.0_dp*dg0fundrho*(g0**2 + rho**2)**2 + &
   4.0_dp*g0*(g0**2+rho**2)*(g0*dg0fundrho + rho))*(2*rho**2 + &
   g0**2*(1.0_dp - sqrt(1.0_dp - g0**2 - &
   rho**2)))/(g0**2 + rho**2)**4
parttwo = - g0*(4.0_dp*rho + 2.0_dp*g0*dg0fundrho*(1.0_dp &
   - sqrt(1.0_dp - g0**2 - rho**2)) + g0**2*(g0*dg0fundrho+ &
   rho)/sqrt(1.0_dp - g0**2 - rho**2))/(g0**2 + rho**2)**2
partoned = - (g0*dg0fundrho + rho)*(2*g0*(1.0_dp - sqrt(1.0_dp - g0**2 &
   - rho**2)) + g0**3/sqrt(1.0_dp - g0**2 - rho**2))/(g0**2 + rho**2)**2
parttwod = (2*dg0fundrho*(1.0_dp - sqrt(1.0_dp - g0**2 - rho**2)) + &
   2*g0*(g0*dg0fundrho + rho)/sqrt(1.0_dp - g0**2 - rho**2) + &
   (3.0_dp*g0**2*dg0fundrho*sqrt(1.0_dp - g0**2 - rho**2) + &
   g0**3*(g0*dg0fundrho + rho)/sqrt(1.0_dp - g0**2 - rho**2))/(1.0 - g0**2 - rho**2))&
   /(2*(g0**2 + rho**2))
dhfundgdrho = partoned + parttwo + partone + parttwod
da21drho = 0.25_dp*(0.5_dp - rho)/sqrt( (1.0_dp - rho)*0.5_dp*rho )
da11drho = da21drho*(1.0_dp + 1.0_dp/rho) - a21fun/(rho**2)
part1 = (1.0_dp - rho)*(1.0_dp + rho)**3*u**2*((3.0_dp*rho/2.0_dp - 1.0_dp + rho*(1.0_dp + & 
     rho)**3*u*a2fun)*a12fun - rho*(1.0_dp + (1.0_dp + rho)**3*u*a1fun)*a22fun)
part2 = 2.0_dp*g0*(1.0_dp + (1.0_dp + rho)**3*u*a2fun)**2
dpart1drho = u**2*(3.0_dp*(1.0_dp - rho)*(1.0_dp + rho)**2 - (1.0_dp + rho)**3)&
     *((3.0_dp*rho/2.0_dp - 1.0_dp + rho*(1.0_dp + rho)**3*u*a2fun)*a12fun - &
     rho*(1.0_dp + (1.0_dp + rho)**3*u*a1fun)*a22fun) + (1.0_dp - rho)*&
     (1.0_dp + rho)**3*u**2*( (3.0_dp/2.0_dp + 3.0_dp*u*(1.0_dp + rho)**2*rho*a2fun + &
     u*(1.0_dp + rho)**3*(a2fun + rho*da2fundrho))*a12fun + (3.0_dp*rho/2.0_dp - 1.0_dp + &
     rho*(1.0_dp + rho)**3*u*a2fun)*(-1.0_dp/2.0_dp) - (rho*(-1.0_dp/4.0_dp) + a22fun)*&
     (1.0_dp + (1.0_dp + rho)**3*u*a1fun) - rho*a22fun*(3.0_dp*(1.0_dp + rho)**2*u*a1fun&
     + (1.0_dp + rho)**3*u*da1fundrho) )
dpart2drho = 2.0_dp*dg0fundrho*(1.0_dp + (1.0_dp + rho)**3*u*a2fun)**2 + &
     4.0_dp*g0*(1.0_dp + (1.0_dp + rho)**3*u*a2fun)*u*(3.0_dp*(1.0_dp + rho)**2*a2fun +&
     (1.0_dp + rho)**3*da2fundrho)
dg0dlbdadrho = (dpart1drho*part2 - dpart2drho*part1)/(part2**2)
dg1fundrho = dg0fundrho + (u*dhfundg0 - 1.0_dp)*dg0dlbdadrho + u*dhfundgdrho*dg0fundlbda
dhfundrho = (4.0_dp*rho + 2.0_dp*g1*dg1fundrho*(1.0_dp - sqrt(1.0_dp - g1**2 - rho**2)) +&
    (g1**3*dg1fundrho + g1**2*rho)/sqrt(1.0_dp - g1**2 - rho**2))/(2.0_dp*(g1**2 + rho**2)) -&
    (g1*dg1fundrho + rho)*(2.0_dp*rho**2 + g1**2*(1.0_dp -&
    sqrt(1.0_dp - g1**2 - rho**2)))/((g1**2 + rho**2)**2)
dfdrho = -2.0_dp*t*dg1fundrho + bigu*dhfundrho
dfdg = -2.0_dp*t + bigu*dhfundg1

sd = sign(1.0_dp,n - 1.0_dp)

! Delta v_Hxc = - dEHxc/dn
!deltavhxc = -sd*dfdrho + dtsdn
deltavhxc = sd*dfdrho

! Ec
ts = -2.0_dp*t*sqrt(n*(2.0_dp-n))
ehx = bigu*(1.0+n*((0.5_dp*n)-1.0_dp))
ffun = -2.0_dp*g1 + bigu*hfun
ec = ffun - ts - ehx 

!!! FINISH !!!

!!! START THE DIFFERENTIATION WITH RESPECT TO U !!!
da2fundu = a22fun
da1fundu = a12fun
gnum = (1.0_dp - rho)*(1.0+rho*(1.0_dp+((1.0_dp+rho)**3)*u*a1fun))
gdenom = 1.0_dp+((1.0_dp+rho)**3)*u*a2fun
g = ((1.0_dp-rho)*(1.0_dp+rho*(1.0_dp+((1.0_dp+rho)**3)*u*a1fun)))/&
  (1.0_dp+((1.0_dp+rho)**3)*u*a2fun)
dgnumdu = (1.0_dp-rho)*rho*((1.0_dp+rho)**3)*(a1fun+u*da1fundu)
dgdenomdu = ((1.0_dp+rho)**3)*(a2fun+u*da2fundu)
dgdu = (dgnumdu*gdenom - gnum*dgdenomdu)/gdenom**2
dg0du = dgdu/(2.0_dp*sqrt(g))
u1 = g0*(g0**4 + 3*g0**2*rho**2 + 2.0_dp*rho**2*(rho**2-1.0_dp-sqrt(1.0_dp-g0**2-rho**2)))
u1prime = dg0du*(g0**4 + 3.0_dp*g0**2*rho**2 + 2.0_dp*rho**2*(rho**2-1.0_dp-&
    sqrt(1.0_dp-g0**2-rho**2))) + g0*(4.0_dp*g0**3*dg0du + 6.0_dp*g0*rho**2*dg0du + &
    2.0_dp*rho**2*g0*dg0du/sqrt(1.0_dp-g0**2-rho**2))
v1 = 2.0_dp*sqrt(1.0_dp-g0**2-rho**2)*(g0**2+rho**2)**2
v1prime = -2.0_dp*g0*dg0du*(g0**2+rho**2)**2/sqrt(1.0_dp-g0**2-rho**2) + &
    8.0_dp*sqrt(1.0_dp-g0**2-rho**2)*g0*dg0du*(g0**2+rho**2)
dhdgdu = (u1prime*v1-u1*v1prime)/v1**2
u2 = (1.0_dp-rho)*(1.0_dp+rho)**3*u**2
u2prime = 2.0_dp*(1.0_dp-rho)*(1.0_dp+rho)**3*u
v2 = (3.0_dp*rho/2.0_dp - 1.0_dp + rho*(1.0_dp+rho)**3*u*a2fun)*a12fun - rho*(1.0_dp+(1.0_dp+rho)**3*u*a1fun)*a22fun
v2prime = a12fun*(rho*(1.0_dp+rho)**3*(a2fun+u*da2fundu)) - &
    a22fun*(rho*(1.0_dp+rho)**3*(a1fun+u*da1fundu))
w2 = 2.0_dp*g0*(1.0_dp+(1.0_dp+rho)**3*u*a2fun)**2
w2prime = 4.0_dp*g0*(1.0_dp+(1.0_dp+rho)**3*u*a2fun)*(1.0_dp+rho)**3*(a2fun+u*da2fundu)&
    + 2.0_dp*dg0du*(1.0_dp+(1.0_dp+rho)**3*u*a2fun)**2
dg0dlbdadu = ((u2prime*v2+u2*v2prime)*w2-u2*v2*w2prime)/w2**2
dg1du = dg0du + (dhfundg0 + u*dhdgdu)*dg0fundlbda + (u*dhfundg0 - 1.0_dp)*dg0dlbdadu
dhdu = dg1du*dhfundg1
dhdbigu = (1.0_dp/(2.0_dp*t))*dhdu
dfdu = -2.0_dp*t*(dg1du - hfun - u*dhdu)
dfdbigu = (1.0_dp/(2.0_dp*t))*dFdu
decdbigu = dfdbigu - (1.0_dp+n*(0.5_dp*n - 1.0_dp))
!!! FINISH !!!

!!! START THE DIFFERENTIATION WITH RESPECT TO t !!!
dtsdt = -2.0_dp*sqrt(n*(2.0_dp-n))
da1dt = -2.0_dp*(u/(2.0_dp*t))*a12fun
da2dt = -2.0_dp*(u/(2.0_dp*t))*a22fun
afunc = (1.0_dp - rho)*(1.0_dp + rho*(1.0_dp + (1.0_dp + rho)**3*u*a1fun))
bfunc = 1.0_dp + (1.0_dp + rho)**3*u*a2fun
cfunc = afunc/bfunc
dafuncdt = (1.0_dp - rho)*rho*(1.0_dp + rho)**3*(-2.0_dp*(u/(2.0_dp*t))*a1fun + u*da1dt)
dbfuncdt = (1.0_dp + rho)**3*(-2.0_dp*(u/(2.0_dp*t))*a2fun + u*da2dt)
dcfuncdt = (dafuncdt*bfunc - afunc*dbfuncdt)/bfunc**2
dg0dt = dcfuncdt/(2.0_dp*sqrt(cfunc))
num = (1.0_dp - rho)*(1.0_dp + rho)**3*u**2*((3.0_dp*rho/2.0_dp - 1.0_dp + &
    rho*(1.0_dp + rho)**3*u*a2fun)*a12fun - rho*(1.0_dp + (1.0_dp + rho)**3*u*a1fun)*a22fun)
denom = 2.0_dp*g0*(1.0_dp + (1.0_dp + rho)**3*u*a2fun)**2
dnumdt = (1.0_dp - rho)*(1.0_dp + rho)**3*(-((bigu**2)/(2.0_dp*t**3))*&
    ((3.0_dp*rho/2.0_dp - 1.0_dp + rho*(1.0_dp + rho)**3*u*a2fun)*a12fun - &
    rho*(1.0_dp + (1.0_dp + rho)**3*u*a1fun)*a22fun) + (u**2)*&
    (rho*(1.0_dp + rho)**3*(-(2.0_dp*u/(2.0_dp*t))*a2fun + u*da2dt)*a12fun - &
    rho*(1.0_dp + rho)**3*(-(2.0_dp*u/(2.0_dp*t))*a1fun + u*da1dt)*a22fun))
ddenomdt = 2.0_dp*dg0dt*(1.0_dp + (1.0_dp + rho)**3*u*a2fun)**2 + &
    4.0_dp*g0*(1.0_dp + (1.0_dp + rho)**3*u*a2fun)*&
    ((1.0_dp + rho)**3*(-(2.0_dp*u/(2.0_dp*t))*a2fun + u*da2dt))
dg0dlbdadt = (dnumdt*denom - num*ddenomdt)/denom**2
dhdgdt = dg0dt*( - 2.0_dp*g0*(2.0_dp*g0*(1.0_dp - &
     sqrt(1.0_dp - g0**2 - rho**2)) + g0**3/sqrt(1.0_dp - &
     g0**2 - rho**2))/((g0**2 + rho**2)**2) - ((g0**2 + rho**2)**2 - &
     4.0_dp*g0**2*(g0**2 + rho**2))*(2.0_dp*rho**2 + g0**2*(1.0_dp - sqrt(1.0_dp - g0**2 - &
     rho**2)))/((g0**2 + rho**2)**4) + (2.0_dp*(1.0_dp - sqrt(1.0_dp - g0**2 - rho**2)) + &
     2.0_dp*g0**2/sqrt(1.0_dp - g0**2 - rho**2) + (3.0_dp*g0**2*sqrt(1.0_dp - g0**2 - rho**2) + &
     g0**4/sqrt(1.0_dp - g0**2 - rho**2))/(1.0_dp - g0**2 - rho**2))/(2.0_dp*(g0**2 + rho**2)))
dg1dt = dg0dt - ((2.0_dp*u/(2.0_dp*t))*dhfundg0 - &
     u*dhdgdt)*dg0fundlbda + (u*dhfundg0 - 1.0_dp)*dg0dlbdadt
dhdt = -(dg1dt*g1/((g1**2 + rho**2)**2))*(2.0_dp*rho**2 + g1**2*(1.0_dp - sqrt(1.0_dp -&
     g1**2 - rho**2))) + (1.0_dp/(2.0_dp*(g1**2 + rho**2)))*(2.0_dp*g1*dg1dt*(1.0_dp -&
     sqrt(1.0_dp - g1**2 - rho**2)) + g1**3*dg1dt/sqrt(1.0_dp - g1**2 - rho**2))
dfdt =- 2.0_dp*g1 - 2.0_dp*t*dg1dt + bigu*dhdt
decdt = dfdt - dtsdt

return

end subroutine

subroutine dimer_param(bigu,t,deltav,e,ean0,ean1,ecat0,ecat1,sii,sij,sjj,n,d)

implicit none

real(dp),intent(in) :: bigu
real(dp),intent(in) :: t
real(dp),intent(in) :: deltav

real(dp) :: twot
real(dp) :: bigu0,bigu1,v0,v1
real(dp) :: a0,a1,a2,q,r,theta
real(dp) :: a,b
real(dp) :: c(3)
real(dp) :: can0(2),can1(2)
real(dp) :: ccat0(2),ccat1(2)

real(dp) :: e,ean0,ean1,ecat0,ecat1
real(dp),intent(out) :: sii(4),sij(4),sjj(4)
real(dp),intent(out) :: n
real(dp),intent(out) :: d

twot = 2.0_dp*t

v0 = -deltav
v1 =  0.0_dp

bigu0 = bigu
bigu1 = 0.0_dp

a0 = (v0 + v1)*(4.0_dp*t**2 - bigu0*bigu1 - 4.0_dp*v0*v1) + &
(bigu0 + bigu1)*2.0_dp*(t**2 - v0*v1) - 2.0_dp*(bigu0*v1**2 + bigu1*v0**2)

a1 = bigu0*bigu1 + 8.0_dp*v0*v1 - 4.0_dp*t**2 + 2.0_dp*(v0**2 + v1**2) + &
bigu0*(v0 + 3.0_dp*v1) + bigu1*(v1 + 3.0_dp*v0)

a2 = -(bigu0 + bigu1) - 3.0_dp*(v0 + v1)

q = (3.0_dp*a1 - a2**2)/9.0_dp
r = (9.0_dp*a2*a1 - 27.0_dp*a0 - 2.0_dp*a2**3)/54.0_dp 

!q = -bigu0**2/9.0_dp + (3.0_dp*deltav + bigu1)*bigu0/9.0_dp - 4.0_dp*t**2/3.0_dp - & 
! deltav**2/3.0_dp - bigu1*deltav/3.0_dp - bigu1**2/9.0_dp
!r = -(bigu0 + bigu1)*(18.0_dp*t**2 - 2.0_dp*bigu0**2 + 5.0_dp*bigu0*bigu1 + &
! 9.0_dp*bigu0*deltav - 2.0_dp*bigu1**2 - 9.0_dp*bigu1*deltav - 9.0_dp*deltav**2)/54.0_dp

theta = acos(r/sqrt(-q)**3)

e = 2.0_dp*sqrt(-q)*cos(theta/3.0_dp + 2.0_dp*pi/3.0_dp) + (bigu0 + bigu1)/3.0_dp + v0 + v1

a = bigu0 + 2.0_dp*v0 - e
b = bigu1 + 2.0_dp*v1 - e

c(1) = twot*b/sqrt( (a**2 + b**2)*twot**2 + 2.0_dp*a**2*b**2 ) 
c(2) = twot*a/sqrt( (a**2 + b**2)*twot**2 + 2.0_dp*a**2*b**2 )
c(3) = a*b/sqrt( (a**2 + b**2)*twot**2 + 2.0_dp*a**2*b**2 )

! c(4) = c(3) 

!n = 2.0_dp*c(1)**2 + 2.0_dp*c(3)**2
d = c(1)**2

ean0 = 0.5_dp*(bigu0 + bigu1) + 1.5_dp*(v0 + v1) - &
0.5_dp*sqrt(twot**2 + ((v1 - v0) + (bigu1 - bigu0))**2)
ean1 = 0.5_dp*(bigu0 + bigu1) + 1.5_dp*(v0 + v1) + &
0.5_dp*sqrt(twot**2 + ((v1 - v0) + (bigu1 - bigu0))**2)

can0(1) = t/sqrt(t**2 + (bigu0 + 2.0_dp*v0 + v1 - ean0)**2 )
can0(2) = (bigu0 + 2.0_dp*v0 + v1 - ean0)/sqrt(t**2 + (bigu0 + 2.0_dp*v0 + v1 - ean0)**2)
can1(1) = t/sqrt(t**2 + (bigu0 + 2.0_dp*v0 + v1 - ean1)**2 )
can1(2) = (bigu0 + 2.0_dp*v0 + v1 - ean1)/sqrt(t**2 + (bigu0 + 2.0_dp*v0 + v1 - ean1)**2)

ecat0 = 0.5_dp*(v0 + v1) - 0.5_dp*sqrt(twot**2 + (v1 - v0)**2)
ecat1 = 0.5_dp*(v0 + v1) + 0.5_dp*sqrt(twot**2 + (v1 - v0)**2)

ccat0(1) = t/sqrt(t**2 + (v0 - ecat0)**2)
ccat0(2) = (v0 - ecat0)/sqrt(t**2 + (v0 - ecat0)**2)
ccat1(1) = t/sqrt(t**2 + (v0 - ecat1)**2)
ccat1(2) = (v0 - ecat1)/sqrt(t**2 + (v0 - ecat1)**2)

sii(1) = ( c(3)*can0(1) + c(2)*can0(2) )**2
sii(2) = ( c(3)*can1(1) + c(2)*can1(2) )**2
sii(3) = ( -c(1)*ccat0(1) - c(3)*ccat0(2) )**2
sii(4) = ( -c(1)*ccat1(1) - c(3)*ccat1(2) )**2

sjj(1) = ( -c(1)*can0(1) - c(3)*can0(2) )**2
sjj(2) = ( -c(1)*can1(1) - c(3)*can1(2) )**2
sjj(3) = ( -c(3)*ccat0(1) - c(2)*ccat0(2) )**2
sjj(4) = ( -c(3)*ccat1(1) - c(2)*ccat1(2) )**2

sij(1) = ( c(3)*can0(1) + c(2)*can0(2) )*( -c(1)*can0(1) - c(3)*can0(2) )
sij(2) = ( c(3)*can1(1) + c(2)*can1(2) )*( -c(1)*can1(1) - c(3)*can1(2) )
sij(3) = ( -c(1)*ccat0(1) - c(3)*ccat0(2) )*( -c(3)*ccat0(1) - c(2)*ccat0(2) )
sij(4) = ( -c(1)*ccat1(1) - c(3)*ccat1(2) )*( -c(3)*ccat1(1) - c(2)*ccat1(2) )

n = 2.0_dp*(sii(3) + sii(4))

return

end subroutine

subroutine gf_dimer(n,t,bigu,deltav,deltavhxcimp,mu,e,ean0,ean1,ecat0,ecat1,sii,sij,sjj,eta,wmin,dw,nw, & 
shift,gfdimer,sfdimer)

implicit none

integer,intent(in) :: nw
real(dp),intent(in) :: n,t,bigu,deltav,deltavhxcimp,mu
real(dp),intent(in) :: eta,wmin,dw
real(dp),intent(in) :: e,ean0,ean1,ecat0,ecat1
real(dp),intent(in) :: sii(4),sij(4),sjj(4)

integer :: i,j,iw
real(dp) :: w,v0,v1
complex(dp) :: detgf,detgf0,gfweiss
complex(dp),dimension(0:1,0:1) :: gf,gf0,sf

real(dp),intent(out) :: shift
complex(dp),dimension(nw),intent(out) :: sfdimer
complex(dp),dimension(nw),intent(out) :: gfdimer
complex(dp),dimension(nw) :: gfdimer0

v0 = -deltav
v1 =  0.0_dp


w = 0.0_dp

	gf(0,0) = sii(1)/(w + e - ean0 + im*eta) + sii(2)/(w + e - ean1 + im*eta) + &
sii(3)/(w - e + ecat0 + im*eta) + sii(4)/(w - e + ecat1 + im*eta)
	
	gf(0,1) = sij(1)/(w + e - ean0 + im*eta) + sij(2)/(w + e - ean1 + im*eta) + &
sij(3)/(w - e + ecat0 + im*eta) + sij(4)/(w - e + ecat1 + im*eta)
	
	gf(1,0) = sij(1)/(w + e - ean0 + im*eta) + sij(2)/(w + e - ean1 + im*eta) + &
sij(3)/(w - e + ecat0 + im*eta) + sij(4)/(w - e + ecat1 + im*eta)
	
	gf(1,1) = sjj(1)/(w + e - ean0 + im*eta) + sjj(2)/(w + e - ean1 + im*eta) + &
sjj(3)/(w - e + ecat0 + im*eta) + sjj(4)/(w - e + ecat1 + im*eta)

	gf0(0,0) = 0.5_dp*( (1.0_dp + deltav/sqrt(4.0_dp*t**2 + deltav**2)) / &
(w + 0.5_dp*sqrt(4.0_dp*t**2 + deltav**2) + im*eta) + &
(1.0_dp - deltav/sqrt(4.0_dp*t**2 + deltav**2)) / &
(w - 0.5_dp*sqrt(4.0_dp*t**2 + deltav**2) + im*eta) ) 
	
	gf0(0,1) = 1.0_dp*( (t/sqrt(4.0_dp*t**2 + deltav**2)) / &
(w + 0.5_dp*sqrt(4.0_dp*t**2 + deltav**2) + im*eta) - &
(t/sqrt(4.0_dp*t**2 + deltav**2)) / &
(w - 0.5_dp*sqrt(4.0_dp*t**2 + deltav**2) + im*eta) ) 
	
	gf0(1,0) = 1.0_dp*( (t/sqrt(4.0_dp*t**2 + deltav**2)) / &
(w + 0.5_dp*sqrt(4.0_dp*t**2 + deltav**2) + im*eta) - &
(t/sqrt(4.0_dp*t**2 + deltav**2)) / &
(w - 0.5_dp*sqrt(4.0_dp*t**2 + deltav**2) + im*eta) ) 

	gf0(1,1) = 0.5_dp*( (1.0_dp - deltav/sqrt(4.0_dp*t**2 + deltav**2)) / &
(w + 0.5_dp*sqrt(4.0_dp*t**2 + deltav**2) + im*eta) + &
(1.0_dp + deltav/sqrt(4.0_dp*t**2 + deltav**2)) / &
(w - 0.5_dp*sqrt(4.0_dp*t**2 + deltav**2) + im*eta) ) 

	detgf0 = gf0(0,0)*gf0(1,1) - gf0(0,1)*gf0(1,0)
	detgf  = gf(0,0)*gf(1,1) - gf(0,1)*gf(1,0)

	gfweiss = 1.0_dp/(w + im*eta - v0 - t**2/(w + im*eta - v1))
	shift = real(1.0_dp/gfweiss - 1.0_dp/gf(0,0))

do iw=1,nw
	w = wmin + dw*float(iw-1) + mu


	gf(0,0) = sii(1)/(w + e - ean0 + im*eta) + sii(2)/(w + e - ean1 + im*eta) + &
sii(3)/(w - e + ecat0 + im*eta) + sii(4)/(w - e + ecat1 + im*eta)
	
	gf(0,1) = sij(1)/(w + e - ean0 + im*eta) + sij(2)/(w + e - ean1 + im*eta) + &
sij(3)/(w - e + ecat0 + im*eta) + sij(4)/(w - e + ecat1 + im*eta)
	
	gf(1,0) = sij(1)/(w + e - ean0 + im*eta) + sij(2)/(w + e - ean1 + im*eta) + &
sij(3)/(w - e + ecat0 + im*eta) + sij(4)/(w - e + ecat1 + im*eta)
	
	gf(1,1) = sjj(1)/(w + e - ean0 + im*eta) + sjj(2)/(w + e - ean1 + im*eta) + &
sjj(3)/(w - e + ecat0 + im*eta) + sjj(4)/(w - e + ecat1 + im*eta)

	gf0(0,0) = 0.5_dp*( (1.0_dp + (v1 - v0)/sqrt(4.0_dp*t**2 + (v1 - v0)**2)) / &
(w - 0.5_dp*(v0 + v1) + 0.5_dp*sqrt(4.0_dp*t**2 + (v1 - v0)**2) + im*eta) - &
(1.0_dp - (v1 - v0)/sqrt(4.0_dp*t**2 + (v1 - v0)**2)) / &
(w - 0.5_dp*(v0 + v1) - 0.5_dp*sqrt(4.0_dp*t**2 + (v1 - v0)**2) + im*eta) ) 
	
	gf0(0,1) = 1.0_dp*( (t/sqrt(4.0_dp*t**2 + (v1 - v0)**2)) / &
(w - 0.5_dp*(v0 + v1) + 0.5_dp*sqrt(4.0_dp*t**2 + (v1 - v0)**2) + im*eta) - &
(t/sqrt(4.0_dp*t**2 + (v1 - v0)**2)) / &
(w - 0.5_dp*(v0 + v1) - 0.5_dp*sqrt(4.0_dp*t**2 + (v1 - v0)**2) + im*eta) ) 
	
	gf0(1,0) = 1.0_dp*( (t/sqrt(4.0_dp*t**2 + (v1 - v0)**2)) / &
(w - 0.5_dp*(v0 + v1) + 0.5_dp*sqrt(4.0_dp*t**2 + (v1 - v0)**2) + im*eta) - &
(t/sqrt(4.0_dp*t**2 + (v1 - v0)**2)) / &
(w - 0.5_dp*(v0 + v1) - 0.5_dp*sqrt(4.0_dp*t**2 + (v1 - v0)**2) + im*eta) ) 

	gf0(1,1) = 0.5_dp*( (1.0_dp - (v1 - v0)/sqrt(4.0_dp*t**2 + (v1 - v0)**2)) / &
(w - 0.5_dp*(v0 + v1) + 0.5_dp*sqrt(4.0_dp*t**2 + (v1 - v0)**2) + im*eta) - &
(1.0_dp + (v1 - v0)/sqrt(4.0_dp*t**2 + (v1 - v0)**2)) / &
(w - 0.5_dp*(v0 + v1) - 0.5_dp*sqrt(4.0_dp*t**2 + (v1 - v0)**2) + im*eta) ) 

	detgf0 = gf0(0,0)*gf0(1,1) - gf0(0,1)*gf0(1,0)
	detgf  = gf(0,0)*gf(1,1) - gf(0,1)*gf(1,0)

	!sf(0,0) =  gf0(1,1)/detgf0 - gf(1,1)/detgf
	!sf(1,1) =  gf0(0,0)/detgf0 - gf(0,0)/detgf
	!sf(0,1) =  -gf0(1,0)/detgf0 + gf(1,0)/detgf
	
	gfweiss = 1.0_dp/(w + im*eta - v0 - t**2/(w + im*eta - v1))
	sf(0,0) = 1.0_dp/gfweiss - 1.0_dp/gf(0,0) - deltavhxcimp

	! some tests
	!write(*,'(f28.20,f28.20,f28.20,f28.20,f28.20)') w,-aimag(sf(1,1)),-aimag(sf(0,1)), &
	!real(sf(1,1)),real(sf(0,1))

	sfdimer(iw) = sf(0,0)
	gfdimer(iw) = gf(0,0) 
	gfdimer0(iw) = gf0(0,0)

end do

return

end subroutine

end module
