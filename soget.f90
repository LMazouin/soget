program main

use global
use balda
use dimer

implicit none

integer :: ns, nel, nw
integer :: i, j , ik, iw, bc
real(dp) :: t, bigu, k
real(dp) :: mu0, mu
real(dp) :: beta, dbetadbigu, dbetadt
real(dp) :: eta, wmin, wmax, dw, w
real(dp) :: n0, n, d, dimp
real(dp) :: ts_balda, ec_balda, decdbigu_balda, decdt_balda, vhxc_balda
real(dp),allocatable :: vemb(:)
real(dp),allocatable :: vk(:), epsbath(:), cbath(:,:)
real(dp) :: gfimp_dimer, gfweiss_dimer
character*2 :: opt

real(dp) :: tsimp, ehximp, ecimp, deltavimp, decdbiguimp, decdtimp
real(dp) :: deltavks, deltavemb, deltavhxcimp, v0, v1
real(dp) :: eimp, ean0, ean1, ecat0, ecat1, sii(4), sij(4), sjj(4), ndimer, ddimer
real(dp) :: ehx, ec, etot

integer :: maxiter, iter
real(dp) :: alpha, conv
real(dp) :: res1, res2


integer, parameter :: limit = 1000
integer, parameter :: lenw = limit*4
real(dp) :: abserr
real(dp), parameter :: epsabs = 0.0_dp
real(dp), parameter :: epsrel = 1.0e-6
integer :: ier
integer :: iwork(limit)
integer :: last
integer :: neval
real(dp) :: work(lenw)

! read the input data

read(*,*) ns, nel
read(*,*) bigu, t
read(*,*) mu
read(*,*) wmin,wmax
read(*,*) nw
read(*,*) eta
read(*,*) maxiter
read(*,*) alpha
read(*,*) opt

! calculate the beta parameter of the BALDA

call beta_parameter(t, bigu, beta, dbetadt, dbetadbigu)

if (ns .eq. 2) t = 0.5_dp*t

dw = (wmax - wmin)/(float(nw) - 1.0_dp)

conv = 1.0e-6

! print input parameters

write(*,'(/a25,i14)')  'sites        : ', ns
write(*,'(a25,i14)')   'nelec        : ', nel
write(*,'(a25,f14.6)') 't (hopping)  : ', t
write(*,'(a25,f14.6)') 'U (Coulomb)  : ', bigu
write(*,'(a25,i14)')   'grid size    : ', nw
write(*,'(a25,f14.6)') 'grid spacing : ', dw
write(*,'(a25,f14.6)') 'w_min        : ', wmin
write(*,'(a25,f14.6)') 'w_max        : ', wmax
write(*,'(a25,f14.6)') 'eta          : ', eta
write(*,'(a25,i6)')    'maxiter      : ', maxiter
write(*,'(a25,f14.6)') 'damping      : ', alpha
write(*,'(a25,a14)')   'self-energy  : ', opt

! allocation

allocate(vemb(ns))
allocate(vk(ns-1))
allocate(epsbath(ns-1))
allocate(cbath(ns-1,ns-1))

! initialization

do i = 1,ns
vemb(i) = 0.0_dp
end do

! boundary conditions

if (mod(nel, 4) .eq. 0) then
bc = -1
else
bc = 1
end if

! initialization

do ik = 1, ns - 1
k = pi*ik/float(ns)
epsbath(ik) = - 2.0_dp*t*cos(k)
print*,epsbath(ik)
end do

do ik = 1, ns - 1
do i = 1, ns - 1
k = pi*ik*i/float(ns)
cbath(i,ik) = sqrt(2.0_dp/float(ns))*sin(k)
end do
end do

do ik = 1, ns - 1
vk(ik) = -t*cbath(1,ik) - bc*t*cbath(ns-1,ik)
print*,vk(ik)
end do

n0 = float(nel)/float(ns)

tsimp       = -2.0_dp*t*sqrt(n0*(2.0_dp - n0))
ehximp      = 0.0_dp 
ecimp       = 0.0_dp
decdtimp    = 0.0_dp
decdbiguimp = 0.0_dp
deltavimp   = 0.0_dp

ts_balda       = -4.0_dp*t*sin(0.5_dp*n0*pi)/pi
ec_balda       = 0.0_dp
decdbigu_balda = 0.25_dp*( sin(0.5_dp*pi*n0) - n0**2 ) - &
                 0.125_dp*pi*n0*cos(0.5_dp*pi*n0)
decdt_balda    = 0.0_dp
vhxc_balda     = 0.0_dp

iter = 0

! calculate the dimer parameters

call vhxc(bigu*0.5_dp, t, n0, tsimp, ehximp, ecimp, deltavimp, decdbiguimp, decdtimp)

deltavks     = -2.0_dp*t*(1.0_dp - n0)/sqrt(n0*(2.0_dp - n0))
deltavhxcimp = deltavks - deltavimp - 0.5_dp*bigu
deltavemb    = deltavimp + 0.5_dp*bigu

v0 = -deltavks
v1 = 0.0_dp

call dimer_param(bigu, t, deltavemb, eimp, ean0, ean1, ecat0, ecat1, sii, sij, sjj, ndimer, ddimer)

! print the dimer parameters

write(*,'(/a15, a15, a15, a15, a15, a15, a15, a15, a15, a15)') &
'iteration','Deltav^KS', 'Deltav^imp_Hxc', 'Deltav^emb', 'n^dimer', 'd^dimer', 'n_0', 'n', '|n0 - n|', 'alpha'

write(*,'(i15, f15.6, f15.6, f15.6, f15.6, f15.6, f15.6, f15.6, f15.6, f15.6)') &
iter, deltavks, deltavhxcimp, deltavemb, ndimer, ddimer, n0, n0, 0.0, alpha 

! calculate the BALDA potential

call v_balda(n0, t, bigu, beta, dbetadbigu, dbetadt, &
						 ec_balda, decdbigu_balda, decdt_balda, vhxc_balda)

! calculate the non-interacting chemical potential

mu0 = -2.0_dp*t*cos(0.5_dp*pi*n0)

if (opt .eq. '2L') then

! **************************************************************
! start the self-consistent loop on the impurity site occupation
! **************************************************************

do iter=1,maxiter

call vhxc(bigu*0.5_dp, t, n0, tsimp, ehximp, ecimp, deltavimp, decdbiguimp, decdtimp)

deltavks = -2.0_dp*t*(1.0_dp - n0)/sqrt(n0*(2.0_dp - n0))
deltavhxcimp  = deltavks - deltavimp - 0.5_dp*bigu
deltavemb = deltavimp + 0.5_dp*bigu

v0 = -deltavks
v1 = 0.0_dp

call dimer_param(bigu, t, deltavemb, eimp, ean0, ean1, ecat0, ecat1, sii, sij, sjj, ndimer, ddimer)

call v_balda(n0, t, bigu, beta, dbetadbigu, dbetadt, &
             ec_balda, decdbigu_balda, decdt_balda, vhxc_balda)

! calculate the bath orbital eigenvalues \eps_{k}

do ik=1,ns-1
k = pi*ik/float(ns)
epsbath(ik) = - 2.0_dp*t*cos(k)
end do

! calculated the bath orbitals (eigenvectors)
! it is actually the unitary transformation that
! diagonalizes the bath

do ik=1,ns-1
do i=1,ns-1
k = pi*ik*i/float(ns)
cbath(i,ik) = sqrt(2.0_dp/float(ns))*sin(k)
end do
end do

! calculate the impurity-bath coupling parameters V_{k}
! for the hybridization function \Delta(\omega)

do ik=1,ns-1
vk(ik) = -t*cbath(1,ik) - bc*t*cbath(ns-1,ik)
end do

n = 0.0_dp

!call dxintany(-1000.0_dp, 0.0_dp, gf_imp, n, 1.0e-12_dp)
call dqagi(gf_imp, 0.0_dp, -1, epsabs, epsrel, n, abserr, neval, &
					 ier, limit, lenw, last, iwork, work)

n = 2.0_dp*n/pi

if (abs(n0 - n) .lt. 1.0e-3_dp)  alpha = alpha/1.50_dp

if (abs(n0 - n) .lt. conv) exit

write(*,'(i15, f15.6, f15.6, f15.6, f15.6, f15.6, f15.6, f15.6, f15.6, f15.6)') &
iter, deltavks, deltavhxcimp, deltavemb, ndimer, ddimer, n0, n0, 0.0, alpha 

n0 = (1.0_dp - alpha)*n + alpha*n0

end do

if (iter .ge. maxiter .and. iter .gt. 2) then
write(*,'(a, f15.6, a, i15)') 'Convergence failed for U = ', bigu, ' nel = ', nel 
stop
end if

else if (opt .eq. 'HI') then

! ################################################
! works only at half-filling (just one interation)
! ################################################

! calculate the bath orbital eigenvalues \eps_{k}

do ik=1,ns-1
k = pi*ik/float(ns)
epsbath(ik) = - 2.0_dp*t*cos(k)
end do

! calculated the bath orbitals (eigenvectors)
! it is actually the unitary transformation that
! diagonalizes the bath

do ik=1,ns-1
do i=1,ns-1
k = pi*ik*i/float(ns)
cbath(i,ik) = sqrt(2.0_dp/float(ns))*sin(k)
end do
end do

! calculate the impurity-bath coupling parameters V_{k}
! for the hybridization function \Delta(\omega)

do ik=1,ns-1
vk(ik) = -t*cbath(1,ik) - bc*t*cbath(ns-1,ik)
end do

n = 0.0_dp
!call dxintany(-1000.0_dp,0.0_dp,gf_hi,n,1.0e-12_dp)
call dqagi(gf_hi, 0.0_dp, -1, epsabs, epsrel, n, abserr, neval, &
					 ier, limit, lenw, last, iwork, work)
n = 2.0_dp*n/pi

write(*,'(i15, f15.6, f15.6, f15.6, f15.6, f15.6, f15.6, f15.6, f15.6, f15.6)') &
1, deltavks, deltavhxcimp, deltavemb, ndimer, ddimer, n0, n, 0.0, alpha 

end if


if (opt .eq. '2L') then
!call dxintany(-1000.0_dp,0.0_dp,sfgf_imp,dimp,1.0e-12_dp)
call dqagi(sfgf_imp, 0.0_dp, -1, epsabs, epsrel, dimp, abserr, neval, &
					 ier, limit, lenw, last, iwork, work)
else if (opt .eq. 'HI') then
!call dxintany(-1000.0_dp,0.0_dp,sfgf_hi,dimp,1.0e-12_dp)
call dqagi(sfgf_hi, 0.0_dp, -1, epsabs, epsrel, dimp, abserr, neval, &
					 ier, limit, lenw, last, iwork, work)
end if
dimp = dimp/pi/bigu


write(*, '(/a15, a15, a15, a15, a15)') '# POTENTIALS', 'de^BA_Hxc/dn', 'Delta^imp_Hxc', 'mu_0', 'mu - U/2'
write(*,'(15x, f15.6, f15.6, f15.6, f15.6)') vhxc_balda, deltavhxcimp, mu0, mu - 0.5_dp*bigu

write(*, '(/a15, a15, a15)') '# OCCUPATION', 'N/L', 'n_0'

write(*,'(15x, f15.6, f15.6)')  float(nel)/float(ns), n0

write(*,'(/a20, a20, a20)') '# DOUBLE OCCUPATION', 'd^imp', 'd'
write(*,'(20x, f20.6, f20.6)')  dimp, dimp + decdbigu_balda - 0.5_dp*decdbiguimp 

ts_balda = -4.0_dp*t*sin(0.5_dp*n0*pi)/pi
ehx = 0.25*bigu*n0**2
etot     = ts_balda + bigu*dimp + t*decdt_balda + bigu*decdbigu_balda - 0.5_dp*bigu*decdbiguimp
write(*, '(/a20, a15, a15, a15, a20, a15)') '# PER-SITE ENERGY', 'ts', 'U*d^imp', 't*de^BA_c/dt', 'U*de^bath_c/dU', 'e'
write(*,'(20x, f15.6, f15.6, f15.6, f20.6, f15.6)') ts_balda, bigu*dimp, t*decdt_balda, &
bigu*decdbigu_balda - 0.5_dp*bigu*decdbiguimp, etot

ec = etot - ehx - ts_balda
write(*, '(/a20, a15, a15, a15, a15)') '# CORRELATION ENERGY', 'n0', 'U', 'ec', 'etot'
write(*,'(20x, f15.6, f15.6, f15.6, f15.6)') float(nel)/float(ns), bigu, ec, etot

write(*, '(/a35, a15, a15, a15, a25, a25)') '# GRAND-CANONICAL PER-SITE ENERGY', 'mu', 'N/L', 'n0', 'e(n0,N) - mu*n0', &
'eimp(n0,N) - mu*n0'
!write(*,'(35x, f15.6, f15.6, f25.6)')  mu, float(nel)/float(ns), etot - mu*float(nel)/float(ns)
write(*,'(35x, f15.6, f15.6, f15.6, f25.6, f25.6)')  mu, float(nel)/float(ns), n0, etot - mu*n0, &
ts_balda + bigu*dimp + t*decdt_balda - mu*n0

!open(unit=40,file='gfimp.cnt',status='replace')
!write(40,'(a)') '# INTERACTING IMPURITY GREEN FUNCTION'
!do iw=1,nw
!w = wmin + dw*float((iw-1))
!write(40,'(f14.6,f14.6)') w,gf_imp(w)/pi
!end do
!close(40)

!open(unit=45,file='sfimp.cnt',status='replace')
!write(45,'(a)') '# IMPURITY SELF-ENERGY'
!do iw=1,nw
!w = wmin + dw*float((iw-1))
!write(45,'(f14.6,f14.6,f14.6)') w,real(sf_imp(w)),-aimag(sf_imp(w)/pi)
!end do
!close(45)

contains

! function for the impurity Green's function with
! the Anderson dimer self-energy

function gf_imp(z) result(res)

implicit none

real(dp),intent(in) :: z
real(dp) :: res

complex(dp) :: w
complex(dp) :: gfweiss_dimer,gfimp_dimer
complex(dp) :: delta,gfimp,sfimp

w = z + im*eta

delta = cmplx(0.0_dp,0.0_dp)
do ik=1,ns-1
delta = delta + vk(ik)**2/(w + mu0 - epsbath(ik))
end do

gfimp_dimer = sii(1)/(w + eimp - ean0 ) + sii(2)/(w + eimp - ean1 ) + &
         sii(3)/(w - eimp + ecat0) + sii(4)/(w - eimp + ecat1)
gfweiss_dimer = 1.0_dp/(w - v0 - t**2/(w - v1))
sfimp = 1.0_dp/gfweiss_dimer - 1.0_dp/gfimp_dimer
gfimp = 1.0_dp/(w + mu0 - delta - sfimp)

res = -aimag(gfimp)

end function

! function for the impurity Green's function with
! the Hubbard-I self-energy

function gf_hi(z) result(res)

implicit none

real(dp),intent(in) :: z
real(dp) :: res

complex(dp) :: w
complex(dp) :: gfweiss_dimer,gfimp_dimer
complex(dp) :: delta,gfimp,sfimp

w = z + im*eta

delta = cmplx(0.0_dp,0.0_dp)
do ik=1,ns-1
delta = delta + vk(ik)**2/(w + mu0 - epsbath(ik))
end do

sfimp = 0.25_dp*bigu**2*(n0*(2.0_dp - n0)/(w + mu0 + 0.5_dp*bigu - 0.5_dp*bigu*(2.0_dp - n0))) 
gfimp = 1.0_dp/(w + mu0 - delta - sfimp)

res = -aimag(gfimp)

end function

! function for self-energy impurity Green's function product for
! the impurity double occupation (Anderson dimer)

function sfgf_imp(z) result(res)

implicit none

real(dp),intent(in) :: z
real(dp) :: res

complex(dp) :: w
complex(dp) :: gfweiss_dimer,gfimp_dimer
complex(dp) :: delta,gfimp,sfimp

w = z + im*eta

delta = cmplx(0.0_dp,0.0_dp)
do ik=1,ns-1
delta = delta + vk(ik)**2/(w + mu0 - epsbath(ik))
end do

gfimp_dimer = sii(1)/(w + eimp - ean0 ) + sii(2)/(w + eimp - ean1 ) + &
         sii(3)/(w - eimp + ecat0) + sii(4)/(w - eimp + ecat1)

gfweiss_dimer = 1.0_dp/(w - v0 - t**2/(w - v1))

sfimp = 1.0_dp/gfweiss_dimer - 1.0_dp/gfimp_dimer

gfimp = 1.0_dp/(w + mu0 - delta - sfimp)

res = -aimag(gfimp*sfimp) + aimag(gfimp)*deltavhxcimp

end function

! function for self-energy impurity Green's function product for
! the impurity double occupation (Hubbard-I)

function sfgf_hi(z) result(res)

implicit none

real(dp),intent(in) :: z
real(dp) :: res

complex(dp) :: w
complex(dp) :: gfweiss_dimer,gfimp_dimer
complex(dp) :: delta,gfimp,sfimp

w = z + im*eta

delta = cmplx(0.0_dp,0.0_dp)
do ik=1,ns-1
delta = delta + vk(ik)**2/(w + mu0 - epsbath(ik))
end do

sfimp = 0.25_dp*bigu**2*(n0*(2.0_dp - n0)/(w + mu0  + 0.5_dp*bigu- 0.5_dp*bigu*(2.0_dp - n0))) 
gfimp = 1.0_dp/(w + mu0 - delta - sfimp)
res = -aimag(gfimp*sfimp) - aimag(gfimp)*0.5_dp*bigu*n0

end function

function sf_imp(z) result(res)

implicit none

real(dp),intent(in) :: z
complex(dp) :: res

complex(dp) :: w
complex(dp) :: gfweiss_dimer,gfimp_dimer
complex(dp) :: delta,gfimp,sfimp

w = z + im*eta

delta = cmplx(0.0_dp,0.0_dp)
do ik=1,ns-1
delta = delta + vk(ik)**2/(w + mu - epsbath(ik))
end do

gfimp_dimer = sii(1)/(w + eimp - ean0 ) + sii(2)/(w + eimp - ean1 ) + &
         sii(3)/(w - eimp + ecat0) + sii(4)/(w - eimp + ecat1)
gfweiss_dimer = 1.0_dp/(w - v0 - t**2/(w - v1))
sfimp = 1.0_dp/gfweiss_dimer - 1.0_dp/gfimp_dimer - deltavhxcimp
res = sfimp

end function

end program
