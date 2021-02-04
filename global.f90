!******************!
! global variables !
!******************!

module global

! useful constants

integer,parameter :: dp = kind(0.d0)
complex(dp),parameter :: im = cmplx(0.0_dp,1.0_dp)
real(dp),parameter :: pi = 4.0_dp*atan(1.0_dp)

end module
