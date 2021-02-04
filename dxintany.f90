			subroutine dxintany (xl, xu, g, res, eps)
			use global
!
! 		this was the old sintax:
! 		subroutine dintany(xl,xu,g,res,eps,n,a,naux)
!
! 		integration of function with singularities or peaks 
! 		(intany) Neil H. Callwood (Oct. 1980)
!
! 		Intany is a general purpose integration routine,
! 		particularly useful for integrating functions with peaks
! 		or singularities. It uses binary splitting to isolate
! 		the peaks and further subdivide them, giving a result
! 		guaranteed to an absolute precision epsilon.
! 		A suggested method for integrating around a singularity
! 		x=c, is to divide in two intervals in the main program,
! 		one up to c-delta and the other from c+delta onwards.
! 		Then decrease delta until a satisfactory result is obtained.
!
! 		input parameters...  (double precision)
! 		                    xl   lower bound
! 		                    xu   upper bound
! 		                    g    function name of integrand (external)
! 		                    eps  absolute precision required
! 		                    a    auxiliary array, 50 elements is usuall
! 		                    naux dimension of a
! 		output parameters...  (double precision)
! 		                    res  result
! 		                    n    no. of splits performed
! 		suggestion call the subroutine with eps smaller (in an order of 
! 		magnitude) than the accuracy you need, because the
! 		result comes from the summation of n+1 integrals, 
! 		each one with an error < eps. (n = 70-150)      
! 		routines called...
! 		dqg8(ibm), Gaussian quadrature with 8 points
!
! 		error conditions...
!    			               message = 'auxiliary table full.'
!    			               this message indicates that less
!      		               precision must be specified, or that
!      		               the auxiliary table must be made
!      		               larger.  the contents of the table
!    			               are printed, each subinterval is
!     		               represented by 3 entries, as follows,
!     		               lower limit,
!     		               upper limit,
!      		               and best known value of integral to date
!
			implicit none
!     ---- parameters -----
			real(dp),intent(in) :: xl,xu,eps
			real(dp),intent(out) :: res
!     ---- local data -----
			integer :: n, k, i, kk
			integer, parameter :: naux = 500
			real(dp) :: a(naux),w,r,rr,av,dif
			real(dp), external :: g
!
!     ----- code section -----
!
! 		write(*,*) '  dxintany'
			w = 0.5_dp*(xu - xl)
			n = 0
! 		Note that k is a pointer to the last row of the aux. table, 
! 		which i ''lifo'' stack.
			k=1
			a(1) = xl
			a(2) = xu
			res = 0.0_dp
! 		Put first estimate of result into a(3).
			call dqg8(xl,xu,g,a(3))
! 		Start main loop.
100   continue
      if ( k+5 .gt. naux) go to 900
200   n = n + 1
!     Divide the last interval of stack in two and put rr = sum 
!     of result each half.
!	    r = best known result for this interval.
	    r = a(k+2)
      a(k+4) = a(k+1)
      a(k+1) = 0.5_dp*(a(k) + a(k+1))
      a(k+3) = a(k+1)
      call dqg8(a(k),a(k+1),g,a(k+2))
      k = k + 3
      call dqg8(a(k),a(k+1),g,a(k+2))
      rr = a(k-1) + a(k+2)
!     if result is bad go to subdivide more.
      if (abs(r-rr) .gt. eps) go to 100
!     the result is good, add to total and cross off stack.
      res = res + rr
!     print*,'dx...', n, res
      k = k - 6
      if (k .gt. 0) go to 200
!     Stack empty, return.
!     write(*,*) 'n =', n
      return
!     Error message.
900   kk = k + 2
      write(6,*) 'problems at   ',a(k),a(k+1)
      dif = a(k+1) - a(k)
      av = a(k+2)/dif
      write(6,*) 'av. function  ',av
      write(6,910) (a(i),i=1,kk)
  910 format(//' intany... aux.table full'//50(/6e20.10))
      stop
      end
!
!        subroutine dqg8
!        purpose
!           to compute integral(fct(x), summed over x from xl to xu)
!        usage
!           call dqg8 (xl,xu,fct,y)
!           parameter fct requires an external statement
!
!        description of parameters
!           xl     - double precision lower bound of the interval.
!           xu     - double precision upper bound of the interval.
!           fct    - the name of an external double precision function
!                    subprogram used.
!           y      - the resulting double precision integral value.
!
!           The external double precision function subprogram fct(x)
!           must be furnished by the user.
!
!        method
!           Evaluation is done by means of 8-point gauss quadrature
!           formula, which integrates polynomials up to degree 15
!           exactly. For reference, see
!           V. I. Krylov, Approximate calculation of integrals,
!           MacMillan, New york/London, 1962, pp. 100-111 and 337-340.
!     ................................................................
!
      subroutine dqg8(xl,xu,fct,y)
      use global
			implicit none
      real(dp),intent(in) ::  xl,xu
			real(dp),intent(out) :: y
			real(dp) :: a,b,c
      real(dp),external :: fct
!     write(*,*) '  dqg8'
      a = 0.5_dp*(xu + xl)
      b = xu - xl
      c = 0.48014492824877_dp*b
      y = 0.50614268145188e-1_dp*(fct(a+c) + fct(a-c))
      c = 0.39833323870681d0*b
      y = y + 0.11119051722668_dp*(fct(a+c) + fct(a-c))
      c = 0.26276620495816_dp*b
      y = y + 0.15685332293894_dp*(fct(a+c) + fct(a-c))
      c = 0.9171732124782e-1_dp*b
      y = b*(y + 0.18134189168918_dp*(fct(a+c) + fct(a-c)))
      return
      end
