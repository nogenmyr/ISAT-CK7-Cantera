!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the ISAT-CK7     !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

	subroutine gecp( nr, nc, a, lda, b, x, info, ipiv )

!  solve the (possibly underdetermined) system  a x = b using Gaussian
!  elimination with column pivoting (and row pivoting if required).

!  input:
!	nr	- number of rows of a
!	nc	- number of columns of a
!	a	- matrix
!	lda	- leading dimension of a
!	b	- right-hand-side
!  output:
!	x	- solution
!	info	=  0 successful
!		= -1 failed due to zero pivot
!  workspace
!	ipiv

	use isat_abort_m
	implicit real(kind(1.d0)) (a-h, o-z), integer (i-n)

	integer, parameter :: k_dp = kind(1.d0)
	real(k_dp) :: a(lda,nc), b(nr), x(nc)
	integer    :: ipiv(nc)

	data epslon/ 0.d0 /
	data nwgecp, iwarn/ 5, 0 /

!  check input

	if( nr < 1  .or.  nc < nr ) then
	   write(lu_err,*)'gecp: must have 0 < nr <= nc', nr, nc
	   call isat_abort('gecp',1)
	endif

!  set pointer to elements of x

	do i  = 1, nc
	   ipiv(i)  = i
	end do

	info     = 0

!-------  loop over rows to perform elimination  -----------------------------

	do 100 j = 1, nr

!  look for column pivot along row j

110	continue
	pivot    = -1.

	do i = j, nc
	   absa     = abs( a(j,i) )
	   if( absa > pivot ) then
	      pivot = absa
	      ip    = i
	   endif
	end do

!  pivot found?

	jlast = j
	if( pivot .gt. epslon ) go to 138

!  look for any non-zero element

	do jj = j+1, nr
	   jp        = jj
	   do i = j,  nc
	      if( a(jj,i) .gt. epslon ) go to 126
	   end do
	end do

!  all remaining elements are zero: matrix is singular.
!  report largest rhs element.

	jlast     = j - 1
	bmax      = 0.d0
	do jj = j, nr
	   bmax      = max( bmax, b(jj) )
	end do

	if( iwarn < nwgecp ) then
	   iwarn = iwarn + 1
	   write(lu_err,125) nr-jlast, bmax
125	   format( 'gecp: warning, rank def= ',i3, ' max res = ',1p,e13.4)
	endif
	go to 170

!  commute rows j and jp

126	continue
	bjp      = b(jp)
	b(jp)    = b(j)
	b(j)     = bjp

	do i = j, nc
	   ajp      = a(jp,i)
	   a(jp,i)  = a(j,i)
	   a(j,i)   = ajp
	end do 

	go to 110

!  commute columns j and ip

138	continue
	ipj       = ipiv(j)
	ipiv(j)   = ipiv(ip)
	ipiv(ip)  = ipj

	do jj = 1, nr
	   ajj       = a(jj,j)
	   a(jj,j)   = a(jj,ip)
	   a(jj,ip)  = ajj
	end do

	if( j == nr ) go to 170

!  perform elimination

	do jj = j+1, nr
	   amult     = a(jj,j) / a(j,j)
	   b(jj)     = b(jj) - amult * b(j)

	   do i  = j, nc
	      a(jj,i)   = a(jj,i) - amult * a(j,i)
	   end do
	end do

100	continue

!--------  back-substitution  ---------------------------------------------

170	continue
	a(jlast,1) = b(jlast) / a(jlast,jlast)

	do j  = jlast-1, 1, -1
	   do i   = j+1, jlast
	      b(j)       = b(j) - a(j,i) * a(i,1)
	   end do
	   a(j,1)     = b(j) / a(j,j)
	end do

!  set x

	do i  = 1, nc
	   x(i)      = 0.
	end do

	do j  = 1, jlast
	   x( ipiv(j) ) = a(j,1)
	end do

	return
	end
