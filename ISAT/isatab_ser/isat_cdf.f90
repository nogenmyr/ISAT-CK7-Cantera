!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

module isat_cdf

!  Routines for forming and outputting cumulative distribution functions (CFDs).
!  F(x) denotes the CDF of the random variable X:  F(x) = Prob{ X < x }.

!  Usage:
!    call isat_cdf_init  to initialize the data structure (of type cdf_type).
!    call isat_cdf_add  and/or  isat_cdf_add1  to add samples to the CDF.
!    call isat_cdf_op  to output the CDF
!    call isat_cdf_reset  to reset the variables in the data structure to their initial values.
!    call isat_cdf_kill  to delete the data structure 

!  Method:
!    A histogram, consisting of nbin bins in the interval (x_lower, x_upper),
!    is constructed from the samples of X (through calls to  isat_cdf_add/1).
!    The bins have a logarithmic spacing if if_log=.true.; otherwise they have
!    equal spacing.  A sample consists of a value of X and a non-negative 
!    numerical weight.  Sample values X < x_lower and X > x_upper are placed in the 
!    first and last bins, respectively.  The CDF is formed from this histogram
!    when required for output.

!    If op_inc > 0, then the CDF is output to the file
!    op_file each time (approximately) op_inc samples are added.
!    isat_cdf_op can also be called to output the CDF.

!  Output:
!    The output file  op_file  consists of two numbers per line, 
!    with the following lines:

!      x_min    x_max
!      samples  wt_sum
!      xr(i)    cdf(i)

!    where: 
!       x_min  is the smallest sample value
!       x_max  is the largest sample value
!       sample  is the total number of samples
!       wt_sum  is the sum of the numerical weights of the samples
!       xr(:)  is a vector of values of x (in increasing order)
!       cdf(:) ia a vector containing the corresponding values of F( x(:) ) 

!    (The number of lines may be less than nbin+2.)

use isat_abort_m
use isat_prec

implicit none

type :: cdf_type    !-------------------------------------

  logical             :: initialized          ! indicate whether initialized
  logical             :: on                   ! =.true. to add to CDF
  integer             :: nbin                 ! number of bins
  integer             :: op_inc               ! output every op_inc samples 
  real(k_xf)          :: samples              ! number of samples
  real(k_xf)          :: op_next              ! output after op_next samples
  real(k_xf)          :: x_lower, x_upper     ! range is (x_lower,x_upper)
  real(k_xf)          :: a, b                 ! bin number is: j = a + b * x
  real(k_xf)          :: x_min, x_max         ! smallest and largest sample values
  real(k_xf)          :: wt_sum               ! sum of numerical weights
  real(k_xf), pointer :: cdf(:), pdf(:), xr(:)! cdf, pdf, and x at right of bin
  logical             :: log                  ! = .true. for log scale
  character(30)       :: op_file              ! output file name

end type cdf_type


contains  !--------------------------------------------------------------------------------------

  subroutine isat_cdf_init( nbin, x_lower, x_upper, if_log, op_file, op_inc, cdf )  !------------

! initialize data structure cdf

! input:
!   nbin    - number of bins  (nbin >= 2)
!   x_lower - lower limit of range  (most be strictly positive if if_log=.true.)
!   x_upper - upper limit of range  (x_upper > x_lower)
!   if_log  - = .true. if logarithmic bin spacing to be used (otherwise equal spacing)
!   op_file - file name for output
!   op_inc  - output every op_inc samples (set op_inc = -1 to suppress output)
! output:
!   cdf     - initialized data structure for CDF

  integer, intent(in)     :: nbin, op_inc
  real(k_xf), intent(in)  :: x_lower, x_upper
  logical                 :: if_log
  character(*)            :: op_file
  type(cdf_type), pointer :: cdf

  integer :: i

!  check input
  if( nbin < 2 ) call isat_abort('isat_cdf_init', 1, &
     mess = 'nbin must be at least 2', isv=nbin )

  if( x_lower >= x_upper ) call isat_abort('isat_cdf_init', 2, &
     mess = 'x_lower must be less than x_upper', rvar = (/ x_lower, x_upper /) )

  if( if_log  .and.  x_lower <= 0.d0 ) call isat_abort('isat_cdf_init', 3, &
     mess = 'x_lower must be positive when using log scale', rsv=x_lower ) 

  if( .not.associated(cdf) ) allocate( cdf )
  allocate( cdf%cdf(nbin) )
  allocate( cdf%pdf(nbin) )
  allocate( cdf%xr(nbin) )

  cdf%initialized = .true.
  cdf%on      = .true.
  cdf%cdf     = 0.d0
  cdf%pdf     = 0.d0
  cdf%samples = 0.d0
  cdf%wt_sum  = 0.d0
  cdf%op_inc  = op_inc
  cdf%op_next = op_inc
  if( op_inc < 0 ) cdf%op_next = huge( cdf%op_next )
  cdf%op_file = op_file
  cdf%nbin    = nbin
  cdf%x_lower = x_lower
  cdf%x_upper = x_upper
  cdf%log     = if_log
  cdf%x_max   = -huge(1.d0)
  cdf%x_min   =  huge(1.d0)


  if( .not.if_log ) then   ! linear bin spacing
 !  given sample X, the bin number j is: j = cdf%a + cdf%b * X
     cdf%b = float(nbin)/(x_upper-x_lower)
     cdf%a = 1. - x_lower * cdf%b

	 do i = 1, nbin	! right limit of bin
        cdf%xr(i) = x_lower + float(i) / cdf%b
	 end do
  else
 !  given sample X, the bin number j is: j = cdf%a + cdf%b * log(X)
     cdf%b = float(nbin)/( log(x_upper) - log(x_lower) )
     cdf%a = 1. - log(x_lower) * cdf%b

     do i = 1, nbin	! right limit of bin
        cdf%xr(i) = exp( log(x_lower) + float(i) / cdf%b )
     end do
  endif

  return
  end subroutine isat_cdf_init

  subroutine isat_cdf_add1( cdf, x )	!-------------------------------------

!  add a single unit-weight sample to the cdf

! input:
!   cdf - data structure for CDF
!   x   - sample value

! output:
!   cdf - incremented data structure

  type(cdf_type), pointer :: cdf
  real(k_xf), intent(in)  :: x

  real(k_xf) :: xx(1), wt(1)

  xx(1) = x
  wt(1) = 1.d0

  call isat_cdf_add( cdf, 1, xx, wt )

  return
  end subroutine isat_cdf_add1

  subroutine isat_cdf_add( cdf, n, x, wt )	!-------------------------------------

!  add n samples to the cdf

! input:
!   cdf - data structure for CDF
!   n   - number of samples to be added
!   x   - vector of n samples
!   wt  - vector of corresponding numerical weights

! output:
!   cdf - incremented data structure

  type(cdf_type), pointer :: cdf
  integer, intent(in)     :: n
  real(k_xf), intent(in)  :: x(n), wt(n)

  integer i, j

  if( .not.cdf%log ) then   ! linear bin spacing
     do i = 1, n
        j = cdf%a + cdf%b * x(i)
	    j = max( 1, min( j, cdf%nbin ) )
	    cdf%pdf(j) = cdf%pdf(j) + wt(i)
		cdf%x_max  = max( x(i), cdf%x_max ) 
		cdf%x_min  = min( x(i), cdf%x_min ) 
     end do
  else
     do i = 1, n
        j = cdf%a + cdf%b * log( max( x(i), cdf%x_lower ) )
	    j = max( 1, min( j, cdf%nbin ) )
	    cdf%pdf(j) = cdf%pdf(j) + wt(i)
		cdf%x_max  = max( x(i), cdf%x_max ) 
		cdf%x_min  = min( x(i), cdf%x_min ) 
     end do
  endif

  cdf%samples = cdf%samples + n

  if( cdf%samples > cdf%op_next ) call isat_cdf_op( cdf )
 
  return
  end subroutine isat_cdf_add

  subroutine isat_cdf_f2x( cdf, F, x )	!-------------------------------------

!  For the CDF F(x), given a value of F, return the corresponding value of x. 

  type(cdf_type), pointer :: cdf  ! given CDF F(x)
  real(k_xf), intent(in)  :: F    ! given value of F (0 < F < 1)
  real(k_xf), intent(out) :: x    ! corresponding value of x

  integer    :: j
  real(k_xf) :: jx

  if( cdf%samples < 1.d0 ) then
     x = 0.d0  !  no samples, no valid output
     return
  endif
  
  call isat_cdf_form( cdf )  ! form CDF

  do j = 1, cdf%nbin    !  find smallest j: F(j) >=F
     if( cdf%cdf(j) >= F ) exit
  end do

!  return imposed limits if x is in first or last bin

  if( j == 1 ) then
     x = cdf%x_lower
	 return
  elseif( j >= cdf%nbin ) then
     x = cdf%x_upper
	 return
  endif

  jx = j-1 + (F-cdf%cdf(j-1))/(cdf%cdf(j)-cdf%cdf(j-1))
  x  = (jx - cdf%a) / cdf%b

  if( cdf%log ) x = exp( x )

  return
  end subroutine isat_cdf_f2x

  subroutine isat_cdf_op( cdf )	!-------------------------------------

  !  output cdf

  !  The output file  cdf%op_file  consists of two numbers per line, 
  !  with the following lines:

  !      x_min    x_max
  !      samples  wt_sum
  !      xr(i)    cdf(i)  - for all i at which the CDF increases

  type(cdf_type), pointer :: cdf

  real(k_xf) :: cdf_last
  integer    :: i, lu

  call isat_cdf_form( cdf )

  call isat_lu( lu )
  open( lu, file = trim(cdf%op_file), err=100 )

200     format(1p,2e25.12)
  write(lu,200) cdf%x_min, cdf%x_max
  write(lu,200) cdf%samples, cdf%wt_sum

!  write out CDF only where it is greater than the last value
  cdf_last = 0.d0
  do i = 1, cdf%nbin
     if( cdf%cdf(i)-cdf_last > 1.d-11 *( cdf%cdf(i) + cdf_last ) ) then
	    write(lu,200) cdf%xr(i), cdf%cdf(i)
		cdf_last = cdf%cdf(i)
	 endif
  end do

  close( lu )

  if( cdf%op_inc > 0 ) cdf%op_next = cdf%samples + cdf%op_inc

  return

100     call isat_abort('isat_cdf_op', 1, mess='error opening lu', isv = lu )

  end subroutine isat_cdf_op

  subroutine isat_cdf_form( cdf )	!-------------------------------------

!  form cdf from pdf

! input/output:
!   cdf - data structure for CDF

  type(cdf_type), pointer :: cdf

  integer    :: i

  cdf%cdf(1) = cdf%pdf(1)
  do i = 2, cdf%nbin
     cdf%cdf(i) = cdf%cdf(i-1) + cdf%pdf(i)
  end do

  cdf%wt_sum = cdf%cdf( cdf%nbin )

  cdf%cdf = cdf%cdf / cdf%wt_sum

  return
  end subroutine isat_cdf_form

  subroutine isat_cdf_reset( cdf )	!-------------------------------------

! restore data structure to its initial state

  type(cdf_type), pointer :: cdf

  cdf%samples = 0.d0
  cdf%op_next = cdf%op_inc
  if( cdf%op_inc < 0 ) cdf%op_next = huge( cdf%op_next )
  cdf%wt_sum = 0.d0

  cdf%cdf     = 0.
  cdf%pdf     = 0.
  cdf%x_max   = -huge(1.d0)
  cdf%x_min   =  huge(1.d0)

  return
  end subroutine isat_cdf_reset

  subroutine isat_cdf_kill( cdf )	!-------------------------------------

!  delete the data structure  cdf

  type(cdf_type), pointer :: cdf

  if( .not.associated(cdf) ) return
  if( associated(cdf%cdf) ) deallocate( cdf%cdf )
  if( associated(cdf%pdf) ) deallocate( cdf%pdf )
  if( associated(cdf%xr ) ) deallocate( cdf%xr  )
  deallocate( cdf )

  return
  end subroutine isat_cdf_kill

  subroutine isat_cdf_write( lu, cdf )   !--------------------------------

!  write CDF for checkpointing

  integer, intent(in)     :: lu  ! logical unit for writing
  type(cdf_type), pointer :: cdf ! CDF to be written

  write(lu) cdf%initialized
  if( .not.cdf%initialized ) return

  write(lu) cdf%nbin, cdf%x_lower, cdf%x_upper, cdf%log, cdf%op_file, cdf%op_inc
  write(lu) cdf%on
  write(lu) cdf%samples, cdf%x_min, cdf%x_max, cdf%wt_sum
  write(lu) cdf%cdf(:), cdf%pdf(:), cdf%xr(:)

  return
  end subroutine isat_cdf_write

  subroutine isat_cdf_read( lu, cdf )   !--------------------------------

!  Read checkpointed CDF.
!  If CDF exists on entry, and there is a non-trivial CDF to read, 
!  then the CDF is killed and then re-initialized.

  integer, intent(in)     :: lu  ! logical unit for writing
  type(cdf_type), pointer :: cdf ! CDF to be written

  logical       :: initialized, log
  integer       :: nbin, op_inc
  real(k_xf)    :: x_lower, x_upper
  character(30) :: op_file 

  read(lu,err=100,end=100) initialized
  if( .not.initialized ) return

  if( associated( cdf ) ) call isat_cdf_kill( cdf )

  read(lu,err=110,end=110) nbin, x_lower, x_upper, log, op_file, op_inc
  
  call isat_cdf_init( nbin, x_lower, x_upper, log, op_file, op_inc, cdf )
  
  read(lu,err=120,end=120) cdf%on
  read(lu,err=130,end=130) cdf%samples, cdf%x_min, cdf%x_max, cdf%wt_sum
  read(lu,err=140,end=140) cdf%cdf(:), cdf%pdf(:), cdf%xr(:)

  return

100  call isat_abort('isat_cdf_read',1, mess='error reading initialized')
110  call isat_abort('isat_cdf_read',1, mess='error reading nbin')
120  call isat_abort('isat_cdf_read',1, mess='error reading on')
130  call isat_abort('isat_cdf_read',1, mess='error reading samples')
140  call isat_abort('isat_cdf_read',1, mess='error reading cdf')

  end subroutine isat_cdf_read

end module isat_cdf