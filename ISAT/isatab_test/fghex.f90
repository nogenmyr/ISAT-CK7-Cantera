subroutine fghex( need, nx, x, nf, nh, iusr, rusr, f, g, h )

!==========================================================================
!-------  INPUT
!	need	- integer array indicating function values to be returned
!	need(1)	= 1, return  f	
!	need(2)	= 1, return  g = df/dx
!	nx	- number of components of x  ( nx >= 1 )
!	x	- components of x
!	nf	- number of components of f  ( nf >= 1 )
!	nh	- number of components of h  ( nh >= 0 )
!	iusr	- integer array: iusr(1) = idtab, iusr(2) = ifsin
!	rusr	- not used
!-------  OUTPUT
!	f	- f(x)
!	g	- g(x)
!	h	- not used
!==========================================================================

   use isat_rnu
   implicit none
   integer, intent(in) :: need(3), nx, nf, nh, iusr(*)
   double precision, intent(in)  :: x(nx), rusr(*)
   double precision, intent(out) :: f(nf), g(nf,nx), h(nh)

   integer, parameter :: ntab = 3
   integer, save :: ifst(ntab) = 0 
   integer :: idtab, ifsin, ngh, j, is1, is2
   double precision, allocatable, save, target :: a1(:,:), a2(:,:), a3(:,:)
   double precision, pointer :: a(:,:)

   idtab = iusr(1)
   ifsin = iusr(2)

   if( ifst(idtab) == 0 ) then
      call rnuget( is1, is2 )
      call rnuput( 9753, 8642 )
      ifst(idtab) = 1
      if( idtab == 1 ) then
         allocate( a1(nf,nx) )
         call rnu( a1 )
      elseif( idtab == 2 ) then
         allocate( a2(nf,nx) )
         call rnu( a2 )
      elseif( idtab == 3 ) then
         allocate( a3(nf,nx) )
         call rnu( a3 )
      else
         write(0,*)'fghex: bad idtab=', idtab
         stop
      endif
      call rnuput( is1, is2 )
   endif

!-------  all calls:  evaluate f(x)  -----------------------------------

   if( idtab == 1 ) a => a1(1:nf,1:nx)
   if( idtab == 2 ) a => a2(1:nf,1:nx)
   if( idtab == 3 ) a => a3(1:nf,1:nx)

   f = matmul( a, x )
   g = a
   h = 0.

   if( ifsin == 1 ) then
      do j = 1, nf
         g(j,:) = a(j,:) * cos( f(j) )
      enddo
      f = sin( f )
   endif

   ngh = min(nf,nh)
   h(1:ngh) = f(1:ngh)

   if( need(1) == 0 ) f = -1.
   if( need(2) == 0 ) g = -1.
   if( need(3) == 0 ) h = -1.

! build in delays if so desired

    is1 = 0
	if( need(1) > 0 ) is1 = 1000
	if( need(2) > 0 ) is1 = 10000
	is1 = 0 !!!!  over-ridden
	do j = 1, is1
	   is2 = exp(1. + 1./j)
	end do

   return
   end
