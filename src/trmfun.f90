module radbelt_module 

   ! trmfun.for	1987

   implicit none 

   private 

   public :: trara1, trara2

   contains

!***********************************************************************
!*** trara1 finds particle fluxes for given energies, magnetic field ***
!*** strength and l-value. function trara2 is used to interpolate in ***
!*** b-l-space.                                                      ***
!***   input: descr(8)   header of specified trapped radition model  ***
!***          map(...)   map of trapped radition model               ***
!***                     (descr and map are explained at the begin   ***
!***                     of the main program model)                  ***
!***          n          number of energies                          ***
!***          e(n)       array of energies in mev                    ***
!***          fl         l-value                                     ***
!***          bb0        =b/b0  magnetic field strength normalized   ***
!***                     to field strength at magnetic equator       ***
!***  output: f(n)       decadic logarithm of integral fluxes in     ***
!***                     particles/(cm*cm*sec)                       ***
!***********************************************************************
subroutine trara1(descr,map,fl,bb0,e,f,n)
 
   real bb0 , e , e0 , e1 , e2 , escale , f , f0 , f1 , f2 , fistep , fl , fscale , xnl
   integer i0 , i1 , i2 , i3 , ie , l3 , map , n , nb , nl
   logical s0 , s1 , s2
   dimension e(n) , f(n) , map(*)
   integer descr(8)
   common /tra2  / fistep
   data f1 , f2/1.001 , 1.002/
!
   fistep = descr(7)/descr(2)
   escale = descr(4)
   fscale = descr(7)
   xnl = amin1(15.6,abs(fl))
   nl = xnl*descr(5)
   if ( bb0<1. ) bb0 = 1.
   nb = (bb0-1.)*descr(6)
!
! i2 is the number of elements in the flux map for the first energy.
! i3 is the index of the last element of the second energy map.
! l3 is the length of the map for the third energy.
! e1 is the energy of the first energy map (unscaled)
! e2 is the energy of the second energy map (unscaled)
!
   i1 = 0
   i2 = map(1)
   i3 = i2 + map(i2+1)
   l3 = map(i3+1)
   e1 = map(i1+2)/escale
   e2 = map(i2+2)/escale
!
! s0, s1, s2 are logical variables which indicate whether the flux for
! a particular e, b, l point has already been found in a previous call
! to function trara2. if not, s.. =.true.
!
   s1 = .true.
   s2 = .true.
!
!			energy loop
!
   do ie = 1 , n
!
! for each energy e(i) find the successive energies e0,e1,e2 in
! model map, which obey  e0 < e1 < e(i) < e2 .
!
      do while ( (e(ie)>e2) .and. (l3/=0) )
         i0 = i1
         i1 = i2
         i2 = i3
         i3 = i3 + l3
         l3 = map(i3+1)
         e0 = e1
         e1 = e2
         e2 = map(i2+2)/escale
         s0 = s1
         s1 = s2
         s2 = .true.
         f0 = f1
         f1 = f2
      enddo
!
! call trara2 to interpolate the flux-maps for e1,e2 in l-b/b0-
! space to find fluxes f1,f2 [if they have not already been
! calculated for a previous e(i)].
!
      if ( s1 ) f1 = trara2(map(i1+3),nl,nb)/fscale
      if ( s2 ) f2 = trara2(map(i2+3),nl,nb)/fscale
      s1 = .false.
      s2 = .false.
!
! finally, interpolate in energy.
!
      f(ie) = f1 + (f2-f1)*(e(ie)-e1)/(e2-e1)
      if ( f2<=0.0 ) then
         if ( i1/=0 ) then
!
! --------- special interpolation ---------------------------------
! if the flux for the second energy cannot be found (i.e. f2=0.0),
! and the zeroth energy map has been defined (i.e. i1 not equal 0),
! then interpolate using the flux maps for the zeroth and first
! energy and choose the minimum of this interpolations and the
! interpolation that was done with f2=0.
!
            if ( s0 ) f0 = trara2(map(i0+3),nl,nb)/fscale
            s0 = .false.
            f(ie) = amin1(f(ie),f0+(f1-f0)*(e(ie)-e0)/(e1-e0))
         endif
      endif
!
! the logarithmic flux is always kept greater or equal zero.
!
      f(ie) = amax1(f(ie),0.)
   enddo
end subroutine trara1

!*****************************************************************
!***  trara2 interpolates linearly in l-b/b0-map to obtain     ***
!***  the logarithm of integral flux at given l and b/b0.      ***
!***    input: map(..) is sub-map (for specific energy) of     ***
!***                   trapped radiation model map             ***
!***           il      scaled l-value                          ***
!***           ib      scaled b/b0-1                           ***
!***   output: trara2  scaled logarithm of particle flux       ***
!*****************************************************************
!***  see main program 'model' for explanation of map format   ***
!***  scaling factors.                                         ***
!***  the stepsize for the parameterization of the logarithm   ***
!***  of flux is obtained from 'common/tra2/'.                 ***
!*****************************************************************
function trara2(map,il,ib)

   real dfl , fincr1 , fincr2 , fistep , fkb , fkb1 , fkb2 , fkbj1 , fkbj2 , &
        fkbm , fll1 , fll2 , flog , flog1 , flog2 , flogm ,   &
        fnb , fnl , sl1 , sl2
   real trara2
   integer i1 , i2 , ib , il , itime , j1 , j2 , kt , l1 , l2 , map(*)

   common /tra2  / fistep
   integer :: spag_nextblock_1

   spag_nextblock_1 = 1
   main: do
      select case (spag_nextblock_1)
      case (1)
         fnl = il
         fnb = ib
         itime = 0
         i2 = 0
         do
            !
            ! find consecutive sub-sub-maps for scaled l-values ls1,ls2,
            ! with il less or equal ls2.  l1,l2 are lengths of sub-sub-maps.
            ! i1,i2 are indeces of first elements minus 1.
            !
            l2 = map(i2+1)
            if ( map(i2+2)<=il ) then
               i1 = i2
               l1 = l2
               i2 = i2 + l2
               !
               ! if sub-sub-maps are empty, i. e. length less 4, than trara2=0
               !
            elseif ( (l1<4) .and. (l2<4) ) then
               trara2 = 0.
               return
            else
               !
               ! if flog2 less flog1, than ls2 first map and ls1 second map
               !
               if ( map(i2+3)<=map(i1+3) ) exit
               spag_nextblock_1 = 3
               cycle main
            endif
         enddo
         spag_nextblock_1 = 2
      case (2)
         kt = i1
         i1 = i2
         i2 = kt
         kt = l1
         l1 = l2
         l2 = kt
         spag_nextblock_1 = 3
      case (3)
         !
         ! determine interpolate in scaled l-value
         !
         fll1 = map(i1+2)
         fll2 = map(i2+2)
         dfl = (fnl-fll1)/(fll2-fll1)
         flog1 = map(i1+3)
         flog2 = map(i2+3)
         fkb1 = 0.
         fkb2 = 0.
         if ( l1>=4 ) then
            !
            ! b/b0 loop
            !
            do j2 = 4 , l2
               fincr2 = map(i2+j2)
               if ( fkb2+fincr2>fnb ) goto 10
               fkb2 = fkb2 + fincr2
               flog2 = flog2 - fistep
            enddo
            itime = itime + 1
            if ( itime==1 ) then
               spag_nextblock_1 = 2
               cycle main
            endif
            trara2 = 0.
            return
 10         if ( itime/=1 ) then
               if ( j2==4 ) then
                  spag_nextblock_1 = 4
                  cycle main
               endif
               sl2 = flog2/fkb2
               do j1 = 4 , l1
                  fincr1 = map(i1+j1)
                  fkb1 = fkb1 + fincr1
                  flog1 = flog1 - fistep
                  fkbj1 = ((flog1/fistep)*fincr1+fkb1)/((fincr1/fistep)*sl2+1.)
                  if ( fkbj1<=fkb1 ) goto 15
               enddo
               if ( fkbj1<=fkb2 ) then
                  trara2 = 0.
                  return
               endif
 15            if ( fkbj1<=fkb2 ) then
                  fkbm = fkbj1 + (fkb2-fkbj1)*dfl
                  flogm = fkbm*sl2
                  flog2 = flog2 - fistep
                  fkb2 = fkb2 + fincr2
                  sl1 = flog1/fkb1
                  sl2 = flog2/fkb2
                  spag_nextblock_1 = 5
                  cycle main
               else
                  fkb1 = 0.
               endif
            endif
            fkb2 = 0.
         endif
         j2 = 4
         fincr2 = map(i2+j2)
         flog2 = map(i2+3)
         flog1 = map(i1+3)
         spag_nextblock_1 = 4
      case (4)
         flogm = flog1 + (flog2-flog1)*dfl
         fkbm = 0.
         fkb2 = fkb2 + fincr2
         flog2 = flog2 - fistep
         sl2 = flog2/fkb2
         if ( l1<4 ) then
            fincr1 = 0.
            sl1 = -900000.
            spag_nextblock_1 = 6
            cycle main
         else
            j1 = 4
            fincr1 = map(i1+j1)
            fkb1 = fkb1 + fincr1
            flog1 = flog1 - fistep
            sl1 = flog1/fkb1
         endif
         spag_nextblock_1 = 5
      case (5)
         do while ( sl1>=sl2 )
            fkbj2 = ((flog2/fistep)*fincr2+fkb2)/((fincr2/fistep)*sl1+1.)
            fkb = fkb1 + (fkbj2-fkb1)*dfl
            flog = fkb*sl1
            if ( fkb>=fnb ) then
               spag_nextblock_1 = 7
               cycle main
            endif
            fkbm = fkb
            flogm = flog
            if ( j1>=l1 ) then
               trara2 = 0.
               return
            else
               j1 = j1 + 1
               fincr1 = map(i1+j1)
               flog1 = flog1 - fistep
               fkb1 = fkb1 + fincr1
               sl1 = flog1/fkb1
            endif
         enddo
         spag_nextblock_1 = 6
      case (6)
         fkbj1 = ((flog1/fistep)*fincr1+fkb1)/((fincr1/fistep)*sl2+1.)
         fkb = fkbj1 + (fkb2-fkbj1)*dfl
         flog = fkb*sl2
         if ( fkb<fnb ) then
            fkbm = fkb
            flogm = flog
            if ( j2>=l2 ) then
               trara2 = 0.
               return
            else
               j2 = j2 + 1
               fincr2 = map(i2+j2)
               flog2 = flog2 - fistep
               fkb2 = fkb2 + fincr2
               sl2 = flog2/fkb2
               spag_nextblock_1 = 5
               cycle main
            endif
         endif
         spag_nextblock_1 = 7
      case (7)
         if ( fkb<fkbm+1.e-10 ) then
            trara2 = 0.
         else
            trara2 = flogm + (flog-flogm)*((fnb-fkbm)/(fkb-fkbm))
            trara2 = amax1(trara2,0.)
            return
         endif
         exit main
      end select
   enddo main
end function trara2

end module radbelt_module