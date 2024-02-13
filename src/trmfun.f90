!*****************************************************************************************
!>
!  Trapped radiation model.
!
!### History
!  * Based on: `trmfun.for` 1987

module trmfun_module

   use radbelt_kinds_module

   implicit none

   private

   character(len=10),dimension(4),parameter :: mname = [ 'ae8min.asc' , &
                                                         'ae8max.asc' , &
                                                         'ap8min.asc' , &
                                                         'ap8max.asc'] !! data files available

   type,public :: trm_type
      !! main class for the `aep8` model
      private

      character(len=:),allocatable :: aep8_dir !! directory containing the data files

      ! data read from the files:
      character(len=:),allocatable :: file_loaded !! the file that has been loaded
      integer,dimension(8) :: ihead = 0
      integer,dimension(:),allocatable :: map

      real(wp) :: fistep = 0.0_wp !! the stepsize for the parameterization of the logarithm of flux.
                                  !! formerly stored in common block `tra2`

      ! formerly saved variables in trara1:
      real(wp) :: f1 = 1.001_wp
      real(wp) :: f2 = 1.002_wp

      contains
      private
      procedure,public :: aep8 !! main routine
      procedure,public :: trara1, trara2 !! low-level routine
      procedure,public :: set_data_file_dir, get_data_file_dir
   end type trm_type

   contains

!*****************************************************************************************
!>
!  Set the directory containing the data files.

   subroutine set_data_file_dir(me,dir)
      class(trm_type),intent(inout) :: me
      character(len=*),intent(in) :: dir
      me%aep8_dir = trim(dir)
   end subroutine set_data_file_dir

!*****************************************************************************************
!>
!  Get the directory containing the data files.

   function get_data_file_dir(me) result(dir)
      class(trm_type),intent(in) :: me
      character(len=:),allocatable :: dir
      if (allocated(me%aep8_dir)) then
         dir = trim(me%aep8_dir) // '/'
      else
         dir = 'data/aep8/' ! default
      end if
   end function get_data_file_dir

!*****************************************************************************************
!>
!  Main wrapper for the radiation model.
!  Reads the coefficient file and calls the low-level routine.

   subroutine aep8(me,e,l,bb0,imname,flux)

      class(trm_type),intent(inout) :: me

      real(wp),intent(in) :: e
      real(wp),intent(in) :: l
      real(wp),intent(in) :: bb0
      integer,intent(in) :: imname !! which model to load (index in `mname` array)
      real(wp),intent(out) :: flux

      real(wp) :: ee(1), f(1) !! temp variables
      integer :: i , ierr, iuaeap , nmap
      character(len=:),allocatable :: name
      logical :: load_file

      name = me%get_data_file_dir() // trim(mname(Imname)) ! the file to load

      ! JW : do we need to reset some or all of these ?
      me%fistep = 0.0_wp
      me%f1 = 1.001_wp
      me%f2 = 1.002_wp

      ! check to see if this file has already been loaded
      ! [the class can store one file at a time]
      load_file = .true.
      if (allocated(me%file_loaded)) then
         if (name == me%file_loaded) load_file = .false.
      end if

      if (load_file) then
         open (newunit = iuaeap,file=name,status='OLD',iostat=ierr,form='FORMATTED')
         if ( ierr/=0 ) then
            error stop 'error reading '//trim(name)
         end if
         read (iuaeap,'(1X,12I6)') me%ihead
         nmap = me%ihead(8)
         allocate(me%map(nmap))
         read (iuaeap,'(1X,12I6)') (me%map(i),i=1,nmap)
         close (iuaeap)
         me%file_loaded = trim(name)
      end if

      ee(1) = e
      call me%trara1(me%ihead,me%map,L,Bb0,ee,f,1)
      flux = f(1)
      IF ( Flux>0.0_wp ) Flux = 10.0_wp**Flux

   end subroutine aep8
!*****************************************************************************************

!*****************************************************************************************
!>
!  [[trara1]] finds particle fluxes for given energies, magnetic field
!  strength and l-value. function [[trara2]] is used to interpolate in
!  b-l-space.

   subroutine trara1(me,descr,map,fl,bb0,e,f,n)

   class(trm_type),intent(inout) :: me
   integer,intent(in) :: n !! number of energies
   integer,intent(in) :: descr(8) !! header of specified trapped radition model
   real(wp),intent(in) :: e(n) !! array of energies in mev
   real(wp),intent(in) :: fl !! l-value
   real(wp),intent(in) :: bb0 !! =b/b0  magnetic field strength normalized
                                 !! to field strength at magnetic equator
   integer,intent(in) :: map(*) !! map of trapped radition model
                                !! (descr and map are explained at the begin
                                !! of the main program model)
   real(wp),intent(out) :: f(n) !! decadic logarithm of integral fluxes in
                                !! particles/(cm*cm*sec)

   real(wp) :: e0 , e1 , e2 , escale , f0 , fscale , xnl
   real(wp) :: bb0_ !! local copy of `bb0`. in the original code
                    !! this was modified by this routine.
                    !! added this so `bb0` could be `intent(in)`
   integer :: i0 , i1 , i2 , i3 , ie , l3 , nb , nl
   logical :: s0 , s1 , s2

   e0 = 0.0_wp  ! to avoid -Wmaybe-uninitialized warnings
   f0 = 0.0_wp  ! to avoid -Wmaybe-uninitialized warnings
   i0 = 0       ! to avoid -Wmaybe-uninitialized warnings
   s0 = .false. ! to avoid -Wmaybe-uninitialized warnings  -- but not sure what default value here should be !  -JW

   bb0_ = bb0
   me%fistep = descr(7)/descr(2)
   escale = descr(4)
   fscale = descr(7)
   xnl = min(15.6_wp,abs(fl))
   nl = xnl*descr(5)
   if ( bb0_<1.0_wp ) bb0_ = 1.0_wp
   nb = (bb0_-1.0_wp)*descr(6)

   ! i2 is the number of elements in the flux map for the first energy.
   ! i3 is the index of the last element of the second energy map.
   ! l3 is the length of the map for the third energy.
   ! e1 is the energy of the first energy map (unscaled)
   ! e2 is the energy of the second energy map (unscaled)
   i1 = 0
   i2 = map(1)
   i3 = i2 + map(i2+1)
   l3 = map(i3+1)
   e1 = map(i1+2)/escale
   e2 = map(i2+2)/escale

   ! s0, s1, s2 are logical variables which indicate whether the flux for
   ! a particular e, b, l point has already been found in a previous call
   ! to function trara2. if not, s.. =.true.
   s1 = .true.
   s2 = .true.

   ! energy loop

   do ie = 1 , n

      ! for each energy e(i) find the successive energies e0,e1,e2 in
      ! model map, which obey  e0 < e1 < e(i) < e2 .

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
         f0 = me%f1
         me%f1 = me%f2
      enddo

      ! call trara2 to interpolate the flux-maps for e1,e2 in l-b/b0-
      ! space to find fluxes f1,f2 [if they have not already been
      ! calculated for a previous e(i)].

      if ( s1 ) me%f1 = me%trara2(map(i1+3),nl,nb)/fscale
      if ( s2 ) me%f2 = me%trara2(map(i2+3),nl,nb)/fscale
      s1 = .false.
      s2 = .false.

      ! finally, interpolate in energy.

      f(ie) = me%f1 + (me%f2-me%f1)*(e(ie)-e1)/(e2-e1)
      if ( me%f2<=0.0_wp ) then
         if ( i1/=0 ) then
            ! --------- special interpolation ---------------------------------
            ! if the flux for the second energy cannot be found (i.e. f2=0.0),
            ! and the zeroth energy map has been defined (i.e. i1 not equal 0),
            ! then interpolate using the flux maps for the zeroth and first
            ! energy and choose the minimum of this interpolations and the
            ! interpolation that was done with f2=0.
            if ( s0 ) f0 = me%trara2(map(i0+3),nl,nb)/fscale
            s0 = .false.
            f(ie) = min(f(ie),f0+(me%f1-f0)*(e(ie)-e0)/(e1-e0))
         endif
      endif

      ! the logarithmic flux is always kept greater or equal zero.

      f(ie) = max(f(ie),0.0_wp)
   enddo
end subroutine trara1

!*****************************************************************************************
!>
!  [[trara2]] interpolates linearly in l-b/b0-map to obtain
!  the logarithm of integral flux at given l and b/b0.
!
!### Note
!  see main program 'model' for explanation of map format
!  scaling factors.

function trara2(me,map,il,ib)

   class(trm_type),intent(inout) :: me
   integer,intent(in) :: map(*) !! is sub-map (for specific energy) of
                                !! trapped radiation model map
   integer,intent(in) :: il !! scaled l-value
   integer,intent(in) :: ib !! scaled b/b0-1
   real(wp) :: trara2 !! scaled logarithm of particle flux

   real(wp) :: dfl , fincr1 , fincr2 , fistep , fkb , fkb1 , fkb2 , fkbj1 , fkbj2 , &
               fkbm , fll1 , fll2 , flog , flog1 , flog2 , flogm , &
               fnb , fnl , sl1 , sl2
   integer :: i1 , i2 , itime , j1 , j2 , kt , l1 , l2
   integer :: itask

   fistep = me%fistep

   itask = 1
   main: do
      select case (itask)
      case (1)
         fnl = il
         fnb = ib
         itime = 0
         i2 = 0
         do

            ! find consecutive sub-sub-maps for scaled l-values ls1,ls2,
            ! with il less or equal ls2.  l1,l2 are lengths of sub-sub-maps.
            ! i1,i2 are indeces of first elements minus 1.

            l2 = map(i2+1)
            if ( map(i2+2)<=il ) then
               i1 = i2
               l1 = l2
               i2 = i2 + l2

               ! if sub-sub-maps are empty, i. e. length less 4, than trara2=0

            elseif ( (l1<4) .and. (l2<4) ) then
               trara2 = 0.0_wp
               return
            else

               ! if flog2 less flog1, than ls2 first map and ls1 second map

               if ( map(i2+3)<=map(i1+3) ) exit
               itask = 3
               cycle main
            endif
         enddo
         itask = 2
      case (2)
         kt = i1
         i1 = i2
         i2 = kt
         kt = l1
         l1 = l2
         l2 = kt
         itask = 3
      case (3)

         ! determine interpolate in scaled l-value

         fll1 = map(i1+2)
         fll2 = map(i2+2)
         dfl = (fnl-fll1)/(fll2-fll1)
         flog1 = map(i1+3)
         flog2 = map(i2+3)
         fkb1 = 0.0_wp
         fkb2 = 0.0_wp
         if ( l1>=4 ) then

            ! b/b0 loop

            do j2 = 4 , l2
               fincr2 = map(i2+j2)
               if ( fkb2+fincr2>fnb ) goto 10
               fkb2 = fkb2 + fincr2
               flog2 = flog2 - fistep
            enddo
            itime = itime + 1
            if ( itime==1 ) then
               itask = 2
               cycle main
            endif
            trara2 = 0.0_wp
            return
 10         if ( itime/=1 ) then
               if ( j2==4 ) then
                  itask = 4
                  cycle main
               endif
               sl2 = flog2/fkb2
               do j1 = 4 , l1
                  fincr1 = map(i1+j1)
                  fkb1 = fkb1 + fincr1
                  flog1 = flog1 - fistep
                  fkbj1 = ((flog1/fistep)*fincr1+fkb1)/((fincr1/fistep)*sl2+1.0_wp)
                  if ( fkbj1<=fkb1 ) goto 15
               enddo
               if ( fkbj1<=fkb2 ) then
                  trara2 = 0.0_wp
                  return
               endif
 15            if ( fkbj1<=fkb2 ) then
                  fkbm = fkbj1 + (fkb2-fkbj1)*dfl
                  flogm = fkbm*sl2
                  flog2 = flog2 - fistep
                  fkb2 = fkb2 + fincr2
                  sl1 = flog1/fkb1
                  sl2 = flog2/fkb2
                  itask = 5
                  cycle main
               else
                  fkb1 = 0.0_wp
               endif
            endif
            fkb2 = 0.0_wp
         endif
         j2 = 4
         fincr2 = map(i2+j2)
         flog2 = map(i2+3)
         flog1 = map(i1+3)
         itask = 4
      case (4)
         flogm = flog1 + (flog2-flog1)*dfl
         fkbm = 0.0_wp
         fkb2 = fkb2 + fincr2
         flog2 = flog2 - fistep
         sl2 = flog2/fkb2
         if ( l1<4 ) then
            fincr1 = 0.0_wp
            sl1 = -900000.0_wp
            itask = 6
            cycle main
         else
            j1 = 4
            fincr1 = map(i1+j1)
            fkb1 = fkb1 + fincr1
            flog1 = flog1 - fistep
            sl1 = flog1/fkb1
         endif
         itask = 5
      case (5)
         do while ( sl1>=sl2 )
            fkbj2 = ((flog2/fistep)*fincr2+fkb2)/((fincr2/fistep)*sl1+1.0_wp)
            fkb = fkb1 + (fkbj2-fkb1)*dfl
            flog = fkb*sl1
            if ( fkb>=fnb ) then
               itask = 7
               cycle main
            endif
            fkbm = fkb
            flogm = flog
            if ( j1>=l1 ) then
               trara2 = 0.0_wp
               return
            else
               j1 = j1 + 1
               fincr1 = map(i1+j1)
               flog1 = flog1 - fistep
               fkb1 = fkb1 + fincr1
               sl1 = flog1/fkb1
            endif
         enddo
         itask = 6
      case (6)
         fkbj1 = ((flog1/fistep)*fincr1+fkb1)/((fincr1/fistep)*sl2+1.0_wp)
         fkb = fkbj1 + (fkb2-fkbj1)*dfl
         flog = fkb*sl2
         if ( fkb<fnb ) then
            fkbm = fkb
            flogm = flog
            if ( j2>=l2 ) then
               trara2 = 0.0_wp
               return
            else
               j2 = j2 + 1
               fincr2 = map(i2+j2)
               flog2 = flog2 - fistep
               fkb2 = fkb2 + fincr2
               sl2 = flog2/fkb2
               itask = 5
               cycle main
            endif
         endif
         itask = 7
      case (7)
         if ( fkb<fkbm+1.0e-10_wp ) then
            trara2 = 0.0_wp
         else
            trara2 = flogm + (flog-flogm)*((fnb-fkbm)/(fkb-fkbm))
            trara2 = max(trara2,0.0_wp)
            return
         endif
         exit main
      end select
   enddo main

end function trara2

end module trmfun_module