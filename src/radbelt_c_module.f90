!*****************************************************************************************
!>
!  Experimental C interface to the radbelt module.
!
!@note Currently, this makes use of a global variable to store the state,
!      so should not be considered threadsafe.

module radbelt_c_module

use iso_c_binding, only: c_double, c_int, c_char, c_null_char, &
                         c_intptr_t, c_ptr, c_loc, c_f_pointer, c_null_ptr, c_associated
use radbelt_module, only: radbelt_type

implicit none

integer,parameter :: STR_LEN = 256 !! string length for paths

contains

! function c2f_str(cstr) result(fstr)
!     character(kind=c_char,len=1),dimension(*),intent(in) :: cstr
!     character(len=:),allocatable :: fstr
!     integer :: i
!     fstr = ''
!     i = 0
!     do
!         i = i + 1
!         if (cstr(i)==c_null_char) exit
!         fstr = fstr // cstr(i)
!     end do
! end function c2f_str

subroutine initialize_c(ipointer) bind(C, name="initialize_c")
    !! create a [[radbelt_type]] from C
    integer(c_intptr_t),intent(out) :: ipointer
    type(radbelt_type),pointer :: p
    type(c_ptr) :: cp

    allocate(p)
    cp = c_loc(p)
    ipointer = transfer(cp, 0_c_intptr_t)

end subroutine initialize_c

subroutine destroy_c(ipointer) bind(C, name="destroy_c")
    !! destroy a [[radbelt_type]] from C
    integer(c_intptr_t),intent(in) :: ipointer
    type(radbelt_type),pointer :: p
    type(c_ptr) :: cp

    call int_pointer_to_f_pointer(ipointer,p)
    if (associated(p)) deallocate(p)

end subroutine destroy_c

subroutine int_pointer_to_f_pointer(ipointer, p)

    integer(c_intptr_t),intent(in) :: ipointer
    type(radbelt_type),pointer :: p

    type(c_ptr) :: cp

    cp = transfer(ipointer, c_null_ptr)
    if (c_associated(cp)) then
        call c_f_pointer(cp, p)
    else
        p => null()
    end if

end subroutine int_pointer_to_f_pointer

!*****************************************************************************************
!>
!  C interface for setting the data file paths

subroutine set_data_files_paths_c(ipointer, aep8_dir, igrf_dir) bind(C, name="set_data_files_paths_c")

    integer(c_intptr_t),intent(in) :: ipointer
    character(kind=c_char,len=1),dimension(STR_LEN),intent(in) :: aep8_dir
    character(kind=c_char,len=1),dimension(STR_LEN),intent(in) :: igrf_dir

    character(len=:),allocatable :: aep8_dir_, igrf_dir_
    type(radbelt_type),pointer :: p

    c2f : block
        ! convert to c strings to fortran
        integer :: i
        aep8_dir_ = ''
        igrf_dir_ = ''
        do i = 1, STR_LEN
            aep8_dir_ = aep8_dir_//aep8_dir(i)
            igrf_dir_ = igrf_dir_//igrf_dir(i)
        end do
        aep8_dir_ = trim(aep8_dir_)
        igrf_dir_ = trim(igrf_dir_)
    end block c2f

    call int_pointer_to_f_pointer(ipointer, p)

    ! !... doesn't work this way... ??
    ! aep8_dir_ = c2f_str(aep8_dir)
    ! igrf_dir_ = c2f_str(igrf_dir)
    ! write(*,*) trim(aep8_dir_)
    ! write(*,*) trim(igrf_dir_)

    call p%set_data_files_paths(aep8_dir_, igrf_dir_)

 end subroutine set_data_files_paths_c
!*****************************************************************************************

!*****************************************************************************************
!>
!  C interface to [[get_flux_g]].

subroutine get_flux_g_c(ipointer,lon,lat,height,year,e,imname,flux) bind(C, name="get_flux_g_c")

    integer(c_intptr_t),intent(in) :: ipointer
    real(c_double),intent(in) :: lon !! geodetic longitude in degrees (east)
    real(c_double),intent(in) :: lat !! geodetic latitude in degrees (north)
    real(c_double),intent(in) :: height !! altitude in km above sea level
    real(c_double),intent(in) :: year !! decimal year for which geomagnetic field is to
                                      !! be calculated (e.g.:1995.5 for day 185 of 1995)
    real(c_double),intent(in) :: e !! minimum energy
    integer(c_int),intent(in) :: imname !! which method to use:
                                 !!
                                 !! * 1 -- particle species: electrons, solar activity: min
                                 !! * 2 -- particle species: electrons, solar activity: max
                                 !! * 3 -- particle species: protons, solar activity: min
                                 !! * 4 -- particle species: protons, solar activity: max
    real(c_double),intent(out) :: flux !! The flux of particles above the given energy, in units of cm^-2 s^-1.

    type(radbelt_type),pointer :: p

    call int_pointer_to_f_pointer(ipointer, p)

    flux = p%get_flux(lon,lat,height,year,e,imname)

end subroutine get_flux_g_c

end module radbelt_c_module
