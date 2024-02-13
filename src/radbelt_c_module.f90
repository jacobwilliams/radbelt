!*****************************************************************************************
!>
!  Experimental C interface to the radbelt module.
!
!@note Currently, this makes use of a global variable to store the state,
!      so should not be considered threadsafe.

module radbelt_c_module

use iso_c_binding, only: c_double, c_int, c_char, c_null_char
use radbelt_module, only: radbelt_type

implicit none

type(radbelt_type) :: radbelt !! global variable

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

!*****************************************************************************************
!>
!  C interface for setting the data file paths

subroutine set_data_files_paths_c(aep8_dir, igrf_dir) bind(C, name="set_data_files_paths_c")

    character(kind=c_char,len=1),dimension(STR_LEN),intent(in) :: aep8_dir
    character(kind=c_char,len=1),dimension(STR_LEN),intent(in) :: igrf_dir

    character(len=:),allocatable :: aep8_dir_, igrf_dir_

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

    ! !... doesn't work this way... ??
    ! aep8_dir_ = c2f_str(aep8_dir)
    ! igrf_dir_ = c2f_str(igrf_dir)
    ! write(*,*) trim(aep8_dir_)
    ! write(*,*) trim(igrf_dir_)

    call radbelt%set_data_files_paths(aep8_dir_, igrf_dir_)

 end subroutine set_data_files_paths_c
!*****************************************************************************************

!*****************************************************************************************
!>
!  C interface to [[get_flux_g]].

subroutine get_flux_g_c(lon,lat,height,year,e,imname,flux) bind(C, name="get_flux_g_c")

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

    flux = radbelt%get_flux(lon,lat,height,year,e,imname)

end subroutine get_flux_g_c

end module radbelt_c_module
