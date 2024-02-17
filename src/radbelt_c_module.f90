!*****************************************************************************************
!>
!  Experimental C interface to the radbelt module.

    module radbelt_c_module

    use iso_c_binding, only: c_double, c_int, c_char, c_null_char, &
                             c_intptr_t, c_ptr, c_loc, c_f_pointer, &
                             c_null_ptr, c_associated
    use radbelt_module, only: radbelt_type

    implicit none

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Convert C string to Fortran

function c2f_str(cstr) result(fstr)

    character(kind=c_char,len=1),dimension(:),intent(in) :: cstr !! string from C
    character(len=:),allocatable :: fstr !! fortran string

    integer :: i !! counter

    fstr = ''
    do i = 1, size(cstr)
        fstr = fstr//cstr(i)
    end do
    fstr = trim(fstr)

end function c2f_str

!*****************************************************************************************
!>
!  Convert an integer pointer to a [[radbelt_type]] pointer.

subroutine int_pointer_to_f_pointer(ipointer, p)

    integer(c_intptr_t),intent(in) :: ipointer !! integer pointer from C
    type(radbelt_type),pointer :: p !! fortran pointer

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
!  create a [[radbelt_type]] from C

subroutine initialize_c(ipointer) bind(C, name="initialize_c")

    integer(c_intptr_t),intent(out) :: ipointer
    type(radbelt_type),pointer :: p
    type(c_ptr) :: cp

    allocate(p)
    cp = c_loc(p)
    ipointer = transfer(cp, 0_c_intptr_t)

end subroutine initialize_c

!*****************************************************************************************
!>
!  destroy a [[radbelt_type]] from C

subroutine destroy_c(ipointer) bind(C, name="destroy_c")

    integer(c_intptr_t),intent(in) :: ipointer
    type(radbelt_type),pointer :: p

    call int_pointer_to_f_pointer(ipointer,p)
    if (associated(p)) deallocate(p)

end subroutine destroy_c

!*****************************************************************************************
!>
!  C interface for setting the `trm` data file path

subroutine set_trm_file_path_c(ipointer, aep8_dir, n) bind(C, name="set_trm_file_path_c")

    integer(c_intptr_t),intent(in) :: ipointer
    integer(c_int),intent(in) :: n !! size of `aep8_dir`
    character(kind=c_char,len=1),dimension(n),intent(in) :: aep8_dir

    character(len=:),allocatable :: aep8_dir_
    type(radbelt_type),pointer :: p

    call int_pointer_to_f_pointer(ipointer, p)

    if (associated(p)) then
        aep8_dir_ = c2f_str(aep8_dir)
        call p%set_trm_file_path(aep8_dir_)
    else
        error stop 'error in set_trm_file_path_c: class is not associated'
    end if

 end subroutine set_trm_file_path_c
!*****************************************************************************************

!*****************************************************************************************
!>
!  C interface for setting the `igrf` data file path

 subroutine set_igrf_file_path_c(ipointer, igrf_dir, n) bind(C, name="set_igrf_file_path")

    integer(c_intptr_t),intent(in) :: ipointer
    integer(c_int),intent(in) :: n !! size of `igrf_dir`
    character(kind=c_char,len=1),dimension(n),intent(in) :: igrf_dir

    character(len=:),allocatable :: igrf_dir_
    type(radbelt_type),pointer :: p

    call int_pointer_to_f_pointer(ipointer, p)

    if (associated(p)) then
        igrf_dir_ = c2f_str(igrf_dir)
        call p%set_igrf_file_path(igrf_dir_)
    else
        error stop 'error in set_igrf_file_path: class is not associated'
    end if

 end subroutine set_igrf_file_path_c
!*****************************************************************************************

!*****************************************************************************************
!>
!  C interface for setting the data file paths

 subroutine set_data_files_paths_c(ipointer, aep8_dir, igrf_dir, n, m) bind(C, name="set_data_files_paths_c")

    integer(c_intptr_t),intent(in) :: ipointer
    integer(c_int),intent(in) :: n !! size of `aep8_dir`
    character(kind=c_char,len=1),dimension(n),intent(in) :: aep8_dir
    integer(c_int),intent(in) :: m !! size of `igrf_dir`
    character(kind=c_char,len=1),dimension(m),intent(in) :: igrf_dir

    character(len=:),allocatable :: aep8_dir_, igrf_dir_
    type(radbelt_type),pointer :: p

    call int_pointer_to_f_pointer(ipointer, p)

    if (associated(p)) then

        aep8_dir_ = c2f_str(aep8_dir)
        igrf_dir_ = c2f_str(igrf_dir)

        call p%set_data_files_paths(aep8_dir_, igrf_dir_)

    else
        error stop 'error in set_data_files_paths_c: class is not associated'
    end if

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

    if (associated(p)) then

        flux = p%get_flux(lon,lat,height,year,e,imname)

    else
        error stop 'error in get_flux_g_c: class is not associated'
    end if

end subroutine get_flux_g_c

!*****************************************************************************************
    end module radbelt_c_module
!*****************************************************************************************