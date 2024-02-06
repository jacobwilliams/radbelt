program radbelt_test

    !! comparison to the python radbelt example

    use core
    use radbelt_kinds_module

    implicit none

    real(wp) :: lon, lat, height, year, e, flux, error, relerror
    integer :: imname, i

    lon = -45.0_wp
    lat = -30.0_wp
    height = 500.0_wp

    ! >>> from astropy.time import Time
    ! >>> time = Time('2021-03-01')
    ! >>> time.utc.decimalyear

    year = 2021.1616438356164_wp  ! decimal year
    Imname = 4 ! 'p', 'max'
    e = 20.0_wp

    do i = 1, 3
    flux = get_flux(lon,lat,height,year,e,imname)

    ! error = Flux - 2642.50268555_wp  ! difference from python wrapper version (radbelt)
    error = Flux - 2642.50370051985726336559603128948869_wp ! difference from real128 version
    relerror = abs(error/flux)

    write(*,*) 'Flux      = ', flux
    write(*,*) 'Error     = ', error
    write(*,*) 'Rel Error = ', relerror

    if (relerror>10*epsilon(1.0_wp)) error stop 'error'

    end do

end program radbelt_test
