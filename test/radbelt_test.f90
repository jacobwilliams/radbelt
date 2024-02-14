program radbelt_test

    !! comparison to the python radbelt example

    use radbelt_module
    use radbelt_kinds_module

    implicit none

    real(wp) :: lon, lat, height, year, e, flux, error, relerror
    integer :: imname, i, ilat, ilon, ialt, iunit, istat
    real(wp) :: tstart, tend
    type(radbelt_type) :: radbelt
    integer :: n_cases

    lon = -45.0_wp
    lat = -30.0_wp
    height = 500.0_wp

    ! >>> from astropy.time import Time
    ! >>> time = Time('2021-03-01')
    ! >>> time.utc.decimalyear

    year = 2021.1616438356164_wp  ! decimal year
    Imname = 4 ! 'p', 'max'
    e = 20.0_wp

    ! simple example:
    flux = get_flux(lon,lat,height,year,e,imname)
    ! error = Flux - 2642.50268555_wp  ! difference from python wrapper version (radbelt)
    error = Flux - 2642.50370051985726336559603128948869_wp ! difference from real128 version
    relerror = abs(error/flux)
    write(*,*) 'Flux      = ', flux
    write(*,*) 'Error     = ', error
    write(*,*) 'Rel Error = ', relerror
    if (relerror>10*epsilon(1.0_wp)) error stop 'error'
    write(*,*) ''

    ! cartesian test:
    flux = get_flux([    0.66161289996712280_wp,&
                        -0.66161289996712269_wp,&
                        -0.53685097954130001_wp],year,e,imname)
    write(*,*) 'Flux      = ', flux
    write(*,*) 'Error     = ', error
    write(*,*) 'Rel Error = ', relerror
    if (relerror>10*epsilon(1.0_wp)) error stop 'error'
    write(*,*) ''

    ! speed tests:
    open(newunit=iunit, file = 'results.txt', iostat=istat)
    call cpu_time(tstart)
    n_cases = 0
    do ilat = -89, 90, 5
        do ilon = -180, 180, 45
            do ialt = 500, 1000, 100
                n_cases = n_cases + 1
                flux = get_flux(real(ilon,wp),real(ilat,wp),real(ialt,wp),year,e,imname)
                write(iunit,*) year, ',', ilat, ',', ilon, ',', ialt, ',', flux
            end do
        end do
    end do
    call cpu_time(tend)
    close(iunit)
    write(*,'(a25,f6.3,a,i7,a)') 'Function version runtime: ', (tend - tstart), ' sec. ', &
                                 int(n_cases/(tend - tstart)), ' (cases/sec)'

    call cpu_time(tstart)
    n_cases = 0
    do ilat = -90, 90, 5
        do ilon = -180, 180, 45
            do ialt = 500, 1000, 100
                n_cases = n_cases + 1
                flux = radbelt%get_flux(real(ilon,wp),real(ilat,wp),real(ialt,wp),year,e,imname)
            end do
        end do
    end do
    call cpu_time(tend)
    write(*,'(a25,f6.3,a,i7,a)') 'Class version runtime: ', (tend - tstart), ' sec. ', &
                                 int(n_cases/(tend - tstart)), ' (cases/sec)'

    write(*,*) 'n_cases = ', n_cases

end program radbelt_test
