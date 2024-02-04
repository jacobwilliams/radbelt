program radbelt_test

    !! comparison to the python radbelt example
    
    use core 

    implicit none 

    real :: lon, lat, height, year, xl,bbx 
    REAL, DIMENSION(1) :: E
    INTEGER :: Imname
    REAL, DIMENSION(1) :: Flux, error
   
    lon = -45.0
    lat = -30.0 
    height = 500.0 

    ! >>> from astropy.time import Time
    ! >>> time = Time('2021-03-01')
    ! >>> time.utc.decimalyear

    year = 2021.1616438356164  ! decimal year

    Imname = 4 ! 'p', 'max'
    e = 20.0

    call igrf(Lon,Lat,Height,Year,Xl,Bbx)
    call aep8(E,Xl,Bbx,Imname,Flux)

    error = Flux - 2642.50268555  ! difference from python wrapper version (radbelt)

    write(*,*) 'Flux  = ', Flux
    write(*,*) 'Error = ', error

    if (abs(error(1))>1.0e-9) error stop 'error'

end program radbelt_test
