program astrodynamics_lib

    ! Include modules
    use numbers_module
    use kepler_module
    use newton_module
    use lambert_module
    use vector_module

    implicit none
    
    ! Auxiliary Variables
    real(kind=8)            :: x0(6), xf(6), mu, dt, r1(3), r2(3), tof
    real(kind=8),dimension(:,:),allocatable :: v1, v2
    logical                 :: long_way, status_ok
    integer                 :: multi_revs


    WRITE(*,*) 'Testing Kepler: '

    mu = 3.986004418e5
    dt = 56.483821570873260

    x0 = [40799.082147252739, 8352.4309740588305, -74.782959949448980,&
                                    -0.60876318325705059, 2.9746097973556176, -1.1277840825478561E-003]

    CALL kepler_classical(x0,dt,mu,xf)

    WRITE(*,*) xf(1),xf(2),xf(3)
    WRITE(*,*) xf(4),xf(5),xf(6)

    r1 = [40763.630678738205,8557.8020389363464,-52.756128870937815]
    r2 = [40763.436536803158,8560.3934649368257,-53.651459920556590]
    tof = 86165.171018242836
    long_way = .false.
    multi_revs = 2

    WRITE(*,*) 'Testing Lamberts'
    call solve_lambert_izzo(r1,r2,tof,mu,long_way,multi_revs,v1,v2,status_ok)
    
    WRITE(*,*) 'v1', v1(:,1)
    WRITE(*,*) 'v2', v2(:,1)

end program astrodynamics_lib