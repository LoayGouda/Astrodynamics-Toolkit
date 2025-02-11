module numbers_module


    private

    real(kind=8),parameter,public :: zero       = 0.0
    real(kind=8),parameter,public :: one        = 1.0
    real(kind=8),parameter,public :: two        = 2.0
    real(kind=8),parameter,public :: three      = 3.0
    real(kind=8),parameter,public :: four       = 4.0
    real(kind=8),parameter,public :: five       = 5.0
    real(kind=8),parameter,public :: six        = 6.0
    real(kind=8),parameter,public :: seven      = 7.0
    real(kind=8),parameter,public :: eight      = 8.0
    real(kind=8),parameter,public :: nine       = 9.0
    real(kind=8),parameter,public :: ten        = 10.0

    real(kind=8),parameter,public :: pi         = acos(-one)
    real(kind=8),parameter,public :: twopi      = two*pi
    real(kind=8),parameter,public :: fourpi     = four*pi
    real(kind=8),parameter,public :: halfpi     = 0.5*pi

    real(kind=8),parameter,public :: universal_grav_constant = 6.67408e-20 !! CODATA-recommended universal gravitational
                                                                          !! constant \( km^3/kg-s^2  \)

    !> 3x3 identity matrix:
    real(kind=8),dimension(3,3),parameter,public :: identity_3x3 = reshape(&
                                                    [[one,zero,zero],&
                                                     [zero,one,zero],&
                                                     [zero,zero,one]],[3,3])

    !> 6x6 identity matrix:
    real(kind=8),dimension(6,6),parameter,public :: identity_6x6= reshape(&
                                                    [[one,zero,zero,zero,zero,zero],&
                                                     [zero,one,zero,zero,zero,zero],&
                                                     [zero,zero,one,zero,zero,zero],&
                                                     [zero,zero,zero,one,zero,zero],&
                                                     [zero,zero,zero,zero,one,zero],&
                                                     [zero,zero,zero,zero,zero,one] ],[6,6])

    end module numbers_module