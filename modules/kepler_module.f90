module kepler_module

    use numbers_module
    use iso_fortran_env, only: error_unit

    implicit none

    private

    !public :: kepler_shepperd
    !public :: kepler_goodyear_stienon_klumpp
    public :: kepler_classical

    contains
    !*******************************************************************************

    !*******************************************************************************
    !>
    !  Classical Kepler propagator for elliptical and hyperbolic orbits.
    !  Uses Lagrange formulations from Battin & Newton's method.
    !
    !### See also
    !  * Dario Izzo: pykep/src/core_functions/propagate_lagrangian.h

    subroutine kepler_classical(x0, dt, mu, xf)

        use newton_module, only: newton

        implicit none

        real(kind=8),dimension(6),intent(in)  :: x0  !! initial position,velocity vector
        real(kind=8),intent(in)               :: dt  !! propagation time
        real(kind=8),intent(in)               :: mu  !! central body gravitational parameter
        real(kind=8),dimension(6),intent(out) :: xf  !! final position,velocity vector

        real(kind=8),dimension(3) :: r0  !! initial position vector
        real(kind=8),dimension(3) :: v0  !! initial velocity vector
        real(kind=8) :: de     !! eccentric anomaly difference
        real(kind=8) :: dh     !! hyperbolic anomaly difference
        integer  :: iflag  !! newton status flag
        real(kind=8) :: r,rmag,vmag,energy,a,sqrta,f,&
                    g,ft,gt,sigma0,dm,dn,xs,fx,b,p,x,z,term

        real(kind=8),parameter :: parabolic_tol = 1.0e-12 !! zero tol for parabolic orbits (energy)
        real(kind=8),parameter :: ftol = 1.0e-12          !! function tol for root finding
        real(kind=8),parameter :: xtol = 1.0e-12          !! indep variable tol for root finding
        integer,parameter  :: max_iter = 1000             !! maximum number of iterations in newton

        ! check trivial case:
        if (dt==zero) then
            xf = x0
            return
        end if

        r0     = x0(1:3)
        v0     = x0(4:6)
        rmag   = norm2(r0)
        vmag   = norm2(v0)
        energy = (vmag * vmag / two - mu / rmag)
        sigma0 = dot_product(r0,v0) / sqrt(mu)

        ! if not parabolic, then compute semimajor axis
        if (abs(energy) > parabolic_tol) then
            a = -mu / two / energy
        end if

        if (energy < -parabolic_tol) then ! elliptical case

            sqrta = sqrt(a)
            dm = sqrt(mu / a**3) * dt
            de = dm

            call newton(de,kepde_,d_kepde_,ftol,xtol,max_iter,xs,fx,iflag)
            if (iflag<0) then
                write(error_unit,'(A)') 'Error in kepler_classical [elliptical]: newton did not converge'
                write(*,*) xs,fx,iflag
            end if

            de = xs
            r = a + (rmag - a) * cos(de) + sigma0 * sqrta * sin(de) ! eqn 4.42

            ! lagrange coefficients (Battin eqn 4.41)
            f  = one - a / rmag * (one - cos(de))
            g  = a * sigma0 / sqrt(mu) * (one - cos(de)) + rmag * sqrt(a / mu) * sin(de)
            ft = -sqrt(mu * a) / (r * rmag) * sin(de)
            gt = one - a / r * (one - cos(de))

        else if (energy > parabolic_tol) then ! hyperbolic case

            sqrta = sqrt(-a)
            dn = sqrt(-mu / a**3) * dt
            dh = sign(one,dt)      ! todo: need a better initial guess

            call newton(dh,kepdh_,d_kepdh_,ftol,xtol,max_iter,xs,fx,iflag)
            if (iflag<0) then
                write(error_unit,'(A)') 'Error in kepler_classical [hyperbola]: newton did not converge'
                write(*,*) xs,fx,iflag
            end if

            dh = xs
            r = a + (rmag - a) * cosh(dh) + sigma0 * sqrta * sinh(dh)

            ! lagrange coefficients (Battin eqn 4.62)
            f  = one - a / rmag * (one - cosh(dh))
            g  = a * sigma0 / sqrt(mu) * (one - cosh(dh)) + rmag * sqrt(-a / mu) * sinh(dh)
            ft = -sqrt(-mu * a) / (r * rmag) * sinh(dh)
            gt = one - a / r * (one - cosh(dh))

        else  ! parbolic case

            ! See Battin Section 4.2

            p    = two * rmag - sigma0**2
            B    = one / p**1.5 * ( sigma0 * (rmag + p) + three * sqrt(mu) * dt )
            term = B + sqrt(one+B*B)
            z    = (term)**(one/three) - (term)**(-one/three) ! Battin eqn 4.12 where z = tan(f/2)
            x    = sqrt(p) * z - sigma0
            r    = rmag + sigma0 * x + one/two * x**2

            f  = one - x**2 / two / rmag
            g  = x / two / sqrt(mu) * ( two * rmag + sigma0 * x )
            ft = - sqrt(mu) * x / rmag / r
            gt = one - x**2 / two / r

        end if

        ! results:
        xf(1:3) = f  * r0 + g  * v0
        xf(4:6) = ft * r0 + gt * v0

        contains

            ! function wrappers for newtons method:

            subroutine kepde_(x,f)
            implicit none
            real(kind=8),intent(in)  :: x  !! de
            real(kind=8),intent(out) :: f
            f = kepde(x, dm, sigma0, sqrta, a, rmag)
            end subroutine kepde_

            subroutine d_kepde_(x,f)
            implicit none
            real(kind=8),intent(in)  :: x  !! de
            real(kind=8),intent(out) :: f
            f = d_kepde(x, sigma0, sqrta, a, rmag)
            end subroutine d_kepde_

            subroutine kepdh_(x,f)
            implicit none
            real(kind=8),intent(in)  :: x  !! dh
            real(kind=8),intent(out) :: f
            f = kepdh(x, dn, sigma0, sqrta, a, rmag)
            end subroutine kepdh_

            subroutine d_kepdh_(x,f)
            implicit none
            real(kind=8),intent(in)  :: x  !! dh
            real(kind=8),intent(out) :: f
            f = d_kepdh(x, sigma0, sqrta, a, rmag)
            end subroutine d_kepdh_

    end subroutine kepler_classical

    !*******************************************************************************

    !*******************************************************************************
    !>
    !  Elliptic Kepler's equation

    pure function kepe(e, m, ecc)

        implicit none

        real(kind=8),intent(in) :: e   !! eccentric anomaly
        real(kind=8),intent(in) :: m   !! mean anomaly
        real(kind=8),intent(in) :: ecc !! eccentricity
        real(kind=8) :: kepe

        kepe = e - ecc * sin(e) - m

    end function kepe
    !*******************************************************************************

    !*******************************************************************************
    !>
    !  Derivative of [[kepe]] w.r.t. `e`

    pure function d_kepe(e, ecc)

        implicit none

        real(kind=8),intent(in) :: e   !! eccentric anomaly
        real(kind=8),intent(in) :: ecc !! eccentricity
        real(kind=8) :: d_kepe

        d_kepe = one - ecc * cos(e)

    end function d_kepe
    !*******************************************************************************

    !*******************************************************************************
    !>
    !  Elliptic Kepler's equation written in terms of the
    !  eccentric anomaly difference.  See Battin, eqn 4.43.

    pure function kepde(de, dm, sigma0, sqrta, a, r)

        implicit none

        real(kind=8),intent(in) :: de      !! eccentric anomaly difference
        real(kind=8),intent(in) :: dm      !! mean anomaly difference
        real(kind=8),intent(in) :: sigma0
        real(kind=8),intent(in) :: sqrta
        real(kind=8),intent(in) :: a
        real(kind=8),intent(in) :: r
        real(kind=8) :: kepde

        kepde = -dm + de + sigma0 / sqrta * (one - cos(de)) - &
                (one - r / a) * sin(de)

    end function kepde
    !*******************************************************************************

    !*******************************************************************************
    !>
    !  Derivative of [[kepde]] w.r.t `de`.

    pure function d_kepde(de, sigma0, sqrta, a, r)

        implicit none

        real(kind=8),intent(in) :: de      !! eccentric anomaly difference
        real(kind=8),intent(in) :: sigma0
        real(kind=8),intent(in) :: sqrta
        real(kind=8),intent(in) :: a
        real(kind=8),intent(in) :: r
        real(kind=8) :: d_kepde

        d_kepde = one + sigma0 / sqrta * sin(de) - (one - r / a) * cos(de)

    end function d_kepde
    !*******************************************************************************

    !*******************************************************************************
    !>
    !  Battin, eqn. 4.64.

    pure function kepdh(dh, dn, sigma0, sqrta, a, r)

        implicit none

        real(kind=8),intent(in) :: dh      !! hyperbolic anomaly difference
        real(kind=8),intent(in) :: dn
        real(kind=8),intent(in) :: sigma0
        real(kind=8),intent(in) :: sqrta
        real(kind=8),intent(in) :: a
        real(kind=8),intent(in) :: r
        real(kind=8) :: kepdh

        kepdh = -dn - dh + sigma0 / sqrta * (cosh(dh) - one) + &
                (one - r / a) * sinh(dh)

    end function kepdh
    !*******************************************************************************

    !*******************************************************************************
    !>
    !  Derivative of [[kepdh]] w.r.t `dh`.

    pure function d_kepdh(dh, sigma0, sqrta, a, r)

        implicit none

        real(kind=8),intent(in) :: dh      !! hyperbolic anomaly difference
        real(kind=8),intent(in) :: sigma0
        real(kind=8),intent(in) :: sqrta
        real(kind=8),intent(in) :: a
        real(kind=8),intent(in) :: r
        real(kind=8) :: d_kepdh

        d_kepdh = -one + sigma0 / sqrta * sinh(dh) + (one - r / a) * cosh(dh)

    end function d_kepdh
    !*******************************************************************************

    !*******************************************************************************
    !>
    !  Barker time of flight equation

    pure function barker(r1, r2, mu)

        implicit none

        real(kind=8),dimension(3),intent(in) :: r1
        real(kind=8),dimension(3),intent(in) :: r2
        real(kind=8),intent(in) :: mu
        real(kind=8) :: barker

        real(kind=8) :: x, r1mag, r2mag, r21mag, sigma, r1pr2mag
        real(kind=8),dimension(3) :: r21

        r1mag    = norm2(r1)
        r2mag    = norm2(r2)
        r21      = r2 - r1
        r21mag   = norm2(r21)
        x        = r1(1) * r2(2) - r1(2) * r2(1)
        sigma    = sign(one,x)
        r1pr2mag = r1mag + r2mag

        barker = (r1pr2mag + r21mag)**1.5 - &
                sigma * (r1pr2mag - r21mag)**1.5 / (six * sqrt(mu))

    end function barker
    !*******************************************************************************

    !*******************************************************************************
    !>
    !
    pure function kepds(ds, dt, r0, vr0, alpha, mu)

        implicit none

        real(kind=8),intent(in) :: ds  !! universal anomaly difference
        real(kind=8),intent(in) :: dt
        real(kind=8),intent(in) :: r0
        real(kind=8),intent(in) :: vr0
        real(kind=8),intent(in) :: alpha
        real(kind=8),intent(in) :: mu
        real(kind=8) :: kepds

        real(kind=8) :: c,s,ads2,ds2

        ds2  = ds * ds
        ads2 = alpha * ds2
        s    = stumpff_s(ads2)
        c    = stumpff_c(ads2)

        kepds = -sqrt(mu) * dt + r0 * vr0 * ds2 * c / sqrt(mu) + &
                (one - alpha * r0) * ds2 * ds * s + r0 * ds

    end function kepds
    !*******************************************************************************

    !*******************************************************************************
    !>
    !
    pure function d_kepds(ds, r0, vr0, alpha, mu)

        implicit none

        real(kind=8),intent(in) :: ds
        real(kind=8),intent(in) :: r0
        real(kind=8),intent(in) :: vr0
        real(kind=8),intent(in) :: alpha
        real(kind=8),intent(in) :: mu
        real(kind=8) :: d_kepds

        real(kind=8) :: c,s,ads2,ds2

        ds2  = ds * ds
        ads2 = alpha * ds2
        s    = stumpff_s(ads2)
        c    = stumpff_c(ads2)

        d_kepds = r0 * vr0 / sqrt(mu) * ds * (one - ads2 * s) + &
                (one - alpha * r0) * ds2 * c + r0

    end function d_kepds
    !*******************************************************************************

    !*******************************************************************************
    !>
    !  Stumpff function S(z)

    pure function stumpff_s(z) result(s)

        implicit none

        real(kind=8),intent(in) :: z
        real(kind=8) :: s

        if (z > zero) then
            s = (sqrt(z) - sin(sqrt(z))) / sqrt(z)**3
        else if (z < zero) then
            s = (sinh(sqrt(-z)) - sqrt(-z)) / (-z)**(three/two)
        else
            s = one/six
        end if

    end function stumpff_s
    !*******************************************************************************

    !*******************************************************************************
    !>
    !  Stumpff function C(z)

    pure function stumpff_c(z) result(c)

        implicit none

        real(kind=8),intent(in) :: z
        real(kind=8) :: c

        if (z > zero) then
            c = (one - cos(sqrt(z))) / z
        else if (z < zero) then
            c = (cosh(sqrt(-z)) - one) / (-z)
        else
            c = 0.5
        end if

    end function stumpff_c


end module kepler_module