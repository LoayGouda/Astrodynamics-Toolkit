module lambert_module

    use numbers_module
    use vector_module,    only: cross, unit, ucross

    implicit none

    private

    !constants:
    real(kind=8),parameter :: log2       = log(two)
    real(kind=8),parameter :: two_third  = two/three
    real(kind=8),parameter :: four_third = four/three
    real(kind=8),parameter :: five_half  = five/two
    real(kind=8),parameter :: three_half = three/two

    abstract interface
        function func(t) result(f)  !! interface to the [[zeroin]] input function
        implicit none
        real(kind=8),intent(in)  :: t  !! Independant variable for the function.
        real(kind=8)             :: f  !! The function evaluated at `t`.
        end function func
    end interface

    !public routines:
    public :: solve_lambert_izzo
    !public :: solve_lambert_gooding
    !public :: solve_lambert_arorarussell

    !public :: lambert_test

    contains
    !*****************************************************************************************

    !*****************************************************************************************
    !>
    !  Solve Lambert's problem using Izzo's method.
    !
    !# References
    !
    !  1. D. Izzo, [Revisiting Lambert's Problem](http://arxiv-web3.library.cornell.edu/abs/1403.2705)
    !     [v2] Tue, 24 Jun 2014 13:08:37 GMT (606kb,D)
    !  2. [PyKEP](https://github.com/esa/pykep)
    !  3. R. A. Battin, "An Introduction to the Mathematics and Methods of
    !     Astrodynamics (Revised Edition)", AIAA Education Series, 1999.

    subroutine solve_lambert_izzo(r1,r2,tof,mu,long_way,multi_revs,v1,v2,status_ok)

        implicit none

        real(kind=8),dimension(3),intent(in)                :: r1            !! first cartesian position [km]
        real(kind=8),dimension(3),intent(in)                :: r2            !! second cartesian position [km]
        real(kind=8),intent(in)                             :: tof           !! time of flight [sec]
        real(kind=8),intent(in)                             :: mu            !! gravity parameter [km^3/s^2]
        logical,intent(in)                              :: long_way      !! when true, do "long way" (>pi) transfers
        integer,intent(in)                              :: multi_revs    !! maximum number of multi-rev solutions to compute
        real(kind=8),dimension(:,:),allocatable,intent(out) :: v1            !! vector containing 3d arrays with the cartesian components of the velocities at r1
        real(kind=8),dimension(:,:),allocatable,intent(out) :: v2            !! vector containing 3d arrays with the cartesian components of the velocities at r2
        logical,intent(out)                             :: status_ok     !! true if everything is OK

        !local variables:
        real(kind=8),dimension(:),allocatable :: x
        real(kind=8),dimension(3) :: r1_hat,r2_hat,h_hat,it1,it2,c
        real(kind=8) :: s,cmag,lambda2,lambda3,lambda5,t,t00,t0,t1,r1mag,r2mag,&
                    d3t,d2t,dt,err,t_min,x_old,x_new,term,lambda,&
                    gamma,rho,sigma,vr1,vt1,vr2,vt2,y,vt,ly
        integer :: n_solutions,it,m_nmax,i,iter

        !tolerances are from [2]
        integer,parameter   :: max_halley_iters = 12         !! for halley iterations
        real(kind=8),parameter  :: halley_tol       = 1.0e-13   !! for halley iterations
        real(kind=8),parameter  :: htol_singlerev   = 1.0e-5    !! for householder iterations
        real(kind=8),parameter  :: htol_multirev    = 1.0e-8    !! for householder iterations

        !======= Begin Algorithm 1 in [1] =======

        r1mag = norm2(r1)
        r2mag = norm2(r2)

        !check for valid inputs:
        if (tof<=zero .or. mu<=zero .or. r1mag==zero .or. r2mag==zero) then
            write(*,*) 'Error in solve_lambert_izzo: invalid input'
            write(*,*) tof
            status_ok = .false.
            return
        end if

        status_ok = .true.

        c       = r2 - r1
        cmag    = norm2(c)
        s       = (cmag + r1mag + r2mag) / two
        t       = sqrt(two*mu/(s*s*s)) * tof
        r1_hat  = unit(r1)
        r2_hat  = unit(r2)
        h_hat   = ucross(r1_hat,r2_hat)
        lambda2 = one - cmag/s
        lambda  = sqrt(lambda2)

        if ( all(h_hat == zero) ) then
            write(*,*) 'Warning: pi transfer in solve_lambert_izzo'
            !arbitrarily choose the transfer plane:
            h_hat = [zero,zero,one]
        end if

        it1 = ucross(h_hat,r1_hat)
        it2 = ucross(h_hat,r2_hat)
        if (long_way) then
            lambda = -lambda
            it1 = -it1
            it2 = -it2
        end if

        lambda3 = lambda*lambda2
        lambda5 = lambda2*lambda3
        t1 = two_third * (one - lambda3)

        !======= Begin Algorithm 2 in [1] =======
        ![xlist, ylist] = findxy(lambda, tof)

        ! maximum number of revolutions for which a solution exists:
        m_nmax = floor(t/pi)

        t00 = acos(lambda) + lambda*sqrt(one-lambda2)
        t0 = t00 + m_nmax*pi

        if (t < t0 .and. m_nmax > 0) then    ! Compute xm and tm using Halley

            dt = zero
            d2t = zero
            d3t = zero
            it = 0
            err = one
            t_min = t0
            x_old = zero
            x_new = zero

            do

                call dtdx(dt,d2t,d3t,x_old,t_min)
                if (dt /= zero) x_new = x_old - dt*d2t/(d2t * d2t - dt * d3t / two)
                err = abs(x_old-x_new)
                if ( (err<halley_tol) .or. (it>max_halley_iters) ) exit
                call compute_tof(x_new,m_nmax,t_min)
                x_old = x_new
                it = it + 1

            end do

            if (t_min > t) m_nmax = m_nmax - 1

        end if
        !======= End Algorithm 2 =======

        !mmax is the maximum number of revolutions.
        !Truncate to user-input multi_revs value if it is larger.
        m_nmax = min(multi_revs,m_nmax)

        !the number of solutions to the problem:
        n_solutions = m_nmax*2 + 1

        !allocate output arrays:
        allocate(v1(3,n_solutions))
        allocate(v2(3,n_solutions))
        allocate(x(n_solutions))

        ! Find the x value for each solution:

        ! initial guess for 0 rev solution:

        if (t>=t00) then

            x(1) = -(t-t00)/(t-t00+four)                        !from [2]

        elseif (t<=t1) then

            x(1) = five_half * (t1*(t1-t))/(t*(one-lambda5)) + one

        else

            x(1) = (t/t00) ** ( log2 / log(t1/t00) ) - one      !from [2]

        end if

        ! 0 rev solution:

        iter = householder(t, x(1), 0, htol_singlerev)

        ! multi-rev solutions:

        do i = 1,m_nmax

            !Eqn 31:

            ! left solution:
            term     = ((i*pi+pi)/(eight*t)) ** two_third
            x(2*i)   = (term-one)/(term+one)
            iter     = householder(t, x(2*i), i, htol_multirev)

            ! right solution:
            term     = ((eight*t)/(i*pi)) ** two_third
            x(2*i+1) = (term-one)/(term+one)
            iter     = householder(t, x(2*i+1), i, htol_multirev)

        end do

        ! construct terminal velocity vectors using each x:

        gamma = sqrt(mu*s/two)
        rho   = (r1mag-r2mag) / cmag
        sigma = sqrt(one-rho*rho)

        do i=1,n_solutions

            y   = sqrt(one - lambda2 + lambda2*x(i)*x(i))
            ly  = lambda*y
            vr1 = gamma*((ly-x(i))-rho*(ly+x(i)))/r1mag
            vr2 = -gamma*((ly-x(i))+rho*(ly+x(i)))/r2mag
            vt  = gamma*sigma*(y+lambda*x(i))
            vt1 = vt/r1mag
            vt2 = vt/r2mag

            v1(:,i) = vr1*r1_hat + vt1*it1   !terminal velocity vectors
            v2(:,i) = vr2*r2_hat + vt2*it2   !

        end do

        deallocate(x)

        contains
        !*****************************************************************************************

            !*************************************************************************************
                function householder(t,x,n,eps) result(it)

                !! Householder root solver for x.

                implicit none

                integer                 :: it
                real(kind=8),intent(in)     :: t
                real(kind=8),intent(inout)  :: x    !! input is initial guess
                integer,intent(in)      :: n
                real(kind=8),intent(in)     :: eps

                real(kind=8) :: xnew,tof,delta,dt,d2t,d3t,dt2,term

                integer,parameter :: max_iters = 15

                do it = 1, max_iters

                    call compute_tof(x,n,tof)
                    call dtdx(dt,d2t,d3t,x,tof)

                    delta = tof-t
                    dt2   = dt*dt
                    term  = delta*(dt2-delta*d2t/two)/&
                            (dt*(dt2-delta*d2t)+d3t*delta*delta/six)
                    xnew  = x - term    ! Ref. [1], p. 12.
                    x     = xnew

                    if (abs(term)<=eps) exit

                end do

                end function householder
            !*************************************************************************************

            !*************************************************************************************
                subroutine dtdx(dt,d2t,d3t,x,t)

                !! Compute 1st-3rd derivatives for the Householder iterations.

                implicit none

                real(kind=8),intent(out)  :: dt
                real(kind=8),intent(out)  :: d2t
                real(kind=8),intent(out)  :: d3t
                real(kind=8),intent(in)   :: x
                real(kind=8),intent(in)   :: t

                real(kind=8) :: umx2,y,y2,y3,y5,umx2_inv

                umx2     = one-x*x
                umx2_inv = one/umx2
                y        = sqrt(one-lambda2*umx2)    !Ref [1], p. 6
                y2       = y*y
                y3       = y2*y
                y5       = y3*y2

                !Ref [1], Eqn. 22:

                dt       = umx2_inv * (three*t*x - two + two*lambda3*x/y)
                d2t      = umx2_inv * (three*t + five*x*dt + two*(one-lambda2)*lambda3/y3)
                d3t      = umx2_inv * (seven*x*d2t + eight*dt - six*(one-lambda2)*lambda5*x/y5)

                end subroutine dtdx
            !*************************************************************************************

            !*************************************************************************************
                subroutine compute_tof(x,n,tof)

                !!  Compute time of flight from x

                implicit none

                real(kind=8),intent(in)   :: x
                integer,intent(in)    :: n
                real(kind=8),intent(out)  :: tof

                real(kind=8),parameter :: battin   = 0.01
                real(kind=8),parameter :: lagrange = 0.2

                real(kind=8) :: dist,k,e,rho,z,eta,s1,q,y,g,d,l,f,a,alpha,beta

                dist = abs(x-one)

                if (dist < lagrange .and. dist > battin) then    !use lagrange tof expression

                    !See Ref. [1], Eqn. 9

                    a = one / (one-x*x)

                    if (a>zero) then   !ellipse

                        alpha = two * acos(x)
                        beta = two * asin(sqrt(lambda2/a))
                        if (lambda<zero) beta = -beta
                        tof = ((a*sqrt(a)*((alpha-sin(alpha))-(beta-sin(beta))+two*pi*n))/two)

                    else   !hyperbola

                        alpha = two * acosh(x)
                        beta = two * asinh(sqrt(-lambda2/a))
                        if (lambda<zero) beta = -beta
                        tof = (-a*sqrt(-a)*((beta-sinh(beta))-(alpha-sinh(alpha)))/two)

                    end if

                else

                    k    = lambda2
                    e    = x*x - one
                    rho  = abs(e)
                    z    = sqrt(one+k*e)

                    if (dist < battin) then  ! use battin series tof expression

                        !Equation 20 in [1]:
                        !
                        !See also: Ref. [3], Eqn. 7.30, p. 304.

                        eta = z - lambda*x
                        s1  = (one - lambda - x*eta)/two
                        q   = four_third*hypergeo(s1)
                        tof = (eta*eta*eta*q + four*lambda*eta)/two + n*pi / (rho**three_half)

                    else  ! use lancaster tof expresion

                        y = sqrt(rho)
                        g = x*z - lambda*e
                        d = zero
                        if (e < zero) then
                            l = acos(g)
                            d = n*pi+l
                        else
                            f = y*(z-lambda*x)
                            d = log(f+g)
                        end if
                        tof = (x-lambda*z-d/y)/e

                    end if

                end if

                end subroutine compute_tof
            !*************************************************************************************

            !*************************************************************************************
                pure function hypergeo(x) result(f)

                !!  Evaluate the Gaussian (or ordinary) hypergeometric function: F(3,1,5/2,x)
                !!  See Ref. [3], p. 34.

                implicit none

                real(kind=8)             :: f
                real(kind=8),intent(in)  :: x

                real(kind=8) :: term
                integer  :: i

                real(kind=8),parameter :: tol = 1.0e-11
                integer,parameter  :: max_iters = 10000

                !initialize:
                f    = one
                term = one

                !compute the series until the last term is within convergence tolerance:
                do i = 0, max_iters

                    term = term*(three+i)*(one+i) / (five_half+i)*x / (i+one)
                    f = f + term
                    if (abs(term)<=tol) exit

                end do

                end function hypergeo
            !*************************************************************************************

    end subroutine solve_lambert_izzo

end module lambert_module