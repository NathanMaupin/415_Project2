subroutine eulerSolver(x, y)
    use common_parameters
    implicit none

    ! From main program
    double precision, dimension(iMax, jMax), intent(in):: x, y

    ! Upstream calcs (FreeStream)
    double precision :: riem1Free, riem2Free, velFree

    ! Inlet BC Variables
    double precision, dimension(iCell, jCell) :: alpha
    double precision, dimension(iCell, jCell) :: c, p, rho
    double precision, dimension(iCell, jCell) :: riem1, riem2 
    double precision, dimension(iCell, jCell) :: m, u, v, vel, s, E, iE
    double precision, dimension(iCell, jCell) :: aCell
    double precision, dimension(iCell, jCell, 4) :: res

    ! Define State Vector Dimensions
    double precision, dimension(iCell, jCell, 4) :: qCurrent, qNext ! rho, rho*u, rho*v, rho*ener (3D Array)
    double precision, dimension(4) :: qFree
    double precision :: dqMax
    real :: alphaFree = 0.0 ! Free stream flow is horizontal

    integer :: ma

    print *, "TEST"
    
    do ma = 1, size(mach) ! Iterate through mach values 0.3, 0.5, 0.7

        qNext = 0

        ! Define Free Stream state vector
        qFree(1) = 1 ! rho
        qFree(2) = ma*cos(alphaFree) ! x-mom
        qFree(3) = ma*sin(alphaFree) ! y-mom
        qFree(4) = (1/(gam*(gam-1))) + 0.5*(ma**2)

        iterations = 1
        iterate = .TRUE.
        do while (iterate)

            ! =====================
            ! Inlet BCs -> Option 1
            ! =====================

            riem1Free = velFree + (2*sound)/(gam - 1)
            riem2Free = velFree + (2*sound)/(gam - 1)
            velFree = ma*sound
            
            ! (i), (ii) Riemman free-stream conditions for inlet
            do j = 1, jCell
                    riem1(1, j) = riem1Free
                    riem2(1, j) = riem2Free
            end do

            ! (iii) -> V(1, j) velocity vector @ inlet
            do j = 1, jCell
                vel(1, j) = 0.5*(riem1(1, j) + riem2(1, j)) ! velocity vector
            end do

            ! (iv) -> u(1, j) and v(1, j)
            do j = 1, jCell
                u(1, j) = vel(1, j)*cos(alpha(1, j))
                v(1, j) = vel(1, j)*sin(alpha(1, j))
            end do

            ! (v) -> Speed of sound (c)
            do j = 1, jCell
                c(1, j) = 0.25*(gam - 1)*(riem1(1, j) - riem2(1, j))
            end do

            ! (vi) -> mach number
            do j = 1, jCell
                m(1, j) = vel(1, j)/c(1, j)
            end do

            ! (vii) -> static pressure
            do j = 1, jCell
                p(1, j) = pFree/((1 + 0.5*(gam - 1)*(m(1, j)**2)))**(gam/(gam - 1))
            end do

            ! (viii) -> density
            do j = 1, jCell
                rho(1, j) = gam*p(1, j)/(c(1, j)**2)
            end do

            ! ==========================================
            ! Exit BCs -> Using characteristic variables
            ! ==========================================

            ! (i) -> Static Pressure
            do j = 1, jCell
                if (ma < 1) then
                    p(iCell, j) = pFree
                else if (ma >= 1) then
                    p(iCell, j) = p(iCell-1, j)
                end if
            end do

            ! (ii) -> y-component of velocity
            do j = 1, jCell
                v(iCell, j) = v(iCell-1, j)
            end do

            ! (iii) -> riemman
            do j = 1, jCell
                riem1(iCell, j) = riem1(iCell-1, j)
            end do

            ! (iv) -> entropy
            do j = 1, jCell
                s(iCell, j) = s(iCell-1, j)
            end do

            ! Internal Energy
            do i = 1, iCell
                do j = 1, jCell
                    iE(i, j) = cV*(p(i, j))/(rho(i, j)*rS)
                end do
            end do

            ! ===================================================================================
            ! Uniform Initialization -> Set all values between inlet and exit to the inlet values
            ! ===================================================================================
            if (iterations == 1) then
                do i = 1, iCell
                    do j = 1, jCell
                        riem1(i, j) = riem1(1, j)
                        riem2(i, j) = riem2(1, j)
                        rho(i, j) = rho(1, j)
                        u(i, j) = u(1, j)
                        v(i, j) = v(1, j)
                        E(i, j) = iE(i, j) + 0.5*vel(i, j)**2 ! Total Energy
                        qCurrent(i, j, 1) = qFree(1)
                        qCurrent(i, j, 2) = qFree(2)
                        qCurrent(i, j, 3) = qFree(3)
                        qCurrent(i, j, 4) = qFree(4)
                    end do
                end do
            end if

            call cellDimensions(x, y, aCell)
            call temporalDiscretization(aCell, res, qCurrent, x, y, u, v)

            ! Stores largest change of any element within the state vector across all cells
            dqMax = 0
            do i = 1, iCell
                do j = 1, jCell
                    do k = 1, 4
                        if (res(i, j, k) > dqMax) then
                            dqMax = res(i, j, k)
                        end if
                    end do
                end do
            end do

            ! Checks if maximum change in state vector is within tolerance requirements for solver
            if (dqMax < solTol) then
                iterate = .FALSE.
            end if

            print *, "Iteration:", iterations
            iterations = iterations + 1

        end do

    end do

end subroutine eulerSolver
