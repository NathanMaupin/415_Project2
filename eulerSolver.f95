subroutine eulerSolver(x, y)
    use common_parameters
    implicit none

    ! From main program
    double precision, dimension(iMax, jMax), intent(in):: x, y

    ! Upstream calcs (FreeStream)
    double precision :: velFree

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
    double precision, dimension(iCell, jCell) :: qNext1, qCurrent1

    integer :: ma

    do ma = 1, size(mach) ! Iterate through mach values 0.3, 0.5, 0.7

        ! Define Free Stream state vector
        qFree(1) = 1 ! rho
        qFree(2) = ma*cos(alphaFree) ! x-mom
        qFree(3) = ma*sin(alphaFree) ! y-mom
        qFree(4) = (1/(gam*(gam-1))) + 0.5*(ma**2)

        iterations = 1
        iterate = .TRUE.

        do i = 1, iCell
            do j = 1, jCell
                qNext(i, j, 1) = qFree(1)
                qNext(i, j, 2) = qFree(2)
                qNext(i, j, 3) = qFree(3)
                qNext(i, j, 4) = qFree(4)
            end do
        end do

        do while (iterate)
            ! =====================
            ! Inlet BCs -> Option 1
            ! =====================
            velFree = ma*cFree

            ! Find Pressure
            do j = 1, jCell
                p(1, j) = (gam-1)*(qCurrent(1, j, 4) - 0.5*(qCurrent(1, j, 2)**2 + qCurrent(1, j, 3)**2)/qCurrent(1, j, 1))
            end do

            ! (i), (ii) Riemman free-stream conditions for inlet
            do j = 1, jCell
                riem1 = velFree - (2/(gam-1))*sqrt(gam*(cFree))
                riem2 = sqrt(qCurrent(1, j, 2)**2 + qCurrent(1, j, 3)**2) - (2/(gam-1))*(sqrt(gam*(p(1, j))))
            end do

            ! Copy 
            do i = 1, iCell
                do j = 1, jCell
                    do k = 1, 4
                        qCurrent(i, j, k) = qNext(i, j, k)
                    end do
                end do
            end do

            call cellDimensions(x, y, aCell)
            call temporalDiscretization(aCell, res, qNext, x, y, u, v)

            ! Stores largest change of any element within the state vector across all cells
            dqMax = 0
            do i = 1, iCell
                do j = 1, jCell
                    do k = 1, 4
                        if (res(i, j, k) > dqMax) then
                            print *, 'res:', res(i, j, k)
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


    do i = 1, iCell
        do j = 1, jCell
            qNext1(i, j) = qNext(i, j, 1)
            qCurrent1(i, j) =  qCurrent(i, j, 1)
        end do
    end do
    

    call csvExport(qCurrent1, qNext1, iCell, jCell, 'q1_q2_final')

end subroutine eulerSolver
