module common_parameters
    implicit none

    ! ==========
    ! Input Data
    ! ==========
    integer, parameter :: iMax = 10, jMax = 4 ! iMax/jMax = 5
    integer, parameter :: iCell = iMax - 1, jCell = jMax - 1
    integer, parameter :: cellCount = iCell*jCell
    ! Simulation conditions
    real, dimension(3) :: mach = [0.3, 0.5, 0.7] ! Initial mach conditions
    double precision :: heatRatio = 1.4 ! Specific heat ratio for earth's atmosphere
    real :: sound = 340.29 ! Speed of sound at Sea Level (m/s)
    real :: aoa = 0 ! Angle of flow (deg)
    real :: cV = 0.816 ! Specific Volume of Air -> (m3/kg) 
    real :: rS = 0.287 ! Specific Gas Constant ->(J/(kg*K))
    double precision, parameter :: PI = 4.0D0*ATAN(1.0D0)



    ! Convergence Criteria
    double precision, parameter :: meshTol = 1e-7
    double precision, parameter :: solTol = 1e-10
    logical :: iterate 
    integer :: iterations

    ! Boundary Conditions for Grid Division
    integer :: x_l = 0, x_r = 5 ! Left and Right BCs

    ! Loop variables
    integer :: i, j, k

end module common_parameters
program main
    use common_parameters
    implicit none


    ! Subroutine arrays
    double precision, dimension(iMax, jMax) :: ep_array, eta_array, x_array, y_array
    double precision, dimension(iMax, jMax) :: x_next, y_next
    double precision, dimension(iMax, jMax) :: x_current, y_current
    double precision, dimension(iCell, jCell) :: xCell, yCell, aCell


    ! Values for fluid computations

    ! =========================================
    ! Subroutine for generating meshes
    ! =========================================
    call meshGeneration(ep_array, eta_array, x_array, y_array, x_next, y_next, x_current, y_current)

    call cellDimensions(x_next, y_next, xCell, yCell, aCell)
    ! =========================================
    ! Exporting Meshes
    ! =========================================
    call csvExport(ep_array, eta_array, iMax, jMax, 'ep_eta_grid')
    call csvExport(x_array, y_array, iMax, jMax, 'x_y_preliminary_Grid')
    call csvExport(x_next, y_next, iMax, jMax, 'x_y_refined_Grid')
    call csvExport(xCell, yCell, iCell, jCell, 'cellCoordinates')
    call csvExport(aCell, aCell, iCell, jCell, 'cellArea')


    call cellDimensions(x_next, y_next, xCell, yCell, aCell)
    print *, y_next

    ! =========================================
    ! Euler Solver
    ! =========================================
    do i = 1, size(mach)
        ! call eulerSolver(iMax, jMax, x_next, y_next, eulerTol, mach, gamma, sound, aoa)
    end do

end program main

subroutine meshGeneration(ep_array, eta_array, x_array, y_array, x_next, y_next, x_current, y_current)
    
    use common_parameters

    implicit none
    
    double precision :: y_l, y_u

    ! ep vs. eta and x vs. y simple meshes
    double precision, dimension(iMax, jMax), intent(out) :: ep_array, eta_array, x_array, y_array
    double precision :: ep, eta, x, y

    ! Variables for iterations for refining Mesh
    double precision :: dep, deta
    double precision, dimension(iMax, jMax), intent(out) :: x_next, y_next
    double precision, dimension(iMax, jMax), intent(out) :: x_current, y_current
    double precision, dimension(iMax, jMax) :: dx, dy
    integer :: row
    double precision :: alpha
    double precision :: beta
    double precision :: gamma
    double precision :: step
    double precision :: a1, a2, a3, a4, a5, a6, a7, a8

    ! Mesh Tolerance Check
    double precision :: dx_max, dy_max

    ! Loops for collecting ep and eta values and converting them to x-y grid
    do j = 1, jMax
        do i = 1, iMax
            ep = dble(i - 1) / dble(iMax - 1) ! Convert to double precision to ensure floating-point division
            x = x_l + ep * (x_r - x_l)
            
            if(0.4*iMax < i .and. i <= 0.6*iMax) then
                y_l = 0.1*sin((x - 2)*PI)
                y_u = 1 - 0.1*sin((x - 2)*PI)
            else
                y_l = 0.0
                y_u = 1.0
            end if

            eta = dble(j - 1) / dble(jMax - 1)
            y = y_l + eta*(y_u - y_l)

            ep_array(i, j) = ep
            eta_array(i, j) = eta
            x_array(i, j) = x
            y_array(i, j) = y

        end do
    end do

       ! =========================================
    ! Iterations - Refining Mesh
    ! =========================================

    iterate = .TRUE.
    iterations = 0
    
    dep = ep_array(3, 1) - ep_array(2, 1)
    deta = eta_array(1, 3) - eta_array(1, 2)

    x_next = (x_array)
    y_next = y_array

    do while (iterate)
        iterations = iterations + 1
        
        ! Perform a deep copy from x_next to x_current and y_next to y_current
        do i = 1, iMax
            do j = 1, jMax
                x_current(i, j) = x_next(i, j)
                y_current(i, j) = y_next(i, j)
            end do
        end do

        ! Dirichlet BCs
        x_next(1, 1) = x_current(1, 1)
        y_next(1, 1) = y_current(1, 1)
        x_next(1, jMax-1) = x_current(1, jMax-1)
        y_next(1, jMax-1) = y_current(1, jMax-1)
        x_next(iMax-1, 1) = x_current(iMax-1, 1)
        y_next(iMax-1, 1) = y_current(iMax-1, 1)
        x_next(iMax - 1, jMax - 1) = x_current(iMax - 1, jMax - 1)
        y_next(iMax - 1, jMax - 1) = y_current(iMax - 1, jMax - 1)
        do row = 2, (jMax - 1)
            x_next(1, row) = x_current(1, row)
            y_next(1, row) = y_current(3, row)
            x_next(iMax - 1, row) = x_current(iMax - 1, row)
            y_next(iMax - 1, row) = y_current(iMax - 3, row)
        end do

        do i = 2, (iMax - 1)
            if (2 <= x_current(i, 1) .and. x_current(i, 1) <= 3) then
                x_next(i, 1) = x_current(i, 2) + 0.1*pi*cos((x_current(i, 1)-2)*pi)*(y_current(i, 2) - y_current(i, 1))
                y_next(i, 1) = 0.1*sin((x_next(i, 1)-2)*pi)
                x_next(i, jMax) = x_current(i, jMax-1) + 0.1*pi*cos((x_current(i, jMax)-2)*pi)*(y_current(i, jMax) - &
                 y_current(i, jMax-1))
                y_next(i, jMax-1) = 1 - 0.1*sin((x_next(i, jMax)-2)*pi)
            else
                x_next(i, 1) = x_current(i, 1)
                y_next(i, 1) = y_current(i, 1)
                x_next(i, jMax-2) = x_current(i, jMax-3)
                y_next(i, jMax-2) = y_current(i, jMax-2)
            end if

            do j = 2, (jMax - 1)
                ! Gauss-Seidel
                alpha = (1/(4*deta**2))*((x_current(i, j+1) - x_current(i, j-1))**2 + (y_current(i, j+1) - y_current(i, j-1))**2)
                beta = (1/(4*dep*deta))*((x_current(i+1, j) - x_current(i-1, j))*(x_current(i, j+1) - x_current(i, j-1)) + &
                (y_current(i+1, j) - y_current(i-1, j))*(y_current(i, j+1) - y_current(i, j-1)))
                gamma = (1/(4*dep**2))*((x_current(i+1, j) - x_current(i-1, j))**2 + (y_current(i+1, j) - y_current(i-1, j))**2)

                step = (2*((alpha/(dep**2)) + (gamma/(deta**2))))**(-1)

                a1 = (-beta*step)/(2*dep*deta)
                a2 = (gamma*step)/(deta**2)
                a3 = (beta*step)/(2*dep*deta)
                a4 = (alpha*step)/(dep**2)
                a5 = a4
                a6 = a3
                a7 = a2
                a8 = a1

                ! print *, 'x_next before:', x_next(i, j)

                x_next(i, j) = a1*x_next(i-1, j-1) + a2*x_next(i, j-1) + a3*x_next(i+1, j-1) + a4*x_next(i-1, j) + &
                a5*x_next(i+1, j) + a6*x_next(i-1, j+1) + a7*x_next(i, j+1) + a8*x_next(i+1, j+1)
                y_next(i, j) = a1*y_next(i-1, j-1) + a2*y_next(i, j-1) + a3*y_next(i+1, j-1) + a4*y_next(i-1, j) + &
                a5*y_next(i+1, j) + a6*y_next(i-1, j+1) + a7*y_next(i, j+1) + a8*y_next(i+1, j+1)

                ! print *, 'x_next after :', x_next(i, j)
            end do
        end do

        ! Verifying whether the tolerance requirement is satisfied.
        dx_max = 0
        dy_max = 0
        do j = 2, (jMax - 1)
            do i = 2, (iMax - 1)
                dx(i, j) = x_next(i, j) - x_current(i, j)
                dy(i, j) = y_next(i, j) - y_current(i, j)
            end do
        end do

        ! Updating tolerance check
        do j = 2, (jMax - 1)
            do i = 2, (iMax - 1)
                if (ABS(dx(i, j)) > dx_max) then
                    dx_max = dx(i, j)
                end if
                if (ABS(dy(i, j)) > dy_max) then
                    dy_max = dy(i, j)
                end if
            end do
        end do
        if (dx_max < meshTol .and. dy_max < meshTol) then
            iterate = .FALSE.
        end if
        if (iterations == 1000) then
            iterate = .TRUE.
        end if

        do i = 1, iMax
            do j = 1, jMax
                ! print *, 'x_current:', x_current(i, j)
                ! print *, 'x_next:   ', x_next(i, j)
            end do
        end do
        print *, 'ieration:', iterations
        print '(A, F5.3, A, F5.3, A)', '(dx, dy):(', dx_max, ', ', dy_max, ')'
    end do

end subroutine meshGeneration

subroutine cellDimensions(cell, xMesh, yMesh, xCell, yCell, aCell)
    
    use common_parameters
    implicit none
    
    integer, intent(in) :: cell ! Cell of Interest
    double precision, dimension(iMax, jMax), intent(in) :: xMesh, yMesh
    double precision, dimension(iCell, jCell), intent(out):: xCell, yCell, aCell
    double precision :: dx, dy
    double precision, dimension(iMax, jMax) :: alpha

    do i = 1, (iMax - 1)
        do j = 1, (jMax - 1)
            ! Cell Centroid Coordinates
            xCell(i, j) = (xMesh(i, j) + xMesh(i+1, j) + xMesh(i, j+1) + xMesh(i+1, j+1))/4
            yCell(i, j) = (yMesh(i, j) + yMesh(i+1, j) + yMesh(i, j+1) + yMesh(i+1, j+1))/4

            ! Cell Area
            aCell(i, j) = 0.5*((xMesh(i+1, j+1) - xMesh(i, j))*(yMesh(i, j+1) - yMesh(i+1, j)) - &
            (yMesh(i+1, j+1) - yMesh(i, j))*(xMesh(i, j+1) - xMesh(i+1, j)))

            ! Normal Vectors
            dx = xMesh(i + 1, j) - xMesh(i, j)
            dy = yMesh(i, j+1) - yMesh(i, j)
            alpha(i, j) = ATAN(dy/dx)
        end do
    end do

    do i = 1, cellCount
        ! Cell Centroid Coordinates
        xCell(i, j) = (xMesh(i, j) + xMesh(i+1, j) + xMesh(i, j+1) + xMesh(i+1, j+1))/4
        yCell(i, j) = (yMesh(i, j) + yMesh(i+1, j) + yMesh(i, j+1) + yMesh(i+1, j+1))/4
    end do
    
end subroutine cellDimensions

subroutine ghost(cell, qC, u, v)
    
    use common_parameters
    implicit none

    integer, intent(in) :: cell ! Cell of interest 
    double precision, intent(in) :: u, v ! Cell Velocity components
    double precision :: ug, vg ! Ghost Cell Velocities
    double precision, dimension(cellCount, 4) :: alphaC ! Stores alpha values for each face of cells
    double precision, dimension(cellCount, 4), intent(in) :: qC ! Referencing state vector for individual cell

    ! State Vector Components
    
end subroutine ghost

subroutine residual(xMesh, yMesh)
    
    use common_parameters
    implicit none

    double precision, dimension(iCell, jCell, 4) :: R ! rho*u, rho*u^2+p, rho*u*v, rho*H*u
    double precision, dimension(iMax, jMax, 4):: f, g ! dimension 4 refers to N, S, E, & W
    double precision, dimension(iMax, jMax), intent(in) :: xMesh, yMesh
    double precision, dimension(cellCount, 4) :: cellArray ! 1st dimenison is # of cells, 2nd dimension is N, S, E, & W faces
    
    double precision :: dx, dy
    integer :: A ! Cell Area


    ! Call 
    do i = 1, cellCount

    end do

    



    do k = 1, 4
        do i = 1, iMax-1
            do j = 1, jMax-1
                R(i, j, k) = 0.5*dx*(g(i, j+1, k) - g(i, j-1, k)) + 0.5*dy*(f(i+1, j, k) - f(i-1, j, k))
            end do
        end do
    end do

end subroutine residual

subroutine eulerSolver(x, y)
    use common_parameters
    implicit none

    ! From main program
    double precision, dimension(iMax, jMax), intent(in):: x, y

    ! Inlet BC Variables
    double precision, dimension(iCell, jCell) :: alpha
    double precision, dimension(iCell, jCell) :: c, p, rho
    double precision, dimension(iCell, jCell) :: riem1, riem2 
    double precision, dimension(iCell, jCell) :: xCell, yCell, aCell

    ! Define State Vector Dimensions
    double precision, dimension(iCell, jCell, 4) :: qCurrent, qNext, dq ! rho, rho*u, rho*v, rho*ener (3D Array)
    double precision :: dqMax

    ! Get cell dimensions
    call cellDimensions(x, y, xCell, yCell, aCell)

    ! =====================
    ! Inlet BCs -> Option 1
    ! =====================

    iterate = .TRUE.

    do while (iterate)

        do i = 1, cellCount
            ! Update one cell at a time
        end do

        if (dqMax < solTol) then
            iterate = .FALSE.
        end if
    end do


end subroutine eulerSolver

subroutine jst()

end subroutine jst
