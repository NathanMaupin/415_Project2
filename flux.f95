subroutine flux(xMesh, yMesh, u, v, q)
    
    use common_parameters
    implicit none

    double precision, dimension(iCell, jCell, 4) :: R ! rho*u, rho*u^2+p, rho*u*v, rho*H*u
    double precision, dimension(giCell, gjCell, 4):: f, g ! dimension 4 refers to N, S, E, & W
    double precision, dimension(iMax, jMax), intent(in) :: xMesh, yMesh
    double precision, dimension(iCell, jCell) :: u, v, rho, e
    double precision, dimension(iCell, jCell, 4) :: q
    double precision, dimension(iCell, jCell) :: xCell, yCell
    double precision, dimension(iCell, jCell, 4) :: xnode, ynode
    double precision, dimension(iCell, jCell, 4) :: xFace, yFace

    ! Get A, B, C, & D of cell
    do i = 1, iCell
        do j = 1, jCell
            ! Cell Centroid Coordinates
            xCell(i, j) = (xMesh(i, j) + xMesh(i+1, j) + xMesh(i, j+1) + xMesh(i+1, j+1))/4
            yCell(i, j) = (yMesh(i, j) + yMesh(i+1, j) + yMesh(i, j+1) + yMesh(i+1, j+1))/4
        end do
    end do
    do i = 1, iCell
        do j = 1, jCell
                xnode(i, j, 1) =  xMesh(i, j) ! A
                ynode(i, j, 1) =  yMesh(i, j)
                xnode(i, j, 2) = xMesh(i+1, j) ! B
                ynode(i, j, 2) = yMesh(i+1, j)
                xnode(i, j, 3) = xMesh(i, j+1) ! C
                ynode(i, j, 3) = yMesh(i, j+1)
                xnode(i, j, 4) = xMesh(i+1, j+1) ! D
                ynode(i, j, 4) = yMesh(i+1, j+1)

                xFace(i, j, 1) = xMesh(i+1, j+1) - xMesh(i, j+1) ! N -> C-D
                xFace(i, j, 2) = xMesh(i+1, j) - xMesh(i, j) ! S -> B-A
                xFace(i, j, 3) = xMesh(i+1, j+1) - xMesh(i+1, j) ! E -> C-B
                xFace(i, j, 4) = xMesh(i, j+1) - xMesh(i, j) ! W -> D-A
        end do
    end do

    do i = 1, iCell
        do j = 1, jCell
            rho(i, j) = q(i, j, 1)
            e(i, j) = q(i, j, 4)/rho(i,j)
        end do
    end do

    call ghostCell(u, v, rho, e)

    do i = 3, (giCell-2)
        do j = 3, (gjCell-2)
                f(i, j, 1) = 0.5*(f(i, j, 1) + f(i, j+1, 1)) ! North
                f(i, j, 2) = 0.5*(f(i, j, 2) + f(i, j-1, 2)) ! South
                f(i, j, 3) = 0.5*(f(i, j, 3) + f(i+1, j, 3)) ! East
                f(i, j, 4) = 0.5*(f(i, j, 4) + f(i-1, j, 4)) ! West
                g(i, j, 1) = 0.5*(g(i, j, 1) + g(i, j+1, 1)) ! North
                g(i, j, 2) = 0.5*(g(i, j, 2) + g(i, j-1, 2)) ! South
                g(i, j, 3) = 0.5*(g(i, j, 3) + g(i+1, j, 3)) ! East
                g(i, j, 4) = 0.5*(g(i, j, 4) + g(i-1, j, 4)) ! West
        end do
    end do

    do i = 3, (giCell-2)
        do j = 3, (gjCell-2)
            R = 0.5*(f(i, j, 3)+f(i+1, j, 3))*(yMesh(i+1, j+1)-yMesh(i+1,j)) - &
                0.5*(g(i, j, 3)+g(i+1, j, 3))*(xMesh(i+1, j+1)-xMesh(i+1, j)) + &
                0.5*(f(i, j, 1)+f(i, j+1, 1))*(yMesh(i+1, j+1)-yMesh(i,j+1)) - &
                0.5*(g(i, j, 1)+g(i, j+1, 1))*(xMesh(i+1, j+1)-xMesh(i,j+1)) + &
                0.5*(f(i, j, 4)+f(i-1, j, 4))*(yMesh(i, j+1)-yMesh(i, j)) - &
                0.5*(g(i, j, 4)+g(i-1, j, 4))*(xMesh(i, j+1)-xMesh(i, j)) + &
                0.5*(f(i, j, 2)+f(i, j-1, 2))*(yMesh(i+1, j)-yMesh(i, j)) - &
                0.5*(g(i, j, 2)+g(i, j-1, 2))*(xMesh(i+1, j)-xMesh(i, j))
        end do
    end do

end subroutine flux
