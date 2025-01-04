program schrodinger_gaussian_2d
    implicit none

    ! Define a kind for double precision
    integer, parameter :: dp = selected_real_kind(15, 307)

    ! Constants and parameters
    integer, parameter :: nx = 100, ny = 100    ! Number of grid points in x and y
    integer, parameter :: nt = 1500             ! Number of time steps
    real(dp), parameter :: Lx = 2.0_dp, Ly = 2.0_dp  ! Dimensions of the box
    real(dp), parameter :: dx = Lx / nx, dy = Ly / ny ! Grid spacings
    real(dp), parameter :: dt = 0.00001_dp       ! Time step size
    real(dp), parameter :: pi = 3.141592653589793_dp
    real(dp), parameter :: hbar = 1.0_dp, m = 1.0_dp

    ! Parameters for the double-slit potential
    real(dp), parameter :: slit_width = 0.1_dp       ! Width of each slit
    real(dp), parameter :: slit_separation = 0.2_dp  ! Separation between the centers of the slits
    real(dp), parameter :: barrier_height = 5000.0_dp! Height of the barrier

    ! Arrays for the wavefunction and potential
    complex(dp), dimension(nx, ny) :: psi, psi_temp, psi_new
    real(dp), dimension(nx, ny) :: x, y, V
    real(dp), dimension(nx, ny) :: density

    ! Variables
    integer :: i, j, k
    real(dp) :: x0, y0, sigma, kx, ky, norm
    complex(dp) :: imag_unit

    ! Initialization
    imag_unit = (0.0_dp, 1.0_dp)               ! Imaginary unit
    x0 = 0.0_dp                          ! Initial x position of the wavepacket
    y0 = -Ly / 4.0_dp                    ! Initial y position of the wavepacket
    sigma = 0.2_dp                             ! Width of the Gaussian wavepacket
    kx = 0.0_dp * pi                          ! Initial wave vector in x
    ky = 30.0_dp * pi                         ! Initial wave vector in y

    ! Set up spatial grid and potential
    do i = 1, nx
        do j = 1, ny
            x(i, j) = -Lx / 2.0_dp + (i - 1) * dx
            y(i, j) = -Ly / 2.0_dp + (j - 1) * dy

            ! Initialize the potential to zero
            V(i, j) = 0.0_dp

            ! Add double-slit potential at x = 0
            if (abs(y(i, j)) < dy) then
                if ((abs(x(i, j) - slit_separation / 4.0_dp) > slit_width / 8.0_dp) .and. &
                    (abs(x(i, j) + slit_separation / 4.0_dp) > slit_width / 8.0_dp)) then
                    V(i, j) = barrier_height
                end if
            end if
        end do
    end do

    ! Initialize wavefunction (Gaussian wavepacket)
    do i = 1, nx
        do j = 1, ny
            psi(i, j) = exp(-0.5_dp * (((x(i, j) - x0) / sigma)**2 + ((y(i, j) - y0) / sigma)**2)) * &
                        exp(imag_unit * (kx * x(i, j) + ky * y(i, j)))
        end do
    end do

    ! Normalize wavefunction
    norm = 0.0_dp
    do i = 1, nx
        do j = 1, ny
            norm = norm + abs(psi(i, j))**2 * dx * dy
        end do
    end do
    psi = psi / sqrt(norm)

    ! Time evolution loop
    open(unit=10, file="wavepacket_animation_2d.dat", status="unknown")
    do k = 1, nt
        ! Step 1: Solve in x-direction
        call crank_nicolson_x(psi, psi_temp, V, nx, ny, dx, dt, hbar, m)

        ! Step 2: Solve in y-direction
        call crank_nicolson_y(psi_temp, psi_new, V, nx, ny, dy, dt, hbar, m)

        ! Update wavefunction
        psi = psi_new

        ! Calculate probability density and write to file
        write(10, '(F10.5)', advance="no") k * dt
        do i = 1, nx
            do j = 1, ny
                density(i, j) = abs(psi(i, j))**2
                write(10, '(F10.5)', advance="no") density(i, j)
            end do
        end do
        write(10, *)
    end do
    close(10)

contains

    subroutine crank_nicolson_x(psi_in, psi_out, V, nx, ny, dx, dt, hbar, m)
        implicit none
        complex(dp), intent(in) :: psi_in(nx, ny)
        complex(dp), intent(out) :: psi_out(nx, ny)
        real(dp), intent(in) :: V(nx, ny), dx, dt, hbar, m
        integer, intent(in) :: nx, ny
        complex(dp), dimension(nx) :: rhs, diag, lower, upper, psi_temp
        complex(dp) :: imag_unit, factor
        integer :: i, j

        imag_unit = (0.0_dp, 1.0_dp)
        factor = imag_unit * hbar * dt / (2.0_dp * m * dx**2)

        do j = 1, ny
            ! Set up tridiagonal system for x-direction
            diag = (1.0_dp + 2.0_dp * factor)
            lower = -factor
            upper = -factor

            ! Include potential term
            do i = 1, nx
                diag(i) = diag(i) + imag_unit * dt / (2.0_dp * hbar) * V(i, j)
            end do

            ! Right-hand side
            do i = 2, nx-1
                rhs(i) = (1.0_dp - 2.0_dp * factor) * psi_in(i, j) + factor * (psi_in(i-1, j) + psi_in(i+1, j))
                rhs(i) = rhs(i) - imag_unit * dt / (2.0_dp * hbar) * V(i, j) * psi_in(i, j)
            end do
            rhs(1) = psi_in(1, j)
            rhs(nx) = psi_in(nx, j)

            ! Solve tridiagonal system
            call tridiagonal_solver(lower, diag, upper, rhs, psi_temp, nx)

            ! Store results
            psi_out(:, j) = psi_temp
        end do
    end subroutine crank_nicolson_x

    subroutine crank_nicolson_y(psi_in, psi_out, V, nx, ny, dy, dt, hbar, m)
        implicit none
        complex(dp), intent(in) :: psi_in(nx, ny)
        complex(dp), intent(out) :: psi_out(nx, ny)
        real(dp), intent(in) :: V(nx, ny), dy, dt, hbar, m
        integer, intent(in) :: nx, ny
        complex(dp), dimension(ny) :: rhs, diag, lower, upper, psi_temp
        complex(dp) :: imag_unit, factor
        integer :: i, j

        imag_unit = (0.0_dp, 1.0_dp)
        factor = imag_unit * hbar * dt / (2.0_dp * m * dy**2)

        do i = 1, nx
            ! Set up tridiagonal system for y-direction
            diag = (1.0_dp + 2.0_dp * factor)
            lower = -factor
            upper = -factor

            ! Include potential term
            do j = 1, ny
                diag(j) = diag(j) + imag_unit * dt / (2.0_dp * hbar) * V(i, j)
            end do

            ! Right-hand side
            do j = 2, ny-1
                rhs(j) = (1.0_dp - 2.0_dp * factor) * psi_in(i, j) + factor * (psi_in(i, j-1) + psi_in(i, j+1))
                rhs(j) = rhs(j) - imag_unit * dt / (2.0_dp * hbar) * V(i, j) * psi_in(i, j)
            end do
            rhs(1) = psi_in(i, 1)
            rhs(ny) = psi_in(i, ny)

            ! Solve tridiagonal system
            call tridiagonal_solver(lower, diag, upper, rhs, psi_temp, ny)

            ! Store results
            psi_out(i, :) = psi_temp
        end do
    end subroutine crank_nicolson_y

    subroutine tridiagonal_solver(lower, diag, upper, rhs, solution, n)
        implicit none
        integer, intent(in) :: n
        complex(dp), intent(in) :: lower(n), diag(n), upper(n), rhs(n)
        complex(dp), intent(out) :: solution(n)
        complex(dp) :: temp(n)
        integer :: i

        ! Forward elimination
        temp(1) = diag(1)
        solution(1) = rhs(1)
        do i = 2, n
            temp(i) = diag(i) - lower(i) * upper(i-1) / temp(i-1)
            solution(i) = rhs(i) - lower(i) * solution(i-1) / temp(i-1)
        end do

        ! Back substitution
        solution(n) = solution(n) / temp(n)
        do i = n-1, 1, -1
            solution(i) = (solution(i) - upper(i) * solution(i+1)) / temp(i)
        end do
    end subroutine tridiagonal_solver

end program schrodinger_gaussian_2d
