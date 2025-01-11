program schrodinger_gaussian
    implicit none

    ! Define a kind for double precision
    integer, parameter :: dp = selected_real_kind(15, 307)

    ! Constants and parameters
    integer, parameter :: nx = 200          ! Number of spatial grid points
    integer, parameter :: nt = 1000         ! Number of time steps
    real(dp), parameter :: L = 3.0_dp       ! Length of the box
    real(dp), parameter :: dx = L / nx      ! Spatial step size
    real(dp), parameter :: dt = 0.001_dp    ! Time step size
    real(dp), parameter :: pi = 3.141592653589793_dp
    real(dp), parameter :: hbar = 1.0_dp, m = 1.0_dp
    real(dp), parameter :: omega = 10.0_dp  ! Angular frequency of the harmonic oscillator

    ! Arrays for the wavefunction and potential
    complex(dp), dimension(nx) :: psi, psi_new
    complex(dp), dimension(nx, nx) :: A, B
    real(dp), dimension(nx) :: x, V

    ! Variables
    integer :: i, j, k
    real(dp) :: x0, sigma, k0, norm
    complex(dp) :: imag_unit
    real(dp), dimension(nx) :: density

    ! Initialization
    imag_unit = (0.0_dp, 1.0_dp)            ! Imaginary unit
    x0 = 0.0_dp                             ! Initial position of the wavepacket
    sigma = 0.2_dp                          ! Width of the Gaussian wavepacket
    k0 = 5.0_dp                             ! Initial wave vector

    ! Set up spatial grid and potential
    do i = 1, nx
        x(i) = -L / 2.0_dp + (i - 1) * dx
        V(i) = 0.5_dp * m * omega**2 * x(i)**2  ! Harmonic oscillator potential
    end do

    ! Initialize wavefunction (Gaussian wavepacket)
    do i = 1, nx
        psi(i) = exp(-0.5_dp * ((x(i) - x0) / sigma)**2) * exp(imag_unit * k0 * x(i))
    end do

    ! Normalize wavefunction
    norm = 0.0_dp
    do i = 1, nx
        norm = norm + abs(psi(i))**2 * dx
    end do
    psi = psi / sqrt(norm)

    ! Construct matrices A and B for Crank-Nicolson scheme
    do i = 1, nx
        do j = 1, nx
            A(i, j) = (0.0_dp, 0.0_dp)
            B(i, j) = (0.0_dp, 0.0_dp)
        end do
    end do

    do i = 2, nx-1
        A(i, i-1) = -imag_unit * hbar * dt / (4.0_dp * m * dx**2)
        A(i, i) = 1.0_dp + imag_unit * hbar * dt / (2.0_dp * m * dx**2) + imag_unit * dt * V(i) / (2.0_dp * hbar)
        A(i, i+1) = -imag_unit * hbar * dt / (4.0_dp * m * dx**2)

        B(i, i-1) = imag_unit * hbar * dt / (4.0_dp * m * dx**2)
        B(i, i) = 1.0_dp - imag_unit * hbar * dt / (2.0_dp * m * dx**2) - imag_unit * dt * V(i) / (2.0_dp * hbar)
        B(i, i+1) = imag_unit * hbar * dt / (4.0_dp * m * dx**2)
    end do

    A(1, 1) = 1.0_dp
    A(nx, nx) = 1.0_dp
    B(1, 1) = 1.0_dp
    B(nx, nx) = 1.0_dp

    ! Time evolution loop
    open(unit=10, file="wavepacket_harmonic.dat", status="unknown")
    do k = 1, nt
        ! Solve the linear system A * psi_new = B * psi
        call solve_linear_system(A, B, psi, psi_new, nx)

        ! Update wavefunction
        psi = psi_new

        ! Calculate probability density and write to file
        write(10, '(F10.5)', advance="no") k * dt
        do i = 1, nx
            density(i) = abs(psi(i))**2
            write(10, '(F10.5)', advance="no") density(i)
        end do
        write(10, *)
    end do
    close(10)

contains

    subroutine solve_linear_system(A, B, psi, psi_new, n)
        complex(dp), dimension(n, n), intent(in) :: A, B
        complex(dp), dimension(n), intent(in) :: psi
        complex(dp), dimension(n), intent(out) :: psi_new
        integer, intent(in) :: n
        complex(dp), dimension(n) :: rhs, diag, lower, upper
        complex(dp) :: temp
        integer :: i

        ! Compute the right-hand side: B * psi
        do i = 1, n
            rhs(i) = (0.0_dp, 0.0_dp)
            do j = 1, n
                rhs(i) = rhs(i) + B(i, j) * psi(j)
            end do
        end do

        ! Decompose A into tridiagonal components
        lower = (0.0_dp, 0.0_dp)
        diag = (0.0_dp, 0.0_dp)
        upper = (0.0_dp, 0.0_dp)
        do i = 2, n
            lower(i) = A(i, i-1)
        end do
        do i = 1, n
            diag(i) = A(i, i)
        end do
        do i = 1, n-1
            upper(i) = A(i, i+1)
        end do

        ! Forward elimination
        do i = 2, n
            temp = lower(i) / diag(i-1)
            diag(i) = diag(i) - temp * upper(i-1)
            rhs(i) = rhs(i) - temp * rhs(i-1)
        end do

        ! Back substitution
        psi_new(n) = rhs(n) / diag(n)
        do i = n-1, 1, -1
            psi_new(i) = (rhs(i) - upper(i) * psi_new(i+1)) / diag(i)
        end do
    end subroutine solve_linear_system

end program schrodinger_gaussian
