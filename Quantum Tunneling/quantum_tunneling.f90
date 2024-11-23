program quantum_tunneling
    implicit none
    integer, parameter :: nx = 500          ! Number of spatial points
    integer, parameter :: nt = 1000         ! Number of time steps
    real(8), parameter :: L = 10.0          ! Spatial domain (-L, L)
    real(8), parameter :: dx = 2*L/nx       ! Spatial step
    real(8), parameter :: dt = 0.01         ! Time step
    real(8), parameter :: k0 = -5.0          ! Initial wave vector
    real(8), parameter :: x0 = -5.0         ! Initial position of the wave packet
    real(8), parameter :: sigma = 1.0       ! Width of the wave packet
    real(8), parameter :: V0 = 13.0         ! Height of the potential barrier
    real(8), parameter :: hbar = 1.0        ! Planck's constant
    real(8), parameter :: mass = 1.0        ! Particle mass
    real(8), parameter :: C = 20.0          ! Amplitude of the wavepacket


    complex(8), dimension(nx) :: psi        ! Wave function
    real(8), dimension(nx) :: x, V          ! Spatial points and potential
    real(8), dimension(nx) :: probability
    complex(8), dimension(nx, nx) :: H      ! Hamiltonian matrix
    complex(8), dimension(nx) :: b          ! Right-hand side for Crank-Nicholson
    complex(8), dimension(nx) :: temp       ! Temporary array for solving
    integer :: i, j, t

    complex(8), parameter :: I_complex = (0.0, 1.0) ! Imaginary unit
    real(8) :: a                                ! Coefficient for Hamiltonian

    ! Initialize spatial points, potential, and coefficient
    do i = 1, nx
        x(i) = -L + (i - 1) * dx
        if (abs(x(i)) < 0.5) then
            V(i) = V0
        else
            V(i) = 0.0
        end if
    end do
    a = hbar / (2.0 * mass * dx**2)

    ! Initialize the wave function (Gaussian wave packet)
    do i = 1, nx
        psi(i) = C * exp(-((x(i) - x0)**2) / (2.0 * sigma**2)) * exp(I_complex * k0 * x(i))
    end do

    ! Construct the Hamiltonian matrix for Crank-Nicholson
    do i = 1, nx
        do j = 1, nx
            if (i == j) then
                H(i, j) = 1.0 - I_complex * dt / hbar * (2.0 * a + V(i))
            else if (abs(i - j) == 1) then
                H(i, j) = I_complex * dt / hbar * a
            else
                H(i, j) = 0.0
            end if
        end do
    end do

    ! Time evolution loop
    do t = 1, nt
        ! Save probability density at each time step
        probability = abs(psi)**2
        call save_data(t, x, probability)

        ! Calculate the right-hand side (b) for Crank-Nicholson
        do i = 1, nx
            b(i) = (1.0 + I_complex * dt / hbar * (2.0 * a + V(i))) * psi(i)
            if (i > 1) b(i) = b(i) - (I_complex * dt / hbar * a) * psi(i - 1)
            if (i < nx) b(i) = b(i) - (I_complex * dt / hbar * a) * psi(i + 1)
        end do

        ! Solve the linear system H * psi = b using Thomas algorithm
        call thomas_algorithm(H, b, temp, nx)
        psi = temp
    end do

contains

    subroutine thomas_algorithm(H, b, x, n)
        ! Solves H * x = b for a tridiagonal matrix H using Thomas algorithm
        complex(8), intent(in) :: H(n, n)
        complex(8), intent(in) :: b(n)
        complex(8), intent(out) :: x(n)
        integer, intent(in) :: n
        complex(8), dimension(n) :: c_prime, d_prime
        integer :: i

        ! Forward sweep
        c_prime(1) = H(1, 2) / H(1, 1)
        d_prime(1) = b(1) / H(1, 1)
        do i = 2, n - 1
            c_prime(i) = H(i, i + 1) / (H(i, i) - H(i, i - 1) * c_prime(i - 1))
            d_prime(i) = (b(i) - H(i, i - 1) * d_prime(i - 1)) / &
                         (H(i, i) - H(i, i - 1) * c_prime(i - 1))
        end do
        d_prime(n) = (b(n) - H(n, n - 1) * d_prime(n - 1)) / &
                     (H(n, n) - H(n, n - 1) * c_prime(n - 1))

        ! Back substitution
        x(n) = d_prime(n)
        do i = n - 1, 1, -1
            x(i) = d_prime(i) - c_prime(i) * x(i + 1)
        end do
    end subroutine thomas_algorithm

    subroutine save_data(step, x, prob)
        integer, intent(in) :: step
        real(8), intent(in) :: x(nx), prob(nx)
        character(len=20) :: filename
        integer :: i
        write(filename, '(A,I4.4,A)') 'wave_packet_', step, '.dat'
        open(unit=10, file=trim(filename), status='unknown')
        do i = 1, nx
            write(10, '(F10.5, F10.5)') x(i), prob(i)
        end do
        close(10)
    end subroutine save_data

end program quantum_tunneling
