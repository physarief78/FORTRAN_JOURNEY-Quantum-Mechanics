program quantum_harmonic_oscillator
  implicit none
  integer, parameter :: N = 200              ! Number of grid points
  real(8), parameter :: xmax = 5.0d0        ! Left boundary of the grid
  real(8), parameter :: xmin = -xmax         ! Right boundary of the grid
  real(8), parameter :: dx = (xmax - xmin) / (N - 1)  ! Grid spacing
  real(8), dimension(N) :: x                 ! Position grid
  real(8), dimension(N, N) :: H              ! Hamiltonian matrix
  real(8), dimension(N) :: eigenvalues       ! Eigenvalues
  real(8), dimension(N, N) :: eigenvectors   ! Eigenvectors
  integer :: i, j, nlevel
  character(len=20) :: filename

  filename = 'energy_levels_data.txt'

  ! Initialize position grid
  do i = 1, N
    x(i) = 0.0d0 + (i - N/2) * dx
  end do

  ! Construct Hamiltonian matrix for the harmonic oscillator
  do i = 1, N
    do j = 1, N
      if (i == j) then
        ! Diagonal terms: kinetic + potential energy
        H(i, j) = -2.0d0 / (dx * dx) + 0.5d0 * x(i)**2
      else if (abs(i - j) == 1) then
        ! Off-diagonal terms for finite difference approximation
        H(i, j) = 1.0d0 / (dx * dx)
      else
        H(i, j) = 0.0d0
      end if
    end do
  end do

  ! Compute eigenvalues and eigenvectors using the QR algorithm
  call qr_algorithm(H, N, eigenvalues, eigenvectors)

  ! Open file to write data
  open(unit=10, file=filename, status='unknown')

  ! Write data for ground state to fourth energy level
  do nlevel = 1, 5
    write(10, *) "Energy level ", nlevel-1, ": ", eigenvalues(nlevel)
    write(10, *) "Position  Probability Density for level ", nlevel-1
    do i = 1, N
      write(10, '(F12.6, F12.6)') x(i), eigenvectors(i, nlevel)**2
    end do
    write(10, *)
  end do

  ! Close the file
  close(10)

  print *, "Data for energy levels 0 to 4 saved to ", filename

contains

  subroutine qr_algorithm(A, n, eigenvalues, eigenvectors)
    ! QR algorithm for eigenvalue decomposition
    integer, intent(in) :: n
    real(8), intent(inout) :: A(n, n)
    real(8), intent(out) :: eigenvalues(n)
    real(8), intent(out) :: eigenvectors(n, n)
    real(8), dimension(n, n) :: Q, R, temp
    integer :: k, i, j, max_iter
    real(8) :: tol

    max_iter = 1000
    tol = 1.0d-8

    ! Initialize eigenvectors as identity matrix
    eigenvectors = 0.0d0
    do i = 1, n
      eigenvectors(i, i) = 1.0d0
    end do

    do k = 1, max_iter
      ! Step 1: Perform QR decomposition of A into Q and R
      call qr_decompose(A, Q, R, n)

      ! Step 2: Update A as R * Q, which shifts A towards upper triangular
      A = matmul(R, Q)

      ! Step 3: Update eigenvectors as eigenvectors * Q
      eigenvectors = matmul(eigenvectors, Q)

      ! Check for convergence by examining off-diagonal elements
      if (converged(A, n, tol)) exit
    end do

    ! Extract eigenvalues from the diagonal of A
    do i = 1, n
      eigenvalues(i) = A(i, i)
    end do
  end subroutine qr_algorithm

  subroutine qr_decompose(A, Q, R, n)
    ! QR decomposition using Gram-Schmidt orthogonalization
    integer, intent(in) :: n
    real(8), intent(in) :: A(n, n)
    real(8), intent(out) :: Q(n, n), R(n, n)
    integer :: i, j, k
    real(8) :: norm

    Q = 0.0d0
    R = 0.0d0

    ! Gram-Schmidt process
    do i = 1, n
      Q(:, i) = A(:, i)
      do j = 1, i - 1
        R(j, i) = dot_product(Q(:, j), A(:, i))
        Q(:, i) = Q(:, i) - R(j, i) * Q(:, j)
      end do
      norm = sqrt(dot_product(Q(:, i), Q(:, i)))
      R(i, i) = norm
      Q(:, i) = Q(:, i) / norm
    end do
  end subroutine qr_decompose

  logical function converged(A, n, tol)
    ! Check for convergence of the QR algorithm
    integer, intent(in) :: n
    real(8), intent(in) :: A(n, n)
    real(8), intent(in) :: tol
    integer :: i, j
    converged = .true.

    do i = 1, n - 1
      if (abs(A(i + 1, i)) > tol) then
        converged = .false.
        exit
      end if
    end do
  end function converged

end program quantum_harmonic_oscillator