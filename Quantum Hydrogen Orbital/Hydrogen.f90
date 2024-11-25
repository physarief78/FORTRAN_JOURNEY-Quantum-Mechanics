program hydrogen_orbital_animation
    implicit none
    ! Parameters for grid
    integer, parameter :: n_points = 100
    real(8), parameter :: xmin = -20.0d0, xmax = 20.0d0
    real(8), parameter :: ymin = -20.0d0, ymax = 20.0d0
    real(8), parameter :: pi = 3.141592653589793d0

    ! Variables
    real(8) :: x, y, r, phi, psi, normalization_factor
    real(8) :: dx, dy, radial_part, angular_part, sum_psi_squared
    integer :: i, j, l, m
    real(8), allocatable :: orbital_data(:,:)
    character(len=50) :: filename

    ! Grid spacing
    dx = (xmax - xmin) / (n_points - 1)
    dy = (ymax - ymin) / (n_points - 1)

    ! Allocate 2D array for orbital data
    allocate(orbital_data(n_points, n_points))

    ! Loop over quantum states
    do l = 0, 3
        do m = 0, l
            ! Initialize the sum for normalization
            sum_psi_squared = 0.0d0

            ! Loop over grid points
            do i = 1, n_points
                x = xmin + (i - 1) * dx
                do j = 1, n_points
                    y = ymin + (j - 1) * dy
                    r = sqrt(x**2 + y**2)       ! Radial distance
                    phi = atan2(y, x)          ! Azimuthal angle

                    ! Compute radial part (for n=4, l varies)
                    select case (l)
                    case (0)  ! l=0
                        radial_part = exp(-r / 4.0d0) * (1.0d0 - r / 4.0d0)
                    case (1)  ! l=1
                        radial_part = r * exp(-r / 4.0d0) * (1.0d0 - r / 8.0d0)
                    case (2)  ! l=2
                        radial_part = r**2 * exp(-r / 4.0d0) * (1.0d0 - r / 4.0d0 + r**2 / 72.0d0)
                    case (3)  ! l=3
                        radial_part = r**3 * exp(-r / 4.0d0) * (1.0d0 - 3.0d0*r / 4.0d0 + r**2 / 32.0d0)
                    end select

                    ! Compute angular part (spherical harmonics, θ=π/2 for 2D)
                    select case (m)
                    case (0)
                        angular_part = sqrt((2.0d0*l + 1.0d0) / (4.0d0*pi))  ! m=0
                    case (1)
                        angular_part = sqrt((2.0d0*l + 1.0d0) / (4.0d0*pi)) * sin(phi)  ! m=1
                    case (2)
                        angular_part = sqrt((2.0d0*l + 1.0d0) / (4.0d0*pi)) * cos(2.0d0 * phi)  ! m=2
                    case (3)
                        angular_part = sqrt((2.0d0*l + 1.0d0) / (4.0d0*pi)) * sin(3.0d0 * phi)  ! m=3
                    end select

                    ! Combine radial and angular parts
                    psi = radial_part * angular_part

                    ! Accumulate sum of squared wavefunction values for normalization
                    sum_psi_squared = sum_psi_squared + psi**2 * dx * dy

                    ! Store in data array
                    orbital_data(i, j) = psi
                end do
            end do

            ! Compute the normalization factor
            normalization_factor = sqrt(sum_psi_squared)

            ! Normalize the wavefunction
            do i = 1, n_points
                do j = 1, n_points
                    orbital_data(i, j) = orbital_data(i, j) / normalization_factor
                end do
            end do

            ! Write normalized orbital data to a file
            write(filename, '(A,I1,A,I1,A,I1,A)') "orbital_data_4_", l, "_", m, "_normalized.dat"
            open(unit=10, file=filename, status="replace")
            do i = 1, n_points
                x = xmin + (i - 1) * dx
                do j = 1, n_points
                    y = ymin + (j - 1) * dy
                    write(10, '(3F12.6)') x, y, orbital_data(i, j)
                end do
            end do
            close(10)
            print *, "Saved data for l =", l, "m =", m
        end do
    end do

    ! Deallocate the array
    deallocate(orbital_data)

    print *, "All orbital data saved."

end program hydrogen_orbital_animation

