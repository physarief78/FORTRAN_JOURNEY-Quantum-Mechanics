program schrodinger_gaussian_2d
    use omp_lib
    implicit none

    ! Define a kind for double precision
    integer, parameter :: dp = selected_real_kind(15, 307)

    ! Grid/time parameters
    integer, parameter :: nx = 200, ny = 200
    integer, parameter :: nt = 2000, nsteps = 2
    real(dp), parameter :: Lx = 2.5_dp, Ly = 2.5_dp
    real(dp), parameter :: dx = Lx/nx, dy = Ly/ny
    real(dp), parameter :: dt = 1.0e-5_dp, pi = 3.141592653589793_dp
    real(dp), parameter :: hbar = 1.0_dp, m = 1.0_dp

    ! Gaussian‐pillar lattice parameters
    real(dp), parameter :: V0         = 5000.0_dp
    real(dp), parameter :: a_latt     = 0.15_dp
    real(dp), parameter :: x_lat_max  = 0.75_dp
    real(dp), parameter :: sigma_latt = 0.03_dp

    ! Arrays
    complex(dp), dimension(nx, ny) :: psi, psi_temp, psi_new
    real(dp),    dimension(nx, ny) :: x, y, V, density

    ! Loop indices and other locals
    integer :: i, j, k, save_count
    integer :: nx_cells, ny_cells, ix, iy
    real(dp) :: x0, y0, sigma, kx, ky, norm
    real(dp) :: cx, cy, dx2, dy2
    complex(dp) :: imag_unit
    character(len=100) :: filename, folder_name
    logical :: dir_exists

    ! Create output folder if missing
    folder_name = 'output_data'
    inquire(file=trim(folder_name)//'/.', exist=dir_exists)
    if (.not. dir_exists) call system('mkdir '//trim(folder_name))

    ! Initial wavepacket
    imag_unit = (0.0_dp, 1.0_dp)
    x0 = -Lx/5.0_dp;  y0 = -Ly/5.0_dp
    sigma = 0.2_dp
    kx = 30.0_dp*pi; ky = 30.0_dp

    ! Precompute how many pillars fit in each direction
    nx_cells = int(x_lat_max / a_latt)
    ny_cells = int((Ly/2.0_dp) / a_latt)

    ! Build grid and Gaussian‐pillar potential
    !$omp parallel do private(i,j,ix,iy,cx,cy,dx2,dy2) shared(x,y,V)
    do i = 1, nx
        do j = 1, ny
            x(i,j) = -Lx/2.0_dp + (i-1)*dx
            y(i,j) = -Ly/2.0_dp + (j-1)*dy
            V(i,j) = 0.0_dp

            if (x(i,j) >= 0.0_dp .and. x(i,j) <= x_lat_max) then
                do ix = 0, nx_cells
                    cx = ix * a_latt
                    do iy = -ny_cells, ny_cells
                        cy = iy * a_latt
                        dx2 = (x(i,j)-cx)**2
                        dy2 = (y(i,j)-cy)**2
                        V(i,j) = V(i,j) + V0 * exp(-(dx2 + dy2)/(2.0_dp*sigma_latt**2))
                    end do
                end do
            end if
        end do
    end do
    !$omp end parallel do

    ! Initialize wavefunction
    !$omp parallel do private(i,j) shared(psi)
    do i = 1, nx
        do j = 1, ny
            psi(i,j) = exp(-0.5_dp*(((x(i,j)-x0)/sigma)**2 + ((y(i,j)-y0)/sigma)**2)) * &
                       exp(imag_unit*(kx*x(i,j) + ky*y(i,j)))
        end do
    end do
    !$omp end parallel do

    ! Normalize
    norm = 0.0_dp
    !$omp parallel do private(i,j) reduction(+:norm) shared(psi)
    do i = 1, nx
        do j = 1, ny
            norm = norm + abs(psi(i,j))**2 * dx * dy
        end do
    end do
    !$omp end parallel do
    psi = psi / sqrt(norm)

    ! Open main output file
    open(unit=10, file=trim(folder_name)//'/wavepacket_animation_2d.dat', status="unknown")
    save_count = 0

    ! Time evolution
    do k = 1, nt
        call crank_nicolson_x(psi, psi_temp, V, nx, ny, dx, dt, hbar, m)
        call crank_nicolson_y(psi_temp, psi_new, V, nx, ny, dy, dt, hbar, m)
        psi = psi_new

        !$omp parallel do private(i,j) shared(density,psi)
        do i = 1, nx
            do j = 1, ny
                density(i,j) = abs(psi(i,j))**2
            end do
        end do
        !$omp end parallel do

        if (mod(k, nsteps) == 0) then
            save_count = save_count + 1
            write(10,'(F10.5)',advance='no') k*dt
            do i = 1, nx
                do j = 1, ny
                    write(10,'(F10.5)',advance='no') density(i,j)
                end do
            end do
            write(10,*)

            write(filename,'(A,I0.5,A)') trim(folder_name)//'/wavepacket_', save_count, '.dat'
            open(unit=20, file=filename, status='unknown')
            write(20,'(A,F10.5)') '# Time step: ', k*dt
            do i = 1, nx
                do j = 1, ny
                    write(20,'(3F15.8)') x(i,j), y(i,j), density(i,j)
                end do
                write(20,*)
            end do
            close(20)
        end if
    end do

    close(10)

contains

    subroutine crank_nicolson_x(psi_in, psi_out, V, nx, ny, dx, dt, hbar, m)
        implicit none
        complex(dp), intent(in)  :: psi_in(nx,ny)
        complex(dp), intent(out) :: psi_out(nx,ny)
        real(dp),    intent(in)  :: V(nx,ny), dx, dt, hbar, m
        integer,     intent(in)  :: nx, ny
        complex(dp) :: rhs(nx), diag(nx), lower(nx), upper(nx), sol(nx)
        complex(dp) :: imag_unit, factor
        integer :: i,j

        imag_unit = (0.0_dp,1.0_dp)
        factor = imag_unit*hbar*dt/(2.0_dp*m*dx**2)

        !$omp parallel do private(i,j,rhs,diag,lower,upper,sol)
        do j = 1, ny
            diag  = 1.0_dp + 2.0_dp*factor
            lower = -factor
            upper = -factor
            do i = 1, nx
                diag(i) = diag(i) + imag_unit*dt/(2.0_dp*hbar)*V(i,j)
            end do
            rhs(1) = psi_in(1,j)
            rhs(nx)= psi_in(nx,j)
            do i = 2, nx-1
                rhs(i) = (1.0_dp-2.0_dp*factor)*psi_in(i,j) + &
                         factor*(psi_in(i-1,j)+psi_in(i+1,j)) - &
                         imag_unit*dt/(2.0_dp*hbar)*V(i,j)*psi_in(i,j)
            end do
            call tridiagonal_solver(lower,diag,upper,rhs,sol,nx)
            psi_out(:,j) = sol
        end do
        !$omp end parallel do
    end subroutine

    subroutine crank_nicolson_y(psi_in, psi_out, V, nx, ny, dy, dt, hbar, m)
        implicit none
        complex(dp), intent(in)  :: psi_in(nx,ny)
        complex(dp), intent(out) :: psi_out(nx,ny)
        real(dp),    intent(in)  :: V(nx,ny), dy, dt, hbar, m
        integer,     intent(in)  :: nx, ny
        complex(dp) :: rhs(ny), diag(ny), lower(ny), upper(ny), sol(ny)
        complex(dp) :: imag_unit, factor
        integer :: i,j

        imag_unit = (0.0_dp,1.0_dp)
        factor = imag_unit*hbar*dt/(2.0_dp*m*dy**2)

        !$omp parallel do private(i,j,rhs,diag,lower,upper,sol)
        do i = 1, nx
            diag  = 1.0_dp + 2.0_dp*factor
            lower = -factor
            upper = -factor
            do j = 1, ny
                diag(j) = diag(j) + imag_unit*dt/(2.0_dp*hbar)*V(i,j)
            end do
            rhs(1) = psi_in(i,1)
            rhs(ny)= psi_in(i,ny)
            do j = 2, ny-1
                rhs(j) = (1.0_dp-2.0_dp*factor)*psi_in(i,j) + &
                         factor*(psi_in(i,j-1)+psi_in(i,j+1)) - &
                         imag_unit*dt/(2.0_dp*hbar)*V(i,j)*psi_in(i,j)
            end do
            call tridiagonal_solver(lower,diag,upper,rhs,sol,ny)
            psi_out(i,:) = sol
        end do
        !$omp end parallel do
    end subroutine

    subroutine tridiagonal_solver(lower, diag, upper, rhs, sol, n)
        implicit none
        integer, intent(in)           :: n
        complex(dp), intent(in)       :: lower(n), diag(n), upper(n), rhs(n)
        complex(dp), intent(out)      :: sol(n)
        complex(dp) :: tmp(n)
        integer :: i

        tmp(1)   = diag(1)
        sol(1)   = rhs(1)
        do i = 2, n
            tmp(i) = diag(i) - lower(i)*upper(i-1)/tmp(i-1)
            sol(i) = rhs(i)  - lower(i)*sol(i-1)/tmp(i-1)
        end do

        sol(n) = sol(n)/tmp(n)
        do i = n-1,1,-1
            sol(i) = (sol(i) - upper(i)*sol(i+1)) / tmp(i)
        end do
    end subroutine

end program schrodinger_gaussian_2d
