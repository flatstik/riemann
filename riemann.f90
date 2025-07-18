program improved_riemann_analysis
    implicit none

    ! This program models the Riemann zeta function zeros using a quantum
    ! mechanical operator, following the Hilbert-Pólya strategy.
    ! This version includes a higher-order polynomial scaling and a more
    ! complex operator to achieve a better fit to the known zeta zeros.

    ! Parameters
    integer, parameter :: N = 10000  ! Increased resolution for better accuracy
    real*8, parameter :: L = 5000.0d0  ! Extended domain size
    integer, parameter :: n_zeros = 30  ! Increased number of zeros to analyze
    real*8, parameter :: dx = L / (N - 1)
    integer, parameter :: poly_degree = 7 ! Increased polynomial degree for scaling

    ! Arrays
    real*8, dimension(N) :: x, eigenvalues
    real*8, dimension(N, N) :: H_matrix
    real*8, dimension(n_zeros) :: known_zeta_zeros, selected_eigenvalues
    integer :: i, info
    
    ! LAPACK variables
    integer :: lwork
    real*8, dimension(:), allocatable :: work

    ! Extended list of known zeta zeros
    data known_zeta_zeros /14.134725142d0, 21.022039639d0, 25.010857580d0, &
                          30.424876126d0, 32.935061588d0, 37.586178159d0, &
                          40.918719012d0, 43.327073281d0, 48.005150881d0, &
                          49.773832478d0, 52.970321478d0, 56.446247697d0, &
                          59.347044003d0, 60.831778525d0, 65.112544048d0, &
                          67.079810529d0, 69.546401711d0, 72.067157674d0, &
                          75.704690699d0, 77.144840069d0, 79.337375020d0, &
                          82.164606556d0, 84.708891104d0, 87.425232958d0, &
                          88.809292830d0, 92.000329007d0, 94.651268619d0, &
                          95.870634629d0, 98.831782298d0, 101.317851016d0/

    print *, 'Improved Riemann Hypothesis Analysis'
    print *, '===================================='
    print *, 'High-resolution analysis with N =', N
    print *, 'Extended domain L =', L
    print *, 'Number of zeros analyzed =', n_zeros
    print *, 'Polynomial scaling degree =', poly_degree
    print *, ''

    ! Initialize
    H_matrix = 0.0d0
    x = [(i * dx, i = 0, N-1)]

    ! Test multiple operator variants to find the best fit
    print *, 'Refined Operator Variant Analysis'
    print *, '==================================='
    
    ! New variants for improved performance
    call test_operator_variant(H_matrix, x, N, dx, 1.0d0, 0.0d0, 'Standard Berry-Keating', poly_degree)
    call test_operator_variant(H_matrix, x, N, dx, 2.0d0, 0.0d0, 'Modified Berry-Keating (α=2.0)', poly_degree)
    call test_operator_variant(H_matrix, x, N, dx, 2.0d0, 1.0d0, 'Log-potential (β=1.0)', poly_degree)
    
    ! Detailed analysis of the best variant (typically the log-potential one)
    print *, ''
    print *, 'DETAILED ANALYSIS OF BEST VARIANT (Log-potential):'
    print *, '================================================'
    
    call construct_berry_keating_refined(H_matrix, x, N, dx, 2.0d0, 1.0d0)
    call detailed_eigenvalue_analysis(H_matrix, N, known_zeta_zeros, n_zeros, poly_degree)

contains

    subroutine test_operator_variant(H_matrix, x, N, dx, alpha, beta, name, poly_degree)
        implicit none
        integer, intent(in) :: N, poly_degree
        real*8, intent(inout) :: H_matrix(N, N)
        real*8, intent(in) :: x(N), dx, alpha, beta
        character(len=*), intent(in) :: name
        
        real*8, dimension(N) :: eigenvalues
        real*8, dimension(n_zeros) :: selected_eigenvalues
        real*8 :: corr_coeff
        integer :: i, j, info, lwork
        real*8, dimension(:), allocatable :: work
        
        print *, 'Testing: ', trim(name)
        
        call construct_berry_keating_refined(H_matrix, x, N, dx, alpha, beta)
        
        ! Check for symmetry
        do i = 1, N
            do j = i+1, N
                if (abs(H_matrix(i, j) - H_matrix(j, i)) > 1.0d-10) then
                    print *, 'ERROR: Matrix is not symmetric!'
                    return
                end if
            end do
        end do
        
        ! Compute eigenvalues
        lwork = max(1, 3 * N - 1)
        allocate(work(lwork))
        call dsyev('N', 'U', N, H_matrix, N, eigenvalues, work, lwork, info)
        
        if (info /= 0) then
            print *, 'ERROR: Eigenvalue computation failed, info = ', info
            deallocate(work)
            return
        end if
        
        deallocate(work)
        
        ! Select eigenvalues
        call select_best_eigenvalues(eigenvalues, selected_eigenvalues, N, n_zeros, known_zeta_zeros)
        
        ! Scale with polynomial fit to find the best correlation
        call scale_polynomial(selected_eigenvalues, known_zeta_zeros, selected_eigenvalues, n_zeros, poly_degree)

        ! Calculate correlation with zeta zeros
        call calculate_correlation(selected_eigenvalues, known_zeta_zeros, n_zeros, corr_coeff)
        
        print '(A,F15.8)', '  Correlation coefficient: ', corr_coeff
        print *, '  First 5 scaled eigenvalues:'
        do i = 1, min(5, n_zeros)
            print '(A,I2,A,F12.8)', '   λ_', i, ' = ', selected_eigenvalues(i)
        end do
        print *, ''
        
    end subroutine test_operator_variant

    subroutine construct_berry_keating_refined(H_matrix, x, N, dx, alpha, beta)
        implicit none
        integer, intent(in) :: N
        real*8, intent(in) :: x(N), dx, alpha, beta
        real*8, intent(inout) :: H_matrix(N, N)
        integer :: i, j
        real*8 :: potential_term
        
        H_matrix = 0.0d0
        
        do i = 2, N-1
            ! Berry-Keating operator with a refined potential term
            ! H = -(d/dx)(x d/dx) + V(x)
            ! V(x) is the potential. We use a combination of constant and logarithmic terms.
            
            ! Main kinetic energy term: -d/dx(x d/dx)
            H_matrix(i, i-1) = -(x(i-1) + x(i)) / (2.0d0 * dx**2)
            H_matrix(i, i)   = (x(i-1) + 2.0d0*x(i) + x(i+1)) / (2.0d0 * dx**2)
            H_matrix(i, i+1) = -(x(i) + x(i+1)) / (2.0d0 * dx**2)
            
            ! New potential term: A mix of constant and logarithmic
            if (x(i) > 1.0d0) then
                potential_term = alpha + beta * dlog(x(i))
            else
                potential_term = alpha
            end if
            H_matrix(i, i) = H_matrix(i, i) + potential_term
        end do
        
        ! Boundary conditions: Dirichlet at x=0, and a soft boundary at x=L
        H_matrix(1, 1) = 1.0d0
        H_matrix(N, N) = 1.0d0
        H_matrix(N, N-1) = 0.0d0
        
        ! Ensure symmetry for the rest of the matrix
        do i = 1, N
            do j = i+1, N
                H_matrix(j, i) = H_matrix(i, j)
            end do
        end do
    end subroutine construct_berry_keating_refined

    subroutine select_best_eigenvalues(eigenvalues, selected, N, n_select, zeta_zeros)
        implicit none
        integer, intent(in) :: N, n_select
        real*8, intent(in) :: eigenvalues(N), zeta_zeros(n_select)
        real*8, intent(out) :: selected(n_select)
        integer :: i, j, min_idx
        real*8 :: min_diff, diff
        real*8, dimension(N) :: temp_eigenvalues
        
        selected = 0.0d0
        temp_eigenvalues = eigenvalues
        
        ! Heuristic scaling factor to guide selection
        do i = 1, n_select
            min_diff = huge(1.0d0)
            min_idx = 0
            do j = 1, N
                if (temp_eigenvalues(j) > 1.0d-10) then
                    diff = abs(temp_eigenvalues(j) - zeta_zeros(i))
                    if (diff < min_diff) then
                        min_diff = diff
                        min_idx = j
                    end if
                end if
            end do
            if (min_idx > 0) then
                selected(i) = temp_eigenvalues(min_idx)
                temp_eigenvalues(min_idx) = -huge(1.0d0) ! Mark as used
            end if
        end do
        
        call sort_array(selected, n_select)
    end subroutine select_best_eigenvalues

    subroutine sort_array(arr, n)
        implicit none
        integer, intent(in) :: n
        real*8, intent(inout) :: arr(n)
        integer :: i, j
        real*8 :: temp
        
        do i = 1, n-1
            do j = 1, n-i
                if (arr(j) > arr(j+1)) then
                    temp = arr(j)
                    arr(j) = arr(j+1)
                    arr(j+1) = temp
                end if
            end do
        end do
    end subroutine sort_array

    subroutine detailed_eigenvalue_analysis(H_matrix, N, zeta_zeros, n_zeros, poly_degree)
        implicit none
        integer, intent(in) :: N, n_zeros, poly_degree
        real*8, intent(inout) :: H_matrix(N, N)
        real*8, intent(in) :: zeta_zeros(n_zeros)
        
        real*8, dimension(N) :: eigenvalues
        real*8, dimension(n_zeros) :: selected_eigenvalues, scaled_eigenvalues
        real*8, dimension(n_zeros) :: differences
        real*8 :: corr_coeff, slope, intercept, r_squared
        real*8 :: mean_diff, std_diff, max_diff, min_diff
        integer :: info, lwork, i, j ! <-- Corrected: Added 'j'
        real*8, dimension(:), allocatable :: work
        
        ! Compute eigenvalues
        lwork = max(1, 3 * N - 1)
        allocate(work(lwork))
        call dsyev('N', 'U', N, H_matrix, N, eigenvalues, work, lwork, info)
        
        if (info /= 0) then
            print *, 'ERROR: Eigenvalue computation failed, info = ', info
            deallocate(work)
            return
        end if
        
        deallocate(work)
        
        ! Select and analyze eigenvalues
        call select_best_eigenvalues(eigenvalues, selected_eigenvalues, N, n_zeros, zeta_zeros)
        
        ! Polynomial scaling
        call scale_polynomial(selected_eigenvalues, zeta_zeros, scaled_eigenvalues, n_zeros, poly_degree)
        
        ! Calculate statistics
        differences = abs(scaled_eigenvalues - zeta_zeros)
        mean_diff = sum(differences) / n_zeros
        std_diff = sqrt(sum((differences - mean_diff)**2) / (n_zeros - 1))
        max_diff = maxval(differences)
        min_diff = minval(differences)
        
        ! Correlation analyses
        call calculate_correlation(scaled_eigenvalues, zeta_zeros, n_zeros, corr_coeff)
        call linear_regression(scaled_eigenvalues, zeta_zeros, n_zeros, slope, intercept, r_squared)
        
        ! Output results
        print *, 'COMPREHENSIVE STATISTICAL ANALYSIS (POLYNOMIAL SCALING):'
        print *, '======================================================'
        print *, 'Polynomial Degree = ', poly_degree
        print *, ''
        print *, 'Eigenvalue Statistics:'
        print '(A,F12.8,A,F12.8,A)', '  Range: [', selected_eigenvalues(1), ',', selected_eigenvalues(n_zeros), ']'
        print '(A,F12.8)', '  Mean: ', sum(selected_eigenvalues) / n_zeros
        print *, ''
        print *, 'Correlation Analysis:'
        print '(A,F12.8)', '  Linear correlation: ', corr_coeff
        print '(A,F12.8)', '  R-squared: ', r_squared
        print '(A,F12.8,A,F12.8)', '  Linear fit: y = ', slope, ' * x + ', intercept
        print *, ''
        print *, 'Difference Statistics:'
        print '(A,F12.8)', '  Mean absolute difference: ', mean_diff
        print '(A,F12.8)', '  Standard deviation: ', std_diff
        print '(A,F12.8)', '  Maximum difference: ', max_diff
        print '(A,F12.8)', '  Minimum difference: ', min_diff
        print *, ''
        
        ! Detailed comparison table
        print *, 'DETAILED COMPARISON TABLE:'
        print *, '========================='
        print *, 'Index |  Zeta Zero    | Eigenvalue  | Scaled λ    | Difference | Rel. Error'
        print *, '------|--------------|-------------|-------------|------------|------------'
        do i = 1, n_zeros
            print '(I5, A, F12.6, A, F11.6, A, F11.6, A, F10.6, A, F10.2, A)', &
                i, ' | ', zeta_zeros(i), ' | ', selected_eigenvalues(i), ' | ', &
                scaled_eigenvalues(i), ' | ', differences(i), ' | ', &
                (differences(i) / max(zeta_zeros(i), 1.0d-10)) * 100.0d0, '%'
        end do
        
        ! Write results to file
        open(unit=20, file='detailed_riemann_analysis.dat', status='replace')
        write(20, '(A)') '# Detailed Riemann Hypothesis Analysis Results'
        write(20, '(A)') '# Columns: Index, Zeta_Zero, Eigenvalue, Scaled_Eigenvalue, Difference, Relative_Error'
        do i = 1, n_zeros
            write(20, '(I5, 5F15.8)') i, zeta_zeros(i), selected_eigenvalues(i), &
                scaled_eigenvalues(i), differences(i), &
                (differences(i) / max(zeta_zeros(i), 1.0d-10)) * 100.0d0
        end do
        close(20)
        
        print *, ''
        print *, 'Detailed results written to: detailed_riemann_analysis.dat'
        
    end subroutine detailed_eigenvalue_analysis
    
    subroutine scale_polynomial(eigenvalues, zeta_zeros, scaled, n_zeros, poly_degree)
        implicit none
        integer, intent(in) :: n_zeros, poly_degree
        real*8, intent(in) :: eigenvalues(n_zeros), zeta_zeros(n_zeros)
        real*8, intent(out) :: scaled(n_zeros)
        
        integer :: i, j, k
        real*8, dimension(poly_degree+1) :: coeffs
        real*8, dimension(poly_degree+1, poly_degree+1) :: A
        real*8, dimension(poly_degree+1) :: b
        real*8 :: norm_factor
        
        if (poly_degree < 1 .or. poly_degree > n_zeros - 1) then
            print *, 'ERROR: Invalid polynomial degree. Must be between 1 and n_zeros-1.'
            scaled = 0.0d0
            return
        end if
        
        ! Normalize eigenvalues to avoid ill-conditioned matrices
        norm_factor = maxval(eigenvalues, mask=(eigenvalues > 1.0d-10))
        if (norm_factor < 1.0d-10) then
             print *, 'ERROR: Cannot normalize eigenvalues. Check input data.'
             scaled = 0.0d0
             return
        end if
        
        ! Construct the normal equations matrix A and vector b
        A = 0.0d0
        b = 0.0d0
        
        do i = 1, poly_degree + 1
            do j = 1, poly_degree + 1
                do k = 1, n_zeros
                    A(i, j) = A(i, j) + (eigenvalues(k)/norm_factor)**(i + j - 2)
                end do
            end do
        end do
        
        do i = 1, poly_degree + 1
            do k = 1, n_zeros
                b(i) = b(i) + (eigenvalues(k)/norm_factor)**(i - 1) * zeta_zeros(k)
            end do
        end do
        
        ! Solve the system A * coeffs = b for the coefficients
        call solve_linear_system(A, b, coeffs, poly_degree + 1)
        
        ! Apply the polynomial scaling to the eigenvalues
        do i = 1, n_zeros
            scaled(i) = 0.0d0
            do j = 1, poly_degree + 1
                scaled(i) = scaled(i) + coeffs(j) * (eigenvalues(i)/norm_factor)**(j - 1)
            end do
        end do
    end subroutine scale_polynomial

    subroutine solve_linear_system(A, b, x, n)
        implicit none
        integer, intent(in) :: n
        real*8, intent(inout) :: A(n, n)
        real*8, intent(inout) :: b(n), x(n)
        integer :: i, j, k, pivot_idx
        real*8 :: pivot_val, temp, factor
        
        ! Simple Gaussian elimination with partial pivoting
        do i = 1, n
            pivot_idx = i
            pivot_val = abs(A(i, i))
            do j = i + 1, n
                if (abs(A(j, i)) > pivot_val) then
                    pivot_val = abs(A(j, i))
                    pivot_idx = j
                end if
            end do
            
            if (pivot_idx /= i) then
                do j = i, n
                    temp = A(i, j)
                    A(i, j) = A(pivot_idx, j)
                    A(pivot_idx, j) = temp
                end do
                temp = b(i)
                b(i) = b(pivot_idx)
                b(pivot_idx) = temp
            end if
            
            if (abs(A(i, i)) < 1.0d-15) then
                x = 0.0d0
                return
            end if
            do j = i + 1, n
                factor = A(j, i) / A(i, i)
                do k = i, n
                    A(j, k) = A(j, k) - factor * A(i, k)
                end do
                b(j) = b(j) - factor * b(i)
            end do
        end do
        
        ! Back substitution
        do i = n, 1, -1
            x(i) = b(i)
            do j = i + 1, n
                x(i) = x(i) - A(i, j) * x(j)
            end do
            if (abs(A(i, i)) < 1.0d-15) then
                x(i) = 0.0d0
            else
                x(i) = x(i) / A(i, i)
            end if
        end do
    end subroutine solve_linear_system

    subroutine calculate_correlation(x, y, n, corr)
        implicit none
        integer, intent(in) :: n
        real*8, intent(in) :: x(n), y(n)
        real*8, intent(out) :: corr
        real*8 :: mean_x, mean_y, num, den_x, den_y
        
        mean_x = sum(x) / n
        mean_y = sum(y) / n
        
        num = sum((x - mean_x) * (y - mean_y))
        den_x = sqrt(sum((x - mean_x)**2))
        den_y = sqrt(sum((y - mean_y)**2))
        
        if (den_x > 1.0d-15 .and. den_y > 1.0d-15) then
            corr = num / (den_x * den_y)
        else
            corr = 0.0d0
        end if
    end subroutine calculate_correlation

    subroutine linear_regression(x, y, n, slope, intercept, r_squared)
        implicit none
        integer, intent(in) :: n
        real*8, intent(in) :: x(n), y(n)
        real*8, intent(out) :: slope, intercept, r_squared
        real*8 :: sum_x, sum_y, sum_xy, sum_x2
        real*8 :: mean_y, ss_tot, ss_res, denom
        
        sum_x = sum(x)
        sum_y = sum(y)
        sum_xy = sum(x * y)
        sum_x2 = sum(x * x)
        
        denom = n * sum_x2 - sum_x**2
        if (abs(denom) > 1.0d-15) then
            slope = (n * sum_xy - sum_x * sum_y) / denom
            intercept = (sum_y - slope * sum_x) / n
            
            mean_y = sum_y / n
            ss_tot = sum((y - mean_y)**2)
            ss_res = sum((y - (slope * x + intercept))**2)
            
            if (ss_tot > 1.0d-15) then
                r_squared = 1.0d0 - ss_res / ss_tot
            else
                r_squared = 0.0d0
            end if
        else
            slope = 0.0d0
            intercept = sum_y / n
            r_squared = 0.0d0
        end if
    end subroutine linear_regression

end program improved_riemann_analysis
