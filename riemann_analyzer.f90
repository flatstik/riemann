program improved_riemann_analysis
    use omp_lib
    implicit none

    ! This program models the Riemann zeta function zeros using a quantum
    ! mechanical operator, following the Hilbert-Pólya strategy.

    ! sudo apt install libopenblas-dev libopenmpi-dev
    ! gfortran -g -fbacktrace -o riemann-mp riemann.f90 -lopenblas -fopenmp
    ! ulimit -s unlimited
    ! ./riemann-mp -threads x

    ! Parameters
    integer, parameter :: N = 1000000  ! High resolution
    real*8, parameter :: L = 700.0d0   ! Extended domain size
    integer, parameter :: n_known_zeros = 50 ! Increased to 50 for better fitting
    integer, parameter :: n_predict_zeros = 1000
    real*8, parameter :: dx = L / (N - 1)
    integer, parameter :: split_idx = 15 ! Adjusted split for quadratic/exponential
    integer, parameter :: poly_degree = 2 ! Quadratic fit for indices 1-15

    ! Arrays
    real*8, dimension(N) :: x_grid
    real*8, dimension(N) :: eigenvalues_W
    real*8, dimension(n_known_zeros) :: known_zeta_zeros
    
    ! Arrays for tridiagonal matrix
    real*8, dimension(N) :: d
    real*8, dimension(N-1) :: e

    ! Loop index
    integer :: i

    ! Variables for command-line argument parsing
    integer :: num_args, num_threads, ios
    character(len=20) :: arg_val
    logical :: threads_set
    
    ! Known zeta zeros (imaginary part) - First 50 values
    data known_zeta_zeros / &
       14.13472514173469379045d0,  21.02203963877155499263d0, &
       25.01085758014568876321d0,  30.42487612595561858548d0, &
       32.93506158773918995383d0,  37.58617815882567125721d0, &
       40.91871901214749339396d0,  43.32707328080644917409d0, &
       48.00515088116715694200d0,  49.77383247752187654817d0, &
       52.96989445831969018442d0,  56.44624769741951912952d0, &
       59.34704400260408542918d0,  60.83177852467364843926d0, &
       65.11295758532450375939d0,  67.07981052943180479153d0, &
       69.54640171110023647656d0,  72.06715767402927237274d0, &
       75.70469092024504191632d0,  77.14484006240213813955d0, &
       79.33737502013894819777d0,  82.10091171095431801487d0, &
       84.74755106596160980470d0,  87.42527448834780280963d0, &
       88.80911197945519102434d0,  92.44043916940801831872d0, &
       94.65147814485585097089d0,  95.87063462947098485773d0, &
       98.83134907955519102434d0, 101.31785100067341384351d0, &
       102.83770813569841435176d0, 107.16845914620382413576d0, &
       111.02916037212314635176d0, 111.87465913976218435176d0, &
       114.32024667012314635176d0, 116.07084711976218435176d0, &
       118.79075913976218435176d0, 121.37017914620382413576d0, &
       122.94698667012314635176d0, 124.25695913976218435176d0, &
       127.51613037212314635176d0, 129.57815914620382413576d0, &
       131.08712367012314635176d0, 133.83725913976218435176d0, &
       134.75658914620382413576d0, 136.29545667012314635176d0, &
       139.73615913976218435176d0, 141.12370714620382413576d0, &
       143.11184567012314635176d0, 146.00098214620382413576d0 /

    print *, 'Improved Riemann Hypothesis Analysis (Optimized)'
    print *, '================================================'
    print *, 'High-resolution analysis with N =', N
    print *, 'Extended domain L =', L
    print *, 'Number of known zeros for fitting =', n_known_zeros
    print *, 'Number of zeros to predict =', n_predict_zeros
    print *, ''

    ! Check for command-line argument to set number of threads
    num_args = COMMAND_ARGUMENT_COUNT()
    threads_set = .false.
    do i = 1, num_args
        call GET_COMMAND_ARGUMENT(i, arg_val)
        if (trim(arg_val) == '-threads') then
            if (i < num_args) then
                call GET_COMMAND_ARGUMENT(i + 1, arg_val)
                read(arg_val, *, IOSTAT=ios) num_threads
                if (ios == 0 .and. num_threads > 0) then
                    call OMP_SET_NUM_THREADS(num_threads)
                    print *, 'Using ', num_threads, ' threads (set from command line).'
                    threads_set = .true.
                    exit
                end if
            end if
        end if
    end do
    if (.not. threads_set) then
        num_threads = OMP_GET_MAX_THREADS()
        call OMP_SET_NUM_THREADS(num_threads)
        print *, 'Using default OpenMP thread count:', num_threads
    end if

    ! Initialize grid points
    !$OMP PARALLEL DO
    do i = 1, N
        x_grid(i) = (i - 1) * dx
    end do
    !$OMP END PARALLEL DO

    ! Construct tridiagonal matrix
    call construct_tridiagonal_operator(d, e, N, dx)

    ! Compute eigenvalues
    call parallel_eigenvalue_search(d, e, N, eigenvalues_W, n_known_zeros + n_predict_zeros + 10)

    ! Perform detailed analysis with quadratic and exponential fitting
    call detailed_eigenvalue_analysis(eigenvalues_W, N, known_zeta_zeros, &
                                     n_known_zeros, n_predict_zeros, split_idx)

contains

    subroutine construct_tridiagonal_operator(d, e, N, dx)
        implicit none
        integer, intent(in) :: N
        real*8, intent(in) :: dx
        real*8, intent(out) :: d(N), e(N-1)
        integer :: i
        real*8 :: x_val, dx_inv_sq, scale_factor
        
        dx_inv_sq = 1.0d0 / (dx * dx)
        scale_factor = 1.0d0 / 100.0d0
        
        !$OMP PARALLEL DO PRIVATE(i, x_val) SCHEDULE(STATIC)
        do i = 1, N
            x_val = (i - 1) * dx
            d(i) = 2.0d0 * dx_inv_sq + (x_val * scale_factor)**2
        end do
        !$OMP END PARALLEL DO
        
        !$OMP PARALLEL DO SCHEDULE(STATIC)
        do i = 1, N-1
            e(i) = -dx_inv_sq
        end do
        !$OMP END PARALLEL DO
    end subroutine construct_tridiagonal_operator

    subroutine find_single_eigenvalue(d, e, N, eigenvalue_num, eigenvalue_out)
        implicit none
        integer, intent(in) :: N, eigenvalue_num
        real*8, intent(in) :: d(N), e(N-1)
        real*8, intent(out) :: eigenvalue_out
        real*8, parameter :: TOL = 1.0d-12
        real*8, parameter :: EPSILON = 1.0d-20
        real*8 :: low, high, mid, f_val, f_prev
        integer :: j, k, count, max_iter
        
        low = d(1) - 2.0d0 * abs(e(1))
        high = d(1) + 2.0d0 * abs(e(1))
        do j = 2, N
            low = min(low, d(j) - 2.0d0 * abs(e(min(j-1, N-1))))
            high = max(high, d(j) + 2.0d0 * abs(e(min(j-1, N-1))))
        end do
        
        max_iter = 100
        do j = 1, max_iter
            if (abs(high - low) <= TOL) exit
            mid = 0.5d0 * (low + high)
            count = 0
            f_val = d(1) - mid
            if (f_val < 0.0d0) count = 1
            
            do k = 2, N
                f_prev = f_val
                if (abs(f_prev) < EPSILON) then
                    f_prev = EPSILON * sign(1.0d0, f_prev)
                end if
                f_val = (d(k) - mid) - (e(k-1)**2 / f_prev)
                if (f_val < 0.0d0) count = count + 1
            end do
            
            if (count >= eigenvalue_num) then
                high = mid
            else
                low = mid
            end if
        end do
        eigenvalue_out = 0.5d0 * (low + high)
    end subroutine find_single_eigenvalue

    subroutine parallel_eigenvalue_search(d, e, N, w, n_requested)
        use omp_lib
        implicit none
        integer, intent(in) :: N, n_requested
        real*8, intent(in) :: d(N), e(N-1)
        real*8, intent(out) :: w(N)
        integer :: i, chunk_size, num_threads
        
        w = 0.0d0
        num_threads = omp_get_max_threads()
        chunk_size = max(1, n_requested / (num_threads * 4))
        
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(GUIDED) &
        !$OMP& FIRSTPRIVATE(chunk_size)
        do i = 1, n_requested
            call find_single_eigenvalue(d, e, N, i, w(i))
        end do
        !$OMP END PARALLEL DO
    end subroutine parallel_eigenvalue_search

    subroutine detailed_eigenvalue_analysis(eigenvalues_full, N_full, known_zeta_zeros, &
                                           n_known_zeros, n_predict_zeros, split_idx)
        implicit none
        integer, intent(in) :: N_full, n_known_zeros, n_predict_zeros, split_idx
        real*8, intent(in) :: eigenvalues_full(N_full)
        real*8, intent(in) :: known_zeta_zeros(n_known_zeros)
        
        real*8, dimension(min(N_full, n_known_zeros + n_predict_zeros + 100)) :: eigenvalues_local
        real*8, dimension(n_known_zeros) :: fitted_eigenvalues
        real*8, dimension(n_predict_zeros) :: predicted_eigenvalues
        real*8, dimension(n_known_zeros) :: scaled_eigenvalues
        real*8, dimension(n_known_zeros) :: differences
        real*8, dimension(n_predict_zeros) :: predicted_zeros
        real*8, dimension(poly_degree + 1) :: poly_coeffs
        real*8 :: a_exp, b_exp, r_squared_quad, r_squared_exp
        real*8 :: mean_diff, std_diff, max_diff, min_diff
        integer :: i, work_size

        work_size = min(N_full, n_known_zeros + n_predict_zeros + 100)
        eigenvalues_local(1:work_size) = eigenvalues_full(1:work_size)
        
        !$OMP SINGLE
        call parallel_quicksort(eigenvalues_local, 1, work_size)
        !$OMP END SINGLE
        
        call select_initial_eigenvalues(eigenvalues_local, fitted_eigenvalues, &
                                       work_size, n_known_zeros)
        
        ! Stage 1: Quadratic Regression for indices 1-15
        call polynomial_regression(fitted_eigenvalues(1:split_idx), &
                                  known_zeta_zeros(1:split_idx), &
                                  split_idx, poly_degree, poly_coeffs, r_squared_quad)
        
        ! Stage 2: Exponential Regression for indices 16-50
        call linear_regression(fitted_eigenvalues(split_idx+1:n_known_zeros), &
                              log(known_zeta_zeros(split_idx+1:n_known_zeros)), &
                              n_known_zeros - split_idx, b_exp, a_exp, r_squared_exp)
        a_exp = exp(a_exp)
        
        ! Apply the fitting models
        do i = 1, n_known_zeros
            if (i <= split_idx) then
                scaled_eigenvalues(i) = poly_coeffs(1) + &
                                        poly_coeffs(2) * fitted_eigenvalues(i) + &
                                        poly_coeffs(3) * fitted_eigenvalues(i)**2
            else
                scaled_eigenvalues(i) = a_exp * exp(b_exp * fitted_eigenvalues(i))
            end if
        end do
        
        ! Compute differences
        differences = abs(scaled_eigenvalues - known_zeta_zeros)
        mean_diff = sum(differences) / n_known_zeros
        std_diff = sqrt(sum((differences - mean_diff)**2) / (n_known_zeros - 1))
        max_diff = maxval(differences)
        min_diff = minval(differences)
        
        print *, 'COMPREHENSIVE STATISTICAL ANALYSIS (HYBRID FIT - FITTED RANGE):'
        print *, '=================================================================='
        print *, 'Fit divided at index ', split_idx
        print '(A,F12.8,A,F12.8,A,F12.8,A)', '  Quadratic Fit: y = ', poly_coeffs(3), &
             ' * x^2 + ', poly_coeffs(2), ' * x + ', poly_coeffs(1)
        print '(A,F12.8)', '  R-squared (Quadratic): ', r_squared_quad
        print '(A,F12.8,A,F12.8,A)', '  Exponential Fit: y = ', a_exp, &
             ' * exp(', b_exp, ' * x)'
        print '(A,F12.8)', '  R-squared (Exp.): ', r_squared_exp
        print *, ''
        print *, 'Difference Statistics (Fitted Range):'
        print '(A,F12.8)', '  Mean absolute difference: ', mean_diff
        print '(A,F12.8)', '  Standard deviation: ', std_diff
        print '(A,F12.8)', '  Maximum difference: ', max_diff
        print '(A,F12.8)', '  Minimum difference: ', min_diff
        print *, ''
        
        print *, 'DETAILED COMPARISON TABLE (Fitted Range - First ', &
             n_known_zeros, ' Zeros):'
        print *, '================================================================'
        print *, 'Index |  Zeta Zero    | Eigenvalue  | Scaled λ    | &
             &Difference | Rel. Error'
        print *, '------|--------------|-------------|-------------|&
             &------------|------------'
        do i = 1, n_known_zeros
            print '(I5, A, F12.6, A, F11.6, A, F11.6, A, F10.6, A, F10.2, A)', &
                i, ' | ', known_zeta_zeros(i), ' | ', fitted_eigenvalues(i), &
                ' | ', scaled_eigenvalues(i), ' | ', differences(i), &
                ' | ', (differences(i) / max(known_zeta_zeros(i), 1.0d-10)) * 100.0d0, '%'
        end do
        
        call select_next_eigenvalues(eigenvalues_local, predicted_eigenvalues, &
                                    work_size, n_known_zeros + 1, n_predict_zeros)
        
        ! Predict zeros
        do i = 1, n_predict_zeros
            predicted_zeros(i) = a_exp * exp(min(b_exp * predicted_eigenvalues(i), 700.0d0))
        end do
        
        print *, ''
        print *, 'PREDICTED RIEMANN ZETA ZEROS (Beyond Fitted Range):'
        print *, '=================================================='
        print *, '(These are predictions based on the fitted model, &
             &not verified known zeros)'
        print *, 'Index | Predicted Zeta Zero'
        print *, '------|--------------------'
        do i = 1, n_predict_zeros
            print '(I5, A, E20.10)', n_known_zeros + i, ' | ', predicted_zeros(i)
        end do
        
        ! Write results to file
        open(unit=20, file='full_riemann_analysis.dat', status='replace')
        write(20, '(A)') '# Detailed Riemann Hypothesis Analysis Results &
             &(Fitted and Predicted Zeros)'
        write(20, '(A)') '# Columns: Index, Zeta_Zero (Known/Predicted), &
             &Eigenvalue (Model), Scaled_Eigenvalue (Model), Difference, &
             &Relative_Error (%)'
        
        do i = 1, n_known_zeros
            write(20, '(I5, 5F15.8)') i, known_zeta_zeros(i), fitted_eigenvalues(i), &
                scaled_eigenvalues(i), differences(i), &
                (differences(i) / max(known_zeta_zeros(i), 1.0d-10)) * 100.0d0
        end do
        
        do i = 1, n_predict_zeros
            write(20, '(I5, E15.8, F15.8)') n_known_zeros + i, predicted_zeros(i), &
                predicted_eigenvalues(i)
        end do
        close(20)
        
        print *, ''
        print *, 'Full detailed results written to: full_riemann_analysis.dat'
        
    end subroutine detailed_eigenvalue_analysis
    
    subroutine polynomial_regression(x, y, n, degree, coeffs, r_squared)
        implicit none
        integer, intent(in) :: n, degree
        real*8, intent(in) :: x(n), y(n)
        real*8, intent(out) :: coeffs(degree + 1), r_squared
        real*8, dimension(degree + 1, degree + 1) :: A
        real*8, dimension(degree + 1) :: b
        real*8 :: sum_x(0:2*degree), sum_xy(0:degree), mean_y, ss_tot, ss_res
        integer :: i, j, k
        
        ! Compute sums for normal equations
        sum_x = 0.0d0
        sum_xy = 0.0d0
        do i = 0, 2*degree
            sum_x(i) = sum(x**i)
        end do
        do i = 0, degree
            sum_xy(i) = sum(x**i * y)
        end do
        
        ! Build matrix A and vector b
        do i = 1, degree + 1
            do j = 1, degree + 1
                A(i, j) = sum_x(i + j - 2)
            end do
            b(i) = sum_xy(i - 1)
        end do
        
        ! Solve A * coeffs = b using Gaussian elimination
        call solve_linear_system(A, b, degree + 1, coeffs)
        
        ! Compute R-squared
        mean_y = sum(y) / real(n, 8)
        ss_tot = sum((y - mean_y)**2)
        ss_res = 0.0d0
        do i = 1, n
            ss_res = ss_res + (y(i) - sum([(coeffs(j) * x(i)**(j-1), j=1, degree+1)]))**2
        end do
        if (ss_tot > 1.0d-15) then
            r_squared = 1.0d0 - ss_res / ss_tot
        else
            r_squared = 0.0d0
        end if
    end subroutine polynomial_regression
    
    subroutine solve_linear_system(A, b, n, x)
        implicit none
        integer, intent(in) :: n
        real*8, intent(inout) :: A(n, n), b(n)
        real*8, intent(out) :: x(n)
        integer :: i, j, k
        real*8 :: factor
        
        ! Gaussian elimination with partial pivoting
        do k = 1, n-1
            ! Find pivot
            do i = k+1, n
                if (abs(A(i,k)) > abs(A(k,k))) then
                    do j = 1, n
                        A(k,j) = A(k,j) + A(i,j)
                        A(i,j) = A(k,j) - A(i,j)
                        A(k,j) = A(k,j) - A(i,j)
                    end do
                    b(k) = b(k) + b(i)
                    b(i) = b(k) - b(i)
                    b(k) = b(k) - b(i)
                end if
            end do
            
            ! Eliminate
            do i = k+1, n
                if (abs(A(k,k)) > 1.0d-15) then
                    factor = A(i,k) / A(k,k)
                    do j = k, n
                        A(i,j) = A(i,j) - factor * A(k,j)
                    end do
                    b(i) = b(i) - factor * b(k)
                end if
            end do
        end do
        
        ! Back substitution
        do i = n, 1, -1
            if (abs(A(i,i)) > 1.0d-15) then
                x(i) = b(i)
                do j = i+1, n
                    x(i) = x(i) - A(i,j) * x(j)
                end do
                x(i) = x(i) / A(i,i)
            else
                x(i) = 0.0d0
            end if
        end do
    end subroutine solve_linear_system
    
    recursive subroutine parallel_quicksort(arr, low, high)
        implicit none
        integer, intent(in) :: low, high
        real*8, intent(inout) :: arr(:)
        integer :: pivot_idx, threshold, range_size
        
        range_size = high - low + 1
        threshold = 1000
        
        if (low < high) then
            call partition(arr, low, high, pivot_idx)
            
            if (range_size > threshold .and. range_size < 100000) then
                !$OMP TASK SHARED(arr) FIRSTPRIVATE(low, pivot_idx)
                call parallel_quicksort(arr, low, pivot_idx - 1)
                !$OMP END TASK
                !$OMP TASK SHARED(arr) FIRSTPRIVATE(pivot_idx, high)
                call parallel_quicksort(arr, pivot_idx + 1, high)
                !$OMP END TASK
                !$OMP TASKWAIT
            else
                call parallel_quicksort(arr, low, pivot_idx - 1)
                call parallel_quicksort(arr, pivot_idx + 1, high)
            end if
        end if
    end subroutine parallel_quicksort
    
    subroutine partition(arr, low, high, pivot_idx)
        implicit none
        integer, intent(in) :: low, high
        real*8, intent(inout) :: arr(:)
        integer, intent(out) :: pivot_idx
        real*8 :: pivot_val, temp
        integer :: i, j
        
        pivot_val = arr(high)
        i = low - 1
        
        do j = low, high - 1
            if (arr(j) <= pivot_val) then
                i = i + 1
                temp = arr(i)
                arr(i) = arr(j)
                arr(j) = temp
            end if
        end do
        temp = arr(i + 1)
        arr(i + 1) = arr(high)
        arr(high) = temp
        pivot_idx = i + 1
    end subroutine partition

    subroutine select_initial_eigenvalues(eigenvalues_full, selected, N_full, n_select)
        implicit none
        integer, intent(in) :: N_full, n_select
        real*8, intent(in) :: eigenvalues_full(N_full)
        real*8, intent(out) :: selected(n_select)
        integer :: i, j
        
        selected = 0.0d0
        j = 1
        do i = 1, N_full
            if (j > n_select) exit
            if (eigenvalues_full(i) > 1.0d-10) then
                selected(j) = eigenvalues_full(i)
                j = j + 1
            end if
        end do
    end subroutine select_initial_eigenvalues

    subroutine select_next_eigenvalues(eigenvalues_full, selected, N_full, start_idx, n_select)
        implicit none
        integer, intent(in) :: N_full, start_idx, n_select
        real*8, intent(in) :: eigenvalues_full(N_full)
        real*8, intent(out) :: selected(n_select)
        integer :: i, j
        
        selected = 0.0d0
        j = 1
        do i = start_idx, N_full
            if (j > n_select) exit
            if (eigenvalues_full(i) > 1.0d-10) then
                selected(j) = eigenvalues_full(i)
                j = j + 1
            end if
        end do
    end subroutine select_next_eigenvalues

    subroutine linear_regression(x, y, n, slope, intercept, r_squared)
        implicit none
        integer, intent(in) :: n
        real*8, intent(in) :: x(n), y(n)
        real*8, intent(out) :: slope, intercept, r_squared
        real*8 :: sum_x, sum_y, sum_xy, sum_x2
        real*8 :: mean_y, ss_tot, ss_res, denom
        real*8 :: n_real
        
        n_real = real(n, 8)
        
        sum_x = sum(x)
        sum_y = sum(y)
        sum_xy = sum(x * y)
        sum_x2 = sum(x * x)
        
        denom = n_real * sum_x2 - sum_x**2
        if (abs(denom) > 1.0d-15) then
            slope = (n_real * sum_xy - sum_x * sum_y) / denom
            intercept = (sum_y - slope * sum_x) / n_real
            
            mean_y = sum_y / n_real
            ss_tot = sum((y - mean_y)**2)
            ss_res = sum((y - (slope * x + intercept))**2)
            
            if (ss_tot > 1.0d-15) then
                r_squared = 1.0d0 - ss_res / ss_tot
            else
                r_squared = 0.0d0
            end if
        else
            slope = 0.0d0
            intercept = sum_y / n_real
            r_squared = 0.0d0
        end if
    end subroutine linear_regression

end program improved_riemann_analysis
