subroutine corr_fft_full(d_1d, w_1d, ne,ns,nt, n_thread, c_1d)
    ! =====================================================
    ! HSR correlate with the wavetet
    ! =====================================================
    use, intrinsic :: iso_c_binding 
    include 'fftw3.f03'
    
    integer*8, intent(in)  :: ne,ns,nt, n_thread
    real*8,  intent(in)  :: d_1d(nt*ns*ne)
    real*8,  intent(in)  :: w_1d(nt*ns*ne)
    real*8,  intent(out) :: c_1d(nt*ns*ne)
    real*8 :: d(nt,ns,ne)
    real*8 :: w(nt,ns,ne)
    real*8 :: c(nt,ns,ne)
    integer*8 :: i,j,k,l
    integer :: nt_f
    
    ! FFTW
    type(C_PTR) :: forward, backward
    real*8, dimension(nt) :: w_ij, d_ij, c_ij, in
    double complex, dimension(nt/2+1) :: w_ij_fft, d_ij_fft, c_ij_fft, out
    nt_f = nt
    
    CALL omp_set_num_threads(n_thread)
    ! init
    !$omp parallel do private(k,j,i)
    do i = 1,ne,1
        do j=1, ns,1
            do k=1,nt,1
                d(k,j,i) = d_1d( ((i-1)*ns+j-1)*nt+k)
                w(k,j,i) = w_1d( ((i-1)*ns+j-1)*nt+k)
            end do
        end do
    end do
    !$omp end parallel do
    
    forward  = fftw_plan_dft_r2c_1d(nt_f,in,out,FFTW_ESTIMATE)
    backward = fftw_plan_dft_c2r_1d(nt_f,out,in,FFTW_ESTIMATE)

    !$omp parallel do private(w_ij, d_ij, c_ij,w_ij_fft,d_ij_fft,c_ij_fft, j,i)
    do i = 1,ne,1
        do j=1, ns,1
            w_ij = w(:,j,i)
            d_ij = d(:,j,i)

            call fftw_execute_dft_r2c(forward, w_ij, w_ij_fft)
            call fftw_execute_dft_r2c(forward, d_ij, d_ij_fft)

            c_ij_fft = conjg(w_ij_fft)*d_ij_fft
 
            call fftw_execute_dft_c2r(backward,c_ij_fft, c_ij)

            c(:,j,i)=c_ij
            
        end do
    end do
    ! $omp end parallel do

    !$omp parallel do private(k,j,i)
    do i = 1,ne,1
        do j=1, ns,1
            do k=1,nt,1
                c_1d(((i-1)*ns+j-1)*nt+k) = c(k,j,i)
            end do
        end do
    end do
    !$omp end parallel do

    call fftw_destroy_plan(forward)
    call fftw_destroy_plan(backward)
return
end subroutine corr_fft_full

subroutine corr_fft_padding_f64(d_1d, w_1d, ne,ns,nt, n_thread, c_1d)
    ! =====================================================
    ! HSR correlate with the wavetet
    ! Thanks for https://github.com/lbovard/fftw-tests/blob/master/1d_real_fft.f90
    ! =====================================================
    use, intrinsic :: iso_c_binding 
    include 'fftw3.f03'
    
    integer*8, intent(in)  :: ne,ns,nt,n_thread
    real*8,  intent(in)  :: d_1d(nt*ns*ne)
    real*8,  intent(in)  :: w_1d(nt*ne)
    real*8,  intent(out) :: c_1d(nt*2*ns*ne)
    real*8 :: d(nt*2,ns,ne)
    real*8 :: w(nt*2,ne)
    real*8 :: c(nt*2,ns,ne)
    integer*8 :: i,j,k, STATUS,t1,t2,t3
    integer :: nt_f
    ! FFTW
    type(C_PTR) :: forward, backward
    real*8, dimension(nt*2) :: w_ij, d_ij, c_ij, in
    double complex, dimension(nt+1) :: w_ij_fft, d_ij_fft, c_ij_fft, out
    nt_f=nt
    CALL omp_set_num_threads(n_thread)
    ! init
    call system_clock(t1)
    !$omp parallel do private(k,j,i)
    do i = 1,ne,1
        do j=1, ns,1
            do k=1,nt*2,1
                if (k<=nt) then
                    d(k,j,i) = d_1d( ((i-1)*ns+j-1)*nt+k)
                else
                    d(k,j,i) = 0
                end if
            end do
        end do
    end do
    !$omp end parallel do

    !$omp parallel do private(j,i)
    do i = 1,ne,1
        do j=1,nt*2,1
            if (j<=nt) then
                w(j,i) = 0
            else
                w(j,i) = w_1d( (i-1)*nt+j-nt)
            end if
        end do
    end do
    !$omp end parallel do

    call system_clock(t2)
    forward  = fftw_plan_dft_r2c_1d(nt_f*2,in,out,FFTW_ESTIMATE)
    backward = fftw_plan_dft_c2r_1d(nt_f*2,out,in,FFTW_ESTIMATE)

    !$omp parallel do private(STATUS,w_ij_fft,d_ij_fft,c_ij_fft, j,i)
    do i = 1,ne,1
        ! call omp_get_thread_num(STATUS)
        ! write(*,*) i, STATUS
        do j=1, ns,1

            call fftw_execute_dft_r2c(forward, w(:,i), w_ij_fft)
            call fftw_execute_dft_r2c(forward, d(:,j,i), d_ij_fft)

            c_ij_fft = conjg(w_ij_fft)*d_ij_fft

            call fftw_execute_dft_c2r(backward,c_ij_fft, c(:,j,i))
        
            ! write(*,*) w(:,i)
            ! write(*,*) d(:,j,i)
            ! write(*,*) c(:,j,i)
        end do
    end do
    !$omp end parallel do


    !$omp parallel do private(k,j,i)
    do i = 1,ne,1
        do j=1, ns,1
            do k=1,nt*2,1
                c_1d(((i-1)*ns+j-1)*nt*2+k) = c(k,j,i)/(2.0*nt)
            end do
        end do
    end do
    !$omp end parallel do

    call system_clock(t3)
    write(*,*) "init time:",(t2-t1)/1e9,"corr time:",(t3-t2)/1e9

    call fftw_destroy_plan(forward)
    call fftw_destroy_plan(backward)

return
end subroutine corr_fft_padding_f64

subroutine corr_fft_padding_f32(d_1d, w_1d, ne,ns,nt, n_thread, c_1d)
    ! =====================================================
    ! HSR correlate with the wavetet
    ! Thanks for https://github.com/lbovard/fftw-tests/blob/master/1d_real_fft.f90
    ! =====================================================
    use, intrinsic :: iso_c_binding 
    include 'fftw3.f03'
    
    integer*8, intent(in)  :: ne,ns,nt,n_thread
    real*4,  intent(in)  :: d_1d(:)
    real*4,  intent(in)  :: w_1d(:)
    real*4,  intent(out) :: c_1d(nt*2*ns*ne)
    real*4 :: d(nt*2,ns,ne)
    real*4 :: w(nt*2,ne)
    real*4 :: c(nt*2,ns,ne)
    integer*8 :: i,j,k, STATUS,t1,t2,t3
    integer :: nt_f
    ! FFTW
    type(C_PTR) :: forward, backward
    real*4, dimension(nt*2) :: w_ij, d_ij, c_ij, in
    complex, dimension(nt+1) :: w_ij_fft, d_ij_fft, c_ij_fft, out

    nt_f=nt

    CALL omp_set_num_threads(n_thread)
    ! init
    call system_clock(t1)
    !$omp parallel do private(k,j,i)
    do i = 1,ne,1
        do j=1, ns,1
            do k=1,nt*2,1
                if (k<=nt) then
                    d(k,j,i) = d_1d( ((i-1)*ns+j-1)*nt+k)
                else
                    d(k,j,i) = 0
                end if
            end do
        end do
    end do
    !$omp end parallel do

    !$omp parallel do private(j,i)
    do i = 1,ne,1
        do j=1,nt*2,1
            if (j<=nt) then
                w(j,i) = 0
            else
                w(j,i) = w_1d( (i-1)*nt+j-nt)
            end if
        end do
    end do
    !$omp end parallel do

    call system_clock(t2)
    forward  = fftwf_plan_dft_r2c_1d(nt_f*2,in,out,FFTW_ESTIMATE)
    backward = fftwf_plan_dft_c2r_1d(nt_f*2,out,in,FFTW_ESTIMATE)

    !$omp parallel do private(w_ij_fft,d_ij_fft,c_ij_fft, j,i)
    do i = 1,ne,1
        do j=1, ns,1

            call fftwf_execute_dft_r2c(forward, w(:,i), w_ij_fft)
            call fftwf_execute_dft_r2c(forward, d(:,j,i), d_ij_fft)

            c_ij_fft = conjg(w_ij_fft)*d_ij_fft

            call fftwf_execute_dft_c2r(backward,c_ij_fft, c(:,j,i))
        
            ! write(*,*) w(:,i)
            ! write(*,*) d(:,j,i)
            ! write(*,*) c(:,j,i)
        end do
    end do
    !$omp end parallel do
    call fftwf_destroy_plan(forward)
    call fftwf_destroy_plan(backward)

    !$omp parallel do private(k,j,i)
    do i = 1,ne,1
        do j=1, ns,1
            do k=1,nt*2,1
                c_1d(((i-1)*ns+j-1)*nt*2+k) = c(k,j,i)/(2.0*nt)
            end do
        end do
    end do
    !$omp end parallel do

    call system_clock(t3)
    write(*,*) "init time:",(t2-t1)/1e9,"corr time:",(t3-t2)/1e9

return
end subroutine corr_fft_padding_f32

subroutine corr_fft_valid_f64(d_1d, w_1d, ne,ns,nt,sp,ep, n_thread, c_1d)
    ! =====================================================
    ! HSR correlate with the wavetet, 
    ! return valid corr result limited in sp*dt-ep*dt
    ! Thanks for https://github.com/lbovard/fftw-tests/blob/master/1d_real_fft.f90
    ! =====================================================

    interface
        subroutine corr_fft_padding_f64(d_1d, w_1d, ne,ns,nt,n_thread, c_1d)
            integer*8, intent(in)  :: ne,ns,nt,n_thread
            real*8,  intent(in)  :: d_1d(nt*ns*ne)
            real*8,  intent(in)  :: w_1d(nt*ne)
            real*8,  intent(out) :: c_1d(nt*ns*ne*2)
        end subroutine corr_fft_padding_f64
    end interface

    
    integer*8, intent(in)  :: ne,ns,nt,sp,ep,n_thread
    real*8,  intent(in)  :: d_1d(nt*ns*ne)
    real*8,  intent(in)  :: w_1d(nt*ne)
    real*8,  intent(out) :: c_1d((ep-sp)*ns*ne)

    real*8  :: c_full(2*nt*ns*ne)
    integer*8 :: i,j,k, idx_c1d,idx_cfull

    call corr_fft_padding_f64(d_1d, w_1d, ne,ns,nt, n_thread, c_full)

    !$omp parallel do private(idx_c1d, idx_cfull, k,j,i)
    do i = 1,ne,1
        do j=1, ns,1
            do k=sp,ep-1,1
                idx_c1d   = ((i-1)*ns+j-1)*(ep-sp) + k-sp+1
                idx_cfull = ((i-1)*ns+j-1)*( 2*nt) + k+nt+1
                c_1d(idx_c1d) = c_full(idx_cfull)
            end do
        end do
    end do
    !$omp end parallel do
return
end subroutine corr_fft_valid_f64

subroutine corr_fft_valid_f32(d_1d, w_1d, ne,ns,nt,sp,ep, n_thread, c_1d)
    ! =====================================================
    ! HSR correlate with the wavetet, 
    ! return valid corr result limited in sp*dt-ep*dt
    ! Thanks for https://github.com/lbovard/fftw-tests/blob/master/1d_real_fft.f90
    ! =====================================================

    interface
        subroutine corr_fft_padding_f32(d_1d, w_1d, ne,ns,nt,n_thread, c_1d)
            integer*8, intent(in)  :: ne,ns,nt,n_thread
            real*4,  intent(in)  :: d_1d(:)
            real*4,  intent(in)  :: w_1d(:)
            real*4,  intent(out) :: c_1d(nt*2*ns*ne)
        end subroutine corr_fft_padding_f32
    end interface

    
    integer*8, intent(in)  :: ne,ns,nt,sp,ep,n_thread
    real*4,  intent(in)  :: d_1d(:)
    real*4,  intent(in)  :: w_1d(:)
    real*4,  intent(out) :: c_1d((ep-sp)*ns*ne)

    real*4  :: c_full(2*nt*ns*ne)
    integer*8 :: i,j,k, idx_c1d,idx_cfull

    call corr_fft_padding_f32(d_1d, w_1d, ne,ns,nt, n_thread, c_full)

    !$omp parallel do private(idx_c1d, idx_cfull, k,j,i)
    do i = 1,ne,1
        do j=1, ns,1
            do k=sp,ep-1,1
                idx_c1d   = ((i-1)*ns+j-1)*(ep-sp) + k-sp+1
                idx_cfull = ((i-1)*ns+j-1)*( 2*nt) + k+nt+1
                c_1d(idx_c1d) = c_full(idx_cfull)
            end do
        end do
    end do
    !$omp end parallel do
return
end subroutine corr_fft_valid_f32

program hello
    interface
        subroutine corr_fft_full(d_1d, w_1d, ne,ns,nt, n_thread, c_1d)
            integer*8, intent(in)  :: ne,ns,nt, n_thread
            real*8,  intent(in)  :: d_1d(nt*ns*ne)
            real*8,  intent(in)  :: w_1d(nt*ns*ne)
            real*8,  intent(out) :: c_1d(nt*ns*ne)
        end subroutine corr_fft_full
        
        subroutine corr_fft_padding_f64(d_1d, w_1d, ne,ns,nt,n_thread, c_1d)
            integer*8, intent(in)  :: ne,ns,nt,n_thread
            real*8,  intent(in)  :: d_1d(nt*ns*ne)
            real*8,  intent(in)  :: w_1d(nt*ne)
            real*8,  intent(out) :: c_1d(nt*2*ns*ne)
        end subroutine corr_fft_padding_f64

        subroutine corr_fft_padding_f32(d_1d, w_1d, ne,ns,nt,n_thread, c_1d)
            integer*8, intent(in)  :: ne,ns,nt,n_thread
            real*4,  intent(in)  :: d_1d(nt*ns*ne)
            real*4,  intent(in)  :: w_1d(nt*ne)
            real*4,  intent(out) :: c_1d(nt*2*ns*ne)
        end subroutine corr_fft_padding_f32
    end interface
    
    real*8, allocatable:: d(:,:,:), w(:,:),c(:)
    real*4, allocatable:: d1(:,:,:), w1(:,:),c1(:)
    integer*8 :: i,j,k
    integer*8 :: ne,ns,nt
    integer*8 :: n_thread=10

    ne=3000
    ns=100
    nt=10000
    allocate(d(nt, ns, ne))
    allocate(w(nt, ne))
    allocate(c((nt*2)*ns*ne))

    call corr_fft_padding_f64(reshape(d,(/nt*ns*ne/)), reshape(w,(/nt*ne/)), ne,ns,nt,n_thread, c)

    deallocate(d)
    deallocate(w)
    deallocate(c)
    
    ne=3000
    ns=100
    nt=10000
    allocate(d1(nt, ns, ne))
    allocate(w1(nt, ne))
    allocate(c1((nt*2)*ns*ne))
    
    call corr_fft_padding_f32(reshape(d1,(/nt*ns*ne/)), reshape(w1,(/nt*ne/)), ne,ns,nt,n_thread, c1)

    deallocate(d1)
    deallocate(w1)
    deallocate(c1)

end