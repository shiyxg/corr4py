subroutine do_correlate2(d, wavelet, ne,ns,nt_d, nt_w,start_coor, end_coor, result)
    ! =====================================================
    ! HSR correlate with the wavetet
    ! =====================================================
    integer, intent(in)   :: ne,ns,nt_d, nt_w, start_coor, end_coor
    real*8, intent(in)  :: d(ne,ns,nt), wavelet(:,:)
    real*8, intent(out)  :: result(end_coor-start_coor, ns, ne)
    integer*8 :: i,j,k,l

    write(*,*) shape(d), shape(wavelet)
    
    do i = 1,ne,1
        do j=1, ns,1
            do k=1,end_coor-start_coor,1
                result(k,j,i) = 0
            end do
        end do
    end do
    CALL omp_set_num_threads(24)
    !$omp parallel do private(l,k,j,i)
    do i = 1,ne,1
        do j=1, ns,1
            do k=1,end_coor-start_coor,1
                do l=1,nt_w,1
                    result(k,j,i) = result(k,j,i) + d(k+l+start_coor-1,j,i)*wavelet(l,i)
                end do
            end do
        end do
    end do
    !$omp end parallel do
return
end subroutine do_correlate2

program hello
    integer :: ne,ns,nt_d, nt_w,start_coor, end_coor
    interface
        subroutine do_correlate2(d, wavelet, ne,ns,nt_d, nt_w,start_coor, end_coor, result)
            integer, intent(in)   :: ne,ns,nt_d, nt_w, start_coor, end_coor
            real*8, intent(in)  :: d(:,:,:), wavelet(:,:)
            real*8, intent(out)  :: result(end_coor-start_coor, ns, ne)
            integer :: i,j,k,l
        end subroutine do_correlate2
    end interface
    
    real*8, allocatable:: d(:,:,:), d1(:,:)
    real*8, allocatable:: wavelet(:,:), w1(:)
    real*8, allocatable:: result(:,:,:), r1(:,:)
    integer :: i,j,k

    ne=100
    ns=2000
    nt_d=20000
    nt_w = 500
    start_coor=10000
    end_coor=10500
    allocate(d(nt_d, ns, ne))
    allocate(wavelet(nt_w, ne))
    allocate(result(end_coor-start_coor, ns, ne))

    do i = 1,ne,1
        do j=1, ns,1
            do k=1,nt_d
                d(k,j,i)=i+0.1
            end do
        end do
    end do
    do i = 1,ne,1
        do j=1, nt_w
            wavelet(j,i)=i+0.1
        end do
    end do
    call do_correlate2(d, wavelet, ne,ns,nt_d, nt_w,start_coor, end_coor, result)
    write(*,*) result(2,2,1)
    deallocate(d)
    deallocate(wavelet)
    deallocate(result)

    ns=10
    nt=8*100*3600
    start_coor=-3000
    end_coor=3000
    allocate(d1(nt, ns))
    allocate(w1(nt))
    allocate(r1(end_coor-start_coor, ns))

    do j=1, ns,1
        do k=1,nt
            d1(k,j)=j-1
        end do
    end do
    do j=1, nt,1
        w1(j)=j
    end do
    call do_correlate1(d1, w1, ns,nt,start_coor, end_coor, r1)
    write(*,*) r1(1,1),d1(1,1),d1(1,2)
    write(*,*) r1(1,2)
    deallocate(d1)
    deallocate(w1)
    deallocate(r1)
end