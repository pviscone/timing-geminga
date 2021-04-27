module mod
    implicit none 
    contains   
    subroutine matrix(ph_time,f0_arr,f1_arr)
    ! ph_time= array tempi fotoni, N=# fotoni, n_harm=#armoniche, f0_arr=array_f0, size_f0=len(f0_arr), f1_arr=array_f1, size_f1=len(f1_arr), result=matrice degli z test OUTPUT
        real, dimension(30957), intent(in) :: ph_time
        real, dimension(500), intent(in) :: f0_arr
        real, dimension(50), intent(in) :: f1_arr
        !real, dimension(500,50), intent(out)  :: result
        !real, dimension(30957) :: phi_arr
        real, dimension(:), allocatable :: phi_arr       
        integer :: size_f0, size_f1, i, j
        real (kind=16):: t0
      
        allocate(phi_arr(30957))
        size_f0=500
        size_f1=50        
        t0=250992001.      

        do i=1, 500
            do j=1, 50
                phi_arr=modulo((f0_arr(i)*(ph_time-t0)+(f1_arr(j)*(ph_time-t0)**2)/2),1.)
                print *, phi_arr   
            
            end do 
        end do            
    end subroutine matrix
end module mod

