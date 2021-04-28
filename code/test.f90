module mod
    implicit none 
    contains   
    subroutine matrix(ph_time,f0_arr,f1_arr,result)
    ! ph_time= array tempi fotoni, N=# fotoni, n_harm=#armoniche, f0_arr=array_f0, size_f0=len(f0_arr), f1_arr=array_f1, size_f1=len(f1_arr), result=matrice degli z test OUTPUT
        real*8, dimension(30957), intent(in) :: ph_time
        real*8, dimension(500), intent(in) :: f0_arr
        real*8, dimension(50), intent(in) :: f1_arr
        real*8, dimension(500,50), intent(out)  :: result
        !real, dimension(30957) :: phi_arr
        real*8, dimension(:), allocatable :: phi_arr       
        !integer*8 :: N, n_harm, size_f0, size_f1, i, j, k, p
        integer*8 :: i,j,k,p 
        !real*16 ::  t0, cs, sn, ztest
        real*8 :: ztest, cs, sn
        allocate(phi_arr(30957))
        !N=30957
        !n_harm=10
        !size_f0=50
        !size_f1=50        
        !t0=250992001.      
        
        do i=1, 500
            do j=1, 50
                ztest=0.
                phi_arr=modulo((f0_arr(i)*(ph_time)+(f1_arr(j)*(ph_time)**2)/2),1.)
                
                phi_arr=phi_arr*2*3.1415927
                
                do k=1, 10
                    cs=0.
                    sn=0.  

                    do p=1, 30957             

                        cs=cs+cos(k*phi_arr(p))
                        sn=sn+sin(k*phi_arr(p))
                        !print *,cs        
                        
                    end do
                
                ztest=ztest+cs**2+sn**2     
                !print *, ztest 
                
                end do     
            result(i,j)=ztest
            end do 
        end do            
    end subroutine matrix
end module mod




