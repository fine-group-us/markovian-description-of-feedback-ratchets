  !     The subroutine compute a frequency vector depending on the number
  !     of bins used in the discretization of the position 
  !---------------------------------------------------

      subroutine joint_prob(x,L, nbins,Px)

      integer :: i,k,n
      integer, intent(in) :: nbins
      real*8, intent(in) :: x,L
      double precision, intent(out):: Px(nbins)
      real*8 ::  delta_L
      
  
      delta_L=L/dble(nbins)
      n=floor(dble(x)/dble(delta_L))
      Px(n+1)=Px(n+1)+1.d0
     
     
      
     
      end subroutine joint_prob

      
