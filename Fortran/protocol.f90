
      subroutine protocol(x, deltax,L,p)

      real*8, intent(in) :: x,L, deltax
      real*8, intent(out)::p
      real*8:: a, ll, lu
      a=1.d0
      lu=modulo(x+deltax,L)
      ll=modulo(x-deltax,L)
      if (ll.LE. lu) then
         if (lu .LE.a) then
            p=(lu-ll)/(2.d0*deltax)
         else
            if (ll .LE. a) then
               p=(a-ll)/(2.d0*deltax)
            else
               p=0.d0
            endif
         endif
      else
         if (ll .LE. a) then
            p=(a-ll+lu)/(2.d0*deltax)
         else
            if (lu .LE. a) then
               p=lu/(2.d0*deltax)
            else
               p=a/(2.d0*deltax)
            endif
         endif
      endif
      
     
      end subroutine protocol
