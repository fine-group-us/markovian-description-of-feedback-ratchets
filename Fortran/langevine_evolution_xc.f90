program langevine_evolution_xc
  implicit none

  real ran3
  integer*8 :: Nsim, Nstep, nummed,nt,n,i,j, iseed,k,c,r, nbins
  real*8 :: gamma, beta, D, deltatm, deltat,t1,t2,t3,Fext,V0,a,L, Fin, Feff, u1, u2, norm, pi,p, dice, deltax
  parameter(Nsim=1000000, Nstep=int(10**2), deltat=10.d0**(-4), deltatm=deltat*dble(Nstep), Fext=0.d0, V0=5.d0, a=1.d0, nummed=50, nt=Nstep*nummed, nbins=300)
  parameter(gamma=1.d0, beta=1.d0, D=1.d0,L=2.d0, deltax=0.5d0) 
  logical :: control(Nsim)
  real*8 :: x(Nsim), PC1x(nbins),PC0x(nbins), PC

  
  iseed = -2345678
  pi=4.d0*atan(1.d0)
  ! Integration time
  t1=L**2/(2*D)
  t2=gamma*L*a/V0
  t3=gamma*L*(L-a)/V0
 

  
  ! Output 
 
  open(24,file='datos_de_P_C_1.dat',status="unknown")
  open(25,file='datos_de_P_C_0.dat',status="unknown")
  open(26,file='datos_de_P_C.dat',status="unknown")

 ! Initial condition: an uniform distribution
       do k=1,Nsim
          x(k)=L*ran3(iseed)
          control(k)=.FALSE.     
       enddo
       Fin=0.d0
  
      temporal: do i=1,nt
         do k=1,Nbins
         PC1x(k)=0.d0
         PC0x(k)=0.d0 
      enddo
      PC=0.d0
      simulation: do n=1,Nsim
          measure: if  (modulo(i,Nstep)==0) then ! perform a measure 
             ! Apply the protocol
             call protocol(x(n), deltax,L,p)
             dice=ran3(iseed)
           if(dice .GE. p) then !switch off 
              control(n) = .FALSE.
              Fin=0.d0    
           else ! switch on
              control(n)= .TRUE.
              Fin=V0/(L-a)
           endif
           ! Build the joint probability distribution for the process (x,c) and for (c)
           if (control(n)==.TRUE.)then
               call joint_prob(x(n), L, nbins, PC1x)
               PC=PC+1.d0
            else
              call joint_prob(x(n), L, nbins, PC0x)
              endif
            u1 = ran3(iseed)
            u2 = ran3(iseed)
           do while (u1==0.d0) 
               u1 = ran3(iseed)
            end do
            norm = sqrt(-2.d0*log(u1))*cos(2.d0*pi*u2)
           Feff=Fin-Fext
           x(n)=x(n)+Feff*deltat+norm*sqrt(2.d0)*(deltat)**(0.5d0)
           x(n)=modulo(x(n),L)
           
         else !no measure is performed

             if (control(n)==.FALSE.) then
                Fin=0.d0
              else 
                 potential:if(x(n)>a) then  
                 Fin=-V0/a
              else
                  Fin=V0/(L-a)
               endif potential
            endif 
            u1 = ran3(iseed)
            u2 = ran3(iseed)
            do while(u1==0.d0) 
               u1 = ran3(iseed)
            end do
              ! Build the joint probability distribution for the process (x,c) and for (c)
            if (control(n)==.TRUE.)then
               call joint_prob(x(n), L, nbins, PC1x)
               PC=PC+1.d0
            else
               call joint_prob(x(n), L, nbins, PC0x)
               endif
            norm = sqrt(-2.d0*log(u1))*cos(2.d0*pi*u2)
           Feff=Fin-Fext
           x(n)=x(n)+Feff*deltat+norm*sqrt(2.d0)*(deltat)**(0.5d0)
           x(n)=modulo(x(n),L)
                
        endif measure
      
      enddo simulation
      
   ! Save the results
      write(24,199) PC1x/Nsim
      write(25,199) PC0x/Nsim
      write(26,199) PC/Nsim
     enddo temporal

   199   format(1000000000(e14.7,2x))  
  
   end program langevine_evolution_xc
