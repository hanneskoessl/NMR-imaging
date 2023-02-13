subroutine random_init(mag,Nimax,Njmax)
  integer Nimax,Njmax
  integer mag(Nimax,Njmax)
 
  integer i,j
  
  do i=1,Nimax
     do j=1,Njmax
	mag(i,j) = INT(rand()*5.d0)+1
     enddo
  enddo

  return
end subroutine random_init


program simulatedAnnealing
  external CalcE


!dimension of the image
  integer Nimax,Njmax,NN,dummyi,dummyj
  parameter (Nimax=254,Njmax=333,NN=4)
  real*8 JJ                    !/* Coupling constant */
  parameter (JJ=1.45d0)


!mean values and standard deviation
  integer MW(5)
!                 BG  WM  GM  CSF  SB
  parameter (MW=(/30,426,602,1223,167/))
  integer SIGMA(5)
  parameter (SIGMA=(/30,59,102,307,69/))

  integer MR(Nimax,Njmax), mag(Nimax,Njmax), z(Nimax,Njmax)           	!/* 2D Ising Lattice */
  integer i, j, k, m, n, sweeps, accept, move      				!/* Loop counters */
  integer Inn(NN)                    						!/* Nearest neighbor array I */
  parameter (Inn=(/1,-1,0,0/))
  integer Jnn(NN)                       						!/* Nearest neighbor array J */
  parameter (Jnn=(/0,0,1,-1/)) 
  integer Inew, Jnew               							!/* Nearest neighbot indices */ 
  real*8 T, Tf, lambda, En
  integer seed, correct(0:5), gesamt(0:5)
  CALL SYSTEM_CLOCK(seed)
  call srand(seed)
  T=27.
  Tf=0.01
  lambda=1.1

  open(unit=11, file='./CorrectSegImage.dat', status='old')
  !read in image
  do i=1,Nimax
     do j=1,Njmax
        read(11,'(3I5)') dummyi,dummyj,MR(i,j)
     enddo
     read (11,*)
  enddo
  close(11)

  open(unit=12, file='./SimMRimage.dat', status='old')
  !read in image
  do i=1,Nimax
     do j=1,Njmax
        read(12,'(3I5)') dummyi,dummyj,z(i,j)
     enddo
     read (12,*)
  enddo
  close(12)

  call random_init(mag,Nimax,Njmax) 
  write (*,*) CalcE(mag,z,Jnn,Inn,NN,JJ,Nimax,Njmax,MW,SIGMA)



  do while (T > Tf)
    sweeps=200000
    do m=1,sweeps,1
      call sweep(mag,z,Jnn,Inn,NN,JJ,Nimax,Njmax,T,accept,MW,SIGMA)
    end do
    T=T/lambda
  end do

  correct=0
  gesamt=0
  do k=1,5
    do i=1,Nimax
        do j=1,Njmax
          if ((mag(i,j)==k) .and. (MR(i,j)==k)) then
            correct(k)=correct(k)+1
          endif
          if (MR(i,j)==k) then
            gesamt(k)=gesamt(k)+1
	  endif  
        enddo
     enddo
   enddo
open(unit=12, file='SegSA1.dat')
!write out image
   do i=1,Nimax
      do j=1,Njmax
            write(12,'(3I5)') i,j,mag(i,j)
      enddo
      write (12,*)
   enddo
  close(12)

write (*,*) 'Gewebe 1: richtig=' , correct(1), ', gesamt=',gesamt(1), ', Fehler=', 1.-DBLE(correct(1))/gesamt(1)
write (*,*) 'Gewebe 2: richtig=' , correct(2), ', gesamt=',gesamt(2), ', Fehler=', 1.-DBLE(correct(2))/gesamt(2)
write (*,*) 'Gewebe 3: richtig=' , correct(3), ', gesamt=',gesamt(3), ', Fehler=', 1.-DBLE(correct(3))/gesamt(3)
write (*,*) 'Gewebe 4: richtig=' , correct(4), ', gesamt=',gesamt(4), ', Fehler=', 1.-DBLE(correct(4))/gesamt(4)
write (*,*) 'Gewebe 5: richtig=' , correct(5), ', gesamt=',gesamt(5), ', Fehler=', 1.-DBLE(correct(5))/gesamt(5)

correct(0)=0.
gesamt(0)=0.
do i=1,5
  correct(0)=correct(0)+correct(i)
  gesamt(0)=gesamt(0)+gesamt(i)
enddo
write (*,*) 'gesamt richtig=', correct(0), 'gesamt=',gesamt(0), 'Fehler gesamt=', 1-DBLE(correct(0))/gesamt(0)

write (*,*) CalcE(mag,z,Jnn,Inn,NN,JJ,Nimax,Njmax,MW,SIGMA)

end




subroutine sweep(mag,z,Jnn,Inn,NN,JJ,Nimax,Njmax,T,accept,MW,SIGMA)
  integer Nimax, Njmax,MW(5),SIGMA(5)
  integer mag(Nimax,Njmax),z(Nimax,Njmax),Jnn(NN),Inn(NN),L,NN, accept, magalt
  
  real*8 JJ, E, Etemp, deltaE, r, T
  integer i,j,k,Inew,Jnew
  integer seed
  

    i=int(rand()*Nimax)+1			!random number between (0 and L) +1 so that all are equally likely, then convert to integer
    j=int(rand()*Njmax)+1
  
    E=EnergyNN(i,j,mag,z,Jnn,Inn,NN,JJ,Nimax,Njmax,MW,SIGMA)
    magalt=mag(i,j)
    mag(i,j)=INT(rand()*5.d0)+1 		!randomly change tissue type

    Etemp=EnergyNN(i,j,mag,z,Jnn,Inn,NN,JJ,Nimax,Njmax,MW,SIGMA)

    deltaE=Etemp-E
!    write (*,*) E, Etemp, deltaE
    !check for floating point overflow
    if ((deltaE)/T.gt.20.d0) then
      r=0.d0
    else
      if ((deltaE)/T.lt.0.05d0) then
        r=1.d0
      else
        r=exp(-(deltaE)/T)
      end if
    end if

    if(rand() < r) then
      mag(i,j)=mag(i,j)
      accept=accept+1
    else
      mag(i,j)=magalt
    end if

  return 
end subroutine sweep


real function EnergyNN(i,j,mag,z,Jnn,Inn,NN,JJ,Nimax,Njmax,MW,SIGMA)
  integer Nimax, Njmax,MW(5),SIGMA(5),mittel,standardabweichung
  integer mag(Nimax,Njmax),z(Nimax,Njmax),Jnn(NN),Inn(NN),NN
  
  real*8 JJ
  integer i,j,k,Inew,Jnew
  real*8 Energy

  !/* Determine the initial energy */
  Energy = 0.0
  
	!/* Loop over nearest neighbors */
	do k=1, NN  
            Inew = i + Inn(k)       
            Jnew = j + Jnn(k)
    
         !new boundary condiition (period possible since BG->BG)
         !/* Check periodic boundary conditions */
              if (Inew .le. 0) then
                 Inew = Nimax
              else 
                 if(Inew .gt. Nimax ) then 
                    Inew = 1
                 endif
              endif
              if (Jnew .le. 0) then
                 Jnew = Njmax
              else 
                 if(Jnew .gt. Njmax) then 
                    Jnew = 1
                 endif
              endif

	    
	      !/* Update the energy */
            if (mag(i,j)==mag(Inew,Jnew)) then
 	       Energy = Energy-JJ
            end if
        enddo
        mittel=MW(mag(i,j))
        standardabweichung=SIGMA(mag(i,j))
        Energy=Energy+(z(i,j)-mittel)**2/(2.*standardabweichung**2)+log(DBLE(standardabweichung))
   EnergyNN=Energy
   return
end function EnergyNN

real function CalcE(mag,z,Jnn,Inn,NN,JJ,Nimax,Njmax,MW,SIGMA)
  integer Nimax, Njmax,MW(5),SIGMA(5),mittel,standardabweichung
  integer mag(Nimax,Njmax),z(Nimax,Njmax),Jnn(NN),Inn(NN),NN
  
  real*8 JJ
  integer i,j,k,Inew,Jnew
  real*8 Energy

  ! Determine the initial energy 
  Energy = 0.0
  
  do i=1,Nimax
     do j=1,Njmax
        
	! Loop over nearest neighbors
	do k=1, NN  
            Inew = i + Inn(k)       
            Jnew = j + Jnn(k)
    
!new boundary condiition (period possible since BG->BG)
         ! Check periodic boundary conditions
              if (Inew .le. 0) then
                 Inew = Nimax
              else 
                 if(Inew .gt. Nimax ) then 
                    Inew = 1
                 endif
              endif
              if (Jnew .le. 0) then
                 Jnew = Njmax
              else 
                 if(Jnew .gt. Njmax) then 
                    Jnew = 1
                 endif
              endif

	    
	    ! Update the energy
            if (mag(i,j)==mag(Inew,Jnew)) then
 	       Energy = Energy-JJ
            end if
        enddo
        mittel=MW(mag(i,j))
        standardabweichung=SIGMA(mag(i,j))
        Energy=Energy+2.0*((z(i,j)-mittel)**2/(2.*standardabweichung**2)+log(DBLE(standardabweichung)))
     enddo
   enddo

   ! Account for double counting
   Energy = Energy/2.d0;

   CalcE=Energy
   return
 end function CalcE

