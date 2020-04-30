subroutine f(s, E, mu, rho1, H1, km, GSK,n)
implicit none 
complex*16, intent(in) :: s
integer, intent(in) :: N
complex*16,intent(out) :: GSK(n*2,n*2)
complex*16 :: fn, disp(4,4), T(4,4)
real*8, parameter :: pi = 4.D0*DATAN(1.D0)
integer ::  mode, NL, IS1, j, i,pa
real*8, intent(in) ::  E(N), mu(N), rho1(N),H1(N), km
real*8 :: uniload(N*2), rho, H, unit1(4,4), d0, d1, mm
complex*16 :: CI, T3inv(2,2), D3inv(2,2),uhat(5),kp, ks, S3INV(4,4), stress(2),lamda
complex*16 :: D3(2,2),T3(2,2), SK3(2,2), lmu,cp2,cs2, S2INV(4,4),sk(4,4)


!f2py intent(in) s, E, mu, rho1, H1, km,fm, r,load,m,n
!f2py intent(out) GSK
!f2py depend(n) E , mu , rho1 ,H1,GSK, s4inv, uniload 

   
    
   CI = (0.0 , 1.0) 
   pa= n*2
   unit1 = 0 
   unit1(1,1) = -1.0 
   unit1(2,2) = -1.0 
   unit1(3,3) = 1.0
   unit1(4,4) = 1.0

   GSK = 0 

  DO NL = 1,n                           !for different layers
                 
	H   = H1(NL)
      rho = rho1(NL)
             
      lamda = E(NL)*mu(NL)/ ((1-2*mu(NL))*(1+mu(NL)))  !*(1+damp*s*CI) 
	lmu  = (E(NL)/(2*(1+mu(NL))))                    !*(1+damp*s*CI) 
   
	CP2= (lamda+2.0*lmu)/rho         
      CS2 =(lmu)/rho 

      ! kp = sqrt( s*ci/CP2 + km**2)    
      ! ks = sqrt( s*ci/CS2 + km**2)

      kp = sqrt( (s*ci)**2/CP2 + km**2)    !s represents the frequency in this case, not the laplace  variable
      ks = sqrt( (s*ci)**2/CS2 + km**2)    
        
     if (kp .NE. 0 .and. ks .NE. 0 ) then 
 
               
      IF (NL<=n-1) THEN
          
      DISP(1,1) = -km
      DISP(1,2) = -km*EXP(-kp*H)
      DISP(1,3) =  ks
      DISP(1,4) = -ks*EXP(-ks*H)

      DISP(2,1) = -kp
      DISP(2,2) =  kp*EXP(-kp*H)
      DISP(2,3) =  km 
      DISP(2,4) =  km*EXP(-ks*H)
   

      DISP(3,1) = -km*EXP(-kp*H)
      DISP(3,2) = -km
      DISP(3,3) =  ks*EXP(-ks*H)
      DISP(3,4) = - ks

      DISP(4,1) = -kp*EXP(-kp*H)
      DISP(4,2) =  kp 
      DISP(4,3) = km*EXP(-ks*H)
      DISP(4,4) = km

!	*******************STRESS MATRIX**********************************
      T(1,1) =  2.0*kp*km*lmu
      T(1,2) =  -2.0*kp*km*lmu*EXP(-kp*H)
      T(1,3) = -(ks**2+km**2)*lmu
      T(1,4) = -(ks**2+km**2)*EXP(-ks*H)*lmu

      T(2,1) =  (kp**2*(lamda+2*lmu)-km**2*lamda) 
      T(2,2) =  (kp**2*(lamda+2*lmu)-km**2*lamda)*EXP(-kp*H)
      T(2,3) = -2*km*ks*lmu
      T(2,4) = 2*km*ks*EXP(-ks*H)*lmu

      T(3,1) = 2.0*kp*EXP(-kp*H)*km*lmu
      T(3,2) = -2.0*kp*km*lmu
      T(3,3) = -(ks**2+km**2)*EXP(-ks*H)*lmu
      T(3,4) = -(ks**2+km**2)*lmu

      T(4,1) =  ((kp**2*(lamda+2*lmu)-km**2 *lamda))*EXP(-kp*H)
      T(4,2) =  kp**2*(lamda+2*lmu)-km**2*lamda
      T(4,3) = -2*km*ks*EXP(-ks*H)*lmu
      T(4,4) =  2*km*ks *lmu

        
!	********FIND INVERSE OF DISPLACEMENT MATRIX**********************

      call zInverse(DISP, S2INV, 4)
           
      SK = matmul(T,S2INV)   !	MULTIPLICATION TO GET STIFFNESS MATRIX OF AN ELEMENT
      
      SK = MATMUL(UNIT1, SK)
     
      IS1=2*(NL-1)  ! NL IS number
      
      DO I = 1,4
          DO J = 1,4
             GSK(I+IS1,J+IS1) = GSK(I+IS1,J+IS1) + SK(I,J) 
          ENDDO
      ENDDO

!!********************FOR LAST LAYER - 2DEGREE OF FREEDOM****************
    ELSE
          
           D3(1,1) = -km
           D3(1,2) = ks 
           D3(2,1) = -kp
           D3(2,2) = km

		T3(1,1) =  2.0*lmu*km*kp
		T3(1,2) = -lmu *(ks**2+km**2)
		T3(2,1) = (kp**2*(lamda+2*lmu)-km**2*lamda) 
		 T3(2,2) = -2.0*lmu *km*ks

           call zInverse(D3, d3inv, 2)
           SK3 = matmul(T3,D3inv)       
           
           SK3 = -SK3
          
	     GSK(pa-1,pa-1) = GSK(pa-1,pa-1)+ SK3(1,1)
           GSK(pa-1,pa) = GSK(pa-1,pa)+ SK3(1,2)
           GSK(pa,pa-1) = GSK(pa,pa-1)+ SK3(2,1)
           GSK(pa,pa) = GSK(pa,pa)+ SK3(2,2)
    END IF   !FOR LAYERS
    
    else  
           GSK = 0 

    END IF 
   
    END DO		 !END NL LOOP 

   
   
                                               
end subroutine 


subroutine zInverse(a,ra,n)
      integer, intent(in)::n
      integer:: lda,ipiv(n),info,lwork
      complex*16, intent(in)::a(n,n)
	complex*16, intent(out)::ra(n,n)
	complex*16 :: work(n)
	  !f2py intent(in) n, a
      !f2py intent(out) ra
      !f2py depend(n) a, ra	, ipiv
      ra=a
      lwork=n
      lda=n
      call zgetrf(n, n, ra, lda, ipiv, info)
      if(info/=0) then 
            ra = 0  
      end if  
      call zgetri(n, ra, lda, ipiv, work, lwork, info)
      if(info/=0) then 
            ra = 0 
      end if  
end subroutine
