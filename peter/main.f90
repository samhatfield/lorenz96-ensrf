PROGRAM LORENZ95

!!!!!!!!!!!!!!!!!!!!!!!  
!Lorenz95 with RK4 time stepping and plenty of different parametrisation schemes.
!Based on a C++ code of Fenwick Cooper.
!Peter Dueben 27.03.13
!!!!!!!!!!!!!!!!!!!!!!!

  implicit none

!Parameter time stepping:
  REAL :: f_param, c_param
  INTEGER :: N_spinup, N_run, N_length, N_var, N_param, N_timestep
  REAL :: delta_T, deltat, sqrtdeltat
  REAL :: noise_amp, noise_cor, sigma0, sigma1, sigmaM, sigmaA, b0, b1, b2, b3 
  INTEGER i,j,k,l,m
  CHARACTER(len=70) :: output1

  REAL, POINTER, dimension(:,:) :: y
  REAL,POINTER, dimension(:,:) :: y_small
  REAL,POINTER, dimension(:) :: y_ini
  REAL,POINTER, dimension(:,:) :: y_ini_small
  REAL, POINTER, dimension(:) :: y0,ek,sigma 
  REAL :: df !No forcing = 0, with forcing = 5

!Set parameter:  
  f_param = 20.0
  c_param = 10.0
  N_spinup = 100000
  N_run = 1
  N_length = 40001
  N_var = 8
  N_param = 10
  N_timestep = 20
  delta_T = 0.01
  deltat = delta_T/real(N_timestep)
  sqrtdeltat = sqrt(deltat)
  df = 0.0

! N_param = 0: No parametrisation / no y's
! N_param = 1: Wilks 2005 deterministic (f=18,20, c=10)
! N_param = 2: Wilks 2005 AR0 (f=18,20, c=10)
! N_param = 3: Wilks 2005 AR1 (f=18,20, c=10)
! N_param = 4: Arnold et al. 2012 deterministic (f=20, c=4,10)
! N_param = 5: Arnold et al. 2012 AR0 (f=20, c=4,10)
! N_param = 6: Arnold et al. 2012 AR1 (f=20, c=4,10)
! N_param = 7: Arnold et al. 2012 state dependent (f=20, c=4,10)
! N_param = 8: Arnold et al. 2012 multiplicative (f=20, c=4,10)
! N_param = 9: Arnold et al. 2012 multiplicative additive (f=20, c=4,10)
! N_param = 10: Full system with y's

  CALL ini_para(N_param,noise_amp,noise_cor,sigma0,sigma1,sigmaM,sigmaA,b0,b1,b2,b3)

  allocate(y(N_length,N_var))
  allocate(y0(N_var))
  allocate(y_ini(N_var))
  allocate(ek(N_var))
  allocate(sigma(N_var))
  IF(N_param==10)THEN
     allocate(y_small(N_var,32))
     allocate(y_ini_small(N_var,32))
  END IF

  write(*,*) 'Start run.'

  DO i=1,N_run
     WRITE(output1,'(a,i0,a,i0,a)')'./Output/out_',N_param,'_',i, '.dat'
     open (60, FILE=(TRIM(output1)), STATUS='REPLACE', ACTION='write')


     !Setup initial conditions, if true all runs have different initial conditions
     IF(.FALSE..or.i==1)THEN
        DO j=1,N_var
           y0(j) = rand_normal(0.0,1.0)
           ek(j) = 0.0
           sigma(j) = sigma1*abs(y0(j))+sigma0
           IF(N_param==10)THEN
              DO k=1,32
                 y_small(j,k) = rand_normal(0.0,1.0)
              END DO
           END IF
        END DO

        !Spin up:
        DO j = 1,N_spinup
           !Calculate the deterministic part:
           CALL timestep(deltat,N_param,noise_amp,noise_cor,sigma,ek,sigma0,sigma1,sigmaM,sigmaA,&
                &b0,b1,b2,b3,y0,f_param,c_param,df,y_small)
        END DO

        IF(i==1)THEN
           y_ini = y0
           IF(N_param==10) y_ini_small = y_small
        END IF
     ELSE
        y0 = y_ini
        IF(N_param==10) y_small = y_ini_small
     END IF

 
     !Run:
     DO j= 1,N_length
        DO k = 1,N_timestep
           CALL timestep(deltat,N_param,noise_amp,noise_cor,sigma,ek,sigma0,sigma1,sigmaM,sigmaA,&
                &b0,b1,b2,b3,y0,f_param,c_param,df,y_small) 
        END DO
        y(j,:) = y0(:)
        write(60,"(9F10.2)") real(j)*delta_T,y0
     END DO
     
     close(60)
  END DO

  deallocate(y)
  deallocate(y_ini)
  deallocate(y0)
  deallocate(ek)
  deallocate(sigma)  
  IF(N_param==10)THEN
     deallocate(y_small)
     deallocate(y_ini_small)
  END IF

  write(*,*) 'Finish run.'

  contains


    SUBROUTINE timestep(deltat,N_param,noise_amp,noise_cor,sigma,ek,sigma0,sigma1,sigmaM,sigmaA,&
         &b0,b1,b2,b3,y0,f_param,c_param,df,y_small)
    
      IMPLICIT NONE

      INTEGER :: N_param,i,j
      REAL :: deltat,noise_amp,noise_cor,sigma0,sigma1,sigmaM,sigmaA,b0,b1,b2,b3,f_param,c_param,df
      REAL :: sigmaNew, up
      REAL,dimension(:):: y0, sigma, ek
      REAL :: yOut(size(y0(:)))
      REAL,POINTER, dimension(:,:) :: y_small
      REAL, dimension(:,:) ::  yOut_small(size(y_small(:,1)),size(y_small(1,:)))

      IF(N_param.lt.10)THEN
         CALL deterministic(deltat,y0,yOut,N_param,f_param,c_param,df)

         !Add the stochastic part:
         IF(N_param==0.or.N_param==1.or.N_param==4)THEN
            DO i=1,size(y0(:))
               y0(i) = yOut(i)
            END DO
         ELSE IF(N_param==2.or.N_param==5)THEN
            DO i=1,size(y0(:))
               y0(i) = yOut(i)+deltat*noise_amp*rand_normal(0.0,1.0) 
            END DO
         ELSE IF(N_param==3.or.N_param==6)THEN
            DO i=1,size(y0(:))
               ek(i) = noise_cor*ek(i)+noise_amp*sqrt(1.0-noise_cor*noise_cor)*rand_normal(0.0,1.0)
               y0(i) = yOut(i)+deltat*ek(i) 
            END DO
         ELSE IF(N_param==7)THEN
            DO i=1,size(y0(:))
               sigmaNew=sigma1*abs(y0(i))+sigma0
               ek(i) = sigmaNew/sigma(i)*noise_cor*ek(i)+sigmaNew &
                    & *sqrt(1.0-noise_cor*noise_cor)*rand_normal(0.0,1.0)
               sigma(i) = sigmaNew
               y0(i) = yOut(i)+deltat*ek(i) 
            END DO
         ELSE IF(N_param==8)THEN
            DO i=1,size(y0(:))
               up = b0+b1*y0(i)+b2*y0(i)*y0(i)+b3*y0(i)*y0(i)*y0(i)
               ek(i) = noise_cor*ek(i)+noise_amp*sqrt(1.0-noise_cor*noise_cor)*rand_normal(0.0,1.0)
               y0(i) = yOut(i)-deltat*up*ek(i) 
            END DO
         ELSE IF(N_param==9)THEN
            DO i=1,size(y0(:))
               up = b0+b1*y0(i)+b2*y0(i)*y0(i)+b3*y0(i)*y0(i)*y0(i)
               ek(i) = noise_cor*ek(i)+sqrt(1.0-noise_cor*noise_cor)*rand_normal(0.0,1.0)
               y0(i) = yOut(i)+deltat*ek(i)*(sigmaM*abs(up)+sigmaA) 
            END DO
         END IF
         
      ELSE IF(N_param==10)THEN
         CALL step_full_L95(deltat,y0,yOut,y_small,yOut_small,f_param,c_param,df)
         DO i=1,size(y0(:))
            y0(i) = yOut(i)
         END DO
         DO i=1,size(y_small(:,1))
            DO j=1,size(y_small(1,:))
               y_small(i,j) = yOut_small(i,j)
            END DO
         END DO
      END IF    

    END SUBROUTINE timestep


    SUBROUTINE deterministic(deltat,y,Out,N_param,f_param,c_param,df)
    
      IMPLICIT NONE

      REAL,dimension(:) :: y,Out
      REAL,dimension(:) :: dy(size(y(:)))
      REAL, dimension(:,:):: work(size(y(:)),3)
      REAL deltat, deltat_2, deltat_6,f_param,c_param,df
      INTEGER :: N_param
      INTEGER :: i

      deltat_2 = 0.5*deltat
      deltat_6 = deltat/6.0

      !RK4 timestepping:
      CALL discrete_L95(y,dy,N_param,f_param,c_param,df)
   
      DO i = 1,size(y(:))
         work(i,1)= y(i) + deltat_2 * dy(i)
      END DO
     
      CALL discrete_L95(work(:,1),work(:,2),N_param,f_param,c_param,df)

      DO i = 1,size(y(:))
         work(i,1)= y(i) + deltat_2 * work(i,2)
      END DO
        
      CALL discrete_L95(work(:,1),work(:,3),N_param,f_param,c_param,df)

      DO i = 1,size(y(:))
         work(i,1)= y(i) + deltat * work(i,3)
         work(i,3)=work(i,3)+work(i,2)
      END DO
    
      CALL discrete_L95(work(:,1),work(:,2),N_param,f_param,c_param,df)

      DO i = 1,size(y(:))
         Out(i)= y(i) + deltat_6 * (dy(i)+work(i,2)+2.0*work(i,3))
      END DO

    END SUBROUTINE deterministic
    
    SUBROUTINE discrete_L95(y1,y2,N_param,f_param,c_param,df)
    
      IMPLICIT NONE

      REAL,dimension(:) :: y1,y2
      REAL :: up(size(y1(:)))
      REAL :: f_param,c_param,df
      INTEGER :: N_param
      REAL :: c0, c1, c2, c3, c4
      INTEGER :: i
      INTEGER :: n
      
      n = size(y1(:))

      IF(N_param==0)THEN
         DO i=1,size(y1(:))
            up(i) = 0.0
         END DO
      ELSE IF(N_param==1.or.N_param==2.or.N_param==3)THEN
         IF(f_param==18.0)THEN
            c0= 0.275
            c1= 1.59
            c2=-0.0190
            c3=-0.0130
            c4= 0.000707
         ELSE IF(f_param==20.0)THEN
            c0= 0.262
            c1= 1.45
            c2=-0.0121
            c3=-0.00713
            c4= 0.000296
         END IF
         DO i=1,n
            up(i) = c0+c1*y1(i)+c2*y1(i)*y1(i)+c3*y1(i)*y1(i)*y1(i)+c4*y1(i)*y1(i)*y1(i)*y1(i)
         END DO
      ELSE IF(N_param==4.or.N_param==5.or.N_param==6.or.N_param==7.or.N_param==8.or.N_param==9)THEN
         IF(c_param==4.0)THEN
            c0=-0.198
            c1= 0.575
            c2=-0.00550
            c3=-0.000223
         ELSE IF(c_param==10.0)THEN
            c0= 0.341
            c1= 1.30
            c2=-0.0136
            c3=-0.00235
         END IF
         DO i=1,n
            up(i) = c0+c1*y1(i)+c2*y1(i)*y1(i)+c3*y1(i)*y1(i)*y1(i)
         END DO
      END IF
         
      y2(1)=y1(n)*(y1(2)-y1(n-1))-y1(1)+f_param-up(1)+df
      y2(2)=y1(1)*(y1(3)-y1(n))  -y1(2)+f_param-up(2)+df
      DO i=3,n-1
         y2(i)=y1(i-1)*(y1(i+1)-y1(i-2))-y1(i)+f_param-up(i)+df
      END DO
      y2(n)=y1(n-1)*(y1(1)-y1(n-2))-y1(n)+f_param-up(n)+df

    END SUBROUTINE discrete_L95


    SUBROUTINE step_full_L95(deltat,y,Out,y_small,Out_small,f_param,c_param,df)
    
      IMPLICIT NONE

      REAL,dimension(:) :: y,Out
      REAL,dimension(:,:) :: y_small,Out_small
      REAL,dimension(:) :: dy(size(y(:)))
      REAL,dimension(:,:) :: dy_small(size(y_small(:,1)),size(y_small(1,:)))
      REAL,dimension(:,:):: work(size(y(:)),3)
      REAL,dimension(:,:,:) :: work_small(size(y_small(:,1)),size(y_small(1,:)),3)
      REAL deltat, deltat_2, deltat_6,f_param,c_param,df
      INTEGER :: i,j

      deltat_2 = 0.5*deltat
      deltat_6 = deltat/6.0

      !RK4 timestepping:
      CALL discrete_L95_full(y,dy,y_small,dy_small,f_param,c_param,df)
   
      DO i = 1,size(y(:))
         work(i,1)= y(i) + deltat_2 * dy(i)
      END DO
      DO i = 1,size(y_small(:,1))
         DO j = 1,size(y_small(1,:))
            work_small(i,j,1)= y_small(i,j) + deltat_2 * dy_small(i,j)
         END DO
      END DO

      CALL discrete_L95_full(work(:,1),work(:,2),work_small(:,:,1),work_small(:,:,2),f_param,c_param,df)

      DO i = 1,size(y(:))
         work(i,1)= y(i) + deltat_2 * work(i,2)
      END DO
      DO i = 1,size(y_small(:,1))
         DO j = 1,size(y_small(1,:))
            work_small(i,j,1)= y_small(i,j) + deltat_2 * work_small(i,j,2)
         END DO
      END DO

      CALL discrete_L95_full(work(:,1),work(:,3),work_small(:,:,1),work_small(:,:,3),f_param,c_param,df)

      DO i = 1,size(y(:))
         work(i,1)= y(i) + deltat * work(i,3)
         work(i,3)=work(i,3)+work(i,2)
      END DO
      DO i = 1,size(y_small(:,1))
         DO j = 1,size(y_small(1,:))
            work_small(i,j,1)= y_small(i,j) + deltat * work_small(i,j,3)
            work_small(i,j,3)=work_small(i,j,3)+work_small(i,j,2)
         END DO
      END DO

      CALL discrete_L95_full(work(:,1),work(:,2),work_small(:,:,1),work_small(:,:,2),f_param,c_param,df)

      DO i = 1,size(y(:))
         Out(i)= y(i) + deltat_6 * (dy(i)+work(i,2)+2.0*work(i,3))
      END DO

      DO i = 1,size(y_small(:,1))
         DO j = 1,size(y_small(1,:))
            Out_small(i,j)= y_small(i,j) + deltat_6 * (dy_small(i,j)+work_small(i,j,2)+2.0*work_small(i,j,3))
         END DO
      END DO

    END SUBROUTINE step_full_L95

    SUBROUTINE discrete_L95_full(y1,y2,y1_small,y2_small,f_param,c_param,df)
    
      IMPLICIT NONE

      REAL,dimension(:) :: y1,y2
      REAL,dimension(:,:) :: y1_small,y2_small
      REAL :: B(size(y1(:)))
      REAL :: f_param,c_param,df
      INTEGER :: i,j,im,ip,ipp
      INTEGER :: n
      
      n = size(y1(:))

      DO i=1,n
         B(i) = 0.0
         DO j=1,32
            B(i) = B(i)+c_param/10.0 * y1_small(i,j)
         END DO
      END DO
  
      y2(1)=y1(n)*(y1(2)-y1(n-1))-y1(1)+f_param+df-B(1)
      y2(2)=y1(1)*(y1(3)-y1(n))  -y1(2)+f_param+df-B(2)
      DO i=3,n-1
         y2(i)=y1(i-1)*(y1(i+1)-y1(i-2))-y1(i)+f_param+df-B(i)
      END DO
      y2(n)=y1(n-1)*(y1(1)-y1(n-2))-y1(n)+f_param+df-B(n)

      DO i=1,n
         im=i-1
         ip=i+1
         ipp=i+2
         IF(i==1) im=n
         IF(i==n) ip=1
         
         y2_small(i,1)= -c_param*10.0*y1_small(i,2)*(y1_small(i,3)-y1_small(im,32)) &
              &-c_param*y1_small(i,1)+c_param/10.0*y1(i)
         DO j=2,32-2
            y2_small(i,j)= -c_param*10.0*y1_small(i,j+1)*(y1_small(i,j+2)-y1_small(i,j-1)) &
                 & -c_param*y1_small(i,j)+c_param/10.0*y1(i)
         END DO
         y2_small(i,31)= -c_param*10.0*y1_small(i,32)*(y1_small(ip,1)-y1_small(i,31-1)) &
              & -c_param*y1_small(i,31)+c_param/10.0*y1(i)
         y2_small(i,32)= -c_param*10.0*y1_small(ip,1)*(y1_small(ip,2)-y1_small(i,32-1)) &
              & -c_param*y1_small(i,32)+c_param/10.0*y1(i)
      END DO

    END SUBROUTINE discrete_L95_full



    SUBROUTINE ini_para(N_param,noise_amp,noise_cor,sigma0,sigma1,sigmaM,sigmaA,b0,b1,b2,b3)

      INTEGER :: N_param
      REAL :: noise_amp,noise_cor,sigma0,sigma1,sigmaM,sigmaA,b0,b1,b2,b3
      
      IF(N_param==0.or.N_param==10)THEN
         noise_amp = 0.0
         write(*,*) 'No parametrisation.'
      ELSE IF(N_param==1)THEN
         IF(f_param==18.0)THEN
            noise_amp = 1.74
         ELSE IF(f_param==20.0)THEN
            noise_amp = 1.99
         END IF
         noise_cor = 0.986
         write(*,*) 'Wilks 2005 deterministic'
      ELSE IF(N_param==2)THEN
         IF(f_param==18.0)THEN
            noise_amp = 1.74
         ELSE IF(f_param==20.0)THEN
            noise_amp = 1.99
         END IF
         noise_cor = 0.986
         write(*,*) 'Wilks 2005 AR0'
      ELSE IF(N_param==3)THEN
         IF(f_param==18.0)THEN
            noise_amp = 1.74
         ELSE IF(f_param==20.0)THEN
            noise_amp = 1.99
         END IF
         noise_cor = 0.986
         write(*,*) 'Wilks 2005 AR1'
      ELSE IF(N_param==4)THEN
         IF(c_param==10.0)THEN
            noise_amp = 1.99
            noise_cor = 0.986
         ELSE IF(c_param==4.0)THEN
            noise_amp = 2.12
            noise_cor = 0.993
         END IF
         write(*,*) 'Arnold et al. 2012 deterministic'
      ELSE IF(N_param==5)THEN
         IF(c_param==10.0)THEN
            noise_amp = 1.99
            noise_cor = 0.986
         ELSE IF(c_param==4.0)THEN
            noise_amp = 2.12
            noise_cor = 0.993
         END IF
         write(*,*) 'Arnold et al. 2012 AR0'
      ELSE IF(N_param==6)THEN
         IF(c_param==10.0)THEN
            noise_amp = 1.99
            noise_cor = 0.986
         ELSE IF(c_param==4.0)THEN
            noise_amp = 2.12
            noise_cor = 0.993
         END IF
         write(*,*) 'Arnold et al. 2012 AR1'
      ELSE IF(N_param==7)THEN
         IF(c_param==10.0)THEN
            sigma0 = 1.47
            sigma1 = 0.0873
            noise_cor = 0.986
         ELSE IF(c_param==4.0)THEN
            sigma0 = 1.62
            sigma1 = 0.0780
            noise_cor = 0.993
         END IF
         write(*,*) 'Arnold et al. 2012 state dependent'
      ELSE IF(N_param==8)THEN
         IF(c_param==10.0)THEN
            noise_amp = 0.469
            noise_cor = 0.94
            b0 = 0.341
            b1 = 1.3
            b2 = -0.0136
            b3 = -0.00235
         ELSE IF(c_param==4.0)THEN
            noise_amp = 0.746
            noise_cor = 0.95
            b0 = -0.198
            b1 = 0.575
            b2 = -0.0055
            b3 = 0.000223
         END IF
         write(*,*) 'Arnold et al. 2012 multiplicative'
      ELSE IF(N_param==9)THEN
         IF(c_param==10.0)THEN
            sigmaM = 0.101
            sigmaA = 1.37
            noise_cor = 0.988
            b0 = 0.341
            b1 = 1.3
            b2 = -0.0136
            b3 = -0.00235
         ELSE IF(c_param==4.0)THEN
            sigmaM = 0.177
            sigmaA = 1.55
            noise_cor = 0.993
            b0 = -0.198
            b1 = 0.575
            b2 = -0.0055
            b3 = -0.000223
         END IF
         write(*,*) 'Arnold et al. 2012 multiplicative additive'
      ELSE 
         write(*,*) "N_param is set wrong!", N_param
         STOP
      END IF

    END SUBROUTINE INI_PARA
    
    
    FUNCTION rand_normal(mean,stdev) RESULT(c)
      IMPLICIT NONE
      
      REAL :: mean,stdev,c,rand(2),dummy1,dummy2
      IF(stdev <= 0.0d0) THEN
         WRITE(*,*) "Standard Deviation must be .gt. zero."
      ELSE
         CALL RANDOM_NUMBER(rand)
         dummy1=(-2.0d0*log(rand(1)))**0.5
         dummy2 = 2.0d0*3.141592654*rand(2)
         c= mean+stdev*dummy1*sin(dummy2)
      END IF
    END FUNCTION rand_normal
    
END PROGRAM LORENZ95

