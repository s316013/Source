
    program Source

    implicit none
    
   PARAMETER (ntens=6, nstatv=125, nprops=155,)
   
   dimension stress(ntens),statev(nstatv),&
     ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),&
     stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),&
     props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)
   
   

   
     CALL UMAT(stress,statev,ddsdde,sse,spd,scd,
     1 rpl, ddsddt, drplde, drpldt,
     2 stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     3 ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     4 celent,dfgrd0,dfgrd1,noel,npt,kspt)   
       
       
   open(unit=125, file='C:\Temp\New Framework\ GAMMA(I).txt')
                   Write (125,*) STATEV(NSLPTL+I)
    
    end program Source    
    
    
      
      SUBROUTINE UMAT(stress,statev,ddsdde,sse,spd,scd,
     1 rpl, ddsddt, drplde, drpldt,
     2 stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     3 ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     4 celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)


      PARAMETER (ND=150)

      include 'aba_param.inc'
!
      CHARACTER*8 CMNAME
      EXTERNAL F

      dimension stress(ntens),statev(nstatv),
     1 ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     2 stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     3 props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)

      DIMENSION ISPDIR(3), ISPNOR(3), NSLIP(3), 
     2          SLPDIR(3,ND), SLPNOR(3,ND), SLPDEF(6,ND), 
     3          SLPSPN(3,ND), DSPDIR(3,ND), DSPNOR(3,ND), 
     4          DLOCAL(6,6), D(6,6), ROTD(6,6), ROTATE(3,3), 
     5          FSLIP(ND), DFDXSP(ND), DDEMSD(6,ND), 
     6          H(ND,ND), DDGDDE(ND,6), 
     7          DSTRES(6), DELATS(6), DSPIN(3), DVGRAD(3,3),
     8          DGAMMA(ND), DTAUSP(ND), DGSLIP(ND), 
     9          WORKST(ND,ND), INDX(ND), TERM(3,3), TRM0(3,3), ITRM(3)

      DIMENSION FSLIP1(ND), STRES1(6), GAMMA1(ND), TAUSP1(ND), 
     2          GSLP1(ND), SPNOR1(3,ND), SPDIR1(3,ND), DDSDE1(6,6),
     3          DSOLD(6), DGAMOD(ND), DTAUOD(ND), DGSPOD(ND), 
     4          DSPNRO(3,ND), DSPDRO(3,ND), 
     5          DHDGDG(ND,ND)
      

!-----  Elastic matrix in local cubic crystal system: DLOCAL
      DO J=1,6
         DO I=1,6
            DLOCAL(I,J)=0.
         END DO
      END DO

      CHECK=0.
      DO J=10,21
         CHECK=CHECK+ABS(PROPS(J))
      END DO

      IF (CHECK.EQ.0.) THEN
         DO J=4,9
            CHECK=CHECK+ABS(PROPS(J))
         END DO

         IF (CHECK.EQ.0.) THEN

            IF (PROPS(3).EQ.0.) THEN

!-----  Isotropic material
               GSHEAR=PROPS(1)/2./(1.+PROPS(2))
               E11=2.*GSHEAR*(1.-PROPS(2))/(1.-2.*PROPS(2))
               E12=2.*GSHEAR*PROPS(2)/(1.-2.*PROPS(2))

               DO J=1,3
                  DLOCAL(J,J)=E11

                  DO I=1,3
                     IF (I.NE.J) DLOCAL(I,J)=E12
                  END DO

                  DLOCAL(J+3,J+3)=GSHEAR
               END DO

            ELSE

!-----  Cubic material
               DO J=1,3
                  DLOCAL(J,J)=PROPS(1)

                  DO I=1,3
                     IF (I.NE.J) DLOCAL(I,J)=PROPS(2)
                  END DO

                  DLOCAL(J+3,J+3)=PROPS(3)
               END DO

            END IF

         ELSE

!-----  Orthotropic metarial
            DLOCAL(1,1)=PROPS(1)
            DLOCAL(1,2)=PROPS(2)
            DLOCAL(2,1)=PROPS(2)
            DLOCAL(2,2)=PROPS(3)

            DLOCAL(1,3)=PROPS(4)
            DLOCAL(3,1)=PROPS(4)
            DLOCAL(2,3)=PROPS(5)
            DLOCAL(3,2)=PROPS(5)
            DLOCAL(3,3)=PROPS(6)

            DLOCAL(4,4)=PROPS(7)
            DLOCAL(5,5)=PROPS(8)
            DLOCAL(6,6)=PROPS(9)

         END IF

      ELSE

!-----  General anisotropic material
         ID=0
         DO J=1,6
            DO I=1,J
               ID=ID+1
               DLOCAL(I,J)=PROPS(ID)
               DLOCAL(J,I)=DLOCAL(I,J)
            END DO
         END DO
      END IF

!!-----  Rotation matrix: ROTATE, i.e. direction cosines of [100], [010]
!     and [001] of a cubic crystal in global system
!
      CALL ROTATION (PROPS(57), ROTATE)
      

!-----  Rotation matrix: ROTD to transform local elastic matrix DLOCAL 
!     to global elastic matrix D
!
      DO J=1,3
         J1=1+J/3
         J2=2+J/2

         DO I=1,3
            I1=1+I/3
            I2=2+I/2

            ROTD(I,J)=ROTATE(I,J)**2
            ROTD(I,J+3)=2.*ROTATE(I,J1)*ROTATE(I,J2)
            ROTD(I+3,J)=ROTATE(I1,J)*ROTATE(I2,J)
            ROTD(I+3,J+3)=ROTATE(I1,J1)*ROTATE(I2,J2)+
     2                    ROTATE(I1,J2)*ROTATE(I2,J1)

         END DO
      END DO

!-----  Elastic matrix in global system: D
!     {D} = {ROTD} * {DLOCAL} * {ROTD}transpose
!
      DO J=1,6
         DO I=1,6
            D(I,J)=0.
         END DO
      END DO

      DO J=1,6
         DO I=1,J

            DO K=1,6
               DO L=1,6
                  D(I,J)=D(I,J)+DLOCAL(K,L)*ROTD(I,K)*ROTD(J,L)
               END DO
            END DO

            D(J,I)=D(I,J)

         END DO
      END DO

!-----  Total number of sets of slip systems: NSET
      NSET=NINT(PROPS(25))
      IF (NSET.LT.1) THEN
         WRITE (6,*) '***ERROR - zero sets of slip systems'
         STOP
      ELSE IF (NSET.GT.3) THEN
         WRITE (6,*) 
     2     '***ERROR - more than three sets of slip systems'
         STOP
      END IF

!-----  Implicit integration parameter: THETA
      THETA=PROPS(145)

!-----  Finite deformation ?
!-----  NLGEOM = 0,   small deformation theory
!       otherwise, theory of finite rotation and finite strain, Users 
!     must declare "NLGEOM" in the input file, at the *STEP card
!
      IF (PROPS(146).EQ.0.) THEN
         NLGEOM=0
      ELSE
         NLGEOM=1
      END IF

!-----  Iteration?
!-----  ITRATN = 0, no iteration
!       otherwise, iteration (solving increments of stresses and 
!     solution dependent state variables)
!
      IF (PROPS(153).EQ.0.) THEN
         ITRATN=0
      ELSE
         ITRATN=1
      END IF

      ITRMAX=NINT(PROPS(154))
      GAMERR=PROPS(155)

      NITRTN=-1

      DO I=1,NTENS
         DSOLD(I)=0.
      END DO

      DO J=1,ND
         DGAMOD(J)=0.
         DTAUOD(J)=0.
         DGSPOD(J)=0.
         DO I=1,3
            DSPNRO(I,J)=0.
            DSPDRO(I,J)=0.
         END DO
      END DO

!-----  Increment of spin associated with the material element: DSPIN
!     (only needed for finite rotation)
!
      IF (NLGEOM.NE.0) THEN
         DO J=1,3
            DO I=1,3
               TERM(I,J)=DROT(J,I)
               TRM0(I,J)=DROT(J,I)
            END DO

            TERM(J,J)=TERM(J,J)+1.D0
            TRM0(J,J)=TRM0(J,J)-1.D0
         END DO

         CALL LUDCMP (TERM, 3, 3, ITRM, DDCMP)

         DO J=1,3
            CALL LUBKSB (TERM, 3, 3, ITRM, TRM0(1,J))
         END DO

         DSPIN(1)=TRM0(2,1)-TRM0(1,2)
         DSPIN(2)=TRM0(1,3)-TRM0(3,1)
         DSPIN(3)=TRM0(3,2)-TRM0(2,3)

      END IF

!-----  Increment of dilatational strain: DEV
      DEV=0.D0
      DO I=1,NDI
         DEV=DEV+DSTRAN(I)
      END DO

!-----  Iteration starts (only when iteration method is used)
1000  CONTINUE

!-----  Parameter NITRTN: number of iterations
!       NITRTN = 0 --- no-iteration solution
!
      NITRTN=NITRTN+1

!-----  Check whether the current stress state is the initial state
      IF (STATEV(1).EQ.0.) THEN

!-----  Initial state
!
!-----  Generating the following parameters and variables at initial 
!     state:
!          Total number of slip systems in all the sets NSLPTL
!          Number of slip systems in each set NSLIP
!          Unit vectors in initial slip directions SLPDIR
!          Unit normals to initial slip planes SLPNOR
!
         NSLPTL=0
         DO I=1,NSET
            ISPNOR(1)=NINT(PROPS(25+8*I))
            ISPNOR(2)=NINT(PROPS(26+8*I))
            ISPNOR(3)=NINT(PROPS(27+8*I))

            ISPDIR(1)=NINT(PROPS(28+8*I))
            ISPDIR(2)=NINT(PROPS(29+8*I))
            ISPDIR(3)=NINT(PROPS(30+8*I))

            CALL SLIPSYS (ISPDIR, ISPNOR, NSLIP(I), SLPDIR(1,NSLPTL+1), 
     2                    SLPNOR(1,NSLPTL+1), ROTATE)

            NSLPTL=NSLPTL+NSLIP(I)
         END DO
         

         IF (ND.LT.NSLPTL) THEN
            WRITE (6,*) 
     2 '***ERROR - parameter ND chosen by the present user is less than
     3             the total number of slip systems NSLPTL'
            STOP
         END IF

!-----  Slip deformation tensor: SLPDEF (Schmid factors)
         DO J=1,NSLPTL
            SLPDEF(1,J)=SLPDIR(1,J)*SLPNOR(1,J)
            SLPDEF(2,J)=SLPDIR(2,J)*SLPNOR(2,J)
            SLPDEF(3,J)=SLPDIR(3,J)*SLPNOR(3,J)
            SLPDEF(4,J)=SLPDIR(1,J)*SLPNOR(2,J)+SLPDIR(2,J)*SLPNOR(1,J)
            SLPDEF(5,J)=SLPDIR(1,J)*SLPNOR(3,J)+SLPDIR(3,J)*SLPNOR(1,J)
            SLPDEF(6,J)=SLPDIR(2,J)*SLPNOR(3,J)+SLPDIR(3,J)*SLPNOR(2,J)
         END DO

!-----  Initial value of state variables: unit normal to a slip plane 
!     and unit vector in a slip direction
!
         STATEV(NSTATV)=FLOAT(NSLPTL)
         
         DO I=1,NSET
            STATEV(NSTATV-4+I)=FLOAT(NSLIP(I))
         END DO

         IDNOR=3*NSLPTL
         IDDIR=6*NSLPTL
         DO J=1,NSLPTL
            DO I=1,3
               IDNOR=IDNOR+1
               STATEV(IDNOR)=SLPNOR(I,J)

               IDDIR=IDDIR+1
               STATEV(IDDIR)=SLPDIR(I,J)
            END DO
         END DO

!-----  Initial value of the current strength for all slip systems
!
         CALL GSLPINIT (STATEV(1), NSLIP, NSLPTL, NSET, PROPS(97))

!-----  Initial value of shear strain in slip systems
!FIX--  Initial value of cumulative shear strain in each slip systems

         DO I=1,NSLPTL
            STATEV(NSLPTL+I)=0.
!FIXA
            STATEV(9*NSLPTL+I)=0.
!FIXB
         END DO

!FIXA
         STATEV(10*NSLPTL+1)=0.
!FIXB

!-----  Initial value of the resolved shear stress in slip systems
              
         DO I=1,NSLPTL
            TERM1=0.

            DO J=1,NTENS
               IF (J.LE.NDI) THEN
                  TERM1=TERM1+SLPDEF(J,I)*STRESS(J)
               ELSE
                  TERM1=TERM1+SLPDEF(J-NDI+3,I)*STRESS(J)
               END IF
            END DO

            STATEV(2*NSLPTL+I)=TERM1
         END DO
             

      ELSE

!-----  Current stress state
!
C-----  Copying from the array of state variables STATVE the following
C          parameters and variables at current stress state:
C          Total number of slip systems in all the sets NSLPTL
!          Number of slip systems in each set NSLIP
!          Current slip directions SLPDIR
!          Normals to current slip planes SLPNOR
C
         NSLPTL=NINT(STATEV(NSTATV))
         DO I=1,NSET
            NSLIP(I)=NINT(STATEV(NSTATV-4+I))
         END DO

         IDNOR=3*NSLPTL
         IDDIR=6*NSLPTL
         DO J=1,NSLPTL
            DO I=1,3
               IDNOR=IDNOR+1
               SLPNOR(I,J)=STATEV(IDNOR)

               IDDIR=IDDIR+1
               SLPDIR(I,J)=STATEV(IDDIR)
            END DO
         END DO

C-----  Slip deformation tensor: SLPDEF (Schmid factors)
         DO J=1,NSLPTL
            SLPDEF(1,J)=SLPDIR(1,J)*SLPNOR(1,J)
            SLPDEF(2,J)=SLPDIR(2,J)*SLPNOR(2,J)
            SLPDEF(3,J)=SLPDIR(3,J)*SLPNOR(3,J)
            SLPDEF(4,J)=SLPDIR(1,J)*SLPNOR(2,J)+SLPDIR(2,J)*SLPNOR(1,J)
            SLPDEF(5,J)=SLPDIR(1,J)*SLPNOR(3,J)+SLPDIR(3,J)*SLPNOR(1,J)
            SLPDEF(6,J)=SLPDIR(2,J)*SLPNOR(3,J)+SLPDIR(3,J)*SLPNOR(2,J)
         END DO

      END IF
      

!-----  Slip spin tensor: SLPSPN (only needed for finite rotation)
      IF (NLGEOM.NE.0) THEN
         DO J=1,NSLPTL
            SLPSPN(1,J)=0.5*(SLPDIR(1,J)*SLPNOR(2,J)-
     2                       SLPDIR(2,J)*SLPNOR(1,J))
            SLPSPN(2,J)=0.5*(SLPDIR(3,J)*SLPNOR(1,J)-
     2                       SLPDIR(1,J)*SLPNOR(3,J))
            SLPSPN(3,J)=0.5*(SLPDIR(2,J)*SLPNOR(3,J)-
     2                       SLPDIR(3,J)*SLPNOR(2,J))
         END DO
      END IF

!-----  Double dot product of elastic moduli tensor with the slip 
!     deformation tensor (Schmid factors) plus, only for finite 
!     rotation, the dot product of slip spin tensor with the stress: 
!     DDEMSD
!
      DO J=1,NSLPTL
         DO I=1,6
            DDEMSD(I,J)=0.
            DO K=1,6
               DDEMSD(I,J)=DDEMSD(I,J)+D(K,I)*SLPDEF(K,J)
            END DO
         END DO
      END DO

      IF (NLGEOM.NE.0) THEN
         DO J=1,NSLPTL

            DDEMSD(4,J)=DDEMSD(4,J)-SLPSPN(1,J)*STRESS(1)
            DDEMSD(5,J)=DDEMSD(5,J)+SLPSPN(2,J)*STRESS(1)

            IF (NDI.GT.1) THEN
               DDEMSD(4,J)=DDEMSD(4,J)+SLPSPN(1,J)*STRESS(2)
               DDEMSD(6,J)=DDEMSD(6,J)-SLPSPN(3,J)*STRESS(2)
            END IF

            IF (NDI.GT.2) THEN
               DDEMSD(5,J)=DDEMSD(5,J)-SLPSPN(2,J)*STRESS(3)
               DDEMSD(6,J)=DDEMSD(6,J)+SLPSPN(3,J)*STRESS(3)
            END IF

            IF (NSHR.GE.1) THEN
               DDEMSD(1,J)=DDEMSD(1,J)+SLPSPN(1,J)*STRESS(NDI+1)
               DDEMSD(2,J)=DDEMSD(2,J)-SLPSPN(1,J)*STRESS(NDI+1)
               DDEMSD(5,J)=DDEMSD(5,J)-SLPSPN(3,J)*STRESS(NDI+1)
               DDEMSD(6,J)=DDEMSD(6,J)+SLPSPN(2,J)*STRESS(NDI+1)
            END IF

            IF (NSHR.GE.2) THEN
               DDEMSD(1,J)=DDEMSD(1,J)-SLPSPN(2,J)*STRESS(NDI+2)
               DDEMSD(3,J)=DDEMSD(3,J)+SLPSPN(2,J)*STRESS(NDI+2)
               DDEMSD(4,J)=DDEMSD(4,J)+SLPSPN(3,J)*STRESS(NDI+2)
               DDEMSD(6,J)=DDEMSD(6,J)-SLPSPN(1,J)*STRESS(NDI+2)
            END IF

            IF (NSHR.EQ.3) THEN
               DDEMSD(2,J)=DDEMSD(2,J)+SLPSPN(3,J)*STRESS(NDI+3)
               DDEMSD(3,J)=DDEMSD(3,J)-SLPSPN(3,J)*STRESS(NDI+3)
               DDEMSD(4,J)=DDEMSD(4,J)-SLPSPN(2,J)*STRESS(NDI+3)
               DDEMSD(5,J)=DDEMSD(5,J)+SLPSPN(1,J)*STRESS(NDI+3)
            END IF

         END DO
      END IF
      

!-----  Shear strain-rate in a slip system at the start of in!rement: 
!     FSLIP, and its derivative: DFDXSP
!
      ID=1
      DO I=1,NSET
         IF (I.GT.1) ID=ID+NSLIP(I-1)
         CALL STRAINRATE (STATEV(NSLPTL+ID), STATEV(2*NSLPTL+ID), 
     2                    STATEV(ID), NSLIP(I), FSLIP(ID), DFDXSP(ID), 
     3                    PROPS(65+8*I))     
      END DO

!-----  Self- and latent-hardening laws
!FIXA  
       CALL LATENTHARDEN (STATEV(NSLPTL+1), STATEV(2*NSLPTL+1), 
     2                   STATEV(1), STATEV(9*NSLPTL+1),
     3                   STATEV(10*NSLPTL+1), NSLIP, NSLPTL, 
     4                   NSET, H(1,1), PROPS(97), ND)
       
!FIXB

!-----  LU de!omposition to solve the in!rement of shear strain in a 
!     slip system
!
      TERM1=THETA*DTIME
      DO I=1,NSLPTL
         TAUSLP=STATEV(2*NSLPTL+I)
         GSLIP=STATEV(I)
         X=TAUSLP/GSLIP
         TERM2=TERM1*DFDXSP(I)/GSLIP
         TERM3=TERM1*X*DFDXSP(I)/GSLIP

         DO J=1,NSLPTL
            TERM4=0.
            DO K=1,6
               TERM4=TERM4+DDEMSD(K,I)*SLPDEF(K,J)
            END DO

            WORKST(I,J)=TERM2*TERM4+H(I,J)*TERM3*DSIGN(1.D0,FSLIP(J))

            IF (NITRTN.GT.0) WORKST(I,J)=WORKST(I,J)+TERM3*DHDGDG(I,J)

         END DO

         WORKST(I,I)=WORKST(I,I)+1.
      END DO

      CALL LUDCMP (WORKST, NSLPTL, ND, INDX, DDCMP)

!-----  Increment of shear strain in a slip system: DGAMMA
      TERM1=THETA*DTIME
      DO I=1,NSLPTL

         IF (NITRTN.EQ.0) THEN
            TAUSLP=STATEV(2*NSLPTL+I)
            GSLIP=STATEV(I)
            X=TAUSLP/GSLIP
            TERM2=TERM1*DFDXSP(I)/GSLIP

            DGAMMA(I)=0.
            DO J=1,NDI
               DGAMMA(I)=DGAMMA(I)+DDEMSD(J,I)*DSTRAN(J)
            END DO

            IF (NSHR.GT.0) THEN
               DO J=1,NSHR
                  DGAMMA(I)=DGAMMA(I)+DDEMSD(J+3,I)*DSTRAN(J+NDI)
               END DO
            END IF

            DGAMMA(I)=DGAMMA(I)*TERM2+FSLIP(I)*DTIME
         ELSE
            DGAMMA(I)=TERM1*(FSLIP(I)-FSLIP1(I))+FSLIP1(I)*DTIME
     2                -DGAMOD(I)

         END IF

      END DO

      CALL LUBKSB (WORKST, NSLPTL, ND, INDX, DGAMMA)

      DO I=1,NSLPTL
         DGAMMA(I)=DGAMMA(I)+DGAMOD(I)
      END DO

!-----  Update the shear strain in a slip system: STATEV(NSLPTL+1) - 
!     STATEV(2*NSLPTL)
!
      DO I=1,NSLPTL
         STATEV(NSLPTL+I)=STATEV(NSLPTL+I)+DGAMMA(I)-DGAMOD(I)
      END DO

!-----  Increment of current strength in a slip system: DGSLIP
      DO I=1,NSLPTL
         DGSLIP(I)=0.
         DO J=1,NSLPTL
            DGSLIP(I)=DGSLIP(I)+H(I,J)*ABS(DGAMMA(J))
         END DO
      END DO

!-----  Update the current strength in a slip system: STATEV(1) - 
!     STATEV(NSLPTL)
!
      DO I=1,NSLPTL
         STATEV(I)=STATEV(I)+DGSLIP(I)-DGSPOD(I)
      END DO

!-----  Increment of strain associated with lattice stretching: DELATS
      DO J=1,6
         DELATS(J)=0.
      END DO

      DO J=1,3
         IF (J.LE.NDI) DELATS(J)=DSTRAN(J)
         DO I=1,NSLPTL
            DELATS(J)=DELATS(J)-SLPDEF(J,I)*DGAMMA(I)
         END DO
      END DO

      DO J=1,3
         IF (J.LE.NSHR) DELATS(J+3)=DSTRAN(J+NDI)
         DO I=1,NSLPTL
            DELATS(J+3)=DELATS(J+3)-SLPDEF(J+3,I)*DGAMMA(I)
         END DO
      END DO

!-----  In!rement of deformation gradient asso!iated with latti!e 
!     stret!hing in the !urrent state, i.e. the velo!ity gradient 
!     (asso!iated with latti!e stret!hing) times the in!rement of time:
!     DVGRAD (only needed for finite rotation)
!
      IF (NLGEOM.NE.0) THEN
         DO J=1,3
            DO I=1,3
               IF (I.EQ.J) THEN
                  DVGRAD(I,J)=DELATS(I)
               ELSE
                  DVGRAD(I,J)=DELATS(I+J+1)
               END IF
            END DO
         END DO

         DO J=1,3
            DO I=1,J
               IF (J.GT.I) THEN
                  IJ2=I+J-2
                  IF (MOD(IJ2,2).EQ.1) THEN
                     TERM1=1.
                  ELSE
                     TERM1=-1.
                  END IF

                  DVGRAD(I,J)=DVGRAD(I,J)+TERM1*DSPIN(IJ2)
                  DVGRAD(J,I)=DVGRAD(J,I)-TERM1*DSPIN(IJ2)

                  DO K=1,NSLPTL
                     DVGRAD(I,J)=DVGRAD(I,J)-TERM1*DGAMMA(K)*
     2                                       SLPSPN(IJ2,K)
                     DVGRAD(J,I)=DVGRAD(J,I)+TERM1*DGAMMA(K)*
     2                                       SLPSPN(IJ2,K)
                  END DO
               END IF

            END DO
         END DO

      END IF

!-----  In!rement of resolved shear stress in a slip system: DTAUSP
      DO I=1,NSLPTL
         DTAUSP(I)=0.
         DO J=1,6
            DTAUSP(I)=DTAUSP(I)+DDEMSD(J,I)*DELATS(J)
         END DO
      END DO

!-----  Update the resolved shear stress in a slip system: 
!     STATEV(2*NSLPTL+1) - STATEV(3*NSLPTL)
!
      DO I=1,NSLPTL
         STATEV(2*NSLPTL+I)=STATEV(2*NSLPTL+I)+DTAUSP(I)-DTAUOD(I)
      END DO

!-----  In!rement of stress: DSTRES
      IF (NLGEOM.EQ.0) THEN
         DO I=1,NTENS
            DSTRES(I)=0.
         END DO
      ELSE
         DO I=1,NTENS
            DSTRES(I)=-STRESS(I)*DEV
         END DO
      END IF

      DO I=1,NDI
         DO J=1,NDI
            DSTRES(I)=DSTRES(I)+D(I,J)*DSTRAN(J)
         END DO

         IF (NSHR.GT.0) THEN
            DO J=1,NSHR
               DSTRES(I)=DSTRES(I)+D(I,J+3)*DSTRAN(J+NDI)
            END DO
         END IF

         DO J=1,NSLPTL
            DSTRES(I)=DSTRES(I)-DDEMSD(I,J)*DGAMMA(J)
         END DO
      END DO

      IF (NSHR.GT.0) THEN
         DO I=1,NSHR

            DO J=1,NDI
               DSTRES(I+NDI)=DSTRES(I+NDI)+D(I+3,J)*DSTRAN(J)
            END DO

            DO J=1,NSHR
               DSTRES(I+NDI)=DSTRES(I+NDI)+D(I+3,J+3)*DSTRAN(J+NDI)
            END DO

            DO J=1,NSLPTL
               DSTRES(I+NDI)=DSTRES(I+NDI)-DDEMSD(I+3,J)*DGAMMA(J)
            END DO

         END DO
      END IF

!-----  Update the stress: STRESS
      DO I=1,NTENS
         STRESS(I)=STRESS(I)+DSTRES(I)-DSOLD(I)
      END DO

!-----  In!rement of normal to a slip plane and a slip dire!tion (only 
!     needed for finite rotation)
!
      IF (NLGEOM.NE.0) THEN
         DO J=1,NSLPTL
            DO I=1,3
               DSPNOR(I,J)=0.
               DSPDIR(I,J)=0.

               DO K=1,3
                  DSPNOR(I,J)=DSPNOR(I,J)-SLPNOR(K,J)*DVGRAD(K,I)
                  DSPDIR(I,J)=DSPDIR(I,J)+SLPDIR(K,J)*DVGRAD(I,K)
               END DO

            END DO
         END DO

!-----  Update the normal to a slip plane and a slip dire!tion (only 
!     needed for finite rotation)
!
         IDNOR=3*NSLPTL
         IDDIR=6*NSLPTL
         DO J=1,NSLPTL
            DO I=1,3
               IDNOR=IDNOR+1
               STATEV(IDNOR)=STATEV(IDNOR)+DSPNOR(I,J)-DSPNRO(I,J)

               IDDIR=IDDIR+1
               STATEV(IDDIR)=STATEV(IDDIR)+DSPDIR(I,J)-DSPDRO(I,J)
            END DO
         END DO

      END IF

!-----  Derivative of shear strain in!rement in a slip system w.r.t. 
!     strain in!rement: DDGDDE
!
      TERM1=THETA*DTIME
      DO I=1,NTENS
         DO J=1,NSLPTL
            TAUSLP=STATEV(2*NSLPTL+J)
            GSLIP=STATEV(J)
            X=TAUSLP/GSLIP
            TERM2=TERM1*DFDXSP(J)/GSLIP
            IF (I.LE.NDI) THEN
               DDGDDE(J,I)=TERM2*DDEMSD(I,J)
            ELSE
               DDGDDE(J,I)=TERM2*DDEMSD(I-NDI+3,J)
            END IF
         END DO

         CALL LUBKSB (WORKST, NSLPTL, ND, INDX, DDGDDE(1,I))

      END DO

!-----  Derivative of stress in!rement w.r.t. strain in!rement, i.e. 
!     Ja!obian matrix
!
!-----  Jacobian matrix: elastic part
      DO J=1,NTENS
         DO I=1,NTENS
            DDSDDE(I,J)=0.
         END DO
      END DO

      DO J=1,NDI
         DO I=1,NDI
            DDSDDE(I,J)=D(I,J)
            IF (NLGEOM.NE.0) DDSDDE(I,J)=DDSDDE(I,J)-STRESS(I)
         END DO
      END DO

      IF (NSHR.GT.0) THEN
         DO J=1,NSHR
            DO I=1,NSHR
               DDSDDE(I+NDI,J+NDI)=D(I+3,J+3)
            END DO

            DO I=1,NDI
               DDSDDE(I,J+NDI)=D(I,J+3)
               DDSDDE(J+NDI,I)=D(J+3,I)
               IF (NLGEOM.NE.0)
     2            DDSDDE(J+NDI,I)=DDSDDE(J+NDI,I)-STRESS(J+NDI)
            END DO
         END DO
      END IF

!-----  Jacobian matrix: plastic part (slip)
      DO J=1,NDI
         DO I=1,NDI
            DO K=1,NSLPTL
               DDSDDE(I,J)=DDSDDE(I,J)-DDEMSD(I,K)*DDGDDE(K,J)
            END DO
         END DO
      END DO

      IF (NSHR.GT.0) THEN
         DO J=1,NSHR

            DO I=1,NSHR
               DO K=1,NSLPTL
                  DDSDDE(I+NDI,J+NDI)=DDSDDE(I+NDI,J+NDI)-
     2                                DDEMSD(I+3,K)*DDGDDE(K,J+NDI)
               END DO
            END DO

            DO I=1,NDI
               DO K=1,NSLPTL
                  DDSDDE(I,J+NDI)=DDSDDE(I,J+NDI)-
     2                            DDEMSD(I,K)*DDGDDE(K,J+NDI)
                  DDSDDE(J+NDI,I)=DDSDDE(J+NDI,I)-
     2                            DDEMSD(J+3,K)*DDGDDE(K,I)
               END DO
            END DO

         END DO
      END IF

      IF (ITRATN.NE.0) THEN
         DO J=1,NTENS
            DO I=1,NTENS
               DDSDDE(I,J)=DDSDDE(I,J)/(1.+DEV)
            END DO
         END DO
      END IF

!-----  Iteration ?
      IF (ITRATN.NE.0) THEN

!-----  Save solutions (without iteration):
!            Shear strain-rate in a slip system FSLIP1
!            Current strength in a slip system GSLP1
!            Shear strain in a slip system GAMMA1
!            Resolved shear stress in a slip system TAUSP1
!            Normal to a slip plane SPNOR1
!            Slip direction SPDIR1
!            Stress STRES1
!            Jacobian matrix DDSDE1
!
         IF (NITRTN.EQ.0) THEN

            IDNOR=3*NSLPTL
            IDDIR=6*NSLPTL
            DO J=1,NSLPTL
               FSLIP1(J)=FSLIP(J)
               GSLP1(J)=STATEV(J)
               GAMMA1(J)=STATEV(NSLPTL+J)
               TAUSP1(J)=STATEV(2*NSLPTL+J)
               DO I=1,3
                  IDNOR=IDNOR+1
                  SPNOR1(I,J)=STATEV(IDNOR)

                  IDDIR=IDDIR+1
                  SPDIR1(I,J)=STATEV(IDDIR)
               END DO
            END DO

            DO J=1,NTENS
               STRES1(J)=STRESS(J)
               DO I=1,NTENS
                  DDSDE1(I,J)=DDSDDE(I,J)
               END DO
            END DO

         END IF

!-----  Increments of stress DSOLD, and solution dependent state 
!     variables DGAMOD, DTAUOD, DGSPOD, DSPNRO, DSPDRO (for the next 
C     iteration)
!
         DO I=1,NTENS
            DSOLD(I)=DSTRES(I)
         END DO

         DO J=1,NSLPTL
            DGAMOD(J)=DGAMMA(J)
            DTAUOD(J)=DTAUSP(J)
            DGSPOD(J)=DGSLIP(J)
            DO I=1,3
               DSPNRO(I,J)=DSPNOR(I,J)
               DSPDRO(I,J)=DSPDIR(I,J)
            END DO
         END DO

!-----  Check if the iteration solution converges
         IDBACK=0
         ID=0
         DO I=1,NSET
            DO J=1,NSLIP(I)
               ID=ID+1
               X=STATEV(2*NSLPTL+ID)/STATEV(ID)
               RESIDU=THETA*DTIME*F(X,PROPS(65+8*I))+DTIME*(1.0-THETA)*
     2                FSLIP1(ID)-DGAMMA(ID)
               IF (ABS(RESIDU).GT.GAMERR) IDBACK=1
            END DO
         END DO

         IF (IDBACK.NE.0.AND.NITRTN.LT.ITRMAX) THEN
!-----  Iteration: arrays for iteration
!FIXA
            CALL ITERATION (STATEV(NSLPTL+1), STATEV(2*NSLPTL+1), 
     2                      STATEV(1), STATEV(9*NSLPTL+1), 
     3                      STATEV(10*NSLPTL+1), NSLPTL, 
     4                      NSET, NSLIP, ND, PROPS(97), DGAMOD,
     5                      DHDGDG)
!FIXB

            GO TO 1000

         ELSE IF (NITRTN.GE.ITRMAX) THEN
!-----  Solution not converge within maximum number of iteration (the 
!     solution without iteration will be used)
!
            DO J=1,NTENS
               STRESS(J)=STRES1(J)
               DO I=1,NTENS
                  DDSDDE(I,J)=DDSDE1(I,J)
               END DO
            END DO

            IDNOR=3*NSLPTL
            IDDIR=6*NSLPTL
            DO J=1,NSLPTL
               STATEV(J)=GSLP1(J)
               STATEV(NSLPTL+J)=GAMMA1(J)
               STATEV(2*NSLPTL+J)=TAUSP1(J)

               DO I=1,3
                  IDNOR=IDNOR+1
                  STATEV(IDNOR)=SPNOR1(I,J)

                  IDDIR=IDDIR+1
                  STATEV(IDDIR)=SPDIR1(I,J)
               END DO
            END DO

         END IF

      END IF

!-----  Total cumulative shear strains on all slip systems (sum of the 
!       absolute values of shear strains in all slip systems)
!FIX--  Total cumulative shear strains on each slip system (sum of the 
!FIX    absolute values of shear strains in each individual slip system)
!
      DO I=1,NSLPTL
!FIXA
         STATEV(10*NSLPTL+1)=STATEV(10*NSLPTL+1)+ABS(DGAMMA(I))
         STATEV(9*NSLPTL+I)=STATEV(9*NSLPTL+I)+ABS(DGAMMA(I))
!FIXB
      END DO

      RETURN
      END


!----------------------------------------------------------------------


      SUBROUTINE ROTATION (PROP, ROTATE)

!-----  This subroutine calculates the rotation matrix, i.e. the 
!     dire!tion !osines of !ubi! !rystal [100], [010] and [001] 
!     dire!tions in global system

!-----  The rotation matrix is stored in the array ROTATE.

!-----  Use single pre!ision on !ray
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PROP(16), ROTATE(3,3), TERM1(3,3), TERM2(3,3), INDX(3) 

!-----  Subroutines:
!
!       !ROSS  -- !ross produ!t of two ve!tors
!
!       LUD!MP -- LU de!omposition
!
!       LUBKSB -- linear equation solver based on LU de!omposition 
!                 method (must call LUDCMP first)


!-----  PROP -- !onstants !hara!terizing the !rystal orientation 
!               (INPUT)
!
!            PROP(1) - PROP(3) -- dire!tion of the first ve!tor in 
!                                 lo!al !ubi! !rystal system
!            PROP(4) - PROP(6) -- dire!tion of the first ve!tor in 
!                                 global system
!
!            PROP(9) - PROP(11)-- dire!tion of the se!ond ve!tor in 
!                                 lo!al !ubi! !rystal system
!            PROP(12)- PROP(14)-- dire!tion of the se!ond ve!tor in 
!                                 global system
!
!-----  ROTATE -- rotation matrix (OUTPUT):
!
!            ROTATE(i,1) -- dire!tion !osines of dire!tion [1 0 0] in 
!                           lo!al !ubi! !rystal system
!            ROTATE(i,2) -- dire!tion !osines of dire!tion [0 1 0] in 
!                           lo!al !ubi! !rystal system
!            ROTATE(i,3) -- dire!tion !osines of dire!tion [0 0 1] in 
!                           local cubic crystal system

!-----  local matrix: TERM1
      CALL CROSS (PROP(1), PROP(9), TERM1, ANGLE1)

!-----  LU decomposition of TERM1
      CALL LUDCMP (TERM1, 3, 3, INDX, DCMP)

!-----  inverse matrix of TERM1: TERM2
      DO J=1,3
         DO I=1,3
            IF (I.EQ.J) THEN
               TERM2(I,J)=1.
            ELSE
               TERM2(I,J)=0.
            END IF
         END DO
      END DO

      DO J=1,3
         CALL LUBKSB (TERM1, 3, 3, INDX, TERM2(1,J))
      END DO

!-----  global matrix: TERM1
      CALL CROSS (PROP(4), PROP(12), TERM1, ANGLE2)

C-----  !he!k: the angle between first and se!ond ve!tor in lo!al and 
!     global systems must be the same.  The relative differen!e must be
!     less than 0.1%.
!
      IF (ABS(ANGLE1/ANGLE2-1.).GT.0.001) THEN 
         WRITE (6,*) 
     2      '***ERROR - angles between two vectors are not the same'
         STOP
      END IF

!-----  rotation matrix: ROTATE
      DO J=1,3
         DO I=1,3
            ROTATE(I,J)=0.
            DO K=1,3
               ROTATE(I,J)=ROTATE(I,J)+TERM1(I,K)*TERM2(K,J)
            END DO
         END DO
      END DO

      RETURN
      END


!-----------------------------------


           SUBROUTINE CROSS (A, B, C, ANGLE)

!-----  (1) normalize vectors A and B to unit vectors
!       (2) store A, B and A*B (cross product) in C

!-----  Use single precision on cray
!
           IMPLICIT REAL*8 (A-H,O-Z)
           DIMENSION A(3), B(3), C(3,3)

           SUM1=SQRT(A(1)**2+A(2)**2+A(3)**2)
           SUM2=SQRT(B(1)**2+B(2)**2+B(3)**2)

           IF (SUM1.EQ.0.) THEN
              WRITE (6,*) '***ERROR - first vector is zero'
              STOP
           ELSE
              DO I=1,3
                 C(I,1)=A(I)/SUM1
              END DO
           END IF

           IF (SUM2.EQ.0.) THEN
              WRITE (6,*) '***ERROR - second vector is zero'
              STOP
           ELSE
              DO I=1,3
                 C(I,2)=B(I)/SUM2
              END DO
           END IF

           ANGLE=0.
           DO I=1,3
              ANGLE=ANGLE+C(I,1)*C(I,2)
           END DO
           ANGLE=ACOS(ANGLE)

           C(1,3)=C(2,1)*C(3,2)-C(3,1)*C(2,2)
           C(2,3)=C(3,1)*C(1,2)-C(1,1)*C(3,2)
           C(3,3)=C(1,1)*C(2,2)-C(2,1)*C(1,2)
           SUM3=SQRT(C(1,3)**2+C(2,3)**2+C(3,3)**2)
           IF (SUM3.LT.1.E-8) THEN
              WRITE (6,*) 
     2           '***ERROR - first and second vectors are parallel'
               STOP
            END IF

           RETURN
           END


C----------------------------------------------------------------------


      SUBROUTINE SLIPSYS (ISPDIR, ISPNOR, NSLIP, SLPDIR, SLPNOR, 
     2                    ROTATE)

!-----  This subroutine generates all slip systems in the same set for 
!     a CUBIC crystal.  For other crystals (e.g., HCP, Tetragonal, 
!     Orthotropic, ...), it has to be modified to include the effect of
!     crystal aspect ratio.

C-----  Denote s as a slip dire!tion and m as normal to a slip plane.  
!     In a cubic crystal, (s,-m), (-s,m) and (-s,-m) are NOT considered
!     independent of (s,m).

!-----  Subroutines:  LINE1 and LINE

!-----  Variables:
!
!     ISPDIR -- a typical slip direction in this set of slip systems 
C               (integer)  (INPUT)
!     ISPNOR -- a typical normal to slip plane in this set of slip 
!               systems (integer)  (INPUT)
C     NSLIP  -- number of independent slip systems in this set 
!               (OUTPUT)
!     SLPDIR -- unit vectors of all slip directions  (OUTPUT)
!     SLPNOR -- unit normals to all slip planes  (OUTPUT)
!     ROTATE -- rotation matrix (INPUT)
!          ROTATE(i,1) -- dire!tion !osines of [100] in global system
!          ROTATE(i,2) -- dire!tion !osines of [010] in global system
!          ROTATE(i,3) -- dire!tion !osines of [001] in global system
!
!     NSPDIR -- number of all possible slip dire!tions in this set
!     NSPNOR -- number of all possible slip planes in this set
!     IWKDIR -- all possible slip dire!tions (integer)
!     IWKNOR -- all possible slip planes (integer)


!-----  Use single pre!ision on !ray
!
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ISPDIR(3), ISPNOR(3), SLPDIR(3,50), SLPNOR(3,50), 
     *          ROTATE(3,3), IWKDIR(3,24), IWKNOR(3,24), TERM(3)

      NSLIP=0
      NSPDIR=0
      NSPNOR=0

!-----  Generating all possible slip dire!tions in this set
!
!       Denote the slip dire!tion by [lmn].  I1 is the minimum of the 
!     absolute value of l, m and n, I3 is the maximum and I2 is the 
!     mode, e.g. (1 -3 2), I1=1, I2=2 and I3=3.  I1<=I2<=I3.

      I1=MIN(IABS(ISPDIR(1)),IABS(ISPDIR(2)),IABS(ISPDIR(3)))
      I3=MAX(IABS(ISPDIR(1)),IABS(ISPDIR(2)),IABS(ISPDIR(3)))
      I2=IABS(ISPDIR(1))+IABS(ISPDIR(2))+IABS(ISPDIR(3))-I1-I3

      RMODIR=SQRT(FLOAT(I1*I1+I2*I2+I3*I3))

!     I1=I2=I3=0
      IF (I3.EQ.0) THEN 
         WRITE (6,*) '***ERROR - slip direction is [000]'
         STOP

!     I1=I2=0, I3>0   ---   [001] type
      ELSE IF (I2.EQ.0) THEN
         NSPDIR=3
         DO J=1,3
            DO I=1,3
               IWKDIR(I,J)=0
               IF (I.EQ.J) IWKDIR(I,J)=I3
            END DO
         END DO

!     I1=0, I3>=I2>0
      ELSE IF (I1.EQ.0) THEN

!        I1=0, I3=I2>0   ---   [011] type
         IF (I2.EQ.I3) THEN
            NSPDIR=6
            DO J=1,6
               DO I=1,3
                  IWKDIR(I,J)=I2
                  IF (I.EQ.J.OR.J-I.EQ.3) IWKDIR(I,J)=0
                  IWKDIR(1,6)=-I2
                  IWKDIR(2,4)=-I2
                  IWKDIR(3,5)=-I2
               END DO
            END DO

!        I1=0, I3>I2>0   ---   [012] type
         ELSE
            NSPDIR=12
            CALL LINE1 (I2, I3, IWKDIR(1,1), 1)
            CALL LINE1 (I3, I2, IWKDIR(1,3), 1)
            CALL LINE1 (I2, I3, IWKDIR(1,5), 2)
            CALL LINE1 (I3, I2, IWKDIR(1,7), 2)
            CALL LINE1 (I2, I3, IWKDIR(1,9), 3)
            CALL LINE1 (I3, I2, IWKDIR(1,11), 3)

         END IF

!     I1=I2=I3>0   ---   [111] type
      ELSE IF (I1.EQ.I3) THEN
         NSPDIR=4
         CALL LINE (I1, I1, I1, IWKDIR)

!     I3>I2=I1>0   ---   [112] type
      ELSE IF (I1.EQ.I2) THEN
         NSPDIR=12
         CALL LINE (I1, I1, I3, IWKDIR(1,1))
         CALL LINE (I1, I3, I1, IWKDIR(1,5))
         CALL LINE (I3, I1, I1, IWKDIR(1,9))

!     I3=I2>I1>0   ---   [122] type
      ELSE IF (I2.EQ.I3) THEN
         NSPDIR=12
         CALL LINE (I1, I2, I2, IWKDIR(1,1))
         CALL LINE (I2, I1, I2, IWKDIR(1,5))
         CALL LINE (I2, I2, I1, IWKDIR(1,9))

!     I3>I2>I1>0   ---   [123] type
      ELSE
         NSPDIR=24
         CALL LINE (I1, I2, I3, IWKDIR(1,1))
         CALL LINE (I3, I1, I2, IWKDIR(1,5))
         CALL LINE (I2, I3, I1, IWKDIR(1,9))
         CALL LINE (I1, I3, I2, IWKDIR(1,13))
         CALL LINE (I2, I1, I3, IWKDIR(1,17))
         CALL LINE (I3, I2, I1, IWKDIR(1,21))

      END IF

C-----  Generating all possible slip planes in this set
!
!       Denote the normal to slip plane by (pqr).  J1 is the minimum of
!     the absolute value of p, q and r, J3 is the maximum and J2 is the
!     mode, e.g. (1 -2 1), J1=1, J2=1 and J3=2.  J1<=J2<=J3.

      J1=MIN(IABS(ISPNOR(1)),IABS(ISPNOR(2)),IABS(ISPNOR(3)))
      J3=MAX(IABS(ISPNOR(1)),IABS(ISPNOR(2)),IABS(ISPNOR(3)))
      J2=IABS(ISPNOR(1))+IABS(ISPNOR(2))+IABS(ISPNOR(3))-J1-J3

      RMONOR=SQRT(FLOAT(J1*J1+J2*J2+J3*J3))

      IF (J3.EQ.0) THEN 
         WRITE (6,*) '***ERROR - slip plane is [000]'
         STOP

!     (001) type
      ELSE IF (J2.EQ.0) THEN
         NSPNOR=3
         DO J=1,3
            DO I=1,3
               IWKNOR(I,J)=0
               IF (I.EQ.J) IWKNOR(I,J)=J3
            END DO
         END DO

      ELSE IF (J1.EQ.0) THEN

!     (011) type
         IF (J2.EQ.J3) THEN
            NSPNOR=6
            DO J=1,6
               DO I=1,3
                  IWKNOR(I,J)=J2
                  IF (I.EQ.J.OR.J-I.EQ.3) IWKNOR(I,J)=0
                  IWKNOR(1,6)=-J2
                  IWKNOR(2,4)=-J2
                  IWKNOR(3,5)=-J2
               END DO
            END DO

!     (012) type
         ELSE
            NSPNOR=12
            CALL LINE1 (J2, J3, IWKNOR(1,1), 1)
            CALL LINE1 (J3, J2, IWKNOR(1,3), 1)
            CALL LINE1 (J2, J3, IWKNOR(1,5), 2)
            CALL LINE1 (J3, J2, IWKNOR(1,7), 2)
            CALL LINE1 (J2, J3, IWKNOR(1,9), 3)
            CALL LINE1 (J3, J2, IWKNOR(1,11), 3)

         END IF

!     (111) type
      ELSE IF (J1.EQ.J3) THEN
         NSPNOR=4
         CALL LINE (J1, J1, J1, IWKNOR)

!     (112) type
      ELSE IF (J1.EQ.J2) THEN
         NSPNOR=12
         CALL LINE (J1, J1, J3, IWKNOR(1,1))
         CALL LINE (J1, J3, J1, IWKNOR(1,5))
         CALL LINE (J3, J1, J1, IWKNOR(1,9))

!     (122) type
      ELSE IF (J2.EQ.J3) THEN
         NSPNOR=12
         CALL LINE (J1, J2, J2, IWKNOR(1,1))
         CALL LINE (J2, J1, J2, IWKNOR(1,5))
         CALL LINE (J2, J2, J1, IWKNOR(1,9))

!     (123) type
      ELSE
         NSPNOR=24
         CALL LINE (J1, J2, J3, IWKNOR(1,1))
         CALL LINE (J3, J1, J2, IWKNOR(1,5))
         CALL LINE (J2, J3, J1, IWKNOR(1,9))
         CALL LINE (J1, J3, J2, IWKNOR(1,13))
         CALL LINE (J2, J1, J3, IWKNOR(1,17))
         CALL LINE (J3, J2, J1, IWKNOR(1,21))

      END IF

!-----  Generating all slip systems in this set
!
!-----  Unit ve!tors in slip dire!tions: SLPDIR, and unit normals to 
!     slip planes: SLPNOR in lo!al !ubi! !rystal system
!
      WRITE (6,*) '          '
      WRITE (6,*) ' #          Slip plane          Slip direction'

!      DO J=1,NSPNOR
         NSLIP=12
         DO I=1,NSLIP
               
               DO K=1,3
                  SLPDIR(K,I)=1/SQRT(2.)
                  SLPNOR(K,I)=1/SQRT(3.)
                  IF (I.EQ.K.OR.I-K.EQ.3) SLPDIR(K,I)=0
                  IF (I-K.EQ.6) SLPDIR(K,I)=0
                  IF (I-K.EQ.9) SLPDIR(K,I)=0

                  SLPDIR(1,2)=-1/SQRT(2.)
                  SLPDIR(1,5)=-1/SQRT(2.)
                  SLPDIR(1,9)=-1/SQRT(2.)
                  SLPDIR(1,12)=-1/SQRT(2.)
                  SLPDIR(2,3)=-1/SQRT(2.)
                  SLPDIR(2,4)=-1/SQRT(2.)
                  SLPDIR(2,7)=-1/SQRT(2.)
                  SLPDIR(2,12)=-1/SQRT(2.)
                  SLPDIR(3,1)=-1/SQRT(2.)
                  SLPDIR(3,5)=-1/SQRT(2.)
                  SLPDIR(3,7)=-1/SQRT(2.)
                  SLPDIR(3,11)=-1/SQRT(2.)
                  
                  SLPNOR(1,7)=-1/SQRT(3.)
                  SLPNOR(1,8)=-1/SQRT(3.)
                  SLPNOR(1,9)=-1/SQRT(3.)
                  SLPNOR(1,10)=-1/SQRT(3.)
                  SLPNOR(1,11)=-1/SQRT(3.)
                  SLPNOR(1,12)=-1/SQRT(3.)
                  SLPNOR(2,4)=-1/SQRT(3.)
                  SLPNOR(2,5)=-1/SQRT(3.)
                  SLPNOR(2,6)=-1/SQRT(3.)
                  SLPNOR(2,7)=-1/SQRT(3.)
                  SLPNOR(2,8)=-1/SQRT(3.)
                  SLPNOR(2,9)=-1/SQRT(3.)
                  SLPNOR(3,4)=-1/SQRT(3.)
                  SLPNOR(3,5)=-1/SQRT(3.)
                  SLPNOR(3,6)=-1/SQRT(3.)
                  SLPNOR(3,10)=-1/SQRT(3.)
                  SLPNOR(3,11)=-1/SQRT(3.)
                  SLPNOR(3,12)=-1/SQRT(3.)
               END DO 
         END DO
         
      
         
         
10    FORMAT(1X,I2,9X,'(',3(1X,I2),1X,')',10X,'[',3(1X,I2),1X,']')

      WRITE (6,*) 'Number of slip systems in this set = ',NSLIP
      WRITE (6,*) '          '

      IF (NSLIP.EQ.0) THEN
         WRITE (6,*) 
     *      'There is no slip direction normal to the slip planes!'
         STOP

      ELSE
          
C-----  Unit vectors in slip directions: SLPDIR, and unit normals to 
!     slip planes: SLPNOR in global system
!
         DO J=1,NSLIP
            DO I=1,3
               TERM(I)=0.
               DO K=1,3
                  TERM(I)=TERM(I)+ROTATE(I,K)*SLPDIR(K,J)
               END DO
            END DO
            DO I=1,3
               SLPDIR(I,J)=TERM(I)
            END DO

            DO I=1,3
               TERM(I)=0.
               DO K=1,3
                  TERM(I)=TERM(I)+ROTATE(I,K)*SLPNOR(K,J)
               END DO
            END DO
            DO I=1,3
               SLPNOR(I,J)=TERM(I)
            END DO
         END DO

      END IF

      RETURN
      END


!----------------------------------


           SUBROUTINE LINE (I1, I2, I3, IARRAY)

!-----  Generating all possible slip directions <lmn> (or slip planes 
!     {lmn}) for a cubic crystal, where l,m,n are not zeros.

C-----  Use single precision on cray
!
           IMPLICIT REAL*8 (A-H,O-Z)
           DIMENSION IARRAY(3,4)

           DO J=1,4
              IARRAY(1,J)=I1
              IARRAY(2,J)=I2
              IARRAY(3,J)=I3
           END DO

           DO I=1,3
              DO J=1,4
                 IF (J.EQ.I+1) IARRAY(I,J)=-IARRAY(I,J)
              END DO
           END DO

           RETURN
           END


C-----------------------------------


           SUBROUTINE LINE1 (J1, J2, IARRAY, ID)

!-----  Generating all possible slip directions <0mn> (or slip planes 
!     {0mn}) for a cubic crystal, where m,n are not zeros and m does 
!     not equal n.

!-----  Use single precision on cray
!
           IMPLICIT REAL*8 (A-H,O-Z)
           DIMENSION IARRAY(3,2)

           IARRAY(ID,1)=0
           IARRAY(ID,2)=0

           ID1=ID+1
           IF (ID1.GT.3) ID1=ID1-3
           IARRAY(ID1,1)=J1
           IARRAY(ID1,2)=J1

           ID2=ID+2
           IF (ID2.GT.3) ID2=ID2-3
           IARRAY(ID2,1)=J2
           IARRAY(ID2,2)=-J2
  
           RETURN
           END


!----------------------------------------------------------------------


      SUBROUTINE GSLPINIT (GSLIP0, NSLIP, NSLPTL, NSET, PROP)

!-----  This subroutine !al!ulates the initial value of !urrent 
!     strength for ea!h slip system in a rate-dependent single !rystal.
!     Two sets of initial values, proposed by Asaro, Pier!e et al, and 
!     by Bassani, respe!tively, are used here.  Both sets assume that 
!     the initial values for all slip systems are the same (initially 
!     isotropi!).

!-----  These initial values are assumed the same for all slip systems 
!     in ea!h set, though they !ould be different from set to set, e.g.
!     <110>{111} and <110>{100}.

!-----  Users who want to use their own initial values may !hange the 
!     fun!tion subprogram GSLP0.  The parameters !hara!terizing these 
!     initial values are passed into GSLP0 through array PROP.

!-----  Use single pre!ision on !ray
!
      IMPLICIT REAL*8 (A-H,O-Z)
      EXTERNAL GSLP0
      DIMENSION GSLIP0(NSLPTL), NSLIP(NSET), PROP(16,NSET)

!-----  Fun!tion subprograms:
!
!       GSLP0 -- User-supplied fun!tion subprogram given the initial 
!                value of !urrent strength at initial state

!-----  Variables:
!
!     GSLIP0 -- initial value of !urrent strength (OUTPUT)
!
!     NSLIP  -- number of slip systems in ea!h set (INPUT)
!     NSLPTL -- total number of slip systems in all the sets (INPUT)
!     NSET   -- number of sets of slip systems (INPUT)
!
!     PROP   -- material !onstants !hara!terizing the initial value of 
!               !urrent strength (INPUT)
!
!               For Asaro, Pier!e et al's law 
!               PROP(1,i) -- initial hardening modulus H0 in the ith 
!                            set of slip systems
!               PROP(2,i) -- saturation stress TAUs in the ith set of  
!                            slip systems
!               PROP(3,i) -- initial !riti!al resolved shear stress 
!                            TAU0 in the ith set of slip systems
!
!               For Bassani's law 
!               PROP(1,i) -- initial hardening modulus H0 in the ith 
!                            set of slip systems
!               PROP(2,i) -- stage I stress TAUI in the ith set of  
!                            slip systems (or the breakthrough stress 
!                            where large plasti! flow initiates)
!               PROP(3,i) -- initial !riti!al resolved shear stress 
!                            TAU0 in the ith set of slip systems
!

      ID=0
      DO I=1,NSET
         ISET=I
         DO J=1,NSLIP(I)
            ID=ID+1
            GSLIP0(ID)=GSLP0(NSLPTL,NSET,NSLIP,PROP(1,I),ID,ISET)
         END DO
      END DO

      RETURN
      END


!----------------------------------


!-----  Use single pre!ision on !ray
!
           REAL*8 FUNCTION GSLP0(NSLPTL,NSET,NSLIP,PROP,ISLIP,ISET)

!-----     User-supplied fun!tion subprogram given the initial value of
!        !urrent strength at initial state

!-----  Use single pre!ision on !ray
!
           IMPLICIT REAL*8 (A-H,O-Z)
           DIMENSION NSLIP(NSET), PROP(16)

           GSLP0=PROP(3)

           RETURN
           END


!----------------------------------------------------------------------


      SUBROUTINE STRAINRATE (GAMMA, TAUSLP, GSLIP, NSLIP, FSLIP, 
     2                       DFDXSP, PROP)
      
!-----  This subroutine !al!ulates the shear strain-rate in ea!h slip 
!     system for a rate-dependent single !rystal.  The POWER LAW 
!     relation between shear strain-rate and resolved shear stress 
!     proposed by Hut!hinson, Pan and Ri!e, is used here.

!-----  The power law exponents are assumed the same for all slip 
!     systems in ea!h set, though they !ould be different from set to 
!     set, e.g. <110>{111} and <110>{100}.  The strain-rate !oeffi!ient
!     in front of the power law form are also assumed the same for all 
!     slip systems in ea!h set. 

!-----  Users who want to use their own !onstitutive relation may 
!     !hange the fun!tion subprograms F and its derivative DFDX, 
!     where F is the strain hardening law, dGAMMA/dt = F(X), 
!     X=TAUSLP/GSLIP.  The parameters !hara!terizing F are passed into 
!     F and DFDX through array PROP.

!-----  Fun!tion subprograms:
!
!       F    -- User-supplied fun!tion subprogram whi!h gives shear 
!               strain-rate for ea!h slip system based on !urrent 
!               values of resolved shear stress and !urrent strength
!
!       DFDX -- User-supplied fun!tion subprogram dF/dX, where x is the
!               ratio of resolved shear stress over !urrent strength

!-----  Variables:
!
!     GAMMA  -- shear strain in ea!h slip system at the start of time 
!               step  (INPUT)
!     TAUSLP -- resolved shear stress in ea!h slip system (INPUT)
!     GSLIP  -- !urrent strength (INPUT)
!     NSLIP  -- number of slip systems in this set (INPUT)
!
!     FSLIP  -- !urrent value of F for ea!h slip system (OUTPUT)
!     DFDXSP -- !urrent value of DFDX for ea!h slip system (OUTPUT)
!
!     PROP   -- material !onstants !hara!terizing the strain hardening 
!               law (INPUT)
!
!               For the !urrent power law strain hardening law 
!               PROP(1) -- power law hardening exponent
!               PROP(1) = infinity !orresponds to a rate-independent 
!               material
!               PROP(2) -- !oeffi!ient in front of power law hardening


!-----  Use single pre!ision on !ray
!
      IMPLICIT REAL*8 (A-H,O-Z)
      EXTERNAL F, DFDX
      DIMENSION GAMMA(NSLIP), TAUSLP(NSLIP), GSLIP(NSLIP), 
     2          FSLIP(NSLIP), DFDXSP(NSLIP), PROP(8)
      do i=1,12
      
      end do
      


      DO I=1,NSLIP
         X=TAUSLP(I)/GSLIP(I)
         FSLIP(I)=F(X,PROP)
         DFDXSP(I)=DFDX(X,PROP)
      END DO

      RETURN
      END


!-----------------------------------


!-----  Use single precision on cray
!
           REAL*8 FUNCTION F(X,PROP)

!-----     User-supplied function subprogram which gives shear 
!        strain-rate for each slip system based on current values of 
!        resolved shear stress and current strength
C
!-----  Use single precision on cray
!
           IMPLICIT REAL*8 (A-H,O-Z)
           DIMENSION PROP(8)

           F=PROP(2)*(ABS(X))**PROP(1)*DSIGN(1.D0,X)
                 
           RETURN
           END


!-----------------------------------


!-----  Use single precision on cray
!
           REAL*8 FUNCTION DFDX(X,PROP)

!-----     User-supplied function subprogram dF/dX, where x is the 
!        ratio of resolved shear stress over current strength

!-----  Use single precision on cray
!
           IMPLICIT REAL*8 (A-H,O-Z)
           DIMENSION PROP(8)

           DFDX=PROP(1)*PROP(2)*(ABS(X))**(PROP(1)-1.)

           RETURN
           END


!----------------------------------------------------------------------

!FIXA
      SUBROUTINE LATENTHARDEN (GAMMA, TAUSLP, GSLIP, GMSLTL, GAMTOL, 
     2                         NSLIP, NSLPTL, NSET, H, PROP, ND)
!FIXB

!-----  This subroutine !al!ulates the !urrent self- and latent-
!     hardening moduli for all slip systems in a rate-dependent single 
!     !rystal.  Two kinds of hardening law are used here.  The first 
!     law, proposed by Asaro, and Pier!e et al, assumes a HYPER SE!ANT 
!     relation between self- and latent-hardening moduli and overall 
!     shear strain.  The Baus!hinger effe!t has been negle!ted.  The 
!     se!ond is Bassani's hardening law, whi!h gives an expli!it 
!     expression of slip intera!tions between slip systems.  The 
!     !lassi!al three stage hardening for F!! single !rystal !ould be 
!     simulated.

!-----  The hardening !oeffi!ients are assumed the same for all slip 
!     systems in ea!h set, though they !ould be different from set to 
!     set, e.g. <110>{111} and <110>{100}.

!-----  Users who want to use their own self- and latent-hardening law 
!     may !hange the fun!tion subprograms HSELF (self hardening) and 
!     HLATNT (latent hardening).  The parameters !hara!terizing these 
!     hardening laws are passed into HSELF and HLATNT through array 
!     PROP.


!-----  Fun!tion subprograms:
!
!       HSELF  -- User-supplied self-hardening fun!tion in a slip 
!                 system
!
!       HLATNT -- User-supplied latent-hardening fun!tion

!-----  Variables:
!
!     GAMMA  -- shear strain in all slip systems at the start of time 
!               step  (INPUT)
!     TAUSLP -- resolved shear stress in all slip systems (INPUT)
!     GSLIP  -- !urrent strength (INPUT)
!FIX  GMSLTL -- total !umulative shear strains on ea!h individual slip system 
!FIX            (INPUT)
!     GAMTOL -- total !umulative shear strains over all slip systems 
!               (INPUT)
!     NSLIP  -- number of slip systems in ea!h set (INPUT)
!     NSLPTL -- total number of slip systems in all the sets (INPUT)
!     NSET   -- number of sets of slip systems (INPUT)
!
!     H      -- !urrent value of self- and latent-hardening moduli 
!               (OUTPUT)
!               H(i,i) -- self-hardening modulus of the ith slip system
!                         (no sum over i)
!               H(i,j) -- latent-hardening molulus of the ith slip 
!                         system due to a slip in the jth slip system 
!                         (i not equal j)
!
!     PROP   -- material !onstants !hara!terizing the self- and latent-
!               hardening law (INPUT)
!
!               For the HYPER SE!ANT hardening law 
!               PROP(1,i) -- initial hardening modulus H0 in the ith 
!                            set of slip systems
!               PROP(2,i) -- saturation stress TAUs in the ith set of  
!                            slip systems
!               PROP(3,i) -- initial !riti!al resolved shear stress 
!                            TAU0 in the ith set of slip systems
!               PROP(9,i) -- ratio of latent to self-hardening Q in the
!                            ith set of slip systems
!               PROP(10,i)-- ratio of latent-hardening from other sets 
!                            of slip systems to self-hardening in the 
!                            ith set of slip systems Q1
!
!               For Bassani's hardening law 
!               PROP(1,i) -- initial hardening modulus H0 in the ith 
!                            set of slip systems
!               PROP(2,i) -- stage I stress TAUI in the ith set of  
!                            slip systems (or the breakthrough stress 
!                            where large plasti! flow initiates)
!               PROP(3,i) -- initial !riti!al resolved shear stress 
!                            TAU0 in the ith set of slip systems
!               PROP(4,i) -- hardening modulus during easy glide Hs in 
!                            the ith set of slip systems
!               PROP(5,i) -- amount of slip Gamma0 after whi!h a given 
C                            interaction between slip systems in the 
C                            ith set reaches peak strength
C               PROP(6,i) -- amount of slip Gamma0 after whi!h a given 
!                            intera!tion between slip systems in the 
!                            ith set and jth set (i not equal j) 
!                            rea!hes peak strength
!               PROP(7,i) -- representing the magnitude of the strength
!                            of intera!tion in the ith set of slip 
!                            system
!               PROP(8,i) -- representing the magnitude of the strength
!                            of intera!tion between the ith set and jth
!                            set of system
!               PROP(9,i) -- ratio of latent to self-hardening Q in the
!                            ith set of slip systems
!               PROP(10,i)-- ratio of latent-hardening from other sets 
!                            of slip systems to self-hardening in the 
!                            ith set of slip systems Q1
!
!     ND     -- leading dimension of arrays defined in subroutine UMAT 
!               (INPUT) 


C-----  Use single precision on cray
C
      IMPLICIT REAL*8 (A-H,O-Z)
      EXTERNAL HSELF, HLATNT
CFIXA
      DIMENSION GAMMA(NSLPTL), TAUSLP(NSLPTL), GMSLTL(NSLPTL),
     2          GSLIP(NSLPTL), NSLIP(NSET), PROP(16,NSET), 
     3          H(ND,NSLPTL)
CFIXB

      CHECK=0.
      DO I=1,NSET
         DO J=4,8
            CHECK=CHECK+ABS(PROP(J,I))
         END DO
      END DO

C-----  CHECK=0   --  HYPER SECANT hardening law
C       otherwise --  Bassani's hardening law

      ISELF=0
      DO I=1,NSET
         ISET=I
         DO J=1,NSLIP(I)
            ISELF=ISELF+1

            DO LATENT=1,NSLPTL
               IF (LATENT.EQ.ISELF) THEN
CFIXA
                  H(LATENT,ISELF)=HSELF(GAMMA,GMSLTL,GAMTOL,NSLPTL,
     2                                  NSET,NSLIP,PROP(1,I),CHECK,
     3                                  ISELF,ISET)
CFIXB
               ELSE
CFIXA
                  H(LATENT,ISELF)=HLATNT(GAMMA,GMSLTL,GAMTOL,NSLPTL,
     2                                   NSET,NSLIP,PROP(1,I),CHECK,
     3                                   ISELF,ISET,LATENT)
CFIXB

               END IF
            END DO

         END DO
      END DO

      RETURN
      END


C-----------------------------------


C-----  Use single precision on cray
CFIXA
           REAL*8 FUNCTION HSELF(GAMMA,GMSLTL,GAMTOL,NSLPTL,NSET,
     2                           NSLIP,PROP,CHECK,ISELF,ISET)
CFIXB

C-----     User-supplied self-hardening function in a slip system

C-----  Use single precision on cray
C
           IMPLICIT REAL*8 (A-H,O-Z)
CFIXA
           DIMENSION GAMMA(NSLPTL), NSLIP(NSET), PROP(16),
     2               GMSLTL(NSLPTL)
CFIXB

           IF (CHECK.EQ.0.) THEN

C-----  HYPER SECANT hardening law by Asaro, Pierce et al
              TERM1=PROP(1)*GAMTOL/(PROP(2)-PROP(3))
              TERM2=2.*EXP(-TERM1)/(1.+EXP(-2.*TERM1))
              HSELF=PROP(1)*TERM2**2

           ELSE

C-----  Bassani's hardening law
CFIXA
              TERM1=(PROP(1)-PROP(4))*GMSLTL(ISELF)/(PROP(2)-PROP(3))
CFIXB
              TERM2=2.*EXP(-TERM1)/(1.+EXP(-2.*TERM1))
              F=(PROP(1)-PROP(4))*TERM2**2+PROP(4)

              ID=0
              G=1.
              DO I=1,NSET
                 IF (I.EQ.ISET) THEN
                    GAMMA0=PROP(5)
                    FAB=PROP(7)
                 ELSE
                    GAMMA0=PROP(6)
                    FAB=PROP(8)
                 END IF

                 DO J=1,NSLIP(I)
                    ID=ID+1
                    IF (ID.NE.ISELF) THEN
CFIXA
		       G=G+FAB*TANH(GMSLTL(ID)/GAMMA0)
CFIXB
		    END IF

                 END DO
              END DO

              HSELF=F*G

           END IF

           RETURN
           END


C-----------------------------------


C-----  Use single precision on cray
CFIXA
           REAL*8 FUNCTION HLATNT(GAMMA,GMSLTL,GAMTOL,NSLPTL,NSET,
     2                            NSLIP,PROP,CHECK,ISELF,ISET,LATENT)
CFIXB

C-----     User-supplied latent-hardening function

C-----  Use single precision on cray
C
           IMPLICIT REAL*8 (A-H,O-Z)
CFIXA
           DIMENSION GAMMA(NSLPTL), NSLIP(NSET), PROP(16),
     2               GMSLTL(NSLPTL)
CFIXB

           ILOWER=0
           IUPPER=NSLIP(1)
           IF (ISET.GT.1) THEN
              DO K=2,ISET
                 ILOWER=ILOWER+NSLIP(K-1)
                 IUPPER=IUPPER+NSLIP(K)
              END DO
           END IF

           IF (LATENT.GT.ILOWER.AND.LATENT.LE.IUPPER) THEN
              Q=PROP(9)
           ELSE
              Q=PROP(10)
           END IF

           IF (CHECK.EQ.0.) THEN

C-----  HYPER SECANT hardening law by Asaro, Pierce et al
              TERM1=PROP(1)*GAMTOL/(PROP(2)-PROP(3))
              TERM2=2.*EXP(-TERM1)/(1.+EXP(-2.*TERM1))
              HLATNT=PROP(1)*TERM2**2*Q

           ELSE

C-----  Bassani's hardening law
CFIXA
              TERM1=(PROP(1)-PROP(4))*GMSLTL(ISELF)/(PROP(2)-PROP(3))
CFIXB
              TERM2=2.*EXP(-TERM1)/(1.+EXP(-2.*TERM1))
              F=(PROP(1)-PROP(4))*TERM2**2+PROP(4)

              ID=0
              G=1.
              DO I=1,NSET
                 IF (I.EQ.ISET) THEN
                    GAMMA0=PROP(5)
                    FAB=PROP(7)
                 ELSE
                    GAMMA0=PROP(6)
                    FAB=PROP(8)
                 END IF

                 DO J=1,NSLIP(I)
                    ID=ID+1
                    IF (ID.NE.ISELF) THEN
CFIXA
		       G=G+FAB*TANH(GMSLTL(ID)/GAMMA0)
CFIXB
		    END IF

                 END DO
              END DO

              HLATNT=F*G*Q

           END IF

           RETURN
           END


C----------------------------------------------------------------------

CFIXA
      SUBROUTINE ITERATION (GAMMA, TAUSLP, GSLIP, GMSLTL, GAMTOL, 
     2                      NSLPTL, NSET, NSLIP, ND, PROP, DGAMOD, 
     3                      DHDGDG)
CFIXB

C-----  This subroutine generates arrays for the Newton-Rhapson 
C     iteration method.

C-----  Users who want to use their own self- and latent-hardening law 
C     may change the function subprograms DHSELF (self hardening) and 
C     DHLATN (latent hardening).  The parameters characterizing these 
C     hardening laws are passed into DHSELF and DHLATN through array 
C     PROP.


C-----  Function subprograms:
C
C       DHSELF -- User-supplied function of the derivative of self-
C                 hardening moduli
C
C       DHLATN -- User-supplied function of the derivative of latent-
C                 hardening moduli

C-----  Variables:
C
C     GAMMA  -- shear strain in all slip systems at the start of time 
C               step  (INPUT)
C     TAUSLP -- resolved shear stress in all slip systems (INPUT)
C     GSLIP  -- current strength (INPUT)
CFIX  GMSLTL -- total cumulative shear strains on each individual slip system 
CFIX            (INPUT)
C     GAMTOL -- total cumulative shear strains over all slip systems 
C               (INPUT)
C     NSLPTL -- total number of slip systems in all the sets (INPUT)
C     NSET   -- number of sets of slip systems (INPUT)
C     NSLIP  -- number of slip systems in each set (INPUT)
C     ND     -- leading dimension of arrays defined in subroutine UMAT 
C               (INPUT) 
C
C     PROP   -- material constants characterizing the self- and latent-
C               hardening law (INPUT)
C
C               For the HYPER SECANT hardening law 
C               PROP(1,i) -- initial hardening modulus H0 in the ith 
C                            set of slip systems
C               PROP(2,i) -- saturation stress TAUs in the ith set of  
C                            slip systems
C               PROP(3,i) -- initial critical resolved shear stress 
C                            TAU0 in the ith set of slip systems
C               PROP(9,i) -- ratio of latent to self-hardening Q in the
C                            ith set of slip systems
C               PROP(10,i)-- ratio of latent-hardening from other sets 
C                            of slip systems to self-hardening in the 
C                            ith set of slip systems Q1
C
C               For Bassani's hardening law 
C               PROP(1,i) -- initial hardening modulus H0 in the ith 
C                            set of slip systems
C               PROP(2,i) -- stage I stress TAUI in the ith set of  
C                            slip systems (or the breakthrough stress 
C                            where large plastic flow initiates)
C               PROP(3,i) -- initial critical resolved shear stress 
C                            TAU0 in the ith set of slip systems
C               PROP(4,i) -- hardening modulus during easy glide Hs in 
C                            the ith set of slip systems
C               PROP(5,i) -- amount of slip Gamma0 after which a given 
C                            interaction between slip systems in the 
C                            ith set reaches peak strength
C               PROP(6,i) -- amount of slip Gamma0 after which a given 
C                            interaction between slip systems in the 
C                            ith set and jth set (i not equal j) 
C                            reaches peak strength
C               PROP(7,i) -- representing the magnitude of the strength
C                            of interaction in the ith set of slip 
C                            system
C               PROP(8,i) -- representing the magnitude of the strength
C                            of interaction between the ith set and jth
C                            set of system
C               PROP(9,i) -- ratio of latent to self-hardening Q in the
C                            ith set of slip systems
C               PROP(10,i)-- ratio of latent-hardening from other sets 
C                            of slip systems to self-hardening in the 
C                            ith set of slip systems Q1
C
C-----  Arrays for iteration:
C
C       DGAMOD (INPUT)
C
C       DHDGDG (OUTPUT)
C

C-----  Use single precision on cray
C
      IMPLICIT REAL*8 (A-H,O-Z)
      EXTERNAL DHSELF, DHLATN
CFIXA
      DIMENSION GAMMA(NSLPTL), TAUSLP(NSLPTL), GMSLTL(NSLPTL),
     2          GSLIP(NSLPTL), NSLIP(NSET), PROP(16,NSET), 
     3          DGAMOD(NSLPTL), DHDGDG(ND,NSLPTL)
CFIXB

      CHECK=0.
      DO I=1,NSET
         DO J=4,8
            CHECK=CHECK+ABS(PROP(J,I))
         END DO
      END DO

C-----  CHECK=0   --  HYPER SECANT hardening law
C       otherwise --  Bassani's hardening law

      ISELF=0
      DO I=1,NSET
         ISET=I
         DO J=1,NSLIP(I)
            ISELF=ISELF+1

            DO KDERIV=1,NSLPTL
               DHDGDG(ISELF,KDERIV)=0.

               DO LATENT=1,NSLPTL
                  IF (LATENT.EQ.ISELF) THEN
CFIXA
                     DHDG=DHSELF(GAMMA,GMSLTL,GAMTOL,NSLPTL,NSET,
     2                           NSLIP,PROP(1,I),CHECK,ISELF,ISET,
     3                           KDERIV)
CFIXB
                  ELSE
CFIXA
                     DHDG=DHLATN(GAMMA,GMSLTL,GAMTOL,NSLPTL,NSET,
     2                           NSLIP,PROP(1,I),CHECK,ISELF,ISET,
     3                           LATENT,KDERIV)
CFIXB
                  END IF

                  DHDGDG(ISELF,KDERIV)=DHDGDG(ISELF,KDERIV)+
     2                                 DHDG*ABS(DGAMOD(LATENT))
               END DO

            END DO
         END DO
      END DO

      RETURN
      END


C-----------------------------------


C-----  Use single precision on cray
CFIXA
           REAL*8 FUNCTION DHSELF(GAMMA,GMSLTL,GAMTOL,NSLPTL,NSET,
     2                            NSLIP,PROP,CHECK,ISELF,ISET,
     3                            KDERIV)
CFIXB

C-----  User-supplied function of the derivative of self-hardening
C     moduli

C-----  Use single precision on cray
C
           IMPLICIT REAL*8 (A-H,O-Z)
CFIXA
           DIMENSION GAMMA(NSLPTL), GMSLTL(NSLPTL), 
     2               NSLIP(NSET), PROP(16)
CFIXB

           IF (CHECK.EQ.0.) THEN

C-----  HYPER SECANT hardening law by Asaro, Pierce et al
              TERM1=PROP(1)*GAMTOL/(PROP(2)-PROP(3))
              TERM2=2.*EXP(-TERM1)/(1.+EXP(-2.*TERM1))
              TERM3=PROP(1)/(PROP(2)-PROP(3))*DSIGN(1.D0,GAMMA(KDERIV))
              DHSELF=-2.*PROP(1)*TERM2**2*TANH(TERM1)*TERM3

           ELSE

C-----  Bassani's hardening law
CFIXA
              TERM1=(PROP(1)-PROP(4))*GMSLTL(ISELF)/(PROP(2)-PROP(3))
CFIXB
              TERM2=2.*EXP(-TERM1)/(1.+EXP(-2.*TERM1))
              TERM3=(PROP(1)-PROP(4))/(PROP(2)-PROP(3))

              IF (KDERIV.EQ.ISELF) THEN
                 F=-2.*(PROP(1)-PROP(4))*TERM2**2*TANH(TERM1)*TERM3
                 ID=0
                 G=1.
                 DO I=1,NSET
                    IF (I.EQ.ISET) THEN
                       GAMMA0=PROP(5)
                       FAB=PROP(7)
                    ELSE
                       GAMMA0=PROP(6)
                       FAB=PROP(8)
                    END IF

                    DO J=1,NSLIP(I)
                       ID=ID+1
CFIXA
                       IF (ID.NE.ISELF) G=G+FAB*TANH(GMSLTL(ID)/GAMMA0)
CFIXB
                    END DO
                 END DO

              ELSE
                 F=(PROP(1)-PROP(4))*TERM2**2+PROP(4)
                 ILOWER=0
                 IUPPER=NSLIP(1)
                 IF (ISET.GT.1) THEN
                    DO K=2,ISET
                       ILOWER=ILOWER+NSLIP(K-1)
                       IUPPER=IUPPER+NSLIP(K)
                    END DO
                 END IF

                 IF (KDERIV.GT.ILOWER.AND.KDERIV.LE.IUPPER) THEN
                    GAMMA0=PROP(5)
                    FAB=PROP(7)
                 ELSE
                    GAMMA0=PROP(6)
                    FAB=PROP(8)
                 END IF

CFIXA
                 TERM4=GMSLTL(KDERIV)/GAMMA0
CFIXB
                 TERM5=2.*EXP(-TERM4)/(1.+EXP(-2.*TERM4))
                 G=FAB/GAMMA0*TERM5**2

              END IF

              DHSELF=F*G

           END IF

           RETURN
           END


C-----------------------------------


C-----  Use single precision on cray
CFIXA
           REAL*8 FUNCTION DHLATN(GAMMA,GMSLTL,GAMTOL,NSLPTL,NSET,
     2                            NSLIP,PROP,CHECK,ISELF,ISET,LATENT,
     3                            KDERIV)
CFIXB

C-----  User-supplied function of the derivative of latent-hardening 
C     moduli

C-----  Use single precision on cray
C
           IMPLICIT REAL*8 (A-H,O-Z)
CFIXA
           DIMENSION GAMMA(NSLPTL), GMSLTL(NSLPTL), NSLIP(NSET), 
     2               PROP(16)
CFIXB

           ILOWER=0
           IUPPER=NSLIP(1)
           IF (ISET.GT.1) THEN
              DO K=2,ISET
                 ILOWER=ILOWER+NSLIP(K-1)
                 IUPPER=IUPPER+NSLIP(K)
              END DO
           END IF

           IF (LATENT.GT.ILOWER.AND.LATENT.LE.IUPPER) THEN
              Q=PROP(9)
           ELSE
              Q=PROP(10)
           END IF

           IF (CHECK.EQ.0.) THEN

C-----  HYPER SECANT hardening law by Asaro, Pierce et al
              TERM1=PROP(1)*GAMTOL/(PROP(2)-PROP(3))
              TERM2=2.*EXP(-TERM1)/(1.+EXP(-2.*TERM1))
              TERM3=PROP(1)/(PROP(2)-PROP(3))*DSIGN(1.D0,GAMMA(KDERIV))
              DHLATN=-2.*PROP(1)*TERM2**2*TANH(TERM1)*TERM3*Q

           ELSE

C-----  Bassani's hardening law
CFIXA
              TERM1=(PROP(1)-PROP(4))*GMSLTL(ISELF)/(PROP(2)-PROP(3))
CFIXB
              TERM2=2.*EXP(-TERM1)/(1.+EXP(-2.*TERM1))
              TERM3=(PROP(1)-PROP(4))/(PROP(2)-PROP(3))

              IF (KDERIV.EQ.ISELF) THEN
                 F=-2.*(PROP(1)-PROP(4))*TERM2**2*TANH(TERM1)*TERM3
                 ID=0
                 G=1.
                 DO I=1,NSET
                    IF (I.EQ.ISET) THEN
                       GAMMA0=PROP(5)
                       FAB=PROP(7)
                    ELSE
                       GAMMA0=PROP(6)
                       FAB=PROP(8)
                    END IF

                    DO J=1,NSLIP(I)
                       ID=ID+1
CFIXA
                       IF (ID.NE.ISELF) G=G+FAB*TANH(GMSLTL(ID)/GAMMA0)
CFIXB
                    END DO
                 END DO

              ELSE
                 F=(PROP(1)-PROP(4))*TERM2**2+PROP(4)
                 ILOWER=0
                 IUPPER=NSLIP(1)
                 IF (ISET.GT.1) THEN
                    DO K=2,ISET
                       ILOWER=ILOWER+NSLIP(K-1)
                       IUPPER=IUPPER+NSLIP(K)
                    END DO
                 END IF

                 IF (KDERIV.GT.ILOWER.AND.KDERIV.LE.IUPPER) THEN
                    GAMMA0=PROP(5)
                    FAB=PROP(7)
                 ELSE
                    GAMMA0=PROP(6)
                    FAB=PROP(8)
                 END IF
CFIXA
                 TERM4=GMSLTL(KDERIV)/GAMMA0
CFIXB
                 TERM5=2.*EXP(-TERM4)/(1.+EXP(-2.*TERM4))
                 G=FAB/GAMMA0*TERM5**2

              END IF

              DHLATN=F*G*Q

           END IF

           RETURN
           END


C----------------------------------------------------------------------


      SUBROUTINE LUDCMP (A, N, NP, INDX, D)

C-----  LU decomposition

C-----  Use single precision on cray
C
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NMAX=200, TINY=1.0E-20)
      DIMENSION A(NP,NP), INDX(N), VV(NMAX)

      D=1.
      DO I=1,N
         AAMAX=0.

         DO J=1,N
            IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
         END DO

         IF (AAMAX.EQ.0.) PAUSE 'Singular matrix.'
         VV(I)=1./AAMAX
      END DO

      DO J=1,N
         DO I=1,J-1
            SUM=A(I,J)

            DO K=1,I-1
               SUM=SUM-A(I,K)*A(K,J)
            END DO

            A(I,J)=SUM
         END DO
         AAMAX=0.

         DO I=J,N
            SUM=A(I,J)

            DO K=1,J-1
               SUM=SUM-A(I,K)*A(K,J)
            END DO

            A(I,J)=SUM
            DUM=VV(I)*ABS(SUM)
            IF (DUM.GE.AAMAX) THEN
               IMAX=I
               AAMAX=DUM
            END IF
         END DO

         IF (J.NE.IMAX) THEN
            DO K=1,N
               DUM=A(IMAX,K)
               A(IMAX,K)=A(J,K)
               A(J,K)=DUM
            END DO

            D=-D
            VV(IMAX)=VV(J)
         END IF

         INDX(J)=IMAX
         IF (A(J,J).EQ.0.) A(J,J)=TINY
         IF (J.NE.N) THEN
            DUM=1./A(J,J)
            DO I=J+1,N
               A(I,J)=A(I,J)*DUM
            END DO
         END IF

      END DO

      RETURN
      END


C----------------------------------------------------------------------


      SUBROUTINE LUBKSB (A, N, NP, INDX, B)

C-----  Linear equation solver based on LU decomposition

C-----  Use single precision on cray
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NP,NP), INDX(N), B(N)

      II=0
      DO I=1,N
         LL=INDX(I)
         SUM=B(LL)
         B(LL)=B(I)

         IF (II.NE.0) THEN
            DO J=II,I-1
               SUM=SUM-A(I,J)*B(J)
            END DO
         ELSE IF (SUM.NE.0.) THEN
            II=I
         END IF

         B(I)=SUM
      END DO

      DO I=N,1,-1
         SUM=B(I)

         IF (I.LT.N) THEN
            DO J=I+1,N
               SUM=SUM-A(I,J)*B(J)
            END DO
         END IF

         B(I)=SUM/A(I,I)
      END DO

      RETURN
      END
