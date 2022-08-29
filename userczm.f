c Copyright Giovanni Florian
*deck,userCZM     parallel       optimize  USERDISTRIB               gal
      subroutine userCZM (matId, elemId, kMatIntPt, ldstep,isubst,
     &                    keycut, ncomp,nProp, nstatev,
     &                    Time, dTime, Temp, dTemp,
     &                    coords, prop, Strain, dStrain, 
     &                    stress, dsdePl, sedEl, sedPl, statev,
     &                    var1, var2, var3, var4, var5)
c
c*************************************************************************
c
c     *** primary function ***
c
c           user cohesive zone model example
c
c           Commands
c             TB,CZM,mat,NTEMP,NPTS,user 
c                TBTEMP if mat. constants are temperature dependent
c                TBDATA define material constants
c
c*************************************************************************
c
c     input arguments
c     ===============
c      matId     (int,sc,in)              material #
c      elemId    (int,sc,in)              element #
c      kMatIntPt (int,sc,in)              material integration point #
c      ldstep    (int,sc,in)              load step number
c      isubst    (int,sc,in)              substep number
c      ncomp     (int,sc,in)              number of stress, strain components
c      nProp     (int,sc,in)              Number of material ocnstants
c      nstatev   (int,sc,in)              Number of state variables
c
c      Temp      (dp ,sc,in)              temperature at beginning of time increment
c      dTemp     (dp ,sc,in)              temperature increment 
c      Time      (dp ,sc,in)              time at beginning of increment (t)
c      dTime     (dp ,sc,in)              time increment (dt)
c
c      prop     (dp,ar(nprop),i)          Material constants defined by TB command 
c      Strain   (dp,ar(ncomp),i)          Interface separation at beginning of time increment
c      dStrain  (dp,ar(ncomp),i)          Interface separation increment
c      coords   (dp,ar(3),i)              current coordinates
c
c     output arguments              
c     ======================             
c      stress   (dp,ar(nTesn),io)         Traction stress
c      sedEl    (dp,sc,io)                elastic work
c      sedPl    (dp,sc,io)                plastic work
c      keycut   (int,sc,io)               loading bisect/cut control
c                                         0 - no bisect/cut
c                                         1 - bisect/cut 
c                                         (factor will be determined by ANSYS solution control)
c      dsdePl   (dp,ar(ncomp,ncomp),io)   consistent tangent jacobian matrix
c
c     input output arguments              
c     ======================             
c      statev   (dp,ar(nstatev,io)        user defined solution state variables
c
c     misc.
c     ======================             
c      var1, var2, var3, var4, var5       currently not used
c
c     local variables
c     ======================             
c
c      debugflag (in,sc, l)                debugflag to print debug information
c
c
c*************************************************************************
c
#include "impcom.inc"
c
      INTEGER          
     &                 matId,elemId,kMatIntPt,ldstep,isubst,keycut,
     &                 ncomp,nProp,nstatev
      DOUBLE PRECISION 
     &                 Time,dTime,Temp,dTemp,sedEl,sedPl

      DOUBLE PRECISION coords(3),prop(nProp),Disp(ncomp), 
     &                 Strain(ncomp),dStrain(ncomp),stress(ncomp), 
     &                 dsdePl(ncomp,ncomp),statev(nstatev), 
     &                 DSEC(ncomp,ncomp),PenK(nComp)
      DOUBLE PRECISION var1, var2, var3, var4, var5
c
c --- parameters
c
      DOUBLE PRECISION ZERO
      PARAMETER       (ZERO       = 0.d0
     &                 )

      INTEGER          debugflag

      DOUBLE PRECISION toler
      PARAMETER (toler = 1.0d-15)      
      
      DOUBLE PRECISION stiffi, sigm, deltaN, deltaT, deltaN1, 
     &                 phiN, dsigCap, kdsigCap, sigRef,
     &                 deltaBN, deltaBX,deltaBY, deltaBT,deltaBT1,
     &                 kdeltaT,qfact, rfact,
     &                 DtDx,DtDy,DtDx2,DtDy2,DnDn,DnDt,DtDn,DtDt,
     &                 con1,con2,con3,con4,tmp1,tmp2,tmp3,sn,st
      INTEGER          i,j
c
      DOUBLE PRECISION GIC,GIIC,Sigc,Tauc,Pen,Eta,IDTFLAG,RT_Old
c
c --- debug includes
#include "locknm.inc"
      INTEGER          wrinqr, iott
      EXTERNAL         pplock, ppunlock, wrinqr
c
c*************************************************************************
c
c *** initialization
c     debugflag
      debugflag = 1
      keycut = 0
      sedEl    = ZERO
      sedPl    = ZERO
* ---------------------------------------------------------------------*
C   Initialization of variables
* ---------------------------------------------------------------------*
      dsdePl = 0.d0
      Stress = 0.d0
*----------------------------------------------------------------------*
C   Properties of cohesive element
*----------------------------------------------------------------------*
      GIC     = prop(1) ! Fracture toughness, mode I
      GIIC    = prop(2) ! Fracture toughness, mode II
      Sigc    = prop(3) ! Strength, mode I
      Tauc    = prop(4) ! Strength, mode II
      Pen     = prop(5) ! Penalty stiffness, mode I
      Eta     = prop(6) ! B-K exponent
*----------------------------------------------------------------------*
C   Compute total displacements
*----------------------------------------------------------------------*
      Disp = Strain !+ dstrain ! does not converge
C
      IDTFLAG = 0
c
      CALL Damage(ncomp,Disp,stress,Pen,Sigc,Tauc,GIC,GIIC,Eta,
     &        statev,PenK,IDTFLAG,RT_Old,nstatev)
* ---------------------------------------------------------------------*
C   Calculation of the consistent tangent stiffness matrix
* ---------------------------------------------------------------------*
      CALL STIFFNUM(ncomp,Disp,stress,Pen,Sigc,Tauc,GIC,GIIC,Eta,
     &        statev,PenK,dsdePl,dStrain,RT_Old,nstatev)
C
      if(debugflag .gt. 0) then
         call pplock(LOCKOT)
         iott = wrinqr(2)
         write(iott,*) 'userCZM debug :',
     &   ' isubst=',isubst!,
c     &   ' coords=',coords,
c     &   ' ncomp=',ncomp,
c     &   ' Penk=',Penk
c         write(iott,*) 'prop      :',prop
c         write(iott,*) 'statev    :',statev
c         write(iott,*) 'dStrain    :',dstrain
c         write(iott,*) 'Disp    :',Disp
c         write(iott,*) 'dsdePl(i,i):',
c     &    (dsdePl(i,i),i=1,ncomp)
c         write(iott,*) 'dsec(i,i):',
c     &    (dsec(i,i),i=1,ncomp)
c         write(iott,*)
c     &                 'statev    :',(statev(i),i=1,nstatev)
         call ppunlock(LOCKOT)
c 1000    format(a/4(a,i5,1x),4(a,i1,1x)/5x,7(a,e12.5,1x),a,3e12.4)
      end if

      RETURN
      END SUBROUTINE userCZM
C
*************************************************************************
**                            SUBROUTINE Damage                        **
*************************************************************************
      SUBROUTINE Damage(ncomp,Disp,stress,Pen,Sigc,Tauc,GIC,GIIC,Eta,
     &            statev,PenK,IDTFLAG,RT_Old,nstatev)
C
#include "impcom.inc"
C
      INTEGER          ncomp,nstatev,iott

      DOUBLE PRECISION stress(ncomp),statev(nstatev),PenK(nComp),
     &                 Disp(ncomp),DSEC(ncomp,ncomp)
      
      DOUBLE PRECISION ZERO, ONE, TWO
      PARAMETER       (ZERO       = 0.d0,
     &                 ONE        = 1.d0,
     &                 TWO        = 2.d0
     &                 )
c
      DOUBLE PRECISION RT,RT_Old,IDTFLAG,Pen,GIC,GIIC,Tauc,Sigc,Dx,DI,
     &                 DII,Beta,DL,PenB,XLo,XLf,R,Damg,DelDam,Dy,eta,
     &                 DIma,KDI,KDII,KDIII,K2DI,K2DII,K2DIII
      INTEGER          I
c        
      RT    = statev(1)
      IF (IDTFLAG.NE.0)  RT = RT_Old
c 
      PenK(1) = Pen   ! Pen :: mode I interface penalty stiffness
      DO I=2,ncomp        
        PenK(I) = Pen*(GIC/GIIC)*(Tauc/Sigc)**2 ! Turon constraint
      ENDDO
* ---------------------------------------------------------------------
C   Solution dependent state variable Damg: damage at the end of the
C   last converged increment
* ---------------------------------------------------------------------
      if (ncomp .eq. 3) then
       DI   = Disp(1)
       Dx   = Disp(2)
       Dy   = Disp(3)
       DII  = SQRT(Dx*Dx + Dy*Dy)
      else
       DI   = Disp(1)
       Dx   = Disp(2)
       DII  = SQRT(Dx*Dx)
      end if
* ---------------------------------------------------------------------
c   determine Macaulay displacement and Kii*Deltai
* ---------------------------------------------------------------------
      DIma = (ONE/TWO)*(DI + abs(DI))
c
      if (ncomp .eq. 3) then
        KDI    = Penk(1)*DIma**2
        KDII   = Penk(2)*Dx**2
        KDIII  = Penk(3)*Dy**2
        K2DI   = (Penk(1)*DIma)**2
        K2DII  = (Penk(2)*Dx)**2
        K2DIII = (Penk(3)*Dy)**2
      else
        KDI    = Penk(1)*DIma**2
        KDII   = Penk(2)*Dx**2
        K2DI   = (Penk(1)*DIma)**2
        K2DII  = (Penk(2)*Dx)**2
      end if
* ---------------------------------------------------------------------
c   determine mixed-mode ratios
* ---------------------------------------------------------------------
      If (DI.LT.1.D-19) THEN ! mode II
        Beta = ONE
        DL   = DII
        PenB = PenK(2)
      Else
        if (ncomp .eq. 3) then
          Beta = (KDII + KDIII) / (KDI + KDII + KDIII)
          DL   = (KDI + KDII + KDIII) / sqrt(K2DI + K2DII + K2DIII)
        else
          Beta = (KDII) / (KDI + KDII)
          DL   = (KDI + KDII) / sqrt(K2DI + K2DII)
        end if
          PenB = Penk(1)*(1-Beta) + Beta*Penk(2)
      EndIf
* ---------------------------------------------------------------------
c    determine mixed-mode onset and final displacement B-K criterion
* ---------------------------------------------------------------------
      XLo = sqrt((Penk(1)*(Sigc/Pen)**2 + (Penk(2)*(Tauc/Pen)**2 
     &  - Penk(1)*(Sigc/Pen)**2)*(Beta**Eta))/(PenB))
c
      XLf = (Penk(1)*2*GIC/Pen + (Penk(2)*2*GIIC/Pen -
     &  Penk(1)*2*GIC/Pen)*(Beta**Eta))/(PenB*XLo)
* ---------------------------------------------------------------------
c    calculate mixed-mode damage threshold
* ---------------------------------------------------------------------
      R = (DL - XLo) / (XLf - XLo)
      If (R.LT.Zero) R=Zero
      If (R.GT.One) R=One
* ---------------------------------------------------------------------
c    update internal variables
* ---------------------------------------------------------------------
      if (RT.lt.R) then
        RT=R
      end if
c
      Damg = RT*XLf/(RT*XLf+(One-RT)*XLo)
      write(iott,*) 'DL  :',DL
      write(iott,*) 'XLo :',XLo
      write(iott,*) 'XLf :',XLf
* ---------------------------------------------------------------------*
C   Compute TAU tensor (tractions)
* ---------------------------------------------------------------------*
      DSEC = 0.D0
      Do I=1,ncomp
        DSEC(I,I)=(1.D0-Damg)*PenK(I)
      EndDO
      If (DI.lt.0.D0) then
        DSEC(1,1)= PenK(1) ! interpenetration
      end if
      Do I=1,ncomp
        stress(I) = DSEC(I,I)*Disp(I)
      EndDo
* ---------------------------------------------------------------------
C   Update state variables
C   (if Damage is called from UEL(IDTFLAG=0) and not from STIFF(IDTFLAG=1)
* ---------------------------------------------------------------------
      IF (IDTFLAG.EQ.0) THEN ! Save State variables only on call from UEL
        RT_Old = STATEV(1)
        DelDam = RT - RT_Old
        STATEV(1) = RT
        STATEV(2) = Beta     ! mixed-mode ratio
      END IF
      RETURN
      END SUBROUTINE Damage
*************************************************************************
**                        SUBROUTINE STIFFNUM                          **
*************************************************************************
      SUBROUTINE STIFFNUM(ncomp,Disp,stress,Pen,Sigc,Tauc,GIC,GIIC,Eta,
     &           statev,PenK,dsdePl,dStrain,RT_Old,nstatev)
C
#include "impcom.inc"
C
      DOUBLE PRECISION dsdePl(ncomp,ncomp),statev(nstatev),
     &                 Disp(ncomp),stresspert(ncomp),
     &                 DispPERT(ncomp),PenK(ncomp),
     &                 dStrain(ncomp),stress(ncomp)
      DOUBLE PRECISION Pen,GIC,GIIC,Tauc,Sigc,Eta,RT_Old,
     &                   IDTFLAG,SIGNPERT,PERT
      INTEGER          IPERT,ISTRE,ncomp,nstatev
C
* ---------------------------------------------------------------------
C      TANGENT STIFFNESS MATRIX   DTAN
* ---------------------------------------------------------------------
      PERT=1.D-8 ! perturbation constant
C
      IDTFLAG=1
      DO IPERT=1,ncomp
      DispPERT = Disp
C
      IF (abs(dStrain(IPERT)).GT.1.D-8) THEN
        SIGNPERT = dStrain(IPERT)/abs(dStrain(IPERT))
      ELSE
        SIGNPERT = 1.D0
      ENDIF
C
      DispPERT(IPERT) = Disp(IPERT) + SIGNPERT*PERT
C
      CALL Damage(ncomp,DispPERT,stresspert,Pen,Sigc,Tauc,
     &     GIc,GIIc,Eta,statev,PenK,IDTFLAG,RT_Old,nstatev)
c
      DO ISTRE=1,ncomp   
        dsdePl(ISTRE,IPERT) = (stresspert(ISTRE) - stress(ISTRE))
     &  /(PERT*SIGNPERT)
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE STIFFNUM
