MODULE GTMP_2_GFCONVPAR
   USE module_mp_thompson_GF, only : thompson_init, mp_gt_driver
   
   PRIVATE
   PUBLIC GTMP_2_GFCONVPAR_interface
   Contains

   subroutine GTMP_2_GFCONVPAR_interface(its,ite,itf, kts,kte,ktf,cumulus &
                        ,ierr,klcl,kbcon,ktop,k22,kpbli            &
                        ,zo_cup,po_cup,tn_cup,qo_cup                 &
                        ,zo,tempco,qo,qco,qrco,hco,pwo,pwavo,rho       &
                        ,ztexec,zqexec,use_excess,zuo,up_massentr,up_massdetr,entr_rate_2d,cd &
                        ,dbyo,GAMMAo_cup,qeso_cup,vvel2d,vvel1d)
   IMPLICIT NONE
   real, parameter :: ctea=1./3. ,cteb=2., visc=2000., eps=0.622
   real, parameter :: f=2., C_d=0.506, gam=0.5, beta=1.875
   character *(*)              ,intent(in) ::  cumulus
   INTEGER                     ,iNTENT (IN)    :: its,ite, itf,kts,kte,ktf,use_excess
   integer, dimension (its:ite),intent (in)    :: kbcon,ktop,k22,klcl,kpbli
   integer, dimension (its:ite),intent (inOUT) :: ierr
   real,    dimension (its:ite),intent (in)    :: ztexec,zqexec
   real,    dimension (its:ite),intent (inOUT) :: pwavo,vvel1d
   real,    dimension (its:ite,kts:kte),intent (in)  ::                        &
                         zo_cup,po_cup,tn_cup,qo_cup                             &
                        ,hco,zo,qo,rho,zuo,up_massentr,up_massdetr,entr_rate_2d,cd &
                        ,dbyo,GAMMAo_cup,qeso_cup
                        
   real,    dimension (its:ite,kts:kte),intent (inOUT)  :: tempco,qco,qrco,pwo,vvel2d
    
   integer,PARAMETER :: mynum = 1

   REAL, DIMENSION( its:ite , kts:kte )  ::   temp     &
                                             ,dz       &
					     ,dz_cup   &
                                             ,p        &
                                             ,w                          

   REAL, DIMENSION( its:ite, kts:kte) ::  &
                  qv_curr , qc_curr , qr_curr, qi_curr, qs_curr, qg_curr &
                 ,qnc_curr, qnr_curr, qni_curr                           &
                 ,qnwfa_curr,qnifa_curr                           
                     
   REAL, DIMENSION( its:ite )  ::         RAINNC &
                                        ,RAINNCV &
                                         ,SNOWNC &
                                        ,SNOWNCV &
                                      ,GRAUPELNC &
                                     ,GRAUPELNCV &
                                             ,SR
                                                  
!----------------------------------------------------------------------
!--  qv         water vapor    mixing ratio (kg/kg)
!--  qc         cloud water    mixing ratio (kg/kg)
!--  qr         rain water     mixing ratio (kg/kg)
!--  qi         cloud ice      mixing ratio (kg/kg)
!--  qs         snow           mixing ratio (kg/kg)
!--  qg         graupel        mixing ratio (kg/kg)
!--  qh         hail           mixing ratio (kg/kg)
!--  qnc        droplet   number mixing ratio  (#/kg)
!--  qni        cloud ice number concentration (#/kg)
!--  qnr        rain      number concentration (#/kg)
!
!----------------------------------------------------------------------
!-- temp          temperature    (K)
!-- pi_phy        exner function           (dimensionless)
!-- p             pressure                 (Pa)
!-- RAINNC        grid scale precipitation (mm)
!-- RAINNCV       one time step grid scale precipitation (mm/step)
!-- SNOWNC        grid scale snow and ice (mm)
!-- SNOWNCV       one time step grid scale snow and ice (mm/step)
!-- GRAUPELNC     grid scale graupel (mm)
!-- GRAUPELNCV    one time step grid scale graupel (mm/step)
!-- HAILNC        grid scale hail (mm)
!-- HAILNCV       one time step grid scale hail (mm/step)
!-- SR            one time step mass ratio of snow to total precip
!-- z             Height above sea level   (m)
!-- dt            Time step              (s)
!-- G             acceleration due to gravity  (m/s^2)
!-- CP            heat capacity at constant pressure for dry air (J/kg/K)
!-- R_d           gas constant for dry air (J/kg/K)
!-- R_v           gas constant for water vapor (J/kg/K)
!-- XLS           latent heat of sublimation   (J/kg)
!-- XLV           latent heat of vaporization  (J/kg)
!-- XLF           latent heat of melting       (J/kg)
!-- precr         rain precipitation rate at all levels (kg/m2/s)
!-- preci         ice precipitation rate at all levels (kg/m2/s)
!-- precs         snow precipitation rate at all levels (kg/m2/s)
!-- precg         graupel precipitation rate at all levels (kg/m2/s)                             &

     REAL, DIMENSION(its:ite,kts:kte) :: HGT
     REAL, DIMENSION(its:ite,kts:kte) :: NWFA , NIFA,refl_10cm , rainprod,evapprod
     REAL, DIMENSION(its:ite)         :: qnwfa2d ,NWFA2D 
     REAL, DIMENSION(its:ite, kts:kte):: re_cloud, re_ice, re_snow
     REAL, DIMENSION(its:ite, kts:kte):: temp_env, qv_env,qnwfa_env,qnifa_env


     INTEGER :: itimestep = 1    !
     INTEGER :: do_radar_ref = 0 ! 
     INTEGER :: has_reqc= 1 , has_reqi= 1 , has_reqs= 1 
     REAL :: dx=10000., dy=10000.! grid spacing (m)
     REAL :: dt = 120.            ! model timestep (s)
     LOGICAL :: start_of_simulation =.true., &  
                diagflag=.false. 
     
     REAL, DIMENSION(kts:kte):: &
                          qv1d, qc1d, qi1d, qr1d, qs1d, qg1d, ni1d,     &
                          nr1d, nc1d, nwfa1d, nifa1d,                   &
                          t1d, p1d, w1d, dz1d, dBZ, pii1,th1

     integer, dimension (its:ite)   ::     start_level 
     REAL :: qaver
     INTEGER :: kr
!------------------------------------------------------------------
!-for gate soundings
     real, parameter ::                    &
        rgas     = 287.               &
    ,   cp       = 1004.              &
    ,   cv       = 717.               &
    ,   rm       = 461.               &
    ,   p00      = 1.e5               &
    ,   t00      = 273.16             &
    ,   g        = 9.80               &
    ,   cpor     = cp / rgas          &
    ,   rocv     = rgas / cv          
      integer, parameter :: gate=1
      integer, parameter :: klon = 2!161 ! number of soundings for gate
      integer, parameter :: klev = 41  ! must be = KTE,KTS, etc
      real, dimension(klon,klev,100):: vargrads
      character (len=50),dimension(100,2) :: gradsname
      integer :: nrec,nvx,nvar,klevgrads,nvar2d,nvar3d
      character (len=50) :: runname, runlabel
      namelist /run/ runname, runlabel
      integer ::jl,jk,nruns,k,kk,i,j
      integer, parameter :: levs=klev  

      OPEN(15,FILE='gf2.inp',STATUS='old',FORM='formatted')         
      read(15,nml=run)
      print*,"x=",runname,runlabel
      close(15)
      NRUNS = 2
 
      DO JL=1,NRUNS
        write(0,*) "#######################################",jl

        IF(JL==1) THEN

         open(10, FILE="cloud2.txt",STATUS='old',FORM='formatted')
          DO k=kts,kte

            read(10,101) kk, t1d(k) , th1(k), pii1(k),&
            p1d(k) , & 
            w1d(k) , &
            dz1d(k) ,&
            qv1d(k) ,&
            qc1d(k) ,&
            qi1d(k) ,&
            qr1d(k) ,&
            qs1d(k) ,&
            qg1d(k) ,&
            ni1d(k) ,&
            nr1d(k) ,&
            nc1d(k) ,&  
            nwfa1d(k) ,&
            nifa1d(k) 
            
            qv_curr(:,k)      = qv1d(k)         
            qc_curr(:,k)      = qc1d(k)         
            qr_curr(:,k)      = qr1d(k)         
            qi_curr(:,k)      = qi1d(k)         
            qs_curr(:,k)      = qs1d(k)         
            qg_curr(:,k)      = qg1d(k)         
            qni_curr(:,k)     = ni1d(k)         
            qnr_curr(:,k)     = nr1d(k)         
            temp(:,k)         =  t1d(k)
                p(:,k)        =  p1d(k)
                w(:,k)        =  w1d(k)


          ENDDO
        101 format(1x,i4,17e14.4)
        close(10)
ELSE
!GOTO 400
        DO i=its,itf
         if(ierr(i) /= 0) CYCLE

         DO k=kts,kte
           kr=k!+1
           qv_curr(i,k)      = qco(i,kr)!qv1d(k)         
           qc_curr(i,k)      = 0.       !qc1d(k)         
           qr_curr(i,k)      = 0.       !qr1d(k)         
           qi_curr(i,k)      = 0.       !qi1d(k)         
           qs_curr(i,k)      = 0.       !qs1d(k)         
           qg_curr(i,k)      = 0.       !qg1d(k)         
           qni_curr(i,k)     = 0.       !ni1d(kr)         
           qnr_curr(i,k)     = 0.       !nr1d(kr)         
               temp(i,k)     = tempco(i,kr)      ! t1d(k)
	   temp_env(i,k)     = tn_cup(i,kr)      ! t1d(k)
	     qv_env(i,k)     = qo_cup(i,kr)      ! t1d(k)
                  p(i,k)     = po_cup(i,kr)*100. ! p1d(k)
                  w(i,k)     = vvel2d(i,kr)      ! w1d(k)
                hgt(i,k)     = zo(i,kr)
         print*,"ENV1=",k,qco(i,kr),qo_cup(i,kr)  ,tempco(i,kr) ,tn_cup(i,kr) 		
         ENDDO
	 
         DO k=kts,kte-1 
	       kr=k!+1
               dz    (i,k) = zo    (i,kr+1)-zo    (i,kr)
	       dz_cup(i,k) = zo_cup(i,kr+1)-zo_cup(i,kr)
         ENDDO
	 dz(i,kte)         = dz(i,kte-1)
	 dz_cup(i,kte)     = dz_cup(i,kte-1)
         !
	 !-- if aerosol-aware
 	 DO k=kts,kte
           kr=k!+1
           qnc_curr  (i,k)   =0.   ! NC=qnc_curr,     
           qnwfa_curr(i,k)   =0.   ! NWFA=qnwfa_curr, air-parcel
           qnifa_curr(i,k)   =0.   ! NIFA=qnifa_curr, air-parcel
           qnwfa_env (i,k)   =0.  ! NWFA=qnwfa_curr, environ
           qnifa_env (i,k)   =0.  ! NIFA=qnifa_curr, environ
         ENDDO
         qnwfa2d   (i)=0.0      
	 !
	 !--output----------------------
         SR        (i)=0.0   
         RAINNC    (i)=0.0        
         RAINNCV   (i)=0.0       
         SNOWNC    (i)=0.0        
         SNOWNCV   (i)=0.0       
         GRAUPELNC (i)=0.0     
         GRAUPELNCV(i)=0.0    
         rainprod  (i,kts:kte) = 0. 
         evapprod  (i,kts:kte) = 0. 
        
         if (diagflag .and. do_radar_ref == 1) refl_10cm(i,:)=0.

         if (has_reqc .ne.0 .and. has_reqi.ne.0 .and. has_reqs.ne.0) then
           re_cloud(i,:)     = 0.
           re_ice  (i,:)     = 0.
           re_snow (i,:)     = 0.
         endif
        ENDDO
400 continue

        IF(start_of_simulation) then 
            NWFA=0.0;NIFA=0.0;NWFA2D=0.0
            CALL thompson_init(hgt,dx, dy, start_of_simulation,its, ite, kts, kte,   &
                               NWFA2D,NWFA, NIFA)
            qnwfa_curr=NWFA
            qnifa_curr=NIFA
            qnwfa2d   =0.
!---temporary comment for GATE soundings (uncomment for normal runs)
!start_of_simulation =.false.
!----------
        ENDIF

        CALL mp_gt_driver(mynum,                &
                     qv_curr,                   &! QV=qv_curr,     
                     qc_curr,                   &! QC=qc_curr,     
                     qr_curr,                   &! QR=qr_curr,     
                     qi_curr,                   &! QI=qi_curr,     
                     qs_curr,                   &! QS=qs_curr,     
                     qg_curr,                   &! QG=qg_curr,     
                     qni_curr,                  &! NI=qni_curr,    
                     qnr_curr,                  &! NR=qnr_curr,    
                     TEMP,                      &!  temperature    (K)
                     P,                         &! pressure(Pa)
                     W,                         &
                     dz,                      &
                     dt,                        &! time step              (s)
                     itimestep,                 &
                     RAINNC,                    &
                     RAINNCV,                   &
                     SNOWNC,                    &
                     SNOWNCV,                   &
                     GRAUPELNC,                 & 
                     GRAUPELNCV,                & 
                     SR,                        &
                     rainprod,                  &
                     evapprod,                  &
                     refl_10cm,                 &
                     diagflag,                  &
                     do_radar_ref,              &
                     re_cloud,                  & 
                     re_ice,                    &
                     re_snow,                   &
                     has_reqc,                  &
                     has_reqi,                  &
                     has_reqs,                  &
                     ITS,ITE, KTS,KTE,KTF,      &! replace KTE by KTF
!
!- environmental state for cumulus parameterization 
                     ierr,kbcon,ktop,klcl,kpbli,&
		     temp_env,                  &
		     qv_env  ,                  & 
                     zuo,up_massentr,up_massdetr,entr_rate_2d,cd,dz_cup, &
!- optional arrays
		     qnwfa_env,                 &
                     qnifa_env,                 &
!- environmental state for cumulus parameterization 
!
!- moving optional arrays to the end of the argument list
                     qnc_curr,                  &! NC=qnc_curr,     
                     qnwfa_curr,                &! NWFA=qnwfa_curr, 
                     qnifa_curr,                &! NIFA=qnifa_curr, 
                     qnwfa2d                    &! NWFA2D=qnwfa2d,  
                                                )
RETURN !<<<<<<< 
        ENDIF
GOTO 500
        DO i=its,itf
         if(ierr(i) /= 0) CYCLE
         DO k=kts,kte
          kr=k+1
          qco (i,kr)   =qv_curr(i,k) !in-cloud water vapor mixing ratio (kg/kg)
          qrco(i,kr)   =qc_curr(i,k) !in-cloud liquid water mixing ratio (kg/kg)     
          !?           =qr_curr(i,k) !rainfall mixing ratio (kg/kg)      
          !?           =qi_curr(i,k) !in-cloud ice  mixing ratio (kg/kg)     
          !?           =qs_curr(i,k) !snow mixing ratio (kg/kg)          
          !?           =qg_curr(i,k) !graupel mixing ratio (kg/kg)        
          !?           =qni_curr(i,k)!in-cloud ice number concentration (#/m3)    
          !?           =qnr_curr(i,k)!in-cloud rainfall number concentration (#/m3)        
          tempco(i,kr) =temp(i,k)    !in-cloud air temperature (K) 
          vvel2d(i,kr) =  w(i,k)     !in-cloud updraft vert velocity(m/s)
          !pwo(i,kr) = ??
         ENDDO
         !pwavo(i)=??

        ENDDO
500 continue
!
!-SRF - OUTPUT================================================================= 
       IF(gate==1) THEN
         do jk=1,kte
           vargrads(jl,jk,1) = qv_curr   (1,jk)*1000.
           vargrads(jl,jk,2) = qc_curr   (1,jk)*1000.
           vargrads(jl,jk,3) = qr_curr   (1,jk)*1000.
           vargrads(jl,jk,4) = qi_curr   (1,jk)*1000.
           vargrads(jl,jk,5) = qs_curr   (1,jk)*1000.
           vargrads(jl,jk,6) = qg_curr   (1,jk)*1000.
           vargrads(jl,jk,7) = qni_curr  (1,jk)
           vargrads(jl,jk,8) = qnr_curr  (1,jk)
           vargrads(jl,jk,9) = qnc_curr  (1,jk)
           vargrads(jl,jk,10)= qnwfa_curr(1,jk)
           vargrads(jl,jk,11)= qnifa_curr(1,jk)
           vargrads(jl,jk,12)= NWFA(1,jk)
           vargrads(jl,jk,13)= NIFA(1,jk)
           vargrads(jl,jk,14)= w(1,jk)
           vargrads(jl,jk,15)= t1d(jk)
           vargrads(jl,jk,16)= temp(1,jk)
          enddo
          nvar3d=16
        !- surface quantities
         vargrads(jl,1,21) = RAINNC(1)*3600.
         vargrads(jl,1,22) = RAINNCV(1)
         vargrads(jl,1,23) = -999.
         vargrads(jl,1,24) = -999.
         vargrads(jl,1,25) = -999.
         vargrads(jl,1,26) = -999.
         vargrads(jl,1,27) = -999.
         nvar2d=7
         !total number of variables out
         nvar = nvar3d+nvar2d ; if(nvar>100) stop "increase nvar"
       ENDIF  
           
 ENDDO  ! loop on nruns (number of GATE Soundings) 
 !
 IF(gate==1) THEN
    IF(gate==1) THEN
      PRINT*,'Writing GrADS control file:',runname//'.ctl'
      gradsname(1,1)  ='qv     '         ;gradsname(1,2)  =  'qv          '
      gradsname(2,1)  ='qc     '         ;gradsname(2,2)  =  'qc          '
      gradsname(3,1)  ='qr     '         ;gradsname(3,2)  =  'qr          '
      gradsname(4,1)  ='qi     '         ;gradsname(4,2)  =  'qi          '
      gradsname(5,1)  ='qs     '         ;gradsname(5,2)  =  'qs          '
      gradsname(6,1)  ='qg     '         ;gradsname(6,2)  =  'qg          '
      gradsname(7,1)  ='qni    '         ;gradsname(7,2)  =  'qni    '
      gradsname(8,1)  ='qnr    '         ;gradsname(8,2)  =  'qnr    '
      gradsname(9,1)  ='qnc    '         ;gradsname(9,2)  =  'qnc    '
      gradsname(10,1) ='qnwfa  '         ;gradsname(10,2) =  'qnwfa  '
      gradsname(11,1) ='qnifa  '         ;gradsname(11,2) =  'qnifa  '
      gradsname(12,1) ='nwfa  '         ;gradsname(12,2) =  'zero  '
      gradsname(13,1) ='nifa  '         ;gradsname(13,2) =  'zero  '
      gradsname(14,1) ='W  '         ;gradsname(14,2) =  'W  '
      gradsname(15,1) ='T  '         ;gradsname(15,2) =  'T   '
      gradsname(16,1) ='T2  '         ;gradsname(16,2) =  'T2  '
! 
      
      gradsname(21,1) ='rainnc'    ;gradsname(21,2) =  'liquid surf prec [mm/h]'
      gradsname(22,1) ='rainncv'   ;gradsname(22,2) =  'rainncv'
      gradsname(23,1) ='ierr  '    ;gradsname(23,2) =  'IERR'
      gradsname(24,1) ='cltop '    ;gradsname(24,2) =  'cloud top index'
      gradsname(25,1) ='clbas '    ;gradsname(25,2) =  'cloud base indec'
      gradsname(26,1) ='dnmf  '    ;gradsname(26,2) =   'downdraft mass flux at initiation level'
!     gradsname(27,1) ='pcape  '   ;gradsname(27,2) =  'pcape'
      gradsname(27,1) ='pcin   '   ;gradsname(27,2) =  'not-working : pcin'
      
       
!       OPEN(20,file=trim(runname)//'.ctl',status='unknown')
       OPEN(20,file='x.ctl',status='unknown')
       write(20,2001) '^'//trim(runname)//'.gra'
       write(20,2002) 'undef -9.99e33'
       write(20,2002) 'options byteswapped' ! zrev'
       !write(20,2002) 'title '//"GF with GATE soundings"
       write(20,2002) 'title '//trim(runlabel)
       write(20,2003) 1,0.,1. ! units m/km
       write(20,2004) klon,1.,1.
       write(20,2005) klev,(0.01*p1d(jk),jk=1,klev)
       write(20,2006) 1,'00:00Z01JAN2000','1mn'
       write(20,2007) nvar3d+nvar2d
       do nvx=1,nvar3d         
         klevgrads=klev
         write(20,2008) gradsname(nvx,1),klevgrads,gradsname(nvx,2)
       enddo
       do nvx=21,20+nvar2d         
         klevgrads=0
         write(20,2008) gradsname(nvx,1),klevgrads,gradsname(nvx,2)
       enddo
       write(20,2002) 'endvars'
       
  2001 format('dset ',a)
  2002 format(a)
  2003 format('xdef ',i4,' linear ',2f15.3)
  2004 format('ydef ',i4,' linear ',2f15.3)
  2005 format('zdef ',i4,' levels ',60f6.0)
  2006 format('tdef ',i4,' linear ',2a15)
  2007 format('vars ',i4)
  2008 format(a10,i4,' 99 ',a40)!'[',a8,']')
! 2008 format(a10,i4,' 99' )
  2055 format(60f7.0)
   133 format (1x,F7.0)
   CLOSE(20)
  
   print*, 'opening GrADS file:',trim(runname)//'.gra'
   OPEN(19,FILE= trim(runname)//'.gra',  &
   FORM='unformatted',ACCESS='direct',STATUS='unknown', RECL=(klon)) !INTEL
   !FORM='unformatted',ACCESS='direct',STATUS='unknown', RECL=4*(klon))!PGI
   NREC=0
   do nvx=1,nvar3d
      klevgrads=klev   
      do jk=1,klevgrads
        nrec=nrec+1
        WRITE(19,REC=nrec) real(vargrads(:,jk,nvx),4)
      enddo    
   enddo
   do nvx=21,20+nvar2d
      klevgrads=1
      do jk=1,klevgrads
        nrec=nrec+1
        WRITE(19,REC=nrec) real(vargrads(:,jk,nvx),4)
      enddo    
   enddo


   close (19)
  ENDIF
!SRF - end

!-- formats
167     format(A12,1x,i2,1x,f7.0,2(1x,e13.4),2(1x,f6.1))
121     format(i3,2x,f8.2,2x,f8.2,2x,e13.6,2x,f10.2,2x,f10.2)
122     format(f5.1,2x,f8.2,2x,f6.2)
124     format(f8.2,11(2x,e14.6))
123     format(i3,2x,i3,2x,e13.6)
125     format(i2,2x,f8.2,11(2x,e13.6))
132     format(1x,i2,f8.0,1x,2(1x,f8.3),5(1x,e12.4))

 ENDIF !gate==1
!-SRF - OUTPUT================================================================= 
end SUBROUTINE GTMP_2_GFCONVPAR_interface
END MODULE GTMP_2_GFCONVPAR
