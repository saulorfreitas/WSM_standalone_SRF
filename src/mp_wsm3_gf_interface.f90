 USE GTMP_2_GFCONVPAR, only : GTMP_2_GFCONVPAR_interface



ELSEIF(CLOUDMP==2) then
 print*,"z1=",zo(1,:)
 print*,"z2=",zo_cup(1,:)

       do i=its,itf
          qco(i,:) = qo_cup(i,:)
          tempco(i,:) = t_cup(i,:)
          vvel2d(i,:) = 0.
          if(ierr(i) /= 0) cycle
          start_level(i)=klcl(i)
          x_add =float(use_excess)*zqexec(i)
          call get_cloud_bc(kts,kte,qo_cup (i,kts:kte),qco(i,kts),k22(i),x_add)

          do k=kts,start_level(i)
            qco(i,k) = qco(i,kts) !- air parcel water vapor
            tempco(i,k) = (1./cp)*(hco(i,k)-g*zo_cup(i,k)-xlv*qco(i,k))
            vvel2d(i,k) = zws(i)
            !---- temporary !<<<<<<<<<<<<<<<<<<<<<
            print*,"Ks=",ktop(i),kbcon(i),start_level(i)
            print*,"CL=",k, qco(i,k),tempco(i,k),hco(i,k),zo_cup(i,k)
            !---- temporary !<<<<<<<<<<<<<<<<<<<<<
          enddo

         enddo
         call GTMP_2_GFCONVPAR_interface(its,ite,itf, kts,kte,ktf
             ,cumulus &
             ,ierr,klcl,kbcon,ktop,k22,kpbl             &
             ,zo_cup,po_cup,tn_cup,qo_cup                 &
             ,zo,tempco,qo,qco,qrco,hco,pwo,pwavo,rho       &
             ,ztexec,zqexec,use_excess
             ,zuo,up_massentr,up_massdetr,entr_rate_2d,cd &
             ,dbyo,GAMMAo_cup,qeso_cup,vvel2d,vvel1d)
         
ENDIF



      call cup_up_moisture(cumulus,start_level,klcl,ierr,ierrc,zo_cup,qco,qrco,pwo,pwavo,hco,tempco,xland   &
         ,po,p_cup,kbcon,ktop,cd,dbyo,clw_all,t_cup,qo,GAMMAo_cup,zuo,qeso_cup &
         ,k22,qo_cup,ZQEXEC,use_excess,ccn,rho,up_massentr,up_massdetr,psum    &
         ,psumh,c1d,x_add_buoy,vvel2d,vvel1d,zws,entr_rate_2d                  &
         ,1,itf,ktf,ipr,jpr,its,ite, kts,kte,z1                                )


!======================================================

   subroutine cup_up_moisture(name,start_level,klcl,ierr,ierrc,z_cup,qc,qrc,pw,pwav,hc,tempc,xland,&
      po,p_cup,kbcon,ktop,cd,dby,clw_all,                   &
      t_cup,q,gamma_cup,zu,qes_cup,k22,qe_cup,           &
      zqexec,use_excess,ccn,rho,                            &
      up_massentr,up_massdetr,psum,psumh,c1d,x_add_buoy,  &
      vvel2d,vvel1d,zws,entr_rate_2d,                         &
      itest,itf,ktf,ipr,jpr,its,ite, kts,kte,z1      )

      implicit none
      real, parameter :: bdispm = 0.366       !berry--size dispersion (maritime)
      real, parameter :: bdispc = 0.146       !berry--size dispersion (continental)
      real, parameter :: T_BF   = 268.16 , T_ice_BF = 235.16
      real, parameter :: rk = 3 ! or 2
      real, parameter :: xexp = 2.
      real, parameter :: FRACT = 1.
      !
      !  on input
      integer  ,intent (in   ) ::  use_excess,itest,itf,ktf    &
         ,its,ite,ipr,jpr, kts,kte
      ! cd= detrainment function
      ! q = environmental q on model levels
      ! qe_cup = environmental q on model cloud levels
      ! qes_cup = saturation q on model cloud levels
      ! dby = buoancy term
      ! cd= detrainment function
      ! zu = normalized updraft mass flux
      ! gamma_cup = gamma on model cloud levels
      !
      character *(*)                    ,intent (in) ::  name
      integer, dimension (its:ite)      ,intent (in) ::  kbcon,ktop,k22,klcl,start_level
      real,  dimension (its:ite,kts:kte),intent (in) ::  t_cup,p_cup,rho,q,zu,gamma_cup       &
         ,qe_cup,hc,po,up_massentr,up_massdetr &
         ,dby,qes_cup,z_cup,cd,c1d

      real,  dimension (its:ite)        ,intent (in) ::  zqexec,xland,x_add_buoy
      real,  dimension (its:ite)        ,intent (in) ::  zws,ccn,z1
      real,  dimension (its:ite,kts:kte),intent (in) ::  entr_rate_2d
      real,  dimension (its:ite,kts:kte),intent (in) ::  vvel2d
      real,  dimension (its:ite        ),intent (in) ::  vvel1d
      !
      ! input and output
      !
      ! ierr error value, maybe modified in this routine
      integer, dimension (its:ite)  ,intent (inout)                   ::  ierr
      ! qc = cloud q (including liquid water) after entrainment
      ! qrch = saturation q in cloud
      ! qrc = liquid water content in cloud after rainout
      ! pw = condensate that will fall out at that level
      ! pwav = totan normalized integrated condensate (I1)
      ! c0 = conversion rate (cloud to rain)

      real,   dimension (its:ite,kts:kte),intent (out)   :: qc,qrc,pw,clw_all,tempc
      real,   dimension (its:ite)        ,intent (out)   :: pwav,psum,psumh
      character*128                      ,intent (inout) :: ierrc(its:ite)
      !
      !  local variables in this routine
      !
      integer                              ::                           &
         iounit,iprop,i,k,k1,k2,n,nsteps
      real                                 ::                           &
         dp,rhoc,dh,qrch,dz,radius,berryc0,q1,berryc
      real :: qaver,denom,aux,cx0,qrci,step,cbf,qrc_crit_BF,min_liq,qavail
      real delt,tem1,qrc_0,cup,hei,qrc_crit_yin,dz1m,vs,Nc,w_upd,Faut
      integer, parameter :: n_smooth=1
      real, parameter :: qrc_crit_yin_ref=-500.*log(2000./9500)*1.e-6/0.93
      real, parameter :: peaut = .55        ! collection efficiency
      real, parameter :: r0 = .8e-5         ! 8 microm  in contrast to 10 micro m
      real, parameter :: xmyu = 1.718e-5    ! the dynamic viscosity kgm-1s-1
      real, parameter :: pi = 4.*atan(1.)
      real, parameter :: rhowater     = 1000.      ! density of liquid water at 0^oc (kg m^-3)
      real, parameter :: rhosnow      = 100.       ! density of snow (kg m^-3)
      real, parameter :: rhoair0      = 1.28       ! density of dry air at 0^oc and 1000mb pressure (kg m^-3)
      !real,   dimension (its:ite,kts:kte) :: qc2

      !--- no precip for small clouds
      !if(name.eq.'shallow')  c0 = c0_shal
      !if(name.eq.'mid'    )  c0 = c0_mid
      !if(name.eq.'deep'   )  c0 = c0_deep
      do i=its,itf
         pwav (i)=0.
         psum (i)=0.
         psumh(i)=0.
      enddo
      do k=kts,ktf
         do i=its,itf
            pw      (i,k)=0.
            clw_all (i,k)=0.
            tempc   (i,k)=t_cup (i,k)
            qrc     (i,k)=0.          !--- liq/ice water
            qc      (i,k)=qe_cup(i,k) !--- total water: liq/ice = vapor water
         !qc2     (i,k)=qe_cup(i,k) !--- total water: liq/ice = vapor water
         enddo
      enddo

      !--- get boundary condition for qc
      do i=its,itf
         if(ierr(i)  /= 0) cycle
         call get_cloud_bc(name,kts,kte,ktf,xland(i),po(i,kts:kte),qe_cup (i,kts:kte),qaver,k22(i))
         qc  (i,kts:start_level(i)) = qaver + zqexec(i) +     FRACT* x_add_buoy(i)/xlv
         qrc (i,kts:start_level(i)) = 0.
      enddo

      !--- option to produce linear fluxes in the sub-cloud layer.
      if(name == 'shallow' .and. use_linear_subcl_mf == 1) then
         do i=its,itf
            if(ierr(i) /= 0) cycle
            call get_delmix(name,kts,kte,ktf,xland(i),start_level(i),po(i,kts:kte) &
               ,qe_cup(i,kts:kte), qc(i,kts:kte))
         enddo
      endif
      do i=its,itf
         if (ierr(i) /= 0) cycle

         do k=start_level(i) + 1,ktop(i) + 1

            DZ=Z_cup(i,K+1)-Z_cup(i,K)
            !
            !--- saturation  in cloud, this is what is allowed to be in it
            !
            QRCH = qes_cup(I,K)+(1./XLV)*(GAMMA_cup(i,k)/(1.+GAMMA_cup(i,k)))*DBY(I,K)

            !-    1. steady state plume equation, for what could
            !-       be in cloud without condensation
            denom =  (zu(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1)) + 1.e-12
            
            !if (denom == 0.)  STOP "DENOM == ZERO"
            !
            !if(denom > 0.) then

            qc (i,k)=  ( qc (i,k-1)*zu(i,k-1)-.5*up_massdetr(i,k-1)* qc(i,k-1) +   &
                                                    up_massentr(i,k-1)* q (i,k-1))/ denom

            if(k==start_level(i)+1) &
            qc(i,k)= qc(i,k) + (zqexec(i) + FRACT*x_add_buoy(i)/xlv) * up_massentr(i,k-1)/denom
               
            !--- assuming no liq/ice water in the environment
            qrc(i,k)=  ( qrc(i,k-1)*zu(i,k-1)-.5*up_massdetr(i,k-1)* qrc(i,k-1))/ denom

            !else
            !   qc (i,k)=    qc (i,k-1)
            !   qrc(i,k)=    qrc(i,k-1)
            !endif

            !            qc2(i,k)= ( (1.-0.5*entr_rate_2d(i,k-1)*dz)*qc2(i,k-1)     &
            !                              + entr_rate_2d(i,k-1)*dz *q  (i,k-1) ) / &
            !                        (1.+0.5*entr_rate_2d(i,k-1)*dz)

            !-- updraft temp
            tempc(i,k) = (1./cp)*(hc(i,k)-g*z_cup(i,k)-xlv*QRCH)

            !--- total condensed water before rainout
            clw_all(i,k)= max(0.,qc(i,k)-qrch)

            qrc   (i,k) = min(clw_all(i,k),qrc(i,k))

            !--- production term => condensation/diffusional growth
            cup         = max(0.,qc(i,k)-qrch-qrc(i,k))/dz

            if(c0 < 1.e-6)then
               qrc (i,k) = clw_all(i,k)
               qc  (i,k) = qrc(i,k)+min(qc(i,k),qrch)
               pwav(i)   = 0.
               psum(i)   = psum(i)+clw_all(i,k)*zu(i,k) *dz
               cycle
            endif

            if (autoconv == 1 ) then
               min_liq  = qrc_crit * ( xland(i)*1. + (1.-xland(i))*0.7 )
               if(trim(name) == 'mid') min_liq = min_liq*0.5 
               
               cx0     = (c1d(i,k)+c0)*DZ
               !cx0  = 1.e-1 * float(JL)/5. *DZ * 1.e-3
               qrc(i,k)= clw_all(i,k)/(1.+cx0)
               pw (i,k)= cx0*max(0.,qrc(i,k) - min_liq)! units kg[rain]/kg[air]

               !--- convert pw to normalized pw
               pw (i,k)=pw(i,k)*zu(i,k)

             elseif (autoconv == 2 ) then
               
               min_liq  = qrc_crit * ( xland(i)*1. + (1.-xland(i))*0.7 )
                              
               !-- Yin, J., et al, 2015. J Meteorol Res.
               !-- includes Yin's factor 
               !hei = max(2000.,min(9000.,z1(i)+0.5*(z_cup(i,k)+z_cup(i,k-1))))/9500.
               !qrc_crit_yin = -500.*log(hei)*1.e-6/rho(i,k)  ! units kg/kg
               !min_liq=(qrc_crit_yin/qrc_crit_yin_ref)*min_liq
               !---
               !
               if(trim(name).eq.'mid') min_liq = min_liq*0.5 

               !-- Eq: qrc + pw = clw_all,  if qrc>min_liq
               !--     pw       = cx0 * (qrc-min_liq) * dz
               !-- =>  qrc      = (clw_all + c0*dz*min_liq) / (1 + c0*dz)

               if(clw_all(i,k) <= min_liq) then !=> more heating at upper levels, more detrained ice
                  qrc(i,k) = clw_all(i,k)
                  pw (i,k) = 0.
               else
                 !-- autoconversion factor
                  cx0     = (c1d(i,k)+c0)*DZ
                 !cx0  = 1.e-2 * float(JL)/1. *DZ * 1.e-3
              
                  qrc(i,k) = (clw_all(i,k)+min_liq*cx0)/(1.+cx0)  ! units kg[cloud/ice]/kg[air]
                  pw (i,k) = cx0*(qrc(i,k) - min_liq)             ! units kg[rain]/kg[air]
                  !--- convert pw to normalized pw
                  pw (i,k) = pw(i,k)*zu(i,k)
          
               endif

            elseif (autoconv == 3 ) then
              !NC=300.e+6 over the land, NC=50.e+6 over the oceans
               Nc = n_cldrop * 1.e+6; Faut= 1.5
               w_upd = min(15.,max(vvel2d(i,k),2.))
               !cx0 = 1350. * (Nc * 1.e-6)**(-1.79) * (dz/w_upd)*Faut
               !pw (i,k)= cx0*(clw_all(i,k))**2.47! units kg[rain]/kg[air]
             
               cx0 = 7.98e+10 * (Nc * 1.e-6)**(-3.01) *Faut
               !-- Yin, J., et al, 2015. J Meteorol Res.
               !-- includes Yin's factor 
               !hei = max(2000.,min(9000.,z1(i)+0.5*(z_cup(i,k)+z_cup(i,k-1))))/9500.
               !qrc_crit_yin = -500.*log(hei)*1.e-6/rho(i,k)  ! units kg/kg
               !cx0 = 1.0*(qrc_crit_yin/qrc_crit_yin_ref)*cx0
               
               pw (i,k)= cx0*exp(log(clw_all(i,k))*4.22)*(dz/w_upd)! units kg[rain]/kg[air]
               
               pw (i,k)= min( pw(i,k),clw_all(i,k) )
               qrc(i,k)=clw_all(i,k)-pw(i,k)     ! units kg[cloud/ice]/kg[air]
               
               if(qrc(i,k) < 0.) stop 'qrc(i,k) < 0. AUTO 3'
               !--- convert pw to normalized pw
               pw (i,k)=pw(i,k)*zu(i,k)

            elseif (autoconv == 4 ) then
              !NC=300.e+6 over the land, NC=50.e+6 over the oceans
               Nc = n_cldrop *1.e+6 
               min_liq  = 4./3.*pi*rhowater*r0**3*Nc/rhoair0  ! 0.419e-3 -- .61e-3
              
              if(clw_all(i,k) <= min_liq) then !=> more heating at upper levels, more detrained ice
              
                  qrc(i,k)= clw_all(i,k)
                  pw(i,k) = 0.
              
              else
                  w_upd = min(15.,max(vvel2d(i,k),2.))
                 
                  cx0 = .104*9.8*peaut/(Nc*rhowater)**(1./3.)/xmyu*rhoair0**(4./3.)! 7.03

                 !pw (i,k)= cx0*(clw_all(i,k)**7./3.)* dz/w_upd ! units kg[rain]/kg[air]
                  pw (i,k)= cx0*exp(log(clw_all(i,k))*((7./3.)))* dz/w_upd ! units kg[rain]/kg[air]

                  pw (i,k)= min( pw(i,k),clw_all(i,k) )
                  qrc(i,k)=clw_all(i,k)-pw(i,k)          ! units kg[cloud/ice]/kg[air]
               
                  if(qrc(i,k) < 0.) stop 'qrc(i,k) < 0. AUTO 4'
                  !--- convert pw to normalized pw
                  pw (i,k)=pw(i,k)*zu(i,k)
               endif

            elseif (autoconv == 9 ) then
               !  c0_deep     = 1.5e-3; c0_mid     = 1.5e-3 ; qrc_crit        = 1.e-4 !(kg/kg)

               min_liq  = qrc_crit * ( xland(i)*0.4 + (1.-xland(i))*1. )

               if(clw_all(i,k) <= min_liq) then !=> more heating at upper levels, more detrained ice

                  qrc(i,k)= clw_all(i,k)
                  pw(i,k) = 0.
               else

                  cx0     = (c1d(i,k)+c0)*(1.+ 0.33*fract_liq_f(tempc(i,k)))
                  !cx0     = (c1d(i,k)+c0)*(1.+ 2.*fract_liq_f(tempc(i,k)))
                  !--- v0
                  qrc(i,k)= qrc(i,k)*exp(-cx0*dz) + (cup/cx0)*(1.-exp(-cx0*dz))
                  qrc(i,k)= max(qrc(i,k),min_liq)
                  pw (i,k)= max(0.,clw_all(i,k)-qrc(i,k)) ! units kg[rain]/kg[air]
                  qrc(i,k)= clw_all(i,k)-pw(i,k)
                  !--- v1
                  !  qrc_0   = qrc(i,k)
                  !  qrc(i,k)= (qrc_0-min_liq)*exp(-cx0*dz) + (cup/cx0)*(1.-exp(-cx0*dz))+min_liq
                  !  qrc(i,k)= max(qrc(i,k),min_liq)
                  !  pw (i,k)= max(0.,clw_all(i,k)-qrc(i,k)) ! units kg[rain]/kg[air]
                  !  qrc(i,k)= clw_all(i,k)-pw(i,k)

                  !  qrc(i,k)= (clw_all(i,k)-min_liq)*exp(-cx0*dz)+min_liq
                  !  pw (i,k)= clw_all(i,k)-qrc(i,k) ! units kg[rain]/kg[air]
                  !--- v3
                  !  qrc(i,k)= (clw_all(i,k)-min_liq) / (1.+cx0*dz)+min_liq
                  !  pw (i,k)= cx0*dz*(qrc(i,k)-min_liq) ! units kg[rain]/kg[air]
                  !  print*,"BG=",k,real(cx0*1.e+3,4),real(pw(i,k),4),real(qrc(i,k),4)&
                  !              ,real(clw_all(i,k)-pw(i,k)-qrc(i,k),4) !==> must be zero

                  !--- convert pw to normalized pw
                  pw (i,k)= pw(i,k)*zu(i,k)
               endif

            elseif (autoconv == 6 ) then
               min_liq  = 0.5* qrc_crit * (xland(i)*1.5+(1.-xland(i))*2.5)

               if(clw_all(i,k) <= min_liq) then
                  qrc(i,k)= clw_all(i,k)
                  pw(i,k) = 0.
               else
                  cx0     = (c1d(i,k)+c0)*dz
                  qrc(i,k)= (clw_all(i,k))*exp(-cx0)
                  pw (i,k)= clw_all(i,k) - qrc(i,k)
                  pw (i,k)= pw(i,k)*zu(i,k)
               endif
                !
                !print*,"6mass=",pw(i,k)/(1.e-12+zu(i,k))+qrc(i,k),clw_all(i,k)
            elseif (autoconv == 7 ) then
               min_liq  = 0.5* qrc_crit * (xland(i)*1.5+(1.-xland(i))*2.5)

               if(clw_all(i,k) <= min_liq) then
                  qrc(i,k)= clw_all(i,k)
                  pw(i,k) = 0.
               else
                  cx0     = c1d(i,k)+c0
                  qrc_0   = qrc(i,k)
                  qrc(i,k)= qrc_0*exp(-cx0*dz) + (cup/cx0)*(1.-exp(-cx0*dz))

                  pw (i,k)= max(clw_all(i,k) - qrc(i,k),0.)
                  qrc(i,k)= clw_all(i,k) - pw (i,k)
                  pw (i,k)= pw(i,k)*zu(i,k)
               endif
                !
                !print*,"6mass=",pw(i,k)/(1.e-12+zu(i,k))+qrc(i,k),clw_all(i,k)

            elseif (autoconv == 8 ) then
               min_liq  =qrc_crit ! * (xland(i)*1.5+(1.-xland(i))*2.5)

               if(clw_all(i,k) <= min_liq) then
                  qrc(i,k)= clw_all(i,k)
                  pw(i,k) = 0.
               else
                  DELT=-5.
                  if(t_cup(i,k) > 273.16 + DELT) then
                     aux = 1.
                  else
                     aux = 1. * exp(0.07* (t_cup(i,k) - (273.16 + DELT)))
                  endif
                  cx0     = aux*c0
                  !                      cx0     = max(cx0,c0)
                  !                      cx0     = max(cx0,0.25*c0)
                  cx0     = max(cx0,0.50*c0)
                  qrc_0   = qrc(i,k)
                  qrc(i,k)= qrc_0*exp(-cx0*dz) + (cup/cx0)*(1.-exp(-cx0*dz))
                  qrc(i,k)= min(clw_all(i,k), qrc(i,k))
                  pw (i,k)= clw_all(i,k) - qrc(i,k)
                  pw (i,k)= pw(i,k)*zu(i,k)
                 !if(pw(i,k)<0.) stop " pw<0 autoc 3"
               endif

            elseif (autoconv == 5 ) then

               min_liq  = (xland(i)*0.3+(1.-xland(i))*0.5)*1.e-3

               if(clw_all(i,k) > min_liq) then

                  tem1 = fract_liq_f(tempc(i,k))
                  cbf  = 1.
                  if(tempc(i,k) < T_BF) cbf=1.+0.5*sqrt(min(max(T_BF-tempc(i,k),0.),T_BF-T_ice_BF))
                  !qrc_crit_BF = 0.5e-3/cbf
                  qrc_crit_BF = 3.e-4/cbf
                  cx0 = c0*cbf*(tem1*1.3+(1.-tem1))/(0.75*min(15.,max(vvel2d(i,k),1.)))

                  !---solution 1 by Runge-Kutta
                  !step = cx0*dz
                  !do n=int(rk),1,-1
                  !  aux     = qrc(i,k)/qrc_crit_BF
                  !  pw (i,k)= auto_rk(n,step,aux,xexp,qrc(i,k))
                  !  qrc(i,k)= max(clw_all(i,k) - pw(i,k), min_liq)
                  !enddo
                  !---

                  !---solution 2 by Runge-Kutta
                  !qrc_0 = qrc(i,k)
                  !step  = cx0*dz
                  !do n = int(rk),1,-1
                  !aux      = qrc(i,k)/qrc_crit_BF
                  !pw (i,k) =-step*qrc(i,k)*(1.0-exp(-aux**xexp))/float(n) + cup*dz/float(n)
                  !pw (i,k) = max(-qrc_0, pw(i,k))
                  !qrc(i,k) = qrc_0 + pw(i,k)
                  !enddo
                  !---

                  !---analytical solution
                  qrc_0   = qrc(i,k)
                  cx0     = cx0* (1.- exp(- (qrc_0/qrc_crit_BF)**2))
                  qrc(i,k)= qrc_0*exp(-cx0*dz) + (cup/cx0)*(1.-exp(-cx0*dz))

                  pw (i,k)= max(clw_all(i,k) - qrc(i,k),0.)
                    !--- convert PW to normalized PW
                  pw (i,k) = pw(i,k)*zu(i,k)

                   !if(pw(i,k)<-1.e-12) stop " pw<0 autoc 4"
               else
                  pw (i,k) = 0.0
                  qrc(i,k) = clw_all(i,k)
               endif
            endif
            !- total water (vapor + condensed) in updraft after the rainout
            qc(i,k)=qrc(i,k)+min(qc(i,k),qrch)

            !--- integrated normalized condensates
            pwav(i)=pwav(i)+pw(i,k)
            psum(i)=psum(i)+clw_all(i,k)*zu(i,k) *dz

         enddo
         if(pwav(i) < 0.) then
               ierr(i)=66
               ierrc(i)="pwav negative"
         endif

      enddo

      !--- get back water vapor qc
      do i=its,itf
         if (ierr(i)  /= 0) cycle
         do k=start_level(i)+1,kbcon(i) + 1
            vs  = 0.
            dz1m= 0.
            do k1 = max(k-n_smooth,kts),min(k+n_smooth,ktf)
               dz   = z_cup(i,k1+1)-z_cup(i,k1)
               vs   =  vs + dz*qrc(i,k1)
               dz1m = dz1m + dz
            enddo
            qrc(i,k) = vs/(1.e-16+dz1m)
               !if(k>ktop(i)-3)print*,"v2=",k,ktop(i),sqrt(vvel2d(i,k)),sqrt(vvel2d(i,ktop(i)))
         enddo
         do k=kts,ktop(i)+1
            qc(i,k)= qc(i,k)-qrc(i,k)
           !print*,"qc=",kbcon(i),qc(i,k)*1000,qrc(i,k)*1000,1000*(qc(i,k)+qrc(i,k))
           !if(qc(i,k) < 0.)stop " qc negative"
         enddo
      enddo

   end subroutine cup_up_moisture
