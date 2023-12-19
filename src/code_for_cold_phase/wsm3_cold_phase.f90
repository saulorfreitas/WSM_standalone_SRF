!
!===============================================================
!
! cold rain processes
!
! - follows the revised ice microphysics processes in HDC
! - the processes same as in RH83 and LFO behave
!   following ice crystal hapits defined in HDC, inclduing
!   intercept parameter for snow (n0s), ice crystal number
!   concentration (ni), ice nuclei number concentration
!   (n0i), ice diameter (d)
!
!===============================================================
!
            supcol = t0c-t(i,k)
            n0sfac(i,k) = max(min(exp(alpha*supcol),n0smax/n0s),1.)
            ifsat = 0
!-------------------------------------------------------------
! Ni: ice crystal number concentraiton   [HDC 5c]
!-------------------------------------------------------------
            xni(i,k) = min(max(5.38e7                                          &
                    *exp(log((den(i,k)*max(qci(i,k),qmin)))*(0.75)),1.e3),1.e6)
            eacrs = exp(0.07*(-supcol))
            if(qrs(i,k).gt.qcrmin.and.qci(i,k).gt.qmin) then
              xmi = den(i,k)*qci(i,k)/xni(i,k)
              diameter  = min(dicon * sqrt(xmi),dimax)
              vt2i = 1.49e4*diameter**1.31
              vt2s = pvts*rslopeb(i,k)*denfac(i,k)
!-------------------------------------------------------------
! praci: Accretion of cloud ice by rain [HL A15] [LFO 25]
!        (T<T0: I->R)
!-------------------------------------------------------------
              acrfac = 2.*rslope3(i,k)+2.*diameter*rslope2(i,k)                &
                      +diameter**2*rslope(i,k)
              pacr(i,k) = min(pi*qci(i,k)*eacrs*n0s*n0sfac(i,k)                &
                       *abs(vt2s-vt2i)*acrfac/4.,qci(i,k)/dtcld_1d(i,k) )
            endif
!-------------------------------------------------------------
! pidep: Deposition/Sublimation rate of ice [HDC 9]
!       (T<T0: V->I or I->V)
!-------------------------------------------------------------
            if(qci(i,k).gt.0.) then
              xmi = den(i,k)*qci(i,k)/xni(i,k)
              diameter = dicon * sqrt(xmi)
              pisd(i,k) = 4.*diameter*xni(i,k)*(rh(i,k)-1.)/work1(i,k)
              if(pisd(i,k).lt.0.) then
                pisd(i,k) = max(pisd(i,k),satdt/2)
                pisd(i,k) = max(pisd(i,k),-qci(i,k)/dtcld_1d(i,k) )
              else
                pisd(i,k) = min(pisd(i,k),satdt/2)
              endif
              if(abs(pisd(i,k)).ge.abs(satdt)) ifsat = 1
            endif
!-------------------------------------------------------------
! psdep: deposition/sublimation rate of snow [HDC 14]
!        (V->S or S->V)
!-------------------------------------------------------------
            if(qrs(i,k).gt.0..and.ifsat.ne.1) then
              coeres = rslope2(i,k)*sqrt(rslope(i,k)*rslopeb(i,k))
              pres(i,k) = (rh(i,k)-1.)*n0sfac(i,k)*(precs1*rslope2(i,k)        &
                        +precs2*work2(i,k)*coeres)/work1(i,k)
              supice = satdt-pisd(i,k)
              if(pres(i,k).lt.0.) then
                pres(i,k) = max(pres(i,k),-qrs(i,k)/dtcld_1d(i,k) )
                pres(i,k) = max(max(pres(i,k),satdt/2),supice)
              else
                pres(i,k) = min(min(pres(i,k),satdt/2),supice)
              endif
              if(abs(pisd(i,k)+pres(i,k)).ge.abs(satdt)) ifsat = 1
            endif
!-------------------------------------------------------------
! pigen: generation(nucleation) of ice from vapor [HDC 7-8]
!       (T<T0: V->I)
!-------------------------------------------------------------
            if(supsat.gt.0.and.ifsat.ne.1) then
              supice = satdt-pisd(i,k)-pres(i,k)
              xni0 = 1.e3*exp(0.1*supcol)
              roqi0 = 4.92e-11*exp(log(xni0)*(1.33))
              pgen(i,k) = max(0.,(roqi0/den(i,k)-max(qci(i,k),0.))/dtcld_1d(i,k) )
              pgen(i,k) = min(min(pgen(i,k),satdt),supice)
            endif
!-------------------------------------------------------------
! psaut: conversion(aggregation) of ice to snow [HDC 12]
!       (T<T0: I->S)
!-------------------------------------------------------------
            if(qci(i,k).gt.0.) then
              qimax = roqimax/den(i,k)
              paut(i,k) = max(0.,(qci(i,k)-qimax)/dtcld_1d(i,k) )
            endif
!