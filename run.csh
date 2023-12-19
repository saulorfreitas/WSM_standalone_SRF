#!/bin/csh

#---------------------------create the executable
rm -f wsm.x
mk


#---------------------------create wsm.inp namelist
cat << Eof1 > wsm.inp

 &run
  prefix   ="check2_cond_evap_before_SedimAfter_check11_ref0",
  mynum = 243,
  mcphys_type=30,
  time=57600,
  nsteps=100,
  sedim = .false. 
  sedim_after =.true.
  cold_phase= .true. 
  warm_phase= .true. 
  cond_evap= .false.
  cond_evap_before=.true.
 &end
Eof1

#--------------------------- run the job
wsm.x >wsm.out




