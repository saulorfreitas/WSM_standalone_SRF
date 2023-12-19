'reinit'
'open SedimAfter_check7_ref0_WSM30_dataOut-57600-00243-00100.ctl'
'open cond_evap_before_SedimAfter_check7_ref0_WSM30_dataOut-57600-00243-00100.ctl'

'set lat -4.65'
*'set lon -55.55'
'set z 1 41'

is1=1 ; is2=2
px=7.8
pxx=8.4
dummy = colors2()
'set display white'
'clear'
'set xlopts 1 7 0.13'
'set ylopts 1 7 0.13'
'set ylint 1000'
"set mproj off"
"set lwid 13 3";"set cthick 13"
"set axlim -0.5 0.5" 
'set vpage 0.5 11 0.5 8'
'set parea 1. 8 1 7'

ipl=1
while(ipl<=4)

   'clear'
   'set  clopts -1 -1 0.12'
   'set grid on'
   'set annot 1 10'
   'set grads off'
   "set cmark 0" 

   if(ipl=1);  var='temp.'is2'-temp.'is1'' ; name=var; icx=1 ;endif
   if(ipl=2);  var='rv.'is2'-rv.'is1'' ; name=var; icx=23 ;endif
   if(ipl=3);  var='rcp.'is2'-rcp.'is1'' ; name=var; icx=22 ;endif
   if(ipl=4);  var='rrp.'is2'-rrp.'is1'' ; name=var; icx=4 ;endif


  
   'set ccolor 'icx
   'd 'var''

*******
px=px-0.3
'set strsiz 0.12 0.18'
'set string 'icx' tl 10 0'
'draw string  'pxx' 'px' 'name''
*******

'q pos'
ipl=ipl+1
endwhile




return
******************************************
function colors2
*karla
*light yellow to dark red
'set rgb 21 255 250 170'
'set rgb 24 255 160   0'
'set rgb 25 255  96   0'
'set rgb 26 255  50   0'
'set rgb 27 225  20   0'
'set rgb 28 192   0   0'
'set rgb 29 165   0   0'
*
*light green to dark green
'set rgb 31 230 255 225'
'set rgb 32 200 255 190'
'set rgb 33 180 250 170'
'set rgb 34 150 245 140'
'set rgb 35 120 245 115'
'set rgb 36  80 240  80'
'set rgb 37  55 210  60'
'set rgb 38  30 180  30'
'set rgb 39  15 160  15'
*set rgb 39   5 150   5
*
*light blue to dark blue
'set rgb 41 225 255 255'
'set rgb 42 180 240 250'
'set rgb 43 150 210 250'
'set rgb 44 120 185 250'
'set rgb 45  80 165 245'
'set rgb 46  60 150 245'
'set rgb 47  40 130 240'
'set rgb 48  30 110 235'
'set rgb 49  20 100 210'

'set rgb  16  0    0  255 80 '
'set rgb  17  55   55  255 '
'set rgb  18  0    160  0   '
'set rgb  19  165  165  255 '
'set rgb  20  244  223  66 '
'set rgb  21  0    100  255 '
'set rgb  22  200  0 0 '
'set rgb  23  34   180 34 '
'set rgb  24  255  0 0 '
'set rgb  25  0    200  255 '
'set rgb  26  34   255 34 '
'set rgb  27  66   244 212 '
'set rgb  28  101 66 244 '
'set rgb  29  135 132 59 '
'set rgb  30 119 115 115'
*brown
'set rgb  50  142 72 72'
*brown2
'set rgb  51  200 168 137'
*green clear
'set rgb  52  102 255 51'
*blue clear
'set rgb  54  102 179 255'
'set rgb 56 204 204 51'
'set rgb 58 200 137 200'

return
