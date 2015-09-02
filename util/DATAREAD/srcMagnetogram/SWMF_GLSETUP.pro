pro SWMF_GLSETUP, DemoMode=DemoMode, PlotRadius=PlotRadius, USEPIL=USEPIL
;--------------------------------------------------------------------------------------------------
; NAME :
;   SWMF_GLSETUP
; PURPOSE :
;   Determine the Gibson-Low flux rope parameters from the input magnetogram and observed CME speed.
; INPUT PARAMETERS :
;   The observed CME speed, the input magnetic field of SWMF, the location of the CME source region
;   (Interactive Selection).
; OUTPUTS :
;   Recommended GL flux rope parameters.
; KEYWORDS:
;   DemoMode = If set, the pre-saved magnetogram data at Rs=1.0 will be loaded.
;   PlotRadius = Set up the layer of the magnetogram. Cannot be used with DemoMode.
;   UsePIL = If set, the orientation of the flux rope will be caluclated according to the PIL direcion.
; CALLS   :
;   PLOT_IMAGE
; RESTRICTIONS:
;   PlotRadius can only be 0.015 increament from 1.0. Please do not use large raduis for the GL setup.   
; MODIFICATION HISTORY:
;   Originally coded by Meng Jin@ AOSS, University of Michigan
;   v0.1 06/02/2014 Demo Version.
;   v0.2 06/05/2014 Added option to read binary data, remove the ASCII template file dependance.
;                   Added DemoMode and PlotRadius keywords.
;   v0.3 06/06/2014 Added option to calculate the GL orientation according to the PIL.                                     
;---------------------------------------------------------------------------------------------------

;Setup the color mode and a better IDL font.
device,decomposed=1
!p.font=1

;Turn on/off demo mode. With the demo mode on, the pre-saved 2D data will be read instead
;of reading from 3D data which is much more time consuming
if not keyword_set(DemoMode) then DemoMode=0


;Read Observed CME speed.
file=''
CMESpeed=0.0
read,prompt='Please Input the Observed CME Speed (km/s): ',CMESpeed

;Setup the magnetogram layer, default is at the solar surface
if not keyword_set(PlotRadius) then  PlotRadius=1.00

if keyword_set(DemoMode) and keyword_set(PlotRadius) then begin
  print,'DemoMode and PlotRadius Cannot be Used at the same time!'
  print,'Play DemoMode...'
  PlotRadius=1.00
endif

;Read the SWMF input magnetic field
if not DemoMode then begin
  read, prompt='Input Magnetic Field of SWMF (Format can be ASCII/Binary): ',file
  data={field1:dblarr(6562980),field2:dblarr(6562980),field3:dblarr(6562980),field4:dblarr(6562980)}
  if QUERY_ASCII(file,info) then begin
    openr,lun,file,/get_lun
    line=''
    linedata=dblarr(6)    
    for i=0,4 do begin
      readf,lun,line
    endfor
    
    for i=0,6562980-1 do begin
       readf,lun,linedata
       data.field1[i]=linedata[0]
       data.field2[i]=linedata[1]
       data.field3[i]=linedata[2]
       data.field4[i]=linedata[3]
    endfor

  endif else begin
    openr, lun, file, /get_lun, /F77_UNFORMATTED
    line1=bytarr(48)
    line2=lonarr(6)
    line3=lonarr(3)
    line4=dblarr(4)
    line5=bytarr(80)
    tmp_line1=dblarr(6562980*3)  
    tmp_line2=dblarr(6562980)
    readu,lun,line1
    readu,lun,line2
    readu,lun,line3
    readu,lun,line4
    readu,lun,line5
    readu,lun,tmp_line1
    data.field1=tmp_line1[0:6562980-1]
    data.field2=tmp_line1[6562980:6562980*2-1]
    data.field3=tmp_line1[6562980*2:6562980*3-1]
    readu,lun,tmp_line2
    data.field4=tmp_line2
    free_lun,lun
  endelse
  
  ;Read out data into 2D arrays
  nn=long(0)
  index=where(abs(data.field1-PlotRadius) lt 0.001)
  Br_field=dblarr(360,180)
  ;Bphi_field=dblarr(360,180)
  ;Btheta_field=dblarr(360,180)
  Longitude=dblarr(360,180)
  Latitude=dblarr(360,180)
  for i=0,179 do begin
    for j=0,360 do begin
      if j eq 360 then begin
        nn=nn+1
        break
      endif else begin
        Longitude[j,i]=data.field2[index[nn]]
        Latitude[j,i]=data.field3[index[nn]]
        Br_field[j,i]=data.field4[index[nn]]
        ;Bphi_field[j,i]=data.field5[index[nn]]
        ;Btheta_field[j,i]=data.field6[index[nn]]
        nn=nn+1
      endelse
    endfor
  endfor  
endif

;Restore pre-saved 2D data for the demo mode
if DemoMode then begin
  restore,'magneticfield_1.0.sav'
endif

;Display the magnetogram and let user interactively select the CME source region. The
;procedure to select is:
; 1. Click the CME source region of positive polarity with 'left' button of mouse
; 2. Click the CME source region of negative polarity with 'right' button of mouse
;
;Note that the user can click anywhere inside the active region. However, click closer
;to the center of the positive/negative patterns is recommended.
;
;Note the solar latitude is expressed in pixel due to the non-uniform spacing. The latitude
;is uniform in sin(latitude). This will be changed in the future to degree. 
 
window,2,xs=1200,ys=800
plot_image,br_field,min=-40,max=40,charsize=3,title='SWMF Input Magnetogram (R ='$
+strtrim(PlotRadius,2)+' Rs)',xtitle='Solar Longitude (Degree)',$
ytitle='Solar Latitude (Pixel)'
print,'Please Select the CME Source Region (POSITIVE)'
while(!MOUSE.button ne 1) do begin
  cursor,xPositiveSelect,yPositiveSelect,/data,/down
  if br_field[xPositiveSelect,yPositiveSelect] lt 0 then begin
    delvarx,xPositiveSelect,yPositiveSelect
    print,'Negative Polarity! Please Select POSITIVE Polarity!'   
    !MOUSE.button=0
  endif else begin
    plots,xPositiveSelect,yPositiveSelect,/data,psym=-2,color='0000FF'XL 
  endelse
endwhile
print,'Positive Source Region Selected: ',round(xPositiveSelect),round(yPositiveSelect)
print,'Please Select the CME Source Region (NEGATIVE)'
while(!MOUSE.button ne 4) do begin
  cursor,xNegativeSelect,yNegativeSelect,/data,/down
  if br_field[xNegativeSelect,yNegativeSelect] gt 0 then begin
    delvarx,xNegativeSelect,yNegativeSelect
    print,'Positive Polarity! Please Select NEGATIVE Polarity!'   
    !MOUSE.button=0
  endif else begin
    plots,xNegativeSelect,yNegativeSelect,/data,psym=-2,color='0000FF'XL
  endelse
endwhile
print,'Negative Source Region Selected: ',round(xNegativeSelect),round(yNegativeSelect)

;Set the box size. This box size is used to search the weight center around the
;selected positive/negative points in the interactive selection. It maybe increased
;for larger active regions.
boxsize=16

;Calculate the weighted center for the positive polarity.
xPositiveWeight=0.
yPositiveWeight=0.
TotalPositiveFlux=0.
for i=xPositiveSelect-boxsize/2,xPositiveSelect+boxsize/2 do begin
  for j=yPositiveSelect-boxsize/2,yPositiveSelect+boxsize/2 do begin
    if br_field[i,j] gt 0 then begin
      xPositiveWeight=xPositiveWeight+br_field[i,j]*i
      yPositiveWeight=yPositiveWeight+br_field[i,j]*j
      TotalPositiveFlux=TotalPositiveFlux+br_field[i,j]
    endif
  endfor
endfor
xPositiveWeight=xPositiveWeight/TotalPositiveFlux
yPositiveWeight=yPositiveWeight/TotalPositiveFlux

;Calculate the weighted center for the negative polarity.
xNegativeWeight=0.
yNegativeWeight=0.
TotalNegativeFlux=0.
for i=xNegativeSelect-boxsize/2,xNegativeSelect+boxsize/2 do begin
  for j=yNegativeSelect-boxsize/2,yNegativeSelect+boxsize/2 do begin
    if br_field[i,j] lt 0 then begin
      xNegativeWeight=xNegativeWeight+br_field[i,j]*i
      yNegativeWeight=yNegativeWeight+br_field[i,j]*j
      TotalNegativeFlux=TotalNegativeFlux+br_field[i,j]
    endif
  endfor
endfor
xNegativeWeight=xNegativeWeight/TotalNegativeFlux
yNegativeWeight=yNegativeWeight/TotalNegativeFlux

;Plot the weighted centers on the magnetogram.
plot_image,br_field,min=-40,max=40,charsize=3,title='SWMF Input Magnetogram (R = '+strtrim(PlotRadius,2)+' Rs)',$
xtitle='Solar Longitude (Degree)',ytitle='Solar Latitude (Pixel)'
plots,xPositiveWeight,yPositiveWeight,/data,psym=-2,color='0000FF'XL
plots,xNegativeWeight,yNegativeWeight,/data,psym=-2,color='FF0000'XL

;Calculate the GL flux rope orientation from the two weighted points.
;Calculate the orientation from the PIL (In development)
r1=[xPositiveWeight-xNegativeWeight,yPositiveWeight-yNegativeWeight]
r1=r1/sqrt(r1[0]^2+r1[1]^2)
r2=[1.0,0.0]
GL_Orientation=acos(r1[0]*r2[0]+r1[1]*r2[1])*180/3.1415926
if r1[1] lt 0 then begin
  GL_Orientation=360-GL_Orientation
endif

;Extract the profile along the two weighted centers in order to determine the 
;center of the flux rope.
aa=(yPositiveWeight-yNegativeWeight)/(xPositiveWeight-xNegativeWeight)
bb=yPositiveWeight-aa*xPositiveWeight
if abs(xPositiveWeight-xNegativeWeight) gt abs(yPositiveWeight-yNegativeWeight) then begin
  xProfile=min([xPositiveWeight,xNegativeWeight])+indgen(round(abs(xPositiveWeight-xNegativeWeight))+1)
  yProfile=round(aa*xProfile+bb)
  xProfile=round(xProfile)
endif else begin
  yProfile=min([yPositiveWeight,yNegativeWeight])+indgen(round(abs(yPositiveWeight-yNegativeWeight))+1)
  xProfile=round((yProfile-bb)/aa)
  yProfile=round(yProfile)
endelse

nProfile=n_elements(xProfile)
magProfile=fltarr(nProfile)

for i=0,nProfile-1 do begin
  magProfile[i]=br_field[xProfile[i],yProfile[i]]
endfor

temp=min(abs(magProfile),index)
plots,xProfile[index],yProfile[index],/data,psym=-2,color='00FF00'XL
GL_Latitude=Latitude[xProfile[index],yProfile[index]]
GL_Longitude=Longitude[xProfile[index],yProfile[index]]
GL_Latitude=GL_Latitude*180./3.1415926
GL_Longitude=GL_Longitude*180./3.1415926

;Calculate the gradient of the Br field
ddx=(shift(br_field,-1,0)-shift(br_field,1,0))/2.
ddy=(shift(br_field,0,-1)-shift(br_field,0,1))/2.
ddx[0,*]=br_field[1,*]-br_field[0,*]
ddx[359,*]=br_field[359,*]-br_field[358,*]
ddy[*,0]=br_field[*,1]-br_field[*,0]
ddy[*,179]=br_field[*,179]-br_field[*,178]
br_field_gradient=sqrt(ddx^2+ddy^2)

;Cell size is used to divide the magnetogram to sub regions in order to determine
;the PIL. 
cell_size=2

;Setup the threshold for selecting cells near the PIL.
flux_threshold=1.0

;Calculate the Bitmap (1/0) for determining the PIL.
M=360/cell_size
N=180/cell_size
bitmap=fltarr(360,180)
bitmap[*,*]=0.0
for i=0,M-2 do begin
  for j=0,N-2 do begin
    index1=where(br_field[i*cell_size:(i+1)*cell_size,j*cell_size:(j+1)*cell_size] lt -flux_threshold)
    index2=where(br_field[i*cell_size:(i+1)*cell_size,j*cell_size:(j+1)*cell_size] gt flux_threshold)
    if index1[0] ne -1 and index2[0] ne -1 then begin
      bitmap[i*cell_size:(i+1)*cell_size,j*cell_size:(j+1)*cell_size]=1.0
    endif
  endfor
endfor  

;Setup Bitmap for magnetic gradient
bitmap_gradient=fltarr(360,180)
bitmap_gradient[*,*]=1.0
bitmap_gradient[where(br_field_gradient lt 0.5)]=0.0

;Distance cut-off for determining the PIL. 
Dis_threshold=6

DisCenter=fltarr(360,180)
for i=0,359 do begin
  for j=0,179 do begin
    DisCenter[i,j]=sqrt((i-xProfile[index])^2+(j-yProfile[index])^2)
  endfor
endfor

dismap=DisCenter
dismap[where(Discenter gt Dis_threshold)]=0
dismap[where(Discenter le Dis_threshold)]=1

;The final weighted map showing the PIL of the CME source region.
wmap=bitmap*br_field*bitmap_gradient*dismap

;At this moment, the PIL length is represented by degree and does not 
;take into account the effect of different latitude. It will be improved later. 
PIL_Length=(n_elements(where(wmap lt 0))+n_elements(where(wmap gt 0)))/2.

;Showing the PIL
showpoints=where(wmap gt 0)
NN=n_elements(showpoints)
for i=0,NN-1 do begin
  y_show=floor(showpoints[i]/360)
  x_show=showpoints[i]-(y_show*360)
  plots,x_show,y_show,psym=-1,color='00FFFF'XL
endfor

;Calculate the orientation of the flux rope according to PIL 
;(make it vertial to PIL).
if keyword_set(USEPIL) then begin
  Dis_threshold_s=4
  dismap[where(Discenter gt Dis_threshold_s)]=0
  wmap_s=bitmap*br_field*bitmap_gradient*dismap
  PILpoints=where(wmap_s gt 0)
  MM=n_elements(PILpoints)
  PIL_x=fltarr(MM)
  PIL_y=fltarr(MM)
  for i=0,MM-1 do begin
    PIL_y[i]=floor(PILpoints[i]/360)
    PIL_x[i]=PILpoints[i]-(PIL_y[i]*360)
  endfor
  PIL_fit=linfit(PIL_x,PIL_y,/double)
  aa_PIL=-1./PIL_fit[1]
  r3=[1.,aa_PIL]
  r3=r3/sqrt(r3[0]^2+r3[1]^2)
  GL_Orientation_s=acos(r3[0]*r2[0]+r3[1]*r2[1])*180/3.1415926
  if GL_Orientation gt 180 then begin
    GL_Orientation_s=360-GL_Orientation_s
  endif
  GL_Orientation=GL_Orientation_s
endif

;Relationship between the observed CME speed and GL Bstrength.
;This factor is now based on the 2011 March 7 CME. More tests
;are needed in order to get a more precise value.
factor_BV=2200./2.25
GL_Bstrength=CMESpeed/factor_BV

;Relationship between the PIL length and the GL flux rope Radius.
;This factor is now based on the 2011 March 7 CME. More tests
;are needed in order to get a more precise value.
factor_RL=35.
GL_Radius=PIL_Length/factor_RL

;Recommended GL flux rope parameters
print,'========================================'
print,'The Recommended GL FLux Rope Parameters'
print,'========================================'
print,FORMAT='(A20,5X,F6.2)','Latitude: ',GL_Latitude
print,FORMAT='(A20,5X,F6.2)','Longitude: ',GL_Longitude
print,FORMAT='(A20,5X,F6.2)','Orientation: ',GL_Orientation
print,FORMAT='(A20,5X,F6.2)','Radius: ', GL_Radius
print,FORMAT='(A20,5X,F6.2)','Bstrength: ',GL_Bstrength
print,FORMAT='(A20,5X,F6.2)','Stretch (FIXED): ',0.6
print,FORMAT='(A20,5X,F6.2)','Distance (FIXED): ',1.8
print,'-----------------------------------------'

;The region size is used to cover the whole area of active region in order to show a zoom-in image.
RegionSize=50

;Display the zoom-in image of the active region with weighted centers and PIL.
window,3,xs=800,ys=800
plot_image,br_field[xProfile[index]-RegionSize/2:xProfile[index]+RegionSize/2,$
yProfile[index]-RegionSize/2:yProfile[index]+RegionSize/2],min=-10,max=10,charsize=3,$
title='CME Source Region (R = '+strtrim(PlotRadius,2)+' Rs)',$
xtitle='Solar Longitude (Pixel)',ytitle='Solar Latitude (Pixel)'
plots,xPositiveWeight-(xProfile[index]-RegionSize/2),yPositiveWeight-(yProfile[index]-RegionSize/2),/data,psym=-2,$
color='0000FF'XL,symsize=3,thick=3
plots,xNegativeWeight-(xProfile[index]-RegionSize/2),yNegativeWeight-(yProfile[index]-RegionSize/2),/data,psym=-2,$
color='FF0000'XL,symsize=3,thick=3
for i=0,NN-1 do begin
  y_show=floor(showpoints[i]/360)
  x_show=showpoints[i]-(y_show*360)
  y_show=y_show-(yProfile[index]-RegionSize/2)
  x_show=x_show-(xProfile[index]-RegionSize/2)
  plots,x_show,y_show,psym=-1,color='00FFFF'XL,symsize=3,thick=3
endfor

;Display the PIL map
window,4,xs=1200,ys=800
plot_image,bitmap*br_field*bitmap_gradient,min=-20,max=20,title='PIL Map (R = '+strtrim(PlotRadius,2)+' Rs)',$
xtitle='Solar Longitude (Degree)',ytitle='Solar Latitude (Pixel)',charsize=3


!mouse.button=0
end