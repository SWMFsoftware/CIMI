pro SWMF_GLSETUP, FILE=FILE, UseBATS=UseBATS, PlotRadius=PlotRadius,   $
                  nSmooth=nSmooth, CMEspeed=CMEspeed, USEPIL=USEPIL,   $
                  CMEGrid=CMEGrid, ARMag=ARMag, ARSize_OFF=ARSize_OFF, $
                  GLRadius=GLRadius, SizeFactor=SizeFactor,            $
                  GLRadiusRange=GLRadiusRange, Help=Help

;-----------------------------------------------------------------------
; NAME:
;   SWMF_GLSETUP
; PURPOSE:
;   Determine the Gibson-Low flux rope parameters from the input
;   magnetogram and observed CME speed.
;
; INPUT PARAMETERS:
;   The observed CME speed, the input magnetic field of SWMF, the
;   location of the CME source region(Interactive Selection).
;
; OUTPUTS:
;   Recommended GL flux rope parameters.
;
; SYSTEM REQUIREMENTS:
; Mouse with left and right button. See GLSETUP_README about
; possible issues with cursor in Mac OS.
;
;
; KEYWORDS: 
;
;   FILE = input magnetogram file (can be FITS or SWMF format).
;   UseBATS = if set, will read BATS-R-US format (2D or 3D). Default
;             will read FITS format.
;   PlotRadius = Set up the layer of the magnetogram for 3D input
;         Cannot be used 2D file, requires /UseBATS. Default is 1.0.
;   nSmooth = If nSmooth is ODD integer larger than 1, apply boxcar smoothing on the 
;             magnetic field. This can help finding the PIL. Default
;             is 5.   
;   CMESpeed = Observed CME speed in km/s.
;
;   UsePIL = If set (default), the orientation of the flux rope will be
;   calculated according to the PIL direction.
;
;   CMEGrid = If set, the grid refinement parameters for CME will be
;   calculated.
;
;   ARMag = 1 or 2, if 1, the AR magnetic parameter is calculated
;   based on the average Br around the weighted center; if 2, the
;   AR magnetic parameter is calculated based on the average total field
;   along the PIL. The corresponding empirical relationships will be
;   used accordingly. The default is 1.
;
;   =================This group of paramaters serves to find GLRadius
;   There are three ways to set GLRadius:
;   1. Explicitly: 
; 
;      /ARSize_OFF , GLRadius = GLRadius
;
;   2. In terms of Polarity Inversion Line length:
;
;      /ARSize_OFF[,SizeFactor=SizeFactor]
;
;      In this case the formula, GLRadius=PILLength/SizeFactor is
;      applied.  Default SizeFactor is 26.25.
;   3.  In terms ofthe Active Region size (do not use /ARSize_OFF): 
;
;      [GLRadiusRange=GLRadiusRange]
;
;       In this case GL is limited, using
;       GLRadiusRange = 2-elements array to specify the range for GL
;           Radius. Default is [0.2,2.0].
;   ==================
;   Help = If set, print all available keywords
;
; RESTRICTIONS:
;   - PlotRadius can only be 0.015 increament from 1.0. Please use
;     default value (1.0) for GL setup since the empirical
;     relationship is derived from that layer.
;   - Current empirical relationship is based on the GONG magnetogram.
;
; MODIFICATION HISTORY:
;   Originally coded by Meng Jin@ AOSS, University of Michigan
;   v0.1 06/02/2014 Demo Version.
;   v0.2 06/05/2014 Added option to read binary data, remove the ASCII
;    template file dependance; added DemoMode and PlotRadius keywords.
;   v0.3 06/06/2014 Added option to calculate the GL orientation
;    according to the PIL.         
;   v0.4 11/24/2014 Added option to calculate the CME grid
;   v0.5 02/28/2015 Calculate the total magnetic field along the PIL
;   v0.6 03/10/2015 Added option for determining the AR parameter 
;    information.
;   v0.7 03/16/2015 Added option to manually set the radius of the
;   flux rope
;   v0.8 03/17/2015 Added two empirical equations to determine the
;   poloroidal flux of the GL flux rope.
;   v0.9 03/20/2015 Improved the coefficients of the empirical
;   equations.
;   v1.0 07/31/2015 Fixed reading new format FDIPS output. Fixed
;   PlotRadius issue. FDIPS output resolution can be flexible now.
;   2D Input option. Run without SSW option.
;   08/04/2015 G.Toth added optional nSmooth parameter for smoothing.
;      Added optional CMESpeed parameter to set the CME speed.
;      Added optional FILE parameter to set the name of the input
;      file. Show centerpoint of CME Source region with green.
;      Initialize !MOUSE.button before starting the cursor loop.
;      Removed DELVARX calls. Failed without solarsoft and not needed.
;   01/27/2016 M.Jin revised the method to calculate the GL flux rope.
;      Added option to read and use FITS file from GONG directly.
;      Added option to use Active Region Size to determine the GL size.
;      Limited the GL size between 0.2 and 2.0 Rs.
;      Updated the empirical relationship based on the GONG
;      magnetogram with nsmooth = 5.
;      Added Warning information when the poloidal flux is negative.
;      Fixed/introduced bugs.
;    01/29/2016 Added option to change GL Radius range.
;    02/18/2016 Deleted DemoMode; Orientation Angle changed according
;    to the new definition; Seperate function for magnetogram reading;
;    determine the flux rope helicity based on hemisphere 
;    (North - Negative); improve code stability.
;    03/16/2016 Igor Sokolov - edit the list of input parameters 
;    03/18/2016 Igor Sokolov: (1) renaming xPositiveSelect=>
;    xPositive,...;(2) when the weighted centers are calculated, again
;    rename xPositiveWeighted=> xPositive...;(3) reordered the lines to
;    separate two visualization sessions (browser side in web app);
;    (4) some if loops (sometimes a couple of hundred lines apart) are
;    merged to if endif else elseif. In a couple of occurences the 
;    non-tranferrable to python where construction is simplified. The
;    algorithm had not been modified even where this is hardly
;    desired. Some comments are added (need Meng's approval)
;    03/19/2016 Igor Sokolov: sorted out calculations: each of them 
;    is now done if and where is neded. Some redundant calculations
;    are commented out. 
;   
;------------------------------------------------------------------------

;Setup the color mode and a better IDL font.
  device,decomposed=0
  !p.font=1

;Setup help option
  if not keyword_set(help) then help=0
  if help then begin
       PRINT,'Use:  '
       PRINT,"IDL> SWMF_GLSETUP, FILE='FILE', /UseBATS, PlotRadius=PlotRadius,      $ "
       PRINT,'                    nSmooth=nSmooth, CMEspeed=CMEspeed, /usepil,      $ '
       PRINT,'                    /CMEGrid, ARMag=ARMag, /ARSize_OFF,               $ '
       PRINT,'                    GLRadius=GLRadius, SizeFactor=SizeFactor,         $ ' 
       PRINT,'                    GLRadiusRange=[Rmin,Rmax], /help                    '
       PRINT,'All parameters are optional, any order is allowed'
       RETURN
  endif
;Input file name is not given
  if not keyword_set(file) then begin
     file=''
     if(not keyword_set(UseBATS))then begin
        read, prompt= $
              'Input magnetogram fits file name: ', file
     endif else begin
        read, prompt= $
              'Input magnetogram BATSRUS file name: ', file
     endelse
  endif

;set default for UseFits
  if not keyword_set(UseBATS) then UseBATS=0

;Setup default PlotRadius
 if not keyword_set(PlotRadius) then PlotRadius=1.0


;set default for nSmooth
  if not keyword_set(nSmooth) then nSmooth=5

;Read Observed CME speed.
  if not keyword_set(CMESpeed) then begin
     CMESpeed=0.0
     read,prompt='Please Input the Observed CME Speed (km/s): ',CMESpeed
  endif

;set default option for ARMag
  if not keyword_set(ARMag) then ARMag=1
  if ARMag ne 1 and ARMag ne 2 then begin
     print, 'ARMag can only be 1 or 2!'
     print, 'Use Default option: ARMag = 1'
     ARMag=1
  endif

;set default for ARSize_OFF
  if not keyword_set(ARSize_OFF) then begin
     ARSize_OFF=0
     if keyword_set(GLRadius) then begin
        print, $
        'To set GLRadius explicitly, use /ARSize_OFF' 
        RETURN
     endif
     if keyword_set(SizeFactor) then begin
        print, $
        'To set SizeFactor  explicitly, use /ARSize_OFF'
        RETURN
     endif
  endif else begin
      if keyword_set(GLRadiusRange) then begin
         print, $
        'To set GLRadiusRange, do not use use /ARSize_OFF'
        RETURN
     endif
   endelse
         
;Setup default range for GL Radius
  if not keyword_set(GLRadiusRange) then GLRadiusRange=[0.2,2.0]

;If GLRadius is not given then it is set to PIL_length/SizeFactor
  if not keyword_set(SizeFactor) then SizeFactor = 26.25


;Read the magnetogram
  mag_info=read_magnetogram(file,PlotRadius,UseBATS)
  nlat=mag_info.nlat
  nlon=mag_info.nlon
  longitude=mag_info.longitude
  latitude=mag_info.latitude
  br_field=mag_info.br_field
  btheta_field=mag_info.btheta_field
  bphi_field=mag_info.bphi_field
  
;Smooth the image, NOTE: to get equivalent result in Python, one can
;use the following function:
;
; smoothed=scipy.ndimage.filters.uniform_filter(data,size=nsmooth)
; 
  if nSmooth gt 1 then begin
     br_field     = smooth(br_field, nsmooth, /edge_truncate)
     btheta_field = smooth(btheta_field, nsmooth, /edge_truncate)
     bphi_field   = smooth(bphi_field, nsmooth, /edge_truncate)
  endif

  bt_field = sqrt(bphi_field^2+btheta_field^2+br_field^2)

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

  br_field_show=br_field
  index=where(br_field lt -20)
  br_field_show[index]=-20
  index=where(br_field gt 20)
  br_field_show[index]=20

  window,2,xs=1200,ys=1200.*float(nlat)/float(nlon)*4./3.
  loadct,0
  contour,br_field_show,min=-20,max=20,charsize=3,$
          title='SWMF Input Magnetogram (R ='$
          +strtrim(PlotRadius,2)+' Rs)',xtitle='Solar Longitude (Pixel)',$
          ytitle='Solar Latitude (Pixel)',/fill,nlevels=60,/iso,xstyle=1,ystyle=1
  
  print,'Please Select the CME Source Region (POSITIVE) with the left button'
  loadct,39
  !MOUSE.button = 0
  while(!MOUSE.button ne 1) do begin
     cursor,xPositive,yPositive,/data,/down
     if br_field[xPositive,yPositive] lt 0 then begin
        print,'Negative Polarity! Please Select POSITIVE Polarity!'   
        !MOUSE.button=0
     endif else begin
        plots,xPositive,yPositive,/data,psym=-2,color=250 
     endelse
  endwhile
  print,'Positive Source Region Selected: ',round(xPositive),round(yPositive)
  print,'Please Select the CME Source Region (NEGATIVE) with the right button'
  while(!MOUSE.button ne 4) do begin
     cursor,xNegative,yNegative,/data,/down
     if br_field[xNegative,yNegative] gt 0 then begin
        print,'Positive Polarity! Please Select NEGATIVE Polarity!'   
        !MOUSE.button=0
     endif else begin
        plots,xNegative,yNegative,/data,psym=-2,color=50
     endelse
  endwhile
  print,'Negative Source Region Selected: ',round(xNegative),round(yNegative)
; We have now initial guesses for centers of positive and negative
; spots. Now, we improve these estimates

;Set the box size. This box size is used to search the weight center around the
;selected positive/negative points in the interactive selection. It maybe increased
;for larger active regions.
  boxsize = round(long(16)*nlon/360)

;Calculate the weighted center for the positive polarity.
  xPositiveWeight=0.
  yPositiveWeight=0.
  TotalPositiveFlux=0.
  for i=xPositive-boxsize/2,xPositive+boxsize/2 do begin
     for j=yPositive-boxsize/2,yPositive+boxsize/2 do begin
        if br_field[i,j] gt 0 then begin
           xPositiveWeight=xPositiveWeight+br_field[i,j]*i
           yPositiveWeight=yPositiveWeight+br_field[i,j]*j
           TotalPositiveFlux=TotalPositiveFlux+br_field[i,j]
        endif
     endfor
  endfor
  xPositive=xPositiveWeight/TotalPositiveFlux
  yPositive=yPositiveWeight/TotalPositiveFlux

;Calculate the weighted center for the negative polarity.
  xNegativeWeight=0.
  yNegativeWeight=0.
  TotalNegativeFlux=0.
  for i=xNegative-boxsize/2,xNegative+boxsize/2 do begin
     for j=yNegative-boxsize/2,yNegative+boxsize/2 do begin
        if br_field[i,j] lt 0 then begin
           xNegativeWeight=xNegativeWeight+br_field[i,j]*i
           yNegativeWeight=yNegativeWeight+br_field[i,j]*j
           TotalNegativeFlux=TotalNegativeFlux+br_field[i,j]
        endif
     endfor
  endfor
  xNegative=xNegativeWeight/TotalNegativeFlux
  yNegative=yNegativeWeight/TotalNegativeFlux

;Distance between the spot centers
  Dis_Weight=sqrt((xNegative-xPositive)^2+(yNegative-yPositive)^2)


;Extract the profile along the two weighted centers in order to determine the 
;find center of the flux rope.
  nProfile = max([round(abs(xPositive-xNegative)), $
                 round(abs(yPositive-yNegative))]) +1
  Profile = indgen(nProfile)
  iXProfile=round(xPositive+(xNegative-xPositive)*Profile/(nProfile-1.))
  iYProfile=round(yPositive+(yNegative-yPositive)*Profile/(nProfile-1.))
  magProfile=fltarr(nProfile)

  for i = 0, nProfile-1 do begin
     magProfile[i]=br_field[iXProfile[i],iYProfile[i]]
  endfor

  temp=min(abs(magProfile),iPIL)
  ar_center=[iXProfile[iPIL],iYProfile[iPIL]]

; Distances from the AR center and spot centers:
  DisCenter=fltarr(nlon,nlat)
  P_Center=fltarr(nlon,nlat)
  N_Center=fltarr(nlon,nlat)
  for i=0,nlon-1 do begin
     for j=0,nlat-1 do begin
        DisCenter[i,j]=sqrt((i-ar_center[0])^2+(j-ar_center[1])^2)
        P_Center[i,j]=sqrt((i-xPositive)^2+(j-yPositive)^2)
        N_Center[i,j]=sqrt((i-xNegative)^2+(j-yNegative)^2)
     endfor
  endfor


;and set GLLatitude and GLLongitude
  GL_Latitude=Latitude[ar_center[0],ar_center[1]]
  GL_Longitude=Longitude[ar_center[0],ar_center[1]]
  GL_Latitude=GL_Latitude*180./!DPI
  GL_Longitude=GL_Longitude*180./!DPI

;Bug: neither ndim nor param is defined
;   if UseBATS and ndim eq 2 then GL_Longitude = GL_Longitude + param


;Calculate Active Region Size for determining the GL size.
  AR_threshold=round(long(12)*nlon/360*Dis_Weight/13.6)
  ARmap_p=fltarr(nlon,nlat)
  ARmap_n=fltarr(nlon,nlat)
  ARmap_p[where(br_field gt 20)]=1
  ARmap_n[where(br_field lt -20)]=1

  Regionmap_p=fltarr(nlon,nlat)
  Regionmap_n=fltarr(nlon,nlat)
  Regionmap_p[where(P_Center le AR_threshold)]=1
  Regionmap_n[where(N_Center le AR_threshold)]=1

  sizemap_p=br_field*ARmap_p*Regionmap_p
  sizemap_n=br_field*ARmap_n*Regionmap_n
 
;Calculate the gradient of the Br field
  ddx=(shift(br_field,-1,0)-shift(br_field,1,0))/2.
  ddy=(shift(br_field,0,-1)-shift(br_field,0,1))/2.
  ddx[0,*]=br_field[1,*]-br_field[0,*]
  ddx[nlon-1,*]=br_field[nlon-1,*]-br_field[nlon-2,*]
  ddy[*,0]=br_field[*,1]-br_field[*,0]
  ddy[*,nlat-1]=br_field[*,nlat-1]-br_field[*,nlat-2]
  br_field_gradient=sqrt(ddx^2+ddy^2)
;Setup Bitmap for magnetic field gradient
  bitmap_gradient=fltarr(nlon,nlat)
  bitmap_gradient[*,*]=1.0
  bitmap_gradient[where(br_field_gradient lt 0.5)]=0.0

;Cell size is used to divide the magnetogram to sub regions in order to determine
;the PIL. 
  cell_size=1

;Setup the threshold for selecting cells near the PIL.
  flux_threshold=2.0

;Calculate the Bitmap (1/0) for determining the PIL.
  M=nlon/cell_size
  N=nlat/cell_size
  bitmap=fltarr(nlon,nlat)
  bitmap[*,*]=0.0
  for i=0,M-2 do begin
     for j=0,N-2 do begin
        if(min(br_field[i*cell_size:(i+1)*cell_size,j*cell_size:(j+1)*cell_size]) $
                     lt -flux_threshold  and                                      $
           max(br_field[i*cell_size:(i+1)*cell_size,j*cell_size:(j+1)*cell_size]) $
                     gt flux_threshold) then begin
           ;PIL indersects this quadrant. Even among these 4 points 
           ;we do not include those in which the graient is too small: 
           bitmap[i*cell_size:(i+1)*cell_size,j*cell_size:(j+1)*cell_size]=       $
              bitmap_gradient[i*cell_size:(i+1)*cell_size,j*cell_size:(j+1)*cell_size]
        endif
     endfor
  endfor  


;Distance cut-off for determining the PIL.
; round(8.*(nLon/360.)*DisWeight/13.6) for calculating PIL_Length
; round(6.*(nLon/360.))                for calculating bt_pil
; 3 for calculating PIL  orientation)
; round(6.*(nLon/360.))                for vizualization
  Dis_threshold=round(long(8)*nlon/360*Dis_Weight/13.6)
  DisMax=round(long(6)*nlon/360)
  Dis_threshold_s=3

  dismap=fltarr(nlon,nlat)
  dismap[where(Discenter le Dis_threshold)]=1

;The final weighted map showing the PIL of the CME source region.
  wmap=bitmap*br_field*dismap

  if keyword_set(USEPIL) then begin
;Calculate the orientation of the flux rope according to PIL 
;(make it perpendicular to PIL).
     dismap[where(Discenter gt Dis_threshold_s)]=0
     wmap_s=bitmap*br_field*dismap
     PILpoints=where(wmap_s gt 0)
     MM=n_elements(PILpoints)
     PIL_x=fltarr(MM)
     PIL_y=fltarr(MM)
     for i=0,MM-1 do begin
        PIL_y[i]=floor(PILpoints[i]/nlon)
        PIL_x[i]=PILpoints[i]-(PIL_y[i]*nlon)
     endfor
     
     PIL_xx=PIL_x[sort(PIL_x)]
     PIL_yy=PIL_y[sort(PIL_x)]
     PIL_fit=ladfit(PIL_xx,PIL_yy,/double)  
        
     if PIL_fit[1] eq 0 then begin
        if br_field(ar_center[0],ar_center[1]-2) lt 0 then r1=[0,-1] else r1=[0,1]
     endif else begin
        aa_PIL=-1./PIL_fit[1]
        if abs(aa_PIL) le 1 then begin
           bb_PIL=ar_center[1]-aa_PIL*ar_center[0]
           xx=[ar_center[0]-2,ar_center[0]-3,ar_center[0]-4]
           yy=aa_PIL*xx+bb_PIL
           ave_field=(br_field(xx[0],yy[0])+br_field(xx[1],yy[1])+br_field(xx[2],yy[2]))/3.
           if ave_field lt 0 then r1=[-1.,-aa_PIL] else r1=[1.,aa_PIL]
        endif else begin
           bb_PIL=ar_center[1]-aa_PIL*ar_center[0]
           yy=[ar_center[1]-2,AR_center[1]-3,AR_Center[1]-4]
           xx=floor((yy-bb_PIL)/aa_PIL)
           ave_field=(br_field(xx[0],yy[0])+br_field(xx[1],yy[1])+br_field(xx[2],yy[2]))/3.
           if ave_field lt 0 then r1=[-signum(aa_PIL),-signum(aa_PIL)*aa_PIL] $
           else r1=[signum(aa_PIL),signum(aa_PIL)*aa_PIL]
        endelse
     endelse   
        
     if array_equal(PIL_xx,PIL_xx[0]) then begin
        if br_field(ar_center[0]-2,ar_center[1]) lt 0 then r1=[-1,0] else r1=[1,0]
     endif
  endif else begin
;Calculate the GL flux rope orientation from the two weighted points.
     r1=[xNegative-xPositive,yNegative-yPositive]
  endelse
  r1=r1/sqrt(r1[0]^2+r1[1]^2)
  r2=[1.0,0.0]
  GL_Orientation=acos(r1[0]*r2[0]+r1[1]*r2[1])*180/!DPI
  if r1[1] lt 0 then begin
     GL_Orientation=360-GL_Orientation
  endif
  
;Calculate the poloidal flux needed for the observed CME velocity.
;These relationships are based on the GONG magnetogram with nsmooth = 5
  if ARMag eq 1 then begin
     RegionSize_ARMag=round(long(4)*nlon/360)
     br_ar=mean(abs(br_field[ar_center[0]-RegionSize_ARMag/2:ar_center[0]+RegionSize_ARMag/2,$
                             ar_center[1]-RegionSize_ARMag/2:ar_center[1]+RegionSize_ARMag/2]))
     GL_poloidal=(CMESpeed*br_ar^0.43989278-3043.9307)/565.05018
  endif else begin
;Calculate the AR magnetic field strength in order to determine the
;right poloidal flux needed.
     showpoints=where(wmap ne 0)
     NN=n_elements(showpoints)
     ;pillines=intarr(NN,2)
     bt_pil=0.0
     for i=0,NN-1 do begin
        y_show=floor(showpoints[i]/nlon)
        x_show=showpoints[i]-(y_show*nlon)
        ;pillines[i,*]=[x_show,y_show]
        if sqrt((ar_center[0]-x_show)^2+(ar_center[1]-y_show)^2) lt DisMax then begin
           bt_pil=bt_pil+bt_field[x_show,y_show]
        endif
     endfor

     bt_pil=bt_pil/nn
     GL_poloidal=(CMESpeed*bt_pil^0.58148678-2814.1030)/507.60065
  endelse
  
;Print WARNING information is GL_Bstrength is negative                                               
  if GL_poloidal le 0 then begin
     print,'*********************************************'
     print,'WARNING: CALCULATION FAILED!USE WITH CAUTION!'
     print,'Either the active region is too weak or the'    
     print,'CME speed is too small!'
     print,'GL Poloidal Flux is set to 0!'
     print,'*********************************************'
     GL_poloidal = 0.0
  endif

;Relationship between the PIL length and the GL flux rope Radius.   
;This factor is now based on the 2011 March 7 CME. More tests  
;are needed in order to get a more precise value.  
  if  ARSize_OFF then begin
     if not keyword_set(GLRadius) then begin
;At this moment, the PIL length is represented by degree and does not 
;take into account the effect of different latitude. It will be improved later. 
        PIL_Length=(n_elements(where(wmap lt 0))+n_elements(where(wmap gt 0)))/2.*360./nlon
        GLRadius=PIL_Length/SizeFactor
     endif

;Use Active Region size to specify the GL flux rope Radius.
  endif else begin
     ARSize=n_elements(where(abs(sizemap_p) ne 0))+n_elements(where(abs(sizemap_n) ne 0))
;Redundant?
;    ARStrength=(total(sizemap_p)+total(abs(sizemap_n)))/float(ARSize)
     GLRadius=0.8/280.*ARSize
     if GLRadius gt GLRadiusRange[1] then GLRadius=GLRadiusRange[1]
     if GLRadius lt GLRadiusRange[0] then GLRadius=GLRadiusRange[0]
  endelse

;Relationship between the GL Poloidal flux and GL Bstrength.
;Flux rope helicity is determined by the hemisphere, northern
;hemisphere - negative helicity.
 
  GL_Bstrength=GL_poloidal/(21.457435*GLRadius^4)
  if GL_Latitude gt 0 then GL_Bstrength = -GL_Bstrength
 

;Recommended GL flux rope parameters
  Distance = 1.8
  Stretch  = 0.6

  print,'========================================'
  print,'The Recommended GL FLux Rope Parameters'
  print,'========================================'
  print,FORMAT='(A25,5X,F6.2)','Latitude: ',GL_Latitude
  print,FORMAT='(A25,5X,F6.2)','Longitude: ',GL_Longitude
  print,FORMAT='(A25,5X,F6.2)','Orientation: ',GL_Orientation
  print,FORMAT='(A25,5X,F6.2)','Radius: ', GLRadius
  print,FORMAT='(A25,5X,F6.2)','Bstrength: ',GL_Bstrength
  print,FORMAT='(A25,5X,F6.2)','Stretch (FIXED): ',Stretch
  print,FORMAT='(A25,5X,F6.2)','Distance (FIXED): ',Distance
  print,FORMAT='(A25,5X,F6.2)','Height [Rs]: ', $
        GLRadius + Distance - Stretch - 1.0
  print,FORMAT='(A25,5X,F6.2)','Angular size [deg]: ', $
        2*GLRadius/Distance/!dtor
  print,FORMAT='(A25,5X,F6.2)','Poloidal flux [1E21 Mx]: ', GL_poloidal
  print,'-----------------------------------------'

  if keyword_set(CMEGrid) then begin
;Calculate the CME grid refinement parameters based on the flux rope
;location and size.                                                

     CMEbox_Start=[1.1,GL_Longitude-40.*GLRadius,GL_Latitude-20.*GLRadius]
     CMEbox_End=[20.0,GL_Longitude+40.*GLRadius,GL_Latitude+20.*GLRadius]

     print,'=========================================='
     print,'The Recommended Grid Refinement Parameters'
     print,'=========================================='
     print,FORMAT='(A25,5X,F6.2)','R_Start: ', CMEbox_start[0]
     print,FORMAT='(A25,5X,F6.2)','R_End: ', CMEbox_end[0]
     print,FORMAT='(A25,5X,F6.2)','Longitude_Start: ', CMEbox_start[1]
     print,FORMAT='(A25,5X,F6.2)','Longitude_End: ', CMEbox_end[1]
     print,FORMAT='(A25,5X,F6.2)','Latitude_Start: ', CMEbox_start[2]
     print,FORMAT='(A25,5X,F6.2)','Latitude_End: ', CMEbox_end[2]
     print,'-----------------------------------------'
  endif
; Final magnetogram (large) 
;Plot the weighted centers on the magnetogram.
  device,decomposed=0
  loadct,0
  contour,br_field_show,min=-20,max=20,charsize=3,title='SWMF Input Magnetogram (R ='$
          +strtrim(PlotRadius,2)+' Rs)',xtitle='Solar Longitude (Pixel)',$
          ytitle='Solar Latitude (Pixel)',/fill,nlevels=60,/iso,xstyle=1,ystyle=1

  loadct,39
  plots,xPositive,yPositive,/data,psym=-2,color=250
  plots,xNegative,yNegative,/data,psym=-2,color=50
;plot center of the flux rope 
  plots,ar_center[0],ar_center[1],/data,psym=-2,color=150
; Showing positive and negative spots
  contour,abs(sizemap_p),/overplot,c_color=100
  contour,abs(sizemap_n),/overplot,c_color=100  
;Showing the PIL
  showpoints=where(wmap gt 0)
  NN=n_elements(showpoints)
  for i=0,NN-1 do begin
     y_show=floor(showpoints[i]/nlon)
     x_show=showpoints[i]-(y_show*nlon)
     if sqrt((ar_center[0]-x_show)^2+(ar_center[1]-y_show)^2) lt DisMax then begin
        plots,x_show,y_show,psym=-1,color=200
     endif
  endfor

;The region size is used to cover the whole area of active region in
;order to show a zoom-in image. Shorter RegionSize for near-Limb
;regions when needed.

  RegionSize=round(long(50)*nlon/360)

;Display the zoom-in image of the active region with weighted centers and PIL.
  window,3,xs=800,ys=800

  device,decomposed=0
  loadct,0
  
  sub_x1=max([ar_center[0]-RegionSize/2,0])
  sub_x2=min([ar_center[0]+RegionSize/2,nlon-1])
  sub_y1=max([ar_center[1]-RegionSize/2,0])
  sub_y2=min([ar_center[1]+RegionSize/2,nlat-1])

  contour,br_field_show[sub_x1:sub_x2,sub_y1:sub_y2],$                        
          min=-20,max=20,charsize=3,title='CME Source Region (R ='$
          +strtrim(PlotRadius,2)+' Rs)',xtitle='Solar Longitude (Degree)',$
          ytitle='Solar Latitude (Pixel)',/fill,nlevels=60,/iso,xstyle=1,ystyle=1

  loadct,39
  plots,xPositive-sub_x1,yPositive-sub_y1,$
        /data,psym=-2,color=250,symsize=3,thick=3
  plots,xNegative-sub_x1,yNegative-sub_y1,$
        /data,psym=-2,color=50,symsize=3,thick=3
  plots,ar_center[0]-sub_x1,ar_center[1]-sub_y1,/data,psym=-2,color=150,symsize=3,thick=3
  for i=0,NN-1 do begin
     y_show=floor(showpoints[i]/nlon)
     x_show=showpoints[i]-(y_show*nlon)
     y_show=y_show-sub_y1
     x_show=x_show-sub_x1
     xcenter=ar_center[0]-sub_x1
     ycenter=ar_center[1]-sub_y1
     if sqrt((xcenter-x_show)^2+(ycenter-y_show)^2) lt DisMax then begin
        plots,x_show,y_show,psym=-1,color=200,symsize=3,thick=3
     endif
  endfor

;------------------------------
;Save data for testing purpose. 
;WARNING!!!Uncomment pillines which is currently commented out!!!
;------------------------------
  ;save,pillines,ar_center,filename='AR6.sav'
  ;save,pillines,ar_center,xPositive,yPositive,xNegative,yNegative,$
  ;     sizemap_p,sizemap_n,filename='AR5_show.sav'

  !mouse.button=0
end
