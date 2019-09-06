;---------------------------------------------------------------
;                constq
;
;  This processor plots constant Q1 and Q2 curves in (a0,E) space
;  
;  Input: constQ.dat
;         constQ2.dat
;
;  Creaded by S.-B. Kang, Code 673, GSFC
;---------------------------------------------------------------
pro constq

iwc=0
ipa=0
openr,1,'constQ.dat'
readf,1,ipc,iwc,ipa
ckeV=fltarr(iwc)
a=fltarr(ipc,iwc,ipa)
a1=fltarr(iwc)
readf,1,ckeV
for j=0,ipc-1 do begin
   for m=0,ipa-1 do begin
       readf,1,a1
       a(j,*,m)=a1
   endfor
endfor
close,1
openr,1,'constQ2.dat'
readf,1,ipc,iq,ipa
Ec=fltarr(iq)
a0=findgen(ipa)+1.
E1=fltarr(ipa)
E=fltarr(ipc,iq,ipa)
readf,1,Ec
for j=0,ipc-1 do begin
    for k=0,iq-1 do begin
        readf,1,E1
        E(j,k,*)=E1
    endfor
endfor
close,1

set_plot,'ps'
device,filename='diag_constQ.ps',/inches,/color,yoffset=0., $
       xoffset=0.1,xsize=6.,ysize=6.
for j=0,ipc-1 do begin
    ymin=min(a(j,*,*))
    ymax=max(a(j,*,*))
    print,ymin,ymax
    plot,a(j,*,0),Ec,/xs,/ys,$
	xrange=[0.1,90.],yrange=[ymin,ymax],$
	/ylog,/noerase,$
	xtitle='Pitch-angle (deg)',$
	ytitle='Energy (keV)',$
	title='Constant Q!D1!N (deg) curve'
    for m=1,ipa-1 do oplot,a(j,*,m),Ec
    erase
endfor
;;plot,a(*,0),ckeVm,/xs,/ys,xrange=[1.e1,90.],yrange=[1.e0,1.e3],/ylog
;;for m=1,ipa-1 do oplot,a(*,m),ckeVm
device,/close_file

pos1=[0.1,0.445,0.688,0.9]
device,filename='diag_constQ2.ps',/inches,/color,yoffset=0., $
       xoffset=0.1,xsize=8.5,ysize=11.
for j=0,ipc-1 do begin
    ymin=min(E)
    ymax=max(E)
    print,ymin,ymax
    ymin=0.1
    plot,a0,E(j,0,*),/xs,/ys,$
	xrange=[0.1,90.],yrange=[ymin,ymax],$
	/ylog,/noerase,pos=pos1,$
	xtitle='Pitch-angle (deg)',$
	ytitle='Energy (keV)',$
	title='Constant Q!D2!N (keV) curve'
    for k=1,iq-1 do oplot,a0,E(j,k,*)
    erase
endfor
;;plot,a0,E(0,*),/xs,/ys,xrange=[0.1,90.],yrange=[1.e-1,1.e4],/ylog
;;for k=1,iwc-1 do oplot, a0m,E(k,*)
device,/close_file

end
