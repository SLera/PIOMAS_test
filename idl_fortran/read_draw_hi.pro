; This procedure reads lon. and lat. for an arctic ice/ocean model in generalized 
;curvilinear coordinates and the model's output of 2-D scalar fields such as sea
; ice thickness, and plots arctic sea ice thickness together with a map of the 
;Arctic. Other 2-D scalar fields can be plotted in the same way.

; written by Jinlun Zhang, PSC/APL/UW, zhang@apl.washington.edu, 3/11/2005

pro hi, ps=ps
ncol = 3 & nrow = 4
!p.multi = [0, ncol, nrow]

; Set up PostScript file or plotting window.
if keyword_set(ps) then begin
   set_plot, 'PS'
   device, filename='foo.ps', /inches, xsize=7, ysize=9.0, $
      xoffset=0.75, yoffset=1.5, $
      /portrait, /color, bits_per_pixel=8
endif else begin
   window, /free, xsize=500, ysize=500
   win = !d.window
endelse

;params for contour plot
level_dif=0.5
lowest=0.0
highest=5.0
highest=6.0

data_source=1
n_month=12
iregion=1      ;iregion=1, arctic; 2, antarctic; 2, global

openfile0=strarr(data_source)
openfile=strarr(data_source) & openfile1=strarr(data_source)
openfile0=['grid.dat']

openfile=['heff.H1985']

filestr2=' NCEP/NCAR Reanalysis Forcing'

nx=360 & ny=120 ;Arctic Parallel Ocean and sea Ice Model (POIM) grid
max_levels=30
desired_level =1
mask_value=99999.0
barwidth=200.0

x1=fltarr(nx,ny, /NOZERO)
y1=fltarr(nx,ny, /NOZERO)
x2=fltarr(nx, /NOZERO)
y2=fltarr(ny, /NOZERO)
data0=fltarr(nx,ny, /NOZERO)
data1=fltarr(nx,ny, /NOZERO)
data3=fltarr(nx,ny, /NOZERO)
data4=fltarr(nx,ny, /NOZERO) & idata=intarr(nx,ny, /NOZERO)
data2=fltarr(2,nx,ny, /NOZERO)
con_data=fltarr(data_source,n_month,nx,ny, /NOZERO)
vec_data=fltarr(data_source,n_month,2,nx,ny, /NOZERO)
kmt=intarr(nx,ny) & kmu=intarr(nx,ny) & kmt0=intarr(4)
lon=fltarr(nx,ny) & lat=fltarr(nx,ny)

; read grid info (lon., and lat.)
       close,9 & openr, 9, openfile0(0)
       readf, 9, lon, format='(10f8.2)'
       readf, 9, lat, format='(10f8.2)'
       close,9

mask_value=99999.0

filestr1=strarr(data_source)
filestr1(0)='Hi, '

month=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','Annual']
;month=string(format='(i3)',indgen(n_month)+1)

; Define color table.

;num_levels=fix((highest-lowest)/level_dif)+1
num_levels=fix((highest-lowest)/level_dif)
nconts=num_levels
print,'num_levels=',num_levels
; now set up color table
          r=intarr(nconts+3)
          g=intarr(nconts+3)
          b=intarr(nconts+3)
con_clr,nconts,r,g,b,white, black, brown, blue
;redblue,nconts,r,g,b,white, black, brown, blue

; Define the extent of the map.  Note that 1000 km is about 9 degrees of lat.
xmin = -2600
xmax =  2700
ymin = -2300
ymax =  2200
dx=40.0
cells=3.0

; read grid mask (ocean levels)
openr,9,'io.dat_360_120.output'
  readf, 9, kmt,format='(360i2)'
close,9
;print,kmt

for i=0,nx-2 do begin
for j=0,ny-2 do begin
    kmt0(0)=kmt(i,j)
    kmt0(1)=kmt(i+1,j)
    kmt0(2)=kmt(i,j+1)
    kmt0(3)=kmt(i+1,j+1)
    kmu(i,j)=min(kmt0)
endfor
endfor

;print,kmt

; define user symbol, use psym=8 for user defined symbol
aa=findgen(16)*(!pi*2/16.)
;usersym, cos(aa),sin(aa), /fill
usersym, 0.2*cos(aa),0.2*sin(aa)

window=fltarr(ncol*nrow,4)
ymod=0
xstart=0.25 & xend=0.95 & ystart=0.83 & yend=0.88
for i=0,ncol*nrow-1 do begin
  xmod=i mod ncol
  if xmod eq 0 then ymod=ymod+1
  window(i,0)=xmod*1.0/ncol + xstart/ncol
  window(i,1)=(nrow-ymod)*1.0/nrow + ystart/nrow
  window(i,2)=xmod*1.0/ncol + xend/ncol
  window(i,3)=(nrow-ymod)*1.0/nrow + yend/nrow
;  print, i,window(i,0),window(i,1),window(i,2),window(i,3)
endfor

; start big loop
ip=-1
for flag=0,data_source-1 do begin
       openr, flag+1, openfile(flag), /swap_if_big_endian
for count=0, n_month-1 do begin
       print, count

; read 2-D fields
       readu, flag+1, data4
       
    filestr=filestr1(flag) + month(count)

    data4(where (kmt lt 1))=99999.9

; Set up the plotting space with no margins or axes; draw a map
xmargin=0.5 & ymargin=0.3

if iregion eq 1 then begin

map_set, 90, 0, -90, limit=[49,-180,90,180],latdel=20, londel=90, /grid, /continents, /lamb, /advance,XMARGIN=xmargin, YMARGIN=ymargin

endif else if iregion eq 2 then begin

map_set, -90, 0, 0, limit=[-90,-180,-50,180],latdel=5, londel=90, /grid, /continents, /lamb, /advance

endif else begin

map_set, limit=[-90,-180,90,180],latdel=20, londel=20, /grid, /continents, /cyl, title=' ', /advance

endelse

; plot data in color
   idata=fix((nconts-1)*(data4-lowest)/(highest-lowest))
   for i=0,nx-1 do begin
   for j=0,ny-1 do begin
     if idata(i,j) lt 0 then idata(i,j)=0
     if idata(i,j) ge nconts then idata(i,j)=nconts-1
   endfor
   endfor
   for i=0,nx-1 do begin
   for j=0,ny-1 do begin
     if kmt(i,j) ge 1 then begin
       o=(i+1) mod nx
       p=(j+1) mod ny
       polyfill, color=idata(i,j), $
       [lon(i,j),lon(o,j),lon(o,p),lon(i,p)],[lat(i,j),lat(o,j),lat(o,p),lat(i,p)]
     endif
   endfor
   endfor

; draw a color bar
ip=ip+1 & ipp=ip mod (ncol*nrow) & if ipp eq 0 then ip=0
color_size=nconts & div=nconts/2
;colorbar, maxrange=highest,minrange=lowest,ncolors=color_size-1, $
colorbar, maxrange=highest,minrange=lowest,ncolors=color_size, $
;divisions=color_size, position=[0.05,0.40,0.08,0.95], /vertical
;divisions=color_size, position=[0.05,0.40,0.08,0.95], /vertical
;divisions=color_size, position=[0.05,0.40,0.08,0.95]
;divisions=color_size, position=[window(ip,*)]
divisions=div, position=[window(ip,0),window(ip,1),window(ip,2),window(ip,3)],charsize=1.4
;divisions=div, position=[window(ip,0),window(ip,1),window(ip,2),window(ip,3)],charsize=1.4, color=black, format='(f3.1)'

; Draw some characters
if iregion eq 1 then begin
xyouts, 90, 56, filestr, color=white, charsize=1.0, charthick=8, alignment=0.5
xyouts, 90, 56, filestr, color=black, charsize=1.0, charthick=1, alignment=0.5
endif else if iregion eq 2 then begin
endif

endfor
endfor
; end of big loops

close,1,2,3,4,5,6,7,8,9,10,99


; Close up.
if keyword_set(ps) then begin
   device, /close_file
   set_plot, 'X'
;  set_plot, 'CGM' ; if current graphics device is CGM by typing: help,/device
endif else begin
   hak, /mesg
   wdelete, win
endelse

!p.multi = 0
end

