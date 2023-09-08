;* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
;*                                                                 *
;* Copyright (C) 2005 Raul Morales-Juberias and T.E. Dowling       *
;*                                                                 *
;* This program is free software; you can redistribute it and/or   *
;* modify it under the terms of the GNU General Public License     *
;* as published by the Free Software Foundation; either version 2  *
;* of the License, or (at your option) any later version.          *
;* A copy of this License is in the file:                          *
;*   $EPIC4_PATH/License.txt                                       *
;*                                                                 *
;* This program is distributed in the hope that it will be useful, *
;* but WITHOUT ANY WARRANTY; without even the implied warranty of  *
;* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.            *
;*                                                                 *
;* You should have received a copy of the GNU General Public       *
;* License along with this program; if not, write to the Free      *
;* Software Foundation, Inc., 51 Franklin Street, Fifth Floor,     *
;* Boston, MA 02110-1301, USA.                                     *
;*                                                                 *
;* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

pro nc_held

path = getenv('IDL_EPIC_PATH')

filename = strcompress(path+'/extract.nc')
print,strcompress('Reading '+filename)
ncid = NCDF_OPEN(filename)

varid = NCDF_VARID(ncid,'time')
ncdf_diminq, ncid, varid,time_dim,trange
NCDF_VARGET, ncid,  varid, time
   NCDF_ATTGET, ncid,  varid, 'units', time_units
   time_units = STRING(time_units)

varid = NCDF_VARID(ncid,'lon_u')
NCDF_VARGET, ncid,  varid, lon_u      ; Read in variable 'lon_u'
   NCDF_ATTGET, ncid,  varid, 'units', lon_u_units
   lon_u_units = STRING(lon_u_units)

varid = NCDF_VARID(ncid,'lat_u')
NCDF_VARGET, ncid,  varid, lat_u      ; Read in variable 'lat_u'
   NCDF_ATTGET, ncid,  varid, 'units', lat_u_units
   lat_u_units = STRING(lat_u_units)
   NCDF_ATTGET, ncid,  varid, 'mapping', lat_u_mapping
   lat_u_mapping = STRING(lat_u_mapping)

NCDF_ATTGET, ncid, /GLOBAL, 'planet_cp', planet_cp
planet_cp = STRING(planet_cp)
NCDF_ATTGET, ncid, /GLOBAL, 'planet_rgas', planet_rgas
planet_rgas = STRING(planet_rgas)

NCDF_ATTGET, ncid, /GLOBAL, 'grid_nk', grid_nk
grid_nk = STRING(grid_nk)
NCDF_ATTGET, ncid, /GLOBAL, 'grid_nj', grid_nj
grid_nj = STRING(grid_nj)
NCDF_ATTGET, ncid, /GLOBAL, 'grid_ni', grid_ni
grid_ni = STRING(grid_ni)

x  = fltarr(grid_nk)
y  = fltarr(grid_nk)
p2 = fltarr(grid_nk)

kappa     = float(planet_rgas)/float(planet_cp)
kappa1    = kappa+1.
inv_kappa = 1./kappa

; Interpolate variables onto p/pbot vertical coordinate.
print,'Interpolate variables onto p/pbot vertical coordinate.'

n_interp = 24
sigma_trad    = fltarr(n_interp)
for k=0,n_interp-1 do begin
  sigma_trad[k] = float(k+1)/(n_interp+1) 
endfor

varid = NCDF_VARID(ncid,'p3')
NCDF_VARGET, ncid,  varid, p3      ; Read in variable 'p3'
   NCDF_ATTGET, ncid,  varid, 'units', p3_units
   p3_units = STRING(p3_units)

varid = NCDF_VARID(ncid,'ptop')
NCDF_VARGET, ncid,  varid, ptop      ; Read in variable 'ptop'
   NCDF_ATTGET, ncid,  varid, 'units', ptop_units
   ptop_units = STRING(ptop_units)

; A Held-Suarez extract file can be very large
; so we have to free memory whenever possible.

varid = NCDF_VARID(ncid,'u')
NCDF_VARGET, ncid,  varid, u      ; Read in variable 'u'
   NCDF_ATTGET, ncid,  varid, 'units', u_units
   u_units = STRING(u_units)

u_interp = fltarr(grid_ni,grid_nj,n_interp,trange)

for t=0,trange-1 do begin
  percent = 100.*float(t)/(trange-1)
  print, format='(%"U: %5.1f")',percent
  for i=0,grid_ni-1 do begin
    for j=0,grid_nj-1 do begin
      ; Form p2 from p3
      k = 0
      p2[k]=(((p3[i,j,k,t]^kappa1)-(ptop[i,j]^kappa1))/(kappa1*(p3[i,j,k,t]-ptop[i,j])))^inv_kappa
      for k=1,grid_nk-1 do begin
        p2[k]=(((p3[i,j,k,t]^kappa1)-(p3[i,j,k-1,t]^kappa1))/(kappa1*(p3[i,j,k,t]-p3[i,j,k-1,t])))^inv_kappa
      endfor

      for k=0,grid_nk-1 do begin
        x[k] = p2[k]/p3[i,j,grid_nk-1,t]
      endfor

      ; u interpolation
      for k = 0,grid_nk-1 do begin
        y[k] = u[i,j,k,t]
      endfor
      for k = 0,n_interp-1 do begin
        u_interp[i,j,k,t] = interpol(y,x,sigma_trad[k],/quadratic)
      endfor
    endfor
  endfor
endfor 

;free memory no longer needed
u = 0

varid = NCDF_VARID(ncid,'t2')
NCDF_VARGET, ncid,  varid, t2      ; Read in variable 't2'
   NCDF_ATTGET, ncid,  varid, 'units', t2_units
   t2_units = STRING(t2_units)

t2_interp = fltarr(grid_ni,grid_nj,n_interp,trange)

for t=0,trange-1 do begin
  percent = 100.*float(t)/(trange-1)
  print, format='(%"T2: %5.1f")',percent
  for i=0,grid_ni-1 do begin
    for j=0,grid_nj-1 do begin
      ; Form p2 from p3
      k = 0
      p2[k]=(((p3[i,j,k,t]^kappa1)-(ptop[i,j]^kappa1))/(kappa1*(p3[i,j,k,t]-ptop[i,j])))^inv_kappa
      for k=1,grid_nk-1 do begin
        p2[k]=(((p3[i,j,k,t]^kappa1)-(p3[i,j,k-1,t]^kappa1))/(kappa1*(p3[i,j,k,t]-p3[i,j,k-1,t])))^inv_kappa
      endfor

      for k=0,grid_nk-1 do begin
        x[k] = p2[k]/p3[i,j,grid_nk-1,t]
      endfor

      ; T2 interpolation
      for k = 0,grid_nk-1 do begin
        y[k] = t2[i,j,k,t]
      endfor
      for k = 0,n_interp-1 do begin
        t2_interp[i,j,k,t] = interpol(y,x,sigma_trad[k],/quadratic)
      endfor
    endfor
  endfor
endfor

;free memory no longer needed
t2 = 0

varid = NCDF_VARID(ncid,'theta2')
NCDF_VARGET, ncid,  varid, theta2      ; Read in variable 'theta2'
   NCDF_ATTGET, ncid,  varid, 'units', theta2_units
   theta2_units = STRING(theta2_units)

NCDF_CLOSE, ncid      ; Close the NetCDF file

theta2_interp = fltarr(grid_ni,grid_nj,n_interp,trange)

for t=0,trange-1 do begin
  percent = 100.*float(t)/(trange-1)
  print, format='(%"THETA2: %5.1f")',percent
  for i=0,grid_ni-1 do begin
    for j=0,grid_nj-1 do begin
      ; Form p2 from p3
      k = 0
      p2[k]=(((p3[i,j,k,t]^kappa1)-(ptop[i,j]^kappa1))/(kappa1*(p3[i,j,k,t]-ptop[i,j])))^inv_kappa
      for k=1,grid_nk-1 do begin
        p2[k]=(((p3[i,j,k,t]^kappa1)-(p3[i,j,k-1,t]^kappa1))/(kappa1*(p3[i,j,k,t]-p3[i,j,k-1,t])))^inv_kappa
      endfor

      for k=0,grid_nk-1 do begin
        x[k] = p2[k]/p3[i,j,grid_nk-1,t]
      endfor

      ; theta interpolation
      for k = 0,grid_nk-1 do begin
        y[k] = theta2[i,j,k,t]
      endfor
      for k = 0,n_interp-1 do begin
        theta2_interp[i,j,k,t] = interpol(y,x,sigma_trad[k],/quadratic)
      endfor
    endfor
  endfor
endfor 

;free memory no longer needed
theta2 = 0

; Calculate mean values.
print,'Calculate mean values.'

umean      = fltarr(grid_nj,n_interp)
theta2mean = fltarr(grid_nj,n_interp)
t2mean     = fltarr(grid_nj,n_interp)

for j=0,grid_nj-1 do begin
  percent = 100.*float(j)/(grid_nj-1)
  print, format='(%"Mean values: %5.1f")',percent
  for k=0,n_interp-1 do begin
    umean[j,k]      = mean(u_interp[*,j,k,*])
    theta2mean[j,k] = mean(theta2_interp[*,j,k,*])
    t2mean[j,k]     = mean(t2_interp[*,j,k,*])
  endfor
endfor

; Calculate eddy values.
print, 'Calculate eddy values.'

ueddy2 = fltarr(grid_ni,grid_nj,n_interp,trange)

for t=0,trange-1 do begin
  percent = 100.*float(t)/(trange-1)
  print, format='(%"UEDDY: %5.1f")',percent
  for k = 0,n_interp-1 do begin
    for j = 0,grid_nj-1 do begin
      for i = 0,grid_ni-1 do begin
        tmp             = u_interp[i,j,k,t]-umean[j,k]
        ueddy2[i,j,k,t] = tmp*tmp
      endfor
    endfor
  endfor
endfor

; Free memory no longer needed.
u_interp = 0

teddy2 = fltarr(grid_ni,grid_nj,n_interp,trange)

for t=0,trange-1 do begin
  percent = 100.*float(t)/(trange-1)
  print, format='(%"TEDDY: %5.1f")',percent
  for k = 0,n_interp-1 do begin
    for j = 0,grid_nj-1 do begin
      for i = 0,grid_ni-1 do begin
        tmp             = t2_interp[i,j,k,t]-t2mean[j,k]
        teddy2[i,j,k,t] = tmp*tmp
      endfor
    endfor
  endfor
endfor

; Free memory no longer needed.
t2_interp = 0

; Calculate Fourier transform of ueddy2.
print,'Calculate Fourier transform of ueddy2.'

spect = fltarr(grid_ni/2+1,grid_nj,n_interp,trange)
ps    = fltarr(grid_ni/2+1)

for t=0,trange-1 do begin
  percent = 100.*float(t)/(trange-1)
  print, format='(%"Spectrum of UEDDY: %5.1f")',percent
  for k = 0,n_interp-1 do begin
     for j = 0,grid_nj-1 do begin
       ps    = abs(FFT(ueddy2[*,j,k,t]))
       ps[0] = 0
       for i = 0,grid_ni/2 do begin
        spect[i,j,k,t] = ps[i]
       endfor
     endfor
  endfor
endfor

; Average eddy values.
print,'Average eddy values.'

teddy2mean = fltarr(grid_nj,n_interp)

for j = 0,grid_nj-1 do begin
  for k = 0,n_interp-1 do begin
    teddy2mean[j,k] = mean(teddy2[*,j,k,*])
  endfor
endfor

; Free memory no longer needed.
teddy2 = 0

spectmean = fltarr(grid_ni/2+1,grid_nj)

for j = 0, grid_nj-1 do begin
  for i=0,grid_ni/2 do begin
    spectmean[i,j] = mean(spect[i,j,*,*])
  endfor
endfor

; Free memory no longer needed.
spect = 0

; Latitude and wavenumber vectors

lat = fltarr(grid_nj)
for j=0,grid_nj-1 do begin
 lat[j] = lat_u[j]
endfor

pmst = lon_u[grid_ni-1]-lon_u[grid_ni-2]
wn   = findgen(grid_ni/2+1)

;----------------------------------------------

print,'Make plots.'

Loadct,39
ncolors=!d.table_size

set_plot,'ps'
device,filename='held_suarez.ps',preview=0
device,/landscape,/inches,xsize=9.0,ysize=7,/helvetica,font_size=12, $
       bits_per_pixel=1

!p.multi = [0,2,2,0,0]

!p.charthick=2. 

!x.omargin=[2,2]
!y.omargin=[5,1]

thickv    = fltarr(20)+1.0
thickv[5] = 2.0

;zonal mean wind

!x.margin = [6,1]
!y.margin = [0,1]

contour,umean,lat,sigma_trad,color=0,/device,/follow,  $
	levels=findgen(20)*4-20,		            $
	c_labels=[1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0], $
	c_linestyle=(findgen(20)*5-25 lt 0.0),	            $
	c_thick =thickv,                                    $ 
	c_charthick=2.0, c_charsize=0.7,                    $
        charthick=2.0,                                      $
	xstyle=1,			ystyle=1,           $
	xrange=[-90,90],		yrange=[1,0],       $
	xticks=6,	                yticks=4,           $
        xtickname=[' ',' ',' ',' ',' ',' ',' '],            $
        xthick=2.0,ythick=2.0,                              $  
	ytitle='p/pbot'
plots,[0,0],[0,1],color=0
xyouts,-86,.95,'a)',charsize=2.,charthick=3.

!x.margin = [1,6]
!y.margin = [0,1]
       
;eddy variance of the temperature

thick2 = fltarr(20)+1.0

contour,teddy2mean,lat,sigma_trad,color=0,/device,/follow,	$
        levels = findgen(20)*5,                                 $
        c_labels=[1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0],     $
        c_thick=thick2,                                         $
	c_charthick=2.0, c_charsize=0.7,                        $
        charthick=2.0,                                          $
	xstyle=1,			ystyle=1,	        $
	xrange=[-90,90], 		yrange=[1,0], 	        $
	xticks=6,	                yticks=4,               $   
        xtickname=[' ',' ',' ',' ',' ',' ',' '],                $
        xthick=2.0,ythick=2.0,                                  $  
        ytickname=[' ',' ',' ',' ',' ']         
plots,[0,0],[0,1],color=0
xyouts,-86,.95,'b)',charsize=2.,charthick=3.


;zonal mean temperature

!x.margin = [6,1]
!y.margin = [0,1]

thick3 = fltarr(24)+1.0
	
contour,t2mean,lat,sigma_trad,color=0,/device,/follow,     $	
	levels=findgen(24)*5+190,	                   $
        c_thick = thick3,                                  $
	c_charthick=2.0, c_charsize=0.7,                   $
        charthick=2.0,                                     $
	xstyle=1,			ystyle=1,	   $
	xrange=[-90,90],		yrange=[1,0], 	   $
	xticks=6,	                yticks=4,          $
        xtickname=['-90 ','-60 ','-30 ','0','30','60',' '],$
        ytickname=['1.00','0.75','0.50','0.25',' '],       $  
        xthick=2.0,ythick=2.0,                             $            
	xtitle='Latitude [deg]',    ytitle='p/pbot'
plots,[0,0],[0,1],color=0
xyouts,-86,.95,'c)',charsize=2.,charthick=3.


;potential temperature

!x.margin = [1,6]
!y.margin = [0,1]

thick4 = fltarr(20)+1.0
	
contour,theta2mean,lat,sigma_trad,color=0,/device,/follow,  $
	levels = [260,265,270,275,280,285,290,295,300,305,310,315,320,325,350,400,500,600,700,800],$
        c_thick = thick4,                                   $
	c_charthick=2.0, c_charsize=0.7,	            $
        charthick=2.0,                                      $
	xstyle=1,			ystyle=1,           $
	xrange=[-90,90],		yrange=[1,0],	    $
	xticks=6,	                yticks=4,           $
        xtickname=['-90 ','-60 ','-30 ','0','30','60','90'],$              
        ytickname=[' ',' ',' ',' ',' '],                    $  
        xthick=2.0,ythick=2.0,                              $             
	xtitle='Latitude [deg]'
plots,[0,0],[0,1],color=0
xyouts,-86,.95,'d)',charsize=2.,charthick=3.


; Fourier spectrum of ueddy2.

!p.multi=0

thick5 = fltarr(20)+1.0

contour,spectmean,wn,lat,color=0,/normal,/follow,                  $
        levels = (findgen(20)+1)*2,                                $
        c_thick = thick5,                                          $
        c_labels = [0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1],      $
	c_charthick=2.0, c_charsize=0.7,                           $
        charthick=2.0,                                             $
	xstyle=1,			ystyle=1,	           $
        xrange=[0,15], 		        yrange=[-90,90],           $
	xticks=3,	                yticks=4,                  $  
        xthick=2.0,ythick=2.0,                                     $                          
	xtitle='Zonal Wavenumber',      ytitle='Latitude [deg]'
plots,[0,15],[0,0],color=0

device,/close_file

print,'Done processing extract frames ',0,' to ',trange-1

return
end
