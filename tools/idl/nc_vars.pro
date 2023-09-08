;* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
;*                                                                 *
;* Copyright (C) 2003 Raul Morales-Juberias                        *
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

pro nc_vars

@script

path = getenv('IDL_EPIC_PATH')+'/'

Loadct,39
ncolors=!d.table_size
device,decomposed=0

!p.background = 255
!p.multi =[0,2,2,0,0]
!y.omargin=[1,15]

legend1 = strcompress('Planet: '+planet_name+'; '+planet_type)
legend2 = strcompress('Grid info:'+grid_data_version+'; '+grid_geometry)
legend3	= strcompress('         '+grid_advection_scheme+'; '+grid_uv_timestep_scheme)
legend4 = strcompress('Longitude: '+grid_globe_lonbot+'->'+grid_globe_lontop+' ('+grid_ni+')')
legend5 = strcompress('Latitude: '+grid_globe_latbot+'->'+grid_globe_lattop+' ('+grid_nj+')')
legend6 = strcompress('Sigmatheta: '+grid_sgth_bot+'->'+grid_sgth_top+' ('+grid_nk+')'+'; Sponge: '+grid_k_sponge)

print,strcompress('grid_ni = '+grid_ni)
iindex=0
print,'Input i index to plot'
read,iindex

moviename=''
print,'Name of the movie to be created:'
read, moviename
mpeg   = obj_new('IDLgrMpeg',filename=strcompress(path+moviename+'.mpg'),frame_rate=2,format=0,Quality=100)

years  = 0.0 & days = 0.0 & hours = 0.0 & min = 0.0 & sec = 0.0
min_u      = -100 & max_u=100
min_v      = -100 & max_v=100
min_theta  = min(theta) & max_theta=max(theta)

p3         = p3/100  ;pa -> mb
min_p3     = min(p3) & max_p3=max(p3)

;########################################################	
window,xsize=600,ysize=600,retain=2
erase
;--------------------------------------------------------		
for t=0,trange-1 do begin  ;TEMPORAL LOOP
;--------------------------------------------------------
	aux = time(t)/31536000.0    & years = floor(aux)
	aux = (aux - years) * 365.  & days  = floor(aux)
	aux = (aux - days) * 24.    & hours = floor(aux)
 	aux = (aux - hours) * 60.   & min   = floor(aux)
	sec = ((aux - min) *60.)               
	if (round(sec) eq 60) then begin
   		min = min+1
	        sec = round(sec)
	endif else begin
   		sec = round(sec)
	endelse   
	timer = strcompress('Time-step:'+string(t)+' = '+ $
		string(years)+' Years;'+string(days)+' Days; '+string(hours)+ $
		':'+string(min)+':'+string(sec))
	!p.charsize=1.0
	;################################################
	plot,lat_u,u(iindex,*,0,t),			$
	/device,color=0,linestyle=1,thick=1.0,		$
	xstyle=1,		ystyle=1,		$
	xrange=[grid_globe_latbot,grid_globe_lattop],	$
	yrange=[min_u,max_u],		 		$
	xcharsize=1.0,		ycharsize=1.0,		$
	xtitle='Latitude',      ytitle='U(m/s)'
	for k = 0,grid_nk-1 do begin
		oplot,lat_u,u(iindex,*,k,t),color=k*(254/(grid_nk-1))
	endfor
	;################################################
	plot,lat_v,v(iindex,*,0,t),			$
	/device,color=0,linestyle=1,thick=1.0,		$
	xstyle=1,		ystyle=1,		$
	xrange=[grid_globe_latbot,grid_globe_lattop],	$
	yrange=[min_v,max_v], 				$
	xcharsize=1.0,		ycharsize=1.0,		$
	xtitle='Latitude',      ytitle='V(m/s)'
	for k = 0,grid_nk-1 do begin
		oplot,lat_v,v(iindex,*,k,t),color=k*(254/(grid_nk-1))
	endfor
	;################################################
	plot,lat_u,theta(iindex,*,0,t),			$
	/ylog,/device,color=0,linestyle=1,thick=1.0,	$
	xstyle=1,		ystyle=1,		$
	xrange=[grid_globe_latbot,grid_globe_lattop],	$

        yrange=[min_theta*0.98,max_theta*1.02],		$
        yticks = 12,					$
	xcharsize=1.0,		ycharsize=1.0,		$
	xtitle='Latitude',      ytitle='theta(K)'
	for k = 0,grid_nk-1 do begin
		oplot,lat_u,theta(iindex,*,k,t),color=k*(254/(grid_nk-1))
		oplot,lat_u,theta(iindex,*,k,t),color=k*(254/(grid_nk-1)),linestyle=1
	endfor
	;################################################
	plot,lat_u,p3(iindex,*,0,t),			$
	/ylog,/device,color=0,linestyle=1,thick=1.0,	$
	xstyle=1,		ystyle=1,		$
	xrange=[grid_globe_latbot,grid_globe_lattop],	$
	yrange=[max_p3*1.1,min_p3*0.9],	                $
	xcharsize=1.0,		ycharsize=1.0,		$
	xtitle='Latitude',      ytitle='Pressure(mb)'
	for k = 0,grid_nk-1 do begin
		oplot,lat_u,p3(iindex,*,k,t),color=k*(254/(grid_nk-1))
		oplot,lat_u,p3(iindex,*,k,t),color=k*(254/(grid_nk-1)),linestyle=1
	endfor
	
	!p.charsize=0.0
	xyouts,0.1,0.98,/normal,legend1,color=0,charsize=2
	xyouts,0.1,0.96,/normal,legend2,color=0
	xyouts,0.1,0.94,/normal,legend3,color=0
	xyouts,0.1,0.92,/normal,legend4,color=0
	xyouts,0.1,0.90,/normal,legend5,color=0
	xyouts,0.1,0.88,/normal,legend6,color=0
	xyouts,0.1,0.86,/normal,timer,color=0

	image=TVRD(order=1,true=1)

	mpeg -> put,image

;	image_name=strcompress(moviename+string(t)+'.tiff')
;	write_tiff,image_name,image
endfor

print,'Saving the movie'
	mpeg -> save
	obj_destroy,mpeg
wdelete		
print,'DONE'
;***************************************************************************

end
