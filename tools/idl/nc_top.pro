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

pro nc_top
  
  @script

  path = getenv('IDL_EPIC_PATH')+'/'

  print,'Choose what variable do you want to render:'
  print,'1 - Potential Vorticity'
  print,'2 - Eddy Potential Vorticity'
  option2 = 3
  read,option2

  moviename='pp'
  print,'Name of the movie to be created:'
  read, moviename
  mpeg = obj_new('IDLgrMpeg',filename=strcompress(path+moviename+'.mpg'),frame_rate=2,format=0,Quality=100)
  
  kplot = 15
  print,'Central layer to plot:'
  read,kplot
  kplot  = kplot-1
  
  device,decomposed = 0
  Loadct,3
  ncolors=!d.table_size
  !P.background = 255
  !P.charsize = 0.0

  legend1 = strcompress('Planet: '+planet_name+'; '+planet_type)
  legend2 = strcompress('Grid info:'+grid_data_version+'; '+grid_geometry)
  legend3 = strcompress(grid_advection_scheme+'; '+grid_uv_timestep_scheme)
  legend4 = strcompress('Longitude: '+grid_globe_lonbot+'->'+grid_globe_lontop+' ('+grid_ni+')')
  legend5 = strcompress('Latitude: '+grid_globe_latbot+'->'+grid_globe_lattop+' ('+grid_nj+')')
  legend6 = strcompress('Sigmatheta: '+grid_sgth_bot+'->'+grid_sgth_top+' ('+grid_nk+')'+'; Sponge: '+grid_k_sponge)
  years   = 0.0 & days = 0.0 & hours = 0.0 & min = 0.0 & sec = 0.0

  device,get_screen_size = dev_size

  wysize  = dev_size[1]*.65 & yBmargin= 180. & yTmargin= 80.
  ysize   = (wysize-yBmargin-yTmargin)/5.
  xsize   = (ysize*grid_ni)/(grid_nj*1.0)
  if xsize lt dev_size[0]/2 then xsize=dev_size[0]/2.5
  xRmargin1= 60.& xLmargin2=30.& xRmargin2=xRmargin1+xsize & xLmargin1=xLmargin2+xsize/4.0
  wxsize   = xRmargin2+xLmargin1

  pim = fltarr(4,10)
  for i = 0,4 do begin
     pim[*,2*i]   = [xRmargin1, yBmargin+i*ysize, xRmargin1+xsize,     yBmargin+(i+1)*ysize]
     pim[*,2*i+1] = [xRmargin2, yBmargin+i*ysize, xRmargin2+xsize/4.0, yBmargin+(i+1)*ysize]
  endfor
  pbar1   = [(xRmargin1)/wxsize,(yBmargin+5*ysize+35)/wysize,(xRmargin1+xsize)/wxsize,(yBmargin+5*ysize+45)/wysize]

  Window,1,xpos=10,ypos=50,Xsize=wxsize,Ysize=wysize

  pvmean = fltarr(grid_nj+1,grid_nk)
  for k = 0, grid_nk-1 do begin
    for j = 0, grid_nj do begin
      pvmean[j,k] = mean(pv2[*,j,k,*])
    endfor
  endfor

  for t = 0,trange-1 do begin
  
    erase  
  
    aux   = time(t)/31536000.0   & years = floor(aux)
    aux   = (aux - years) * 365. & days  = floor(aux)
    aux   = (aux - days) * 24.   & hours = floor(aux)
    aux   = (aux - hours) * 60.  & min   = floor(aux)
    sec   = ((aux - min) *60.)               
    if (round(sec) eq 60) then begin
      min = min+1
      sec = round(sec)
    endif else begin
      sec = round(sec)
    endelse   
    timer = strcompress('Time-step:'+string(t)+' = '+string(years)+ $
            ' Years;'+string(days)+' Days; '+string(hours)+         $
	    ':'+string(min)+':'+string(sec))
 
    zonal   = fltarr(grid_nj+1,grid_nk)
    picture = fltarr(grid_ni,grid_nj+1,grid_nk)

    for k = 0, grid_nk-1 do begin
      for j = 0, grid_nj do begin
        zonal[j,k] = mean(u[*,j,k,t])
      endfor
    endfor

    case option2 of
      1: Begin
        for k = 0, grid_nk-1 do begin
          for j = 0, grid_nj do begin
            picture[*,j,k] = pv2[*,j,k,t]
          endfor
        endfor
      end
    
      2: Begin
        for k = 0, grid_nk-1 do begin
          for j = 0, grid_nj do begin
            picture[*,j,k] = pv2[*,j,k,t]-pvmean[j,k]
          endfor
        endfor
      end
    endcase

    pict_min = fltarr(grid_nk)
    pict_max = fltarr(grid_nk)
    for k = 0,grid_nk-1 do begin
      pict_min[k] = min(picture[*,*,k])
      pict_max[k] = max(picture[*,*,k])
    endfor

    for i = 0,4 do begin
       k2plot = kplot+2-i

       plot,zonal[*,k2plot],lat_u,/noerase,/device,position=pim[*,2*i+1],color=0,  $
       xstyle=1.0,    ystyle=1.0,	                                           $
       xtickformat = '(f1.0)', ytickformat ='(f1.0)',                              $
       xticks=1.0,    yticks=1.0,                                                  $
       xrange=[-100.,200.],  yrange=[min(lat_u),max(lat_u)]
           
       image=BytScl(picture[*,*,k2plot],min=pict_min[k2plot],max=pict_max[k2plot],top=ncolors-2)
       tv,congrid(image,xsize,ysize,interp=1,cubic=-.5),xRmargin1,yBmargin+i*ysize
         
       plot,lon_u,lat_u,/nodata,/noerase,/device,position=pim[*,2*i],color=0,     $
       xstyle=1.0,    ystyle=1.0,                                                 $
       xtickformat = '(f1.0)', ytickformat ='(f1.0)',                             $
       xticks=1.0,    yticks=1.0,                                                 $
       xrange=[min(lon_u),max(lon_u)], yrange=[min(lat_u),max(lat_u)]
       xyouts,xRmargin2+80,(i*ysize)+yBmargin+10,/device,strcompress('K = '+string(k2plot+1)),color=0,charsize=1.0  
    endfor
   
    plot,zonal[*,kplot-2],lat_u,/noerase,/device,/nodata,position=pim[*,1],color=0,  $
       xstyle=1,	                        ystyle=1,	                     $
       xrange=[-100.,200.],  yrange=[min(lat_u),max(lat_u)],                         $
       xcharsize=1.0,                           ycharsize=1.0,		             $
       xthick = 1.0,                            ythick = 1.0,                        $ 
       xtitle= 'U (m/s)'
    
    plot,lon_u,lat_u,/nodata,/noerase,/device,position=pim[*,0],color=0,             $
       xstyle=1,		       ystyle=1,		                     $
       xrange=[min(lon_u),max(lon_u)], yrange=[min(lat_u),max(lat_u)],               $
       xcharsize=1.0,		       ycharsize=1.0,	                             $
       xthick = 1.0,                   ythick = 1.0,                                 $
       xtitle='Longitude',             ytitle='Latitude'
    
    xyouts,0.1,0.14,/normal,legend1,color=0,charsize=1.5
    xyouts,0.1,0.12,/normal,legend2,color=0,charsize=1.0
    xyouts,0.1,0.10,/normal,legend3,color=0,charsize=1.0
    xyouts,0.1,0.08,/normal,legend4,color=0,charsize=1.0
    xyouts,0.1,0.06,/normal,legend5,color=0,charsize=1.0
    xyouts,0.1,0.04,/normal,legend6,color=0,charsize=1.0
    xyouts,0.1,0.02,/normal,timer,  color=0,charsize=1.0

    image=TVRD(order=1,true=1)
    mpeg -> put,image
  endfor

  print,'Saving the movie'
  mpeg -> save
  obj_destroy,mpeg
  wdelete	
	
  print,'DONE'

  return
  end
