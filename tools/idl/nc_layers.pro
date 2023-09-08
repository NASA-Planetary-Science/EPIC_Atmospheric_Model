;* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
;*                                                                 *
;* Copyright (C) 2003 Raul Morales-Juberias and Timothy E. Dowling *
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

;* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
;* This procedure allows you to view one time step of the model      *
;* layers in interactive mode					     *
;* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

pro nc_layers,t

@script

path = getenv('IDL_EPIC_PATH')+'/'

if N_Params() eq 0 then begin
   print,strcompress('Input time step to display [0,'+string(trange-1)+']:')
   read,t
endif

Loadct,39
ncolors=!d.table_size
device,decomposed=0

maxvar = max(var2plot)
minvar = min(var2plot)

model2  = obj_new('IDLgrModel',lighting=2)

zscale  = 50.

grid_globe_lontop = float(grid_globe_lontop)
grid_globe_lonbot = float(grid_globe_lonbot)
grid_globe_lattop = float(grid_globe_lattop)
grid_globe_latbot = float(grid_globe_latbot)

dlon = (grid_globe_lontop-grid_globe_lonbot)/(grid_ni-1)
dlat = (grid_globe_lattop-grid_globe_latbot)/(grid_nj)

image_data = bytarr(4,grid_ni,grid_nj)

;
; NOTE: Need to add the layers to the graphics model from the bottom up, because
;       IDL assumes this is the back-front order for transparency.
;
for k = grid_nk-1,0,-1 do begin
  surf = obj_new('IDLgrSurface',zscale*(var2plot(*,*,k,t)-minvar)/(maxvar-minvar),style=2,color=[255,255,255])
  surf->setProperty,xcoord_conv=[grid_globe_lonbot,dlon]
  surf->setProperty,ycoord_conv=[grid_globe_latbot,dlat]
  surf->setProperty,zcoord_conv=[0,zscale*1.5/grid_nk]
  surf->setProperty,shading=1  ; 1=Gourand shading

  ; Could use a second variable to define color through this image_data array.
  for i = 0,grid_ni-1 do begin
    for j = 0,grid_nj-1 do begin
      image_data(0,i,j) = 255
      image_data(1,i,j) = 200
      image_data(2,i,j) =  10
      ; set alpha channel for opacity
      if (k eq grid_nk-1) then begin
        ; make bottom layer opaque
        image_data(3,i,j) = 255
      endif else begin
        image_data(3,i,j) = 255*.75  
      endelse
    endfor
  endfor
  image = obj_new('IDLgrImage',data=image_data,blend_function=[3,4])
  surf->setProperty,texture_map=image

  model2->add,surf
endfor

xaxis = obj_new('IDLgrAxis',0,/exact,range=[grid_globe_lonbot,grid_globe_lontop])
xaxis->setProperty,location=[0,grid_globe_latbot]
xaxis->setProperty,tickinterval=90
model2->add,xaxis

yaxis = obj_new('IDLgrAxis',1,/exact,range=[grid_globe_latbot,grid_globe_lattop])
yaxis->setProperty,location=[grid_globe_lonbot,0]
yaxis->setProperty,tickinterval=30
model2->add,yaxis
	
;model2 -> rotate,[1,0,0],-90
;model2 -> rotate,[0,1,0],45
;model2 -> rotate,[1,0,0],30
				
xobjview,model2,/double,xsize=600,ysize=600,title=strcompress(planet_name+':  '+'timestep = '+string(t))
	
;***************************************************************************
;***************************************************************************

return
end
