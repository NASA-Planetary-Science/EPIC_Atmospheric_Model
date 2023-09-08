pro plot_YvsX_interpolate
;/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
; *                                                                 *
; * Copyright (C) 2008 Kunio M. Sayanagi                            *
; *                    Timothy E. Dowling                           *
; *                                                                 *
; * This program is free software; you can redistribute it and/or   *
; * modify it under the terms of the GNU General Public License     *
; * as published by the Free Software Foundation; either version 2  *
; * of the License, or (at your option) any later version.          *
; * A copy of this License is in the file:                          *
; *   $EPIC_PATH/License.txt                                        *
; *                                                                 *
; * This program is distributed in the hope that it will be useful, *
; * but WITHOUT ANY WARRANTY; without even the implied warranty of  *
; * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.            *
; *                                                                 *
; * You should have received a copy of the GNU General Public       *
; * License along with this program; if not, write to the Free      *
; * Software Foundation, Inc., 59 Temple Place - Suite 330,         *
; * Boston, MA  02111-1307, USA.                                    *
; *                                                                 *
; * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

; This program interpolates an input data on a specified grid and plots Y(X) ...
; The first line of the input data should have the number of data points in the file
; The subsequent lines have X Y data...

  DEFSYSV,'!YES',1
  DEFSYSV,'!NO', 0
 
  plot_num  = 0
   
  ; set up the PS graphics device
  mydevice = !D.NAME
  set_plot,'ps'
  device,/inches,xsize=6.5,ysize=9.,filename='idl.ps', xoffset=1.0,yoffset=1.0
  !p.multi=[0,1,1,0,0]

  print,      'Input File: '
  inputfile  = '' 
  read, inputfile 
  filename = inputfile
  get_lun, file_lun
  openr, file_lun, filename

  data_size = 0

  readf, file_lun, data_size

  X     = dblarr(data_size)
  Y     = dblarr(data_size)

  Xin = 0.0
  Yin = 0.0

  for kount=0, data_size-1 do begin
    readf, file_lun, Xin, Yin
    X[kount] = Xin
    Y[kount] = Yin
  endfor    

  free_lun, file_lun

  print, "Set max and min of output independent variable (leave empty to use max/min of input data)"
  print, "  Enter min then max in one line:"
  maxmin = "junk"
  read, maxmin
  if ( maxmin eq '') then begin
    Xlo = min(X)
    Xhi = max(X)        
  end else begin
    reads, maxmin, Xlo, Xhi
  endelse

  print, "Set the output grid spacing deltaX:"
  deltaX = 1.0
  read, deltaX

  numX = ceil((Xhi -Xlo)/deltaX)+1
  Xinterp = dindgen(numX)*deltaX +Xlo

  Yinterp = interpol(Y,       $ ; uninterpolated data
		     X,       $ ; grid of uninterpolated data
		     Xinterp, $ ; the grid on which interpolation is done
		     /lsquadratic)

  Xlo = min(Xinterp)
  Xhi = max(Xinterp)
  Ylo = min([Yinterp, Y])
  Yhi = max([Yinterp, Y])

  plot, X, Y, $
      linestyle =  0, $
      xstyle    = 1, $
      xrange    = [Xlo, Xhi], $
      yrange    = [Ylo, Yhi]

  oplot, Xinterp, Yinterp, linestyle = 1

  print, 'Output File: ' 
  outputfile = 'junk'
  read, outputfile
  if (outputfile eq '') then begin
    filename = 'output.dat'
  endif else begin
    filename = outputfile
  endelse

  get_lun, file_lun
  openw, file_lun, filename
  printf, file_lun, numX
  for kount=0, numX-1 do begin
    printf, file_lun, Xinterp[kount], Yinterp[kount]
  endfor
  free_lun, file_lun

  ; close device and finish up
  device, /close
  set_plot, mydevice
  return
end


