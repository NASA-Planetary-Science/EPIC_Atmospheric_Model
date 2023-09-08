;/* * * * * -*-Text-*-* * * * * * * * * * * * * * * * * * * * * * * *
; *                                                                 *
; * Copyright (C) 1999 Adam P. Showman                              *
; *                                                                 *
; * This program is free software; you can redistribute it and/or   *
; * modify it under the terms of the GNU General Public License     *
; * as published by the Free Software Foundation; either version 2  *
; * of the License, or (at your option) any later version.          *
; * A copy of this License is in the file:                          *
; *   $EPIC4_PATH/License.txt                                       *
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

; ncread.pro           Adam Showman
;                      October 1999
;
; ==================================================================
; Modified for use with EPIC 4.x by Kunio M. Sayanagi
;                      June 2008.
; ==================================================================
; 
; Requires IDL version 5.0 or later.
;
; Reads data from either a single-timestep netCDF file created by EPIC4 
; or an extract.nc and allows various means of plotting the data.  
;
; When reading an extract.nc file, the program asks which frame to read.
; Note that this frame number is different from EPIC's timestep number,
; and specifying frame_id = i just means reading the data that happens 
; to be stored in the i-th place in extract.nc data.
;
; Quantities to be plotted can either be 2-d arrays (corresponding to one
; of the 3-d scalar fields plotted at a given vertical layer) or
; 1-d arrays (corresponding to a "slice" through one of the
; 3-d scalar fields with either lon, lat, or theta [or p] as a coordinate;
; for such a plot, the values of the other two indep variables
; are fixed at index values obtained from the user).  The plots of 2-d 
; arrays can be done either with greyscale, contours, overlain 
; greyscale and contour plots, or wiremesh.  Plots of velocity
; can also be made: these plot velocity arrows corresponding to
; the horizontal wind (u,v) on a given vertical layer.
;
; When plotting a 1-d quantity along pressure coordinate, the program,
; of course, assumes that pressure (i.e., p3) is saved in the netCDF file.
; In EPIC, Pressure p3 is calculated on layer interfaces where as u and v
; are calculated at the center (see Dowling et al. 1998 for the details)
; -- this program calculates the layer center pressure to properly place
; the pressure corresponding to each altitude.
;
; After each plot, the user is asked whether additional plots should
; be overlain on the plot just created.  For example, wind vectors can be
; overlain on a greyscale pressure plot.  Overplotting
; with greyscale destroys the pre-existing plot, so for plots of
; two overlain quantities (one with greyscale and one with contours),
; do the greyscale first.  Overplotting can also be done for the
; 1-d plots.  It only makes sense if the type of plots are the same,
; and I have not yet added the relevant error checking machinery.
; (So it will let you overplot a 1-d plot on a 2-d plot right now
; even though that makes no sense.)  I have also not yet added
; overplotting for the 3-d wire mesh plotting routine. 
;
; The program asks the user whether the plot should go to screen or
; to a postscript file. 
; 
; The fields read from the netCDF file do not have pads.  I add pads
; to the variables to make it easier to calculate derivatives etc.
; However, that part of the code is not completely finished, so the
; pads are not done consistently in all cases.  This affects the 
; plots along the edges (for a plot of a 2-d array).

; To execute: from the IDL prompt, type ".run ncread.pro" to compile
;             the program, and then "main" to run it.



;*****************************************************************
;*****************************************************************
; MAIN PROGRAM
;*****************************************************************
;*****************************************************************
pro main
   ;---Setup system variables, default filename, etc.
   DEFSYSV,'!YES',1
   DEFSYSV,'!NO',0
   DEFSYSV,'!TRUE',1
   DEFSYSV,'!FALSE',0
   DEFSYSV,'!NEW',0
   DEFSYSV,'!OVERLAY',1
   DEFSYSV,'!EPIC4',1

   ;----Set the version number-------------
   DEFSYSV,'!EPIC_VERSION',!EPIC4

   ; Ask whether the file being read is an extract file and set the mode
   ans = 'junk'
   print, 'Reading an extract.nc file? [y/n]: '
   read, ans

   if ( (ans eq 'y') or (ans eq 'yes') )then begin
     DEFSYSV, '!EXTRACT_MODE', !YES
   endif else begin
     DEFSYSV, '!EXTRACT_MODE', !NO
   endelse

   filename = ''
   print,'Enter the filename [epic.nc]'
   read,filename
   if (filename eq '') then filename = 'epic.nc'

   ; Ask the time frame number in extract 
   ; (note that it's different from time step number!
   frame_id = -1 ; this is a dummy line to set the data type
   if (!EXTRACT_MODE) then begin
     print, 'Enter the frame number to extract from extract.nc'
     read, frame_id
   endif else begin
     ; if not reading extract.nc, keep frame_id < 0.
     frame_id = -1
   endelse

   ;-----get the file ID and # of dims,vars, etc.-----------

  file_id = NCDF_OPEN(filename, /Nowrite)  ; opens the file and returns its ID
                           ; open for reading only
  print,'Opened file ',filename

  basic = NCDF_INQUIRE(file_id) ; returns a structure containing
                         ; info about the netCDF file.
                         ; The info is a structure containing
			  ; ndims,nvars,ngatts,and recdim:
			  ; the number of
			  ; dimensions, variables, global
			  ; attributes and whether there's
                         ; an unlimited dimension for
                         ; this netCDF file, in that
                         ; order.  (If there is an unlimited
                         ; dimension, it returns its ID #, 
                         ; otherwise returns -1.)

  get_var_info,file_id,basic
  get_gatt_info,file_id,basic

  ;NCDF_CONTROL,id,/ENDEF

  ;---For some reason we're in data mode now already although I
  ; never asked IDL to switch.  Not sure why....



  ;%%%%%%%%%%%%%%%%%% NOW HAVE THE COMMAND LINE %%%%%%%%%%%%%%%%%%%%%%%%


  ;-----Query for output to screen or postscript file---------------
  ans='' & postscript = !NO
  print,'Output plots to ps file ("y" or "yes" for yes, otherwise no)?'
  read,ans
  if ((ans eq 'y') or (ans eq 'yes') or (ans eq 'Y') or (ans eq 'YES')) then begin
     postscript = !YES
     print,'Plots will be written to postscript file'
     old_dname = !D.NAME
     set_plot, 'ps'
     device, xoffset=0.5,yoffset=0.7,xsize=7.,ysize=9.7,/inches
  endif

  nplots=0
  print,'How many plots on one page/screen (maximum of 8)?'
  read,nplots

  if (nplots eq 1) then !p.multi=[0,1,1,0,0]
  if (nplots eq 2) then !p.multi=[0,1,2,0,0]
  if (nplots eq 3) then !p.multi=[0,1,3,0,0]
  if (nplots eq 4) then !p.multi=[0,2,2,0,0]
  if ((nplots eq 5) or (nplots eq 6)) then !p.multi=[0,2,3,0,0]
  if ((nplots eq 7) or (nplots eq 8)) then !p.multi=[0,2,4,0,0]

  if (nplots gt 8) then begin
     print,'Program not currently set up to handle more than 8 plots.'
     stop
  endif

  ;-----Loop through the plots--------------------------------------
  overlay = !NEW
  count = 1
  vert_level = 0
  while count le nplots do begin
     print,'Plot ',count,    $
       ': Quantity to plot (u, v, velocity, [or, any other quantity you know you have in extract.nc])?'
     quantity_name=''
     read,quantity_name

     numcontours = 20   ; number of contours between min and max for
			 ; contour plots

     if (quantity_name eq 'velocity') then begin
        vert_level=0
        print,'What vertical level to plot?'
        read,vert_level
     endif else begin
        plot_type=0
        print,'Type of plot (0=contour, 1=greyscale, 2=overlain contour/greyscale,'
        print,' 3=3-d wire-mesh, 4=one-dimensional slice)?'
        read,plot_type

        if (plot_type eq 0 or plot_type eq 1 or plot_type eq 2 or plot_type eq 3) $
	          then begin
	    plot_axis = 'dummy'
           vert_level=0
           print,'What vertical level to plot?'
           read,vert_level

	 endif else begin ; i.e., if (plot_type eq 4)
	    plot_axis = ''
	    print,'Type the axis to plot [lon, lat, sigmatheta, p]:'
	    read,plot_axis
	    if (plot_axis ne 'lon' and plot_axis ne 'lat' and $
	            plot_axis ne 'sigmatheta' and plot_axis ne 'p') then begin
		print,'Error: You did not enter lon, lat, sigmatheta, or p'
		stop
	    endif
	    I_plot = -1 & J_plot = -1 & K_plot = -1
	    if (plot_axis eq 'lon' or plot_axis eq 'lat') then begin
	       print,'Value of K (vertical index) to hold constant for the plot?'
	       read,K_plot
	    endif
	    if (plot_axis eq 'lon' or plot_axis eq 'sigmatheta' or         $
	             plot_axis eq 'p') then begin
	       print,'Value of J (latitude index) to hold constant for the plot?'
	       read,J_plot
	    endif
	    if (plot_axis eq 'lat' or plot_axis eq 'sigmatheta' or      $
	             plot_axis eq 'p') then begin
	       print,'Value of I (longitude index) to hold constant for the plot?'
	       read,I_plot
	    endif
        endelse
     endelse


     ;----- read data out of netCDF file -----------


     ;--create array of variables to extract from the netCDF file.
     ; deal with the special cases (i.e., velocity. pressure as vertical coordinate etc) 
     ; before specifying what needs to be read out of the netCDF file.
     ; How it's done is admittedly apologetically kludgie
     if (quantity_name eq 'velocity') then begin
       var_list = ['u', 'v']
       get_var_data,file_id, frame_id, var_list, var
       get_av_winds, file_id,var  ; uses var.u and var.v (assumed to exist) and averages
                          ; them onto h grid.  Places them into var.u_av and var.v_av
	 	           ; (which are added to the var structure),
		           ; and creates 1-d coordinate arrays var.lon_( ), var.lat_( ),
		           ; and var.sigmatheta_(  ).



       get_var_data,file_id, frame_id, var_list, var
     endif else if ( plot_axis eq 'p' ) then begin
       var_list = [quantity_name, 'p3']
       get_var_data,file_id, frame_id, var_list, var	
       get_layer_p, file_id, var
     endif else begin
       var_list = [quantity_name]
       get_var_data,file_id, frame_id, var_list, var	
     endelse

     numcontours = 20   ; number of contours between min and max for
			 ; contour plots


     if (quantity_name eq 'velocity') then plot_flat_horiz_vector,var,'u_av','v_av',vert_level,1,overlay,-40.,0.
     if (plot_type eq 0) then plot_flat_horiz_scalar,file_id,var,quantity_name,  $
		 overlay, !YES,!NO,numcontours,vert_level,1,1,postscript
     if (plot_type eq 1) then plot_flat_horiz_scalar,file_id, var,quantity_name, $
	          overlay,!NO,!YES,numcontours,vert_level,1,1,postscript
     if (plot_type eq 2) then plot_flat_horiz_scalar,file_id, var,quantity_name, $
		  overlay,!YES,!YES,numcontours,vert_level,1,1,postscript
     print,'main: vert_level=',vert_level
     if (plot_type eq 3) then plot_surface_horiz_scalar,var,quantity_name, $
	         vert_level,1,1
     if (plot_type eq 4) then plot_1dslice_scalar,file_id, var,quantity_name, $
		 plot_axis, I_plot, J_plot, K_plot,overlay

     over_plot = ''
     print,' Overplot an additional quantity on this plot [''y'' for yes, otherwise no]?'
     read,over_plot
     if over_plot eq 'y' then begin
        overlay = !OVERLAY
	 count = count - 1
     endif else begin
        overlay = !NEW
     endelse
  count = count + 1
  endwhile


  ;----re-set plot type if we send plots to postscript file
  if (postscript eq !YES) then begin
     device, /close
     set_plot,old_dname
  endif

  NCDF_CLOSE,file_id

end   ;----program

;***********************************************************************
;***********************************************************************



;***********************************************************************
; get_var_info: Gets info about the dependent and independent variables
;   from a netCDF file created by EPIC.  Right now nothing is done
;   with this info (I had it printing out the name and number of array
;   elements in the various dependent and independent variables, and
;   also their data type, the number of dimensions they have, and the
;   number of "attributes" they have [e.g. what their units are etc.].
;   This info is useful for understanding what a new netCDF file has in
;   it, but the print statements are currently commented out.)
;***********************************************************************
PRO get_var_info,id,basic
         ; id (input) is the ID # of the open netCDF file
	  ; basic (input) is a structure containing some basic info:
         ;     the number of "dimensions" (i.e. the number of possible
         ;     independent variables--these are specified for each dependent
         ;     variable) and the names of these  dimensions.

  ;-------dimension arrays for reading in stuff-------------
  dimname = strarr(basic.ndims)
  dimsize = intarr(basic.ndims)

  ;-------read in stuff----------------------------------
;   print,'Name and size of dimensions:'
  for i=0, basic.ndims-1 do begin
     NCDF_DIMINQ,id,i,sname,ssize  ; s for scalar--it seems
                                ; to not want an element
                                ; of an array for some reason.
;      print,sname,ssize
     dimname[i]=sname
     dimsize[i]=ssize
  endfor


  ;---Define the structure type to be returned by NCDF_VARINQ,
  ;    with one change: NCDF_VARINQ returns the .dim entry to be
  ;    an array of length .ndims.  This varies, so it would be
  ;    impossible to make an array of such structures.  For our
  ;    application, ndims has a maximum of 4, so I'll just define
  ;    .dim to be a long array of length 4 rather than of length
  ;    ndims.
  template = {name:"",datatype:"",ndims:0L,natts:0L,dim:lonarr(4)}

			    ; i.e. the structure has fields .name,
			    ; .datatype, .ndims, .natts,
			    ; and .dim (which is an array
			    ; giving the IDs for each
			    ; of the dimensions).

  template.dim=[-1,-1,-1,-1]  ; set the .dim array to a flag which cannot
                             ; be used as an ID #. (This will distinguish
			      ; between placefiller spots in the array
			      ; and real dimension ID values.)

  ;...copy the template over to make an array of structures
  varstruct = REPLICATE(template,basic.nvars) ; make an array
                                   ; of structures of type
				    ; result; array is
				    ; basic.nvars-1 long.

;   print,'Name,   datatype,   # dimensions, # attributes of variables'
  for i=0,basic.nvars-1 do begin
     result = NCDF_VARINQ(id,i) ;Get the array of variable info from netCDF file
                             ;The array is of type template except for an
                             ;array of possibly different length.

     varstruct[i].name = result.name
     varstruct[i].datatype = result.datatype
     varstruct[i].ndims = result.ndims
     varstruct[i].natts = result.natts
;      print,result.name,result.datatype,result.ndims,result.natts
     for j=0,result.ndims-1 do begin
         varstruct[i].dim[j]=result.dim[j]
     endfor
  endfor

  return 
end    ;-----end procedure get_var_info

;****************************************************************
; get_gatt_info: Get info about the global attributes (e.g. quantities
;    such as the planetary rotation rate, gravity, etc.).  Originally
;    I had it printing all of them to screen, but the print statements
;    are currently commented out.
;****************************************************************
pro get_gatt_info,id,basic

  template = {datatype:"",length:0L}      ; see IDL online notes for NCDF_ATTINQ
                 ; datatype is e.g. double or long--i.e. attribute data type
		  ; length is the number of values in the attribute.  For
		  ; strings this is number of characters plus 1 (for the cr at
		  ; the end).

  gattstruct = REPLICATE(template,basic.ngatts)       ; make array of structures
  gattname = strarr(basic.ngatts)   ; array of strings (for names of global att's)

;   print,'Name, datatype, and # elements for global attributes'
  for i=0,basic.ngatts-1 do begin
     tname = NCDF_ATTNAME(id,/GLOBAL,i)
     tresult = NCDF_ATTINQ(id,/GLOBAL,tname)
     gattname[i]=tname
;      print,tname,tresult.datatype,tresult.length
     ; for some reason IDL won't let me say gattstruct[i]=tresult, so I have
     ; to set each of the fields individually
     gattstruct[i].datatype=tresult.datatype
     gattstruct[i].length=tresult.length

  endfor

  return 
end    ;----end procedure get_gatt_info

;***********************************************************************
; get_dim_data: Gets the values for the dimensions
;***********************************************************************
PRO get_dim_data,id,basic
         ; id (input) is the ID # of the open netCDF file
	  ; basic (input) is the

  return   
end   ;---end procedure get_dim_data
;***********************************************************************
; get_var_data: Read in variable data and place it in a structure for
;   later use, as follows:
;        (1) Read in the values of selected dependent variables from the
;   netCDF file (which will generally by 3- or 4-d arrays for EPIC),
;   as well as 1-d arrays corresponding to their axes.   (These arrays,
;   e.g. lon_u, lat_u, and sigmatheta_u for dependent variable u, map 
;   indeces i,j, and k to values of longitude, latitude, and sigmatheta.)
;   Each variable has its *own* spatial 1-d arrays, since lon/lat/sigmatheta
;   values for a given index are *different* for u, v, p, and q (because
;   of the staggered C grid).  However, there is only one 1-d time array,
;   since all variables are read out at given times.  
;         (2) Transfer these arrays to a structure var, adding pads as
;   follows: 2 pads are added to the longitude index.  The i index therefore 
;   runs from 0 to ni+1, and the data read from the netCDF file is placed in 
;   elements 1 to ni.  Data from element 1 is transfered to element ni+1 and 
;   data from element ni is transfered to element 0.  (This is because the 
;   domain is always periodic in longitude.)  (The size of the netCDF file's
;   arrays in longitude run from 1 to ni, i.e. they contain ni elements.)
;   Either one or two pads are added to nj, depending on the boundary condition. 
;   For a channel boundary condition (which includes doing the whole globe), 
;   the arrays from the netCDF file run from 0 to nj in latitude, i.e. they 
;   contain nj+1 elements.  For that case we add one pad to the north.  (This pad 
;   is truly appropriate for q and v, because for those cases, the pad is the north 
;   pole or the north channel boundary, and its values might be needed for 
;   computation of other quantities [say divergence at the half-grid spacing south 
;   of that v value].  However, for u and h, the pad really should not be added,
;   because for u and h the pad value is _north_ of the north pole/north
;   channel boundary [AND more importantly, it is never needed for calculation].
;   However, to simplify things, I add a pad even for this case.)
;   For a periodic north-south boundary condition, the arrays from the netCDF file
;   run from 1 to nj (i.e., nj elements), and I add two pads, one to the north
;   and one to the south.  Data from element 1 is transferred to element nj+1 and
;   data from element nj is transferred to element 0.
;
;   Therefore, with the pads, the arrays will ALL contain ni+2 elements in
;   lon and nj+2 elements in lat.  They will ALL run from 0-->ni+1 and 0-->nj+1.
;   This will be true independent of the variable (i.e. q vs. v vs. h vs. u) and 
;   geometry (periodic vs. channel in N-S direction).  This makes it easy
;   to have for loops just run over these quantities without checking specific
;   cases (as would have been necessary had I not added a pad to the u and h
;   variables for the channel/globe case).  To prevent problems in finding mins/
;   maxes etc., for that case, I copy the northernmost real u and h values
;   (i.e., at index nj) into the pad values at nj+1.
;
;       The function is intended to work with both EPIC 3.8 (or 3.9) and
;   EPIC 4.0.  The independent height variable is called "sigmatheta" here,
;   but it could either be actual sigmatheta (from EPIC 4.0) or theta
;   (from EPIC 3.9).  Therefore, fields in structure var include var.p,
;   var.theta, var.u, var.sigmatheta_u, var.sigmatheta_v, var.sigmatheta_theta,
;   etc.  For EPIC 4.0 these are self explanatory.  For EPIC 3.9, we simply
;   treat the indep. variable theta as a special case of sigmatheta.  The
;   var.theta should therefore have precisely the same theta value at a given
;   height index k. 
;***********************************************************************
PRO get_var_data, file_id, frame_id, vars_to_get, var
         ; file_id (input) is the ID # of the open netCDF file
	  ; frame_id (input) is the frame number in extract.nc to extract.
	  ;     Note that it is DIFFERENT from EPIC timestep number,
	  ;     and frame_id = i is just the i-th frame that happens to be
	  ;     saved in extract.nc
	  ; vars_to_get (input) is a string array containing a list of
	  ;     each variable to retrieve from the netCDF file and put into
	  ;     structure var.  E.g. ['u', 'v', 'p'].  It need not include 
	  ;     the axes for each of these variables, which read in and added to 
	  ;     var automatically. 
	  ; var (output) is a structure containing the indpendent and
         ;     dependent variables.  Need not be predefined in the 
	  ;     calling routine--it is defined here and then returned.
;***********************************************************************

  ;---Get array containing # dimensions, size of dimensions, etc. of vars_to_get
  sizearray = SIZE(vars_to_get)

  ;---Make sure vars_to_get is a 1-d array----
  if (sizearray[0] gt 1) then begin
     print,'Error in get_var_data: Inputted Array of variables to get'
     print,' is not a 1-d array'
     stop
  endif

  print,'Reading the following variables from the netCDF file: ',vars_to_get
  ;-----get grid_jlo, which tells us whether the netCDF file data
  ;     starts at j index 0 or 1 (there is no index in the netCDF file,
  ;     but the 0 or 1 refers to the starting index internal to EPIC,
  ;     and I want to stay consistent with that).
  NCDF_ATTGET,file_id,'grid_jlo',grid_jlo,/GLOBAL
  NCDF_ATTGET,file_id,'grid_ni',ni,/GLOBAL
  NCDF_ATTGET,file_id,'grid_nj',nj,/GLOBAL
  NCDF_ATTGET,file_id,'grid_nk',nk,/GLOBAL
  NCDF_ATTGET,file_id,'grid_dlt',grid_dlt,/GLOBAL
  NCDF_ATTGET,file_id,'grid_dln',grid_dln,/GLOBAL
  NCDF_ATTGET,file_id,'grid_globe_latbot',latbot,/GLOBAL ;these are +-90 for a whole
  NCDF_ATTGET,file_id,'grid_globe_lattop',lattop,/GLOBAL ; globe geometry

  ;---Get the correct sizes of the pads to use----------
  lonpad = 2
  if (grid_jlo eq 0) then begin       ; ==> channel or whole globe bdry condition
                                      ; in latitude.  Add pad only to north end of
				       ; domain but not south end (since south end
		        	       ; of domain already contains the south pole
				       ; or channel edge--while north end does not--
				       ; and the channel edge/pole is similar to a pad
				       ; in that its values are completely determined
				       ; by those of "interior" points).  But ONLY
				       ; add a pad to the north end if the variable
				       ; is q or v; otherwise, one isn't needed.
	 latpad = 1

  endif else begin     
     if (grid_jlo eq 1) then begin       ; doubly periodic north-south
        latpad = 2
     endif else begin                             
        print,'Error in get_var_info: grid_jlo did not equal 0 or 1, so'
        print,'the boundary condition at the north/south edges of domain is unknown.'
        print,'(This is needed to determine how to load data into the arrays.)'
        stop
     endelse
  endelse



  var = {dummy:0L}  ; FUDGE: just to create var before the loop.  Think of something
                    ; real to put here.

  ;%%%%%%%%%%%%%%%%%%%%%%% NOW LOOP OVER THE VARIABLES TO GET %%%%%%%%%%%%%%%%%%%
  for count=0, sizearray[1]-1 do begin

     if (!EPIC_VERSION eq !EPIC4) then begin

       if (!EXTRACT_MODE eq !NO) then begin
          NCDF_VARGET,file_id,vars_to_get[count],dummyvar  ; get var named vars_to_get[count]
                                                      ; and put it in dummyvar, which IDL
 	  		           		       ; defines to be the appropriate type
						       ; and length.
	endif else if (!EXTRACT_MODE eq !YES) then begin

	   NCDF_VARGET, file_id, vars_to_get[count], dummyvar, $
	      COUNT   = [ni, nj+1, nk,        1], $
	      STRIDE  = [1,   1,    1,        1], $
	      OFFSET  = [0,   0,    0, frame_id]
	endif else begin
	  print, 'get_var_data error: EXTRACT_MODE has to be either YES or NO'
	  stop
	end

        ;-----Get 1-d axis arrays.  Call them lon, lat, sigmatheta regardless of version:
	 ;     For variables like dvdt and dudt which don't have axes, make them.

	 x_axis_name = 'lon_'+vars_to_get[count] 
	 y_axis_name = 'lat_'+vars_to_get[count]
	 z_axis_name = 'sigmatheta_'+vars_to_get[count]

	 if (vars_to_get[count] eq 'dvdt') then begin
	    x_axis_name = 'lon_v'
	    y_axis_name = 'lat_v'
	    z_axis_name = 'sigmatheta_v'
	 endif
	 if (vars_to_get[count] eq 'dudt') then begin
	    x_axis_name = 'lon_u'
	    y_axis_name = 'lat_u'
	    z_axis_name = 'sigmatheta_u'
	 endif
	 if (vars_to_get[count] eq 'pv2') then begin
	    x_axis_name = 'lon_u'
	    y_axis_name = 'lat_v'
	    z_axis_name = 'sigmatheta_u'     
	 endif


;	 print,'x_axis_name=',x_axis_name
;	 print,'y_axis_name=',y_axis_name
;	 print,'z_axis_name=',z_axis_name

        NCDF_VARGET,file_id,x_axis_name,lon_dummyvar 
        NCDF_VARGET,file_id,y_axis_name,lat_dummyvar

        NCDF_VARGET,file_id, z_axis_name, sigmatheta_dummyvar

        ;---(I should see if I can check for attempts to retrieve nonexisting variables)
     endif else  begin
       print, 'Unrecognized EPIC version'
	stop
     endelse

     ;%%%%%%%%% For rest of this function, version number is irrelevant %%%%%%%%%%%%%%%


     ;---Now get size information about these arrays.
     sd = SIZE(dummyvar)
     slon_d = SIZE(lon_dummyvar)
     slat_d = SIZE(lat_dummyvar)
     ssigmatheta_d = SIZE(sigmatheta_dummyvar)


     ;----Error check that all the dimensions are consistent------
     if ((ni ne sd[1]) or (nk ne sd[3]) or (nj+2 ne sd[2]+latpad) or (ni ne slon_d[1]) $ 
         or (nj+2 ne slat_d[1]+latpad) or (nk ne ssigmatheta_d[1])) then begin
         print,'Error in get_var_info: Lengths of 1-d axis arrays and corresponding'
         print,'  length of that dimension in the 3-d (or 4-d) dependent variables'
         print,'  are inconsistent with each other--or with the values of grid_ni,'
         print,'  grid_nj, or grid_nk--in the netCDF file.'
         stop
     endif

     if ((sd[0] ne 3) and (sd[0] ne 4)) then begin
         print,'Error in get_var_data: Tried to read in a variable which does not'
         print,'have 3 or 4 dimensions.  Only dependent variables (and their'
         print,'tendencies) are allowed.'
         stop
     endif

     ;---Create new dummy variables with pads added.  I assume the vars are doubles-----
     ;   We use elements starting at 0 for lon, lat, and time, but 1 for sigmatheta.
     ;   longitude: ni+2 ==> 2 pads added since original array is ni long
     ;   latitude:  ni+2 ==> 1 or 2 pads added since orig. array is nj or nj+1 long
     ;              (depending on BC)
     ;   sigmatheta: ni+1 ==> I want the array defined from 1-->nk (which equals the
     ;              amount of data obtained from the netCDF file), but I don't know
     ;              how to NOT define the 0 element, so I'll define it from 0
     ;              which implies nk+1 elements.  0 element not used.
     if (sd[0] eq 3) then dummyvar_pad = dblarr(ni+2, nj+2,nk+1)

     if (sd[0] eq 4) then dummyvar_pad = dblarr(ni+2, nj+2,nk+1, sd[4])
     lon_dummyvar_pad = dblarr(ni+2)
     lat_dummyvar_pad = dblarr(nj+2)

     ;---Transfer data from orig. dummy vars to padded dummy vars----------------
     if (sd[0] eq 3) then dummyvar_pad [1:ni, grid_jlo:nj, 1:nk] = dummyvar
     if (sd[0] eq 4) then dummyvar_pad [1:ni, grid_jlo:nj, 1:nk, *] = dummyvar
     lon_dummyvar_pad[1:ni] = lon_dummyvar
     lat_dummyvar_pad[grid_jlo:nj] = lat_dummyvar

     ;---Add pad values.----------------------------------------  
     dummyvar_pad [0,*,*,*] = dummyvar_pad [ni,*,*,*]  ; works for 3 or 4-d dummyvar_pad
     dummyvar_pad [ni+1,*,*,*] = dummyvar_pad [1,*,*,*]
     lon_dummyvar_pad[0]=lon_dummyvar_pad[1] - grid_dln
     lon_dummyvar_pad[ni+1]=lon_dummyvar_pad[ni] + grid_dln
     if (grid_jlo eq 1) then begin   ;....doubly periodic in north-south direction
        dummyvar_pad [*, 0,*,*] = dummyvar_pad [*,nj,*,*]
        dummyvar_pad [*, nj+1,*,*] = dummyvar_pad [*,1,*,*]
        lat_dummyvar_pad[0] = lat_dummyvar_pad[1] - grid_dlt
        lat_dummyvar_pad[nj+1] = lat_dummyvar_pad[nj] + grid_dlt

     endif else begin                ;....channel or whole globe geometry

        ;.....First set the variable itself..............................
        if (vars_to_get[count] eq 'v') then dummyvar_pad [*,nj+1,*,*] = 0.0
                                             ; v is zero at the north pole

	 if (vars_to_get[count] eq 'q') then dummyvar_pad[*,nj+1,*,*] = -1.0 ;FUDGE
				; Need to do a circ. theorem on u[nj] to get real
			        ; value of q[nj+1]
        if ((vars_to_get[count] ne 'v') and (vars_to_get[count] ne 'q')) then begin
	    dummyvar_pad[*,nj+1,*,*] = dummyvar_pad[*,nj,*,*]
	 endif

	 ;.....Now set the axes for the variable...........................
        if (lattop eq 90.0) then                                       $
	    lat_dummyvar_pad[nj+1] = lat_dummyvar_pad[nj]+(1. + 0.5*sqrt(2.))*grid_dlt $
	    
	 else  $
           lat_dummyvar_pad[nj+1] = lat_dummyvar_pad[nj] + grid_dlt           
     endelse
     ; NOTE: The above requires lat/lon to INCREASE with increasing index

     ;---Use the size info to add the appropriate fields (tags) to structure var.
     var = CREATE_STRUCT(var, vars_to_get[count],dummyvar_pad,        $
          'lon_'+vars_to_get[count],lon_dummyvar_pad,                 $
          'lat_'+vars_to_get[count],lat_dummyvar_pad,                 $
          'sigmatheta_'+vars_to_get[count],sigmatheta_dummyvar)
     ; --what about time??? If it's a tendency, there's a time_index_var axis, but
     ; there's not if it's an extracted 4-d dependent variable.  FUDGE
  endfor

  return 

end   ;---end procedure get_var_data
;***********************************************************************
; get_gatt_data: Gets the values for the global attributes
;***********************************************************************
PRO get_gatt_data,id,basic
         ; id (input) is the ID # of the open netCDF file
	  ; basic (input) is the

  return   
end   ;---end procedure get_gatt_data


;*****************************************************************************
; get_av_winds: Assuming u and v exist in structure var, creates *averages*
;   u_av and v_av which are on the h-grid.  The u_av and v_av fields (i.e. tags)
;   are added to var.  The 1-d axis coordinates are also created and added 
;   to var.  Their values are taken from the averages of the u and v coords
;   (but should be the same as the lon/lat values of lon_p and lat_p).
;*****************************************************************************
PRO get_av_winds, id, var
   ; id (input)--ID # of open netCDF file.  Used to get some attribute values
   ; var (input and output)--structure containing fields u, lon_u, lat_u, sigmatheta_u,
   ;     as well as v, lon_v, lat_v, and sigmatheta_v.   Average values of all
   ;     these fields are created and added to structure.
;****************************************************************************

  u_av = var.u   ; create temporary arrays u_av and v_av
  v_av = var.v
  lon_u_av = var.lon_u
  lat_u_av = var.lat_u
  sigmatheta_u_av = var.sigmatheta_u

  lon_v_av = var.lon_v
  lat_v_av = var.lat_v
  sigmatheta_v_av = var.sigmatheta_v

  NCDF_ATTGET,id,'grid_ni',ni,/GLOBAL
  NCDF_ATTGET,id,'grid_nj',nj,/GLOBAL

  for i=1,ni do begin
     u_av[i,*,*,*] = 0.5*(var.u[i,*,*,*] + var.u[i+1,*,*,*])
     lon_u_av[i] = 0.5*(var.lon_u[i] + var.lon_u[i+1])
  endfor

  for j=0,nj do begin
     v_av[*,j,*,*] = 0.5*(var.v[*,j,*,*] + var.v[*,j+1,*,*])
     lat_v_av[j] = 0.5*(var.lat_v[j] + var.lat_v[j+1])
  endfor

  lon_v_av = var.lon_v
  lat_u_av = var.lat_u
  sigmatheta_u_av = var.sigmatheta_u
  sigmatheta_v_av = var.sigmatheta_v

  var = CREATE_STRUCT(var,'u_av',u_av, 'v_av', v_av,     $
        'lon_u_av', lon_u_av, 'lat_u_av', lat_u_av, 'sigmatheta_u_av', sigmatheta_u_av, $
        'lon_v_av', lon_v_av, 'lat_v_av', lat_v_av, 'sigmatheta_v_av', sigmatheta_v_av)


  return

end     ;---------procedure get_av_winds---------------


;*******************************************************************************
; get_layer_p
;*******************************************************************************
   ; id (input)--ID # of open netCDF file.  Used to get some attribute values
   ; var (input and output)--structure containing fields temp, CH_4, ...
PRO get_layer_p,id,var

 if (!EPIC_VERSION eq !EPIC4) then begin

  p_layer = var.p3   ; define layer_p and set its size

  NCDF_ATTGET,id,'grid_nk',nk,/GLOBAL
  NCDF_ATTGET,id,'planet_rgas',rgas,/GLOBAL
  NCDF_ATTGET,id,'planet_cp', cp,/GLOBAL

  kappa = rgas/cp

  ooopk = 1./(1.+ kappa) ; One Over One Plus Kappa
  ook  = 1./kappa        ; One Over Kappa
  opk = 1. + kappa       ; One Plus Kappa

  for k=1,nk do begin
     if (k eq 1) then begin
	 p_layer[*,*,k] = ooopk^ook * var.p3[*,*,k]
     endif else begin
	 p_layer[*,*,k] = ( ooopk * (var.p3[*,*,k]^opk - var.p3[*,*,k-1]^opk)/     $
			   (var.p3[*,*,k] - var.p3[*,*,k-1]))^ook
     endelse
  endfor 


   var = CREATE_STRUCT(var, 'p_layer',p_layer,		 $
			    'lon_p_layer',var.lon_p3,     $
			    'lat_p_layer',var.lat_p3,     $
			    'sigmatheta_p_layer', var.sigmatheta_p3)   

 endif else begin
  print, 'get_layer_p error: unrecognized EPIC_VERSION'
 endelse


end     ;------- get_layer_p--------------------------


;************************************************************************
; plot_flat_horiz_scalar: Plot a scalar quantity (e.g., u, v, vorticity, 
;    PV, divergence, p, or theta) on a "horizontal" surface  (i.e., a surface 
;    of constant sigmatheta), using a contour and/or greyscale plot.
;    Quantity is plotted with either contours, greyscale, or both (overlain
;    on top of each other).  Calling routine specifies the quantity to be
;    plotted, which type of plot to make (greyscale/contour), the number of
;    of contours (the value of which is ignored if only greyscale plotting
;    is to be done), which vertical layer and timestep to plot, and how
;    the max/min values of the contoured/greyscaled range is to be chosen.
;    (The contours are assumed to be evenly spaced within this range.)
;         
;    In theory the routine should not have to be altered when the calling
;    routine wants to plot new (previously nonexistent) quantities such as
;    f_p, water vapor amount, dust loading, etc.; it requires only that
;    the inputted var structure contain the relevant fields and that the
;    inputted name of the quantity to be plotted corresponds to one of those
;    fields.  1-d arrays for the axes, of the form lon_quantity, lat_quantity,
;    and sigmatheta_quantity must also exist.
;************************************************************************
PRO plot_flat_horiz_scalar,id, var,quantity_name,overplot, contour,greyscale, $
     numcontours,vert_level, timestep, range,postscript

    ; id (input) is netCDF file ID # so we can get some global attributes
    ; var (input) is a structure containing the arrays of all the data
    ; quantity_name (input) is a string giving the name of the quantity to be
    ;      plotted.  Must be the same as that of the appropriate field in
    ;      structure var.  E.g. 'u' or 'v'.
    ; overplot (input) is a flag (set to !NEW or !OVERLAY) which tells whether we
    ;      want to overlay the plot on a previously existing plot (e.g. a plot of
    ;      wind vectors or a greyscale temperature map).  That will only work if
    ;      we want a contour (but not greyscale) plot, since greyscale destroys
    ;      pre-existing plots.
    ; contour (input) is a yes-no variable saying whether we want to plot contours
    ; greyscale (input) is a yes-no variable saying whether we want to plot greyscale
    ; numcontours (input) gives the number of evenly spaced contours (btw. min/max values)
    ; vert_level (input) is an integer specifying which vertical level to plot
    ; timestep (input) is an integer specifying which timestep to plot
    ; range (input) is an integer flag specifying over what range the
    ;      contours/greyscale are to be divided.  (I.e. do they extend from the
    ;      min to max value on just that level, or for all times on just that
    ;      level, or on all levels at that time, or globally for all levels/times?) 
    ; postscript (input) flags whether output will be to postscript file (!YES)
    ;      or to screen (!NO)
  ;**************************************************************************
  ;--------Do error checking----------------------------------------------- 
  if ((contour ne !YES) and (greyscale ne !YES)) then begin
      print,'Error in plot_horiz_scalar: neither contours nor greyscale'
      print,' was requested; must have one or both to make a plot.'
      stop
  endif

  print,'plot_flat_horiz_scalar: vert_level=',vert_level
  ; I should make sure vert_level and timestep are within the correct bounds


  ;------Get global attributes---------
  NCDF_ATTGET,id,'grid_jlo',grid_jlo,/GLOBAL
  NCDF_ATTGET,id,'grid_ni',ni,/GLOBAL
  NCDF_ATTGET,id,'grid_nj',nj,/GLOBAL
  NCDF_ATTGET,id,'grid_nk',nk,/GLOBAL


  ;------extract the desired arrays from the structure---------------------
   result1 = EXECUTE('quantity = var.'+quantity_name) ; Execute the command in the string.
                                                     ; (Now quantity is an array
						      ; containing what we want to plot.)

   result2 = EXECUTE('lon_quantity = var.lon_'+quantity_name)
   result3 = EXECUTE('lat_quantity = var.lat_'+quantity_name)
   result4 = EXECUTE('sigmatheta_quantity = var.sigmatheta_'+quantity_name)
   if (result1 eq 0) then begin
       print,'Error in plot_horiz_scalar: You tried to plot a quantity which'
       print,'  does not exist in the var structure passed to plot_horiz_scalar.'
       stop
   endif

   if ((result2 eq 0) or (result3 eq 0) or (result4 eq 0)) then begin
      print,'Error in plot_horiz_scalar: You wanted to plot '+quantity_name+','
      print,'  but the arrays lon_'+quantity_name+', lat_'+quantity_name+', or'
      print,'  sigmatheta_'+quantity_name+' did not all exist in structure var.'
      stop
   endif
   ;------extract a 2-d array at a particular level and time---------
   print,'plot_flat_horiz_scalar: max(quantity)=',max(quantity)
   quantity_slice = quantity[*,*,vert_level]   ; assumes it's not a 4-d array

   max_quantity_slice = max(quantity_slice)
   min_quantity_slice = min(quantity_slice)
   print,'plot_flat_horiz_scalar: max(quantity_slice)=',max_quantity_slice 
   print,'plot_flat_horiz_scalar: min(quantity_slice)=',min_quantity_slice


   qmax=max(quantity_slice[1:ni,grid_jlo:nj],min=qmin)     ; get min and max value of quantity_slice
   ; the above is not quite right since for q,v it should go to nj+1 while for
   ; u or h it should go to nj.  But I should generalize to whether it mins/maxes
   ; on just this slice or all values etc.

   ;---preparation for contour and greyscale plots-----------------------
   if (quantity_name eq 'CH_4_rh') then begin
      quantity_slice = quantity_slice^3
   endif
   qpix = BYTSCL(quantity_slice,min=qmin,max=qmax)  ; rescale quantity_slice from
                              ; doubles of range (min,max) into
                              ; bytes in the range (0,top), 
                              ; (255 is default for top).

   levelq = qmin + (qmax-qmin)*(FINDGEN(numcontours))/double(numcontours-1)
               ; Make an array to give u value of contour levels, evenly
               ; spaced from umin to umax.  numcontours different levels.
   if (qmin eq qmax) then begin
      levelq = [0.]
   endif
   xbot = -180. & xtop = 180.
   ybot = -90. & ytop = 90. 

   ;--create the window size etc. with a call to contour; this does not plot 
   ;  any data (because of the /NODATA flag) and is necessary even if user 
   ;  only wants a greyscale plot.

   xbot = min(lon_quantity) & xtop = max(lon_quantity)
   ybot = min(lat_quantity) & ytop = max(lat_quantity) 
   if (greyscale eq !NO) then begin      ;------------------------contour plot ONLY
      if (overplot eq !NEW) then                $
         CONTOUR,quantity_slice,lon_quantity,lat_quantity,levels=levelq,$ 
             xstyle=1,ystyle=1,xrange=[xbot,xtop],yrange=[ybot,ytop],/FOLLOW ,$
	      xtitle='longitude [deg]' ,ytitle='latitude [deg]'
                                 ; plot the data w/o overplot keyword

      if (overplot eq !OVERLAY) then            $
         CONTOUR,quantity_slice,lon_quantity,lat_quantity,levels=levelq,$ 
               xstyle=1,ystyle=1,/OVERPLOT,xrange=[xbot,xtop],yrange=[ybot,ytop],$
		/FOLLOW
                                 ; plot the data w/ overplot keyword
   endif




   if (greyscale eq !YES) then begin      ;----greyscale plot, with or without contours
      if (overplot eq !NEW) then begin
         tit=''
         tit=quantity_name+'   max='+string(max_quantity_slice)+   $
	          '  min='+string(min_quantity_slice)
        if quantity_name eq 'p' then tit=$
	 quantity_name+'   max='+string(max_quantity_slice/1.d5)+   $
	          '  min='+string(min_quantity_slice/1.d5)


         CONTOUR,quantity_slice,lon_quantity,lat_quantity,   $
              levels=levelq,xstyle=1,ystyle=1,/NODATA,xrange=[xbot,xtop],  $
              yrange=[ybot,ytop],xtitle='longitude [deg]',ytitle='latitude [deg]', $
	       title=tit
      endif
      if (overplot eq !OVERLAY) then                $
         CONTOUR,quantity_slice,lon_quantity,lat_quantity,      $
               xrange=[xbot,xtop],yrange=[ybot,ytop],           $
               levels=levelq,xstyle=1,ystyle=1,/NODATA,/OVERPLOT

      print,'Plot w/ greyscale'
      if (postscript ne !YES) then begin
         px = !X.WINDOW*!D.X_VSIZE  ; Make two-element arrays PX and PY containing
         py = !Y.WINDOW*!D.Y_VSIZE  ; containing size of plot window in device 
				     ; pixels

         sx = px[1] - px[0] + 1   ; sx and sy are scalars which give size of image in
         sy = py[1] - py[0] + 1   ; pixels for the x and y directions.
         TVSCL, CONGRID(qpix,sx,sy),px[0],py[0]    ; Plot qpix w/ greyshading.

      endif else begin      ; we want to output to ps file
         TV, qpix,!x.window(0),!y.window(0), xsize=!x.window(1)-!x.window(0), $
	      ysize=!y.window(1)-!y.window(0),/norm
      endelse

      if (contour eq !YES) then begin
         CONTOUR,quantity_slice,lon_quantity,lat_quantity, levels=levelq,xstyle=1,$
             ystyle=1,/NOERASE,/OVERPLOT,xrange=[xbot,xtop],yrange=[ybot,ytop], $
	      /FOLLOW,xtitle='longitude [deg]',ytitle='latitude [deg]',$
	      title=''
      endif

   endif


  print,'min and max values:',qmin,qmax

end    ; procedure plot_flat_horiz_scalar


;************************************************************************
; plot_surface_horiz_scalar: Plot a scalar quantity (e.g., u, v, vorticity,
;    PV, divergence, p, or theta) on a "horizontal" surface (i.e., a surface 
;    of constant sigmatheta) using a 3-d wire-mesh plot.
;    Quantity is plotted with either contours, greyscale, or both (overlain
;    on top of each other).  Calling routine specifies the quantity to be
;    plotted, which type of plot to make (greyscale/contour), the number of
;    of contours (the value of which is ignored if only greyscale plotting
;    is to be done), which vertical layer and timestep to plot, and how
;    the max/min values of the contoured/greyscaled range is to be chosen.
;    (The contours are assumed to be evenly spaced within this range.)
;         
;    In theory the routine should not have to be altered when the calling
;    routine wants to plot new (previously nonexistent) quantities such as
;    f_p, water vapor amount, dust loading, etc.; it requires only that
;    the inputted var structure contain the relevant fields and that the
;    inputted name of the quantity to be plotted corresponds to one of those
;    fields.
;************************************************************************
PRO plot_surface_horiz_scalar,var,quantity_name, vert_level, timestep, range

    ; var (input) is a structure containing the arrays of all the data
    ; quantity_name (input) is a string giving the name of the quantity to be
    ;      plotted.  Must be the same as that of the appropriate field in
    ;      structure var.  E.g. 'u' or 'v'.
    ; vert_level (input) is an integer specifying which vertical level to plot
    ; timestep (input) is an integer specifying which timestep to plot
    ; range (input) is an integer flag specifying over what range the
    ;      contours/greyscale are to be divided.  (I.e. do they extend from the
    ;      min to max value on just that level, or for all times on just that
    ;      level, or on all levels at that time, or globally for all levels/times?) 
  ;**************************************************************************
  ;-------Error checking------------------
  ; I should make sure vert_level and timestep are within the correct bounds

  print,'plot_surface_horiz scalar: vert_level=',vert_level
  ;------extract the desired arrays from the structure---------------------
   result1 = EXECUTE('quantity = var.'+quantity_name) ; Execute the command in the string.
                                                     ; (Now quantity is an array
						      ; containing what we want to plot.)

   result2 = EXECUTE('lon_quantity = var.lon_'+quantity_name)
   result3 = EXECUTE('lat_quantity = var.lat_'+quantity_name)
   result4 = EXECUTE('sigmatheta_quantity = var.sigmatheta_'+quantity_name)
   if (result1 eq 0) then begin
       print,'Error in plot_horiz_scalar: You tried to plot a quantity which'
       print,'  does not exist in the var structure passed to plot_horiz_scalar.'
       stop
   endif

   if ((result2 eq 0) or (result3 eq 0) or (result4 eq 0)) then begin
      print,'Error in plot_horiz_scalar: You wanted to plot '+quantity_name+','
      print,'  but the arrays lon_'+quantity_name+', lat_'+quantity_name+', or'
      print,'  sigmatheta_'+quantity_name+' did not all exist in structure var.'
      stop
   endif
   ;------extract a 2-d array at a particular level and time---------
   print,'plot_surface_horiz_scalar: max(quantity)=',max(quantity)  ;FUDGE

   quantity_slice = quantity[*,*,vert_level]   ; assumes it's not a 4-d array FUDGE

   print,'plot_surface_horiz_scalar: max(quantity_slice)=',max(quantity_slice)  ;FUDGE
   print,'plot_surface_horiz_scalar: min(quantity_slice)=',min(quantity_slice)  ;FUDGE

   qmax=max(quantity_slice,min=qmin)     ; get min and max value of quantity_slice


  ;----procedure SURFACE plots wiremesh lines for every array element spacing
   ;    which may be too dense.  Rescale the array and get new axes
   qplot = CONGRID(quantity_slice,128,64)
   lon_qplot = CONGRID(lon_quantity,128)
   lat_qplot = CONGRID(lat_quantity,64)

   SURFACE,qplot,lon_qplot,lat_qplot,ax=30,az=30   


;    print,'min and max values:',qmin,qmax
;    stop   ; FUDGE
end    ; procedure plot_surface_horiz_scalar

;************************************************************************
; plot_flat_horiz_vector: Plot a vector quantity (e.g., horizontal velocity) 
;    on a "horizontal" surface  (i.e., a surface of constant sigmatheta)
;     Calling routine specifies the quantity to be
;    plotted, which type of plot to make (greyscale/contour), the number of
;    of contours (the value of which is ignored if only greyscale plotting
;    is to be done), which vertical layer and timestep to plot, and how
;    the max/min values of the contoured/greyscaled range is to be chosen.
;    (The contours are assumed to be evenly spaced within this range.)
;         
;    In theory the routine should not have to be altered when the calling
;    routine wants to plot new (previously nonexistent) quantities such as
;    f_p, water vapor amount, dust loading, etc.; it requires only that
;    the inputted var structure contain the relevant fields and that the
;    inputted name of the quantity to be plotted corresponds to one of those
;    fields.  1-d arrays for the axes, of the form lon_quantity, lat_quantity,
;    and sigmatheta_quantity must also exist.
;************************************************************************
PRO plot_flat_horiz_vector,var,quantity_namex,quantity_namey, $
         vert_level,range,overlay,xshift,yshift

    ; var (input) is a structure containing the arrays of all the data
    ; quantity_namex (input) is a string giving the name of the scalar 
    ;      quantity of the x vector component to be plotted (e.g., u).
    ;     Must be the same as that of the appropriate field in
    ;      structure var.
    ;  quantity_namey (input). Same but for y component (e.g., v).
    ; vert_level (input) is an integer specifying which vertical level to plot
    ; range (input) is an integer flag specifying over what range the
    ;      contours/greyscale are to be divided.  (I.e. do they extend from the
    ;      min to max value on just that level, or for all times on just that
    ;      level, or on all levels at that time, or globally for all 
    ;	    levels/times?) 
    ; overlay (input): integer flag that tells whether to create a new plot
    ;      or overlay on an existing plot
    ; xshift, yshift (input): doubles giving amount to add to x and y
    ;      quantities before plotting them (e.g. to allow change in ref. frame
    ;	    for wind).
  ;**************************************************************************
  ;--------Do error checking----------------------------------------------- 

       ;-------Error checking------------------
  ; I should make sure vert_level and timestep are within the correct bounds

  ;------extract the desired arrays from the structure---------------------
   result1 = EXECUTE('quantityx = var.'+quantity_namex)
				 ; Execute the command in the string.
                                ; (Now quantity is an array
				 ; containing what we want to plot.)

   result2 = EXECUTE('lon_quantityx = var.lon_'+quantity_namey)
   result3 = EXECUTE('lat_quantityx = var.lat_'+quantity_namey)
   result4 = EXECUTE('sigmatheta_quantityx = var.sigmatheta_'+quantity_namey)
   result5 = EXECUTE('quantityy = var.'+quantity_namey)
				 ; Execute the command in the string.
                                ; (Now quantity is an array
				 ; containing what we want to plot.)

   result6 = EXECUTE('lon_quantityy = var.lon_'+quantity_namey)
   result7 = EXECUTE('lat_quantityy = var.lat_'+quantity_namey)
   result8 = EXECUTE('sigmatheta_quantityy = var.sigmatheta_'+quantity_namey)

   ;----make sure the quantities to plot, and their axes, exist-------
   if ((result1 eq 0) or (result5 eq 0)) then begin
       print,'Error in plot_flat_horiz_vector: You tried to plot a quantity which'
       print,'  does not exist in the var structure passed in as input.'
       stop
   endif

   if ((result2 eq 0) or (result3 eq 0) or (result4 eq 0) or (result6 eq 0) or     $
       (result7 eq 0) or (result8 eq 0)) then begin
      print,'Error in plot_flat_horiz_vector: Input structure var does not contain'
      print,' the 1-d arrays corresponding to the axes of what you want to plot.' 
      stop
   endif

   ;---I should error check to make sure the lon and lat axes for the two
   ;    are the same.  Skip for now.  Would be better to do when the x and y
   ;    plot axes aren't hardwired to be lat and lon.

   ;------extract a 2-d array at a particular level and time---------
   quantity_slicex = quantityx[*,*,vert_level] 
   quantity_slicey = quantityy[*,*,vert_level]   


   ;-----get max/min values-----------------
   qmaxx=max(quantity_slicex,min=qminx)     ; get min and max value of quantity_slice
   qmaxy=max(quantity_slicey,min=qminy)     ; i.e. of east and north winds
   magnitude = sqrt(quantity_slicex^2 + quantity_slicey^2)
   qmax = max(magnitude, min=qmin)   ;  get min/max of TOTAL wind amplitude

  ;----procedure VEL plots wiremesh lines for every array element spacing
   ;    which may be too dense.  Rescale the array and get new axes

   rebinx = 64
   rebiny = 32
   qplotx = congrid(quantity_slicex,rebinx,rebiny)
   qploty = congrid(quantity_slicey,rebinx,rebiny)
   lon_qplot = CONGRID(lon_quantityx,rebinx)
   lat_qplot = CONGRID(lat_quantityx,rebiny)

;     VELOVECT,quantity_slicex,quantity_slicey,lon_quantityx,lat_quantityx,$
;         ytitle='Latitude [deg]',xtitle='Longitude [deg]' 

  if (overlay eq !NEW) then begin
     veltitle = 'Level '+string(vert_level)+' Velocity    max='+string(qmax)+  $
	       ' m/sec   min='+string(qmin) + ' m/sec'
     VELOVECT,qplotx + xshift,qploty + yshift,lon_qplot,lat_qplot,$
              ytitle='Latitude [deg]',xtitle='Longitude [deg]',length=2.0, $
	       thick=5.0,title=veltitle
  endif else begin
     VELOVECT,qplotx + xshift,qploty + yshift,lon_qplot,lat_qplot,$
             ytitle='Latitude [deg]',xtitle='Longitude [deg]',/OVERPLOT, $
	      length=2.0,thick=5.0
  endelse

   print,'min and max values:',qmin,qmax


end    ; procedure plot_flat_horiz_vector


;************************************************************************
; plot_1dslice_scalar: Plot a one-dimensional slice of a scalar quantity 
;    (e.g., u, v, vorticity, PV, divergence, p, or theta) as a function
;    of latitude, longitude, or sigmatheta, with the other two held
;    constant at particular values to be specified.   Calling routine 
;    specifies the quantity to be plotted and which axis to plot against
;    as input strings, and the fixed values of I, J, and K at which
;    the plot is done are provided as input integers.  (Only two of
;    these make sense for the plot; the third -- whichever one we're
;    plotting against -- should have the value -1.)
; 
;    In theory the routine should not have to be altered when the calling
;    routine wants to plot new (previously nonexistent) quantities such as
;    f_p, water vapor amount, dust loading, etc.; it requires only that
;    the inputted var structure contain the relevant fields and that the
;    inputted name of the quantity to be plotted corresponds to one of those
;    fields.  1-d arrays for the axes, of the form lon_quantity, lat_quantity,
;    and sigmatheta_quantity must also exist.
;************************************************************************
PRO plot_1dslice_scalar,id, var,quantity_name,plot_axis, I_plot, J_plot, $
     K_plot,overlay

    ; id (input) is netCDF file ID # so we can get some global attributes
    ; var (input) is a structure containing the arrays of all the data
    ; quantity_name (input) is a string giving the name of the quantity to be
    ;      plotted.  Must be the same as that of the appropriate field in
    ;      structure var.  E.g. 'u' or 'v'.
    ; plot_axis (input) is a string containing 'lon', 'lat', or 'sigmatheta,'
    ;      which tells us which axis to plot against.  If plot_axis = 'p'
    ;	    then we do a vertical plot (at a single x,y) using p as the vert.coord.
    ; I_plot, J_plot, K_plot (input): For the two of these three variables
    ;      to be held fixed, their values give the values to hold constant.
    ;      The third quantity (corresponding to plot_axis) should be -1.
    ; overlay (input) is a flag (set to !NEW or !OVERLAY) which tells whether we
    ;      want to overlay the plot on a previously existing plot (e.g. a plot of
    ;      wind vectors or a greyscale temperature map).  That will only work if
    ;      we want a contour (but not greyscale) plot, since greyscale destroys
    ;      pre-existing plots.
  ;**************************************************************************
  ;--------Do error checking----------------------------------------------- 
  if ((plot_axis ne 'lon') and (plot_axis ne 'lat') and (plot_axis ne       $
       'sigmatheta') and (plot_axis ne 'p')) then begin
      print,'Error in plot_1dslice_scalar: Axis to plot not clearly specified'
      print,'     as either lon, lat, sigmatheta, or p'
      stop
  endif

  ;------Get global attributes---------
  NCDF_ATTGET,id,'grid_jlo',grid_jlo,/GLOBAL
  NCDF_ATTGET,id,'grid_ni',ni,/GLOBAL
  NCDF_ATTGET,id,'grid_nj',nj,/GLOBAL
  NCDF_ATTGET,id,'grid_nk',nk,/GLOBAL


  ;------extract the desired arrays from the structure---------------------
   result1 = EXECUTE('quantity = var.'+quantity_name) ; Execute the command in the string.
                                                     ; (Now quantity is an array
						      ; containing what we want to plot.)
   if ((plot_axis eq 'lon') or (plot_axis eq 'lat') or        $
                   (plot_axis eq 'sigmatheta')) then begin
      result2 = EXECUTE('axis_quantity = var.'+plot_axis+'_'+quantity_name)
      if (result2 eq 0) then begin
         print,'Error in plot_1dslice_scalar: You wanted to plot '+quantity_name+','
         print,'  but the arrays '+plot_axis+'_'+quantity_name+' did not exist'
         print,'  in structure var.'
         stop
      endif
   endif

   if (plot_axis eq 'p') then begin
      axis_quantity = var.p_layer[I_plot, J_plot, *]
      axis_quantity[nk] = 2.* axis_quantity[nk-1] - axis_quantity[nk-2]
      ;   (p in the abyssal layer is zero, so just make a fake p which is
      ;    higher than the last active layer's p by a typical layer spacing of p.)
   endif


   if (result1 eq 0) then begin
       print,'Error in plot_1dslice_scalar: You tried to plot a quantity which'
       print,'  does not exist in the var structure.'
       stop
   endif


   ;------extract the appropriate 1-d array ---------
   print,'plot_1dslice_scalar: max(quantity)=',max(quantity)
   if (plot_axis eq 'lon') then begin
      quantity_slice = quantity[*,J_plot, K_plot]
      nlo = 1
      nhi = ni
   endif

   if (plot_axis eq 'lat') then begin
      quantity_slice = quantity[I_plot, *, K_plot]
      nlo = grid_jlo
      nhi = nj
   endif

   if ( (plot_axis eq 'sigmatheta') or (plot_axis eq 'p') ) then begin
      quantity_slice = quantity[I_plot, J_plot, *]
      nlo = 1
      nhi = nk-1
   endif

   print,'plot_1dslice_scalar: max(quantity_slice)=',max(quantity_slice)  
   print,'plot_flat_horiz_scalar: min(quantity_slice)=',min(quantity_slice)


   qmax=max(quantity_slice[nlo:nhi],min=qmin)   ; get min and max value of quantity_slice
   ; the above is not quite right since for q,v it should go to nj+1 while for
   ; u or h it should go to nj.  But I should generalize to whether it mins/maxes
   ; on just this slice or all values etc.

   ;---Set the plot ranges for the indep variable----------------
   axis_bot = min(axis_quantity[nlo:nhi])
   axis_top = max(axis_quantity[nlo:nhi])

   ;---If we plot against pressure, swap bot and top values so
   ;   so p increases down:
   if (plot_axis eq 'p') then begin
      tmp = axis_bot
      axis_bot = axis_top
      axis_top = tmp
   endif

   ;---Make the plot---------------------
   if (plot_axis eq 'lon') or (plot_axis eq 'lat')  then begin
      if (overlay eq !NEW) then begin
         plot,axis_quantity,quantity_slice,xrange=[axis_bot,axis_top], xstyle=1,$
	        xtitle=plot_axis+'  [deg]', ytitle=quantity_name
      endif else begin
         oplot, axis_quantity, quantity_slice
      endelse
   endif

   if (plot_axis eq 'sigmatheta' or plot_axis eq 'p') then begin
      if (overlay eq !NEW) then begin
         plot,quantity_slice, axis_quantity, yrange=[axis_bot,axis_top], ystyle=1,$
	        ytitle=plot_axis, xtitle=quantity_name
      endif else begin
         oplot, quantity_slice, axis_quantity
      endelse
   endif

  print,'min and max values:',qmin,qmax

end    ; procedure plot_1dslice_scalar







