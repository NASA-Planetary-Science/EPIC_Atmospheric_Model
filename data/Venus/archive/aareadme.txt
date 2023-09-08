CCSD3ZF0000100000001NJPL3IF0PDS200000001 = SFDU_LABEL                         
RECORD_TYPE                     = FIXED_LENGTH                                
RECORD_BYTES                    = 80                                          
SPACECRAFT_NAME                 = MAGELLAN                                    
TARGET_NAME                     = VENUS                                       
OBJECT                          = TEXT                                        
  PUBLICATION_DATE              = 1993-04-01                                  
  BYTES                         = 80                                          
  NOTE                          = "Magellan GxDR CD-ROM"                      
END_OBJECT                      = TEXT                                        
END                                                                           
                                                                              
                            MAGELLAN GxDR CD-ROM                              
                                                                              
                                                                              
1.  Introduction                                                              
                                                                              
    Overall, the collection of Magellan CD-ROMs contains three data           
sets.  One is the GxDR data set, consisting of MIT Altimetric and             
Radiometric Global Data Records: the Global Topography Data Record            
(GTDR), Global Emissivity Data Record (GEDR), Global Slope Data Record        
(GSDR), and Global Reflectivity Data Record (GREDR).  Another data set        
contains the Altimetry and Radiometry Composite Data Record (ARCDR)           
from MIT.  The third data set is the MIDR (Mosaic Image Data Record)          
a collection of high-resolution radar backscatter images from JPL.            
                                                                              
    This CD-ROM contains MIT Altimetric and Radiometric Global Data           
Record and associated documentation.  The complete data set includes          
GTDR, GREDR, GSDR, and GEDR, i.e. global images of Venus in which the         
data values represent topography, radar reflectivity, meter-scale             
surface slope, and microwave emissivity. For each data type, the planet       
is presented in four map projections: sinusoidal, Mercator, and north         
and south polar stereographic. In addition, the sinusoidal topography         
image (GTDR) is accompanied by a "radius error" image.                        
                                                                              
    The volume also contains documentation and index files to support         
access and use of the images. The Mercator and sinusoidal images cover        
the entire planet. The sinusoidal images are centered on each pole and        
extend to +-33 degrees latitude in the corners. The mapping scale is          
the same for each image--each sinusoidal pixel subtends 360/8192 deg.,        
or about 4.641 kilometers, and this is also true for equatorial pixels        
in the Mercator images and for polar pixels in the stereographic              
images. This is explained more fully in the DSMAP.LBL file in the             
[LABEL] directory.                                                            
                                                                              
    All data formats are based on the Planetary Data System Data              
Preparation Workbook ("PDS Data Preparation Workbook, D-7669. May             
1991.) and are similar to the formats used in generating the Voyager          
CD-ROM set ("Voyagers to the Outer Planets", volumes 1-8, Planetary           
Data System, JPL, 1989).                                                      
                                                                              
2.  Disk Format                                                               
                                                                              
    The disk has been formatted so that a variety of computer systems         
(e.g. IBM PC, Macintosh, Sun, VAX) may access the data.  Specifically,        
the disk is formatted according to the ISO 9660 level 1 Interchange           
Standard, and file attributes are specified by Extended Attribute             
Records (XARs).  For computer software that fully supports XARs, access       
to the CD-ROM volume will be straightforward; the disk will appear to         
the user to be identical to a file system of directories, sub-                
directories, and data files.  Some computer systems that do not support       
XARs will ignore them; others will append the XAR to the beginning of         
the file.  In the latter case the user must ignore the first 512 bytes        
of the file.  For further information, refer to the ISO 9660 Standard         
Document: RF# ISO 9660-1988, 15 April 1988.                                   
                                                                              
3.  File Formats                                                              
                                                                              
    Each GxDR is divided into two arrays of 32 framelets (in sinusoidal       
and Mercator map projections) arranged in four rows and eight columns,        
and two arrays of 4 framelets (in north and south polar stereographic         
projections) arranged in two rows and two columns.  The GTDR contains         
a fifth array of 32 framelets holding the "radius error image" (in a          
sinusoidal map projection).                                                   
                                                                              
    The framelets are numbered in increasing order from left to right,        
top to bottom.  Each framelet is stored in a separate file.  A framelet       
is 1024 lines by 1024 samples. The GREDR, GSDR, and the GTDR radius           
error image use one byte per sample, whereas the GEDR and the remaining       
GTDR images, with their greater dynamic range, use two bytes per pixel.       
The framelet files contain embedded VICAR2 labels and have detached PDS       
labels.  The framelet files and supplementary data files for each GxDR        
are stored in a separate subdirectory.                                        
                                                                              
    Additional general information about each product is found in the         
[LABEL] directory.  The *.LBL files in that directory contain the PDS         
catalog information for the Magellan Mission, the spacecraft, the radar       
instrument, and the data sets.  Further information on GxDR products          
can also be obtained from the Magellan Software Interface Specification       
(SIS) document "Global Altimetry and Radiometry Data Record", MIT-MGN-        
GxDR, version 2.1, October 21, 1991. A copy of this document is located       
in the [DOCUMENT] directory of this disk, in the file named GXDR.TXT.         
                                                                              
    Browse versions of the GxDRs are provided for quick inspection of         
the images.  The browse version of each GxDR is found in the sub-             
directory for that GxDR.  The file name for the browse image is               
BROWSE.IMG.  Each browse image of a 32-framelet image is 512 lines by         
1024 samples, and was created by averaging groups of 8x8 pixels in the        
original 4096x8192 GxDR frame.  Each browse image of a 4-framelet image       
is 1024 lines by 1024 samples and was created by averaging groups of          
2x2 pixels in the original 2048x2048 GxDR frame. Browse images contain        
embedded VICAR2 labels and have detached PDS labels.                          
                                                                              
    Files in the [SOFTWARE] directory are intended to be used in              
specific environments -- MS-DOS, MacOS, and UNIX/X11. They are not            
necessarily useful, or even readable, within other operating systems.         
                                                                              
    All document files and detached label files contain 80-byte fixed-        
length records, with a carriage return character (ASCII 13) in the 79th       
byte and a line feed character (ASCII 10) in the 80th byte.  This             
allows the files to be read by the MacOS, DOS, Unix, and VMS operating        
systems.  All tabular files are also described by PDS labels, either          
embedded at the beginning of the file or detached.  If detached, the          
PDS label file has the same name as the data file it describes, with          
the extension .LBL; for example, the file GCONTENT.TAB is accompanied         
by the detached label file GCONTENT.LBL in the same directory.  The           
detached labels for MIDRs and GxDRs also contain PDS-defined map              
projection keywords that provide information needed to extract latitude       
and longitude values from given line and sample locations.                    
                                                                              
    Tabular files are formatted so that they may be read directly into        
many database management systems on various computers.  All fields are        
separated by commas, and character fields are enclosed in double              
quotation marks (").  Character fields are left justified, and numeric        
fields are right justified.  The "start byte" and "bytes" values listed       
in the labels do not include the commas between fields or the quotation       
marks surrounding character fields.  The records are of fixed length,         
and the last two bytes of each record contain the ASCII carriage return       
and line feed characters.  This allows a table to be treated as a fixed       
length record file on computers that support this file type and as a          
normal text file on other computers.                                          
                                                                              
    PDS labels are object-oriented.  The object to which the label            
refers (e.g. IMAGE, TABLE, etc.) is denoted by a statement of the form        
                                                                              
    ^object = location                                                        
                                                                              
in which the carat character (^, also called a pointer in this context)       
indicates that the object starts at the given location.  In an embedded       
label, the location is an integer representing the starting record            
number of the object (the first record in the file is record 1).  In a        
detached label, the location denotes the name of the file containing          
the object, along with the starting record or byte number if there is         
more than one object.  For example,                                           
                                                                              
    ^IMAGE_HEADER = ("F02.IMG",1)                                             
    ^IMAGE = ("F02.IMG",3)                                                    
                                                                              
indicates that the IMAGE object begins at record 3 of the file                
F02.IMG, in the same directory as the detached label file.  Below             
is a list of the possible formats for the ^object definition.                 
                                                                              
    ^object = n                                                               
    ^object = n<BYTES>                                                        
    ^object = ("filename.ext")                                                
    ^object = ("filename.ext",n)                                              
    ^object = ("[dirlist]filename.ext",n)                                     
    ^object = ("filename.ext",n<BYTES>)                                       
    ^object = ("[dirlist]filename.ext",n<BYTES>)                              
                                                                              
where                                                                         
                                                                              
    n        is the starting record or byte number of the object,             
             counting from the beginning of the file (record 1, byte 1).      
                                                                              
    <BYTES>  indicates that the number given is in units of bytes.            
                                                                              
    filename is the upper-case file name, before the decimal point.           
                                                                              
    ext      is the upper-case file name, after the decimal point.            
                                                                              
    dirlist  is a period-delimited path-list of parent directories, in        
             upper case, that specifies the object file directory.  The       
             list begins on directory level below the root directory of       
             the CD-ROM.  '[dirlist]' may be omitted when the object          
             being described is located either in the same directory as       
             the detached label or in a sub-directory named 'LABEL',          
             located one level below the CD-ROM root directory.               
                                                                              
4.  CD-ROM Contents                                                           
                                                                              
    The files on this CD-ROM are organized in one top-level directory         
with several subdirectories.  The following table shows the structure         
and content of these directories.  In the table, directory names are          
enclosed in square brackets ([]), upper-case letters indicate an actual       
directory or file name, and lower-case letters indicate the general           
form of a set of directory or file names.                                     
                                                                              
    FILE                 CONTENTS                                             
                                                                              
|- AAREADME.TXT       This text file.                                         
|                                                                             
|- GCUMCOMM.TXT       A description of known errors and anomalies in          
|                     this and ALL previous Magellan GxDR CD-ROM              
|                     disks.                                                  
|                                                                             
|- VOLDESC.SFD        A description of the contents of this CD-ROM            
|                     volume in a format readable by both human and           
|                     computers.                                              
|                                                                             
|- [INDEX]            A directory containing index files which describe       
|    |                products published so far or on this volume.            
|    |                                                                        
|    |- GCONTENT.LBL  A PDS label describing GCONTENT.TAB.                    
|    |                                                                        
|    |- GCONTENT.TAB  A table of GxDR frames in terms of their directory      
|    |                name and extent in latitude and longitude.              
|    |                                                                        
|    |- GCUMDIR.TAB   The contents of all previous GxDR products              
|    |                published so far, including the GxDR on this            
|    |                volume.                                                 
|    |                                                                        
|    |- GCUMDIR.LBL   A PDS label describing GCUMDIR.TAB.                     
|    |                                                                        
|    |- GEO.TAB       A table of IAU-defined Venus geologic features.         
|    |                                                                        
|    |- GEO.LBL       A PDS label describing GEO.TAB.                         
|                                                                             
|                                                                             
|- [LABEL]            A directory containing catalog information              
|     |               describing the GxDR magellan data products that         
|     |               will be submitted to PDS.  The information is in        
|     |               the form of catalog templates that will be              
|     |               submitted to PDS for loading into the PDS central       
|     |               node catalog.                                           
|     |                                                                       
|     |- CATALOG.LBL  PDS high-level experiment-description templates.        
|     |                                                                       
|     |- DSMAP.LBL    PDS high-level data set templates describing            
|     |               cartographic projections and references to them.        
|     |                                                                       
|     |- GEDRDS.LBL   PDS high-level dataset templates for GEDRs.             
|     |                                                                       
|     |- GREDRDS.LBL  PDS high-level dataset templates for GREDRs.            
|     |                                                                       
|     |- GSDRDS.LBL   PDS high-level dataset templates for GSDRs.             
|     |                                                                       
|     |- GTDRDS.LBL   PDS high-level dataset templates for GTDRs.             
|                                                                             
|- [DOCUMENT]         On production discs, this directory will contain        
|     |               document files (extension name "TXT") describing        
|     |               products, missions, organization, etc.                  
|     |                                                                       
|     |- GXDR.TXT     A machine-readable version of the MIT-MGN-GxDR          
|                     SIS document describing the format and contents         
|                     of the original data files supplied from MIT to         
|                     the Magellan project on magnetic tape.                  
|                                                                             
|- [SOFTWARE]         This directory contains useful software for             
|   |                 displaying the images on this and other Magellan        
|   |                 image disks.                                            
|   |                                                                         
|   |- [PC]           This directory contains software for MS-DOS             
|   |   |                                                                     
|   |   |- IMDISP.EXE   An executable MS-DOS binary command file for          
|   |   |               enhancing and displaying images on an IBM-PC          
|   |   |               computer (or compatible) with color display.          
|   |   |                                                                     
|   |   |- IMDISP.TXT   An ASCII file describing the MS-DOS IMDISP command.   
|   |   |                                                                     
|   |   |- SOFTWARE.TXT A short description of the IMDISP program.            
|   |                                                                         
|   |- [MAC]            This directory contains software for MAC-OS           
|   |   |                                                                     
|   |   |- IMAGE.HQX    A compressed Apple Macintosh archive file that        
|   |   |               contains the IMAGE application and its                
|   |   |               documentation. Use the public-domain "BINHEX" and     
|   |   |               "Stuffit" Macintosh program to unpack the archive.    
|   |   |                                                                     
|   |   |- SOFTWARE.TXT A short description of the IMAGE program.             
|   |                                                                         
|   |- [UNIX]                                                                 
|   |   |                                                                     
|   |   |- SOFTWARE.TXT A short description of the XMGNDISP program.          
|   |   |                                                                     
|   |   |- XMGNDISP.C   C-language source code for the xmgndisp application,  
|   |   |               a UNIX command that displays Magellan images in       
|   |   |               version 11 of the X window system. It uses only the   
|   |   |               standard Xlib routines and has been tested in the     
|   |   |               following environments:                               
|   |   |                 Hardware: MicroVAX, DecStation, Sun-SPARC           
|   |   |                 OS:       Ultrix 4.1, SunOS 4.0.3, SunOS 4.1.1      
|   |   |                 GUI:      X11R3, X11R4, OSF/Motif, OpenWindows 2    
|   |   |                 Window managers: olwm, mwm, twm                     
|   |   |                                                                     
|   |   |- XMGNDISP.MAN  UNIX manual file that documents MGNDISP.             
|   |   |                                                                     
|   |   |- XMGNDISP.TXT  ASCII version of XMGNDISP.MAN.                       
|                                                                             
|- [xxxxx]            Four sub-directories corresponding to the four          
|   |                 GxDR data types:                                        
|   |                                                                         
|   |                   xxxxx = GEDR  (global emissivity data record)         
|   |                   xxxxx = GREDR (global reflectivity data record)       
|   |                   xxxxx = GSDR  (global slope data record)              
|   |                   xxxxx = GTDR  (global topographic data record)        
|   |                                                                         
|   |- [yyy]          four sub-directories (5 for GTDR) containing data       
|   |   |             in various map projections:                             
|   |   |                                                                     
|   |   |               yyy = MERC   Mercator projection                      
|   |   |               yyy = NORTH  North Polar Stereographic projection     
|   |   |               yyy = SINUS  Sinusoidal projection                    
|   |   |               yyy = SOUTH  South Polar Stereographic projection     
|   |   |               yyy = ERROR  Radius error data in sinusoidal          
|   |   |                            projection (for GTDR only).              
|   |   |                                                                     
|   |   |- BROWSE.IMG  A 1/8 scaled (for SINUS, MERC, and ERROR) or           
|   |   |              a 1/2 scaled (for NORTH and SOUTH)                     
|   |   |              version of the original image.                         
|   |   |                                                                     
|   |   |- BROWSE.LBL  A PDS label describing BROWSE.IMG.                     
|   |   |                                                                     
|   |   |- PALETTE.LBL A PDS label describing PALETTE.LBL.                    
|   |   |                                                                     
|   |   |- PALETTE.TAB A table of RGB color values to assist in viewing       
|   |   |              or printing the image in BROWSE.IMG. For 'BYTE'        
|   |   |              data (i.e. GREDR, and GSDR, and the radius error       
|   |   |              GTDR image, this table will also be appropriate        
|   |   |              for the image framelets themselves.)                   
|   |   |                                                                     
|   |   |- FRAME.LBL   A PDS label describing FRAME.TAB.                      
|   |   |                                                                     
|   |   |- FRAME.TAB   A table describing the range of latitude and           
|   |   |              longitude within each image framelet.                  
|   |   |                                                                     
|   |   |- Fff.IMG     A 1024x1024-pixel image framelet file, where           
|   |   |              ff is the framelet number.                             
|   |   |                                                                     
|   |   |- Fff.LBL     PDS labels describing the corresponding Fff.IMG.       
|   |   |                                                                     
|   |   |- HIST.LBL    A PDS label describing HIST.TAB.                       
|   |   |                                                                     
|   |   |- HIST.TAB    A binary histogram of pixel values in the complete     
|   |   |              image frame.                                           
                                                                              
                                                                              
5.  Recommended CD-ROM Drives and Driver Software                             
                                                                              
    VAX/VMS                                                                   
       Drive:  Digital Equipment Corporation (DEC) RRD40 or RRD50             
       Driver: DEC VFS CD-ROM driver V4.7 or V5.2 and up.                     
                                                                              
       Note:   The driver software may be obtained from Jason Hyon at JPL.    
               It is necessary to use this driver to access the XARs on the   
               CD-ROM.                                                        
                                                                              
    VAX/Ultrix                                                                
       Drive:  DEC RRD40 or RRD50                                             
       Driver: Supplied with Ultrix 3.1                                       
                                                                              
       Note:   Internet users can obtain a copy of the "cdio" software        
               package via anonymous ftp from the "space.mit.edu"             
               server in the file named "src/cdio.shar".  Contact Peter       
               Ford for details (see below).                                  
                                                                              
    IBM PC                                                                    
       Drive:  Toshiba, Hitachi, Sony, or compatible.                         
       Driver: Microsoft MSCDEX version 2.2.                                  
                                                                              
       Note:   The latest version of MSCDEX (released in February 1990) is    
               generally available.  Contact Jason Hyon for assistance in     
               locating a copy.                                               
                                                                              
    Apple Macintosh                                                           
       Drive:  Apple CD SC (Sony) or Toshiba                                  
       Driver: Apple CD-ROM driver                                            
                                                                              
       Note:   The Toshiba drive requires a separate driver, which may be     
               obtained from Toshiba.                                         
                                                                              
    Sun Micro (SunOS 4.0.x and earlier)                                       
       Drive:  Delta Microsystems SS-660 (Sony)                               
       Driver: Delta Microsystems driver or SUN sr.o Driver                   
                                                                              
       Note:   For questions concerning this driver, contact Denis Down at    
               Delta Microsystems, 415-449-6881.                              
                                                                              
    Sun Micro (SunOS 4.0.x and later)                                         
       Drive:  Sun Microsystems                                               
       Driver: SunOS sr.o driver                                              
                                                                              
       Note:   A patch must be made to SunOS before the Sun driver can        
               access any CD-ROM files containing Extended Attribute          
               Records, e.g. most of the files on this disk. A copy of        
               this patch is available to Internet users via anonymous        
               ftp from the "space.mit.edu" server in the file named          
               "src/SunOS.4.x.CD-ROM.patch".                                  
                                                                              
                                                                              
6.  Whom to Contact for Information                                           
                                                                              
    For questions about how to read the CD-ROM:                               
                                                                              
              Jason J. Hyon                                                   
              MS 525-3610                                                     
              Jet Propulsion Laboratory                                       
              4800 Oak Grove Drive                                            
              Pasadena, CA  91109                                             
              818-306-6054                                                    
                                                                              
              Electronic mail addresses:                                      
              Internet:  jjh345@mipl3.jpl.nasa.gov                            
              NASAmail:  JHYON                                                
              SPAN:      MIPL3::JJH345                                        
              X.400:     (ID:JHYON,PRMD:NASAMAIL,ADMD:TELEMAIL,C:USA)         
                                                                              
    For questions concerning the MIT products and documentation:              
                                                                              
              Peter G. Ford                                                   
              Center for Space Research                                       
              Building 37-601                                                 
              Massachusetts Institute of Technology                           
              Cambridge, MA  02139                                            
              617-253-6485                                                    
                                                                              
              Electronic mail addresses:                                      
              Internet:   pgf@space.mit.edu                                   
              NASAmail:   PFORD                                               
              SPAN:       JPLPDS::PFORD                                       
              X.400:     (ID:PFORD,PRMD:NASAMAIL,ADMD:TELEMAIL,C:USA)         
                                                                              
                                                                              
7.  Cognizant Persons                                                         
                                                                              
  Mike Martin (JPL), Gail Woodward (JPL), Peter Ford (MIT), Robert            
  Mehlman (UCLA), and Jason Hyon (JPL) designed the PDS labels and the        
  tables.  Florence Moss (JPL) and Jason Hyon developed the software to       
  generate the PDS labels and tables.  Raymond Arvidson and Gail              
  Woodward created the catalog templates, with help from Pam Woncik           
  (JPL), Mary Dale-Bannister (Washington University), and Susan Slavney       
  (Washington University).  The data and documentation were generated         
  at MIT by Peter Ford, Gordon Pettengill, Fang Liu, and Joan Quigley         
  from data supplied by the Magellan Project.  Judy Mankin and Bob            
  Canada at DADC designed the cover of this CD-ROM.                           
                                                                              
  This disk was produced by Jason Hyon.                                       
                                                                              
