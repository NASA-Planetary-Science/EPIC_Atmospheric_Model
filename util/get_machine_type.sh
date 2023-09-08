#!/bin/sh
#
# This shell is hacked from pvmgetarch.sh
#
# Generate hosttype string
#
# This is a heuristic thing that may need to be tuned from time
# to time.  I don't know of a real solution to determining the
# machine type.  Claims to pick one of:
#   afx8, 
#   bfly,bsd386, 
#   convex,convexn,cm2,cm5,cray,craysmp, 
#   decmips,decalpha,dgav,
#   hp300, hpux, 
#   i860, 
#   ksr1, 
#   LINUX,
#   ncube2,next, 
#   paragon,decmips, 
#   rs6000, 
#   sgi,sun2,sun3,sun4,symmetry,
#   titan, 
#   uvax,unknown,
#   vax
# Need to do:
#   ipsc2
#
# T. Dowling, 5/16/02: added Darwin
#
# Notes:
#   1. Local people mess with things.
#   2. It's good to try a few things for robustness.
#
# 08 Apr 1993  Robert Manchek  manchek@CS.UTK.EDU.
#

#
# begin section that may need to be tuned.
#
ARCH=unknown

if [ -f /bin/machine ]; then
	ht="`/bin/machine`"

	case "$ht" in
	alpha )			ARCH=decalpha ;;
	mips )			ARCH=decmips ;;
	esac
fi

if [ -f /bin/uname ]; then
	os="`/bin/uname -s`"
	ht="`/bin/uname -m`"
	osht="$os,$ht"

	case "$osht" in
	SunOS,sun3* )           ARCH=sun3 ;;
	SunOS,sun4* )           ARCH=sun4 ;;
	ULTRIX,RISC )           ARCH=decmips ;;
	ULTRIX,vax )            ARCH=uvax ;;
	AIX,* )                 ARCH=rs6000 ;;
	*HP*,9000/[2345]* )     ARCH=hp300 ;;
	*HP*,9000/[78]* )       ARCH=hpux ;;
	IRIX,* )                ARCH=sgi ;;
	*,alpha )               ARCH=decalpha ;;
	CRSOS,smp )             ARCH=craysmp ;;
	*,paragon )             ARCH=paragon ;;
	dgux,AViiON )           ARCH=dgav ;;
	Linux,* )               ARCH=LINUX ;;
        Darwin,* )              ARCH=Darwin ;;
	esac
fi

if [ "$ARCH" = unknown ]; then
	if [ -f /usr/bin/uname ]; then
		os="`/usr/bin/uname -s`"
		ht="`/usr/bin/uname -m`"
		osht="$os,$ht"

		case "$osht" in
		SunOS,sun3* )           ARCH=sun3 ;;
		SunOS,sun4* )           ARCH=sun4 ;;
		ULTRIX,RISC )           ARCH=decmips ;;
		ULTRIX,vax )            ARCH=uvax ;;
		AIX,* )                 ARCH=rs6000 ;;
		*HP*,9000/[2345]* )     ARCH=hp300 ;;
		*HP*,9000/[78]* )       ARCH=hpux ;;
		IRIX,* )                ARCH=sgi ;;
		*,alpha )               ARCH=decalpha ;;
		CRSOS,smp )             ARCH=craysmp ;;
		*,paragon )             ARCH=paragon ;;
		dgux,AViiON )           ARCH=dgav ;;
		Linux,* )               ARCH=LINUX ;;
        	Darwin,* )              ARCH=Darwin ;;
		esac
	fi
fi

if [ "$ARCH" = unknown ]; then
	if [ -f /bin/arch ]; then
		case "`/bin/arch`" in
		ksr1 ) ARCH=ksr1 ;;
		sun2 ) ARCH=sun2 ;;
		sun3 ) ARCH=sun3 ;;
		sun4 ) ARCH=sun4 ;;
		esac
	fi
fi

if [ "$ARCH" = unknown ]; then

	if [ -f /ultrixboot ]; then
		if [ -f /pcs750.bin ]; then
			ARCH=uvax
		else
			ARCH=decmips
		fi
	else
		if [ -f /pcs750.bin ]; then ARCH=vax; fi
	fi

	if [ -d /usr/alliant ]; then ARCH=afx8; fi
	if [ -f /usr/bin/cluster ]; then ARCH=bfly; fi
	if [ -d /usr/convex ]; then ARCH=convex; fi
	if [ -f /unicos ]; then ARCH=cray; fi
	if [ -f /hp-ux ]; then ARCH=hp300; fi
	if [ -f /usr/bin/getcube ]; then ARCH=i860; fi
	if [ -f /usr/bin/asm56000 ]; then ARCH=next; fi
	if [ -f /etc/vg ]; then ARCH=rs6000; fi
	if [ -d /usr/include/caif ]; then ARCH=rt; fi
	if [ -f /bin/4d ]; then ARCH=sgi; fi
	if [ -f /dynix ]; then ARCH=symmetry; fi
	if [ -f /bin/titan ]; then ARCH=titan; fi

	if [ -f /usr/bin/machine ]; then
		case "`/usr/bin/machine`" in
		i386 ) ARCH=bsd386 ;;
		esac
	fi
fi

if [ "$ARCH" = sun4 -a -f /dev/cm ]; then ARCH=cm2; fi
if [ "$ARCH" = sun4 -a -f /dev/cmni ]; then ARCH=cm5; fi

# T. Dowling 14 Dec 1994:
if [ "$ARCH" = sun4 -a -d /usr/ncube ]; then ARCH=ncube2; fi

# T. Dowling 24 Jul 2000:
if [ "$ARCH" = rs6000 -a -d /usr/lpp/ppe.poe ] ; then ARCH=sp2; fi

if [ "$ARCH" = convex ]; then
	if getsysinfo -f native_default; then
		ARCH=convexn
	fi
fi

#
# Ugh, done.
#

echo $ARCH
exit


