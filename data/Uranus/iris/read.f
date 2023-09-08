      real lat(37),pln(37)
      real xlat(37),xp(37),t(37,37),fp(37,37)
      real pi

      pi = 3.141592653589793238

c Read input fields of T,fp.
c NOTE: URANUS .ura file has pressure running opposite to the
c       other files, so switch it to match.

      ilatmax = 37
      ipmax   = 15
      open(unit=1,file='t_smth.ura',status='old')
      read(1,*)
     +  ((xlat(ilat),xp(ip), t(ilat,ip),ilat=1,ilatmax),ip=ipmax,1,-1)
      close(1)
      open(unit=1,file='fp_smth.ura',status='old')
      read(1,*)
     +  ((xlat(ilat),xp(ip),fp(ilat,ip),ilat=1,ilatmax),ip=ipmax,1,-1)
      close(1)

c Switch to ln p[bars] coordinate, and to latitude[deg].
      do i=1,ipmax
        pln(i)=(2.d0-xp(i))*log(10.d0)
      enddo
      do i=1,ilatmax
        lat(i)=(90.d0-xlat(i))
      enddo

c Write output p[mbar], lat[deg]:
      open(unit=3,file='fpara.dat',status='unknown')
      write(3,'(1x,i3,1x,i3)') ilatmax,ipmax
      write(3,'(1x,a9,15(1x,f6.2))') ' p[mbar]:',
     +         (1000.*exp(pln(i)),i=1,ipmax)
      write(3,*) 'lat[deg]'
      do ilat=1,ilatmax
        write(3,'(1x,f7.1,2x,15(1x,f6.4))') 
     +        lat(ilat),(fp(ilat,ip),ip=1,ipmax)
      enddo
      close(3)

      open(unit=3,file='temperature.dat',status='unknown')
      write(3,'(1x,i3,1x,i3)') ilatmax,ipmax
      write(3,'(1x,a9,15(1x,f6.2))') ' p[mbar]:',
     +         (1000.*exp(pln(i)),i=1,ipmax)
      write(3,*) 'lat[deg]'
      do ilat=1,ilatmax
        write(3,'(1x,f7.1,2x,15(1x,f6.2))') 
     +        lat(ilat),(t(ilat,ip),ip=1,ipmax)
      enddo
      close(3)

      end


