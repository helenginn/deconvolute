c
c This will (hopefully) write the FT in rho to a MTZ file.
c
      subroutine write_mtz(rho,smin,smax)
c
      use mapdefs
      type (map) rho
      real smin,smax
c
c its based on f2mtz, which contained this header
c
C*************************************************************
C
C     This code is distributed under the terms and conditions of the
C     CCP4 licence agreement as `Part ii)' software.  See the conditions
C     in the CCP4 manual for a copyright statement.
C
C**************************************************************
c
c
      integer   nlprgo
      parameter	(nlprgo=5)
c
      character*30 lsprgo(nlprgo)
      character*1  ctprgo(nlprgo)
      real resol
      integer h,k,l,sortx(5)
c
c Stuff for parser
c
      integer maxtok
c
c Max number of tokens for 'parser'
c
      parameter (maxtok = 2)
      character str*500
      character*4 key, cvalue(maxtok)
      integer ibeg(maxtok), iend(maxtok), ityp(maxtok), idec(maxtok)
      real fvalue(maxtok)      
      integer ntok, dummy
c
c stuff for rdsymm
c
      character spgrnx*10, pgnamx*10
      integer nsymx, nsympx, nspgrx
      real rsymx(4,4,192)
c
c other assorted variables
c
      integer ct
      real adata(nlprgo), rtod
c
c Set up constant to convert radians to degrees
c
      rtod = 180.0/acos(-1.0)
c
c Set up column labels
c
      lsprgo(1)='H'
      lsprgo(2)='K'
      lsprgo(3)='L'
      lsprgo(4)='FC'
      lsprgo(5)='PHIC'
c
c Set up column types
c
      ctprgo(1)='H'
      ctprgo(2)='H'
      ctprgo(3)='H'
      ctprgo(4)='F'
      ctprgo(5)='P'
c
      ntok = maxtok
      nsymx = 0
c
c Fudge to symmetry using the parser and rdsymm
c Note that the true spacegroup is writen to the MTZ file.
c
c For reasons that are less than clear to me, I3 must be specified in the
c format for it to run on a unix machine it wasn't compiled on. JMD
c So
c     write(str,'(A,1x,I)') 'SYMMETRY', rho%true_spgrp
c becomes
c
      write(str,'(A,1x,I3)') 'SYMMETRY', rho%true_spgrp
      call parser (key, str, ibeg, iend, ityp, fvalue, cvalue,
     $     idec, ntok, dummy, .true.)
c reads symmetry from input
         call rdsymm (2, str, ibeg, iend, ityp, fvalue, ntok,
     $        spgrnx, nspgrx, pgnamx, nsymx, nsympx, rsymx)

      call mtzini
c
c Open mtz file 
c
      call lwopen(1,'HKLOUT')
      call lwtitl (1,rho%title,0)   
      call lwcell(1,rho%cell)
c
c Store the column labels in the mtz header:
c
      call lwclab(1, lsprgo, nlprgo, ctprgo, 0)
c
c mtz file number
c lsprgo char*30(nlprgo)  column labels
c nlprgo no. of column labels
c CTPRGO  (nlprgo)  CHARACTER*1  column types 
c IAPPND    INTEGER  =0     replace all existing labels and types
c                    =1     append to the existing labels and types
c
c Store the symmetry in the mtz header:
c
      call lwsymm (1, nsymx, nsympx, rsymx, spgrnx(1:1), nspgrx,
     $     spgrnx, pgnamx)
c
c Store a history header:
c
      call lwhstl (1, ' ')
c
c write all reflections to mtz file...
c
      ct = 0
c
c errr, we should create another rho%array for hkl limits.
c and I'm sure we can do something clever with pointer aliases to sections.
c
      SELECT CASE (rho%lspgrp)
c
      CASE (16,18)
c
c Update sort order
c
        sortx(1)=1
        sortx(2)=2
        sortx(3)=3
        sortx(4)=0
        sortx(5)=0
        call lwsort(1,sortx)
c
c Get HKL
c
        do h = 0,rho%mxyz(1)/2
          do k = 0,rho%mxyz(2)/2
            do l = 0,rho%mxyz(3),2
              adata(1) = h
              adata(2) = k
              adata(3) = l/2
c
c This is a test version for sg 18; the Re and Im are ajacent in l
c Test for systematic abscence.
c
              if(((rho%true_spgrp.eq.23).or.(rho%true_spgrp.eq.197))
     1	           .and.(modulo(h+k+l/2,2).eq.1))cycle
              if((rho%true_spgrp.eq.197) .and.
     1	           ((l/2.lt.h).or.(k.lt.h).or. 
     2	             ((k.eq.h).and.(k.lt.l/2))))cycle
c
              adata(4)=abs(cmplx(rho%sfs(h,k,l),rho%sfs(h,k,l+1)))
              if ( adata(4) .eq. 0.0 ) then
                adata(5)=0.0
              else
                adata(5)=rtod*atan2(rho%sfs(h,k,l+1),rho%sfs(h,k,l))
              end if
              resol=dot_product(
     1          matmul((/h,k,l/2/),rho%recip_metric),(/h,k,l/2/))
c
!!      if( h+k+l/2 .le. 6 )
!!     1  write(6,'(1x,3i3,1x,f11.3,f7.1,f11.8)')
!!     1  h,k,l/2,adata(4),adata(5),resol
c
              if((resol.ge.smin).and.(resol.le.smax))
     1          call lwrefl(1,adata)
              ct=ct+1 
            end do
          end do
        end do
c
      CASE(1)
        do l = 0,rho%mxyz(3)-1
          do k = 0,rho%mxyz(2)-1
            do h = 0,rho%mxyz(1)+1,2
              adata(1) = h/2
              if ( k .le. rho%mxyz(2)/2 ) then
                adata(2) = k
              else
                adata(2) = k - rho%mxyz(2)
              end if
              if ( l .le. rho%mxyz(3)/2 ) then
                adata(3) = l
              else
                adata(3) = l - rho%mxyz(3)
              end if
c
c This is a test version for sg 1; the Re and Im are ajacent in h
c
              adata(4)=abs(cmplx(rho%sfs(h,k,l),rho%sfs(h+1,k,l)))
              if(adata(4).eq.0.0)then
                adata(5) = 0.0
              else
                adata(5) = rtod*atan2(rho%sfs(h+1,k,l), rho%sfs(h,k,l))
              end if
              resol = dot_product(
     1          matmul(adata(1:3),rho%recip_metric),adata(1:3))
              if((resol.ge.smin).and.(resol.le.smax))
     1          call lwrefl(1,adata)
              ct=ct+1 
            end do
          end do
        end do
c
      CASE DEFAULT
        write(6,*) 
     1   'I don''t know how to put this spacegroup into an MTZ file.'
c
      END SELECT
c
c Close mtz file, print full header:
c
      call lwclos(1,3)
      write (6,'(I10, '' input reflections processed'')')ct
c
      return
      end
