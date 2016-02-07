      implicit none
C
C
C     ************************************************************
C     *                                                          *
C     *                                                          *
C     *       PROGRAM TO MERGE  MTZ FILES                        *
C     *                                                          *
C     *                                                          *
C     ************************************************************
C
C     MTZ version, EYJ/DIS/SL July-1993.
C
C     This version, DIS, December 1993, to deal with files from R-free.
C       In addition numerous small corrections/changes made.
C
C     Control input uses R.Esnouf's parser.
C
C     12/11/99 JMD Changes for UNIX, and phase extension in the low res 
C       direction
C
C     INPUT stream 11 IS THE FIRST SET OF DATA
C     INPUT stream 12 IS THE SECOND .. ETC..
C     OUTPUT is to stream 2
C
C
      INTEGER MFILES,MCOLS,MBATCH
      PARAMETER (MFILES=3,MCOLS=200,MBATCH=240)
c
      INTEGER NTYP,LTYP
      PARAMETER (NTYP=14,LTYP=8)
c
      INTEGER MBLENG,CBLENG
      PARAMETER (MBLENG=90+56,CBLENG=70+3*8)
c
      INTEGER NHISLM
      PARAMETER (NHISLM=30)
c
      INTEGER SIZE1
      PARAMETER (SIZE1=20)
c
      CHARACTER*10 VERSN
      PARAMETER (VERSN='MTZ:V1.1')
C     ..
C     .. Arrays in Common ..
      REAL CELL,CRANGE,RSYM,WRANGE,RBATR,RBATW
      INTEGER ISORT,NCOLS,NCOLW,NPLABS,NREFR,NREFS,NREFW,NSPGRP,NSYM,
     &        NSYMP,RLUN,RPOINT,WLUN,NBATCH,BATNUM,NBATW,WOMBAT,
     &        NBATR,HDRST,NHISTL
      LOGICAL SORTB
      CHARACTER CTYPE*1,LTYPE*1,SPGNAM*10,CLABEL*30,PLABS*30,TITLE*70,
     &          CBATR*1,CBATW*1,HSCR*80,PGNAM*6
      CHARACTER LSUSRI(MCOLS)*30,LSUSRO(MCOLS)*30
C     ..
C     .. Common blocks ..
      COMMON /MTZCHR/TITLE(MFILES),CLABEL(MCOLS,MFILES),
     &       CTYPE(MCOLS,MFILES),SPGNAM(MFILES),LTYPE(MFILES),
     &       CBATR(CBLENG,MBATCH,MFILES),PGNAM(MFILES)
      COMMON /MTZHDR/CELL(6,MFILES),NSYM(MFILES),NSYMP(MFILES),
     &       RSYM(4,4,96,MFILES),NCOLS(MFILES),NREFS(MFILES),
     &       NBATCH(MFILES),BATNUM(MBATCH,MFILES),
     &       ISORT(5,MFILES),CRANGE(2,MCOLS,MFILES),NSPGRP(MFILES),
     &       RBATR(MBLENG,MBATCH,MFILES)
      COMMON /MTZWRC/PLABS(MCOLS,MFILES),HSCR(NHISLM,MFILES),
     &       CBATW(CBLENG,MBATCH,MFILES)
      COMMON /MTZWRK/NCOLW(MFILES),RLUN(MFILES),WLUN(MFILES),
     &       RPOINT(MCOLS,MFILES),WRANGE(2,MCOLS,MFILES),NREFW(MFILES),
     &       NREFR(MFILES),NPLABS(MFILES),NBATW(MFILES),NBATR(MFILES),
     &       WOMBAT(MBATCH,MFILES),HDRST(MFILES),SORTB(MFILES),
     &       NHISTL(MFILES),RBATW(MBLENG,MBATCH,MFILES)
      COMMON /MTZLBC/LSUSRI,LSUSRO
      COMMON /MTZLAB/NLUSRI,NLUSRO
C     ..

      CHARACTER LSPRGO(MCOLS)*30, CTPRGO(MCOLS)
      INTEGER NLPRGO, ifil, n_col1, n_col2, n_col, npoint1(20)
      INTEGER npoint2(20)
C      DATA LSPRGO/'H','K','L','FP','SIGFP','FCALC','PHIC'/
C      DATA CTPRGO/'H','H','H','F','Q','F','P'/
      LOGICAL  END,ZERO,EOF1,EOF2,LOGMSS(MCOLS)
      DATA  END/.FALSE./,ZERO/.FALSE./
C
      character*1 i_abs
      logical sysabs,out
c
      integer dunit,IPOINT,ncols_read(2)
      character*30 flagco,labco1(20),labco2(20),charval,r_free
      character*30 G2CHRU,extend_file,col_lab(mcols,mfiles),
     &             col_typ(mcols,mfiles)
      logical      flag,rf,terse
      character*7 direction
C
      LOGICAL ADDIN,REFIN,EXTEND
c
      integer nrf,nrfc,nref1,nref2,int1,int2,ilaue,nopt,keyword
      integer intval,j,iprint,ifail,k,nflag,ioff,nmatch
      integer jj,n_down,nadd,nweight,ih,ik,il,nlusri,nlusro
      real    dnow,dmin,dstep,scale1,realval,smin,snow,sstep,sstep3
      real    snew,smin_new,sinner,sininv,dnew,dmin_new,sininv2
      real    dstep_out,resol,resol2,weight,wsum
      real    adata1(mcols),adata2(mcols),adata3(mcols)
      logical auto_res
c
c---- open stream for R.Esnouf's parser
c
      dunit=70
      call getopen('sys$input',dunit)
c
c---- parse command line
c
      call ccpfyp
c
c---- defaults
c
      flagco=' '
      terse=.false.
      flag=.false.
      rf=.false.
      nrf=0
      nrfc=0
      NREF2=0
      INT1=0
      INT2=0
      DNOW=1.0
      DMIN=1000.0
      auto_res=.false.
      DSTEP=50000.0
      EXTEND=.FALSE.
      ADDIN =.FALSE.
      scale1  =1.
      ilaue = 14
      i_abs='P'
      N_COL1=0
      N_COL2=0
      weight=0.0
c
c Added RME 4/12/12: Uninitialized variable
c
      nadd=0
c
      call header
c
c---- accept = character as a blank in the parsing
      call G2EQ
c
c---- READ THE INPUT CARDS
c
1111  nopt=keyword('FSMELT',
     &  'FLAG!INCO!SCAL!ADDI!NOAD!POIN!ABSE!DNOW!DSTE'//
     &  '!RFRE!TERSE!GO  !EXIT!DMIN',0)
c
c---- flag
      if(nopt.eq.1) then
        flag=.true.
        flagco=charval('FSMELT:FLAGCOL',1,30,0)
        write(6,*)'New refls flagged, column: ',flagco
      endif
c
c---- input columns
      if(nopt.eq.2) then
        ifil=intval('FSMELT:FILE_NO',1,2,0)
        n_col=intval('FSMELT:NO_COLS',1,20,0)
        do j=1,n_col
          if(ifil.eq.1)then
            n_col1=n_col
            labco1(j)=charval('FSMELT:COL_LAB',1,30,0)
            write(6,*)'Input Column file 1: ',labco1(j)
          else
            n_col2=n_col
            labco2(j)=charval('FSMELT:COL_LAB',1,30,0)
            write(6,*)'Input Column file 2: ',labco2(j)
          endif
        enddo
      endif
c
c---- Scale factor
      if(nopt.eq.3) then
        scale1=realval('FSMELT:SCALE1',0.0,1000.0,0)
        write(6,*)'Scale factor for file 1:',scale1
      endif
c
c---- add unobserved data
      if(nopt.eq.4) then
        write(6,*)'Addin unobserved data'
        addin=.true.
      endif
c
c---- don't add unobserved data
      if(nopt.eq.5) then
        write(6,*)'Don''t addin unobserved data'
        addin=.false.
      endif
c
c---- point group 
      if(nopt.eq.6) then
        ipoint=INTVAl('FSMELT:POINTGRP',1,9999,0)
        write(6,*)'Point group: ',ipoint
      endif
c
c---- absences
      if(nopt.eq.7)then
c------ note use of G2CHRU to force upper case.
        i_abs=G2CHRU('FSMELT:ABSENCES',1,1,0)
        write(6,*)'Absences ',i_abs
      endif
c
c-----  dnow
        if(nopt.eq.8)then
        dnow=realval('FSMELT:DNOW',-1.0,1000.0,0)
        write(6,*)'DNOW ',dnow
        endif
c
c----   dstep
        if(nopt.eq.9)then
        dstep=realval('FSMELT:DSTEP',-1000000.0,1000000.0,0)
        if(dstep.ne.0.0)then
          extend_file=charval('FSMELT:EXTEND_FILE',1,30,0)
          extend=.true.
          write(6,*)'Extending resolution dstep', dstep
          write(6,*)'Next round''s input written to ',extend_file
          else
          write(6,*)'Resolution not extended'
          endif
        endif
c
c-----  
        if(nopt.eq.10)then
             r_free=charval('FSMELT:R_FREE_COL_LAB',1,30,0)
             write(6,*)'Column for R-free flag: ',r_free
             rf=.true.
        endif
c
c----     terse
        if(nopt.eq.11) then
          terse=.true.
          write(6,*)'Ouput will be terse'
          endif
c
c----     exit
        if(nopt.eq.13) stop
c
c-----  dmin
        if(nopt.eq.14)then
        dmin=realval('FSMELT:DMIN',0.0,1000000.0,0)
        write(6,*)'DMIN ',dmin
        endif
c
c----     GO
        if(nopt.ne.12) goto 1111
c
c--- Now we get on with it
c
c
c--- Call mtz initialisation 
c
        CALL MTZINI
        IPRINT=1
        if(terse)iprint=0
c
c--- Open 1st input file...usually this is the master Fobs file. 
c
        CALL LROPEN(1,'HKLIN1',IPRINT,IFAIL)
        IF(IFAIL.EQ.1) GOTO 3000
c
c--- Sort out resolution range
c
c--- if dnow < 0, set range to that of file 1
c--- else
c
        if(dnow.lt.0.0)then
          auto_res=.true.
          call lrrsol(1,smin,snow)
          snow=sqrt(snow)
          dnow=1/snow
          smin=sqrt(smin)
          dmin=1/smin
        else
          snow=1/dnow
          smin=1/dmin
        endif
        print*,'DNOW =',dnow
        print*,'DMIN =',dmin
c
c---- Phase extension
c       
c---- sstep is 2sintheta/lambda increment
c---- snow is current maximum 2sintheta/lambda limit
c---- smin is current minimum 2sintheta/lambda limit
c
        if(extend)then
                sstep=1/dstep
                sstep3=sstep*3.0
c
                if(dstep.gt.0.0)then
                        snew=sstep+snow
                        smin_new=smin
                        sinner=snow-sstep3
                        direction='higher'
                else
                        snew=snow
                        smin_new=sstep+smin
                        sinner=smin-sstep3
                        direction=' lower'
                endif
c
                sininv=1.0/sinner
                dnew=1.0/snew
                if(smin_new.gt.0.0)then
                        dmin_new=1.0/smin_new
                else
                        dmin_new=10000
                endif
                print *, 'New resolution limits: ',dnew,' to ',dmin_new
c
                snow=snow**2
                snew=snew**2
                smin=smin**2
                smin_new=smin_new**2
                sinner=sinner**2
                sstep=sstep**2
                if(dstep.gt.0)then
                        sininv2=sqrt(1/(snow-sstep))
                else
                        sininv2=sqrt(1/(smin+sstep))
                endif
c
c---- turn off phase extension if its going to fall off the low end
c
                dstep_out=dstep
                if(dstep.lt.0.0)then
                        if(dmin_new.gt.abs(dstep))then
                                dmin_new=abs(dstep)
                                dstep_out=0.0000000000
                        endif
                endif
c
                print*,'Reflections flagged from old limit'
                print*,'Fcalcs at',sininv,' resolution and ',direction,
     &                  ' attenuated'
                print*,'Fcalcs at',sininv2,' resolution and ',direction,
     &                  ' ignored'
c
c---- write file of new params, for cyclic averaging with phase extension.
c
                open(4,file=extend_file,status='UNKNOWN')
                if(scale1.ne.1.0)    write(4,8761)  scale1
                if(addin) then
                  write(4,8762)
                else
                  write(4,8765)
                endif
                write(4,8763)        dstep_out,extend_file
                write(4,8764)        dnew
                write(4,8771)        dmin_new
                if(ipoint.gt.0)      write(4,8766)      ipoint
                write(4,8767)        n_col1,(labco1(j),j=1,n_col1)
                write(4,8768)        n_col2,(labco2(j),j=1,n_col2)
                if(flag)             write(4,8769)      flagco
                if(i_abs.ne.'P')     write(4,878)       i_abs
                if(rf)               write(4,8781)      r_free
                if(terse)            write(4,8782)
                write(4,8770)
                close(4)
8761    format(' SCALE',t20,f10.3)
8762    format(' ADDIN')
8763    format(' DSTEP',t20,f14.3,1x,a30)
8764    format(' DNOW ',t20,f10.5)
8771    format(' DMIN ',t20,f14.5)
8765    format(' NOADDIN ')
8766    format(' POINT ',t20,i10)
8767    format(' INCO  1 ', i4,t20,a30,a30,49(/,t20,a30,a30))
8768    format(' INCO  2 ', i4,t20,a30,a30,49(/,t20,a30,a30))
8769    format(' FLAGCO ',t20,a30)
8770    format(' GO ')
c876     format(i10,2f15.3,f12.5,i10,' /')
c877     format(a1,'            / absences: A,B,C,F,I,R,P are allowed')
878     format(' ABSENCE ',t20,a1,
     &    ' {absences: A,B,C,F,I,R,P are allowed}')
8781    format(' RFREE ',t20,a30)
8782    format(' TERSE ')
c
c---- else if not extending
c
        else
                snow=snow**2
                smin=smin**2
                snew=snow
                smin_new=smin
                sstep=0.0
                sinner=0.0
        endif
c
c###############################################################################
c
c--- Read columns from Fobs file. 
c
d     write(6,*)'Reading column names for HKLIN1'
      call lrclab(1,lsprgo,ctprgo,nlprgo)
d     write(6,*)'Found ',nlprgo,' columns'
      ncols_read(1)=nlprgo
      do j=1,nlprgo
        col_lab(j,1)=lsprgo(j)
        col_typ(j,1)=ctprgo(j)
d       write(6,'(I4,X,A,X,A20)')j,ctprgo(j),lsprgo(j)
      enddo
c
c--- Open 2nd file of F's - this is usually the file of Fcalcs.
c
d     write(6,*)'Opening HKLIN2'
      call lropen(2,'HKLIN2',iprint,ifail)
      if(ifail.eq.1) goto 3000
c
c--- Read columns from HKLIN2 
c
d     write(6,*)'Reading column names for HKLIN2'
      call lrclab(2,lsprgo,ctprgo,nlprgo)
d     write(6,*)'Found ',nlprgo,' columns'
      ncols_read(2)=nlprgo
      do j=1,nlprgo
        col_lab(j,2)=lsprgo(j)
        col_typ(j,2)=ctprgo(j)
d       write(6,'(I4,X,A,X,A20)')j,ctprgo(j),lsprgo(j)
      enddo
c
c--- write out extra column labels for output file
c--- find out which column the r-free label is in:
c
        if(rf)then
          do j=1,ncols_read(1)
            if(col_lab(j,1).eq.r_free) then
c               write(6,*) 'label found for r_free in file 1 ', j
                nrfc=j
                goto 2221
            endif
          enddo
          write(6,*)'%FSMELT_MTZ-ERR: Column not found for r-free'
          stop
2221      continue
        endif
c
c--- Find out which columns we are after for HKLIN1 
c
d     write(6,*)'User asked for ',n_col1,' columns from HKLIN1'
      do k=1,n_col1
d       write(6,*)'Searching for col ',k,' ',labco1(k)
        do j=1,ncols_read(1)
d         write(6,*)'Checking ',col_lab(j,1),' against ',labco1(k)
          if(col_lab(j,1).eq.labco1(k)) then
            if(.not.terse)write(6,'(A8,A17,I4,A2,A8,A11)')labco1(k),
     &        ' found as column ',j,' (',col_lab(j,1),') in file 1'
d           write(6,'(A8,I2,A2,I4)')'npoint1(',k,')=',j
            npoint1(k)=j
            goto 2222
          endif
        enddo
c
        write(6,*)'%FSMELT_MTZ-ERR: Column not found for file 1'
        stop
c
2222    continue
      enddo
c
c--- Find out which columns we are after for HKLIN2
c
d     write(6,*)'User asked for ',n_col2,' columns from HKLIN2'
      do k=1,n_col2
d       write(6,*)'Searching for col ',k,' ',labco2(k)
        do j=1,ncols_read(2)
d         write(6,*)'Checking ',col_lab(j,2),' against ',labco2(k)
          if(col_lab(j,2).eq.labco2(k)) then
            if(.not.terse)write(6,'(A8,A17,I4,A2,A8,A11)')labco2(k),
     &        ' found as column ',j,' (',col_lab(j,2),') in file 2'
d           write(6,'(A8,I2,A2,I4)')'npoint2(',k,')=',j
            npoint2(k)=j
            goto 3333
          endif
        enddo
c
        write(6,*)'%FSMELT_MTZ-ERR: Column not found for file 2'
        stop
c
3333    continue
      enddo
c
c--- Check what we got.
c
d     write(6,'(A6,I4,X,8A8)')'HKLIN1',n_col1,
d    &  (col_lab(npoint1(j),1),j=1,n_col1)
c
c--- Check what we got.
c
d     write(6,'(A6,I4,X,8A8)')'HKLIN2',n_col2,
d    &  col_lab(npoint2(1),2),col_lab(npoint2(2),2)
c
      if(.not.terse)write(6,*)'All column labels assigned o.k.'
c
c--- Number of columns on output file----------------
c   H + K + L + ( no. of cols taken from file one ) + (no. from file 2)
c
c--- Open output file
c
      call lwopen(1,'HKLOUT')
c
c--- If we are having an add-in flag column, find out where it will be
c
      if(flag) then
        nlprgo=(3+n_col1+n_col2+1)
        nflag=3+n_col1+n_col2+1
      else
        nlprgo=(3+n_col1+n_col2)
        nflag=3+n_col1+n_col2
      endif
c
c--- Dump some info
c
        ioff=n_col1+3
        write(6,678)
        write(6,881)npoint1(1),npoint1(2),npoint2(1),npoint2(2),ioff+1
881     format(' Fobs will come from col',i5,' of file1',/,
     &         ' Sig(Fobs)      from col',i5,' of file1',/,
     &         ' Fcalc from col         ',i5,' of file2',/,
     &         ' Phases from col        ',i5,' of file2',/,
     &         ' Fcalc will be go to col',i5,' of output')
        if(rf)write(6,882)nrfc
        if(flag)write(6,883)nflag
882     format(' R-free data from col   ',i5,' of file1')
883     format(' Flag for new data col  ',i5,' of output')
        write(6,678)
c
c--- Set up column names and types for HKLOUT
c
      lsprgo(1)='H'
      ctprgo(1)='H'
      lsprgo(2)='K'
      ctprgo(2)='H'
      lsprgo(3)='L'
      ctprgo(3)='H'
c
c--- Copy over required columns from HKLIN1
c
      do j=1, n_col1
        lsprgo(j+3)=(col_lab(npoint1(j),1))
        ctprgo(j+3)=(col_typ(npoint1(j),1))
      enddo
c
c--- Copy over required columns from HKLIN2
c
      if(n_col2.gt.0)then
        do j=1, n_col2
          lsprgo(j+3+N_col1)=(col_lab(npoint2(j),2))
          ctprgo(j+3+n_col1)=(col_typ(npoint2(j),2))
        enddo
      endif
C
C RME DEBUG
C
	DO J = 1, 3+N_Col1+N_col2
	  print *,lsprgo(j), ctprgo(j)
	ENDDO
c
c--- Add add-in flag column if required
c
      if(flag) then
        ctprgo(nlprgo)='I'
        lsprgo(nlprgo)=flagco
      endif
c
c--- Check to see what we got
c
d      write(6,*)'NLPRGO = ',nlprgo
d      do j=1,nlprgo
d        write(6,'(A7,I4,A4,A8)')'LSPRGO(',j,') = ',lsprgo(j)
d        write(6,'(A7,I4,A4,A8)')'CTPRGO(',j,') = ',ctprgo(j)
d      enddo
c
c--- Write columns to HKLOUT
c
      call lwclab(1,lsprgo,nlprgo,ctprgo,0)
c
c--- Initialise counters and flags
c
      nmatch=0
      nref1=0
      nref2=0
      eof1=.false.
      eof2=.false.
c
c###############################################################################
c
c--- Check we got the right columns
c
d     write(6,'(A8,I4,X,8A8)')'Z HKLIN1',n_col1,
d    &  (col_lab(npoint1(j),1),j=1,n_col1)
d     write(6,'(8I4)')(npoint1(j),j=1,n_col1)
c
c--- Read a reflection from HKLIN1
c
190   CALL LRREFL(1,RESOL,ADATA1,EOF1)
d     write(6,*)'At 190'
      IF (EOF1) goto 2000
d     write(6,'(A6,X,8F10.3)')'Read1>',(adata1(j),j=1,3),
d    &  (adata1(npoint1(j)),j=1,n_col1)
      CALL LRREFM(1,LOGMSS)
c
c--- assume adata1(npoint(1)) will contain Fobs and
c--- adata1(npoint1(2)) will contain sig Fobs
c
c--- skip if Fobs is missing
c
d     if(logmss(npoint1(1)))write(6,*)'HKLIN1 with missing number flag!'
      if(logmss(npoint1(1)))goto190
c
c--- ignore stuff outside resolution ranges
c
      if((resol.gt.snew).or.(resol.lt.smin_new))then
        if(auto_res)write(6,195)ADATA1(1),ADATA1(2),ADATA1(3),RESOL,
     &    SMIN,SNOW
d       write(6,195)ADATA1(1),ADATA1(2),ADATA1(3),RESOL,SMIN,SNOW
195     format('Rejected: ',3f5.0,f10.8,' on cuts of ',2f12.8)
        goto190
      endif
c
      if(rf)adata1(nrfc)=adata1(nrfc)*scale1
c
c--- Increment HKLIN1 reflection counter
c
      nref1=nref1+1
d     write(6,*)'NREF1 =',nref1
c
c--- Read a reflection from HKLIN2
c
200   call lrrefl(2,resol2,adata2,eof2)
d     write(6,*)'At 200'
      if (eof2) goto 1000
d     write(6,'(A6,X,8F10.3)')'Read2>',(adata2(j),j=1,3),
d    &  (adata2(npoint1(j)),j=1,n_col1)
c
c--- Next HKLIN2 reflection if outside the resolution range
c
      if((resol2.gt.snew).or.(resol2.lt.smin_new))goto200
c
c--- Increment HKLIN2 reflection counter
c
      refin=.true.
      nref2=nref2+1
c
c--- Check for match - if HKL from HKLIN1 is larger than from HKLIN2, process unmatched
c--- HKLIN2 reflection - if HKL from HKLIN2 is larger than from HKLIN1, get
c--- another HKLIN1
c
300   continue
d     write(6,*)'At 300'
      do  j=1,3
        if(adata1(j).gt.adata2(j))goto 600
        if(adata1(j).lt.adata2(j))goto 500
      enddo
c
c--- Matched - Copy desired HKLIN1 columns into adata3 for output 
c--- Apply scale factor to HKLIN1 column 1 and 2 (FOBS and SIGF) 
c
      do j=1,3
        adata3(j)=adata1(j)
      enddo
c
      jj=3
      do j=1, n_col1
        jj=jj+1
        if((j.eq.1).or.(j.eq.2)) then
          adata3(jj)=adata1(npoint1(j))*scale1
        else
          adata3(jj)=adata1(npoint1(j))
        endif
      enddo
c
c--- Copy HKLIN2 data into adata3 for output
c
      nmatch=nmatch+1
      do j=1,n_col2
        adata3(j+ioff)=adata2(npoint2(j))
      enddo
c
c--- Here we insert some code to deal properly with the rFree refls - replace FOBS
c--- with FCALC
c
      if(rf)then
c
c---- Check if this reflection is to be excluded.
c
        if(adata1(nrfc).ne.0.0)then
c
c----- Don't treat reflections in phase extension shells
c
          if((resol.gt.SNOW-sstep).or.(resol.lt.smin+sstep))goto450
c
c----- Increment rFree counter
c
          nrf=nrf+1
c
c----- Replace FOBS with FCALC and SIGF with -100
c
          adata3(4)=adata3(1+ioff)
          adata3(5)=-100.0
c
c----- If phase extension is on, downweight extended reflections
c
          if(extend)then
            if(dstep.gt.0.0)then
              if(resol.gt.sinner)then
                adata3(4)=adata3(4)*0.5
                adata3(1+ioff)=adata3(1+ioff)*0.5
                n_down=n_down+1
              endif
            else
              if(resol.lt.sinner)then
                adata3(4)=adata3(4)*0.5
                adata3(1+ioff)=adata3(1+ioff)*0.5
                n_down=n_down+1
              endif
            endif
          endif
c
c----- Increment added-in counter and write reflection 
c
          nadd=nadd+1
          call lwrefl(1,adata3)
          goto 450
c
        endif
      endif
c
c--- Calculate Rayment weights for downweighted region
c
      if(extend)then
        if(dstep.gt.0.0)then
          if(resol.gt.sinner)then
c
c----- Assume first column taken from file two are F's - calculate weight, increment
c----- counter and keep total
c
            if(adata3(4).gt.0)weight=
     &        exp(-abs((adata3(4)-adata3(ioff+1))/adata3(4)))
            nweight=nweight+1
            wsum=wsum+weight
c
c----- Flag extended data by a negative value in flag column.
c
            if(flag)then
              if(abs(resol).gt.snow)then
                adata3(nflag)=-10
              else
                adata3(nflag)=0
              endif
            endif
          endif
c
c---- Else if extending inwards
c
        else
          if(resol.lt.sinner)then
c
c----- Assume first column taken from file two are F's - calculate weight, increment
c----- counter and keep total
c
            if(adata3(4).gt.0)weight=
     &        exp(-abs((adata3(4)-adata3(ioff+1))/adata1(4)))
            nweight=nweight+1
            wsum=wsum+weight
c
c----- Flag extended data by a negative value in flag column.
c
            if(flag)then
              if(abs(resol).lt.smin)then
                adata3(nflag)=-10
              else
                adata3(nflag)=0
              endif
            endif
          endif
        endif
      endif
c
c--- Write this matched, processed reflection to HKLOUT
c
d400  write(6,*)'At 400'
d     if(adata3(5).eq.-100)write(6,'(A21)')'Dodgy real reflection'
d     write(6,'(3F4.0,4F10.3)')(adata3(j),j=1,7)
      call lwrefl(1,adata3)
c
c--- Set the refin flag to indicate we no longer have a reflection from HKLIN2
c
450   refin=.false.
d     write(6,*)'At 450'
c
c--- Read and process another reflection from HKLIN1
c
500   continue
d     write(6,*)'At 500'
      call lrrefl(1,resol,adata1,eof1)
      if (eof1) goto 2000
d     write(6,'(A6,X,8F10.3)')'Read1>',(adata1(j),j=1,3),
d    &  (adata1(npoint1(j)),j=1,n_col1)
c
c--- Check for missing number flags - next reflection if FOBS missing
c
      call lrrefm(1,logmss)
d     if(logmss(npoint1(1)))write(6,*)'HKLIN1 with missing number flag!'
      if(logmss(npoint1(1)))goto500
c
c--- Next reflection if outside resolution range
c
      if((resol.gt.snew).or.(resol.lt.smin_new))goto500
c
c--- Increment counter
c
      nref1=nref1+1
c
c--- Copy desired HKLIN1 columns into adata3 for output 
c--- Apply scale factor to HKLIN1 column 1 and 2 (FOBS and SIGF) 
c
      do j=1,3
        adata3(j)=adata1(j)
      enddo
c
      jj=3
      do j=1, n_col1
        jj=jj+1
        if((j.eq.1.).or.(j.eq.2)) then
          adata3(jj)=adata1(npoint1(j))*scale1
        else
          adata3(jj)=adata1(npoint1(j))
        endif
      enddo
c
      goto 300
c
c--- Process unmatched reflections from HKLIN2 - if we are not adding in, just get
c--- another one - if the refin flag indicates we no longer hold an unprocessed
c--- HKLIN2 reflection, just get another one
c
600   continue
d     write(6,*)'At 600'
      if(.not.addin)goto 200
      if(.not.refin)goto 200
c
c--- Don't add in in extension bands
c
      if(extend)then
        if(dstep.gt.0.0)then
          if(resol2.gt.snow-sstep)goto200
        else
          if(resol2.lt.smin+sstep)goto200
        endif
      else
        if(resol.gt.dnow)goto200
      endif
c
c--- Unmatched HKLIN2 reflection for adding in - check quadrant
c
      ih=adata2(1)
      ik=adata2(2)
      il=adata2(3)
c
c--- If point group 3 the h0l reflections are not right - should be 0 h -l
c
      if(ipoint.eq.3) then
        if(ik.eq.0)then
          if(ih.ne.0)then
            ik=ih
            ih=0
            il=-il
          endif
        endif
      endif
c
c--- Check to see if this reflection is a systematic absence
c
      if(sysabs(ih,ik,il,i_abs))goto 200
      if(out(ih,ik,il,ipoint))goto 200
c
c--- On outer range scale fcalcs by 0.5
c
      if(extend)then
        if(dstep.gt.0.0)then
          if(resol2.gt.sinner)then
            adata2(npoint2(1))=adata2(npoint2(1))*0.5
            n_down=n_down+1
          endif
        else
          if(resol2.lt.sinner)then
            adata2(npoint2(1))=adata2(npoint2(1))*0.5
            n_down=n_down+1
          endif
        endif
      endif
c
c---------------------------------
c           aDATA2(7)=aDATA2(4)
c           aDATA2(8)=aDATA2(5)
c           aDATA2(6)= 0.0
c           aDATA2(5)=-100.0
c           aDATA2(6)=aDATA2(4)
c           aDATA2(7)=aDATA2(5)
c           aDATA2(5)=-100.0
c           adata2(8)=0
c---------------------------------
c
c--- Prepare output record, first copy over h,k,l:
c
      do j=1,3
        adata3(j)=adata2(j)
      enddo
c
c--- No information from file one so empty those cols in output record
c
      jj=3
      do j=1,n_col1
        jj=jj+1
        adata3(jj)=0.0
      enddo
c
c--- Copy over stuff from HKLIN2
c
      do j=1,n_col2
        adata3(j+3+n_col1)=adata2(npoint2(j))
      enddo
c
c--- Reset dummy FOBS and SIGF:
c
      adata3(4)=adata2(npoint2(1))
      adata3(5)=-100.0
      if(flag)adata3(nflag)=0.0
c
c--- Increment counter and write reflection
c
      nadd=nadd+1
      call lwrefl(1,adata3)
d     write(6,*)'Written Added-In Refl ',nadd
c
c--- Get another reflection from HKLIN2
c
      Goto 200
c
c--- End of HKLIN1 reached - if adding in, set H to an unobtainable
c---  value to get all remaining reflections from HKLIN2
c
2000  continue
      if(.not.addin)goto 1000
      end=.true.
      adata1(1)=10000
      goto 200
c
c--- End of HKLIN2 reached - close HKLOUT
c
1000  continue
      call lwclos(1,iprint)
C
C RME DEBUG
C
	DO J = 1, 3+N_Col1+N_col2
	  print *,lsprgo(j), ctprgo(j)
	ENDDO
c
c--- Dump some stats
c
      write(6,678)
      print*,  'NUMBER OF REFLECTIONS READ 1 :',nref1
      print*,  'NUMBER OF REFLECTIONS READ 2 :',nref2
      print*,  'NUMBER OF REFLECTIONS MATCHED:',nmatch
      write(6,1500)100.0*float(nmatch)/float(nref2)
1500  format (' MATCHED AS PERCENT OF READ 2 :         ',f6.2)
      print*,  'NUMBER OF SYNTHETIC FS ADDED :',NADD
      if(extend)then
        print*,'NUMBER OF FOBS EXTENDED      :',nweight
        if(nweight.ne.0)wsum=wsum/nweight
        print*,'AV. RAYMENT WEIGHT WOULD BE  :',wsum
        print*,'NUMBER OF FC"S DOWNWEIGHTED  :',n_down
      endif
      if(rf)print*,  'NUMBER OF R_FREE FO"S LOST   :',nrf
      write(6,678)
c
c--- The end
c
      stop
c
c--- This only happens if there is a file read error
c
3000  stop '%FSMELT_MTZ-ERR error in reading mtz file'
c
c--- That's all folks!
c
c
678   format(t20,'======================================')
      end
c
c#######################################################################
c
      FUNCTION sysabs(h,k,l,i_abs)
      IMPLICIT NONE
      integer h,k,l,isum
      character*1 i_abs
      logical sysabs
      sysabs=.false.
c
      if(i_abs.eq.'P')return
      if(i_abs.eq.'I')then
        isum=h+k+l
        sysabs=isum/2*2.ne.isum
        return
      endif
      if(i_abs.eq.'A')then
        sysabs=(k+l)/2*2.ne.(k+l)
        return
      endif
      if(i_abs.eq.'B')then
        sysabs=(h+l)/2*2.ne.(h+l)
        return
      endif
      if(i_abs.eq.'C')then
        sysabs=(h+k)/2*2.ne.(h+k)
        return
      endif
      if(i_abs.eq.'R')then
        sysabs=(-h+k+l)/3*3.ne.(-h+k+l)
        return
      endif
      if(i_abs.eq.'F')then
        sysabs=.true.
        if((h+k)/2*2.ne.(h+k))return
        if((h+l)/2*2.ne.(h+l))return
        if((k+l)/2*2.ne.(k+l))return
        sysabs=.false.
        return
      endif
      return
      end
c
c#######################################################################
c
      FUNCTION out(h,k,l,ipoint)
      IMPLICIT NONE
      integer h,k,l,ipoint
      logical out
      out=.false.
      if(ipoint.le.1)return
      if(ipoint.eq.2)then
        out=(k.lt.0).or.(l.lt.0)
        return
      endif
      if(ipoint.eq.222)then
        out=(h.lt.0.or.k.lt.0.or.l.lt.0)
        return
      endif
      if(ipoint.eq.23)then
        out=(h.lt.0.or.k.lt.0.or.l.lt.0)
        if(out)return
        out=(k.gt.l).or.(h.gt.l)
        return
      endif
      if(ipoint.eq.3)then
        out=(h.lt.0).or.(k.le.0)
        if(h.eq.0.and.k.eq.0)then
          out=(l.lt.0)
        endif
        return
      endif
      if(ipoint.eq.32)then
        out=(h.lt.0.or.k.lt.0)
        if(out)return
        out=(l.ge.0).and.(k.gt.h)
        if(out)return
        out=(l.lt.0).and.(k.ge.h)
        return
      endif
      if(ipoint.eq.432)then
        out=(h.lt.k.or.k.lt.l.or.l.lt.0)
        return
      endif
      end
c
c#######################################################################
c
      SUBROUTINE header
c     =================
      IMPLICIT NONE
c
      character*20 date
      character*18 time
      integer      ilen
c
c--- Get date and time from Robert's parser
c
      call p2time(time,ilen)
      call p2date(date,ilen)
c
c--- Print jolly header
c
      write(6,10) time,date
10    format(/,
     &  ' ====================================================',/,
     &  ' Fsmelt_mtz,                 version 1.3, 12-Dec-2002',/,
     &  ' ====================================================',/,
     &  ' Run at    ',a,'  on  ',a,/,
     &  ' ====================================================',/)
c
c--- That's all folks!
c
      return
      end
