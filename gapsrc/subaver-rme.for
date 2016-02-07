c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
      SUBROUTINE subaver(map1,map2,env1,env2,nsym1,nsym2,
     &    i_check,o_check,use_out,no_trans,ITYPE,rot,r_check)
c     ================================================================
      IMPLICIT NONE
c
c--- for refinement routine
c
      INCLUDE 'average.fcm'
c     =====================
c
c--- local variables
c
      integer     map1,map2,env1,env2
      integer     iindex(3),new_index(3)
      real        aindex(3),x(3),x1(3)
      integer     ngadd,nxt,nout,itype,np,nc
      logical     i_check,o_check,use_in,no_trans,check,use_out
      integer*2   irho
      external    gp3,gp3_a
      integer*2   gp3,gp3_a
      integer     iz,jz,iy,jy,ix,jx,jjz,jjy,jjx,j,nsym1,nsym2
      logical     r_check,iiinout,iinout
      real        rot(3,3)
c
	INTEGER   ONXT
	EXTERNAL IINOUT
	ONXT = -987654
c
      use_in=.true.
      check=.true.
      iiinout=.false.
c
c---- the averaging process is very simple, we are driven by the
c---- output map, for each pixel in that map (accessed in a sequential
c---- fashion, but since map is bricked we move slowly through space)
c---- we pick up all the symm related pixels in the input map and take
c---- the average.
c
c---- firstly set up pointers for the maps that drive the process (backwards!)
c
      call init_map(map2)
c
c---- initialize r-factors, corr coefs, etc
c
      ngadd=0
c
      call i_stat
c
      if(o_check)call init_map(env2)
c
c
c---- off to work...
c---- =============
c
c---- loop over bricks
c
      do iz=1,n_brick(3,map2)
c
        jz=(iz-1)*brick(3,map2)
c
        do iy=1,n_brick(2,map2)
c
          jy=(iy-1)*brick(2,map2)
c
          do ix=1,n_brick(1,map2)
c
            jx=(ix-1)*brick(1,map2)
c
c---- now loop within each brick
c
            do jjz=jz+1,jz+brick(3,map2)
c
              if(jjz.le.nx(3,map2))then
c
                do jjy=jy+1,jy+brick(2,map2)
c
                  if(jjy.le.nx(2,map2))then
c
                    do jjx=jx+1,jx+brick(1,map2)
c               
                      if(jjx.le.nx(1,map2))then
c
c---- check that this pixel is in required region of the output map
c
	IF (IY.EQ.1 .AND. IX.EQ.1 .AND. JJZ.EQ.JZ+1 .AND.
     +      JJY.EQ.JY+1 .AND. JJX.EQ. JX+1) THEN
	  print *,'IZ: ',IZ
	ENDIF
c
      if(o_check)then
        if(gp3(env2,jjx,jjy,jjz).eq.0)goto 200
      endif
c
      iindex(1)=jjx
      iindex(2)=jjy
      iindex(3)=jjz
c
      call j3_to_f(iindex,x1,map2)
      call f_to_o(x1,x,map2)
c
      in_temp=0
      sum_temp=zero
c
c
c---- loop over all the required NCS operations:
c
      do j=nsym1,nsym2
C
C RME DEBUG: NXT is never set but, perhaps, should really be J
C
	NXT = J
c
c---- get the appropriate postion for this copy:
c----jons note need rot_new(1,1,j) and vec_new(1,j)
        if((.not.no_trans).and.(.not.r_check))
C          write(6,*) 'am in right place'
     &    call spin (x,x1,ops(1,1,j),vecs(1,j))
c
        if(no_trans.and.(.not.r_check))call vmtply(ops(1,1,j),x,x1)
c
        if(r_check.and.cg_check)call spinr(x,x1,ops(1,1,j),
     &                         vecs(1,j),cg(1),rot(1,1))
c
        if(r_check.and.(.not.cg_check))call spinn(x,x1,ops(1,1,j),
     &                         vecs(1,j),rot(1,1))

        call o_to_j3(x1,aindex,iindex,map1)
c
c---- if required check that this is in a permitted region of the Input map:
c
        if(i_check)then
c
C RME DEBUG
C
	IF (nxt.ne.onxt) then
	  print *,'IINOUT:',env1,iindex(1),iindex(2),iindex(3),
     +                      new_index,check,nxt
	  onxt = nxt
	ENDIF
          iiinout=iinout(env1,iindex(1),iindex(2),iindex(3)
     &                            ,new_index,check,nxt)
c          write(6,*)' env check, on input map...,flag:',
c     &                             iiinout
c          (use the envelope pixel that is closest to non-integral grid point)
          if(.not.iiinout.or.gp3_a(env1,iindex).eq.0)then
            nout=nout+1
c---- pixels outside envelope for i/p counted
            goto 300
          endif
c          (use the envelope pixel that is closest to non-integral grid point)
        endif
c
c---- now get the interpolated electron density:
c
c--- if not using xtal symm....much quicker
        if(xr_con(map1))then
c
          call g_p_i_i_x(map1,irho,aindex,iindex,itype
     &                     ,j,use_in,np,nc)
        else
          call g_p_i_i(map1,irho,aindex,iindex,itype,use_in,np,nc)
        endif
c                     ( itype specifies interp: 1, 8 or 11 point. )
c
c---- if the pixel was found ok add it into the sums:
c
        if(.not.use_in) then
          nout=nout+1
        else
          sum_temp=sum_temp+irho
          in_temp=in_temp+1
          ir_temp(in_temp)=irho
          ngadd=ngadd+1
        endif
300     continue
      enddo
c
c---- thats ended the loop over the ncs relationships in map1, we may
c---- want to add in a component from map2, the output map, for present
c---- assume scale factor of 1 between maps
c
      if(use_out) then
c---- check the sym ops for map 2
c
        irho=gp3(map2,jjx,jjy,jjz)
        sum_temp=sum_temp+irho
        in_temp=in_temp+1
        ir_temp(in_temp)=irho
        ngadd=ngadd+1
      endif
c
c---- now we have gathered all the contributions.
c---- a_stat evaluates mean rho and accumualtes stats for this output pixel:
c
      call a_stat
c
      irho=nint(sum_temp)
      nout=nout+1
c
c----no storage of any pixel
c
c
200                     continue
                      endif
                    enddo
                  endif
                enddo
              endif
            enddo
          enddo
        enddo
      enddo
c
c----phew, averaged, no stats at this level..we're refined..
c
c
      return
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
