!=====================================================================
! AUTHORS:
!  Gary Egbert & Lana Erofeeva
!  College of Atmospheric and Oceanic Sciences
!  104 COAS Admin. Bldg.
!  Oregon State University
!  Corvallis, OR 97331-5503
!  
!  E-mail:  egbert@coas.oregonstate.edu                                      
!  Fax:     (541) 737-2064
!  Ph.:     (541) 737-2947                                        
!  https://www.tpxo.net/
!
! COPYRIGHT: OREGON STATE UNIVERSITY, 2010
! (see the file COPYRIGHT for lisence agreement)
!=====================================================================
      program extract_HC_da
!cc   LANA, 2020 remake (for netcdf files)
!cc   needs min RAM and time to extract/predict by accessing
!cc   4 model nodes corresponding to given lat/lon cell 
!cc   reads OTIS netcdf model file
!cc   (elevations OR transports), reads a list of locations 
!cc   and outputs ASCII file with the complex amplitudes/amp,phases of
!cc   elevations/transports/currents interpolated at the locations
! 
      implicit none
      include 'netcdf.inc'
      include 'constit.h'
      complex, allocatable:: z1(:),dtmp(:,:)
      complex, allocatable:: zl1(:)
      complex d1
      real, allocatable:: lat(:),lon(:),depth(:,:),x(:),y(:),lon0(:)
      real, allocatable:: phase(:)
      real latp,lonp,xt,yt
      real th_lim(2),ph_lim(2),dum,lth_lim(2),lph_lim(2)
      integer, allocatable:: cind(:),lcind(:),mz(:,:)
!
      character*4 c_id(ncmx),c_id_mod(ncmx),lc_id(ncmx)
      character*80 modname,coo_file,outname,gname,ctmp,lname
      character*80 hname,uname,fname
      character*2000 fmt
      character*80 rmCom
      character*1 zuv,c1
      character*80 xy_ll_sub
      logical APRI,geo,interp_micon,ll_km
      integer ncon,nc,n,m,ndat,k,ierr,ierr1,ic,n0,m0
      integer ncl,nl,ml,ibl,nca,l
      integer*4 status,ncid(ncmx),ncidl(ncmx),ncid_new(ncmx)
      integer narg,dimid,len_coo,len_xb,len_yb,varid,coo_idx,row,col
      real, allocatable:: nav_lon(:,:),nav_lat(:,:),z1_re(:,:),z1_imm(:,:)
      real, allocatable:: tides_reout(:,:,:),tides_immout(:,:,:)
      integer, allocatable:: coo_array(:),x_idx(:),y_idx(:)
      integer varid_field1,varid_field2
      character*2 newfield_re,newfield_imm
      character*1 lower
      integer ul_diff
      character*2 upper_ic
      character*80 outname_ic
!
      ll_km=.false.
!
      narg=iargc()
      WRITE(6,*) "narg",narg
      ! NEW: read as line args the var z,u or v,
      ! the coordinates input file (T, U or V 2d bdy NEMO
      ! files) and the out file name
      ! If there is the arg (should be the date in format: yyyymmdd)
      if(narg.eq.3)then
        ! Read the args, namely the var and the coo file 
        call getarg(1,zuv)
        WRITE(6,*) "zuv var: ",zuv
        call getarg(2,coo_file)
        WRITE(6,*) "coo file: ",coo_file
        call getarg(3,outname)
        WRITE(6,*) "Out file: ",outname
        !call getarg(4,mesh_mask)
        !WRITE(6,*) "mesh_mask file: ",mesh_mask
      else
       WRITE(6,*) "zuv var, coo file and outfile &
                   MUST be given as arguments!"
       stop
      endif
!
      ibl=0
      lname='DATA/load_file.nc'
      call rd_inp(modname,c_id,ncon,APRI,geo,interp_micon)
      call rd_mod_file(modname,hname,uname,gname,xy_ll_sub,nca,c_id_mod)
      write(*,*)
      write(*,*)'Lat/Lon file:',trim(coo_file)
      if(ncon.gt.0)write(*,*)'Constituents to include: ',c_id(1:ncon)
      if(geo.and.zuv.eq.'z')then
         write(*,*)'Extract GEOCENTRIC tide HC'
      else
         write(*,*)'Extract OCEAN tide HC'
      endif
!     Open outfile
      open(unit=11,file=outname,status='unknown')
!
!
      write(*,*)
      call rd_mod_header_nc(modname,zuv,n,m,th_lim,ph_lim,nc,c_id_mod,&
                            xy_ll_sub)
      write(*,*)'Lat/Lon file:',trim(coo_file)
      write(*,*)'Model:        ',trim(modname(12:80))
      write(11,*)'Model:        ',trim(modname(12:80))
      if(trim(xy_ll_sub).eq.'')then
       write(*,*)'Lat limits:   ',th_lim
       write(*,*)'Lon limits:   ',ph_lim
      else
       ll_km=.true.
       if(trim(xy_ll_sub).ne.'xy_ll_N'.and.&
          trim(xy_ll_sub).ne.'xy_ll_S'.and.&
          trim(xy_ll_sub).ne.'xy_ll_CATs')then
        write(*,*)'No converting function ', trim(xy_ll_sub),&
                  'in the OTPS'
        stop 
       endif
       if(trim(xy_ll_sub).eq.'xy_ll_N')then
        call xy_ll_N(ph_lim(1),th_lim(1),xt,yt)
        if(xt.gt.180)xt=xt-360
        write(*,*)'Lat,Lon lower left corner:',yt,xt
        call xy_ll_N(ph_lim(1),th_lim(2),xt,yt)
        if(xt.gt.180)xt=xt-360
        write(*,*)'Lat,Lon upper left corner:',yt,xt
        call xy_ll_N(ph_lim(2),th_lim(1),xt,yt)
        if(xt.gt.180)xt=xt-360
        write(*,*)'Lat,Lon lower right corner:',yt,xt
        call xy_ll_N(ph_lim(2),th_lim(2),xt,yt)
        if(xt.gt.180)xt=xt-360
        write(*,*)'Lat,Lon upper right corner:',yt,xt
       elseif(trim(xy_ll_sub).eq.'xy_ll_S')then
        call xy_ll_S(ph_lim(1),th_lim(1),xt,yt)
        if(xt.gt.180)xt=xt-360
        write(*,*)'Lat,Lon lower left corner:',yt,xt
        call xy_ll_S(ph_lim(1),th_lim(2),xt,yt)
        if(xt.gt.180)xt=xt-360
        write(*,*)'Lat,Lon upper left corner:',yt,xt
        call xy_ll_S(ph_lim(2),th_lim(1),xt,yt)
        if(xt.gt.180)xt=xt-360
        write(*,*)'Lat,Lon lower right corner:',yt,xt
        call xy_ll_S(ph_lim(2),th_lim(2),xt,yt)
        if(xt.gt.180)xt=xt-360
        write(*,*)'Lat,Lon upper right corner:',yt,xt
       elseif(trim(xy_ll_sub).eq.'xy_ll_CATs')then
        call xy_ll_CATs(ph_lim(1),th_lim(1),xt,yt)
        if(xt.gt.180)xt=xt-360
        write(*,*)'Lat,Lon lower left corner:',yt,xt
        call xy_ll_CATs(ph_lim(1),th_lim(2),xt,yt)
        if(xt.gt.180)xt=xt-360
        write(*,*)'Lat,Lon upper left corner:',yt,xt
        call xy_ll_CATs(ph_lim(2),th_lim(1),xt,yt)
        if(xt.gt.180)xt=xt-360
        write(*,*)'Lat,Lon lower right corner:',yt,xt
        call xy_ll_CATs(ph_lim(2),th_lim(2),xt,yt)
        if(xt.gt.180)xt=xt-360
        write(*,*)'Lat,Lon upper right corner:',yt,xt
       endif
      endif
      write(*,*)'Constituents: ',c_id_mod(1:nc)
      if(trim(xy_ll_sub).ne.'')then
           write(*,*)'Model is on uniform grid in km'
           write(*,*)'Function to convert x,y to lat,lon:',&
                      trim(xy_ll_sub)
      endif 
!
      if(zuv.eq.'z')write(*,*)'Output elevations (m)'
      if(zuv.eq.'U')write(*,*)'Output WE transport (m^2/s)'
      if(zuv.eq.'u')write(*,*)'Output WE velocity  (cm/s)'
      if(zuv.eq.'V')write(*,*)'Output SN transport (m^2/s)'
      if(zuv.eq.'v')write(*,*)'Output SN velocity  (cm/s)'
!
      if(zuv.eq.'z')write(11,*)'Elevations (m)'
      if(zuv.eq.'U')write(11,*)'WE transport (m^2/s)'
      if(zuv.eq.'u')write(11,*)'WE velocity  (cm/s)'
      if(zuv.eq.'V')write(11,*)'SN transport (m^2/s)'
      if(zuv.eq.'v')write(11,*)'SN velocity  (cm/s)'
!
      if(ncon.eq.0)then
       ibl=1
       ncon=nc
       c_id=c_id_mod
       write(*,*)'Constituents to include: ',c_id(1:ncon)
      endif
!
      allocate(cind(ncon))
      call def_con_ind(c_id,ncon,c_id_mod,nc,cind)
!     Read lat/lon i coo_file
      write(6,*)'Open and read coo_file: ',coo_file
      status=nf_open(trim(coo_file),nf_nowrite,ncid)
      if(status.eq.0)then
       ! Read longitude len
       if(zuv.eq.'z')then
          write(6,*)'Reading xbT dimension in the input file.. '
          status=nf_inq_dimid(ncid,'xbT',dimid)
       else if(zuv.eq.'u')then
          write(6,*)'Reading xbU dimension in the input file.. '
          status=nf_inq_dimid(ncid,'xbU',dimid)
       else if(zuv.eq.'v')then
          write(6,*)'Reading xbV dimension in the input file.. '
          status=nf_inq_dimid(ncid,'xbV',dimid)
       end if
       status=nf_inq_dimlen(ncid,dimid,len_coo)
       len_xb=len_coo
       ! Read latitude len
       status=nf_inq_dimid(ncid,'yb',dimid)
       status=nf_inq_dimlen(ncid,dimid,len_coo)
       len_yb=len_coo
      else
       write(6,*)'ERROR: Something wrong with the input NEMO bdy file..',status
       stop
      endif
      ! Open mesh_mask file
      !write(6,*)'Open mesh_mask: ',mesh_mask
      !status=nf_open(trim(mesh_mask),nf_nowrite,ncid2)
      !if(status.ne.0)then
      !  write(6,*)'ERROR: Something wrong with mesh_mask file..',status
      !  stop
      !end if
      ! Grid points number
      ndat=len_xb*len_yb
      write(6,*)'Grid points number: ',len_xb,' x ',len_yb,' = ',ndat
      ! Read coo values from bdy NEMO files
      allocate(lat(ndat),lon(ndat),nav_lon(len_xb,len_yb),&
              nav_lat(len_xb,len_yb),lon0(len_xb),coo_array(ndat),&
              x_idx(len_xb),y_idx(len_yb))
              !mask(ndat),tmask(len_xb,len_yb,len_lev,1))
      status=nf_inq_varid(ncid,'nav_lon',varid)
      status=nf_get_var(ncid,varid,nav_lon)
      status=nf_inq_varid(ncid,'nav_lat',varid)
      status=nf_get_var(ncid,varid,nav_lat)
      !status=nf_inq_varid(ncid2,'tmask',varid)
      !status=nf_get_var_int(ncid2,varid,tmask)
      ! 
      coo_idx = 1
      do col = 1, len_yb
       do row = 1, len_xb
         lon(coo_idx)=nav_lon(row,col)
         lat(coo_idx)=nav_lat(row,col)
         !mask(coo_idx)=tmask(row,col,1,1)
         coo_array(coo_idx)=coo_idx
         coo_idx=coo_idx+1
       enddo
      enddo
      write(6,*)'coo tot: ',coo_idx-1
      ! Redefn of coordinates..
      write(6,*)'Start interpolation..  '
      if(trim(xy_ll_sub).ne.'')allocate(x(ndat),y(ndat))
      !lon0=lon
      do k=1,ndat
       !write(6,*)'Point 2 be interp:  ',lon(k),lat(k)
       if(trim(xy_ll_sub).eq.'xy_ll_N')then
         call ll_xy_N(lon(k),lat(k),x(k),y(k))
       elseif(trim(xy_ll_sub).eq.'xy_ll_S')then
         call ll_xy_S(lon(k),lat(k),x(k),y(k))
       elseif(trim(xy_ll_sub).eq.'xy_ll_CATs')then
         call ll_xy_CATs(lon(k),lat(k),x(k),y(k))
       endif
       lon0(k)=lon(k)
       if(trim(xy_ll_sub).eq.'')then ! check on lon convention
        if(lon(k).gt.ph_lim(2))lon(k)=lon(k)-360
        if(lon(k).lt.ph_lim(1))lon(k)=lon(k)+360
       endif
      enddo
      write(6,*)'..End interpolation'
      allocate(z1(ncon))
      allocate(z1_re(ncon,k))
      allocate(z1_imm(ncon,k))
!
      if(zuv.eq.'z'.and.geo)then ! NOT our case..
       write(*,'(a,$)')'Reading load correction header...'
       call rd_mod_header1_nc(lname,nl,ml,ncl,lth_lim,lph_lim,lc_id)
       allocate(lcind(ncon))
       call def_con_ind(c_id,ncon,lc_id,ncl,lcind)
       allocate(zl1(ncon))
       write(*,*)'done'
       status=nf_open(trim(lname),nf_nowrite,ncidl(1))
       if(status.ne.0)then
        write(*,*)'Failed to open file:',trim(lname)
        stop
       endif
      endif
      if(zuv.eq.'u'.or.zuv.eq.'v')then
       write(*,'(a,$)')'Reading grid file...'
       allocate(depth(n,m),dtmp(n,m),mz(n,m))
! currents case: need to read grid
       call rd_grd_nc(gname,n,m,depth,mz)
       dtmp=depth
       deallocate(depth)
       write(*,*)'done'
      endif
      allocate(phase(ncon))
! output format
      fmt='(f9.3,x,f9.3,x'
      c1=','
      do ic=1,ncon
       if(ic.eq.ncon)c1=')'
       if(APRI)then
        fmt=trim(fmt)//'f8.3,x,f8.1,x'//c1
       else
        fmt=trim(fmt)//'f8.3,x,f8.3,x'//c1
       endif
      enddo
!
      if(APRI)then
        write(11,*)'  Lat     Lon       ',&
                    (trim(c_id(ic)),'_amp  ',&
                     trim(c_id(ic)),'_ph   ',ic=1,ncon)
      else
        write(11,*)'  Lat       Lon     ', &
                   (trim(c_id(ic)),'_Re   ', &
                    trim(c_id(ic)),'_Im   ',ic=1,ncon)
      endif
      c1=zuv ! since interp change zuv (U->u, V->v)
      latp=0.
      lonp=0.
!
      fname=hname
      if(zuv.ne.'z')fname=uname
      if(nca.eq.0)then
       status=nf_open(trim(fname),nf_nowrite,ncid(1))
       if(status.ne.0)go to 4
      else
       write(*,'(a,$)')'Opening atlas files:'
       do ic=1,ncon
        if(ic.gt.1)then
          k=index(fname,trim(c_id(ic-1)))
          l=len(trim(c_id(ic-1)))
          fname=fname(1:k-1)//trim(c_id(ic))//fname(k+l:80)
        endif
        write(*,'(a,$)')c_id(ic)
        status=nf_open(trim(fname),nf_nowrite,ncid(ic))
        if(status.ne.0)go to 4
       enddo
       write(*,*)'done'
      endif  
!
      do k=1,ndat
       if(lat(k).ne.latp.or.lon(k).ne.lonp)then
        if(trim(xy_ll_sub).eq.'')then
         call interp_da_nc(ncid,n,m,th_lim,ph_lim, &
                        lat(k),lon(k),z1,ncon,cind,ierr,c1,nca)
        else
         z1(1)=-1
         call interp_da_nc(ncid,n,m,th_lim,ph_lim, &
                        y(k),x(k),z1,ncon,cind,ierr,c1,nca)
        endif
        if(ierr.eq.0)then
         if(zuv.eq.'u'.or.zuv.eq.'v')then
          if(trim(xy_ll_sub).eq.'')then
           call interp(dtmp,1,n,m,mz,th_lim,ph_lim, &
                       lat(k),lon(k),d1,ierr1,'z')
          else
           d1=-1
           call interp(dtmp,1,n,m,mz,th_lim,ph_lim, &
                       y(k),x(k),d1,ierr1,'z')
          endif
          z1=z1/real(d1)*100. ! currents cm/s
         elseif(zuv.eq.'z'.and.geo)then       
           call interp_da_nc(ncidl,nl,ml,lth_lim,lph_lim, &
                    lat(k),lon(k),zl1,ncon,lcind,ierr1,'z',0)
           z1=z1+zl1     ! apply load correction to get geocentric tide
         endif  
        endif
        if(ierr.eq.0)then
         if(APRI)then ! NOT our case
          phase=atan2(-imag(z1),real(z1))*180/3.141593
          write(11,fmt)lat(k),lon0(k), &
                (abs(z1(ic)),phase(ic),ic=1,ncon)
         else
           do ic=1,ncon
              z1_re(ic,k)=real(z1(ic))
              z1_imm(ic,k)=imag(z1(ic))
           end do
           write(11,fmt)lat(k),lon0(k), &
                    (real(z1(ic)),imag(z1(ic)),ic=1,ncon)
         endif
        else
          ! If on the land values must be 0.000
           z1=0.000
           do ic=1,ncon
              z1_re(ic,k)=real(z1(ic))
              z1_imm(ic,k)=imag(z1(ic))
           end do
          write(11,fmt)lat(k),lon0(k), &
                    (real(z1(ic)),imag(z1(ic)),ic=1,ncon)
          !write(11,'(1x,f8.3,x,f8.3,a70)')lat(k),lon0(k), &
       !'************* Site is out of model grid OR land ***************'
        endif
        latp=lat(k)
        lonp=lon(k)
       endif  
      enddo
      ! Conversion from coo to nav_lat nav_lon
      write(*,*)'2D Conversion..'
      allocate(tides_reout(len_xb,len_yb,ncon))
      allocate(tides_immout(len_xb,len_yb,ncon))
      do ic=1,ncon
       coo_idx = 1
       do col = 1, len_yb
        do row = 1, len_xb
         tides_reout(row,col,ic)=z1_re(ic,coo_idx)
         tides_immout(row,col,ic)=z1_imm(ic,coo_idx)
         coo_idx=coo_idx+1
        enddo
       enddo
      enddo
!      ! OUT NC FILE
!      do ic=1,ncon
!         !lower = ichar(trim(c_id(ic)(1:1)))
!         !ul_diff=ichar('A')-ichar('a') 
!         !upper_ic=char(lower+ul_diff)//trim(c_id(ic)(2:2))
!         !write(6,*)'Prova: ',lower,ul_diff,upper_ic
!         ! Open the outfiles and add the fields
!         outname_ic=outname ! Da costruire il nome dell'outfile+cambia arg + funzione per convertire in upper le tidal comp
!         write(6,*)'Open outfile: ',outname_ic
!         ! Open  outfile in writing mode
!         status=nf_open(trim(outname_ic),nf_netcdf4,ncid_new)
!         status=nf_redef(ncid_new)
!         ! Define new fields
!         newfield_re=trim(zuv)//'1'
!         status = nf_def_var(ncid_new, newfield_re, nf_float,2,[len_xb,len_yb], varid_field1)
!         newfield_imm=trim(zuv)//'2'
!         status = nf_def_var(ncid_new, newfield_imm, nf_float,2,[len_xb,len_yb], varid_field2)
!         ! Write values in the new fields
!         status=nf_put_vara(ncid_new,varid_field1,[1,1],[len_xb,len_yb],tides_reout(:,:,ic))
!         status=nf_put_vara(ncid_new,varid_field2,[1,1],[len_xb,len_yb],tides_immout(:,:,ic))
!         ! Close the outfile
!         status = nf_close(ncid_new)
!      end do


      deallocate(z1,cind,lat,lon,lon0)
      if(zuv.eq.'u'.or.zuv.eq.'v')deallocate(dtmp,mz)
      if(trim(xy_ll_sub).ne.'')deallocate(x,y)
      if(ibl.eq.1.and.geo)then
         ncon=0
         deallocate(zl1,lcind)
      endif
      close(11)
      status=nf_close(ncid(1))
      if(nca.ne.0.and.ncon.gt.1)then
       do ic=2,ncon
        status=nf_close(ncid(ic))
       enddo
      endif
      if(geo)status=nf_close(ncidl(1))
      write(*,*)'Results are in ',trim(outname)
      if(geo.and.ibl.eq.0)deallocate(zl1,lcind)
      stop
1     write(*,*)'Lat Lon file ',trim(coo_file),' not found'
      write(*,*)'Check setup file, line 2.'
      stop
4     write(*,*)'Can not open ',trim(fname)
      stop 
      end
