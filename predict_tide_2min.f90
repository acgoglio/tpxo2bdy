!=====================================================================
! AUTHORS:
!  AC Goglio CMCC Bologna
!  from
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
      program predict_tide_da
!cc 
!cc   LANA, 2020 remake (for netcdf files)
!cc   needs min RAM and time to extract/predict by directly accessing
!cc   4 model nodes corresponding to given lat/lon cell 
!
!cc   modified March 2011 to optimize for obtaining
!cc   time series at open boundaries
!cc
!cc   reads OTIS netcdf model file
!cc   (elevations OR transports), reads a list of locations,
!cc   reads list of times  
!cc   and outputs ASCII file with the tidal predictions of tidal 
!cc   elevations/transports/currents at the locations and times
!cc   
      implicit none
      include 'netcdf.inc'
      include 'constit.h'
      complex, allocatable:: z1(:),dtmp(:,:)
      complex, allocatable:: zl1(:)
      complex, allocatable:: u1(:),v1(:)
      complex d1
      real, allocatable:: lat(:),lon(:),depth(:,:),x(:),y(:),&
                          lon0(:),zpred(:),upred(:),vpred(:),&
                          z_field(:,:),tides_z(:,:,:)!,bathy(:),mbathy(:,:)
      integer, allocatable:: tmask(:,:,:,:),mask(:),x_idx(:),y_idx(:),&
                             time_mjs(:)
      real xt,yt
      real*8, allocatable:: time_mjd(:) 
      real th_lim(2),ph_lim(2),dum,lth_lim(2),lph_lim(2)
      integer, allocatable:: cind(:),lcind(:),ccind(:),mz(:,:)
!
      character*4 c_id(ncmx),c_id_mod(ncmx),lc_id(ncmx),tcon(ncmx)
      character*80 modname,lltname,outname,ctmp,lname
      character*80 hname,uname,gname,fname,ErrSt
      character*2000 fmt
      character*80 rmCom
      character*1 zuv,c1,c2
      character*80 xy_ll_sub,arg,mesh_mask
      character*10 cdate
      character*8 ctime,day_date
      character*10 deblank 
      logical APRI,geo,interp_micon,ll_km
      integer ncon,nc,n,m,ndat,i,j,k,k1,ierr,ierr1,ic,n0,m0,it
      integer ncl,nl,ml,nmod,imod,ibl,ntime,idum,mjd,julian
      integer yyyy1,mm1,dd1,iargc,narg
      integer nca,l
      integer*4 status,ncid(ncmx),ncidl(ncmx),ncid_new(ncmx)
      integer, allocatable:: yyyy(:),mm(:),dd(:),hh(:),&
                             mi(:),ss(:),coo_array(:)!,time_array(:)
      integer dimid,len_coo,len_lon,len_lat,row,col,len_lev, &
              coo_idx, varid,dimid_lat,dimid_lon,dimid_field,dimid_coo,&
              varid_coo,varid_lat,varid_lon,varid_field,varid_time,&
              dimid_time,varid_navlat,varid_navlon,dimid_vtime !len_mask,
      real, allocatable:: nav_lon(:,:), nav_lat(:,:)!!!,time_array(:)
      !!!character, allocatable:: time_array(:)
      character*4 syyyy !,udm_time
      character*2 smm,sdd,smi,shh,sss
      character*19, allocatable:: time_array(:)
! Loop to read time --> Loop to compute time from ini date given as line arg
      ll_km=.false.
      WRITE(6,*) "Loop on date and times"
      narg=iargc()
      WRITE(6,*) "narg",narg
      ! OLD: This loop counts the num of lines to determine the num of time
      ! steps and calls the routine to read time values
      ! NEW: the num of time steps is fixed=24*30 (dt=2min in 1 day), the
      ! date is fixed from line arg and the hh mm ss runs on dt= 2 minutes
      !
      ! If there is the arg (should be the date in format: yyyymmdd)
      if(narg.eq.3)then
        ! Read the args, namely the date and mesh_mask.nc file  
        call getarg(1,day_date)
        WRITE(6,*) "day_date",day_date
        call getarg(2,mesh_mask)
        WRITE(6,*) "mesh_mask file: ",mesh_mask
        call getarg(3,outname)
        WRITE(6,*) "Out file: ",outname
        ! TIME 
        ! Set the number of time steps = 24*30 (minuts per day/2)
        ntime=720
        ! Allocate the date time values
        allocate(yyyy(ntime),mm(ntime),dd(ntime), &
               hh(ntime),mi(ntime),ss(ntime),time_mjd(ntime),time_mjs(ntime))
        ! Compute single date time values
        call read_time(ntime,day_date,yyyy,mm,dd,hh,mi,ss,&
             time_mjd,time_mjs)
        !do it=1,ntime
        !   write(6,*)"DATE,TIME",yyyy(it),mm(it),dd(it),hh(it),mi(it),ss(it)
        !enddo
      else
       WRITE(6,*) "Ini date, mesh_mask and outfile &
                   MUST be given as arguments!"
       stop
      endif       
      WRITE(6,*) "End loop on date and times"
!
      ibl=0
      lname='DATA/load_file.nc'
      ! Set values previously given in setup.inp file
      call rd_inp(modname,zuv,c_id,ncon,APRI,geo, &
                  interp_micon)
      ! ASCII OUTFILE:
      !! Compose the outname string and open the new file
      !outname='tide_tpxo9_'//day_date
      !outname=trim(outname)
      !WRITE(6,*) "Outname: ",outname
      !open(unit=11,file=outname,status='unknown')
      !
      ! NETCDF OUTFILE: the file is built when time lat and lon dims 
      ! have been defined: look forward!
!      
      ! Set mesh_mask file name from 2^ arg
      lltname=mesh_mask
!
      ! Read infos on tpxo grid, etc
      call rd_mod_file(modname,hname,uname,gname,xy_ll_sub,nca,c_id_mod)
      write(*,*)
      ! Print on stdout setup infos
      write(*,*)'Lat/Lon file:',trim(lltname)
      if(ncon.gt.0)write(*,*)'Constituents to include: ',c_id(1:ncon)
      if(zuv.eq.'z')then
       if(geo)then
         write(*,*)'Predict GEOCENTRIC tide'
       else
         write(*,*)'Predict OCEAN tide'
       endif
      endif
      if(interp_micon)write(*,*)'Interpolate minor constituents'
!
      ! Infos on tpxo9 input files..
      call rd_mod_header_nc(modname,zuv,n,m,th_lim,ph_lim,nc,c_id_mod,&
                          xy_ll_sub)
      write(*,*)'Model:        ',trim(modname(12:80))
      write(11,'(60a1)')('-',i=1,60)
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
                  ' in the OTPS'
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
       else
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
      if(zuv.eq.'z')then
        write(*,*)'Predict elevations (m)'
      else
        write(*,*)'Predict transport (m^2/s) and currents (cm/s)'
      endif
!
      k1=1
      tcon=''  
      if(ncon.eq.0)then
       ibl=1
       ncon=nc
       c_id=c_id_mod
      else
! check if all required constituents are in the model
       do ic=1,ncon
        do k=1,nc
         if(c_id(ic).eq.c_id_mod(k))then
          tcon(k1)=c_id(ic)
          k1=k1+1
          go to 14
         endif
        enddo
        write(*,*)'Constituent ',c_id(ic), ' is NOT in the model'
14      continue
       enddo 
       ncon=k1-1
       c_id=tcon
      endif
 
      write(*,*)'Constituents to include: ',c_id(1:ncon)
      write(11,*)'Constituents included: ',c_id(1:ncon)
!
      allocate(cind(ncon),ccind(ncon))
      call def_con_ind(c_id,ncon,c_id_mod,nc,cind)
! find corresponding indices in constit.h
      call def_cid(ncon,c_id,ccind)
! Read latitudes and longitudes from mesh_mask.nc
      write(6,*)'Open and read NEMO mesh_mask.nc file: ',lltname
      status=nf_open(trim(lltname),nf_nowrite,ncid)
      if(status.eq.0)then
       ! Read longitude len
       status=nf_inq_dimid(ncid,'x',dimid)
       status=nf_inq_dimlen(ncid,dimid,len_coo)
       len_lon=len_coo
       ! Read latitude len
       status=nf_inq_dimid(ncid,'y',dimid)
       status=nf_inq_dimlen(ncid,dimid,len_coo)
       len_lat=len_coo
       ! Read vertical level num 4 tmask
       status=nf_inq_dimid(ncid,'z',dimid)
       status=nf_inq_dimlen(ncid,dimid,len_lev)
      else
       write(6,*)'ERROR: Something wrong with the mesh_mask file..', status
       stop
      endif
      ! Grid points number
      ndat=len_lon*len_lat
      write(6,*)'Grid points number: ',len_lon,' x ',len_lat,' = ',ndat
      write(6,*)'Vertical levels: ',len_lev
      ! Read coo values from mesh_mask.nc
      allocate(lat(ndat),lon(ndat),nav_lon(len_lon,len_lat),&
              nav_lat(len_lon,len_lat),lon0(ndat),coo_array(ndat),& 
              x_idx(len_lon),y_idx(len_lat),&
              mask(ndat),tmask(len_lon,len_lat,len_lev,1))
              !bathy(ndat),mbathy(len_lon,len_lat))
      status=nf_inq_varid(ncid,'nav_lon',varid)
      status=nf_get_var(ncid,varid,nav_lon)
      status=nf_inq_varid(ncid,'nav_lat',varid)
      status=nf_get_var(ncid,varid,nav_lat)
      !ErrSt=nf_strerror(status)
      !status=nf_inq_varid(ncid,'x',varid)
      !status=nf_get_var(ncid,varid,x_idx)
      !status=nf_inq_varid(ncid,'y',varid)
      !status=nf_get_var(ncid,varid,y_idx)
      status=nf_inq_varid(ncid,'tmask',varid)
      !ErrSt=nf_strerror(status)
      status=nf_get_var_int(ncid,varid,tmask)
      !status=nf_inq_varid(ncid,'mbathy',varid)
      !status=nf_get_var(ncid,varid,bathy)
      !ndat=10 ! TMP
      ! 
      coo_idx = 1
      do col = 1, len_lat
       do row = 1, len_lon
         lon(coo_idx)=nav_lon(row,col)
         lat(coo_idx)=nav_lat(row,col)
         !bathy(coo_idx)=mbathy(row,col)
         mask(coo_idx)=tmask(row,col,1,1)
         !write(6,*)'tMask: ',mask(coo_idx)
         coo_array(coo_idx)=coo_idx
         coo_idx=coo_idx+1
         !if (col.eq.1) then
         !   x_idx(row)=row
         !endif
         !   y_idx(col)=col
       enddo
      enddo
      write(6,*)'coo tot: ',coo_idx-1
      ! Redefn of coordinates..
      write(6,*)'Start interpolation..  '
      if(zuv.eq.'z')then
       allocate(zpred(ntime)) ! TMP ndat
      else
       allocate(upred(ndat),vpred(ndat))
      endif
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
!     
      if(zuv.eq.'z')then      
       allocate(z1(ncon))
       if(geo) then ! NOT our case..
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
      else ! OUR CASE
       allocate(u1(ncon),v1(ncon))
      endif
!
       allocate(depth(n,m),dtmp(n,m),mz(n,m))
       write(*,'(a,$)')'Reading grid file...'
       call rd_grd_nc(gname,n,m,depth,mz)
       dtmp=depth
       deallocate(depth)
       write(*,*)'done'
!     ASCII OUTFILE HEADER:
      !ctmp='    Lat       Lon        mm.dd.yyyy hh:mm:ss'
      !write(11,*)'' 
      !if(zuv.eq.'z')then
      !  write(11,*)trim(ctmp),'     z(m)   Depth(m)'
      !else
      !  write(11,*)trim(ctmp),'   U(m^2/s)  V(m^2/s)',&
      !                        '   u(cm/s)   v(cm/s) Depth(m)'
      !endif
      !write(11,*)''
      ! NETCDF OUTFILE: create the file with lat,lon and time
      write(6,*)'Open outfile: ',outname
      ! Open new file in writing mode
      status=nf_create(trim(outname),nf_netcdf4,ncid_new) !nf_write
      ! Set the dimensions
      !status = nf_def_dim(ncid_new, 'coo', ndat, dimid_coo)
      status = nf_def_dim(ncid_new, 'x', len_lon, dimid_lon)
      status = nf_def_dim(ncid_new, 'y', len_lat, dimid_lat)
      status = nf_def_dim(ncid_new, 'time_counter', nf_unlimited, dimid_time) !ntime
      write(6,*)'Time status def dim.. ',status
      ! Set vars
      status = nf_def_var(ncid_new, 'time_counter', nf_int, 1, dimid_time, &
               varid_time) !nf_char
      write(6,*)'Time status def var.. ',status
      status = nf_put_att_text(ncid_new, varid_time,'units',33,&
               'seconds since 1994-12-24 06:28:16') ! 1994-12-24 06:28:16 1858-11-17
      write(6,*)'Time status put att.. ',status
     ! status = nf_def_var(ncid_new, 'coo', nf_int, 1,dimid_coo, &
     !          varid_coo)
     ! status = nf_def_var(ncid_new, 'longitude', nf_float,1,dimid_coo, &
     !          varid_lon)
     ! status = nf_def_var(ncid_new, 'latitude', nf_float,1,dimid_coo, &
     !          varid_lat)
      status = nf_def_var(ncid_new, 'nav_lat', nf_float,2,[dimid_lon,&
               dimid_lat],varid_navlat)
      status = nf_def_var(ncid_new, 'nav_lon', nf_float,2,[dimid_lon, &
               dimid_lat],varid_navlon)
      status = nf_def_var(ncid_new, 'tide_z', nf_float,3,[dimid_lon,&
               dimid_lat,dimid_time], varid_field)
      ! Compression: (WARNING: it increases the run time a lot!!)
      !status = nf_def_var_deflate(ncid_new,varid_field,0,1,1)
      write(6,*)'Time status field def var.. ',status
      ! Compress field 
      !status = nf_def_var_chunking(ncid_new, varid_field, NF_CHUNKED, &
      !         [10,101])
      ! Add attributes
      !status = nf_put_att(ncid_new, NF_GLOBAL, 'note', 'tide z tpxo9')
      !udm_time='minutes since 1858-11-17'
      !status = nf_put_att_text(ncid_new, dimid_time,'units',25,&
      !         'minutes since 1858-11-17')
      !write(6,*)'Att status ',status
      !status = nf_put_att(ncid_new, varid_lat, 'units', 'degree_north')
      !status = nf_put_att(ncid_new, varid_field, '_FillValue', 0)
      !
!
      c1=zuv ! since interp change zuv (U->u, V->v)
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
! Interpolation from tpxo9 to given lat/lon points and time evolution
      allocate(z_field(ntime,ndat),time_array(ntime)) ! BOH time_array(ndat)
      do k=1,ndat
       !write(*,*)'IDX: ',k,mask(k)
       if( mask(k).ne.0 ) then
        xt=lon(k)
        yt=lat(k)
        if(trim(xy_ll_sub).ne.'') then
         xt=x(k)
         yt=y(k)
        endif 
        if(zuv.eq.'z') then ! Our case
          if(ll_km)z1(1)=-1
          !write(*,*)'Start interpolation..'
          call interp_da_nc(ncid,n,m,th_lim,ph_lim, &
                         yt,xt,z1,ncon,cind,ierr,c1,nca)
          !write(*,*)'End interpolation..'
         else ! NOT our case
          if(ll_km)u1(1)=-1
          if(ll_km)v1(1)=-1
          c2='u'                           
          if(c1.eq.'U'.or.c1.eq.'V')c2='U'
          call interp_da_nc(ncid,n,m,th_lim,ph_lim, &
                        yt,xt,u1,ncon,cind,ierr,c2,nca)     
          c2='v'                           
          if(c1.eq.'U'.or.c1.eq.'V')c2='V'
          call interp_da_nc(ncid,n,m,th_lim,ph_lim, &
                        yt,xt,v1,ncon,cind,ierr,c2,nca)
       endif
       endif
       if(ierr.eq.0 .and. mask(k).ne.0 ) then !mask(k).eq.1 bathy(k).eq.0 
          !write(*,*)'Inside and on sea'
          !write(*,*)'PROVA IN:',mask(k) !bathy(k)
          if(ll_km)d1=-1
          call interp(dtmp,1,n,m,mz,th_lim,ph_lim, &
                       yt,xt,d1,ierr1,'z')
         if(zuv.eq.'z'.and.geo)then
          call interp_da_nc(ncidl,nl,ml,lth_lim,lph_lim, &
                lat(k),lon(k),zl1,ncon,lcind,ierr1,'z',0)
          z1=z1+zl1    ! apply load correction to get geocentric tide
         endif  
! predict tide 
! OB usage March 2011 (the other case has been removed)
         !do it=1,ntime
          if(zuv.eq.'z')then
           ! This is our case:
           !write(*,*)'Start time evolution..'
           call ptide(z1,c_id,ncon,ccind,lat(k),time_mjd,ntime,&
                     interp_micon,zpred,)
           ! Write in field!
           do it=1,ntime
            z_field(it,k)=zpred(it)
           enddo
           !write(*,*)'End time evolution..',z_field(k,it)
          !else
          ! call ptide(u1,c_id,ncon,ccind,lat(k),time_mjd(it),1,&
          !           interp_micon,upred(k))
          ! call ptide(v1,c_id,ncon,ccind,lat(k),time_mjd(it),1,&
          !          interp_micon,vpred(k))
          endif
          ! WRITE IN ASCII FILE
          !write(*,*)'Start writing values..'
          !write(cdate,'(i2,a1,i2,a1,i2)')hh(it),':',mi(it),':',ss(it)
          !ctime=deblank(cdate)
          !write(cdate,'(i2,a1,i2,a1,i4)')mm(it),'.',dd(it),'.',yyyy(it)
          !cdate=deblank(cdate)
          !if(it.eq.1)then
           !write(11,'(1x,f10.4,f10.4)')lat(k),lon0(k) 
          !endif
          !if(zuv.eq.'z')then ! OUR case
           !write(11,'(26x,a10,1x,a8,f10.3,f10.3)')&
                 !cdate,ctime,zpred(k),real(d1)
          !else ! NOT our case
           !write(11,'(26x,a10,1x,a8,5(f10.3))')&
                 !cdate,ctime,upred(k),vpred(k),&
            !upred(k)/real(d1)*100,vpred(k)/real(d1)*100,real(d1)
          !endif
         !enddo
       else 
          !write(*,*)'PROVA OUT:',mask(k) !bathy(k)
          ! WRITE IN ASCII FILE
          !write(11,'(1x,f10.4,f10.4,a)')lat(k),lon(k),&
          !'***** Site is out of model grid OR land *****'
          do it=1,ntime
           z_field(it,k)=0.0
          enddo
       endif  
      enddo
      ! Conversion from coo to nav_lat nav_lon
      write(*,*)'2D Conversion..'
      allocate(tides_z(len_lon,len_lat,ntime))
      do it=1,ntime
       coo_idx = 1
       do col = 1, len_lat
        do row = 1, len_lon
         tides_z(row,col,it)=z_field(it,coo_idx)
         coo_idx=coo_idx+1
        enddo
       enddo
      enddo
      ! WRITE IN THE NETCDF
      write(*,*)'Start netcdf writing..'
      ! time in char (if needed!)
      !write (syyyy,'(i4.4)') yyyy(1)
      !write (smm,'(i2.2)') mm(1)
      !write (sdd,'(i2.2)') dd(1)
      !do it=1,ntime
      !!   write(syyyy,'(i4)') yyyy(it)
      !!   write (smm,'(i2)') mm(it)
      !!   write (sdd,'(i2)') dd(it)
      !    write (shh,'(i2.2)') hh(it)
      !    write (smi,'(i2.2)') mi(it)
      !    write (sss,'(i2.2)') ss(it)
      !    time_array(it)=syyyy//'-'//smm//'-'//sdd//' '//shh//':'//smi//':'//sss !hh
      !enddo
      status = nf_enddef(ncid_new)
      status = nf_put_var(ncid_new, varid_navlat,nav_lat)
      status = nf_put_var(ncid_new, varid_navlon,nav_lon)
      status = nf_put_vara(ncid_new, varid_time,1,ntime,time_mjs)
      status = nf_put_vara(ncid_new,varid_field,[1,1,1],[len_lon,len_lat,ntime],tides_z)
      ! Write records with loop
      !do it=1,ntime
      !  status = nf_put_vara(ncid_new, varid_time,it,ntime,time_mjs(it)) ! time_array
      !  write(6,*)'Time status put var..',status,time_mjs(1),time_mjs(720)
      !  status = nf_put_vara(ncid_new,varid_field,&
      !           [1,1,it],[len_lon,len_lat,ntime],tides_z)
      !  write(6,*)'Time status field put var..',status
      !enddo
      ! Close the outfile
      status = nf_close(ncid_new)
      !
      deallocate(cind,ccind,lat,lon,lon0,dtmp,mz)
      if(zuv.eq.'z')then
       deallocate(z1,zpred)
      else
       deallocate(u1,v1,upred,vpred)
      endif
      if(trim(xy_ll_sub).ne.'')deallocate(x,y)
      if(zuv.eq.'z'.and.ibl.eq.1.and.geo)then
         ncon=0
         deallocate(zl1,lcind)
      endif
      close(11)
      close(12)
      status=nf_close(ncid(1))
      if(nca.ne.0.and.ncon.gt.1)then
       do ic=2,ncon
        status=nf_close(ncid(ic))
       enddo
      endif
      if(geo)status=nf_close(ncidl(1))
      write(*,*)'Results are in ',trim(outname)
      stop
1     write(*,*)'Lat lon file ',trim(lltname),' not found'
      write(*,*)'Check setup file, line 2.'
      stop
4     write(*,*)'Grid file ',trim(gname),' not found'
      write(*,*)'Check file ',trim(modname),', line 3'
      stop
6     write(*,*)'File ',trim(lltname),' not found'
      write(*,*)'Check setup file, line 7.'
      stop
16    write(*,*)'File ',trim(day_date),' not found'
      write(*,*)'Check spelling in command line'
      call usage()
      stop
11    write(*,*)'File ''model.list'' was NOT found...'
      write(*,*)'TO CREATE please do:'
      write(*,*)'ls -1 DATA/Model_*>model.list'
      stop
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine usage()
      write(*,*)'Usage:'
      write(*,*)'predict_tide [-t<time_file>]<setup.inp'
      write(*,*)'Default: lat_lon_time file is used for input'
      write(*,*)'          if option -t is given, then'
      write(*,*)'          lats/lons are read from lat_lon file'
      write(*,*)'          set in setup.inp, snd times are read'
      write(*,*)'          from <time_file>.'
      write(*,*)'          Use, when need output for time series'
      write(*,*)'          for open boundaries, i.e. times are'
      write(*,*)'          the same in all nodes'
      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine read_time(ntime,day_date,yyyy,mm,dd,hh,&
                           mi,ss,time_mjd,time_mjs)
      implicit none
      integer k,k_hh,k_mi,julian,mjd,mm1,dd1,yyyy1
      integer ntime,yyyy(ntime),mm(ntime),dd(ntime),yyyy_o,mm_o,dd_o
      integer hh(ntime),mi(ntime),ss(ntime),time_mjs(ntime)
      real*8 time_mjd(ntime)
      character*10 cdate,deblank
      character*8 day_date
!
      ! Loop to set date and compute time
      ! Compute date arrays
      read(day_date(1:4),'(i4)') yyyy_o
      read(day_date(5:6),'(i2)') mm_o
      read(day_date(7:8),'(i2)') dd_o
      WRITE(6,*) "year ",yyyy_o
      WRITE(6,*) "month ",mm_o
      WRITE(6,*) "day ",dd_o
      ! convert to mjd
      call date_mjd(mm_o,dd_o,yyyy_o,mjd)
      WRITE(6,*) "mjd ",mjd
      ! check if exists such a date
      julian=mjd+2400001
      WRITE(6,*) "julian ",julian
      call CALDAT (julian,MM1,DD1,YYYY1)
      if(mm_o.ne.mm1.or.dd_o.ne.dd1.or.yyyy_o.ne.yyyy1)then
       write(cdate,'(i2,a1,i2,a1,i4)')mm_o,'.',dd_o,'.',yyyy_o
       cdate=deblank(cdate)
       write(*,*)'Wrong date in line arg:',cdate
       stop
      endif
      yyyy=yyyy_o
      mm=mm_o
      dd=dd_o
      ss=0
      ! Loop to build hh and min arrays 
      k=1
      k_hh=0
      do while (k_hh<24)
       k_mi=0 
       do while (k_mi<60)
        hh(k)=k_hh
        mi(k)=k_mi
        !WRITE(6,*) "k hh min ",k,hh(k),mi(k)
        ! Julian days
        time_mjd(k)=dble(mjd)+dble(hh(k))/24.D0+ &
                    dble(mi(k))/(24.D0*60.D0)  + &
                    dble(ss(k))/(24.D0*60.D0*60.D0)
        ! Julian seconds
        time_mjs(k)=mjd*(24*60*60)+ &
                    hh(k)*(60*60)+ &
                    mi(k)*(60)+ &
                    ss(k)*1
        k=k+1
        k_mi=k_mi+2
       enddo
       k_hh=k_hh+1
      enddo
     return
     end
