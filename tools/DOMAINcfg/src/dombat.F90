MODULE dombat

   USE dom_oce           ! ocean domain
!   USE closea            ! closed seas
   !
   USE in_out_manager    ! I/O manager
   USE iom               ! I/O library
   USE lbclnk            ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp           ! distributed memory computing library
   USE agrif_modutil
   USE bilinear_interp

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dom_bat        ! called by dom_zgr.F90

CONTAINS

   SUBROUTINE dom_bat

      INTEGER :: inum, isize, jsize, id, ji, jj
      INTEGER :: tabdim1, tabdim2, nxhr, nyhr, nxyhr
      INTEGER, DIMENSION(2) :: ddims
      INTEGER, DIMENSION(3) :: status
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: trouble_points ,  vardep,mask_oce
      INTEGER :: iimin,iimax,jjmin,jjmax
      INTEGER, DIMENSION(1) :: i_min,i_max
      INTEGER, DIMENSION(1) ::j_min,j_max
      REAL(wp), DIMENSION(jpi,jpj) :: myglamf
      INTEGER,DIMENSION(:)  ,POINTER     ::   src_add,dst_add => NULL()
      REAL(wp), DIMENSION(:)  ,ALLOCATABLE ::   vardep1d, lon_new1D,lat_new1D 
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: bathy_new, lat_new, lon_new, bathy_test
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: coarselon, coarselat, coarsebathy
      REAL(wp) :: Cell_lonmin, Cell_lonmax, Cell_latmin, Cell_latmax
      LOGICAL :: identical_grids
      LOGICAL, DIMENSION(:,:), ALLOCATABLE     ::   masksrc
      REAL*8, DIMENSION(:,:),POINTER     ::   matrix,interpdata => NULL() 
      LOGICAL :: lonlat_2D

      CHARACTER(32) :: bathyfile, bathyname, lonname,latname       

     lonlat_2D=.false.

      bathyfile=TRIM(cn_topo)
      bathyname=TRIM(cn_bath)
      lonname=TRIM(cn_lon)
      latname=TRIM(cn_lat)
   
      CALL iom_open( bathyfile, inum, lagrif=.FALSE. )
      id = iom_varid( inum, bathyname, ddims )
      
      status=-1
      ALLOCATE(lon_new  (ddims(1),ddims(2)), STAT=status(1)) 
      ALLOCATE(lat_new  (ddims(1),ddims(2)), STAT=status(2)) 
      ALLOCATE(bathy_new(ddims(1),ddims(2)), STAT=status(3)) 
      IF( sum(status) /= 0 )   CALL ctl_stop( 'STOP', 'dom_bat : unable to allocate arrays' )

      IF (lonlat_2D) THEN
        CALL iom_get  ( inum, jpdom_unknown, lonname, lon_new )
        CALL iom_get  ( inum, jpdom_unknown, latname, lat_new )
      ELSE
        ALLOCATE(lon_new1D(ddims(1)), lat_new1D(ddims(2)))
        CALL iom_get  ( inum, jpdom_unknown, lonname, lon_new1D )
        CALL iom_get  ( inum, jpdom_unknown, latname, lat_new1D )
        DO ji=1, ddims(1)
           lon_new(ji,:)=lon_new1D(ji)
        ENDDO              
        DO ji=1, ddims(2)
           lat_new(:,ji)=lat_new1D(ji)
        ENDDO              
      ENDIF  
      CALL iom_get  ( inum, jpdom_unknown, bathyname, bathy_new )
      CALL iom_close (inum)
      WHERE (bathy_new > 0.)        
        bathy_new=0.
      ENDWHERE 
      bathy_new=-bathy_new 

       ! Eventually add here a pre-selection of the area (coarselon/lat)
       i_min=10  ; j_min=10
       i_max= ddims(1)-10 ; j_max=ddims(2)-10 

       tabdim1 = i_max(1) -  i_min(1) + 1 
       tabdim2 = j_max(1) - j_min(1) + 1
       !

       ALLOCATE(coarselon(tabdim1,tabdim2))
       ALLOCATE(coarselat(tabdim1,tabdim2))
       ALLOCATE(coarsebathy(tabdim1,tabdim2))          
       !
       WHERE( lon_new < 0. )   lon_new = lon_new + 360.
       myglamf=glamf
       WHERE( myglamf < 0. )   myglamf = myglamf + 360.

       coarselat(:,:)   =  lat_new  (i_min(1):i_max(1), j_min(1):j_max(1))
       coarselon  (:,:) =  lon_new  (i_min(1):i_max(1), j_min(1):j_max(1))
       coarsebathy(:,:) =  bathy_new(i_min(1):i_max(1), j_min(1):j_max(1))

       IF( nn_interp == 0 .OR. nn_interp == 1 ) THEN   ! arithmetic or median averages
        !                                                            ! ----------------------------- 
        ALLOCATE(trouble_points(jpi,jpj))
        trouble_points(:,:) = 0
        !                       
        DO jj = 2, jpj
           DO ji = 2, jpi
              !       
              ! fine grid cell extension               
              Cell_lonmin = MIN(myglamf(ji-1,jj-1),myglamf(ji,jj-1),myglamf(ji,jj),myglamf(ji-1,jj))
              Cell_lonmax = MAX(myglamf(ji-1,jj-1),myglamf(ji,jj-1),myglamf(ji,jj),myglamf(ji-1,jj))
              Cell_latmin = MIN(gphif(ji-1,jj-1),gphif(ji,jj-1),gphif(ji,jj),gphif(ji-1,jj))
              Cell_latmax = MAX(gphif(ji-1,jj-1),gphif(ji,jj-1),gphif(ji,jj),gphif(ji-1,jj)) 
              !               
              ! look for points in G0 (bathy dataset) contained in the fine grid cells  
              iimin = 1
              DO WHILE( coarselon(iimin,1) < Cell_lonmin ) 
                 iimin = iimin + 1
              ENDDO
              !               
              jjmin = 1
              DO WHILE( coarselat(iimin,jjmin) < Cell_latmin ) 
                 jjmin = jjmin + 1
              ENDDO
             !                
              iimax = iimin 
              DO WHILE( coarselon(iimax,1) <= Cell_lonmax ) 
                 iimax = iimax + 1
                 iimax = MIN( iimax,SIZE(coarsebathy,1))
              ENDDO
              !                               
              jjmax = jjmin 
              DO WHILE( coarselat(iimax,jjmax) <= Cell_latmax ) 
                 jjmax = jjmax + 1
                 jjmax = MIN( jjmax,SIZE(coarsebathy,2))
              ENDDO
              !
              IF( .NOT. Agrif_Root() ) THEN
                 iimax = iimax-1
                 jjmax = jjmax-1
              ELSE
                 iimax = MAX(iimin,iimax-1)
                 jjmax = MAX(jjmin,jjmax-1)
              ENDIF
              !               
              iimin = MAX( iimin,1 )
              jjmin = MAX( jjmin,1 )
              iimax = MIN( iimax,SIZE(coarsebathy,1))
              jjmax = MIN( jjmax,SIZE(coarsebathy,2))

              nxhr = iimax - iimin + 1
              nyhr = jjmax - jjmin + 1                    

                
              IF( nxhr == 0 .OR. nyhr == 0 ) THEN
                 !
                 trouble_points(ji,jj) = 1
                 !
              ELSE
                 !
                 ALLOCATE( vardep(nxhr,nyhr), mask_oce(nxhr,nyhr) )
                 vardep(:,:) = coarsebathy(iimin:iimax,jjmin:jjmax)
                 !
                 WHERE( vardep(:,:) .GT. 0. )   ;   mask_oce = 1 ;
                 ELSEWHERE                      ;   mask_oce = 0 ;
                 ENDWHERE
                 !
                 nxyhr = nxhr*nyhr
                 IF( SUM(mask_oce) < 0.5*(nxyhr) ) THEN ! if more than half of the points are on land then bathy fine = 0
                    bathy(ji,jj) = 0.
                 ELSE
                    IF( nn_interp == 0 ) THEN     ! Arithmetic average
                      bathy(ji,jj) = SUM( vardep(:,:) * mask_oce(:,:) ) / SUM( mask_oce(:,:) )
                    ELSE                                  ! Median average
                       ALLOCATE(vardep1d(nxyhr))
                       vardep1d = RESHAPE(vardep,(/ nxyhr /) )
                       !!CALL ssort(vardep1d,nxyhr)
                       CALL quicksort(vardep1d,1,nxyhr)
                       !
                       ! Calculate median
                       IF (MOD(nxyhr,2) .NE. 0) THEN
                          bathy(ji,jj) = vardep1d( nxyhr/2 + 1 )
                       ELSE
                          bathy(ji,jj) = 0.5 * ( vardep1d(nxyhr/2) + vardep1d(nxyhr/2+1) )
                       END IF
                       DEALLOCATE(vardep1d)   
                    ENDIF
                 ENDIF
                 DEALLOCATE (mask_oce,vardep)
                 !
              ENDIF
           ENDDO
        ENDDO

        IF( SUM( trouble_points ) > 0 ) THEN
           PRINT*,'too much empty cells, proceed to bilinear interpolation'
           nn_interp = 2
           stop
        ENDIF           
     ENDIF

#undef MYTEST
#ifdef MYTEST
         IF( nn_interp == 2) THEN                        ! Bilinear interpolation
        !                                                    ! ----------------------------- 
        identical_grids = .FALSE.

        IF( SIZE(coarselat,1) == jpi .AND. SIZE(coarselat,2) == jpj  .AND.   &
            SIZE(coarselon,1) == jpj .AND. SIZE(coarselon,2) == jpj ) THEN
           IF( MAXVAL( ABS(coarselat(:,:)- gphit(:,:)) ) < 0.0001 .AND.   &
               MAXVAL( ABS(coarselon(:,:)- glamt(:,:)) ) < 0.0001 ) THEN
              PRINT*,'same grid between ', cn_topo,' and child domain'    
              bathy = bathy_new
              identical_grids = .TRUE.                          
           ENDIF
        ENDIF

        IF( .NOT. identical_grids ) THEN 

           ALLOCATE(masksrc(SIZE(bathy_new,1),SIZE(bathy_new,2)))
           ALLOCATE(bathy_test(jpi,jpj))
           !
           !Where(G0%bathy_meter.le.0.00001) 
           !  masksrc = .false.
           !ElseWhere
              masksrc = .TRUE.
           !End where                       
           !            
           ! compute remapping matrix thanks to SCRIP package            
           CALL get_remap_matrix(coarselat,gphit,coarselon,glamt,masksrc,matrix,src_add,dst_add)
           CALL make_remap(bathy_new,bathy_test,jpi,jpj,matrix,src_add,dst_add)  
           !                                  
           bathy = bathy_test               
           !            
           DEALLOCATE(masksrc)
           DEALLOCATE(bathy_test) 

        ENDIF
        !            
     ENDIF
#endif 
   END SUBROUTINE dom_bat

END MODULE dombat
