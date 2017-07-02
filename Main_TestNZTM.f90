Program Main_NZTM
    use MOD_NZTM
    implicit none
    
    double precision:: east, north, lat, lon, east1, north1;
    
    east  = 1576041.150  
    north = 6188574.240
    
    write(*, *) 'The east/nonth is:   ',  east, north
    
    call nztm_geod( north, east, lat, lon )
    write(*, *) 'The lat/lon is:      ', lat, lon
    
    call geod_nztm( lat, lon, north1, east1 ) 
    write(*, *) 'The east1/nonth1 is: ', east1, north1
    
    write(*,*) 'Difference:           ', int( (east1-east)*1000) /1000.0,  int( (north1-north)*1000) /1000.0;
    
    write(*,*) 'Done'
End Program 