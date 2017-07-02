# Projection-Conversion-Between-NZTM-and-LatLon
Fortran 90 module for transforming the projections between NZTM and LatLon. It is simply translated from C source code from LINZ at http://www.linz.govt.nz/geodetic/software-downloads#nztm2000

Only provide two public functions
public nztm_geod, geod_nztm 

nztm_geod                                                                                                                                       !*   Routine to convert from Tranverse Mercator to latitude and longitude. */
    !*   Method based on Redfearn's formulation as expressed in GDA technical  */
    !*   manual at http://www.anzlic.org.au/icsm/gdatm/index.html  
    
geod_nztm
    !*   Routine to convert from latitude and longitude to Transverse Mercator.*/
    !*   Method based on Redfearn's formulation as expressed in GDA technical  */
    !*   manual at http://www.anzlic.org.au/icsm/gdatm/index.html              */
    !*   Loosely based on FORTRAN source code by J.Hannah and A.Broadhurst.    */
    
Example

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
