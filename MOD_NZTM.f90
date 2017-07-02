Module MOD_NZTM
    implicit none
    ! NZGD Datum 2000 Parameters
    double precision:: NZTM_A, NZTM_RF
    parameter( NZTM_A=6378137, NZTM_RF=298.257222101)    

    ! Projection: Transverse Mercator Projection
    double precision:: NZTM_CM, NZTM_OLAT, NZTM_SF, NZTM_FE, NZTM_FN
    parameter(NZTM_CM = 173.0, NZTM_OLAT = 0.0, NZTM_SF= 0.9996, NZTM_FE=1600000.0, NZTM_FN=10000000.0)
    
    double precision:: PI, TWOPI, rad2deg
    Parameter(PI=3.1415926535898, TWOPI=2.0*PI, rad2deg=180.0/PI)
    
    Type tmprojection
	    double precision:: meridian;          !* Central meridian */
	    double precision:: scalef;            !* Scale factor */
	    double precision:: orglat;            !* Origin latitude */
	    double precision:: falsee;            !* False easting */
	    double precision:: falsen;            !* False northing */
	    double precision:: utom;              !* Unit to metre conversion */

	    double precision:: a, rf, f, e2, ep2; !* Ellipsoid parameters */
	    double precision:: om;                !* Intermediate calculation */
    End Type tmprojection

    
    private
    public nztm_geod, geod_nztm    
    
    Contains
    
    !************************************************************************************
    !Initiallize the TM structure  */

    Subroutine define_tmprojection( tm, a, rf, cm, sf, lto, fe, fn, utom ) 
        type(tmprojection) tm
        double precision:: a, rf, cm, sf, lto, fe, fn, utom    
        double precision:: f;

        tm%meridian = cm;
        tm%scalef = sf;
        tm%orglat = lto;
        tm%falsee = fe;
        tm%falsen = fn;
        tm%utom = utom;
        if( rf /= 0.0 ) then
            f = 1.0/rf; 
        else 
            f = 0.0;
        endif
        
        tm%a = a;
        tm%rf = rf;
        tm%f = f;
        tm%e2 = 2.0*f - f*f;
        tm%ep2 = tm%e2/( 1.0 - tm%e2 );

        tm%om = meridian_arc( tm, tm%orglat );
    End Subroutine define_tmprojection
    
    !************************************************************************************
    !***************************************************************************/
    !*                                                                         */
    !*  meridian_arc                                                           */
    !*                                                                         */
    !*  Returns the length of meridional arc (Helmert formula)                 */
    !*  Method based on Redfearn's formulation as expressed in GDA technical   */
    !*  manual at http://www.anzlic.org.au/icsm/gdatm/index.html               */
    !*                                                                         */
    !*  Parameters are                                                         */
    !*    projection                                                           */
    !*    latitude (radians)                                                   */
    !*                                                                         */
    !*  Return value is the arc length in metres                               */
    !*                                                                         */
    !***************************************************************************/


    Function meridian_arc( tm, lt ) 
        type(tmprojection) tm
        double precision:: lt, meridian_arc
        
        double precision:: e2, a;
        double precision:: e4, e6;
        double precision:: A0, A2, A4, A6;

        e2 = tm%e2;
        a  = tm%a;
        e4 = e2*e2;
        e6 = e4*e2;
 
        A0 = 1 - (e2/4.0) - (3.0*e4/64.0) - (5.0*e6/256.0);
        A2 = (3.0/8.0) * (e2+e4/4.0+15.0*e6/128.0);
        A4 = (15.0/256.0) * (e4 + 3.0*e6/4.0);
        A6 = 35.0*e6/3072.0;

        meridian_arc = a*(A0*lt-A2*sin(2*lt)+A4*sin(4*lt)-A6*sin(6*lt));
    
    End Function meridian_arc

    !*************************************************************************/
    !*                                                                       */
    !*   foot_point_lat                                                      */
    !*                                                                       */
    !*   Calculates the foot point latitude from the meridional arc          */
    !*   Method based on Redfearn's formulation as expressed in GDA technical*/
    !*   manual at http://www.anzlic.org.au/icsm/gdatm/index.html            */
    !*                                                                       */
    !*   Takes parameters                                                    */
    !*      tm definition (for scale factor)                                 */
    !*      meridional arc (metres)                                          */
    !*                                                                       */
    !*   Returns the foot point latitude (radians)                           */ /*                                                                       */
    !*************************************************************************/


    Function foot_point_lat( tm, m ) 
        type(tmprojection) tm
        double precision:: m, foot_point_lat
        
        double precision:: f, a;
        double precision:: n, n2, n3, n4, g, sig, phio;
 
        f  = tm%f
        a  = tm%a;
        n  = f/(2.0-f);
        n2 = n*n;
        n3 = n2*n;
        n4 = n2*n2;
 
        g = a*(1.0-n)*(1.0-n2)*(1+9.0*n2/4.0+225.0*n4/64.0);
        sig = m/g;
 
        phio = sig + (3.0*n/2.0 - 27.0*n3/32.0)*sin(2.0*sig) &
                        + (21.0*n2/16.0 - 55.0*n4/32.0)*sin(4.0*sig) &
                        + (151.0*n3/96.0) * sin(6.0*sig) &
                        + (1097.0*n4/512.0) * sin(8.0*sig);
 
        foot_point_lat = phio;
    End Function foot_point_lat
    
    !***************************************************************************/
    !*                                                                         */
    !*   tmgeod                                                                */
    !*                                                                         */
    !*   Routine to convert from Tranverse Mercator to latitude and longitude. */
    !*   Method based on Redfearn's formulation as expressed in GDA technical  */
    !*   manual at http://www.anzlic.org.au/icsm/gdatm/index.html              */
    !*                                                                         */
    !*   Takes parameters                                                      */
    !*      input easting (metres)                                             */
    !*      input northing (metres)                                            */
    !*      output latitude (radians)                                          */
    !*      output longitude (radians)                                         */
    !*                                                                         */
    !***************************************************************************/

    Subroutine tm_geod( tm, ce, cn, ln, lt ) 
        type(tmprojection) tm
        double precision:: ce, cn, ln, lt          
                           
        double precision:: fn, fe, sf, e2, a, cm, om, utom;
        double precision:: cn1;
        double precision:: fphi;
        double precision:: slt;
        double precision:: clt;
        double precision:: eslt;
        double precision:: eta;
        double precision:: rho;
        double precision:: psi;
        double precision:: E;
        double precision:: x, x2;
        double precision:: t, t2, t4;
        double precision:: trm1, trm2, trm3, trm4;
        
        fn   = tm%falsen;
        fe   = tm%falsee;
        sf   = tm%scalef;
        e2   = tm%e2;
        a    = tm%a;
        cm   = tm%meridian;
        om   = tm%om;
        utom = tm%utom;
 
        cn1  =  (cn - fn)*utom/sf + om;
        fphi = foot_point_lat(tm, cn1);
        slt = sin(fphi);
        clt = cos(fphi);
 
        eslt = (1.0-e2*slt*slt);
        eta = a/sqrt(eslt);
        rho = eta * (1.0-e2) / eslt;
        psi = eta/rho;
 
        E = (ce-fe)*utom;
        x = E/(eta*sf);
        x2 = x*x;
 
 
        t = slt/clt;
        t2 = t*t;
        t4 = t2*t2;
 
        trm1 = 1.0/2.0;
 
        trm2 = ((-4.0*psi +9.0*(1-t2))*psi+12.0*t2)/24.0;
 
        trm3 = ((((8.0*(11.0-24.0*t2)*psi-12.0*(21.0-71.0*t2))*psi &
                      +15.0*((15.0*t2-98.0)*t2+15))*psi &
                      +180.0*((-3.0*t2+5.0)*t2))*psi + 360.0*t4)/720.0;
 
        trm4 = (((1575.0*t2+4095.0)*t2+3633.0)*t2+1385.0)/40320.0;
 
        lt = fphi+(t*x*E/(sf*rho))*(((trm4*x2-trm3)*x2+trm2)*x2-trm1);
        !********************************************************************
        ! Modified by Chonghua, convert lat from rad to degree
        lt = lt*rad2deg
        !********************************************************************
 
        trm1 = 1.0;
 
        trm2 = (psi+2.0*t2)/6.0;
 
        trm3 = (((-4.0*(1.0-6.0*t2)*psi &
                   +(9.0-68.0*t2))*psi &
                   +72.0*t2)*psi &
                   +24.0*t4)/120.0;
 
        trm4 = (((720.0*t2+1320.0)*t2+662.0)*t2+61.0)/5040.0;
 
        ln = cm - (x/clt)*(((trm4*x2-trm3)*x2+trm2)*x2-trm1);
        !********************************************************************
        ! Modified by Chonghua, convert lon from rad to degree
        ln = ln*rad2deg
        !********************************************************************
        
    End Subroutine tm_geod
    
    
    !***************************************************************************/
    !*                                                                         */
    !*   geodtm                                                                */
    !*                                                                         */
    !*   Routine to convert from latitude and longitude to Transverse Mercator.*/
    !*   Method based on Redfearn's formulation as expressed in GDA technical  */
    !*   manual at http://www.anzlic.org.au/icsm/gdatm/index.html              */
    !*   Loosely based on FORTRAN source code by J.Hannah and A.Broadhurst.    */
    !*                                                                         */
    !*   Takes parameters                                                      */
    !*      input latitude (radians)                                           */
    !*      input longitude (radians)                                          */
    !*      output easting  (metres)                                           */
    !*      output northing (metres)                                           */
    !*                                                                         */
    !***************************************************************************/

    Subroutine geod_tm( tm, ln, lt, ce, cn) 
        type(tmprojection) tm
        double precision:: ln, lt, ce, cn
        
        double precision:: fn, fe, sf, e2, a, cm, om, utom;
        double precision:: dlon;
        double precision:: m;
        double precision:: slt;
        double precision:: eslt;
        double precision:: eta;
        double precision:: rho;
        double precision:: psi;
        double precision:: clt;
        double precision:: w, wc, wc2;
        double precision:: t, t2, t4, t6;
        double precision:: trm1, trm2, trm3, gce, trm4, gcn;
 
        ! modified by chonghuaconvert lat, lon from degree to rad
        ln   = ln/rad2deg;
        lt   = lt/rad2deg;
        !********************************************************
        
        fn   = tm%falsen;
        fe   = tm%falsee;
        sf   = tm%scalef;
        e2   = tm%e2;
        a    = tm%a;
        cm   = tm%meridian;
        om   = tm%om;
        utom = tm%utom;
        dlon =  ln - cm;
        
        do while ( dlon > PI ) 
            dlon = dlon -TWOPI;
        enddo
        
        do while ( dlon < -PI ) 
            dlon = dlon + TWOPI;
        enddo
 
        m = meridian_arc(tm,lt);
 
        slt = sin(lt);
 
        eslt = (1.0-e2*slt*slt);
        eta = a/sqrt(eslt);
        rho = eta * (1.0-e2) / eslt;
        psi = eta/rho;
 
        clt = cos(lt);
        w = dlon;
 
        wc = clt*w;
        wc2 = wc*wc;
 
        t = slt/clt;
        t2 = t*t;
        t4 = t2*t2;
        t6 = t2*t4;
 
        trm1 = (psi-t2)/6.0;
 
        trm2 = (((4.0*(1.0-6.0*t2)*psi &
                      + (1.0+8.0*t2))*psi &
                      - 2.0*t2)*psi+t4)/120.0;
 
        trm3 = (61 - 479.0*t2 + 179.0*t4 - t6)/5040.0;
 
        gce = (sf*eta*dlon*clt)*(((trm3*wc2+trm2)*wc2+trm1)*wc2+1.0);
        ce = gce/utom +fe;
 
        trm1 = 1.0/2.0;
 
        trm2 = ((4.0*psi+1)*psi-t2)/24.0;
 
        trm3 = ((((8.0*(11.0-24.0*t2)*psi &
                    -28.0*(1.0-6.0*t2))*psi &
                    +(1.0-32.0*t2))*psi &
                    -2.0*t2)*psi &
                    +t4)/720.0;
 
        trm4 = (1385.0-3111.0*t2+543.0*t4-t6)/40320.0;
 
        gcn = (eta*t)*((((trm4*wc2+trm3)*wc2+trm2)*wc2+trm1)*wc2);
        cn = (gcn+m-om)*sf/utom+fn;

    End Subroutine geod_tm
    
    !**************************************************************************************
    ! Initialize NZTM Projection
    Subroutine init_nztm_projection( nztm_projection )
        type(tmprojection) nztm_projection
     
        call define_tmprojection( nztm_projection, NZTM_A, NZTM_RF, NZTM_CM/rad2deg, &
                                  NZTM_SF, NZTM_OLAT/rad2deg, NZTM_FE, NZTM_FN, dble(1.0) );   
       
    End Subroutine init_nztm_projection      
    
    !************************************************************************************
    ! converts a northing north and easting east to a latitude lat and longitude lon
    Subroutine nztm_geod( north, east, lat, lon )
        double precision:: north, east
        double precision, intent(out):: lat, lon
    
        type(tmprojection) nztm 
        call init_nztm_projection(nztm);
        call tm_geod( nztm, east, north, lon, lat );
    
    End Subroutine nztm_geod
    
    !************************************************************************************
    ! converts a latitude lat and longitude lon to a northing north and easting east
    Subroutine geod_nztm( lat, lon, north, east ) 
        double precision:: lat, lon
        double precision, intent(out):: north, east
        
        type(tmprojection) nztm 
        
        call init_nztm_projection(nztm);
        call geod_tm( nztm, lon, lat, east, north );
    End Subroutine geod_nztm 

End Module MOD_NZTM
    