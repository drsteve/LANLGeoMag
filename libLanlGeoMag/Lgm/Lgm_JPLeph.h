/*
 * The JPL DEXXX ephemerides are stored as Chebychev position coefficients for each
 * component of each body.  The Chebychev coefficients for the planets represent 
 * the solar system barycentric positions of the centers of the planetary systems.
 * That is, the coordinates are in the ICRF (relative to the solar system barycenter).
 *
 * The solar system bodies stored are: Sun, Mercury, Venus, EarthMoon (the Earth-Moon
 * system barycenter), Moon (geocentric), Mars, Jupiter, Saturn, Uranus, Neptune,
 * Pluto.
 * 
 * There are three Cartesian components (x, y, z), for each of the 11 solar system bodies; 
 * there are two components for the 12th item, nutations: d(psi) and d(epsilon);
 * there are three components for the 13th item, librations: three Euler angles.
 * 
 * Planetary positions are stored in units of kilometers (TDB-compatible).
 * The nutations and librations are stored in units of radians.
 */

#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

/*
 *  Bitmap for solar system bodies, etc.
 */
#define LGM_DE_SUN            1
#define LGM_DE_EARTHMOON      2
#define LGM_DE_INNERPLANETS   4
#define LGM_DE_OUTERPLANETS   8
#define LGM_DE_LIBR_NUT       16

#define LGM_DE_ALLBODIES (LGM_DE_SUN | LGM_DE_EARTHMOON | LGM_DE_INNERPLANETS | LGM_DE_OUTERPLANETS | LGM_DE_LIBR_NUT)

//And object IDs for calculating positions
#define LGM_DE_EARTH          101
#define LGM_DE_MOON           102
#define LGM_DE_MERCURY        103
#define LGM_DE_VENUS          104
#define LGM_DE_MARS           105
#define LGM_DE_JUPITER        106
#define LGM_DE_SATURN         107
#define LGM_DE_URANUS         108
#define LGM_DE_NEPTUNE        109
#define LGM_DE_PLUTO          110
#define LGM_DE_LIBRATION      111
#define LGM_DE_NUTATION       112


/*
 *  Constants for DE421 - how many of these do we need???
 */
#define LGM_DE421_DENUM   421.0
#define LGM_DE421_LENUM   421.0
#define LGM_DE421_TDATEF   0.0
#define LGM_DE421_TDATEB   12008021118111700.0
#define LGM_DE421_CENTER   0.0
#define LGM_DE421_CLIGHT   299792.45799999998
#define LGM_DE421_AU   149597870.69962621
#define LGM_DE421_EMRAT   81.300569069915298
#define LGM_DE421_GM1   4.9125495718679402e-11
#define LGM_DE421_GM2   7.2434523326984407e-10
#define LGM_DE421_GMB   8.9970114082680488e-10
#define LGM_DE421_GM4   9.5495486956223901e-11
#define LGM_DE421_GM5   2.8253458408550499e-07
#define LGM_DE421_GM6   8.4597060733084774e-08
#define LGM_DE421_GM7   1.2920248257926499e-08
#define LGM_DE421_GM8   1.52435910924974e-08
#define LGM_DE421_GM9   2.17844105199052e-12
#define LGM_DE421_GMS   0.00029591220828559109
#define LGM_DE421_RAD1   2439.8762509925318
#define LGM_DE421_RAD2   6058.8491732307048
#define LGM_DE421_RAD4   3397.5149999999999
#define LGM_DE421_JDEPOC   2440400.5
#define LGM_DE421_X1   0.36176271710692681
#define LGM_DE421_Y1   -0.090781962796313481
#define LGM_DE421_Z1   -0.085714981975454393
#define LGM_DE421_XD1   0.003367493903543759
#define LGM_DE421_YD1   0.024894520087168508
#define LGM_DE421_ZD1   0.01294630135160613
#define LGM_DE421_X2   0.61275194206322925
#define LGM_DE421_Y2   -0.34836536774854598
#define LGM_DE421_Z2   -0.1952782855265712
#define LGM_DE421_XD2   0.01095206841767497
#define LGM_DE421_YD2   0.01561768426034194
#define LGM_DE421_ZD2   0.0063311057275951032
#define LGM_DE421_XB   0.1205174164491314
#define LGM_DE421_YB   -0.92583847312554313
#define LGM_DE421_ZB   -0.40154022623082858
#define LGM_DE421_XDB   0.016811268307633291
#define LGM_DE421_YDB   0.0017483092419644511
#define LGM_DE421_ZDB   0.00075820292874323567
#define LGM_DE421_X4   -0.11018607268761769
#define LGM_DE421_Y4   -1.327599448796144
#define LGM_DE421_Z4   -0.6058891411186913
#define LGM_DE421_XD4   0.014481653057587021
#define LGM_DE421_YD4   0.0002424630984941193
#define LGM_DE421_ZD4   -0.00028152069232697711
#define LGM_DE421_X5   -5.379706768297444
#define LGM_DE421_Y5   -0.8304814016092148
#define LGM_DE421_Z5   -0.22482876556029499
#define LGM_DE421_XD5   0.001092012625257437
#define LGM_DE421_YD5   -0.0065181166775549687
#define LGM_DE421_ZD5   -0.0028207825773521162
#define LGM_DE421_X6   7.894390677484485
#define LGM_DE421_Y6   4.596478075225197
#define LGM_DE421_Z6   1.558695818747458
#define LGM_DE421_XD6   -0.0032175565229011481
#define LGM_DE421_YD6   0.0043358103439020239
#define LGM_DE421_ZD6   0.0019286463028107301
#define LGM_DE421_X7   -18.265400343728221
#define LGM_DE421_Y7   -1.161957232884
#define LGM_DE421_Z7   -0.25010443219098077
#define LGM_DE421_XD7   0.00022119027863849819
#define LGM_DE421_YD7   -0.0037624755734531391
#define LGM_DE421_ZD7   -0.0016510142513223489
#define LGM_DE421_X8   -16.05504055116263
#define LGM_DE421_Y8   -23.942189592251651
#define LGM_DE421_Z8   -9.4001568033834051
#define LGM_DE421_XD8   0.0026427703618792432
#define LGM_DE421_YD8   -0.001498312896391299
#define LGM_DE421_ZD8   -0.00067904184056411464
#define LGM_DE421_X9   -30.483317419840731
#define LGM_DE421_Y9   -0.87242008708120944
#define LGM_DE421_Z9   8.9115845075577749
#define LGM_DE421_XD9   0.000322208527104429
#define LGM_DE421_YD9   -0.0031435756754003069
#define LGM_DE421_ZD9   -0.0010779506297515331
#define LGM_DE421_XM   -0.00080817735456250673
#define LGM_DE421_YM   -0.0019946299870205989
#define LGM_DE421_ZM   -0.0010872626812385891
#define LGM_DE421_XDM   0.00060108481586158369
#define LGM_DE421_YDM   -0.00016744547005835111
#define LGM_DE421_ZDM   -8.556214120977451e-05
#define LGM_DE421_XC   0.0
#define LGM_DE421_YC   0.0
#define LGM_DE421_ZC   0.0
#define LGM_DE421_XDC   0.0
#define LGM_DE421_YDC   0.0
#define LGM_DE421_ZDC   0.0
#define LGM_DE421_BETA   1.0
#define LGM_DE421_GAMMA   1.0
#define LGM_DE421_J2SUN   1.9999999999999999e-07
#define LGM_DE421_GDOT   0.0
#define LGM_DE421_MA0001   1.3863904478558459e-13
#define LGM_DE421_MA0002   2.9882165103302159e-14
#define LGM_DE421_MA0004   3.931009658107358e-14
#define LGM_DE421_MAD1   1.0933017992969569
#define LGM_DE421_MAD2   3.4524557832816911
#define LGM_DE421_MAD3   4.2214632589401244
#define LGM_DE421_RE   6378.1363000000001
#define LGM_DE421_ASUN   696000.0
#define LGM_DE421_PHI   0.0051281320587143629
#define LGM_DE421_THT   0.38239320052300668
#define LGM_DE421_PSI   1.2941680560570821
#define LGM_DE421_OMEGAX   4.5829226788604672e-05
#define LGM_DE421_OMEGAY   -2.1856340816763818e-06
#define LGM_DE421_OMEGAZ   0.2299448610003971
#define LGM_DE421_S31M   5.9008000000000003e-06
#define LGM_DE421_TAUE1   0.01114245096880191
#define LGM_DE421_TAUE2   0.0065742922459716352
#define LGM_DE421_ROTEX   0.0057810799467625162
#define LGM_DE421_ROTEY   -0.01691560136657734
#define LGM_DE421_DROTEX   0.00025662334534043939
#define LGM_DE421_COBLAT   0.00037977129369143867
#define LGM_DE421_EDOT   0.0
#define LGM_DE421_IFAC   0.00069999999999999999
#define LGM_DE421_KVC   1.4937602493056611e-08
#define LGM_DE421_PHIC   -0.0042595183
#define LGM_DE421_THTC   0.40884429999999999
#define LGM_DE421_PSIC   -1.7145090000000001
#define LGM_DE421_OMGCX   -0.006777016229159276
#define LGM_DE421_OMGCY   -0.0011497514682134349
#define LGM_DE421_OMGCZ   0.229470176943708
#define LGM_DE421_AM   1738.0
#define LGM_DE421_J2M   0.00020327325763707239
#define LGM_DE421_J3M   8.4047015259410008e-06
#define LGM_DE421_J4M   -9.6422859999999992e-06
#define LGM_DE421_C22M   2.2389767096524131e-05
#define LGM_DE421_C31M   2.8452435000000001e-05
#define LGM_DE421_C32M   4.8463872400090304e-06
#define LGM_DE421_C33M   1.6740475300391421e-06
#define LGM_DE421_S32M   1.6841984741476e-06
#define LGM_DE421_S33M   -2.4855259999999999e-07
#define LGM_DE421_C41M   -5.692687e-06
#define LGM_DE421_S41M   1.5743929999999999e-06
#define LGM_DE421_C42M   -1.5861999999999999e-06
#define LGM_DE421_S42M   -1.517312e-06
#define LGM_DE421_C43M   -8.1204099999999996e-08
#define LGM_DE421_S43M   -8.0279070000000011e-07
#define LGM_DE421_C44M   -1.273941e-07
#define LGM_DE421_S44M   8.3147499999999997e-08
#define LGM_DE421_LBET   0.00063100220253646294
#define LGM_DE421_LGAM   0.00022773053141991419
#define LGM_DE421_K2M   0.021633683863607409
#define LGM_DE421_TAUM   0.1078586819856588
#define LGM_DE421_AE   6378.1363000000001
#define LGM_DE421_J2E   0.0010826253049999999
#define LGM_DE421_J3E   -2.532474e-06
#define LGM_DE421_J4E   1.6199740000000001e-06
#define LGM_DE421_K2E0   0.33500000000000002
#define LGM_DE421_K2E1   0.32000000000000001
#define LGM_DE421_K2E2   0.32000000000000001
#define LGM_DE421_TAUE0   0.064000000000000001
#define LGM_DE421_DROTEY   -0.0011481721760368501
#define LGM_DE421_GMAST1   3.8038482424406552e-14
#define LGM_DE421_GMAST2   1.13994252599966e-14
#define LGM_DE421_GMAST3   3.1494923361568481e-15
#define LGM_DE421_PSIDOT   0.0
#define LGM_DE421_MGMIS   1.0
#define LGM_DE421_MA0007   1.7744824515429809e-15
#define LGM_DE421_MA0324   1.4733481315551009e-15
#define LGM_DE421_MA0003   3.4242783009416692e-15
#define LGM_DE421_MA0006   1.35001440499976e-15
#define LGM_DE421_MA0009   1.2642013509650081e-15
#define LGM_DE421_MA0010   1.1959347789583871e-14
#define LGM_DE421_MA0019   1.0333648795561431e-15
#define LGM_DE421_MA0020   6.4848099228059792e-16
#define LGM_DE421_MA0024   8.9752675707191536e-16
#define LGM_DE421_MA0031   2.5404535483183989e-15
#define LGM_DE421_MA0041   1.1754838470754729e-15
#define LGM_DE421_MA0052   3.0183251043572349e-15
#define LGM_DE421_MA0139   4.1915762334793281e-16
#define LGM_DE421_MA0354   7.2842240607496357e-16
#define LGM_DE421_MA0511   3.6522758570194068e-15
#define LGM_DE421_MA0532   1.97490211916245e-15
#define LGM_DE421_MA0654   1.9996156724272161e-16
#define LGM_DE421_MA0005   3.5471586289505639e-16
#define LGM_DE421_MA0008   5.2647085053381123e-16
#define LGM_DE421_MA0013   9.193237222276462e-16
#define LGM_DE421_MA0014   7.7599147020627208e-16
#define LGM_DE421_MA0015   3.6525305443719564e-15
#define LGM_DE421_MA0016   4.9792973121502139e-15
#define LGM_DE421_MA0018   5.9442605141587072e-16
#define LGM_DE421_MA0022   1.0943589036506291e-15
#define LGM_DE421_MA0023   2.8719336010791748e-16
#define LGM_DE421_MA0027   1.8778104806675771e-16
#define LGM_DE421_MA0029   2.0208476918505491e-15
#define LGM_DE421_MA0045   8.8528706142174071e-16
#define LGM_DE421_MA0051   3.2010806116771229e-16
#define LGM_DE421_MA0065   1.5477275183826419e-15
#define LGM_DE421_MA0078   1.8907462127462091e-16
#define LGM_DE421_MA0097   1.9819501612500869e-16
#define LGM_DE421_MA0105   1.96597317770212e-16
#define LGM_DE421_MA0111   2.5908997910519999e-16
#define LGM_DE421_MA0344   2.5315613274938212e-16
#define LGM_DE421_MA0372   7.9190973295434792e-16
#define LGM_DE421_MA0405   2.058483140216775e-16
#define LGM_DE421_MA0409   4.8270616906988074e-16
#define LGM_DE421_MA0451   1.3595913621623679e-15
#define LGM_DE421_MA0704   5.495015030752055e-15
#define LGM_DE421_MA0747   4.3595751009390861e-16
#define LGM_DE421_MA0011   7.9395248351137861e-16
#define LGM_DE421_MA0021   3.1048649761980128e-16
#define LGM_DE421_MA0025   8.9461791612460559e-17
#define LGM_DE421_MA0028   3.6782501044471532e-16
#define LGM_DE421_MA0030   2.1104943845115819e-16
#define LGM_DE421_MA0042   2.042153926450012e-16
#define LGM_DE421_MA0060   4.6675023611284533e-17
#define LGM_DE421_MA0063   2.2832139456143961e-16
#define LGM_DE421_MA0069   9.2403534023231557e-16
#define LGM_DE421_MA0094   9.2401948463490617e-16
#define LGM_DE421_MA0098   1.2283796755031901e-16
#define LGM_DE421_MA0135   1.7436068022199109e-16
#define LGM_DE421_MA0145   3.3672012925053062e-16
#define LGM_DE421_MA0187   2.335168388376332e-16
#define LGM_DE421_MA0192   2.3774305146738432e-16
#define LGM_DE421_MA0194   4.0556072782435618e-16
#define LGM_DE421_MA0216   6.6737353354914072e-16
#define LGM_DE421_MA0230   2.802342422426607e-16
#define LGM_DE421_MA0337   7.2719617012796851e-17
#define LGM_DE421_MA0419   2.2735474822040491e-16
#define LGM_DE421_MA0488   3.645968026955162e-16
#define LGM_DE421_MA0554   9.8655294326978136e-17
#define LGM_DE421_XS   0.0045025097686452508
#define LGM_DE421_YS   0.00076707783929749328
#define LGM_DE421_ZS   0.00026605823283250572
#define LGM_DE421_XDS   -3.517495982409655e-07
#define LGM_DE421_YDS   5.1776262656339178e-06
#define LGM_DE421_ZDS   2.2291017721979059e-06
#define LGM_DE421_jalpha   2414992.5
#define LGM_DE421_jomega   2524624.5
#define LGM_DE421_jdelta   32.0



typedef struct Lgm_JPLephemInfo {
    //set DE version
    int         DEnum;              //!< e.g. 421 for DE421 ephemerides

    //Sun
    int         getSun;
    int         SunAlloced;
    double      ***sun;             //!< [ SetIndex ][ AxisIndex ][ CoefficientIndex ]
    int         sun_nvals;
    int         sun_naxes;
    int         sun_ncoeffs;

    //Earth-Moon system
    double      earthshare;         //!< 1.0 / (1.0 + LGM_DEXXX_EMRAT)
    double      moonshare;          //!< EMRAT / (1.0 + LGM_DEXXX_EMRAT)
    int         getEarth;
    int         EarthMoonAlloced;
    double      ***earthmoon;    //!< [ SetIndex ][ AxisIndex ][ CoefficientIndex ]
    int         earthmoon_nvals;
    int         earthmoon_naxes;
    int         earthmoon_ncoeffs;
    double      ***moon_wrt_earth;  //!< [ SetIndex ][ AxisIndex ][ CoefficientIndex ]
    int         moon_wrt_earth_nvals;
    int         moon_wrt_earth_naxes;
    int         moon_wrt_earth_ncoeffs;

    //Planets
    int         getInnerPlanets;
    int         InnerPlanetsAlloced;
    double      ***mercury;         //!< [ SetIndex ][ AxisIndex ][ CoefficientIndex ]
    int         mercury_nvals;
    int         mercury_naxes;
    int         mercury_ncoeffs;
    double      ***venus;           //!< [ SetIndex ][ AxisIndex ][ CoefficientIndex ]
    int         venus_nvals;
    int         venus_naxes;
    int         venus_ncoeffs;
    double      ***mars;            //!< [ SetIndex ][ AxisIndex ][ CoefficientIndex ]
    int         mars_nvals;
    int         mars_naxes;
    int         mars_ncoeffs;
    int         getOuterPlanets;
    int         OuterPlanetsAlloced;
    double      ***jupiter;         //!< [ SetIndex ][ AxisIndex ][ CoefficientIndex ]
    int         jupiter_nvals;
    int         jupiter_naxes;
    int         jupiter_ncoeffs;
    double      ***saturn;          //!< [ SetIndex ][ AxisIndex ][ CoefficientIndex ]
    int         saturn_nvals;
    int         saturn_naxes;
    int         saturn_ncoeffs;
    double      ***uranus;          //!< [ SetIndex ][ AxisIndex ][ CoefficientIndex ]
    int         uranus_nvals;
    int         uranus_naxes;
    int         uranus_ncoeffs;
    double      ***neptune;         //!< [ SetIndex ][ AxisIndex ][ CoefficientIndex ]
    int         neptune_nvals;
    int         neptune_naxes;
    int         neptune_ncoeffs;
    double      ***pluto;           //!< [ SetIndex ][ AxisIndex ][ CoefficientIndex ]
    int         pluto_nvals;
    int         pluto_naxes;
    int         pluto_ncoeffs;

    //other datasets
    int         getLibrationNutation;
    int         LibrNutAlloced;
    double      ***nutation;         //!< [ SetIndex ][ AxisIndex ][ CoefficientIndex ]
    int         nutation_nvals;
    int         nutation_naxes;
    int         nutation_ncoeffs;
    double      ***libration;         //!< [ SetIndex ][ AxisIndex ][ CoefficientIndex ]
    int         libration_nvals;
    int         libration_naxes;
    int         libration_ncoeffs;

    //ancillary information
    double      jalpha;         //!< Start Julian Date for DE
    double      jomega;         //!< End Julian Date for DE
    double      jdelta;         //!< Julian Date delta for DE
    int         verbosity;

} Lgm_JPLephemInfo;


typedef struct Lgm_JPLephemBundle {

    double      tdb;
    double      **coeffs;
    int         coeffsAlloced;
    int         coefficient_count;
    double      days_per_set;
    double      *T;
    int         TAlloced;
    double      twot1;

} Lgm_JPLephemBundle;


/*
 *  Function prototypes
 */
Lgm_JPLephemInfo *Lgm_InitJPLephemInfo( int DEnum, int getBodies, int Verbosity );
void        Lgm_InitJPLephDefaults (int DEnum, int getBodies, int verbosity, Lgm_JPLephemInfo *jpl );
void        Lgm_FreeJPLephemInfo( Lgm_JPLephemInfo  *jpl );
void        Lgm_ReadJPLephem( Lgm_JPLephemInfo *jpl );

Lgm_JPLephemBundle *Lgm_InitJPLephemBundle( double tdb );
void        *Lgm_FreeJPLephemBundle( Lgm_JPLephemBundle *bundle );
void        Lgm_JPLephem_setup_object( int objName, Lgm_JPLephemInfo *jpl, Lgm_JPLephemBundle *bundle );
double      ***Lgm_JPL_getCoeffSet(int objName, Lgm_JPLephemInfo *jpl);
int         Lgm_JPL_getNCoeffs( int objName, Lgm_JPLephemInfo *jpl);
int         Lgm_JPL_getNAxes( int objName, Lgm_JPLephemInfo *jpl);
int         Lgm_JPL_getNSets( int objName, Lgm_JPLephemInfo *jpl);

void        Lgm_JPLephem_position( double tdb, int objName, Lgm_JPLephemInfo *jpl, Lgm_Vector *position);
