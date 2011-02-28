#ifndef LGM_WGS84_H
#define LGM_WGS84_H

// World Geodetic System 1984 (WGS84) parameters
#define WGS84_A         6378.1370               // semi-major axis of earth (equatorial radius) in km
#define WGS84_B         6356.7523142            // semi-minor axis of earth (polar radius) in km ( derived from b=a(1-f) )
#define WGS84_F         0.0033528106718309896   // f Flattening factor
#define WGS84_FINV      298.2572229328697       // 1/f where f is flattening factor
//#define WGS84_E2        0.006694380004260827    // first eccentricity squared ( e^2 = 1-b^2/a^2 = 2f-f^2)
#define WGS84_E2        0.00669437999014    // first eccentricity squared ( e^2 = 1-b^2/a^2 = 2f-f^2)
#define WGS84_E         0.08181919092890624     // first eccentricity ( e )
#define WGS84_EP2       0.006739496756586903    // second eccentricity squared ( ep^2 = a^2/b^2-1 = f*(2-f)/(1-f)^2 )
#define WGS84_EP        0.08209443803685366     // second eccentricity ep
#define WGS84_A2        40680631.59076899       // WGS84_A * WGS84_A (km^2)
#define WGS84_B2        40408299.98408706       // WGS84_B * WGS84_B (km^2)
#define WGS84_A2mB2     272331.6066819355       // WGS84_A2 - WGS84_B2 (km^2)
#define WGS84_E4        4.481472364144719e-05   // WGS84_E2*WGS84_E2
#define WGS84_1mE2      0.993305619995739       // 1 - WGS84_E2


#endif

/*
 *    $Id: Lgm_WGS84.h 46 2010-10-01 20:46:59Z mgh $
 */

