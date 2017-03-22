/*
 * ViewDriftShell.c:
 *
 *   Written by Mike Henderson
 *   (parts based on demo code written by written by Naofumi Yasufuku)
 *
 *   IMPORTANT NOTE: OpenGL uses column-major format for matrices.
 *                       Cg uses row-major format for matrices.
 *                        C also uses row-major format for matrices.
 *
 *
 *      In row-major the rows are stored contiguously.
 *      In col-major the cols are stored contiguously.
 *
 *
 */
#ifndef MAIN
#define MAIN
#endif

#ifndef LGM_IMAGE_DATA_DIR
#warning "hard-coding LGM_IMAGE_DATA_DIR because it was not in config.h"
#define LGM_IMAGE_DATA_DIR    "/usr/local/share/LanlGeoMag/Images"
#endif

#ifndef LGM_BLUEMARBLE_DATA_DIR
#warning "hard-coding LGM_BLUEMARBLE_DATA_DIR because it was not in config.h"
#define LGM_BLUEMARBLE_DATA_DIR    "/usr/local/share/LanlGeoMag/Images/BlueMarble/5400x2700"
#endif


#include "ViewDriftShell.h"
#include <Lgm_DynamicMemory.h>
#include "Vds_DriftShell.h"
#include <Lgm_Objects.h>

typedef struct TleInfo {

    char        Name[80];       // 0

    long int    SatNum;         // 1
    char        IntDesig[10];   // 2
    double      Epoch;          // 3
    double      d1MeanMotion;   // 4
    double      d2MeanMotion;   // 5
    double      BSTAR;          // 6
    int         ElemNum;        // 7

    double      Inc;            // 8
    double      RAAN;           // 9
    double      Ecc;            // 10
    double      AoP;            // 11
    double      MeanAnomaly;    // 12
    double      MeanMotion;     // 13
    int         RevNum;         // 14

    char        Line0[256];
    char        Line1[256];
    char        Line2[256];

    GtkWidget   *Line0Label;
    GtkWidget   *Line1Label;
    GtkWidget   *Line2Label;

} TleInfo;

TleInfo tle;


void BSTAR_to_STR( double bs, char *Str3 ) {

    int sgn;
    char *p, Str[80], Str2[80];
    long int a;
    double val, man;
    int ex;

    sgn = ( bs  < 0.0 ) ? -1 : 1;
    val = fabs( bs );

    //printf( "val = %.5e\n", val );
    sprintf( Str, "%.5e", val );
    //printf( "Str = %s\n", Str );
    
    // get mantissa
    p = strstr( Str, "." );
    strncpy( Str2, p-1, 7 );
    Str2[7] = '\0';
    //printf("Str2 = %s\n", Str2);
    sscanf( Str2, "%lf", &man );
    //printf( "man = %g\n", man );
    
    // get exponent
    p = strstr( Str, "e" );
    sscanf( p+1, "%d", &ex );
    //printf( "ex = %d\n", ex );
    
    // make man < 1
    man /= 10.0; ++ex;
    //printf( "man = %g\n", man );
    //printf( "ex = %d\n", ex );
    
    a = (long int)(man*1e5);
    //printf("STRING:  %05ld%d", a, ex );
    
    
    if ( a == 0 ) {
        sprintf( Str3, " 00000-0" );
    } else {
        if ( sgn < 0 ) {
            if ( ex < 0 ) {
                sprintf( Str3, "-%05ld-%d", a, abs(ex) );
            } else {
                sprintf( Str3, "-%05ld+%d", a, abs(ex) );
            }
        } else {
            if ( ex < 0 ) {
                sprintf( Str3, " %05ld-%d", a, abs(ex) );
            } else {
                sprintf( Str3, " %05ld+%d", a, abs(ex) );
            }
        }
    }

}


GtkWidget *PUKE_SATSEL_VBOX;


Vds_ObjectInfo  *ObjInfo;

void SolidCone( Lgm_Vector *u, double Fov, GLint slices, GLint stacks);
void CreateCutEllipsoid( double ra, double rb, int n);

/*
 *  Menu Bar stuff
 */
static GtkItemFactoryEntry MenuItems[] = {

    {"/File",                           NULL,           0,                              0,      "<Branch>"                              },
    {"/File/tearoff1",                  NULL,           0,                              0,      "<Tearoff>"                             },
//    {"/File/Print ...",               NULL,           PrintRasterFile,                0,      "<StockItem>",  GTK_STOCK_PRINT         },
    {"/File/Save View As ...",          "<CTL>S",       raise_SaveRasterFileDialog,     0,      "<StockItem>",  GTK_STOCK_SAVE_AS       },
    {"/File/Add MagEphem File ...",     NULL,           raise_OpenMagEphemFileDialog,   0,      "<StockItem>",  GTK_STOCK_OPEN          },
    {"/File/Add Earth-Fixed Sites ...", NULL,           raise_OpenMagEphemFileDialog,   0,      "<StockItem>",  GTK_STOCK_OPEN          },
    {"/File/sep1",                      NULL,           0,                              0,      "<Separator>"                           },
    {"/File/_Quit",                     "<CTL>Q",       gtk_main_quit,                  0,      "<StockItem>",  GTK_STOCK_QUIT          },

    {"/View",                           NULL,           0,                              0,      "<Branch>"                              },
    {"/View/tearoff2",                  NULL,           0,                              0,      "<Tearoff>"                             },
    {"/View/Full Screen",               "<CTL>L",       GoFullScreen,                   0,      "<StockItem>",  GTK_STOCK_FULLSCREEN    },
    {"/View/Undo Full Screen",          "<ALT>L",       GoNormalScreen,                 0,      "<StockItem>",  GTK_STOCK_LEAVE_FULLSCREEN    },

    {"/Tools",                          NULL,           0,                              0,      "<Branch>"                              },
    {"/Tools/tearoff3",                 NULL,           0,                              0,      "<Tearoff>"                             },
//    {"/Tools/Modify Color Table ...",   NULL,           raise_xBow,                     0,      "<StockItem>",  GTK_STOCK_SELECT_COLOR  },
//    {"/Tools/Fonts and Colors ...",     NULL,           raise_fontselectiondialog1,     0,      "<StockItem>",  GTK_STOCK_SELECT_FONT   },
//    {"/Tools/Line Plot Legend ...",     NULL,           raise_spacecraftlegenddialog1,  0,      "<StockItem>",  GTK_STOCK_PROPERTIES    },
//    {"/Tools/PSD Window ...",           NULL,           raise_PsdWindow,                0,      "<StockItem>",  GTK_STOCK_PROPERTIES    },
//    {"/Tools/View Drift Shell ...",     NULL,           raise_ViewDriftShell,           0,      "<StockItem>",  GTK_STOCK_PROPERTIES    },

    {"/Help",                           "<CTL>H",       0,                              0,      "<Branch>"                              },
    {"/Help/tearoff4",                  NULL,           0,                              0,      "<Tearoff>"                             },
    {"/Help/Help Contents ...",         NULL,           0,                              0,      "<StockItem>",  GTK_STOCK_HELP          },
    {"/Help/Release Notes ...",         NULL,           0,                              0,      "<StockItem>",  GTK_STOCK_INFO          },
//    {"/Help/About ...",                 "<CTL>A",       create_aboutdialog,             0,      "<StockItem>",  GTK_STOCK_ABOUT         }

};




GLfloat LightPosition[]        = {0.0, 3.0, 3.0, 0.0};

GtkFileFilter *PngFilter;
GtkObject   *StartYearSpinbutton_adj, *StartMonthSpinbutton_adj, *StartDaySpinbutton_adj;
GtkWidget   *StartYearSpinbutton, *StartMonthSpinbutton, *StartDaySpinbutton;
GtkWidget   *StartTimeLabel;
long int    StartDate;
int         StartYear, StartMonth, StartDay;
int         StartHour, StartMin, StartSec, StartMilliSec;
double      StartUT, StartJD;
GtkObject   *EndDaySpinbutton_adj;
GtkWidget   *EndDaySpinbutton;
GtkWidget   *EndTimeLabel;
long int    EndDate;
int         EndYear, EndMonth, EndDay;
int         EndHour, EndMin, EndSec, EndMilliSec;
double      EndUT, EndJD;
GtkWidget   *CurrentTimeLabel;
long int    CurrentDate;
int         CurrentYear, CurrentMonth, CurrentDay;
int         CurrentHour, CurrentMin, CurrentSec, CurrentMilliSec;
double      CurrentUT, CurrentJD;
double      TimeInc;
int         RunTime;
int         StepTime;
int         DumpFrames;
GtkWidget   *DumpFramesCheckbutton;
GtkWidget   *cFramesLabel;
int         MonthDays[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
char        *MonthStr[] = {"Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec" };
double      MaxStarMagnitude;
char        MapImageFilename[1024];


double qradius=0.0;
int    qcount=0;
double qalpha=1.0;

long int    QuadsLoaded[9];         // Id for quads that are load
GLuint      QuadsLoadedDL[9];       // Their corresponding Display lists
GLuint      QuadsLoadedTexId[9];    // Their corresponding Tex ids

/*
 * Conversion Flags
 */
Lgm_MagModelInfo   *mInfo;
Lgm_Vector          Sun;
double              TiltAngle;
double              DipoleTiltAngle;
double              RotAngle, RotAngle2, RotAngle3, RotAngle4, RotAngle5;
Lgm_Vector          RotAxis, RotAxis2, RotAxis3, RotAxis4, RotAxis5;
Lgm_Vector          DipoleOffset_sm;
double              ED_x0, ED_y0, ED_z0;
char                *Coords[] = { "", "ICRF2000", "MOD", "TOD", "TEME", "PEF", "WGS84", "GSE", "GSM", "SM", "ED"};
int                 ObserverCoords    = GSM_COORDS;
int                 StarsConvertFlag  = GEI2000_TO_GSM;
int                 SatsConvertFlag   = GEI2000_TO_GSM;
int                 AtmosConvertFlag  = GSM_TO_GSM;
int                 ObsToGeoConvertFlag  = GSM_TO_GEO;
int                 GeoToObsConvertFlag  = GEO_TO_GSM;
GtkWidget           *RadioCoordSystem[5];
GtkWidget           *RadioLightingStyle[3];

double IdentityMatrix[3][3] = { {1.0, 0.0, 0.0 }, {0.0, 1.0, 0.0 }, {0.0, 0.0, 1.0} };



void DrawScene( );
void CreateSphere( double r, int n);
void CreateEllipsoid( double ra, double rb, int n);
void IridiumFlare( double JD, _SgpTLE *TLE, Lgm_Vector *EarthToSun, int *Flag1, int *Flag2, int *Flag3, Lgm_Vector *g1, Lgm_Vector *g2, Lgm_Vector *g3, Lgm_Vector *P ) ;


/*
 * Texture names
 */
GLuint  Texture_Earth;
GLuint  Texture_TopSide;
GLuint  Texture_Moon;
GLuint  Texture_Logo;
GLuint  Texture_EqPlane;
GLuint  Texture_MeridPlane1;
GLuint  Texture_MeridPlane2;
GLuint  Texture_Debris;
GLuint  Texture_RocketBody;
GLuint  Texture_Spacecraft;
GLuint  Texture_Sun;
GLuint  Texture_TEMP;;
GLuint  Texture_HiResEarthQuad;

int     Logo_Width, Logo_Height;

int FrameNumber = 0;


static float        LightModelTwoSide[] = {GL_TRUE};
static float        LightModelOneSide[] = {GL_FALSE};

//static float        myModelMatrix[16];
static float        myCgViewMatrix[16];
static float        myGlViewMatrix[16];
//static float        myModelViewMatrix[16];
static float        myGlProjectionMatrix[16];
static float        myCgProjectionMatrix[16];
static float        myCgModelViewProjectionMatrix[16];

/*
 * Define names for Shaders
 */
GLuint g_shaderMyTest;
GLuint g_shaderFrontInit;
GLuint g_shaderFrontPeel;
GLuint g_shaderFrontBlend;
GLuint g_shaderFrontFinal;
GLuint FrontVertShaderInit;
GLuint FrontFragShaderInit;
GLuint FrontVertShaderPeel;
GLuint FrontFragShaderPeel;
GLuint FrontVertShaderBlend;
GLuint FrontFragShaderBlend;
GLuint FrontVertShaderFinal;
GLuint FrontFragShaderFinal;
GLuint ShadeVertex;
GLuint ShadeFragment;

GLuint g_quadDisplayList;

/*
 * Uniform var locations
 */
GLint g_shaderFrontInit_AlphaLoc;
GLint g_shaderFrontPeel_AlphaLoc;
GLint g_shaderFrontFinal_BackgroundColorLoc;
GLint g_shaderFrontPeel_DepthTexLoc;
GLint g_shaderFrontBlend_TempTexLoc;
GLint g_shaderFrontFinal_ColorTexLoc;












void raise_OpenMagEphemFileDialog( gpointer callback_data, guint callback_action, GtkWidget *menu_item ) {
    GtkWidget *w;
    w = CreateOpenMagEphemFileDialog( );
    gtk_widget_show( w );
    return;
}

void raise_SaveRasterFileDialog( gpointer callback_data, guint callback_action, GtkWidget *menu_item ) {
    GtkWidget *w;
    w = CreateSaveRasterFileDialog( );
    gtk_widget_show( w );
    return;
}

void GoFullScreen( gpointer callback_data, guint callback_action, GtkWidget *menu_item ) {
    gtk_window_fullscreen( GTK_WINDOW( ViewDriftShellWindow) );
    return;
}

void GoNormalScreen( gpointer callback_data, guint callback_action, GtkWidget *menu_item ) {
    gtk_window_unfullscreen( GTK_WINDOW( ViewDriftShellWindow) );
    return;
}



















void UpdateTimeDepQuants( long int CurrentDate, double CurrentUT ) {

    Lgm_Vector  u, v;
    double      Q[4], Q2[4], Q3[4], A[3][3];

    // Set transforms, get Tilt
    Lgm_Set_Coord_Transforms( CurrentDate, CurrentUT, mInfo->c );
    DipoleTiltAngle = mInfo->c->psi*DegPerRad;;
printf("MIKE MIKE:   CurrentDate, CurrentUT = %ld %g\n", CurrentDate, CurrentUT);
printf("DipoleTiltAngle = %g\n", DipoleTiltAngle);



    /*
     *  Convert the OBS_COORD->GEO rot matrix to a quaternion rotation
     *  And get matrix Asm_to_obs and quat Q2 needed to rot dipole axis properly
     *
     *  Quaternions:
     *               Q: for rotating the earth properly
     *              Q2: for rotating the dipole axis properly
     *              Q3: for rotating the sun vector properly
     *
     */

// we need to make sure the various matrices are defined no that we have
// overhauled the Ctrans lib...
    if        ( ObserverCoords == GEO_COORDS ) {

        Lgm_MatrixToQuat( IdentityMatrix, Q );


        //Lgm_MatTimesMat( mInfo->c->Agsm_to_sm, mInfo->c->Awgs84_to_gsm, A ); // A = Awgs84_to_sm
        Lgm_MatTimesMat( mInfo->c->Agsm_to_wgs84, mInfo->c->Asm_to_gsm, A ); // A = Asm_to_wgs84
        Lgm_MatrixToQuat( A, Q2 );

        //Lgm_MatrixToQuat( mInfo->c->Awgs84_to_gsm, Q3 );
        Lgm_MatrixToQuat( mInfo->c->Agsm_to_wgs84, Q3 );

    } else if ( ObserverCoords == GSM_COORDS ) {

//Lgm_Set_Coord_Transforms( CurrentDate, CurrentUT, mInfo->c );
//FIX ME
        //Lgm_MatrixToQuat( mInfo->c->Agsm_to_wgs84, Q );
        Lgm_MatrixToQuat( mInfo->c->Awgs84_to_gsm, Q );


        //Lgm_MatrixToQuat( mInfo->c->Agsm_to_sm, Q2 );
        Lgm_MatrixToQuat( mInfo->c->Asm_to_gsm, Q2 );

        Lgm_MatrixToQuat( IdentityMatrix, Q3 );

    } else if ( ObserverCoords == GEI2000_COORDS ) {

        //Lgm_MatrixToQuat( mInfo->c->Agei_to_wgs84, Q );
        Lgm_MatrixToQuat( mInfo->c->Awgs84_to_gei, Q );

        //Lgm_MatTimesMat( mInfo->c->Agsm_to_sm, mInfo->c->Amod_to_gsm, A ); // A = Amod_to_sm
        Lgm_MatTimesMat( mInfo->c->Agsm_to_mod, mInfo->c->Asm_to_gsm, A ); // A = Asm_to_mod
        Lgm_MatrixToQuat( A, Q2 );

        //Lgm_MatrixToQuat( mInfo->c->Amod_to_gsm, Q3 );
        Lgm_MatrixToQuat( mInfo->c->Agsm_to_mod, Q3 );

    } else if ( ObserverCoords == GSE_COORDS ) {

        //Lgm_MatTimesMat( mInfo->c->Agsm_to_wgs84, mInfo->c->Agse_to_gsm, A ); // A = Agse_to_wgs84
        Lgm_MatTimesMat( mInfo->c->Agsm_to_gse, mInfo->c->Awgs84_to_gsm, A ); // A = Awgs84_to_gse
        Lgm_MatrixToQuat( A, Q );

        //Lgm_MatTimesMat( mInfo->c->Agsm_to_sm, mInfo->c->Agse_to_gsm, A ); // A = Agse_to_sm
        Lgm_MatTimesMat( mInfo->c->Agsm_to_gse, mInfo->c->Asm_to_gsm, A ); // A = Asm_to_gse
        Lgm_MatrixToQuat( A, Q2 );

        //Lgm_MatrixToQuat( mInfo->c->Agse_to_gsm, Q3 );
        Lgm_MatrixToQuat( mInfo->c->Agsm_to_gse, Q3 );

    } else if ( ObserverCoords == SM_COORDS ) {

        //Lgm_MatTimesMat( mInfo->c->Agsm_to_wgs84, mInfo->c->Asm_to_gsm, A ); // A = Asm_to_wgs84
        Lgm_MatTimesMat( mInfo->c->Agsm_to_sm, mInfo->c->Awgs84_to_gsm, A ); // A = Awgs84_to_sm
        Lgm_MatrixToQuat( A, Q );

        Lgm_MatrixToQuat( IdentityMatrix, Q2 );

        //Lgm_MatrixToQuat( mInfo->c->Asm_to_gsm, Q3 );
        Lgm_MatrixToQuat( mInfo->c->Agsm_to_sm, Q3 );

    }


    Lgm_QuatToAxisAngle(  Q, &RotAngle,  &RotAxis );
    Lgm_QuatToAxisAngle( Q2, &RotAngle2, &RotAxis2 );
    Lgm_QuatToAxisAngle( Q3, &RotAngle3, &RotAxis3 );
printf("RotAxis = %g %g %g   RotAngle = %g\n", RotAxis.x, RotAxis.y, RotAxis.z, RotAngle);
printf("RotAxis2 = %g %g %g   RotAngle2 = %g\n", RotAxis2.x, RotAxis2.y, RotAxis2.z, RotAngle2);
printf("RotAxis3 = %g %g %g   RotAngle3 = %g\n", RotAxis3.x, RotAxis3.y, RotAxis3.z, RotAngle3);


    /*
     *  Convert the dipole offset components to SM coords.
     *  We may not wnat this in there....
     */
    u.x = mInfo->c->ED_x0; u.y = mInfo->c->ED_y0; u.z = mInfo->c->ED_z0;
    Lgm_Convert_Coords( &u, &DipoleOffset_sm, GEO_TO_SM, mInfo->c );



    // Set Position of sun
    u.x = 1.0, u.y = u.z = 0.0;
    Lgm_Convert_Coords( &u, &v, AtmosConvertFlag, mInfo->c );
    Sun        = v;
    aInfo->Sun = v; // This is the struct for the Atmosphere stuff
printf(" Sun = %g %g %g\n", v.x, v.y, v.z);


}

// transpose
static void TransposeMatrix(float src[16], float dest[16] ) {
    int i, j;
    for (i=0; i<4; i++) {
        for (j=0; j<4; j++) {
            *(dest  +i*4 +j ) = *(src + j*4+i);
        }
    }
}

// Simple 4x4 matrix by 4x4 matrix multiply.
static void multMatrix(float dst[16], const float src1[16], const float src2[16]) {

  float tmp[16];
  int i, j;

  for (i=0; i<4; i++) {
    for (j=0; j<4; j++) {
      tmp[i*4+j] = src1[i*4+0] * src2[0*4+j] +
                   src1[i*4+1] * src2[1*4+j] +
                   src1[i*4+2] * src2[2*4+j] +
                   src1[i*4+3] * src2[3*4+j];
    }
  }
  /* Copy result to dst (so dst can also be src1 or src2). */
  for (i=0; i<16; i++) dst[i] = tmp[i];
}


/*
 * Set up matrix to transform from:
 *
 *     Eye Space -> Clip Space
 *
 *  This implements a Perspective projection.
 *
 *
 */
void SetProjectionMatrix( double FieldOfView, double Aspect, double zNear, double zFar, float m[16] ){

    double SinFieldOfView, dZ, dZinv, f;
    double Ang = 0.5*FieldOfView*M_PI/180.0;

    dZ = zNear - zFar;
    dZinv = 1.0/dZ;
    SinFieldOfView = sin( Ang );
    f = cos( Ang ) / SinFieldOfView; // f = cotangent( FieldOfView/2.0 )

    /*
     * Set up in row-major form (i.e. successive rows are contiguous in memory)
     */
    m[0*4+0] = f/Aspect;    m[0*4+1] = 0.0;   m[0*4+2] =  0.0;                 m[0*4+3] = 0.0;
    m[1*4+0] = 0.0;         m[1*4+1] = f;     m[1*4+2] =  0.0;                 m[1*4+3] = 0.0;
    m[2*4+0] = 0.0;         m[2*4+1] = 0.0;   m[2*4+2] = (zFar+zNear)*dZinv;   m[2*4+3] = 2.0*zNear*zFar*dZinv;
    m[3*4+0] = 0.0;         m[3*4+1] = 0.0;   m[3*4+2] = -1.0;                 m[3*4+3] = 0.0;

}

void SetTranslationScaleMatrix( double x, double y, double z, double s, float m[16] ){

    m[0*4+0] = s;      m[0*4+1] = 0.0;   m[0*4+2] = 0.0;  m[0*4+3] = x;
    m[1*4+0] = 0.0;    m[1*4+1] = s;     m[1*4+2] = 0.0;  m[1*4+3] = y;
    m[2*4+0] = 0.0;    m[2*4+1] = 0.0;   m[2*4+2] = s;    m[2*4+3] = z;
    m[3*4+0] = 0.0;    m[3*4+1] = 0.0;   m[3*4+2] = 0.0;  m[3*4+3] = 1.0;

}

/*
 * Set up matrix to transform from:
 *
 *     World Space -> Eye Space
 *
 *
 */
void SetViewMatrix( double Camera_x, double Camera_y, double Camera_z,
                    double Lookat_x, double Lookat_y, double Lookat_z,
                    double Up_x, double Up_y, double Up_z, float m[16] ){

    double  X[3], Y[3], Z[3], mag, mag_inv, Tx, Ty, Tz;

    /*
     * Camera -  location of the Camera in World Space.
     * Lookat -  the location in World Space that the Camera should look at.
     * Up     -  a vector specifying which direction in (world coords) will
     *           be up in the final view.
     *
     * Construct X, Y, and Z unit vectors. These are the world space coord unit
     * vectors in Eye space
     *
     *  Also compute the translation needed to position the camera properly
     */

    /*
     * Z unit vector points from Camera to Lookat
     */
    Z[0] = Camera_x - Lookat_x;
    Z[1] = Camera_x - Lookat_x;
    Z[2] = Camera_x - Lookat_x;

    /*
     * make Z a unit vector
     */
    mag = sqrt(Z[0]*Z[0] + Z[1]*Z[1] + Z[2]*Z[2]);
    if (mag > 0.0) {
        mag_inv = 1.0/mag;
        Z[0] *= mag_inv; Z[1] *= mag_inv; Z[2] *= mag_inv;
    } else {
        printf("Camera and Lookat are the same point!\n");
        exit(1);
    }

    /*
     *  Yhat is in Up direction (Up need not be normalized)
     */
    Y[0] = Up_x; Y[1] = Up_y; Y[2] = Up_z;


    /*
     * X = Y cross Z
     */
    X[0] =  Y[1]*Z[2] - Y[2]*Z[1]; X[1] = -Y[0]*Z[2] + Y[2]*Z[0]; X[2] =  Y[0]*Z[1] - Y[1]*Z[0];

    /*
     * recompute Y to make an orthog. system. I.e. Y = Z cross X
     */
    Y[0] =  Z[1]*X[2] - Z[2]*X[1]; Y[1] = -Z[0]*X[2] + Z[2]*X[0]; Y[2] =  Z[0]*X[1] - Z[1]*X[0];

    /*
     * Normlaize X
     */
    mag = sqrt(X[0]*X[0] + X[1]*X[1] + X[2]*X[2]);
    if (mag > 0.0) {
        mag_inv = 1.0/mag;
        X[0] *= mag_inv; X[1] *= mag_inv; X[2] *= mag_inv;
    } else {
        printf("Magnitude of X is <= 0.0 !\n");
        exit(1);
    }

    /*
     * Normlaize Y
     */
    mag = sqrt(Y[0]*Y[0] + Y[1]*Y[1] + Y[2]*Y[2]);
    if (mag > 0.0) {
        mag_inv = 1.0/mag;
        Y[0] *= mag_inv; Y[1] *= mag_inv; Y[2] *= mag_inv;
    } else {
        printf("Magnitude of Y is <= 0.0 !\n");
        exit(1);
    }

    /*
     *  Compute the translation to position the camera properly
     *
     *      Tx = -( Xhat dot Camera )
     *      Ty = -( Yhat dot Camera )
     *      Tz = -( Zhat dot Camera )
     */
    Tx = -( X[0]*Camera_x + X[1]*Camera_y + X[2]*Camera_z );
    Ty = -( Y[0]*Camera_x + Y[1]*Camera_y + Y[2]*Camera_z );
    Tz = -( Z[0]*Camera_x + Z[1]*Camera_y + Z[2]*Camera_z );




    m[0*4+0] = X[0];    m[0*4+1] = X[1];  m[0*4+2] = X[2];  m[0*4+3] = Tx;
    m[1*4+0] = Y[0];    m[1*4+1] = Y[1];  m[1*4+2] = Y[2];  m[1*4+3] = Ty;
    m[2*4+0] = Z[0];    m[2*4+1] = Z[1];  m[2*4+2] = Z[2];  m[2*4+3] = Tz;
    m[3*4+0] = 0.0;     m[3*4+1] = 0.0;   m[3*4+2] = 0.0;   m[3*4+3] = 1.0;

}








// from dual depth peeling demo
#define CHECK_GL_ERRORS {}
int g_imageWidth = 1024;
int g_imageHeight = 768;
double g_imageRat;

GLuint g_dualBackBlenderFboId;
GLuint g_dualPeelingSingleFboId;
GLuint g_dualDepthTexId[2];
GLuint g_dualFrontBlenderTexId[2];
GLuint g_dualBackTempTexId[2];
GLuint g_dualBackBlenderTexId;

GLuint g_frontFboId[2];
GLuint g_frontDepthTexId[2];
GLuint g_frontColorTexId[2];
GLuint g_frontColorBlenderTexId;
GLuint g_frontColorBlenderFboId;

GLuint g_accumulationTexId[2];
GLuint g_accumulationFboId;

GLenum g_drawBuffers[] = {GL_COLOR_ATTACHMENT0_EXT,
                           GL_COLOR_ATTACHMENT1_EXT,
                           GL_COLOR_ATTACHMENT2_EXT,
                           GL_COLOR_ATTACHMENT3_EXT,
                           GL_COLOR_ATTACHMENT4_EXT,
                           GL_COLOR_ATTACHMENT5_EXT,
                           GL_COLOR_ATTACHMENT6_EXT };


GLfloat g_opacity = 0.6;
float g_white[3] = {1.0,1.0,1.0};
float g_black[3] = {0.0};
float *g_backgroundColor = g_black;


GLuint LoadShaderFromFile( char *Filename, GLenum ShaderType) {

    GLuint  shader;
    GLint   r;
    int     nLines, len, i;
    char    *Str, *p;
    char    **Lines;
    FILE    *fp;

    // open file
    if ( (fp = fopen( Filename, "r" )) == NULL ) {
        printf("Could not open shader file: %s\n", Filename );
        return( 0 );
    }

    // determine how many Lines there are
    nLines = 0;
    Str = (char *)malloc( 4097*sizeof(char) );
    while ( fgets( Str, 400, fp )  ) {
        /*
         * Get rid of // style comments -- they seem to cause a problem -- dont know why
         * Also get rid of blank lines.
         */
        p = strstr( Str, "//" ); if (p != NULL) *p = '\0'; // get rid of //-style comments
        p = index( Str, '\n' );  if (p != NULL) *p = '\0';
        p = index( Str, '\r' );  if (p != NULL) *p = '\0';
        len = strlen( Str );
        if (len > 0 ) {
            ++nLines;
        }
    }
    rewind( fp ); // rewind fp back to start


    /*
     * Alloc memory for the Lines array
     */
    Lines = (char **)calloc( nLines, sizeof(char *) );
    i = 0;
    //printf( "Shader File: %s\n", Filename );
    while ( fgets( Str, 400, fp )  ) {
        p = strstr( Str, "//" ); if (p != NULL) *p = '\0'; // get rid of //-style comments
        p = index( Str, '\n' );  if (p != NULL) *p = '\0';
        p = index( Str, '\r' );  if (p != NULL) *p = '\0';
        len = strlen( Str );
        if (len > 0 ) {
            Lines[i] = (char *)calloc( len+1, sizeof(char) );
            strcpy( Lines[i], Str );
            //printf("Lines[%2d] = \"%s\"\n", i, Lines[i]);
            ++i;
        }
    }
    //printf("\n");

    fclose(fp);

    shader = glCreateShader( ShaderType );
    glShaderSource( shader, nLines, (const  char **)Lines, (int *)NULL );

    glCompileShader( shader );
    glGetShaderiv( shader, GL_COMPILE_STATUS, &r );
    if (r) {
        printf("Succesfully compiled shader: %s\n", Filename);
    }  else {
        glGetShaderInfoLog( shader, 4096, &len, Str);
        printf("Error compiling shader: %s\n", Filename);
        printf("Log:\n%s\n\n", Str);
    }

    for (i=0; i<nLines; i++) free( Lines[i] );
    free( Lines );
    free( Str );


    return( shader );

}




    void BuildShaders2() {

        printf("\nLoading shaders...\n");
        //ShadeVertex   = LoadShaderFromFile( "Shaders/shade_vertex2.glsl",   GL_VERTEX_SHADER );
        //ShadeFragment = LoadShaderFromFile( "Shaders/shade_fragment2.glsl", GL_FRAGMENT_SHADER );
        ShadeVertex   = LoadShaderFromFile( "Shaders/CH14-WardBRDF.vert",   GL_VERTEX_SHADER );
        ShadeFragment = LoadShaderFromFile( "Shaders/CH14-WardBRDF.frag", GL_FRAGMENT_SHADER );

        g_shaderMyTest = glCreateProgram();
    glAttachShader( g_shaderMyTest, ShadeVertex );
    glAttachShader( g_shaderMyTest, ShadeFragment );
    glLinkProgram( g_shaderMyTest );

}


void BuildShaders() {

	printf("\nLoading shaders...\n");

    FrontVertShaderInit  = LoadShaderFromFile( "Shaders/front_peeling_init_vertex.glsl",    GL_VERTEX_SHADER );
    FrontFragShaderInit  = LoadShaderFromFile( "Shaders/front_peeling_init_fragment.glsl",  GL_FRAGMENT_SHADER );
    FrontVertShaderPeel  = LoadShaderFromFile( "Shaders/front_peeling_peel_vertex.glsl",    GL_VERTEX_SHADER );
    FrontFragShaderPeel  = LoadShaderFromFile( "Shaders/front_peeling_peel_fragment.glsl",  GL_FRAGMENT_SHADER );
    FrontVertShaderBlend = LoadShaderFromFile( "Shaders/front_peeling_blend_vertex.glsl",    GL_VERTEX_SHADER );
    FrontFragShaderBlend = LoadShaderFromFile( "Shaders/front_peeling_blend_fragment.glsl",  GL_FRAGMENT_SHADER );
    FrontVertShaderFinal = LoadShaderFromFile( "Shaders/front_peeling_final_vertex.glsl",   GL_VERTEX_SHADER );
    FrontFragShaderFinal = LoadShaderFromFile( "Shaders/front_peeling_final_fragment.glsl", GL_FRAGMENT_SHADER );


    ShadeVertex   = LoadShaderFromFile( "Shaders/shade_vertex.glsl",   GL_VERTEX_SHADER );
    ShadeFragment = LoadShaderFromFile( "Shaders/shade_fragment.glsl", GL_FRAGMENT_SHADER );



    // Create g_shaderFrontInit program
    g_shaderFrontInit = glCreateProgram();
    glAttachShader( g_shaderFrontInit, ShadeVertex );
    glAttachShader( g_shaderFrontInit, FrontVertShaderInit );
    glAttachShader( g_shaderFrontInit, ShadeFragment );
    glAttachShader( g_shaderFrontInit, FrontFragShaderInit );
    glLinkProgram( g_shaderFrontInit );
    g_shaderFrontInit_AlphaLoc = glGetUniformLocation( g_shaderFrontInit, "Alpha" );


    // Create g_shaderFrontPeel program
    g_shaderFrontPeel = glCreateProgram();
    glAttachShader( g_shaderFrontPeel, ShadeVertex );
    glAttachShader( g_shaderFrontPeel, FrontVertShaderPeel );
    glAttachShader( g_shaderFrontPeel, ShadeFragment );
    glAttachShader( g_shaderFrontPeel, FrontFragShaderPeel );
    glLinkProgram( g_shaderFrontPeel );
    g_shaderFrontPeel_AlphaLoc = glGetUniformLocation( g_shaderFrontPeel, "Alpha" );
    g_shaderFrontPeel_DepthTexLoc = glGetUniformLocation( g_shaderFrontPeel, "DepthTex" );

    // Create g_shaderFrontBlend program
    g_shaderFrontBlend = glCreateProgram();
    glAttachShader( g_shaderFrontBlend, FrontVertShaderBlend );
    glAttachShader( g_shaderFrontBlend, FrontFragShaderBlend );
    glLinkProgram( g_shaderFrontBlend );
    g_shaderFrontBlend_TempTexLoc = glGetUniformLocation( g_shaderFrontBlend, "TempTex" );

    // Create g_shaderFrontFinal program
    g_shaderFrontFinal = glCreateProgram();
    glAttachShader( g_shaderFrontFinal, ShadeVertex );
    glAttachShader( g_shaderFrontFinal, FrontVertShaderFinal );
    glAttachShader( g_shaderFrontFinal, ShadeFragment );
    glAttachShader( g_shaderFrontFinal, FrontFragShaderFinal );
    glLinkProgram( g_shaderFrontFinal );
    g_shaderFrontFinal_BackgroundColorLoc = glGetUniformLocation( g_shaderFrontFinal, "BackgroundColor" );
    g_shaderFrontFinal_ColorTexLoc = glGetUniformLocation( g_shaderFrontFinal, "ColorTex" );


    printf("\n");

}


void InitFrontPeelingRenderTargets()
{
    int i;
    glGenTextures(2, g_frontDepthTexId);
    glGenTextures(2, g_frontColorTexId);
    glGenFramebuffersEXT(2, g_frontFboId);

    for (i = 0; i < 2; i++) {

        glBindTexture(GL_TEXTURE_RECTANGLE_ARB, g_frontDepthTexId[i]);
        glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_S, GL_CLAMP);
        glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_T, GL_CLAMP);
        glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexImage2D(GL_TEXTURE_RECTANGLE_ARB, 0, GL_DEPTH_COMPONENT32F_NV,
                     g_imageWidth, g_imageHeight, 0, GL_DEPTH_COMPONENT, GL_FLOAT, NULL);

        glBindTexture(GL_TEXTURE_RECTANGLE_ARB, g_frontColorTexId[i]);
        glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_S, GL_CLAMP);
        glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_T, GL_CLAMP);
        glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexImage2D(GL_TEXTURE_RECTANGLE_ARB, 0, GL_RGBA, g_imageWidth, g_imageHeight,
                     0, GL_RGBA, GL_FLOAT, 0);

        glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, g_frontFboId[i]);
        glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT,
                                  GL_TEXTURE_RECTANGLE_ARB, g_frontDepthTexId[i], 0);
        glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT,
                                  GL_TEXTURE_RECTANGLE_ARB, g_frontColorTexId[i], 0);
    }

    glGenTextures(1, &g_frontColorBlenderTexId);
    glBindTexture(GL_TEXTURE_RECTANGLE_ARB, g_frontColorBlenderTexId);
    glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_T, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexImage2D(GL_TEXTURE_RECTANGLE_ARB, 0, GL_RGBA, g_imageWidth, g_imageHeight,
                 0, GL_RGBA, GL_FLOAT, 0);

    glGenFramebuffersEXT(1, &g_frontColorBlenderFboId);
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, g_frontColorBlenderFboId);
    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT,
                              GL_TEXTURE_RECTANGLE_ARB, g_frontDepthTexId[0], 0);
    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT,
                              GL_TEXTURE_RECTANGLE_ARB, g_frontColorBlenderTexId, 0);
    CHECK_GL_ERRORS;
}



void InitDualPeelingRenderTargets() {

    int i;

	glGenTextures(2, g_dualDepthTexId);
	glGenTextures(2, g_dualFrontBlenderTexId);
	glGenTextures(2, g_dualBackTempTexId);
	glGenFramebuffersEXT(1, &g_dualPeelingSingleFboId);
	for ( i = 0; i < 2; i++ ) {

		glBindTexture(GL_TEXTURE_RECTANGLE_ARB, g_dualDepthTexId[i]);
		glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_S, GL_CLAMP);
		glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_T, GL_CLAMP);
		glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexImage2D(GL_TEXTURE_RECTANGLE_ARB, 0, GL_FLOAT_RG32_NV, g_imageWidth, g_imageHeight,
					 0, GL_RGB, GL_FLOAT, 0);

		glBindTexture(GL_TEXTURE_RECTANGLE_ARB, g_dualFrontBlenderTexId[i]);
		glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_S, GL_CLAMP);
		glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_T, GL_CLAMP);
		glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexImage2D(GL_TEXTURE_RECTANGLE_ARB, 0, GL_RGBA, g_imageWidth, g_imageHeight,
					 0, GL_RGBA, GL_FLOAT, 0);

		glBindTexture(GL_TEXTURE_RECTANGLE_ARB, g_dualBackTempTexId[i]);
		glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_S, GL_CLAMP);
		glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_T, GL_CLAMP);
		glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexImage2D(GL_TEXTURE_RECTANGLE_ARB, 0, GL_RGBA, g_imageWidth, g_imageHeight,
					 0, GL_RGBA, GL_FLOAT, 0);
	}

	glGenTextures(1, &g_dualBackBlenderTexId);
	glBindTexture(GL_TEXTURE_RECTANGLE_ARB, g_dualBackBlenderTexId);
	glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_T, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexImage2D(GL_TEXTURE_RECTANGLE_ARB, 0, GL_RGB, g_imageWidth, g_imageHeight,
				 0, GL_RGB, GL_FLOAT, 0);

	glGenFramebuffersEXT(1, &g_dualBackBlenderFboId);
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, g_dualBackBlenderFboId);
	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT,
							  GL_TEXTURE_RECTANGLE_ARB, g_dualBackBlenderTexId, 0);

	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, g_dualPeelingSingleFboId);

	int j = 0;
		glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT,
								  GL_TEXTURE_RECTANGLE_ARB, g_dualDepthTexId[j], 0);
		glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT1_EXT,
								  GL_TEXTURE_RECTANGLE_ARB, g_dualFrontBlenderTexId[j], 0);
		glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT2_EXT,
								  GL_TEXTURE_RECTANGLE_ARB, g_dualBackTempTexId[j], 0);

	j = 1;
		glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT3_EXT,
								  GL_TEXTURE_RECTANGLE_ARB, g_dualDepthTexId[j], 0);
		glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT4_EXT,
								  GL_TEXTURE_RECTANGLE_ARB, g_dualFrontBlenderTexId[j], 0);
		glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT5_EXT,
								  GL_TEXTURE_RECTANGLE_ARB, g_dualBackTempTexId[j], 0);

	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT6_EXT,
							  GL_TEXTURE_RECTANGLE_ARB, g_dualBackBlenderTexId, 0);

	CHECK_GL_ERRORS;
}

void RenderFrontToBackPeeling() {

	// ---------------------------------------------------------------------
	// 1. Initialize Min Depth Buffer
	// ---------------------------------------------------------------------

	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, g_frontColorBlenderFboId);
	glDrawBuffer(g_drawBuffers[0]);

	glClearColor(0, 0, 0, 1);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glEnable(GL_DEPTH_TEST);

    glUseProgram( g_shaderFrontInit );
    glUniform1f( g_shaderFrontInit_AlphaLoc, g_opacity );
    DrawScene();
    glUseProgram( 0 );

	CHECK_GL_ERRORS;

	// ---------------------------------------------------------------------
	// 2. Depth Peeling + Blending
	// ---------------------------------------------------------------------

int g_numPasses = 4;
	int numLayers = (g_numPasses - 1) * 2;
    int layer, currId, prevId;
//	for (int layer = 1; g_useOQ || layer < numLayers; layer++) {
	for (layer = 1; layer < numLayers; layer++) {
		currId = layer % 2;
		prevId = 1 - currId;

		glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, g_frontFboId[currId]);
		glDrawBuffer(g_drawBuffers[0]);

		glClearColor(0, 0, 0, 0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glDisable(GL_BLEND);
        glEnable(GL_DEPTH_TEST);

//		if (g_useOQ) {
//			glBeginQuery(GL_SAMPLES_PASSED_ARB, g_queryId);
//		}

        glUseProgram( g_shaderFrontPeel );
        glActiveTexture( GL_TEXTURE0+0 );
        glBindTexture( GL_TEXTURE_RECTANGLE_ARB, g_frontDepthTexId[prevId] );
        glUniform1i( g_shaderFrontPeel_DepthTexLoc, 0 );
        glActiveTexture( GL_TEXTURE0 );
        glUniform1f( g_shaderFrontPeel_AlphaLoc, g_opacity );
        DrawScene();
        glUseProgram( 0 );

//		if (g_useOQ) {
//			glEndQuery(GL_SAMPLES_PASSED_ARB);
//		}

		CHECK_GL_ERRORS;

		glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, g_frontColorBlenderFboId);
		glDrawBuffer(g_drawBuffers[0]);

        glDisable(GL_DEPTH_TEST);
        glEnable(GL_BLEND);

		glBlendEquation(GL_FUNC_ADD);
		glBlendFuncSeparate(GL_DST_ALPHA, GL_ONE, GL_ZERO, GL_ONE_MINUS_SRC_ALPHA);
		
        glUseProgram( g_shaderFrontBlend );
        glActiveTexture( GL_TEXTURE0+0 );
        glBindTexture( GL_TEXTURE_RECTANGLE_ARB, g_frontColorTexId[currId] );
        glUniform1i( g_shaderFrontBlend_TempTexLoc, 0 );
        glActiveTexture( GL_TEXTURE0 );
		glCallList(g_quadDisplayList);
        glUseProgram( 0 );

		glDisable(GL_BLEND);

		CHECK_GL_ERRORS;

//		if (g_useOQ) {
//			GLuint sample_count;
//			glGetQueryObjectuiv(g_queryId, GL_QUERY_RESULT_ARB, &sample_count);
//			if (sample_count == 0) {
//				break;
//			}
//		}
	}

	// ---------------------------------------------------------------------
	// 3. Final Pass
	// ---------------------------------------------------------------------

	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
	glDrawBuffer(GL_BACK);
	glDisable(GL_DEPTH_TEST);

    glUseProgram( g_shaderFrontFinal );
    glUniform3fv( g_shaderFrontFinal_BackgroundColorLoc, 3, (GLfloat *)g_backgroundColor );
    glActiveTexture( GL_TEXTURE0+0 );
    glBindTexture( GL_TEXTURE_RECTANGLE_ARB, g_frontColorBlenderTexId );
    glUniform1i( g_shaderFrontFinal_ColorTexLoc, 0 );
    glActiveTexture( GL_TEXTURE0 );
	glCallList(g_quadDisplayList);
    glUseProgram( 0 );

	CHECK_GL_ERRORS;
}











GLfloat mat_red_diffuse[] = { 0.7, 0.0, 0.1, 1.0 };
GLfloat mat_green_diffuse[] = { 0.0, 0.7, 0.1, 1.0 };
GLfloat mat_blue_diffuse[] = { 0.0, 0.1, 0.7, 1.0 };
GLfloat mat_yellow_diffuse[] = { 0.7, 0.8, 0.1, 1.0 };
GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
GLfloat mat_shininess[] = { 100.0 };
GLfloat knots[8] = { 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0 };
GLfloat pts1[4][4][3], pts2[4][4][3];
GLfloat pts3[4][4][3], pts4[4][4][3];
GLUnurbsObj *nurb;
int u, v;



#include "trackball.h"

#define DIG_2_RAD (G_PI / 180.0)
#define RAD_2_DIG (180.0 / G_PI)

#define ANIMATE_THRESHOLD 1.5

#define VIEW_SCALE_MAX 10.0
#define VIEW_SCALE_MIN 0.1

#define NUM_LISTS 7
int res=20;
int lowres=40;
int hires=140; // nurb resolution
float   Orange[3]   =   {1.0, 0.4, 0.0};
float colors[][3] =    {
                        {1.000000, 0.109804, 0.109804}, // 0
                        {0.835294, 0.341176, 0.141176},
                        {0.713725, 0.470588, 0.156863},
                        {0.705882, 0.647059, 0.160784},
                        {0.772549, 0.909804, 0.156863},
                        {0.603922, 0.792157, 0.109804}, // 5
                        {0.501961, 0.941176, 0.207843},
                        {0.207843, 0.639216, 0.192157},
                        {0.105882, 0.423529, 0.133333},
                        {0.121569, 0.733333, 0.752941},
                        {0.250980, 0.529412, 0.650980}, // 10
                        {0.172549, 0.254902, 0.776471},
                        {0.101961, 0.090196, 0.552941},
                        {0.380392, 0.282353, 0.631373},
                        {0.623529, 0.439216, 0.788235},
                        {0.729411, 0.345098, 0.878431}, // 15
                        {0.517647, 0.133333, 0.447058},
                        {0.325490, 0.203921, 0.098039},
                        {0.784314, 0.435294, 0.921569} };

float   colors2[][3] = {{0.000000, 0.000000, 0.000000},
                        {0.490196, 0.184314, 0.835294},
                        {0.231373, 0.184314, 0.835294},
                        {0.184314, 0.411765, 0.835294},
                        {0.184314, 0.737255, 0.835294},
                        {0.184314, 0.835294, 0.643137},
                        {0.184314, 0.835294, 0.305882},
                        {0.568627, 0.835294, 0.184314},
                        {0.796078, 0.835294, 0.184314},
                        {0.835294, 0.721569, 0.184314},
                        {0.490196, 0.184314, 0.184314},
                        {0.835294, 0.552941, 0.184314},
                        {0.835294, 0.427451, 0.184314},
                        {0.835294, 0.184314, 0.184314},
                        {0.835294, 0.184314, 0.721569},
                        {0.705882, 0.184314, 0.835294},
                        {1.000000, 1.000000, 1.000000},
                        {0.000000, 0.000000, 0.000000} };

//float colors[][3] =   { {0.000000, 0.000000, 0.000000},
//                        {0.497212, 0.186952, 0.847248},
//                        {0.261105, 0.207999, 0.942635},
//                        {0.163216, 0.652864, 0.739681},
//                        {0.172225, 0.780508, 0.600954},
//                        {0.5536, 0.81322, 0.179443},
//                        {0.681273, 0.714834, 0.157733},
//                        {0.882935, 0.331983, 0.331983},
//                        {0.820086, 0.542874, 0.180958},
//                        {0.9546, 0.21064, 0.21064},
//                        {0.746409, 0.164701, 0.644785},
//                        {0.746409, 0.164701, 0.644785},
//                        {0.636484, 0.166193, 0.753172}};











int ReadPng( char *Filename, int *Width, int *Height, GLubyte **pImage );
void MakeTube(double *X, double *Y, double *Z, int NumCurvePoints, int NumCirclePoints, double TubeRadius );




int         nInterpPoints, nFieldPoints[30], nShellPoints;


int         nPnts[30][LGM_LSTARINFO_MAX_FL], nPnts2[30][LGM_LSTARINFO_MAX_FL], gap[1000];
double      s_gsm[30][LGM_LSTARINFO_MAX_FL][1000], x_gsm[30][LGM_LSTARINFO_MAX_FL][1000], y_gsm[30][LGM_LSTARINFO_MAX_FL][1000], z_gsm[30][LGM_LSTARINFO_MAX_FL][1000];
double      s2_gsm[30][LGM_LSTARINFO_MAX_FL][1000], x2_gsm[30][LGM_LSTARINFO_MAX_FL][1000], y2_gsm[30][LGM_LSTARINFO_MAX_FL][1000], z2_gsm[30][LGM_LSTARINFO_MAX_FL][1000];
double      x3_gsm[30][LGM_LSTARINFO_MAX_FL][200], y3_gsm[30][LGM_LSTARINFO_MAX_FL][200], z3_gsm[30][LGM_LSTARINFO_MAX_FL][200];
double      nx3_gsm[30][LGM_LSTARINFO_MAX_FL][200], ny3_gsm[30][LGM_LSTARINFO_MAX_FL][200], nz3_gsm[30][LGM_LSTARINFO_MAX_FL][200];
double      x4_gsm[30][10*LGM_LSTARINFO_MAX_FL][200], y4_gsm[30][10*LGM_LSTARINFO_MAX_FL][200], z4_gsm[30][10*LGM_LSTARINFO_MAX_FL][200];
double      nx4_gsm[30][10*LGM_LSTARINFO_MAX_FL][200], ny4_gsm[30][10*LGM_LSTARINFO_MAX_FL][200], nz4_gsm[30][10*LGM_LSTARINFO_MAX_FL][200];

gulong      PitchAngleAllHandler;
gulong      PitchAngleHandler[90];
GtkWidget   *PitchAngleCheckMenuItem[91];
//int         ShowSun                 =  0;
int         ShowStars               =  1;
int         ShowEarth               =  1;
int         ShowThemisFovs          =  1;
int         ShowAtmosphere          =  1;
int         ShowFullFieldLine       =  0;
int         ShowDriftShellSurface   =  0;
int         ShowAllPitchAngles      =  0;
int         ShowAllPitchAngles2     =  0;
int         ShowIridiumSatToGround  =  0;
int         ShowIridiumSatToSun     =  0;
int         ShowPitchAngle[90];
int         ShowPitchAngle2[90];
int         LightingStyle           =  0;

int         nSatTypes = 12;
int         ShowSatellites[12];
char        *SatTypes[] = { "amateur", "catalog", "geo", "globalstar", "inmarsat", "intelsat", "iridium", "navigation", "orbcomm", "spec_interest", "visible", "weather"};


static GLuint ZPSAxesDL           = 0;  // Display list for xyz axes
static GLuint AxesDL           = 0;  // Display list for xyz axes
//static GLuint SMAxesDL         = 0;  // Display list for GSM xyz axes
static GLuint EarthDL          = 0;  // Display List for Earth
static GLuint TopSideDL        = 0;  // Display List for TopSide Image
static GLuint GeoMarkersDL     = 0;  // Display List for Earth
static GLuint MoonDL           = 0;  // Display List for Moon
static GLuint LogoDL           = 0;  // Display List for Logo Image
static GLuint SunDL            = 0;  // Display List for Sun
static GLuint DipoleAxisDL     = 0;  // Display List for Dipole Axis
static GLuint SunDirectionDL   = 0;  // Display List for Direction to Sun
static GLuint ScPositionDL     = 0;  // Display List for S/C position
static GLuint StarsDL          = 0;  // Display List for Stars
static GLuint SatsDL           = 0;  // Display List for Sats
static GLuint ThemisFovsDL     = 0;  // Display List for Themis ASI FOVs
static GLuint SatOrbitsDL      = 0;  // Display List for Sat Orbits
static GLuint EqPlaneDL        = 0;  // Display List
static GLuint EqPlaneGridDL    = 0;  // Display List
static GLuint MeridPlane1DL     = 0;  // Display List
static GLuint MeridPlane2DL     = 0;  // Display List
static GLuint HiResEarthQuadDL = 0;  // Display List







typedef struct _MaterialProp {
  GLfloat ambient[4];
  GLfloat diffuse[4];
  GLfloat specular[4];
  GLfloat shininess;
} MaterialProp;


typedef struct _GuiInfo {
    MaterialProp *FieldLineMaterial;
    MaterialProp *DriftShellMaterial;

    GtkWidget   **FieldLineShowPitchAngleButton;
    gulong      *FieldLineShowPitchAngleButtonHandler;

    GtkWidget   **DriftShellShowPitchAngleButton;
    gulong      *DriftShellShowPitchAngleButtonHandler;

    GtkWidget   **FieldLineDiffuseColorButton;
    GtkWidget   **FieldLineAmbientColorButton;
    GtkWidget   **FieldLineSpecularColorButton;
    GtkWidget   **DriftShellDiffuseColorButton;
    GtkWidget   **DriftShellAmbientColorButton;
    GtkWidget   **DriftShellSpecularColorButton;
    GtkWidget   **FieldLineShininessButton;
    gulong      *FieldLineShininessButtonHandler;
    GtkWidget   **DriftShellShininessButton;
    gulong      *DriftShellShininessButtonHandler;

    GtkWidget   **FieldLineMaterialButton;
    gulong      *FieldLineMaterialButtonHandler;

    GtkWidget   **DriftShellMaterialButton;
    gulong      *DriftShellMaterialButtonHandler;

    int         nInterpedFieldLinePnts;
    int         nInterpedDriftShellPnts;

} GuiInfo;

GuiInfo *gInfo;


static MaterialProp mat_emerald_glass = {
  {0.0215, 0.1745, 0.0215, 0.5},
  {0.07568, 0.61424, 0.07568, 1.0},
  {0.633, 0.727811, 0.633, 1.0},
  0.6
};

static MaterialProp mat_emerald = {
  {0.0215, 0.1745, 0.0215, 1.0},
  {0.07568, 0.61424, 0.07568, 1.0},
  {0.633, 0.727811, 0.633, 1.0},
  0.6
};

static MaterialProp mat_jade = {
  {0.135, 0.2225, 0.1575, 1.0},
  {0.54, 0.89, 0.63, 1.0},
  {0.316228, 0.316228, 0.316228, 1.0},
  0.1
};

static MaterialProp mat_obsidian = {
  {0.05375, 0.05, 0.06625, 1.0},
  {0.18275, 0.17, 0.22525, 1.0},
  {0.332741, 0.328634, 0.346435, 1.0},
  0.3
};

static MaterialProp mat_pearl = {
  {0.25, 0.20725, 0.20725, 1.0},
  {1.0, 0.829, 0.829, 1.0},
  {0.296648, 0.296648, 0.296648, 1.0},
  0.088
};

static MaterialProp mat_ruby = {
  {0.1745, 0.01175, 0.01175, 1.0},
  {0.61424, 0.04136, 0.04136, 1.0},
  {0.727811, 0.626959, 0.626959, 1.0},
  0.6
};

static MaterialProp mat_turquoise = {
  {0.1, 0.18725, 0.1745, 1.0},
  {0.396, 0.74151, 0.69102, 1.0},
  {0.297254, 0.30829, 0.306678, 1.0},
  0.1
};

static MaterialProp mat_brass = {
  {0.329412, 0.223529, 0.027451, 1.0},
  {0.780392, 0.568627, 0.113725, 1.0},
  {0.992157, 0.941176, 0.807843, 1.0},
  0.21794872
};

static MaterialProp mat_bronze = {
  {0.2125, 0.1275, 0.054, 1.0},
  {0.714, 0.4284, 0.18144, 1.0},
  {0.393548, 0.271906, 0.166721, 1.0},
  0.2
};

static MaterialProp mat_chrome = {
  {0.25, 0.25, 0.25, 1.0},
  {0.4, 0.4, 0.4, 1.0},
  {0.774597, 0.774597, 0.774597, 1.0},
  0.6
};

static MaterialProp mat_copper = {
  {0.19125, 0.0735, 0.0225, 1.0},
  {0.7038, 0.27048, 0.0828, 1.0},
  {0.256777, 0.137622, 0.086014, 1.0},
  0.1
};

static MaterialProp mat_gold = {
  {0.24725, 0.1995, 0.0745, 1.0},
  {0.75164, 0.60648, 0.22648, 1.0},
  {0.628281, 0.555802, 0.366065, 1.0},
  0.4
};

static MaterialProp mat_silver = {
  {0.19225, 0.19225, 0.19225, 1.0},
  {0.50754, 0.50754, 0.50754, 1.0},
  {0.508273, 0.508273, 0.508273, 1.0},
  0.4
};

static MaterialProp mat_earth2 = {
  {0.2, 0.2, 0.2, 1.0},
  {3.0, 3.0, 3.0, 1.0},
  {0.3, 0.3, 0.3, 1.0},
  0.25
};
static MaterialProp mat_EqPlane = {
  {1.0, 1.0, 1.0, 0.7},
  {1.0, 1.0, 1.0, 0.7},
  {0.0, 0.0, 0.0, 1.0},
  0.25
};
static MaterialProp mat_MeridPlane1 = {
  {1.0, 1.0, 1.0, 0.7},
  {1.0, 1.0, 1.0, 0.7},
  {0.0, 0.0, 0.0, 1.0},
  0.25
};


static MaterialProp mat_red_plastic = {
  {0.1, 0.06, 0.0, 1.0},
  {0.7038, 0.12, 0.0, 1.0},
  {0.50196078, 0.50196078, 0.50196078, 1.0},
  0.25
};

static MaterialProp mat_orange_plastic = {
  {0.1, 0.06, 0.0, 1.0},
  { 0.96, 0.46, 0.18, 1.0 },
  {0.50196078, 0.50196078, 0.50196078, 1.0},
  0.25
};

static MaterialProp mat_blue_plastic = {
  {0.0, 0.06, 0.1, 1.0},
  {0.0, 0.34, 0.7038, 1.0},
  {0.50196078, 0.50196078, 0.50196078, 1.0},
  0.25
};

static MaterialProp mat_green_plastic = {
  {0.0, 0.1, 0.06, 1.0},
  {0.0, 0.7038, 0.21, 1.0},
  {0.50196078, 0.50196078, 0.50196078, 1.0},
  0.25
};

static MaterialProp mat_yellow_plastic = {
  {0.1, 0.1, 0.0, 1.0},
  {0.65, 0.65, 0.05, 1.0},
  {0.50196078, 0.50196078, 0.50196078, 1.0},
  0.25
};

static MaterialProp mat_orange_clear_plastic = {
  {0.0, 0.0, 0.0, 1.0},
  {1.0, 0.7, 0.12, 0.4},
  {0.50196078, 0.50196078, 0.50196078, 1.0},
  0.25
};



static MaterialProp *mat_current = &mat_silver;

static float view_quat_diff[4]   = { 0.0, 0.0, 0.0, 1.0 };
static float view_quat_diff3[4]  = { 0.0, 0.0, 0.0, 1.0 };
static float view_quat[4]        = { 0.0, 0.0, 0.0, 1.0 };
static float view_quat3[4]       = { 0.0, 0.0, 0.0, 1.0 };
static float view_scale          = 4.0;


int IdleRunning = FALSE;
int AnimateView = TRUE;
int AnimateTime = TIME_REALTIMEPLAY;

static void toggle_animation (GtkWidget *widget);
static void AdjustIdle( GtkWidget *widget );

void MakeFullScreenQuad() {
    g_quadDisplayList = glGenLists(1);
    glNewList(g_quadDisplayList, GL_COMPILE);

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D(0.0, 1.0, 0.0, 1.0);
    glBegin(GL_QUADS);
    {
        glVertex2f(0.0, 0.0);
        glVertex2f(1.0, 0.0);
        glVertex2f(1.0, 1.0);
        glVertex2f(0.0, 1.0);
    }
    glEnd();
    glPopMatrix();

    glEndList();
}




void LoadTextures(){

    /* The image for the texture */
    char    Filename[1024];
    GLubyte *pImage;
    int     Width;
    int     Height;
    int     PngImageOK;





    /*
     *  Texture for the Earth Image
     */
    glGenTextures( 1, &Texture_Earth );
    glBindTexture( GL_TEXTURE_2D, Texture_Earth );





    /*
     * Read Filename
     */
    sprintf( MapImageFilename, "%s/world.topo.bathy.200406.3x5400x2700.png", LGM_BLUEMARBLE_DATA_DIR );
    //strcpy( MapImageFilename, "/data1/mgh/BlueMarble/5400x2700/world.topo.bathy.200406.3x5400x2700.png");
    PngImageOK = ReadPng( MapImageFilename, &Width, &Height, &pImage );
    printf("PNG image %s: Width, Height = %d %d\n", MapImageFilename, Width, Height );

    if ( PngImageOK ) {

	    /* Put the image into texture memory */
	    glTexImage2D( GL_TEXTURE_2D,    /* texture dim and whether or not a proxy */
		      0,                    /* texture level. Always 0 unless mipmapping */
		      GL_RGBA,              /* suggested internal format */
		      Width, Height,        /* width, height of the source image */
		      0,                    /* border. Not covered in this course */
		      GL_RGBA,              /* source image format */
		      GL_UNSIGNED_BYTE,     /* source image data type */
		      pImage                /* pointer to source image */ );

	    /*
	     * Build a sequence of MIPMAPS to reduce artifatcs in rendering the texture.
	     */
	    gluBuild2DMipmaps( GL_TEXTURE_2D, GL_RGBA, Width, Height, GL_RGBA, GL_UNSIGNED_BYTE, pImage );
	    free( pImage );

	    /*
	     * Set Minification filter to:
	     *      GL_NEAREST                  nearest neighbor when not using mipmaps
	     *      GL_LINEAR                   linear interpolation when not using mipmaps
	     *      GL_NEAREST_MIPMAP_NEAREST   nearest neighbor in nearest mipmap level
	     *      GL_NEAREST_MIPMAP_LINEAR    lin. interp in nearest mipmap level
	     *      GL_LINEAR_MIPMAP_NEAREST    nearest neighbor after lin. interp between mipmap levels
	     *      GL_LINEAR_MIPMAP_LINEAR     lin interp after lin. interp between mipmap levels
	     */
	    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR );
	    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
	    glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE );
	    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);

    }


    /*
     *  Texture for the TopSide Image
     */
//    glGenTextures( 1, &Texture_TopSide );
//    glBindTexture( GL_TEXTURE_2D, Texture_TopSide );
//    strcpy( Filename, "/home/mgh/IMPACT/GITM/Image_Lat_Versus_Lon_Top_035.png");
//    ReadPng( Filename, &Width, &Height, &pImage );
//    printf("PNG image %s: Width, Height = %d %d\n", Filename, Width, Height );
//    gluBuild2DMipmaps(GL_TEXTURE_2D, GL_RGBA, Width, Height, GL_RGBA, GL_UNSIGNED_BYTE, pImage);
//    free( pImage );
//    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR );
//    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
//    //glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE );
//    glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE );
//    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);

    /*
     *  Texture for the Meridional plane 1
     */
//    glGenTextures( 1, &Texture_MeridPlane1 );
//    glBindTexture( GL_TEXTURE_2D, Texture_MeridPlane1 );
//    strcpy( Filename, "/home/mgh/IMPACT/GITM/Image_Merid1_035.png");
//    ReadPng( Filename, &Width, &Height, &pImage );
//    printf("PNG image %s: Width, Height = %d %d\n", Filename, Width, Height );
//    gluBuild2DMipmaps(GL_TEXTURE_2D, GL_RGBA, Width, Height, GL_RGBA, GL_UNSIGNED_BYTE, pImage);
//    free( pImage );
//    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR );
//    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR_MIPMAP_LINEAR );
//    //glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE );
//    glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE );
//    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);

    /*
     *  Texture for the Meridional plane 2
     */
//    glGenTextures( 1, &Texture_MeridPlane2 );
//    glBindTexture( GL_TEXTURE_2D, Texture_MeridPlane2 );
//    strcpy( Filename, "/home/mgh/IMPACT/GITM/Image_Merid2_035.png");
//    ReadPng( Filename, &Width, &Height, &pImage );
//    printf("PNG image %s: Width, Height = %d %d\n", Filename, Width, Height );
//    gluBuild2DMipmaps(GL_TEXTURE_2D, GL_RGBA, Width, Height, GL_RGBA, GL_UNSIGNED_BYTE, pImage);
//    free( pImage );
//    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR );
//    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR_MIPMAP_LINEAR );
//    //glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE );
//    glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE );
//    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);

    /*
     *  Texture for the EQ plane
     */
//    glGenTextures( 1, &Texture_EqPlane );
//    glBindTexture( GL_TEXTURE_2D, Texture_EqPlane );
//    strcpy( Filename, "/home/mgh/DREAM/Dream/Dream/Images/checkerboard_lg.png" );
//    ReadPng( Filename, &Width, &Height, &pImage );
//    printf("PNG image %s: Width, Height = %d %d\n", Filename, Width, Height );
//    gluBuild2DMipmaps(GL_TEXTURE_2D, GL_RGBA, Width, Height, GL_RGBA, GL_UNSIGNED_BYTE, pImage);
//    free( pImage );
//    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR );
//    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
//    //glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE );
//    glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE );
//    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);


    /*
     *  Texture for Debris sprites
     */
    glGenTextures( 1, &Texture_Debris );
    glBindTexture( GL_TEXTURE_2D, Texture_Debris );

    strcpy( Filename, "/home/mgh/DREAM/Dream/Dream/Images/HexNut.png");
    if ( ReadPng( Filename, &Width, &Height, &pImage ) ) {
        printf("PNG image %s: Width, Height = %d %d\n", Filename, Width, Height );
        gluBuild2DMipmaps(GL_TEXTURE_2D, GL_RGBA, Width, Height, GL_RGBA, GL_UNSIGNED_BYTE, pImage);
        free( pImage );
        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR );
        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
        //glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE );
        glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE );
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    }


    glGenTextures( 1, &Texture_RocketBody );
    glBindTexture( GL_TEXTURE_2D, Texture_RocketBody );

    strcpy( Filename, "/home/mgh/DREAM/Dream/Dream/Images/RocketBody4.png");
    if ( ReadPng( Filename, &Width, &Height, &pImage ) ) {
        printf("PNG image %s: Width, Height = %d %d\n", Filename, Width, Height );
        gluBuild2DMipmaps(GL_TEXTURE_2D, GL_RGBA, Width, Height, GL_RGBA, GL_UNSIGNED_BYTE, pImage);
        free( pImage );
        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR );
        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
        //glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE );
        glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE );
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    }

    glGenTextures( 1, &Texture_Spacecraft );
    glBindTexture( GL_TEXTURE_2D, Texture_Spacecraft );
    strcpy( Filename, "/home/mgh/DREAM/Dream/Dream/Images/spacecraft2.png");
    if ( ReadPng( Filename, &Width, &Height, &pImage ) ) {
        printf("PNG image %s: Width, Height = %d %d\n", Filename, Width, Height );
        gluBuild2DMipmaps(GL_TEXTURE_2D, GL_RGBA, Width, Height, GL_RGBA, GL_UNSIGNED_BYTE, pImage);
        free( pImage );
        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR );
        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
        //glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE );
        glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE );
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    }


    /*
     * Sun Texture
     */
    glGenTextures( 1, &Texture_Sun );
    glBindTexture( GL_TEXTURE_2D, Texture_Sun );
    strcpy( Filename, "/home/mgh/DREAM/Dream/Dream/Images/happy.png");
    if ( ReadPng( Filename, &Width, &Height, &pImage ) ) {
        printf("PNG image %s: Width, Height = %d %d\n", Filename, Width, Height );
        gluBuild2DMipmaps(GL_TEXTURE_2D, GL_RGBA, Width, Height, GL_RGBA, GL_UNSIGNED_BYTE, pImage);
        free( pImage );
        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR );
        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
        //glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE );
        glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE );
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    }


    /*
     *  Texture for Moon
     */
    glGenTextures( 1, &Texture_Moon );
    glBindTexture( GL_TEXTURE_2D, Texture_Moon );
    strcpy( Filename, "/home/mgh/DREAM/Dream/Dream/Images/Moon.png");
    if ( ReadPng( Filename, &Width, &Height, &pImage ) ) {
        printf("PNG image %s: Width, Height = %d %d\n", Filename, Width, Height );
        glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA, Width, Height, 0, GL_RGBA, GL_UNSIGNED_BYTE, pImage);
        gluBuild2DMipmaps(GL_TEXTURE_2D, GL_RGBA, Width, Height, GL_RGBA, GL_UNSIGNED_BYTE, pImage);
        free( pImage );
        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR );
        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
        glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE );
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    }


    /*
     *  Texture for Logo
     */
    glGenTextures( 1, &Texture_Logo );
    glBindTexture( GL_TEXTURE_2D, Texture_Logo );
    strcpy( Filename, "/home/mgh/DREAM/Dream/Dream/Images/LANL_LOGO_REV_ALPHA_50_V2.png");
    if ( ReadPng( Filename, &Logo_Width, &Logo_Height, &pImage ) ) {
        printf("PNG image %s: Width, Height = %d %d\n", Filename, Logo_Width, Logo_Height );
        glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA, Logo_Width, Logo_Height, 0, GL_RGBA, GL_UNSIGNED_BYTE, pImage);
        gluBuild2DMipmaps(GL_TEXTURE_2D, GL_RGBA, Logo_Width, Logo_Height, GL_RGBA, GL_UNSIGNED_BYTE, pImage);
        free( pImage );
        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR );
        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
        glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE );
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    }


    /*
     * For Cairo rendering
     */
    glGenTextures( 1, &Texture_TEMP );

}




void StarColor( char Type, float *Red, float *Grn, float *Blu ) {
    switch( Type ) {
        case 'W':
        case 'O': *Red = 170.0/255.0; *Grn = 191.0/255.0; *Blu = 255.0/255.0; break;
        case 'B': *Red = 202.0/255.0; *Grn = 215.0/255.0; *Blu = 255.0/255.0; break;
        case 'A': *Red = 226.0/255.0; *Grn = 223.0/255.0; *Blu = 255.0/255.0; break;
        case 'F': *Red = 255.0/255.0; *Grn = 225.0/255.0; *Blu = 255.0/255.0; break;
        case 'G': *Red = 249.0/255.0; *Grn = 246.0/255.0; *Blu = 191.0/255.0; break;
        case 'K': *Red = 255.0/255.0; *Grn = 228.0/255.0; *Blu = 111.0/255.0; break;
        case 'C':
        case 'R':
        case 'N':
        case 'L':
        case 'T':
        case 'Y':
        case 'M': *Red = 255.0/255.0; *Grn = 160.0/255.0; *Blu =  64.0/255.0; break;
        default:
            *Red = 170.0/255.0; *Grn = 191.0/255.0; *Blu = 255.0/255.0; break;
    }
    return;
}

_GLUfuncptr errorCallback(GLenum errorCode) {
   const GLubyte *estring;

   estring = gluErrorString(errorCode);
   fprintf(stderr, "Quadric Error: %s\n", estring);
   exit(0);
}


static MaterialProp mat_blue_trans = {
  {0.0, 0.06, 0.1, 0.5},
  {0.0, 0.34, 0.7038, 0.5},
  {0.50196078, 0.50196078, 0.50196078, 0.5},
  0.25
};

//void SolidCone(GLdouble base, GLdouble height, GLint slices, GLint stacks) {
void SolidCone( Lgm_Vector *u, double Fov, GLint slices, GLint stacks) {

    double      RotAng, height, base;
    Lgm_Vector  v, z, RotAxis;

printf("u = %g %g %g\n", u->x, u->y, u->z);
    v = *u;

    z.x = 0.0; z.y = 0.0; z.z = 1.0;
printf("z = %g %g %g\n", z.x, z.y, z.z);

    height = Lgm_NormalizeVector( &v );
    base   = height*tan( RadPerDeg*0.5*Fov );
    RotAng = -DegPerRad*acos( v.z );
    Lgm_CrossProduct( &v, &z, &RotAxis );
    Lgm_NormalizeVector( &RotAxis );

printf("height, base, RotAng = %g %g %g\n", height, base, RotAng);

//        glPushMatrix();
//        gluCylinder(qobj, 0.08*AxesScale, 0.08*AxesScale, 2.0, 15, 5);
//        glTranslatef( 0.0, 0.0, 2.0 );
//        Lgm_gl_draw_cone( TRUE, 0.14*AxesScale, 0.5*AxesScale, 15, 5 );

//        glPopMatrix();

//    glBegin(GL_LINE_LOOP);
    GLUquadricObj* quadric = gluNewQuadric();
        gluQuadricCallback(quadric, GLU_ERROR, (_GLUfuncptr) errorCallback);
        gluQuadricDrawStyle(quadric, GLU_FILL); /* smooth shaded */
        gluQuadricNormals(quadric, GLU_SMOOTH);
//    gluQuadricDrawStyle(quadric, GLU_LINE);
        glMaterialfv( GL_FRONT, GL_AMBIENT,   mat_blue_trans.ambient);
        glMaterialfv( GL_FRONT, GL_DIFFUSE,   mat_blue_trans.diffuse);
        glMaterialfv( GL_FRONT, GL_SPECULAR,  mat_blue_trans.specular);
        glMaterialf(  GL_FRONT, GL_SHININESS, mat_blue_trans.shininess * 128.0);



    glPushMatrix();

    glTranslatef( 0.0, 0.0, 1.0 );
    glRotatef( RotAng, RotAxis.x, RotAxis.y, RotAxis.z );

    gluCylinder(quadric, base, 0, height, slices, stacks );

    glPopMatrix();

//    gluDeleteQuadric(quadric);
//    glEnd();
}

void CreateZPSAxes( ) {

    GLUquadricObj   *qobj;

    // ZPS instrument Coord Axes
    ZPSAxesDL = glGenLists( 1 );
    glNewList( ZPSAxesDL, GL_COMPILE );
    {
        qobj = gluNewQuadric();
        gluQuadricCallback(qobj, GLU_ERROR, (_GLUfuncptr) errorCallback);
        gluQuadricDrawStyle(qobj, GLU_FILL); /* smooth shaded */
        gluQuadricNormals(qobj, GLU_SMOOTH);

// add this to options
double AxesScale = 0.5;

        // z-axis
        glMaterialfv( GL_FRONT, GL_AMBIENT,   mat_blue_plastic.ambient);
        glMaterialfv( GL_FRONT, GL_DIFFUSE,   mat_blue_plastic.diffuse);
        glMaterialfv( GL_FRONT, GL_SPECULAR,  mat_blue_plastic.specular);
        glMaterialf(  GL_FRONT, GL_SHININESS, mat_blue_plastic.shininess * 128.0);
        glPushMatrix();
        gluCylinder(qobj, 0.08*AxesScale, 0.08*AxesScale, 2.0, 15, 5);
        glTranslatef( 0.0, 0.0, 2.0 );
        //Lgm_gl_draw_cone( TRUE, 0.14*AxesScale, 0.5*AxesScale, 15, 5 );
        Lgm_gl_draw_cone( TRUE, 0.14*AxesScale, 0.5*AxesScale, 15, 5 );

        glPopMatrix();


        // x-axis
        glMaterialfv( GL_FRONT, GL_AMBIENT,   mat_red_plastic.ambient);
        glMaterialfv( GL_FRONT, GL_DIFFUSE,   mat_red_plastic.diffuse);
        glMaterialfv( GL_FRONT, GL_SPECULAR,  mat_red_plastic.specular);
        glMaterialf(  GL_FRONT, GL_SHININESS, mat_red_plastic.shininess * 128.0);
        glPushMatrix();
        glRotatef( 90.0, 0.0, 1.0, 0.0);
        gluCylinder(qobj, 0.08*AxesScale, 0.08*AxesScale, 2.0, 15, 5);
        glTranslatef( 0.0, 0.0, 2.0 );
        Lgm_gl_draw_cone( TRUE, 0.14*AxesScale, 0.5*AxesScale, 15, 5 );
        glPopMatrix();


        // y-axis
        glMaterialfv( GL_FRONT, GL_AMBIENT,   mat_green_plastic.ambient);
        glMaterialfv( GL_FRONT, GL_DIFFUSE,   mat_green_plastic.diffuse);
        glMaterialfv( GL_FRONT, GL_SPECULAR,  mat_green_plastic.specular);
        glMaterialf(  GL_FRONT, GL_SHININESS, mat_green_plastic.shininess * 128.0);
        glPushMatrix();
        glRotatef( -90.0, 1.0, 0.0, 0.0);
        gluCylinder(qobj, 0.08*AxesScale, 0.08*AxesScale, 2.0, 15, 5);
        glTranslatef( 0.0, 0.0, 2.0 );
        Lgm_gl_draw_cone( TRUE, 0.14*AxesScale, 0.5*AxesScale, 15, 5 );
        glPopMatrix();
    }
    glEndList( );

}




void CreateGSMAxes( ) {

    GLUquadricObj   *qobj;

    // GSM Coord Axes
    AxesDL = glGenLists( 1 );
    glNewList( AxesDL, GL_COMPILE );
    {
        qobj = gluNewQuadric();
        gluQuadricCallback(qobj, GLU_ERROR, (_GLUfuncptr) errorCallback);
        gluQuadricDrawStyle(qobj, GLU_FILL); /* smooth shaded */
        gluQuadricNormals(qobj, GLU_SMOOTH);

// add this to options
double AxesScale = 0.5;

        // z-axis
        glMaterialfv( GL_FRONT, GL_AMBIENT,   mat_blue_plastic.ambient);
        glMaterialfv( GL_FRONT, GL_DIFFUSE,   mat_blue_plastic.diffuse);
        glMaterialfv( GL_FRONT, GL_SPECULAR,  mat_blue_plastic.specular);
        glMaterialf(  GL_FRONT, GL_SHININESS, mat_blue_plastic.shininess * 128.0);
        glPushMatrix();
        gluCylinder(qobj, 0.08*AxesScale, 0.08*AxesScale, 2.0, 15, 5);
        glTranslatef( 0.0, 0.0, 2.0 );
        Lgm_gl_draw_cone( TRUE, 0.14*AxesScale, 0.5*AxesScale, 15, 5 );

        glPopMatrix();


        // x-axis
        glMaterialfv( GL_FRONT, GL_AMBIENT,   mat_red_plastic.ambient);
        glMaterialfv( GL_FRONT, GL_DIFFUSE,   mat_red_plastic.diffuse);
        glMaterialfv( GL_FRONT, GL_SPECULAR,  mat_red_plastic.specular);
        glMaterialf(  GL_FRONT, GL_SHININESS, mat_red_plastic.shininess * 128.0);
        glPushMatrix();
        glRotatef( 90.0, 0.0, 1.0, 0.0);
        gluCylinder(qobj, 0.08*AxesScale, 0.08*AxesScale, 2.0, 15, 5);
        glTranslatef( 0.0, 0.0, 2.0 );
        Lgm_gl_draw_cone( TRUE, 0.14*AxesScale, 0.5*AxesScale, 15, 5 );
        glPopMatrix();


        // y-axis
        glMaterialfv( GL_FRONT, GL_AMBIENT,   mat_green_plastic.ambient);
        glMaterialfv( GL_FRONT, GL_DIFFUSE,   mat_green_plastic.diffuse);
        glMaterialfv( GL_FRONT, GL_SPECULAR,  mat_green_plastic.specular);
        glMaterialf(  GL_FRONT, GL_SHININESS, mat_green_plastic.shininess * 128.0);
        glPushMatrix();
        glRotatef( -90.0, 1.0, 0.0, 0.0);
        gluCylinder(qobj, 0.08*AxesScale, 0.08*AxesScale, 2.0, 15, 5);
        glTranslatef( 0.0, 0.0, 2.0 );
        Lgm_gl_draw_cone( TRUE, 0.14*AxesScale, 0.5*AxesScale, 15, 5 );
        glPopMatrix();
    }
    glEndList( );


}


void CreateEarth( ){

    GLUquadricObj       *qobj;

    // Earth
    EarthDL = glGenLists( 1 );
    glNewList( EarthDL, GL_COMPILE );
    {
        glPushMatrix();
        glMaterialfv( GL_FRONT, GL_AMBIENT,   mat_earth2.ambient);
        glMaterialfv( GL_FRONT, GL_DIFFUSE,   mat_earth2.diffuse);
        glMaterialfv( GL_FRONT, GL_SPECULAR,  mat_earth2.specular);
        glMaterialf(  GL_FRONT, GL_SHININESS, mat_earth2.shininess * 128.0);
        glBindTexture( GL_TEXTURE_2D, Texture_Earth );
        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR );
        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR_MIPMAP_LINEAR );
        glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE );
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP );
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP );
        glEnable( GL_TEXTURE_2D );
        glRotatef( RotAngle, RotAxis.x, RotAxis.y, RotAxis.z );
        glRotatef( 180.0, 0.0, 0.0, 1.0); // rotates image around so that 0deg. glon is in the +x direction
        CreateEllipsoid( 1.0, (double)(WGS84_B/WGS84_A), 80 );
        glDisable( GL_TEXTURE_2D );
        glPopMatrix();
    }
    glEndList( );




    // Dipole Axis
    DipoleAxisDL = glGenLists( 1 );
    glNewList( DipoleAxisDL, GL_COMPILE );
    {
        qobj = gluNewQuadric();
        glMaterialfv( GL_FRONT, GL_AMBIENT,   mat_yellow_plastic.ambient);
        glMaterialfv( GL_FRONT, GL_DIFFUSE,   mat_yellow_plastic.diffuse);
        glMaterialfv( GL_FRONT, GL_SPECULAR,  mat_yellow_plastic.specular);
        glMaterialf(  GL_FRONT, GL_SHININESS, mat_yellow_plastic.shininess * 128.0);

        glPushMatrix();

        // do a rotation so that the current x,y,z coord sys corresponds to the SM system
        glRotatef( RotAngle2, RotAxis2.x, RotAxis2.y, RotAxis2.z );

        // This would inplement an offset dipole if we want to do that
        //glTranslatef( DipoleOffset_sm.x, DipoleOffset_sm.y, DipoleOffset_sm.z );

        // Then draw a cylinder+cone to represent the vector
        glTranslatef( 0.0, 0.0, -1.7 );
        gluCylinder(qobj, 0.041, 0.041, 3.2, 15, 5);
        glTranslatef( 0.0, 0.0, 3.2 );
        Lgm_gl_draw_cone( TRUE, 0.07, 0.4, 15, 5 );

        glPopMatrix();
    }
    glEndList( );



    // Sun vector (from earth to Sun)
    SunDirectionDL = glGenLists( 1 );
    glNewList( SunDirectionDL, GL_COMPILE );
    {
        qobj = gluNewQuadric();
        glMaterialfv( GL_FRONT, GL_AMBIENT,   mat_orange_plastic.ambient);
        glMaterialfv( GL_FRONT, GL_DIFFUSE,   mat_orange_plastic.diffuse);
        glMaterialfv( GL_FRONT, GL_SPECULAR,  mat_orange_plastic.specular);
        glMaterialf(  GL_FRONT, GL_SHININESS, mat_orange_plastic.shininess * 128.0);
        glPushMatrix();
        glRotatef( RotAngle3, RotAxis3.x, RotAxis3.y, RotAxis3.z );
        glRotatef( 90.0, 0.0, 1.0, 0.0 );
        glTranslatef( 0.0, 0.0, -1.7 );
        gluCylinder(qobj, 0.041, 0.041, 3.2, 15, 5);
        glTranslatef( 0.0, 0.0, 3.2 );
        Lgm_gl_draw_cone( TRUE, 0.07, 0.4, 15, 5 );
        glPopMatrix();
    }
    glEndList( );



    // EqPlane grid
    double R, x, y, Phi;
    EqPlaneGridDL = glGenLists( 1 );
    glNewList( EqPlaneGridDL, GL_COMPILE );
        glDisable(GL_LIGHTING);
        glPushMatrix();
        glColor4f( 0.5, 0.5, 0.5, 0.1 );
        glLineWidth( 0.5 );
        glEnable(GL_LINE_SMOOTH);
        glHint( GL_LINE_SMOOTH_HINT, GL_NICEST);
        for (R=1.0; R<=10.0; R += 1.0) {
            glBegin(GL_LINE_LOOP);
                for (Phi=0.0; Phi<2.0*M_PI; Phi += .01) {
                    x = R*cos(Phi);
                    y = R*sin(Phi);
                    glVertex3f( x, y, 0.0 );
                }
            glEnd();
        }
        glBegin(GL_LINES);
            for (Phi=0.0; Phi<360.0; Phi += 15.0) {
                x = cos(Phi*M_PI/180.0);
                y = sin(Phi*M_PI/180.0);
                glVertex3f( x, y, 0.0 );
                glVertex3f( 10.0*x, 10.0*y, 0.0 );
            }
        glEnd();
        glDisable(GL_LINE_SMOOTH);
        glPopMatrix();
        glEnable(GL_LIGHTING);
    glEndList( );




    // EqPlane image
    EqPlaneDL = glGenLists( 1 );
    glNewList( EqPlaneDL, GL_COMPILE );
        glMaterialfv( GL_FRONT, GL_AMBIENT,   mat_earth2.ambient);
        glMaterialfv( GL_FRONT, GL_DIFFUSE,   mat_earth2.diffuse);
        glMaterialfv( GL_FRONT, GL_SPECULAR,  mat_earth2.specular);
        glMaterialf(  GL_FRONT, GL_SHININESS, mat_earth2.shininess * 128.0);
        glBindTexture( GL_TEXTURE_2D, Texture_EqPlane );
        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR );
        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR_MIPMAP_LINEAR );
        glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE );
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP );
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP );
        glEnable( GL_TEXTURE_2D );
//        glEnable( GL_BLEND );
        glDepthMask( GL_FALSE );
//        glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
        glPushMatrix();
        // do a rotation so that the current x,y,z coord sys corresponds to the SM system
        glRotatef( RotAngle2, RotAxis2.x, RotAxis2.y, RotAxis2.z );
//        glRotatef( DipoleTiltAngle, 0.0, 1.0, 0.0);
        glBegin(GL_QUADS);
            glNormal3f( 0.0, 0.0, 1.0 ); glTexCoord2f(0.0, 0.0); glVertex3f( -10.0,  10.0, 0.0 );
            glNormal3f( 0.0, 0.0, 1.0 ); glTexCoord2f(1.0, 0.0); glVertex3f( -10.0, -10.0, 0.0 );
            glNormal3f( 0.0, 0.0, 1.0 ); glTexCoord2f(1.0, 1.0); glVertex3f(  10.0, -10.0, 0.0 );
            glNormal3f( 0.0, 0.0, 1.0 ); glTexCoord2f(0.0, 1.0); glVertex3f(  10.0,  10.0, 0.0 );
        glEnd();
        glPopMatrix();
        glDepthMask( GL_TRUE );
//        glDisable( GL_BLEND );
        glDisable( GL_TEXTURE_2D );
    glEndList( );




if (0==1){
    // TopSide
    TopSideDL = glGenLists( 1 );
    glNewList( TopSideDL, GL_COMPILE );
    {
        glPushMatrix();
        glMaterialfv( GL_FRONT, GL_AMBIENT,   mat_earth2.ambient);
        glMaterialfv( GL_FRONT, GL_DIFFUSE,   mat_earth2.diffuse);
        glMaterialfv( GL_FRONT, GL_SPECULAR,  mat_earth2.specular);
        glMaterialf(  GL_FRONT, GL_SHININESS, mat_earth2.shininess * 128.0);
        glBindTexture( GL_TEXTURE_2D, Texture_TopSide );
        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR );
        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR_MIPMAP_LINEAR );
        glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE );
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP );
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP );
        glEnable( GL_TEXTURE_2D );
        glRotatef( RotAngle, RotAxis.x, RotAxis.y, RotAxis.z );
        glRotatef( 180.0, 0.0, 0.0, 1.0); // rotates image around so that 0deg. glon is in the +x direction
        //CreateEllipsoid( 1.1, (double)(WGS84_B/WGS84_A), 80 );

        CreateCutEllipsoid( (1.0+500.0/WGS84_A), (double)(1.0+500.0/WGS84_A), 76*2 );
        glDisable( GL_TEXTURE_2D );
        glPopMatrix();
    }
    glEndList( );

    // MeridPlane1 image
    MeridPlane1DL = glGenLists( 1 );
    glNewList( MeridPlane1DL, GL_COMPILE );
        glMaterialfv( GL_FRONT, GL_AMBIENT,   mat_earth2.ambient);
        glMaterialfv( GL_FRONT, GL_DIFFUSE,   mat_earth2.diffuse);
        glMaterialfv( GL_FRONT, GL_SPECULAR,  mat_earth2.specular);
        glMaterialf(  GL_FRONT, GL_SHININESS, mat_earth2.shininess * 128.0);
        glBindTexture( GL_TEXTURE_2D, Texture_MeridPlane1 );
        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR );
        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR_MIPMAP_LINEAR );
        glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE );
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP );
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP );
        glShadeModel(GL_SMOOTH);
        glEnable( GL_TEXTURE_2D );
        glEnable( GL_BLEND );
//        glDisable( GL_BLEND );
        //glDepthMask( GL_FALSE );
        //glDepthMask( GL_TRUE );
        glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
        glPushMatrix();
        glRotatef( RotAngle, RotAxis.x, RotAxis.y, RotAxis.z );
//        glRotatef( 180.0, 0.0, 0.0, 1.0); // rotates image around so that 0deg. glon is in the +x direction
        glRotatef( 180.0, 0.0, 0.0, 1.0);
        glRotatef( 90.0, 0.0, 0.0, 1.0);
//        glRotatef( 90.0, 1.0, 0.0, 0.0);
        glBegin(GL_QUADS);
            glNormal3f( 1.0, 0.0, 0.0 ); glTexCoord2f(0.0, 0.0);  glVertex3f(  0.0,  0.0, -1.1 );
            glNormal3f( 1.0, 0.0, 0.0 ); glTexCoord2f(1.0, 0.0); glVertex3f( 1.1, 0.0, -1.1 );
            glNormal3f( 1.0, 0.0, 0.0 ); glTexCoord2f(1.0, 1.0); glVertex3f(  1.1, 0.0, 1.1 );
            glNormal3f( 1.0, 0.0, 0.0 ); glTexCoord2f(0.0, 1.0);  glVertex3f(   0.0,  0.0, 1.1 );
        glEnd();
        glPopMatrix();
        glDepthMask( GL_TRUE );
//        glDisable( GL_BLEND );
        glDisable( GL_TEXTURE_2D );
    glEndList( );

    // MeridPlane1 image 2
    MeridPlane2DL = glGenLists( 1 );
    glNewList( MeridPlane2DL, GL_COMPILE );
        glMaterialfv( GL_FRONT, GL_AMBIENT,   mat_earth2.ambient);
        glMaterialfv( GL_FRONT, GL_DIFFUSE,   mat_earth2.diffuse);
        glMaterialfv( GL_FRONT, GL_SPECULAR,  mat_earth2.specular);
        glMaterialf(  GL_FRONT, GL_SHININESS, mat_earth2.shininess * 128.0);
        glBindTexture( GL_TEXTURE_2D, Texture_MeridPlane2 );
        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR );
        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR_MIPMAP_LINEAR );
        glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE );
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP );
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP );
        glShadeModel(GL_SMOOTH);
        glEnable( GL_TEXTURE_2D );
        glEnable( GL_BLEND );
//        glDisable( GL_BLEND );
        //glDepthMask( GL_FALSE );
        //glDepthMask( GL_TRUE );
        glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
        glPushMatrix();
        glRotatef( RotAngle, RotAxis.x, RotAxis.y, RotAxis.z );
//        glRotatef( 180.0, 0.0, 0.0, 1.0); // rotates image around so that 0deg. glon is in the +x direction
        glRotatef( 180.0, 0.0, 0.0, 1.0);
        glRotatef( 90.0, 0.0, 0.0, 1.0);
        glRotatef( 270.0, 0.0, 0.0, 1.0);
//        glRotatef( 90.0, 1.0, 0.0, 0.0);
        glBegin(GL_QUADS);
            glNormal3f( 1.0, 0.0, 0.0 ); glTexCoord2f(0.0, 0.0);  glVertex3f(  0.0,  0.0, -1.1 );
            glNormal3f( 1.0, 0.0, 0.0 ); glTexCoord2f(1.0, 0.0); glVertex3f( 1.1, 0.0, -1.1 );
            glNormal3f( 1.0, 0.0, 0.0 ); glTexCoord2f(1.0, 1.0); glVertex3f(  1.1, 0.0, 1.1 );
            glNormal3f( 1.0, 0.0, 0.0 ); glTexCoord2f(0.0, 1.0);  glVertex3f(   0.0,  0.0, 1.1 );
        glEnd();
        glPopMatrix();
        glDepthMask( GL_TRUE );
//        glDisable( GL_BLEND );
        glDisable( GL_TEXTURE_2D );
    glEndList( );
}


}

void DestroyHiResEarthQuad( GLuint *Texture ){

    glDeleteTextures( 1, Texture ); // destroy just the selected

    return;

}

/*
 * This shows a single hi res quad
 */
void CreateHiResEarthQuad( int LoadQuad, int nph, int nth, long int QuadId, GLuint *DL, GLuint *TexId ){

    char            Filename[1024];
    GLubyte         *pImage;
    int             Width;
    int             Height, i, j, ii, jj;
    double          th1, th2, r;
    double          st1, st2, ct1, ct2, sp1, cp1;
    double          ths, phs, dth, dph, th, ph;
    double          xi1, yi1, yi2, th_res, ph_res;
    Lgm_Vector      u;
    static GLuint   DisplayList;
    GLuint          Texture;

    r = 1.0+3.5/Re;

    // decode the QuadId
    i = QuadId/100;
    j = QuadId - i*100;
    Texture = *TexId;

    if ( LoadQuad ) {
        glGenTextures( 1, &Texture ); *TexId = Texture;
        glBindTexture( GL_TEXTURE_2D, Texture );
        if ((nph == 80 )&&(nth == 40)) {
            sprintf(Filename, "/data1/mgh/BlueMarble/86400x43200/200406/TextureAtlas_80x40/Quad_%02d_%02d.png", i, j );
        } else {
            sprintf(Filename, "/data1/mgh/BlueMarble/43200x21600/200406/TextureAtlas_40x20/Quad_%02d_%02d.png", i, j );
        }
        if ( ReadPng( Filename, &Width, &Height, &pImage ) >= 0 ){
            printf("PNG image %s: Width, Height = %d %d\n", Filename, Width, Height );
            glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA, Width, Height, 0, GL_RGBA, GL_UNSIGNED_BYTE, pImage);
            free( pImage );
        }
    }

    // Hi Res Earth Quad
    DisplayList = glGenLists( 1 ); *DL = DisplayList;
    glNewList( DisplayList, GL_COMPILE );
    {
        glPushMatrix();
        glRotatef( -RotAngle3, RotAxis3.x, RotAxis3.y, RotAxis3.z ); // This implements coord trans from Observer coords -> GSM
        glRotatef( RotAngle, RotAxis.x, RotAxis.y, RotAxis.z );
        glMaterialfv( GL_FRONT, GL_AMBIENT,   mat_earth2.ambient);
        glMaterialfv( GL_FRONT, GL_DIFFUSE,   mat_earth2.diffuse);
        glMaterialfv( GL_FRONT, GL_SPECULAR,  mat_earth2.specular);
        glMaterialf(  GL_FRONT, GL_SHININESS, mat_earth2.shininess * 128.0);
        glBindTexture( GL_TEXTURE_2D, Texture );
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
        glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE );
        glEnable( GL_TEXTURE_2D );
        th_res = 180.0/(double)nth*RadPerDeg;
        ph_res = 180.0/(double)nth*RadPerDeg;
        ths = j* th_res; 
        phs = i* ph_res-M_PI; 
        dth = th_res/10.0;
        dph = ph_res/10.0;
        for (th = ths, jj=0; jj<10; jj++){
            th1 = th; th2 = th1+dth;
            st1 = sin( th1 ); st2 = sin( th2 );
            ct1 = cos( th1 ); ct2 = cos( th2 );
            yi1 = 1.0-0.1*jj; yi2 = yi1-0.1;
            glBegin(GL_QUAD_STRIP );
                for (ph = phs, ii=0; ii<=10; ii++){
                    sp1 = sin( ph );
                    cp1 = cos( ph );
                    xi1 = 0.1*ii;
                    u.x = st1*cp1; u.y = st1*sp1; u.z = ct1;
                    glNormal3f( u.x, u.y, u.z ); glTexCoord2f(xi1, yi1); glVertex3f( r*u.x, r*u.y, r*u.z );
                    u.x = st2*cp1; u.y = st2*sp1; u.z = ct2;
                    glNormal3f( u.x, u.y, u.z ); glTexCoord2f(xi1, yi2); glVertex3f( r*u.x, r*u.y, r*u.z );
                    ph += dph;
                }
            glEnd();
            th += dth;
        }
        glDisable( GL_TEXTURE_2D );
        glPopMatrix();
    }
    glEndList( );

}


void ReCreateEarth ( ){
    glDeleteLists( EarthDL, 1 );
    glDeleteLists( DipoleAxisDL, 1 );
    glDeleteLists( SunDirectionDL, 1 );
    glDeleteLists( EqPlaneGridDL, 1 );
    glDeleteLists( MeridPlane1DL, 1 );
    glDeleteLists( MeridPlane2DL, 1 );
    glDeleteLists( EqPlaneDL, 1 );
    CreateEarth( );
}

void CreateGeoMarkers( ) {

    double  th, ph, x, y, z;

    // Add geographic markers (cities etc.)
    GeoMarkersDL = glGenLists( 1 );
    glNewList( GeoMarkersDL, GL_COMPILE );
    {
        glPushMatrix();
        glRotatef( -RotAngle3, RotAxis3.x, RotAxis3.y, RotAxis3.z ); // This implements coord trans from Observer coords -> GSM
        glRotatef( RotAngle, RotAxis.x, RotAxis.y, RotAxis.z );
        glMaterialfv( GL_FRONT, GL_AMBIENT,   mat_earth2.ambient);
        glMaterialfv( GL_FRONT, GL_DIFFUSE,   mat_earth2.diffuse);
        glMaterialfv( GL_FRONT, GL_SPECULAR,  mat_earth2.specular);
        glMaterialf(  GL_FRONT, GL_SHININESS, mat_earth2.shininess * 128.0);

        //Add here -- point sprites might be better?
        glPushMatrix();
        glTranslatef( 1.0, 0.0, 0.0);
        CreateSphere( 0.005, 20 );
        glPopMatrix();

        // Los Alamos
        glPushMatrix();
        th = (90.0-35.888)*RadPerDeg;
        ph = (-106.306)*RadPerDeg;
        x = 1.0003394*sin(th)*cos(ph); y = 1.0003394*sin(th)*sin(ph); z = 1.0003394*cos(th);
        glTranslatef( x, y, z );
        CreateSphere( 0.001, 20 );
        glPopMatrix();


        glPopMatrix();
    }
    glEndList( );

}

void CreateLogo( ){

    int     w, h;
    double  dx, dy, x1, x2, y1, y2;

    h = 100;
    w = (int)((double)Logo_Width/(double)Logo_Height * (double)h);



    dx = (double)w/(double)g_imageWidth;
    dy = (double)h/(double)g_imageHeight;

    x1 = 1.0 - 2.0*dx-.01;
    x2 = 1.0-.01;
    y1 = -1.0+0.01;
    y2 = -1.0+2.0*dy+0.01;

    // Logo
    LogoDL = glGenLists( 1 );
    glNewList( LogoDL, GL_COMPILE );
    {
        glColor4f( 1.0, 1.0, 1.0, 1.0 );
        glBindTexture( GL_TEXTURE_2D, Texture_Logo );
        glBegin(GL_QUADS);
        {
            glTexCoord2f( 0.0, 0.0 ); glVertex2f(  x1, y1 );
            glTexCoord2f( 1.0, 0.0 ); glVertex2f(  x2, y1 );
            glTexCoord2f( 1.0, 1.0 ); glVertex2f(  x2, y2 );
            glTexCoord2f( 0.0, 1.0 ); glVertex2f(  x1, y2 );
        }
        glEnd();
    }
    glEndList( );

}

void ReCreateLogo( ){
    if (LogoDL >= 0) glDeleteLists( LogoDL, 1 );
    CreateLogo();
}


void CreateMoon( ){

    Lgm_Vector  P_gei, P;


    // Earth
    MoonDL = glGenLists( 1 );
    glNewList( MoonDL, GL_COMPILE );
    {
        glMaterialfv( GL_FRONT, GL_AMBIENT,   mat_earth2.ambient);
        glMaterialfv( GL_FRONT, GL_DIFFUSE,   mat_earth2.diffuse);
        glMaterialfv( GL_FRONT, GL_SPECULAR,  mat_earth2.specular);
        glMaterialf(  GL_FRONT, GL_SHININESS, mat_earth2.shininess * 128.0);
        glBindTexture( GL_TEXTURE_2D, Texture_Moon );
        glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE );
        glEnable( GL_TEXTURE_2D );
        glPushMatrix();
            Lgm_Radec_to_Cart( mInfo->c->RA_moon, mInfo->c->DEC_moon, &P_gei );
            Lgm_ScaleVector( &P_gei, mInfo->c->EarthMoonDistance );
            Lgm_Convert_Coords( &P_gei, &P, StarsConvertFlag, mInfo->c );
            glTranslatef( P.x, P.y, P.z );

//        glRotatef( RotAngle, RotAxis.x, RotAxis.y, RotAxis.z );
//        glRotatef( 180.0, 0.0, 0.0, 1.0); // rotates image around so that 0deg. glon is in the +x direction
        CreateSphere( 3.0*0.2724, 80 );
        glPopMatrix();
        glDisable( GL_TEXTURE_2D );
    }
    glEndList( );


}

void ReCreateMoon( ){
    glDeleteLists( MoonDL, 1 );
    CreateMoon();
}


void CreateSCPos( ) {

    // S/C position -- since this is transparent we want to call it last.
    ScPositionDL = glGenLists( 1 );
    glNewList( ScPositionDL, GL_COMPILE );

        glDepthMask( GL_FALSE );
        glDisable(GL_LIGHTING);
        glEnable(GL_LINE_SMOOTH);
        glLineWidth( 2.0 );
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);;
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

        glPushMatrix();
        glColor4f( 0.5, 0.5, 0.2, 0.5 );
        glTranslatef(ObjInfo->MagEphemInfo->P.x, ObjInfo->MagEphemInfo->P.y, ObjInfo->MagEphemInfo->P.z );  Lgm_gl_draw_sphere( TRUE, 0.11, 20, 20 );
        glPopMatrix();

        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glDisable( GL_BLEND );
        glDisable(GL_LINE_SMOOTH);
        glEnable(GL_LIGHTING);
        glDepthMask( GL_TRUE );
    glEndList( );

}

void LoadStars( ) {

    // Stars
    FILE        *fp;
    float       val, sRed, sGrn, sBlu;
    double      RA, DEC, MAG;
    Lgm_Vector  P_gei, P;
    char        TYPE;
    GLfloat     quadratic[3] = {1.0f, 0.0f, 0.0f };

    if ( ( fp = fopen("/home/mgh/DREAM/Dream/Dream/Misc/Stars.txt", "r") ) != NULL ) {

        StarsDL = glGenLists( 1 );
        glNewList( StarsDL, GL_COMPILE );
            glDisable(GL_LIGHTING);
            glPointSize( 3.0 );
            glPointParameterfv( GL_POINT_DISTANCE_ATTENUATION, &quadratic[0] );
            glEnable(GL_POINT_SMOOTH);
            glEnable(GL_BLEND);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);;
            glBegin( GL_POINTS );
              while( fscanf(fp, "%lf %lf %lf %c", &RA, &DEC, &MAG, &TYPE) != EOF ){
                if ( MAG < MaxStarMagnitude ) {
                    Lgm_Radec_to_Cart( 15.0*RA, DEC, &P_gei );
                    Lgm_ScaleVector( &P_gei, 70.0 );
                    Lgm_Convert_Coords( &P_gei, &P, StarsConvertFlag, mInfo->c );
                    StarColor( TYPE, &sRed, &sGrn, &sBlu ); // sets the star color based on its spectral type
                    val = (0.1-1.5)/(6.5 - -1.47)*(MAG - -1.47) + 1.5; // controls brightness based on stars magnitude
                    glColor3f( sRed*val, sGrn*val, sBlu*val);
                    glVertex3f( P.x, P.y, P.z );
                }
              }
              fclose(fp);
            glEnd();
            glDisable(GL_BLEND);
            glEnable(GL_LIGHTING);
        glEndList( );



        SunDL = glGenLists( 1 );
        glNewList( SunDL, GL_COMPILE );

            glEnable(GL_POINT_SPRITE);
            glEnable(GL_TEXTURE_2D);
            glDepthMask( GL_FALSE );
            glDisable(GL_LIGHTING);
            glEnable(GL_POINT_SMOOTH);
            glEnable(GL_BLEND);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);;
//        glPointParameterfv( GL_POINT_DISTANCE_ATTENUATION, quadratic );
            glTexEnvf(GL_POINT_SPRITE, GL_COORD_REPLACE, GL_TRUE);

            glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE );
            glBindTexture( GL_TEXTURE_2D, Texture_Sun );

            glPointSize( 32.0 );
            glBegin( GL_POINTS );
                    Lgm_Radec_to_Cart( mInfo->c->RA_sun, mInfo->c->DEC_sun, &P_gei );
                    Lgm_ScaleVector( &P_gei, 69.0 );
                    Lgm_Convert_Coords( &P_gei, &P, StarsConvertFlag, mInfo->c );
                    glColor3f( 1.0, 1.0, 1.0);
                    glVertex3f( P.x, P.y, P.z );
            glEnd();
            glDisable(GL_BLEND);
            glDisable(GL_POINT_SMOOTH);
            glEnable(GL_LIGHTING);
            glDepthMask( GL_TRUE );
            glDisable(GL_TEXTURE_2D);
            glDisable(GL_POINT_SPRITE);
        glEndList( );
    }


}

void ReLoadStars( ) {
    glDeleteLists( StarsDL, 1 );
    glDeleteLists( SunDL, 1 );
    LoadStars( );
}


/*
 * get axis/angle rot required for billboarding
 *
 *      P is the position of the billboard
 *      camera is position of camera
 *
 * Use like this;
 *
 *          glRotatef( Angle, Axis.x, Axis.y, Axis.z );
 *          glTranslatef( p.x, p.y, p.z );
 */
void GetBillboardAxisAngle( Lgm_Vector *p, Lgm_Vector *camera, Lgm_Vector *Axis, double *Angle ){


    double      Q[4], A[3][3];
    Lgm_Vector  X, Y, Z;


    // get u-hat vector -- it points to the camera
    X.x = camera->x - p->x;
    X.y = camera->y - p->y;
    X.z = camera->z - p->z;
    Lgm_NormalizeVector( &X );

    // get w-hat vector (its perp to u x up)
    Lgm_CrossProduct( &X, &aInfo->Up, &Z );
    Lgm_NormalizeVector( &Z );

    // get vhat to complete
    Lgm_CrossProduct( &Z, &X, &Y );

    // now form the trans matrix from original coord system to the local billboard system
    A[0][0] = X.x, A[1][0] = Y.x, A[2][0] = Z.x;
    A[0][1] = X.y, A[1][1] = Y.y, A[2][1] = Z.y;
    A[0][2] = X.z, A[1][2] = Y.z, A[2][2] = Z.z;

    // convert trans matrix to a quat
    Lgm_MatrixToQuat( A, Q );

    // convert quat to axis/angle
    Lgm_QuatToAxisAngle(  Q, Angle, Axis );



}





int LoadTLEs( ){

    int         i;
    char        Filename[1024];
    double      r, g, b, mag;
    _SgpTLE     *TLEs;
    int         nTLEs=0;

    /*
     * Read in TLEs
     */
    TLEs = (_SgpTLE *)calloc( 20000, sizeof(_SgpTLE) );
    for (i=0; i<nSatTypes; i++){
        if ( ShowSatellites[i] ) {
            nTLEs = 0;
            sprintf(Filename, "/home/mgh/DREAM/Dream/Data/TLEs/%s_Latest.txt", SatTypes[i]);
            LgmSgp_ReadTlesFromFile( Filename, &nTLEs, TLEs, 0 );
            printf("File = %s nTLEs found=%d\n", Filename, nTLEs);
            SpaceObjects->nSat = nTLEs;
        }
    }

    SpaceObjects->Sat = (_SpaceObjectItem *)calloc( nTLEs, sizeof(_SpaceObjectItem) );
    SpaceObjects->nSat = nTLEs;

    for (i=0; i<SpaceObjects->nSat; i++){

        // copy over the TLE
        SpaceObjects->Sat[i].TLE = TLEs[i];

        SpaceObjects->Sat[i].DrawStreak              = FALSE;
        SpaceObjects->Sat[i].sPeriodFrac             = 25.0;
        SpaceObjects->Sat[i].DrawGroundPathOfStreak  = FALSE;
        SpaceObjects->Sat[i].DrawStreakToGroundLines = FALSE;

        SpaceObjects->Sat[i].DrawOrbit               = TRUE;
        SpaceObjects->Sat[i].oPeriodFrac             = 100.0;
        SpaceObjects->Sat[i].DrawGroundPathOfOrbit   = FALSE;
        SpaceObjects->Sat[i].DrawOrbitToGroundLines  = FALSE;

        SpaceObjects->Sat[i].DrawSatToGroundLine     = FALSE;

        SpaceObjects->Sat[i].DrawLabel               = TRUE;

        // Streak Colors
        r = (double)rand()/(double)RAND_MAX;
        g = (double)rand()/(double)RAND_MAX;
        b = (double)rand()/(double)RAND_MAX;
        mag = sqrt(r*r+g*g+b*b); r /= mag; g /= mag; b /= mag;
        SpaceObjects->Sat[i].sRed = r;
        SpaceObjects->Sat[i].sGrn = g;
        SpaceObjects->Sat[i].sBlu = b;
        SpaceObjects->Sat[i].sAlf = 0.75;

        SpaceObjects->Sat[i].sgpRed = 0.75*r;
        SpaceObjects->Sat[i].sgpGrn = 0.75*g;
        SpaceObjects->Sat[i].sgpBlu = 0.75*b;
        SpaceObjects->Sat[i].sgpAlf = 0.75*0.75;

        SpaceObjects->Sat[i].sglRed = 0.75*r;
        SpaceObjects->Sat[i].sglGrn = 0.75*g;
        SpaceObjects->Sat[i].sglBlu = 0.75*b;
        SpaceObjects->Sat[i].sglAlf = 0.75*0.75;

        // Orbit Colors
        SpaceObjects->Sat[i].oRed = 0.6;
        SpaceObjects->Sat[i].oGrn = 0.6;
        SpaceObjects->Sat[i].oBlu = 0.6;
        SpaceObjects->Sat[i].oAlf = 0.7;

        SpaceObjects->Sat[i].ogpRed = 1.0;
        SpaceObjects->Sat[i].ogpGrn = 0.0;
        SpaceObjects->Sat[i].ogpBlu = 0.0;
        SpaceObjects->Sat[i].ogpAlf = 1.0;

        SpaceObjects->Sat[i].oglRed = 0.7;
        SpaceObjects->Sat[i].oglGrn = 0.7;
        SpaceObjects->Sat[i].oglBlu = 0.7;
        SpaceObjects->Sat[i].oglAlf = 0.1;

        // Single Line
        SpaceObjects->Sat[i].ssglRed = r;
        SpaceObjects->Sat[i].ssglGrn = g;
        SpaceObjects->Sat[i].ssglBlu = b;
        SpaceObjects->Sat[i].ssglAlf = 0.5*0.75;

    }

    free( TLEs );
    return(1);

}


void ReLoadTLEs( ) {
    if ( SpaceObjects->Sat ) free( SpaceObjects->Sat );
    LoadTLEs( );
}

void DrawSatLabels(){

    int             ii;
    cairo_surface_t *cst;
    cairo_t         *cr;
    GLuint          temp_tex;
    GLdouble        ModelViewMatrix[16];
    GLdouble        ProjMatrix[16];
    GLint           ViewPort[4];
    GLdouble        winx, winy, winz, x1, x2, y1, y2;
    GLfloat         bufferZ;
    _GroupNode      *g;
    _SpaceObjects   *Group;

    g = SatSelectorInfo->SatGroupList;

    if ( g == NULL ) return;

    glGetDoublev( GL_MODELVIEW_MATRIX, ModelViewMatrix );
    glGetDoublev( GL_PROJECTION_MATRIX, ProjMatrix );
    glGetIntegerv( GL_VIEWPORT, ViewPort );
    glGenTextures( 1, &temp_tex );


    while ( g != NULL ) {
        Group = g->Group;

        if ( Group->DrawGroup ) {

            for (ii=0; ii<Group->nSat; ii++){

                if ( Group->Sat[ii].Draw && Group->Sat[ii].DrawLabel ){

                    // get window coords of object
                    gluProject( Group->Sat[ii].x, Group->Sat[ii].y, Group->Sat[ii].z, ModelViewMatrix, ProjMatrix, ViewPort, &winx, &winy, &winz);

                    if (winz >= 0.0){

                        // get latest z-buffer value to see if wee are ocluded...
                        glReadPixels(winx, winy,1,1,GL_DEPTH_COMPONENT, GL_FLOAT, &bufferZ);

                        cst = cairo_image_surface_create( CAIRO_FORMAT_ARGB32, 250, 28 );
                        cr  = cairo_create( cst );
                        cairo_select_font_face( cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL );

                        cairo_set_source_rgba( cr, 0.0, 0.0, 0.0, 0.0 );
                        cairo_rectangle( cr, 0.0, 0.0, 250.0, 28.0 );
                        cairo_fill(cr);

                        if (bufferZ >= winz){
                            cairo_set_source_rgba( cr, 1.0, 1.0, 1.0, 1.0 );
//cairo_set_source_rgba( cr, 0.0, 0.0, 0.0, 1.0 );
                        } else {
                            // it is occluded, draw dimly (or not at all?)
                            //cairo_set_source_rgba( cr, 1.0, 1.0, 1.0, 0.3 );
                        }
                        cairo_set_font_size( cr, 24.0 );
                        cairo_move_to( cr, 0.0, 24.0 ); cairo_show_text( cr, Group->Sat[ii].TLE.Name );

                        glBindTexture( GL_TEXTURE_2D, temp_tex );
                        glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA, 250, 28, 0, GL_BGRA, GL_UNSIGNED_BYTE, (const GLvoid *)cairo_image_surface_get_data( cst ) );
                        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
                        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
                        glEnable( GL_TEXTURE_2D );
                        glDisable( GL_LIGHTING );
                        glEnable(GL_BLEND);
                        glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );

                        winx = 2.0*winx/(double)g_imageWidth - 1.0;
                        winy = 2.0*winy/(double)g_imageHeight - 1.0;
                        x1 = winx; x2 = winx+500.0/(double)g_imageWidth;
                        y1 = winy; y2 = winy+40.0/(double)g_imageHeight;

                        glMatrixMode(GL_MODELVIEW); glPushMatrix(); glLoadIdentity();
                        glMatrixMode( GL_PROJECTION ); glPushMatrix(); glLoadIdentity();
                        gluOrtho2D(-1.0, -1.0, 1.0, 1.0);
                        glEnable( GL_TEXTURE_2D );
                        glDisable( GL_LIGHTING );
                        glEnable(GL_BLEND);
                        glDisable(GL_DEPTH_TEST);
                        glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );

                        glColor4f( 1.0, 1.0, 1.0, 0.5 );
                        glBegin(GL_QUADS);
                            glTexCoord2f( 0.0, 1.0 ); glVertex2f( x1, y1 );
                            glTexCoord2f( 1.0, 1.0 ); glVertex2f( x2, y1 );
                            glTexCoord2f( 1.0, 0.0 ); glVertex2f( x2, y2 );
                            glTexCoord2f( 0.0, 0.0 ); glVertex2f( x1, y2 );
                        glEnd();

                        glDisable( GL_TEXTURE_2D );
                        glEnable( GL_LIGHTING );
                        glDisable(GL_BLEND);
                        glEnable(GL_DEPTH_TEST);
                        glMatrixMode( GL_PROJECTION ); glPopMatrix();
                        glMatrixMode( GL_MODELVIEW ); glPopMatrix();

                        glEnable( GL_LIGHTING );
                        glDisable(GL_BLEND);
                        glDisable( GL_TEXTURE_2D );


                        cairo_destroy( cr );
                        cairo_surface_destroy( cst );

                    }

                }

            }

        }
        g = g->Next;
    }

    glDeleteTextures( 1, &temp_tex );

    return;

}

void CreateSats() {

    double           tsince, JD, tUT;
    long int         tDate;
    Lgm_Vector       Ugsm, Ugei, EarthToSun, EarthToSun_obs, g1, g2, g3, P, Pobs, Uobs;
    int              i, Flag1=0, Flag2=0, Flag3=0, tYear, tMonth, tDay;
    _GroupNode      *g;
    _SpaceObjects   *Group;
    _SgpInfo        *s = (_SgpInfo *)calloc( 1, sizeof(_SgpInfo) );;
    Lgm_CTrans      *c = Lgm_init_ctrans( 0 );

    /*
     * All the TLEs have their own epoch times0 in them. And the propagator (sgp4)
     * uses the "time since (in minutes)". So for a given time of interest, we need to
     * compute the tsince needed.
     */
    JD = CurrentJD;

    Lgm_jd_to_ymdh( JD, &tDate, &tYear, &tMonth, &tDay, &tUT );
    Lgm_Set_Coord_Transforms( tDate, tUT, c );

    SatsDL = glGenLists( 1 );
    glNewList( SatsDL, GL_COMPILE );
//        glPointParameterf( GL_POINT_FADE_THRESHOLD_SIZE, 256.0 );
//        GLfloat sizes[2];
//        glGetFloatv( GL_ALIASED_POINT_SIZE_RANGE, sizes );
//        printf("sizes = %f %f\n", sizes[0], sizes[1] );
//        glPointParameterf( GL_POINT_SIZE_MIN_, sizes[0] );
//        glPointParameterf( GL_POINT_SIZE_MIN, 256.0 );

        g = SatSelectorInfo->SatGroupList;
        while ( g != NULL ) {
            Group = g->Group;
            if ( Group->DrawGroup && Group->DrawPosition) {


                if ( Group->DrawAsPointSprites ) {
                    glEnable(GL_POINT_SPRITE);
                    glEnable(GL_TEXTURE_2D);
                }

                glDepthMask( GL_FALSE );
                glDisable(GL_LIGHTING);
                glEnable(GL_POINT_SMOOTH);
                glEnable(GL_BLEND);
                glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);;




                /*
                 * Debris
                 */
                glColor4f( Group->DebRed, Group->DebGrn, Group->DebBlu, Group->DebAlf );
                if ( Group->DrawAsPointSprites ) {
                    float quadratic[] = {1.0, 0.0, 0.1 };
                    glPointParameterfv( GL_POINT_DISTANCE_ATTENUATION, quadratic );
                    glTexEnvf(GL_POINT_SPRITE, GL_COORD_REPLACE, GL_TRUE);
                    glPointSize( 32.0 );
                    glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE );
                    glBindTexture( GL_TEXTURE_2D, Texture_Debris );
                } else {
                    glPointSize( 10.0 );
                }


                glBegin( GL_POINTS );
                for (i=0; i<Group->nSat; i++){
                    if ( strstr(Group->Sat[i].TLE.Name, " DEB") && Group->DrawDebris ){
                        tsince = (JD - Group->Sat[i].TLE.JD)*1440.0;
                        LgmSgp_SGP4_Init( s, &Group->Sat[i].TLE );
                        LgmSgp_SGP4( tsince, s );
                        Ugei.x = s->X/Re; Ugei.y = s->Y/Re; Ugei.z = s->Z/Re;
                        Lgm_Convert_Coords( &Ugei, &Ugsm, SatsConvertFlag, c );
                        //glColor4f( 0.7, 0.1, 0.0, 0.5 );
                        glVertex3f( Ugsm.x, Ugsm.y, Ugsm.z );
                        Group->Sat[i].x = Ugsm.x; Group->Sat[i].y = Ugsm.y; Group->Sat[i].z = Ugsm.z;
                    }
                }
                glEnd();




                /*
                 * RocketBodies
                 */
                glColor4f( Group->RbRed, Group->RbGrn, Group->RbBlu, Group->RbAlf );
                if ( Group->DrawAsPointSprites ) {
                    float quadratic[] = {1.0, 0.0, 0.1 };
                    glPointParameterfv( GL_POINT_DISTANCE_ATTENUATION, quadratic );
                    glTexEnvf(GL_POINT_SPRITE, GL_COORD_REPLACE, GL_TRUE);
                    glPointSize( 32.0 );
                    glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE );
                    glBindTexture( GL_TEXTURE_2D, Texture_RocketBody );
                } else {
                    glPointSize( 10.0 );
                }

                glBegin( GL_POINTS );
                    for (i=0; i<Group->nSat; i++){
                        if ( strstr(Group->Sat[i].TLE.Name, " R/B") && Group->DrawRocketBodies ){ // RocketBodies
                            tsince = (JD - Group->Sat[i].TLE.JD)*1440.0;
                            LgmSgp_SGP4_Init( s, &Group->Sat[i].TLE );
                            LgmSgp_SGP4( tsince, s );
                            Ugei.x = s->X/Re; Ugei.y = s->Y/Re; Ugei.z = s->Z/Re;
                            Lgm_Convert_Coords( &Ugei, &Ugsm, SatsConvertFlag, c );
                            //glColor4f( 0.5, 0.4, 0.0, 0.5 );
                            glVertex3f( Ugsm.x, Ugsm.y, Ugsm.z );
                            Group->Sat[i].x = Ugsm.x; Group->Sat[i].y = Ugsm.y; Group->Sat[i].z = Ugsm.z;
                        }
                    }
                glEnd();



                /*
                 * Satellites
                 */
                glColor4f( Group->SatRed, Group->SatGrn, Group->SatBlu, Group->SatAlf );
                if ( Group->DrawAsPointSprites ) {
                    float quadratic[] = {1.0, 0.0, 0.1 };
                    glPointParameterfv( GL_POINT_DISTANCE_ATTENUATION, quadratic );
                    glTexEnvf(GL_POINT_SPRITE, GL_COORD_REPLACE, GL_TRUE);
                    glPointSize( 64.0 );
                    glPointSize( 128.0 );
                    glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE );
                    glBindTexture( GL_TEXTURE_2D, Texture_Spacecraft );
                } else {
                    glPointSize( 10.0 );
                }

                glBegin( GL_POINTS );
                    for (i=0; i<Group->nSat; i++){
                        if ( (!strstr(Group->Sat[i].TLE.Name, " R/B") && !strstr(Group->Sat[i].TLE.Name, " DEB")) && Group->DrawSatellites ){ // Satellites
                            tsince = (JD - Group->Sat[i].TLE.JD)*1440.0;
                            LgmSgp_SGP4_Init( s, &Group->Sat[i].TLE );
                            LgmSgp_SGP4( tsince, s );
                            Ugei.x = s->X/Re; Ugei.y = s->Y/Re; Ugei.z = s->Z/Re;
                            Lgm_Convert_Coords( &Ugei, &Ugsm, SatsConvertFlag, c );
//printf("Date, Time: %ld %lf   Sat pnt: Ugei.x = %g Ugei.y = %g Ugei.z = %g\n", Ugei.x, Ugei.y, Ugei.z );
//printf("Sat pnt: Ugsm.x = %g Ugsm.y = %g Ugsm.z = %g\n", Ugsm.x, Ugsm.y, Ugsm.z );
                            //glColor4f( 1.0, 1.0, 1.0, 0.5 );
                            glVertex3f( Ugsm.x, Ugsm.y, Ugsm.z );
                            Group->Sat[i].x = Ugsm.x; Group->Sat[i].y = Ugsm.y; Group->Sat[i].z = Ugsm.z;
                        }
                    }
                glEnd();

/*
for (i=0; i<Group->nSat; i++){
if ( (!strstr(Group->Sat[i].TLE.Name, " R/B") && !strstr(Group->Sat[i].TLE.Name, " DEB")) && Group->DrawSatellites ){ // Satellites
glEnable(GL_LIGHTING);
//glEnable(GL_BLEND);
//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);;
printf("Group->Sat[%d] = %g %g %g\n", i, Group->Sat[i].x, Group->Sat[i].y, Group->Sat[i].z );
Ugsm.x = Group->Sat[i].x; Ugsm.y = Group->Sat[i].y; Ugsm.z = Group->Sat[i].z;
SolidCone( &Ugsm, 30.0, 240, 1 );
//glDisable(GL_BLEND);
glDisable(GL_LIGHTING);
}
}
*/

            }

            g = g->Next;
        }

        glDisable(GL_POINT_SPRITE);
        glDisable(GL_TEXTURE_2D);

        // draw a line from object to ground
        glEnable(GL_LINE_SMOOTH);
        glLineWidth( 3.0 );
        g = SatSelectorInfo->SatGroupList;
        while ( g != NULL ) {
            Group = g->Group;
            if ( Group->DrawGroup ) {
                for (i=0; i<Group->nSat; i++){
                    if (Group->Sat[i].DrawSatToGroundLine) {
                        tsince = (JD - Group->Sat[i].TLE.JD)*1440.0;
                        LgmSgp_SGP4_Init( s, &Group->Sat[i].TLE );
                        glBegin( GL_LINES );
                            LgmSgp_SGP4( tsince, s );
                            Ugei.x = s->X/Re; Ugei.y = s->Y/Re; Ugei.z = s->Z/Re;
                            Lgm_Convert_Coords( &Ugei, &Ugsm, SatsConvertFlag, c );
                            glColor4f( Group->Sat[i].ssglRed, Group->Sat[i].ssglGrn, Group->Sat[i].ssglBlu, Group->Sat[i].ssglAlf );
                            glVertex3f( Ugsm.x, Ugsm.y, Ugsm.z );
                            Lgm_NormalizeVector( &Ugsm );
                            Lgm_ScaleVector( &Ugsm, 1.00 );
                            glVertex3f( Ugsm.x, Ugsm.y, Ugsm.z );
                        glEnd();
                    }
                }
            }
            g = g->Next;
        }


        // draw lines from sat to ground due to iridium flares
        glEnable(GL_LINE_SMOOTH);
        glLineWidth( 3.0 );
double th,ph,r;
Lgm_Vector u, uu;
        th =  35.888;
        ph = -106.306;
        r  = 0.000;
        Lgm_GEOD_to_WGS84( th, ph, r, &u );
        //printf("u= %g %g %g\n", u.x, u.y, u.z );
        //Radius = Lgm_Magnitude( &u );

u.x = r*sin(th)*cos(ph); u.y = r*sin(th)*sin(ph); u.z = r*cos(th);
Lgm_Convert_Coords( &u, &uu, GEO_TO_GEI2000, c );

        // convert sun from gei->obs
        EarthToSun = c->Sun; Lgm_ScaleVector( &EarthToSun, c->earth_sun_dist );
        Lgm_Convert_Coords( &EarthToSun, &EarthToSun_obs, SatsConvertFlag, c );

        glBegin( GL_LINES );
        g = SatSelectorInfo->SatGroupList;
        while ( g != NULL ) {
            Group = g->Group;

            for (i=0; i<Group->nSat; i++){

                if ( strstr( Group->Sat[i].TLE.Name, "IRIDIUM" ) ) {
                    if ( ShowIridiumSatToGround || ShowIridiumSatToSun ) {
                        IridiumFlare( JD, &Group->Sat[i].TLE, &EarthToSun, &Flag1, &Flag2, &Flag3, &g1, &g2, &g3, &P );
                        Lgm_Convert_Coords( &P, &Pobs, SatsConvertFlag, c ); // convert S/C pos from GEI->Obs

                        if ( ShowIridiumSatToGround ) {
                            if (Flag1){
                                glColor4f( 1.0, 0.0, 0.0, 1.0);
                                glVertex3f( Pobs.x, Pobs.y, Pobs.z );
                                Lgm_Convert_Coords( &g1, &Uobs, SatsConvertFlag, c );
                                printf("g1 diff: %g\n", Re*Lgm_VecDiffMag( &uu, &g1 ));
                                glVertex3f( Uobs.x, Uobs.y, Uobs.z );
                            }

                            if (Flag2){
                                glColor4f( 0.0, 1.0, 0.0, 1.0);
                                glVertex3f( Pobs.x, Pobs.y, Pobs.z );
                                Lgm_Convert_Coords( &g2, &Uobs, SatsConvertFlag, c );
                                glVertex3f( Uobs.x, Uobs.y, Uobs.z );
                                printf("g2 diff: %g\n", Re*Lgm_VecDiffMag( &uu, &g2 ));
                            }

                            if (Flag3){
                                glColor4f( 0.0, 0.0, 1.0, 1.0);
                                glVertex3f( Pobs.x, Pobs.y, Pobs.z );
                                Lgm_Convert_Coords( &g3, &Uobs, SatsConvertFlag, c );
                                glVertex3f( Uobs.x, Uobs.y, Uobs.z );
                                printf("g3 diff: %g\n", Re*Lgm_VecDiffMag( &uu, &g3 ));
                            }
                        }
                    }

                    // draw sun to sat line?
                    if ( ShowIridiumSatToSun && (Flag1 || Flag2 || Flag3) ) {
                        glColor4f( 1.0, 1.0, 0.0, 0.1);
                        glVertex3f( Pobs.x, Pobs.y, Pobs.z );
                        glVertex3f( EarthToSun_obs.x, EarthToSun_obs.y, EarthToSun_obs.z );
                    }
                }
            }
            g = g->Next;
        }
        glEnd();




        glDisable(GL_BLEND);
        glEnable(GL_LIGHTING);
        glDepthMask( GL_TRUE );

    glEndList( );

    free(s);
    Lgm_free_ctrans( c ); // free the CTrans structure

}

void CreateThemisFovs() {

    int     i, j, k;

    /*
     *  This section of code adds the Themis FOVs as overlays.
            20 Goose Bay GBAY 53.316 N 299.540 E 61.233 N 23.114E 3:36 03 GBO-14 GMAG-6 (10532/5009) Feb-06
            18 Kuujjuaq KUUJ 58.155 N 291.468 E 67.364 N 13.143E 4:12 13 GBO-13 GMAG-8 (10547/5011) Nov-07
            19 Chibougam au CHBG 49.814 N 285.581 E 60.048 N 3.249E 4:44 16 GBO-17 GMAG-9 (10546/5013) Sep-06
            16 Sanikiluaq SNKQ 56.536 N 280.769 E 66.895 N 3.466W 5:08 09 GBO-22 NRCan w/ GPS-9 Oct-06
            17 Kapuskasin g KAPU 49.392 N 277.680 E 60.183 N 8.605W 5:29 21 GBO-15 GMAG-7 (10545/5012) May- 06
            10 Rankin Inlet RANK 62.828 N 267.887 E 72.748 N 25.138W 6:25 12 GBO-09 CGSM w/ GPS-4 (10528) Sep-05
            13 Gillam GILL 56.354 N 265.344 E 66.502 N 28.043W 6:34 19 GBO-19 CGSM w/ GPS-7 (10516) May- 06
            15 Pinawa (LdB) PINA 50.163 N 263.934 E 60.375 N 29.294W 6:39 18 GBO-16 CGSM w/ GPS-8 May- 06
            14 The Pas TPAS 53.994 N 259.059 E 63.532 N 37.044W 7:05 08 GBO-06 GMAG-1 (10505/4001) May- 05
            11 Fort Smith FSMI 59.984 N 248.158 E 67.551 N 54.383W 8:08 10 GBO-10 CGSM w/ GPS-3 (10527) Jul-05
            12 Athabasca ATHA 54.714 N 246.686 E 62.128 N 54.096W 8:08 02 GBO-02 NRCan w/ GPS-0 Aug-04
            7 / Ekati EKAT 64.717 N 250.667 E 72.474 N 53.555W 8:11 04 GBO-04 GMAG-3 (10503/4003) Dec-04
            9 Prince George PGEO 53.815 N 237.172 E 59.218 N 65.096W 8:53 15 GBO-03 GMAG-2 (10501/4002) Sep-04
            8 Fort Simpson FSIM 61.762 N 238.779 E 67.407 N 67.180W 8:58 05 GBO-21 CGSM w/ GPS-6 (10539) Nov-06
            6 White Horse WHIT 61.010 N 224.777 E 63.690 N 81.709W 10:02 07 GBO-07 GMAG-4 (10533/4015) Jul-05
            5 Inuvik INUV 68.413 N 226.230 E 71.246 N 86.064W 10:19 17 GBO-08 GMAG-11 (10550/5017) Jun-05
            1 Gakona GAKO 62.407 N 214.842 E 63.051 N 91.750W 10:49 20 GBO-18 GI w/ GPS-10 Aug-06
            2 Fort Yukon FYKN 66.560 N 214.786 E 67.218 N 94.757W 11:02 14 GBO-12 GI w/ GPS-5 (10529) Oct-05
            3 Mcgrath MCGR 62.953 N 204.404 E 61.679 N 100.823 W 11:33 11 GBO-11 GMAG-5 (10525/4016) Aug-05
            4 Kiana KIAN 66.971 199.562 65.071 107.234 12:04 22 GBO-20 GMAG-10 Sep-06

            Rabbit Lake 53.1416 252.2326
            Lucky Lake  50.9833 252.8667
    */

    Lgm_Vector v, Z, X, Y, Pole;
    double  GeodLat, GeodLon;
    double  Ageo_to_asi[3][3];
    double  Aasi_to_geo[3][3];

    int     N_THEMIS_ASI = 23;

    double  THEMIS_ASI_LAT[] = { 53.316, 58.155, 49.814, 56.536, 49.392, 
                                 62.828, 56.354, 50.163, 53.994, 59.984, 
                                 54.714, 64.717, 53.815, 61.762, 61.010, 
                                 68.413, 62.407, 66.560, 62.953, 66.971,  58.22, 50.9833, 48.4284 };

    double  THEMIS_ASI_LON[] = { 299.540, 291.468, 285.581, 280.769, 277.680, 
                                 267.887, 265.344, 263.934, 259.059, 248.158, 
                                 246.686, 250.667, 237.172, 238.779, 224.777, 
                                 226.230, 214.842, 214.786, 204.404, 199.562, 256.33, 252.8667, 236.6344 };

    int     TREX_STATION[]   = { 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0 };

    
    double      tUT;
    long int    tDate;
    int         tYear, tMonth, tDay;
    Lgm_CTrans  *c = Lgm_init_ctrans( 0 );
    
    ThemisFovsDL = glGenLists( 1 );
    glNewList( ThemisFovsDL, GL_COMPILE );
    glDepthMask( GL_FALSE );
    glDisable(GL_LIGHTING);
    glEnable(GL_LINE_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);;

    Lgm_jd_to_ymdh( CurrentJD, &tDate, &tYear, &tMonth, &tDay, &tUT );
    Lgm_Set_Coord_Transforms( tDate, tUT, c );

    /*
     *  Geo position of the ASI
     */
    for ( k=0; k<N_THEMIS_ASI; k++ ){
        GeodLat = THEMIS_ASI_LAT[k];
        GeodLon = THEMIS_ASI_LON[k];
        //GeodLat =  61.762;
        //GeodLon = 238.779;
        Lgm_GEOD_to_WGS84( GeodLat, GeodLon, 0.0, &v );
        Lgm_ForceMagnitude( &v, WGS84_A ); //km Lgm_GEOD_to_WGS84() appears to use WGS84_A to reprort units of Re.


        /*
         *  Create a local coord system.
         *    Z - zenith pointing
         *    Y - perp to both Z and rotation axis in east direction
         *    X - completes (points northward)
         */
        Z.x = v.x;
        Z.y = v.y;
        Z.z = v.z;
        Lgm_NormalizeVector( &Z );


        Pole.x = 0.0; Pole.y = 0.0; Pole.z = 1.0;
        Lgm_CrossProduct( &Z, &Pole, &Y );
        Lgm_NormalizeVector( &Y );
        
        Lgm_CrossProduct( &Y, &Z, &X );
        Lgm_NormalizeVector( &X );


        /*
         *  Setup transformation matrix between local system and GEO
         */
        Ageo_to_asi[0][0] = X.x; Ageo_to_asi[1][0] = X.y; Ageo_to_asi[2][0] = X.z;
        Ageo_to_asi[0][1] = Y.x; Ageo_to_asi[1][1] = Y.y; Ageo_to_asi[2][1] = Y.z;
        Ageo_to_asi[0][2] = Z.x; Ageo_to_asi[1][2] = Z.y; Ageo_to_asi[2][2] = Z.z;

        for ( i=0; i<3; i++ ) {
            for ( j=0; j<3; j++ ) {
                Aasi_to_geo[i][j] = Ageo_to_asi[j][i];
            }
        }

        /*
         *   Now create a circular FOV in asi coords and intersect with 110km altitude sphere.
         */
        RayType Ray;
        EllipsoidType Sphere_110;
        Lgm_Vector u_asi, u_geo, u_obs;
        double phi, theta, tmin, tmax, t;
        Sphere_110.Origin.x  = Sphere_110.Origin.y = Sphere_110.Origin.z = 0.0;
        Sphere_110.Radius_a  = WGS84_A+110.0; // km 
        Sphere_110.Radius_b  = WGS84_B+110.0; // km 
        Sphere_110.Radius2_a = Sphere_110.Radius_a*Sphere_110.Radius_a;
        Sphere_110.Radius2_b = Sphere_110.Radius_b*Sphere_110.Radius_b;
        theta = 82.5*RadPerDeg; // 85 deg. half angle (total 170 deg FOV)

        Ray.Origin.x = v.x; // km
        Ray.Origin.y = v.y; // km
        Ray.Origin.z = v.z; // km

        glLineWidth( 1.0 );
        if ( TREX_STATION[k]  ){
            glColor4f( 1.0, 0.2, 0.2, 0.7);
        } else {
            glColor4f( 0.5, 0.5, 0.5, 0.7);
        }
        glBegin( GL_LINE_STRIP );
        for ( phi=0.0; phi <=2.0*M_PI+0.09; phi += 0.1 ){

            u_asi.x = cos( phi ) * sin( theta );
            u_asi.y = sin( phi ) * sin( theta );
            u_asi.z = cos( theta );

            Lgm_MatTimesVec( Aasi_to_geo, &u_asi, &u_geo );

            Ray.Direction.x = u_geo.x;
            Ray.Direction.y = u_geo.y;
            Ray.Direction.z = u_geo.z;

                
            Lgm_EllipsoidIntersect( &Sphere_110, &Ray, &tmin, &tmax, &t );
            //printf("tmin, tmax, t = %g %g %g\n", tmin, tmax, t);


            // intersection point 
            u_geo.x = Ray.Origin.x + tmax*Ray.Direction.x;
            u_geo.y = Ray.Origin.y + tmax*Ray.Direction.y;
            u_geo.z = Ray.Origin.z + tmax*Ray.Direction.z;


            Lgm_Convert_Coords( &u_geo, &u_obs, GeoToObsConvertFlag, c );

            glVertex3f( u_obs.x/Re, u_obs.y/Re, u_obs.z/Re );

        }
        glEnd();

    }
    glDisable(GL_BLEND);
    glEnable(GL_LIGHTING);
    glDepthMask( GL_TRUE );

    glEndList( );

    Lgm_free_ctrans( c ); // free the CTrans structure 

}

void ReCreateThemisFovs( ) {
    glDeleteLists( ThemisFovsDL, 1 );
    CreateThemisFovs( );
}



void ReCreateSats( ) {
    glDeleteLists( SatsDL, 1 );
    CreateSats( );
}

void CreateSatOrbits() {

    double              tsince, tsince0, JD, dt, Height, FLL;
    double              tJD[10001], tUT[10001], period, tinc, OrbitFrac = 0.125, Theta;
    long int            tDate[10001];
    int                 tYear, tMonth, tDay, Flag;
    Lgm_Vector          Ugei[10001], UGSM[10001], Ugsm[10001], aa, bb, uu, v1[10001], v2[10001], v3[10001];
    Lgm_Vector          Wgsm, Wcoord;
    int                 i, j, jj, n, nMax;
    Lgm_CTrans          *c = Lgm_init_ctrans( 0 );
    _GroupNode          *g;
    _SpaceObjects       *Group;
    _SgpInfo            *s = (_SgpInfo *)calloc( 1, sizeof(_SgpInfo) );
    Lgm_MagModelInfo    *mInfo = Lgm_InitMagInfo();

    /*
     * All the TLEs have their own epoch times in them. And the propagator (sgp4)
     * uses the "time since (in minutes)". So for a given time of interest, we need to
     * compute the tsince needed.
     */
    JD = CurrentJD;

    Lgm_MagModelInfo_Set_MagModel( LGM_CDIP, LGM_EXTMODEL_T89, mInfo );

    SatOrbitsDL = glGenLists( 1 );
    glNewList( SatOrbitsDL, GL_COMPILE );
    glDepthMask( GL_FALSE );
    glDisable(GL_LIGHTING);
    glEnable(GL_LINE_SMOOTH);
    glLineWidth( 5.0 );
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);;

    g = SatSelectorInfo->SatGroupList;
    while ( g != NULL ) {
        Group = g->Group;

        if ( Group->DrawGroup ) {

            /*
             *  DRAW ORBIT, ETC
             *  Orbits are drawn from a half orbit before current time to
             *  half orbit after.  (If oPeriodFrac is 100%. lengths of
             *  orbits are scaled appropriately for other values of
             *  oPeriodFrac)
             */
            for (i=0; i<Group->nSat; i++){
                if ( Group->Sat[i].Draw ) {

                    if ( (Group->Sat[i].DrawOrbit) || (Group->Sat[i].DrawGroundPathOfOrbit) || (Group->Sat[i].DrawOrbitToGroundLines) ){
                        LgmSgp_SGP4_Init( s, &Group->Sat[i].TLE );
                        tsince0 = (JD - Group->Sat[i].TLE.JD)*1440.0;
                        period = 1440.0/Group->Sat[i].TLE.MeanMotion; // orbit period in minutes
period *= Group->Sat[i].oPeriodFrac/100.0;

                        // init n to zero and start half an orbit ahead of current time.
                        n=0; dt = period/2.0;

                        // do first one to start
                        tinc = period/1000.0; // start out with very small increment
                        tsince = tsince0+dt;
                        LgmSgp_SGP4( tsince, s );
                        Ugei[n].x = s->X/Re; Ugei[n].y = s->Y/Re; Ugei[n].z = s->Z/Re;
                        tJD[n] = JD + dt/1440.0;
                        Lgm_jd_to_ymdh( tJD[n], &tDate[n], &tYear, &tMonth, &tDay, &tUT[n] );
                        Lgm_Set_Coord_Transforms( tDate[n], tUT[n], c );
                        Lgm_Convert_Coords( &Ugei[n], &Ugsm[n], SatsConvertFlag, c ); // Unfortunately, Ugsm looks like it isnt always GSM!
                        Lgm_Convert_Coords( &Ugei[n], &UGSM[n], GEI2000_TO_GSM, c );  // Save a copy of GMS that is really GSM...
                        aa = Ugsm[n]; Lgm_NormalizeVector( &aa );
                        n++; dt -= tinc;

                        while ( (dt >= -period/2.0) && (n<10000) ) {
                            tsince = tsince0+dt;
                            LgmSgp_SGP4( tsince, s );
                            Ugei[n].x = s->X/Re; Ugei[n].y = s->Y/Re; Ugei[n].z = s->Z/Re;
                            tJD[n] = JD + dt/1440.0;
                            Lgm_jd_to_ymdh( tJD[n], &tDate[n], &tYear, &tMonth, &tDay, &tUT[n] );
                            Lgm_Set_Coord_Transforms( tDate[n], tUT[n], c );
                            Lgm_Convert_Coords( &Ugei[n], &Ugsm[n], SatsConvertFlag, c );
                            Lgm_Convert_Coords( &Ugei[n], &UGSM[n], GEI2000_TO_GSM, c );


                            bb = Ugsm[n]; Lgm_NormalizeVector( &bb );
//Lgm_VecSub( &cc, &bb, &aa );
                            Theta = DegPerRad * acos( fabs(Lgm_DotProduct( &aa, &bb )) );
                            if (Theta > 1.0) tinc *= 0.61803398875;
                            if (Theta < 0.5) tinc *= 1.61803398875;
//                            DS = Lgm_Magnitude( &cc );
//                            if (DS > 1.0)  tinc *= 0.61803398875;
//                            if (DS < 0.5) tinc *= 1.61803398875;


                            n++; dt -= tinc;
                            aa = bb;
                        }
//printf("n = %d\n", n);


                        // ORBIT
                        if ( Group->Sat[i].DrawOrbit ) {
glLineWidth( 5.0 );
                            glColor4f( Group->Sat[i].oRed, Group->Sat[i].oGrn, Group->Sat[i].oBlu, Group->Sat[i].oAlf );
glColor4f( 0.6, 0.6, 0.6, 0.7);
//glColor4f( 0.0, 0.0, 0.0, 0.8 );
                            glBegin( GL_LINE_STRIP );
                                for (j=0; j<n; j++) glVertex3f( Ugsm[j].x, Ugsm[j].y, Ugsm[j].z );
                            glEnd();
                        }
                        // GROUND PATH
                        if ( Group->Sat[i].DrawGroundPathOfOrbit ) {
                            glColor4f( Group->Sat[i].ogpRed, Group->Sat[i].ogpGrn, Group->Sat[i].ogpBlu, Group->Sat[i].ogpAlf );
glColor4f( 0.6, 0.6, 0.6, 0.7);
                            glBegin( GL_LINE_STRIP );
                            for (j=0; j<n; j++) {
                                uu = Ugsm[j];
                                Lgm_ForceMagnitude( &uu, 1.01 );
                                glVertex3f( uu.x, uu.y, uu.z );
                            }
                            glEnd();
                        }
                        // LINES BETWEEN ORBIT AND GROUND PATH
                        if ( Group->Sat[i].DrawOrbitToGroundLines ) {
                            nMax = 25;
                            //glColor4f( Group->Sat[i].oglRed, Group->Sat[i].oglGrn, Group->Sat[i].oglBlu, Group->Sat[i].oglAlf*(double)n/(double)nMax );
                            glColor4f( Group->Sat[i].oglRed, Group->Sat[i].oglGrn, Group->Sat[i].oglBlu, Group->Sat[i].oglAlf );
glLineWidth( 1.0 );
glColor4f( 0.0, 1.0, 1.0, 0.3);
                            glBegin( GL_LINES );
                                for (j=0; j<n; j += 1) {
                                    uu = Ugsm[j];
                                    Lgm_ForceMagnitude( &uu, 1.01 );
                                    glVertex3f( Ugsm[j].x, Ugsm[j].y, Ugsm[j].z );
                                    glVertex3f( uu.x, uu.y, uu.z );
                                }
                            glEnd();
                        }

                        // ORBIT FIELD LINES
                        if ( Group->Sat[i].DrawOrbitFieldLines || Group->Sat[i].DrawOrbitFLFootpoints ) {

                            nMax = 25;
                            //glColor4f( Group->Sat[i].oglRed, Group->Sat[i].oglGrn, Group->Sat[i].oglBlu, Group->Sat[i].oglAlf*(double)n/(double)nMax );
glLineWidth( 2.0 );
                            glColor4f( Group->Sat[i].oglRed, Group->Sat[i].oglGrn, Group->Sat[i].oglBlu, Group->Sat[i].oglAlf );

                            for (j=0; j<n; j += 1) {
                                Height = 120.0;
                                Lgm_Set_Coord_Transforms( tDate[j], tUT[j], mInfo->c );
                                
                                Flag = Lgm_Trace( &UGSM[j], &v1[j], &v2[j], &v3[j], Height, 1e-7, 1e-7, mInfo );
                                FLL = mInfo->Stotal;
                                printf("FLL = %g\n", FLL);

                                if ( FLL > 0.0 ){
                                    Lgm_TraceLine3( &v1[j], FLL, 100, 1.0, 1e-7, 0, mInfo );

                                    if ( Group->Sat[i].DrawOrbitFieldLines ) {
                                        glBegin( GL_LINE_STRIP );
                                            for (jj=0; jj<mInfo->nPnts; jj++) {
                                                Wgsm.x = mInfo->Px[jj]; Wgsm.y = mInfo->Py[jj]; Wgsm.z = mInfo->Pz[jj];
                                                Lgm_Convert_Coords( &Wgsm, &Wcoord, AtmosConvertFlag, mInfo->c );
                                                glVertex3f( Wcoord.x, Wcoord.y, Wcoord.z );
                                            }
                                        glEnd();
                                    }
                                }
                            }

                            if ( Group->Sat[i].DrawOrbitFLFootpoints ) {
glLineWidth( 4.0 );
glColor4f( 0.6, 0.6, 0.6, 0.7);
                                glBegin( GL_LINE_STRIP );
                                    for (j=0; j<n; j++) {
                                        uu = v2[j];
                                        //Lgm_ForceMagnitude( &uu, 1.001 );
                                        //glColor4f( Group->Sat[i].ogpRed, Group->Sat[i].ogpGrn, Group->Sat[i].ogpBlu, Group->Sat[i].ogpAlf*(1.0-(double)j/(double)nMax) );
                                        //glColor4f( Group->Sat[i].ogpRed, Group->Sat[i].ogpGrn, Group->Sat[i].ogpBlu, 1.0 );
                                        Lgm_Set_Coord_Transforms( tDate[j], tUT[j], mInfo->c );
                                        Lgm_Convert_Coords( &uu, &Wcoord, AtmosConvertFlag, mInfo->c );
                                        glVertex3f( Wcoord.x, Wcoord.y, Wcoord.z );
                                    }
                                glEnd();
                                glBegin( GL_LINE_STRIP );
                                    for (j=0; j<n; j++) {
                                        uu = v1[j];
                                        //Lgm_ForceMagnitude( &uu, 1.001 );
                                        //glColor4f( Group->Sat[i].ogpRed, Group->Sat[i].ogpGrn, Group->Sat[i].ogpBlu, 1.0 );
                                        Lgm_Set_Coord_Transforms( tDate[j], tUT[j], mInfo->c );
                                        Lgm_Convert_Coords( &uu, &Wcoord, AtmosConvertFlag, mInfo->c );
                                        glVertex3f( Wcoord.x, Wcoord.y, Wcoord.z );
                                    }
                                glEnd();
                            }
                            
                        }

                    }
                }
            }


            /*
             *  DRAW STREAKS, ETC
             */
            Group = g->Group;
            for (i=0; i<Group->nSat; i++){
                if ( Group->Sat[i].Draw ) {

                    if ( (Group->Sat[i].DrawStreak) || (Group->Sat[i].DrawGroundPathOfStreak) || (Group->Sat[i].DrawStreakToGroundLines) ){
                        LgmSgp_SGP4_Init( s, &Group->Sat[i].TLE );
                        tsince0 = (JD - Group->Sat[i].TLE.JD)*1440.0;
                        period = 1440.0/Group->Sat[i].TLE.MeanMotion; // orbit period in minutes

period *= Group->Sat[i].sPeriodFrac/100.0;

                        // init n to zero and start at current time
                        n=0; dt = 0.0;

                        // do first one to start
                        tinc = period/1000.0; // start out with very small increment
                        tsince = tsince0+dt;
                        LgmSgp_SGP4( tsince, s );
                        Ugei[n].x = s->X/Re; Ugei[n].y = s->Y/Re; Ugei[n].z = s->Z/Re;
                        tJD[n] = JD + dt/1440.0;
                        Lgm_jd_to_ymdh( tJD[n], &tDate[n], &tYear, &tMonth, &tDay, &tUT[n] );
                        Lgm_Set_Coord_Transforms( tDate[n], tUT[n], c );
                        Lgm_Convert_Coords( &Ugei[n], &Ugsm[n], SatsConvertFlag, c );
                        Lgm_Convert_Coords( &Ugei[n], &UGSM[n], GEI2000_TO_GSM, c );
                        aa = Ugsm[n]; Lgm_NormalizeVector( &aa );
                        n++; dt -= tinc;

                        // go until we've gone backwards OrbitFrac of an orbit
                        while ( (dt >= -period/2.0) && (n<10000) ) {
                            tsince = tsince0+dt;
                            LgmSgp_SGP4( tsince, s );
                            Ugei[n].x = s->X/Re; Ugei[n].y = s->Y/Re; Ugei[n].z = s->Z/Re;
                            tJD[n] = JD + dt/1440.0;
                            Lgm_jd_to_ymdh( tJD[n], &tDate[n], &tYear, &tMonth, &tDay, &tUT[n] );
                            Lgm_Set_Coord_Transforms( tDate[n], tUT[n], c );
                            Lgm_Convert_Coords( &Ugei[n], &Ugsm[n], SatsConvertFlag, c );
                            Lgm_Convert_Coords( &Ugei[n], &UGSM[n], GEI2000_TO_GSM, c );


                            bb = Ugsm[n]; Lgm_NormalizeVector( &bb );
                            Theta = DegPerRad * acos( fabs(Lgm_DotProduct( &aa, &bb )) );
                            if (Theta > 1.0) tinc *= 0.61803398875;
                            if (Theta < 0.5) tinc *= 1.61803398875;


                            n++; dt -= tinc;
                            aa = bb;
                        }

                        nMax = n;
                        //printf("nMax = %d\n", nMax);


                        // STREAK
                        if ( Group->Sat[i].DrawStreak ) {
                            glBegin( GL_LINE_STRIP );
                                for (j=0; j<n; j++) {
                                    glColor4f( Group->Sat[i].sRed, Group->Sat[i].sGrn, Group->Sat[i].sBlu, Group->Sat[i].sAlf*(1.0-(double)j/(double)nMax) );
                                    glVertex3f( Ugsm[j].x, Ugsm[j].y, Ugsm[j].z );
                                }
                            glEnd();
                        }
                        // STREAK GROUND PATH
                        if ( Group->Sat[i].DrawGroundPathOfStreak ) {
                            glBegin( GL_LINE_STRIP );
                                for (j=0; j<n; j++) {
                                    uu = Ugsm[j];
                                    Lgm_ForceMagnitude( &uu, 1.001 );
                                    glColor4f( Group->Sat[i].sgpRed, Group->Sat[i].sgpGrn, Group->Sat[i].sgpBlu, Group->Sat[i].sgpAlf*(1.0-(double)j/(double)nMax) );
                                    glVertex3f( uu.x, uu.y, uu.z );
                                }
                            glEnd();
                        }
                        // LINES BETWEEN STREAK AND GROUND PATH
                        if ( Group->Sat[i].DrawStreakToGroundLines ) {
                            glBegin( GL_LINES );
                                for (j=0; j<n; j++) {
                                    uu = Ugsm[j];
                                    Lgm_ForceMagnitude( &uu, 1.001 );
                                    glColor4f( Group->Sat[i].sglRed, Group->Sat[i].sglGrn, Group->Sat[i].sglBlu, Group->Sat[i].sglAlf*(1.0-(double)j/(double)nMax) );
                                    glVertex3f( Ugsm[j].x, Ugsm[j].y, Ugsm[j].z );
                                    glVertex3f( uu.x, uu.y, uu.z );
                                }
                            glEnd();
                        }

                        // STREAK FIELD LINES
                        if ( Group->Sat[i].DrawStreakFieldLines || Group->Sat[i].DrawStreakFLFootpoints ) {

                            //nMax = 25;
                            //glColor4f( Group->Sat[i].oglRed, Group->Sat[i].oglGrn, Group->Sat[i].oglBlu, Group->Sat[i].oglAlf*(double)n/(double)nMax );
                            //glColor4f( Group->Sat[i].sglRed, Group->Sat[i].sglGrn, Group->Sat[i].sglBlu, Group->Sat[i].sglAlf );
                            glColor4f( Group->Sat[i].sglRed, Group->Sat[i].sglGrn, Group->Sat[i].sglBlu, Group->Sat[i].sglAlf );

                            for (j=0; j<n; j += 1) {
                                Height = 120.0;
                                Lgm_Set_Coord_Transforms( tDate[j], tUT[j], mInfo->c );
                                Flag = Lgm_Trace( &UGSM[j], &v1[j], &v2[j], &v3[j], Height, 1e-7, 1e-7, mInfo );
                                FLL = mInfo->Stotal;
                                printf("FLL = %g\n", FLL);

                                if ( FLL > 0.0 ){
                                    Lgm_TraceLine3( &v1[j], FLL, 100, 1.0, 1e-7, 0, mInfo );

                                    if ( Group->Sat[i].DrawStreakFieldLines ) {
                                        glLineWidth( 2.0 );
                                        glColor4f( Group->Sat[i].sglRed, Group->Sat[i].sglGrn, Group->Sat[i].sglBlu, Group->Sat[i].sglAlf*(1.0-(double)j/(double)nMax) );
                                        glBegin( GL_LINE_STRIP );
                                            for (jj=0; jj<mInfo->nPnts; jj++) {

                                                //glVertex3f( mInfo->Px[jj], mInfo->Py[jj], mInfo->Pz[jj] );
                                                Wgsm.x = mInfo->Px[jj]; Wgsm.y = mInfo->Py[jj]; Wgsm.z = mInfo->Pz[jj];
                                                Lgm_Convert_Coords( &Wgsm, &Wcoord, AtmosConvertFlag, mInfo->c );
                                                glVertex3f( Wcoord.x, Wcoord.y, Wcoord.z );
                                            }
                                        glEnd();
                                    }
                                }
                            }

                            if ( Group->Sat[i].DrawStreakFLFootpoints ) {
                                glLineWidth( 4.0 );
                                glBegin( GL_LINE_STRIP );
                                    for (j=0; j<n; j++) {
                                        uu = v2[j];
                                        //Lgm_ForceMagnitude( &uu, 1.001 );
                                        glColor4f( Group->Sat[i].sgpRed, Group->Sat[i].sgpGrn, Group->Sat[i].sgpBlu, Group->Sat[i].sgpAlf*(1.0-(double)j/(double)nMax) );
                                        Lgm_Set_Coord_Transforms( tDate[j], tUT[j], mInfo->c );
                                        Lgm_Convert_Coords( &uu, &Wcoord, AtmosConvertFlag, mInfo->c );
                                        glVertex3f( Wcoord.x, Wcoord.y, Wcoord.z );
                                    }
                                glEnd();
                                glLineWidth( 4.0 );
                                glBegin( GL_LINE_STRIP );
                                    for (j=0; j<n; j++) {
                                        uu = v1[j];
                                        //Lgm_ForceMagnitude( &uu, 1.001 );
                                        glColor4f( Group->Sat[i].sgpRed, Group->Sat[i].sgpGrn, Group->Sat[i].sgpBlu, Group->Sat[i].sgpAlf*(1.0-(double)j/(double)nMax) );
                                        Lgm_Set_Coord_Transforms( tDate[j], tUT[j], mInfo->c );
                                        Lgm_Convert_Coords( &uu, &Wcoord, AtmosConvertFlag, mInfo->c );
                                        glVertex3f( Wcoord.x, Wcoord.y, Wcoord.z );
                                    }
                                glEnd();
                            }
                            
                        }

                    }
                }
            }

        }


        g = g->Next;
    }



    glDisable(GL_BLEND);
    glEnable(GL_LIGHTING);
    glDepthMask( GL_TRUE );

    glEndList( );

    free(s);
    Lgm_free_ctrans( c ); // free the CTrans structure 
    Lgm_FreeMagInfo( mInfo );

}

void ReCreateSatOrbits( ) {
    glDeleteLists( SatOrbitsDL, 1 );
    CreateSatOrbits( );
}






static void init_view (void) {
    float view_quat1[4];
    float view_quat2[4];
    double Alpha;

    view_scale = 1.2;


    /*
     *  Nominally,
     *      z-axis points out of screen
     *      x-axis is to the right
     *      y-axis is up.
     *
     *  Define two rotations to force;
     *      z-axis up
     *      x-axis to left
     *      y-axis out of screen
     */

    // rotate by 180deg. around z-axis
    Alpha = M_PI;
    view_quat1[0] = sin( 0.5 * Alpha )*0.0;
    view_quat1[1] = sin( 0.5 * Alpha )*0.0;
    view_quat1[2] = sin( 0.5 * Alpha )*1.0;
    view_quat1[3] = cos( 0.5 * Alpha );

    // rotate by 270deg. around x-axis
    Alpha = 270.0*M_PI/180.0;
    view_quat2[0] = sin( 0.5 * Alpha )*1.0;
    view_quat2[1] = sin( 0.5 * Alpha )*0.0;
    view_quat2[2] = sin( 0.5 * Alpha )*0.0;
    view_quat2[3] = cos( 0.5 * Alpha );

    // merge the rots into a single quaternion
    add_quats( view_quat1, view_quat2, view_quat );

}



/*
   Create a sphere centered at c, with radius r, and precision n
   Draw a point for zero radius spheres
*/
void CreateSphere( double r, int n) {

    int          i, j;
    double       theta1, theta2, Phi;
    Lgm_Vector   e, p;

    if (r < 0) r = -r;
    if (n < 0) n = -n;
    if ( (n < 4) || (r <= 0) ) {
        glBegin( GL_POINTS );
        glVertex3f( 0.0, 0.0, 0.0 );
        glEnd( );
        return;
    }

    for (j=0;j<n/2;j++) {
        // range in latitude from -PI/2 -> PI/2
        theta1 = j * 2.0*M_PI / n - 0.5*M_PI;
        theta2 = (j + 1) * 2.0*M_PI / n - 0.5*M_PI;

        glBegin(GL_QUAD_STRIP);
        for (i=0;i<=n;i++) {
            Phi = i * 2.0*M_PI / n;

            e.x = cos(theta2) * cos(Phi);
            e.y = cos(theta2) * sin(Phi);
            e.z = sin(theta2);
            p.x = r * e.x;
            p.y = r * e.y;
            p.z = r * e.z;

            glNormal3f(e.x,e.y,e.z);
            glTexCoord2f(i/(double)n,2*(j+1)/(double)n);
            glVertex3f(p.x,p.y,p.z);

            e.x = cos(theta1) * cos(Phi);
            e.y = cos(theta1) * sin(Phi);
            e.z = sin(theta1);
            p.x = r * e.x;
            p.y = r * e.y;
            p.z = r * e.z;

            glNormal3f(e.x,e.y,e.z);
            glTexCoord2f(i/(double)n,2*j/(double)n);
            glVertex3f(p.x,p.y,p.z);
        }
        glEnd();
    }
}


/*
   Create an ellipsoid centered at c, with equatorial radius ra, polar radius rb, and precision n
   Draw a point for zero radius spheres
*/
void CreateEllipsoid( double ra, double rb, int n) {

    int          i, j;
    double       ct1, ct2, st1, st2, cp, sp;
    double       r, theta1, theta2, Phi, g, h;
    Lgm_Vector   e, p;


    if (ra < 0.0) ra = -ra;
    if (rb < 0.0) rb = -rb;
    if (n < 0) n = -n;
    if ( (n < 4) || (ra <= 0.0) || (rb <= 0.0) ) {
        glBegin( GL_POINTS );
        glVertex3f( 0.0, 0.0, 0.0 );
        glEnd( );
        return;
    }

    g = ra*ra/(rb*rb);
    h = 1.0-g;

    for (j=0;j<n/2;j++) {
        // range in latitude from -PI/2 -> PI/2
        theta1 = j * 2.0*M_PI / n - 0.5*M_PI;
        theta2 = (j + 1) * 2.0*M_PI / n - 0.5*M_PI;
        ct1 = cos(theta1); st1 = sin(theta1);
        ct2 = cos(theta2); st2 = sin(theta2);

        glBegin(GL_QUAD_STRIP);
        for (i=0;i<=n;i++) {
            Phi = i * 2.0*M_PI / n;
            cp = cos(Phi); sp = sin(Phi);

            e.x = ct2*cp;
            e.y = ct2*sp;
            e.z = st2;
            r = ra/sqrt( ct2*ct2*h + g );
            p.x = r * e.x;
            p.y = r * e.y;
            p.z = r * e.z;

            glNormal3f(e.x,e.y,g*e.z);
            glTexCoord2f(i/(double)n,2*(j+1)/(double)n);
            glVertex3f(p.x,p.y,p.z);

            e.x = ct1*cp;
            e.y = ct1*sp;
            e.z = st1;
            r = ra/sqrt( ct1*ct1*h + g );
            p.x = r * e.x;
            p.y = r * e.y;
            p.z = r * e.z;

            glNormal3f(e.x,e.y,g*e.z);
            glTexCoord2f(i/(double)n,2*j/(double)n);
            glVertex3f(p.x,p.y,p.z);
        }
        glEnd();
    }
}

/*
   Create an ellipsoid centered at c, with equatorial radius ra, polar radius rb, and precision n
   Draw a point for zero radius spheres
*/
void CreateCutEllipsoid( double ra, double rb, int n) {

    int          i, j;
    double       ct1, ct2, st1, st2, cp, sp, Phi0, Phi1, dPhi, fPhi;
    double       r, theta1, theta2, Phi, g, h;
    Lgm_Vector   e, p;

    Phi0 =  90.0*RadPerDeg;
    Phi1 = 360.0*RadPerDeg;
    dPhi = Phi1-Phi0;


    if (ra < 0.0) ra = -ra;
    if (rb < 0.0) rb = -rb;
    if (n < 0) n = -n;
    if ( (n < 4) || (ra <= 0.0) || (rb <= 0.0) ) {
        glBegin( GL_POINTS );
        glVertex3f( 0.0, 0.0, 0.0 );
        glEnd( );
        return;
    }

    g = ra*ra/(rb*rb);
    h = 1.0-g;

    for (j=0;j<n/2;j++) {
        // range in latitude from -PI/2 -> PI/2
        theta1 = j * 2.0*M_PI / n - 0.5*M_PI;
        theta2 = (j + 1) * 2.0*M_PI / n - 0.5*M_PI;
        ct1 = cos(theta1); st1 = sin(theta1);
        ct2 = cos(theta2); st2 = sin(theta2);

        glBegin(GL_QUAD_STRIP);
        for (i=0;i<=n;i++) {
            // range in longitude from Phi0 to Phi1
            Phi = dPhi*i/(double)n + Phi0;
            fPhi = Phi/(2.0*M_PI);
            cp = cos(Phi); sp = sin(Phi);

            e.x = ct2*cp;
            e.y = ct2*sp;
            e.z = st2;
            r = ra/sqrt( ct2*ct2*h + g );
            p.x = r * e.x;
            p.y = r * e.y;
            p.z = r * e.z;

            glNormal3f(e.x,e.y,g*e.z);
            glTexCoord2f( fPhi, 2*(j+1)/(double)n );
            glVertex3f(p.x,p.y,p.z);

            e.x = ct1*cp;
            e.y = ct1*sp;
            e.z = st1;
            r = ra/sqrt( ct1*ct1*h + g );
            p.x = r * e.x;
            p.y = r * e.y;
            p.z = r * e.z;

            glNormal3f(e.x,e.y,g*e.z);
            glTexCoord2f( fPhi, 2*j/(double)n );
            glVertex3f(p.x,p.y,p.z);
        }
        glEnd();
    }
}






static void realize( GtkWidget *widget, gpointer data) {



    GdkGLContext    *glcontext  = gtk_widget_get_gl_context( widget );
    GdkGLDrawable   *gldrawable = gtk_widget_get_gl_drawable( widget );






    /*
     * Load the textures in
     */
    LoadTextures();

    GLfloat ambient[]         = {0.0, 0.0, 0.0, 1.0};
    GLfloat diffuse[]         = {1.0, 1.0, 1.0, 1.0};
    //GLfloat position[]        = {0.0, 3.0, 3.0, 0.0};
    GLfloat lmodel_ambient[]  = {0.2, 0.2, 0.2, 1.0};
    GLfloat local_view[]      = {0.0};

    /*
     * OpenGL BEGIN
     */
    if (!gdk_gl_drawable_gl_begin (gldrawable, glcontext)) return;

    glClearColor( 0.0, 0.0, 0.0, 0.0 );
//glClearColor( 0.8, 0.8, 0.8, 0.8 );




    glClearDepth( 1.0 );

    glLightfv( GL_LIGHT0, GL_AMBIENT, ambient);
    glLightfv( GL_LIGHT0, GL_DIFFUSE, diffuse);
    //glLightfv( GL_LIGHT0, GL_POSITION, position);
    glLightfv( GL_LIGHT0, GL_POSITION, LightPosition);
    glLightModelfv( GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
    glLightModelfv( GL_LIGHT_MODEL_LOCAL_VIEWER, local_view);

    glFrontFace( GL_CW );
    glEnable( GL_LIGHTING );
    glEnable( GL_LIGHT0 );
    glEnable( GL_AUTO_NORMAL );
    glEnable( GL_NORMALIZE );
    glEnable( GL_DEPTH_TEST );
    glEnable( GL_LINE_SMOOTH );
    glHint( GL_LINE_SMOOTH_HINT, GL_NICEST );
    glDepthFunc( GL_LESS );
    glShadeModel( GL_SMOOTH );


// these are time dep.
    // Create 3D GSM axes
    CreateGSMAxes( );
//    CreateZPSAxes( );


    // Create Earth + dipole field, etc..
    CreateEarth( );

    // Create Earth + dipole field, etc..
    CreateThemisFovs( );

    // Create Earth + dipole field, etc..
    CreateMoon( );

    // Create a DL for stars.
    LoadStars( );

    // Load sats
//    LoadTLEs( );
    CreateSats();
    CreateSatOrbits();
    CreateLogo();




// these are SC specific
// Need to process new mag data files
// for this to update properly...
//    TraceFieldLines( );

    // Create sphere to represent Spacecraft Position
    CreateSCPos( );

    // Create FL tubes
    GenerateFieldLineLists( ObjInfo );

    // Create Drift Shell surfaces
    GenerateDriftShellLists( ObjInfo );


    // Misc Field Lines
//GenerateMiscFieldLineLists( ObjInfo );



    // Create Full screen quad (only for use with depth peeling)
    //MakeFullScreenQuad();

    init_view( );
    gdk_gl_drawable_gl_end( gldrawable );
    /*
     *  OpenGL END
     */

  return;
}




static gboolean configure_event( GtkWidget *widget, GdkEventConfigure *event, gpointer data) {

    GdkGLContext  *glcontext = gtk_widget_get_gl_context (widget);
    GdkGLDrawable *gldrawable = gtk_widget_get_gl_drawable (widget);

    double w = widget->allocation.width;
    double h = widget->allocation.height;
    double Aspect;
    double FieldOfView = 45.0;


   /*
    * OpenGL BEGIN
    */
    if ( !gdk_gl_drawable_gl_begin( gldrawable, glcontext ) ) return FALSE;


    glewInit();

    // initialize depth peeling stuff.
    g_imageWidth  = w;
    g_imageHeight = h;
    g_imageRat = (double)g_imageWidth/(double)g_imageHeight;
    BuildShaders2();
//    BuildShaders();
//    InitFrontPeelingRenderTargets();
//    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);

    // set the projection matrix -- use a Perspective projection
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity( );
    Aspect = ( w > h ) ? w/h : h/w;
    gluPerspective( FieldOfView, Aspect, 0.01, 100.0 );

    // set viewport
    glViewport (0, 0, w, h);
    //glDepthRange(0.0, 0.5);

    // reset HUD texture
    glDeleteTextures( 1, &Texture_TEMP );
    glGenTextures( 1, &Texture_TEMP );

    // reset Logo texture
    ReCreateLogo();

    // reset HUD texture
    glDeleteTextures( 1, &Texture_HiResEarthQuad );
    glGenTextures( 1, &Texture_HiResEarthQuad );

    gdk_gl_drawable_gl_end( gldrawable );
   /*
    * OpenGL END
    */

    return TRUE;

}

void DrawScene( ) {

    char        Str[256];
    int         i, ty, tm, td, tD;
    int         iii, ii, jj, nn, nph, nth;
    Lgm_Vector  Camera_geo;
    double      th, ph;
    int         qj, qi, DestroyQuad, LoadQuad;
    long int    QuadsNeeded[9];
    GLuint      QuadsNeededDL[9];
    GLuint      QuadsNeededTexId[9];

//gluLookAt( aInfo->Camera.x, aInfo->Camera.y, aInfo->Camera.z, 10.0, 0.0, 0.0, 0.0, 0.0, 1.0 );
//GLdouble eyeX, GLdouble eyeY, GLdouble eyeZ, GLdouble centerX, GLdouble centerY, GLdouble centerZ, GLdouble upX, GLdouble upY, GLdouble upZ )


    /* Render Objects */
//    glUseProgram( g_shaderMyTest );




    /*
     * These are the various coordinate axes.
     */
    glCallList( AxesDL );



//    glUseProgram( 0 );

    if ( ShowEarth ) glCallList( EarthDL );

    if ( ShowThemisFovs ) glCallList( ThemisFovsDL );


//    CreateGeoMarkers();
//    glCallList( GeoMarkersDL );




//20100305
//glUseProgram( g_shaderMyTest );
glCallList( DipoleAxisDL );
//glUseProgram( 0 );
//20100305
glCallList( SunDirectionDL );
//20100305    glCallList( EqPlaneGridDL );
//glCallList( EqPlaneDL );


/*
Lgm_Vector Ugeo, Ugsm, RotAxis9;
double Glon, Agsm_to_ins[3][3], Ains_to_gsm[3][3], RotAngle9;
double Qins_to_gsm[4], Qgsm_to_ins[4];
Glon = 0.0*M_PI/180.0;
Ugeo.x=6.6*cos(Glon); Ugeo.y=6.6*sin(Glon); Ugeo.z = 0.0;
Lgm_Convert_Coords( &Ugeo, &Ugsm, GEO_TO_GSM, mInfo->c );
ComputeZPSTransMatrix( &Ugeo, CurrentDate, CurrentUT, Agsm_to_ins, Ains_to_gsm, Qins_to_gsm, Qgsm_to_ins );
Lgm_QuatToAxisAngle( Qins_to_gsm, &RotAngle9, &RotAxis9 );

glPushMatrix();
glRotatef( RotAngle9, RotAxis9.x, RotAxis9.y, RotAxis9.z );
glRotatef( RotAngle3, RotAxis3.x, RotAxis3.y, RotAxis3.z ); // This implements coord trans from Observer coords -> GSM
//glTranslatef( Ugsm.x, Ugsm.y, Ugsm.z );
glCallList( ZPSAxesDL );
glPopMatrix();
*/


    /*
     * Show stars
CHECK COORDS!
     */
    if ( ShowStars ){
        glPushMatrix();
        glTranslatef( aInfo->Camera.x, aInfo->Camera.y, aInfo->Camera.z );
        glCallList( StarsDL );
        glCallList( SunDL );
        glPopMatrix();
    }


    glCallList( MoonDL );


    glPushMatrix();
    glRotatef( RotAngle3, RotAxis3.x, RotAxis3.y, RotAxis3.z ); // This implements coord trans from Observer coords -> GSM






if (LightingStyle == 2){
//  cgGLBindProgram(myCgVertexProgram);
//  checkForCgError("binding vertex program");
//  cgGLEnableProfile(myCgVertexProfile);
//  checkForCgError("enabling vertex profile");
//  cgGLBindProgram(myCgFragmentProgram);
//  checkForCgError("binding fragment program");
//  cgGLEnableProfile(myCgFragmentProfile);
//  checkForCgError("enabling fragment profile");
}
/*
*/

//    glUseProgram( g_shaderMyTest );
//    GLint   SurfaceColorLoc = glGetUniformLocation( g_shaderMyTest, "SurfaceColor" );
//    GLint   PLoc            = glGetUniformLocation( g_shaderMyTest, "P" );
//    GLint   ALoc            = glGetUniformLocation( g_shaderMyTest, "A" );
//    GLint   ScaleLoc        = glGetUniformLocation( g_shaderMyTest, "Scale" );
//    GLint   LightDirLoc     = glGetUniformLocation( g_shaderMyTest, "LightDir" );
//    GLint   ViewPositionLoc = glGetUniformLocation( g_shaderMyTest, "ViewPosition" );
//    GLfloat SurfaceColor[4], ViewPosition[3];
//
//ViewPosition[0] = aInfo->Camera.x;
//ViewPosition[1] = aInfo->Camera.y;
//ViewPosition[2] = aInfo->Camera.z;
//glUniform3fv( ViewPositionLoc, 1, LightPosition );
//glUniform3fv( LightDirLoc, 1, LightPosition );





glMaterialfv( GL_FRONT, GL_AMBIENT,   mat_blue_plastic.ambient);
glMaterialfv( GL_FRONT, GL_DIFFUSE,   mat_blue_plastic.diffuse);
glMaterialfv( GL_FRONT, GL_SPECULAR,  mat_blue_plastic.specular);
glMaterialf(  GL_FRONT, GL_SHININESS, mat_blue_plastic.shininess * 128.0);
glCallList( ObjInfo->MiscFieldLines );




    if ( ShowAllPitchAngles ) {
        /*
         * We want to show all of the PAs
         */

ShowFullFieldLine = 1;
        if ( ShowFullFieldLine ){
            for (i=0; i<ObjInfo->MagEphemInfo->nAlpha; i++ ) {
                glMaterialfv( GL_FRONT, GL_DIFFUSE,   gInfo->FieldLineMaterial[i].diffuse );
                glMaterialfv( GL_FRONT, GL_AMBIENT,   gInfo->FieldLineMaterial[i].ambient );
                glMaterialfv( GL_FRONT, GL_SPECULAR,  gInfo->FieldLineMaterial[i].specular );
                glMaterialf(  GL_FRONT, GL_SHININESS, gInfo->FieldLineMaterial[i].shininess*128.0 );
//SurfaceColor[0] = gInfo->FieldLineMaterial[i].diffuse[0];
//SurfaceColor[1] = gInfo->FieldLineMaterial[i].diffuse[1];
//SurfaceColor[2] = gInfo->FieldLineMaterial[i].diffuse[2];
//SurfaceColor[3] = gInfo->FieldLineMaterial[i].diffuse[3];
//glUniform4fv( SurfaceColorLoc, 1, SurfaceColor );
                glCallList( ObjInfo->DriftShellList2 + i );
            }
        } else {
            for (i=0; i<ObjInfo->MagEphemInfo->nAlpha; i++ ) {
                glMaterialfv( GL_FRONT, GL_DIFFUSE,   gInfo->FieldLineMaterial[i].diffuse );
                glMaterialfv( GL_FRONT, GL_AMBIENT,   gInfo->FieldLineMaterial[i].ambient );
                glMaterialfv( GL_FRONT, GL_SPECULAR,  gInfo->FieldLineMaterial[i].specular );
                glMaterialf(  GL_FRONT, GL_SHININESS, gInfo->FieldLineMaterial[i].shininess*128.0 );
//SurfaceColor[0] = gInfo->FieldLineMaterial[i].diffuse[0];
//SurfaceColor[1] = gInfo->FieldLineMaterial[i].diffuse[1];
//SurfaceColor[2] = gInfo->FieldLineMaterial[i].diffuse[2];
//SurfaceColor[3] = gInfo->FieldLineMaterial[i].diffuse[3];
//glUniform4fv( SurfaceColorLoc, 1, SurfaceColor );
                glCallList( ObjInfo->DriftShellList3 + i );
            }
        }

    } else {

        /*
         * We want to show a subset of Pitch Angles
         */

        for (i=0; i<ObjInfo->MagEphemInfo->nAlpha; i++ ) {
            if ( ShowPitchAngle[i] ) {
                glMaterialfv( GL_FRONT, GL_DIFFUSE,   gInfo->FieldLineMaterial[i].diffuse );
                glMaterialfv( GL_FRONT, GL_AMBIENT,   gInfo->FieldLineMaterial[i].ambient );
                glMaterialfv( GL_FRONT, GL_SPECULAR,  gInfo->FieldLineMaterial[i].specular );
                glMaterialf(  GL_FRONT, GL_SHININESS, gInfo->FieldLineMaterial[i].shininess*128.0 );
//SurfaceColor[0] = gInfo->FieldLineMaterial[i].diffuse[0];
//SurfaceColor[1] = gInfo->FieldLineMaterial[i].diffuse[1];
//SurfaceColor[2] = gInfo->FieldLineMaterial[i].diffuse[2];
//SurfaceColor[3] = gInfo->FieldLineMaterial[i].diffuse[3];
//glUniform4fv( SurfaceColorLoc, 1, SurfaceColor );
                if ( ShowFullFieldLine ){
                    glCallList( ObjInfo->DriftShellList2 + i );
                } else {
                    glCallList( ObjInfo->DriftShellList3 + i);
                }
            }
        }

    }

//    glUseProgram(0);




    /*
     * Convert Camera Position (in Observer Coords) to GEO
     */
    Lgm_Convert_Coords( &aInfo->Camera, &Camera_geo, ObsToGeoConvertFlag, mInfo->c );
//    ra = Lgm_Magnitude( &Camera_geo );                  // radius
    th = DegPerRad*acos( Camera_geo.z/aInfo->Rcam );             // co-latitude
    ph = DegPerRad*atan2( Camera_geo.y, Camera_geo.x ); // longitude


    /*
     * Only do hi-res imagery when we are close (or come up with a bettwer triggering scheme...)
     */
    if ( aInfo->Rcam < 1.4 ) {
// note: we currently arent making sure we load a whole new set each time
// we go between Full res and half res. This usually wont be a problem
// because they are numbered differently (i.e. zooming in/out at one region
// will cause quads with quite different names to be requested between the two sets.)
// Still, we should have a means to ensure there is no interference between the
// sets...

        // figure out what we need in the 3x3 grid centered on qi,qj
        if ( aInfo->Rcam > 1.2 ) {
            // Use Half-Res Quads
            nph = 40; nth = 20;
            qj = (int)(th/180.0*nth);
            qi = (int)((ph+180.0)/360.0*nph);
        } else {
            // Use Full-Res Quads
            nph = 80; nth = 40;
            qj = (int)(th/4.5);
            qi = (int)((ph+180.0)/4.5);
        }

        for ( nn=0, ii=qi-1; ii<=qi+1; ii++ ){
            iii = ii;
            if (iii<0) iii = nph-1;
            if (iii>=nph) iii = 0;
            for ( jj=qj-1; jj<=qj+1; jj++ ){
                if ( ( jj <= 0 ) || (jj >= (nth-1) ) ) QuadsNeeded[nn] = -1; // dont load the quads at the poles
                else QuadsNeeded[nn] = (long int)iii*100+jj;
                nn++;
            }
        }


        // Find all quads in Loaded that arent in Needed. These need to be destroyed.
        for ( ii=0; ii<9; ii++ ) {

            DestroyQuad = TRUE;

            for ( jj=0; jj<9; jj++ ) {
                if ( (QuadsLoaded[ii] == QuadsNeeded[jj]) || (QuadsLoaded[ii] < 0) ) {
                    DestroyQuad = FALSE;
                    break;
                }
            }
            if ( DestroyQuad ) DestroyHiResEarthQuad( &QuadsLoadedTexId[ii] );
            if (QuadsNeededDL[ii]>=0) glDeleteLists( QuadsLoadedDL[ii], 1 ); // destroy all the DLs

        }


        // Load or copy over all the qauds that are needed
        for ( ii=0; ii<9; ii++ ) {

            LoadQuad = TRUE;
            if ( QuadsNeeded[ii] < 0 ) {
                QuadsNeededDL[ii] = -1;
                QuadsNeededTexId[ii] = -1;
            } else {

                for ( jj=0; jj<9; jj++ ) {
                    if ( QuadsNeeded[ii] == QuadsLoaded[jj] ) {
//                        QuadsNeededDL[ii] = QuadsLoadedDL[jj];
                        QuadsNeededTexId[ii] = QuadsLoadedTexId[jj];
                        LoadQuad = FALSE;
                        break;
                    }
                }
                // If LoadQuad is true, then we need to load the texture. But Im thinking
                // we may need to reload all the DLs each time because they (can?) have time-dep
                // coords. So I re-coded CreateHiResEarthQuad() to take the LoadQuad flag and we
                // now call it every time.
                CreateHiResEarthQuad( LoadQuad, nph, nth, QuadsNeeded[ii], &QuadsNeededDL[ii], &QuadsNeededTexId[ii] );

            }

        }

        // finally copy over the new to the loaded
        for (ii=0; ii<9; ii++){
            QuadsLoaded[ii]      = QuadsNeeded[ii];
            QuadsLoadedDL[ii]    = QuadsNeededDL[ii];
            QuadsLoadedTexId[ii] = QuadsNeededTexId[ii];
        }


        for (ii=0; ii<9; ii++){
            if ( QuadsLoaded[ii] >= 0 ) glCallList( QuadsLoadedDL[ii] );
        }

    }





    glPopMatrix();






//glCallList( TopSideDL );
//glCallList( MeridPlane1DL );
//glCallList( MeridPlane2DL );












if (LightingStyle == 2){
//cgGLDisableProfile(myCgVertexProfile);
//checkForCgError("disabling vertex profile");
//cgGLDisableProfile(myCgFragmentProfile);
//checkForCgError("disabling fragment profile");
}





    if ( ShowAtmosphere ){

        /*
         * Compute the location of the Camera. It was originally at +Z, but its
         * been rotated by quat_view.
         */
        glEnable( GL_DEPTH_TEST );
        glEnable( GL_CULL_FACE );
        glEnable( GL_BLEND );
        glBlendFunc( GL_SRC_COLOR, GL_ONE_MINUS_SRC_COLOR );
        glPushMatrix();

            DrawSphereVertices( &aInfo->InnerSphere );

            glFrontFace(GL_CW);
            DrawSphereVertices( &aInfo->OuterSphere );
            glFrontFace(GL_CCW);

        glPopMatrix();

        glDisable( GL_BLEND );
        glDisable( GL_CULL_FACE );
//        glDisable( GL_DEPTH_TEST );

    }







// now that we have gone to a two-sided lighting model for the drift shells,
// we need to pay more attention to the face orientation. Check to see that they are really
// correct. Seems like there might be a problem with the CCW versus CW stuff....?

    glPushMatrix();
    glRotatef( RotAngle3, RotAxis3.x, RotAxis3.y, RotAxis3.z ); // This implements coord trans from Observer coords -> GSM

    //glUseProgram( g_shaderMyTest );

glFrontFace(GL_CW);
    glLightModelfv( GL_LIGHT_MODEL_TWO_SIDE, LightModelTwoSide);
    for (i=0; i<ObjInfo->MagEphemInfo->nAlpha; i++ ) {
        if ( ShowPitchAngle2[i] ) {
            glMaterialfv( GL_FRONT, GL_DIFFUSE,   gInfo->DriftShellMaterial[i].diffuse );
            glMaterialfv( GL_FRONT, GL_AMBIENT,   gInfo->DriftShellMaterial[i].ambient );
            glMaterialfv( GL_FRONT, GL_SPECULAR,  gInfo->DriftShellMaterial[i].specular );
            glMaterialf(  GL_FRONT, GL_SHININESS, gInfo->DriftShellMaterial[i].shininess*128.0 );
glMaterialfv( GL_BACK, GL_DIFFUSE,   gInfo->DriftShellMaterial[i].diffuse );
glMaterialfv( GL_BACK, GL_AMBIENT,   gInfo->DriftShellMaterial[i].ambient );
glMaterialfv( GL_BACK, GL_SPECULAR,  gInfo->DriftShellMaterial[i].specular );
glMaterialf(  GL_BACK, GL_SHININESS, gInfo->DriftShellMaterial[i].shininess*128.0 );
            if ( (int)(gInfo->DriftShellMaterial[i].diffuse[3]*128.0) < 128 ) {
                glEnable( GL_BLEND );
                glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
            }
            glCallList( ObjInfo->DriftShellList4 + i );
            if ( (int)(gInfo->DriftShellMaterial[i+1].diffuse[3]*128.0) < 128 ) {
                glDisable( GL_BLEND );
            }
        }
    }
    glLightModelfv( GL_LIGHT_MODEL_TWO_SIDE, LightModelOneSide);
glFrontFace(GL_CCW);


//    glCallList( ScPositionDL );

    //glUseProgram( 0 );


    glPopMatrix();


    /*
     * The sat positions and orbits are recomputed in the current coord system each time
     * (this is because the orbits need different times at each point in the orbit -- i.e
     *  a simple rotation wont do it.)
     */
    glCallList( SatsDL );
    glCallList( SatOrbitsDL );
    DrawSatLabels();


    int w = g_imageWidth;
    int h = g_imageHeight*0.1;
    if (h<100) h = 100;
    double rat = (double)h/(double)g_imageHeight;
    cairo_surface_t *cst = cairo_image_surface_create( CAIRO_FORMAT_ARGB32, w, h );
    cairo_t *cr  = cairo_create( cst );

    cairo_select_font_face( cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL );
    cairo_set_source_rgba( cr, 0.5, 0.5, 0.5, 0.0 );
    cairo_rectangle( cr, 0.0, 0.0, (double)w, (double)h );
    cairo_fill(cr);

    cairo_set_source_rgba( cr, 1.0, 1.0, 1.0, 0.9 );

    int y = 25.0;
    cairo_set_font_size( cr, 24.0 );
    sprintf(Str, "%s Coordinates", Coords[ObserverCoords] );
    cairo_move_to( cr, 10.0, y ); cairo_show_text( cr, Str);

    cairo_set_font_size( cr, 18.0 );

    y += 20;
    Lgm_Doy( CurrentDate, &ty, &tm, &td, &tD);
    sprintf(Str, "Date: %s %2d, %4d ( %4d / %03d )", MonthStr[CurrentMonth-1], CurrentDay, CurrentYear, CurrentYear, tD);
    cairo_move_to( cr, 15.0, y ); cairo_show_text( cr, Str);

    y += 20;
    Lgm_UT_to_hmsms( CurrentUT, &CurrentHour, &CurrentMin, &CurrentSec, &CurrentMilliSec );
    sprintf(Str, "Time: %02d:%02d:%02d.%03d UTC", CurrentHour, CurrentMin, CurrentSec, CurrentMilliSec );
    cairo_move_to( cr, 15.0, y ); cairo_show_text( cr, Str);

    y += 20;
    sprintf(Str, "Camera: (%6.2lf,%6.2lf,%6.2lf )Re", aInfo->Camera.x, aInfo->Camera.y, aInfo->Camera.z);
    cairo_move_to( cr, 15.0, y ); cairo_show_text( cr, Str);



    cairo_destroy( cr );



    glBindTexture( GL_TEXTURE_2D, Texture_TEMP );
    glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA, w, h, 0, GL_BGRA, GL_UNSIGNED_BYTE, (const GLvoid *)cairo_image_surface_get_data( cst ) );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );

    glMatrixMode(GL_MODELVIEW); glPushMatrix();
    glLoadIdentity();


    glMatrixMode( GL_PROJECTION ); glPushMatrix();
    glLoadIdentity();
    gluOrtho2D(-1.0, -1.0, 1.0, 1.0);
    glEnable( GL_TEXTURE_2D );
    glDisable( GL_LIGHTING );
    glEnable(GL_BLEND);
    glDisable(GL_DEPTH_TEST);
    glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );

    glColor4f( 1.0, 1.0, 1.0, 0.5 );
    glBegin(GL_QUADS);
    {
        glTexCoord2f( 0.0, 1.0 ); glVertex2f( -1.0, -1.0 );
        glTexCoord2f( 1.0, 1.0 ); glVertex2f(  1.0, -1.0 );
        glTexCoord2f( 1.0, 0.0 ); glVertex2f(  1.0,  -1.0+2.0*rat );
        glTexCoord2f( 0.0, 0.0 ); glVertex2f( -1.0,  -1.0+2.0*rat );
    }
    glEnd();

    glCallList( LogoDL );


    glEnable( GL_LIGHTING );
    glDisable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);
    glDisable( GL_TEXTURE_2D );

    glMatrixMode(GL_PROJECTION); glPopMatrix();

    glMatrixMode(GL_MODELVIEW); glPopMatrix();



    cairo_surface_destroy( cst );

    glCallList( EqPlaneGridDL );




}




gboolean expose_event( GtkWidget *widget, GdkEventExpose *event, gpointer data) {

    int             i;
    double          qq[4];
    double          TBQ[4]; // TB stands for trackball
    Lgm_Vector      vv;

    GdkGLContext    *glcontext = gtk_widget_get_gl_context (widget);
    GdkGLDrawable   *gldrawable = gtk_widget_get_gl_drawable (widget);


    if (!gdk_gl_drawable_gl_begin (gldrawable, glcontext)) return( FALSE );



    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );





    /*
     * Update trackball rotations/zooming
     */
    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity( );
    add_quats( view_quat_diff, view_quat, view_quat );
    add_quats( view_quat_diff3, view_quat3, view_quat3 );

//printf("view_quat = %g %g %g %g\n", view_quat[0], view_quat[1], view_quat[2], view_quat[3] );
//printf("view_quat3 = %g %g %g %g\n", view_quat3[0], view_quat3[1], view_quat3[2], view_quat3[3] );

    /*
     * Compute the location of the Camera. It was originally at +Z, but its
     * been rotated by quat_view3 and quat_view. And scaled by view_scale.
     */
    vv.x = vv.y = 0.0; vv.z = 1.0;
    for (i=0; i<4; i++) qq[i] = view_quat[i];
    Lgm_QuatRotateVector( qq, &vv, &aInfo->Camera );
    Lgm_ScaleVector( &aInfo->Camera, 10.0*view_scale );
    aInfo->nCamera = aInfo->Camera;
    aInfo->Rcam =  Lgm_NormalizeVector( &aInfo->nCamera );
    if (aInfo->Rcam < 1.02) {
        // dont go below surface
        aInfo->Camera = aInfo->nCamera;
        Lgm_ScaleVector( &aInfo->Camera, 1.02 );
        aInfo->Rcam = 1.02;
    }



    /*
     * Compute the location of the Up vector. It was originally +Y, but its
     * been rotated by quat_view3 and quat_view.
     */
    vv.x = vv.z = 0.0; vv.y = 1.0;
    for (i=0; i<4; i++) qq[i] = view_quat[i];
    Lgm_QuatRotateVector( qq, &vv, &aInfo->Up );
    Lgm_NormalizeVector( &aInfo->Up );

    /*
     * Compute the location of the Right vector. It was originally +X, but its
     * been rotated by quat_view3 and quat_view.
     */
    vv.x = vv.y = 0.0; vv.x = 1.0;
    for (i=0; i<4; i++) qq[i] = view_quat[i];
    Lgm_QuatRotateVector( qq, &vv, &aInfo->Right );
    Lgm_NormalizeVector( &aInfo->Right );




    TBQ[0] = view_quat[0]; TBQ[1] = view_quat[1]; TBQ[2] = view_quat[2]; TBQ[3] = view_quat[3];
    Lgm_QuatToAxisAngle(  TBQ, &RotAngle4,  &RotAxis4 );
    Lgm_NormalizeVector( &RotAxis4 );

    TBQ[0] = view_quat3[0]; TBQ[1] = view_quat3[1]; TBQ[2] = view_quat3[2]; TBQ[3] = view_quat3[3];
    Lgm_QuatToAxisAngle(  TBQ, &RotAngle5,  &RotAxis5 );

    glRotatef( RotAngle5, RotAxis5.x, RotAxis5.y, RotAxis5.z ); // This is the panning-type rotation
    glRotatef( -RotAngle4, RotAxis4.x, RotAxis4.y, RotAxis4.z ); // This is the spinning type rotation
    glTranslatef( -aInfo->Camera.x, -aInfo->Camera.y, -aInfo->Camera.z ); // Position things properly rel. to "camera" or "eye"


    glClearColor (0.0, 0.0, 0.0, 1.0);
//glClearColor( 1.0, 1.0, 1.0, 1.0 );
//glClearColor( 0.8, 0.8, 0.8, 0.8 );

    if ( LightingStyle == 0 ) {

        /*
         *  Solar illumination (light fixed at position of sun)
         */
        GLfloat ambient[]         = {0.0, 0.0, 0.0, 1.0};
        GLfloat position[]        = {5000.0, 0.0, 0.0, 0.0};
        GLfloat lmodel_ambient[]  = {0.3, 0.3, 0.3, 1.0};
        position[0] = Sun.x; position[1] = Sun.y; position[2] = Sun.z;
        glPushMatrix();
        glLightfv( GL_LIGHT0, GL_AMBIENT, ambient);
        glLightfv (GL_LIGHT0, GL_POSITION, position);
//        LightPosition[0] = 5000.0; LightPosition[1] = LightPosition[2] = LightPosition[3] = 0.0;
//        glLightfv (GL_LIGHT0, GL_POSITION, LightPosition );
        glLightModelfv( GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
        glPopMatrix();

    } else if ( LightingStyle == 1 ) {

        /*
         * Fixed illumination (light moves with viewer)
         */
        GLfloat ambient[]         = {0.0, 0.0, 0.0, 1.0};
        GLfloat position[]        = {0.0, 3000.0, 3000.0, 0.0};
        GLfloat lmodel_ambient[]  = {0.2, 0.2, 0.2, 1.0};
        position[0] = aInfo->Camera.x; position[1] = aInfo->Camera.y; position[2] = aInfo->Camera.z;
        glPushMatrix();
        glLightfv( GL_LIGHT0, GL_AMBIENT, ambient);
        LightPosition[0] = position[0]; LightPosition[1] = position[1]; LightPosition[2] = position[2];
        //glLightfv (GL_LIGHT0, GL_POSITION, position);
        glLightfv (GL_LIGHT0, GL_POSITION, LightPosition);
        glLightModelfv( GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
        glPopMatrix();

    } else {

        /*
         *  Solar illumination 2 (light fixed at position of sun but more ambient)
         */
        GLfloat ambient[]         = {0.7, 0.7, 0.7, 1.0};
        GLfloat position[]        = {5000.0, 0.0, 0.0, 0.0};
        GLfloat lmodel_ambient[]  = {0.7, 0.7, 0.7, 1.0};
        position[0] = Sun.x; position[1] = Sun.y; position[2] = Sun.z;
        glPushMatrix();
        glLightfv( GL_LIGHT0, GL_AMBIENT, ambient);
        glLightfv (GL_LIGHT0, GL_POSITION, position);
//        LightPosition[0] = 5000.0; LightPosition[1] = LightPosition[2] = LightPosition[3] = 0.0;
//        glLightfv (GL_LIGHT0, GL_POSITION, LightPosition );
        glLightModelfv( GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
        glPopMatrix();

    } 


    DrawScene();


    // experimental depth peeling
    // This may be too slow...
    // Could try dual depth peeling instead(?)
    //RenderFrontToBackPeeling();



    /* If we can, swap buffers */
    if ( gdk_gl_drawable_is_double_buffered( gldrawable ) ) {
        gdk_gl_drawable_swap_buffers( gldrawable );
    } else {
        glFlush( );
    }

    gdk_gl_drawable_gl_end( gldrawable );






    return( TRUE );

}




static float begin_x = 0.0;
static float begin_y = 0.0;

static float dx = 0.0;
static float dy = 0.0;

static gboolean button_press_event( GtkWidget *widget, GdkEventButton *event, gpointer data) {

    if (event->button == 1) {
        AnimateView = FALSE;
        view_quat_diff[0] = 0.0;
        view_quat_diff[1] = 0.0;
        view_quat_diff[2] = 0.0;
        view_quat_diff[3] = 1.0;
        AdjustIdle( widget );
    }

    begin_x = event->x;
    begin_y = event->y;

    return FALSE;
}




static gboolean button_release_event( GtkWidget *widget, GdkEventButton *event, gpointer data) {

    int diff=dx*dx + dy*dy;

    if (event->button == 1 && (diff >= ANIMATE_THRESHOLD)){
        AnimateView = TRUE;
    } else if ( event->button == 1 && (diff < ANIMATE_THRESHOLD ) ){
        // view animation rotations to zero
        view_quat_diff[0] = 0.0;
        view_quat_diff[1] = 0.0;
        view_quat_diff[2] = 0.0;
        view_quat_diff[3] = 1.0;
    } else if (event->button == 3 ){
        // dont animate the second trackball
        // that controls side to side viewing.
        view_quat_diff3[0] = 0.0;
        view_quat_diff3[1] = 0.0;
        view_quat_diff3[2] = 0.0;
        view_quat_diff3[3] = 1.0;
    }
    AdjustIdle( widget );

    dx = 0.0;
    dy = 0.0;

    return FALSE;
}

void trackball2( float q[4], Lgm_Vector *u, Lgm_Vector *r, float px, float py, float qx, float qy ) {

    float   a[3]; /* Axis of rotation */
    float   phi;  /* how much to rotate about axis */
    float  mag;

    Lgm_Vector  Up, Right;



    /*
     * Check for zero rotation
     */
    if ( (px == qx) && (py == qy) ) {
        q[0] = q[1] = q[2] = 0.0; q[3] = 1.0;
        return;
    }

    /*
     *  Find 2D vector defined by points
     */
    dx = qx - px; // x-component
    dy = qy - py; // y-component

    /*
     * normalize
     */
    mag = sqrt(dx*dx+dy*dy);
    dx /= mag;
    dy /= mag;

    phi = mag * -90.0*M_PI/180.0; // -45-deg


    /*
     *  Now, we actually want to rotate about an axis thats perp to this.
     *  Above, we have bvec = dx xhat + dy yhat. To get the perp vector
     *  rot by 90deg. This gives, avec = -dy xhat + dx yhat
     *
     *  Here, xhat is the Right vector (r), and yhat if the Up vector (u)
     */
//    Up = *u; Right = *r;
// do the rotation before other rots -- its easier
Up.x = Up.z = 0.0; Up.y = 1.0;
Right.y = Right.z = 0.0; Right.x = 1.0;
    Lgm_ScaleVector( &Right, -dy ); // -dy xhat
    Lgm_ScaleVector( &Up,     dx ); //  dx yhat
    a[0] = Right.x + Up.x; a[1] = Right.y + Up.y; a[2] = Right.z + Up.z;


    axis_to_quat( a, phi, q );
}



static gboolean motion_notify_event( GtkWidget *widget, GdkEventMotion *event, gpointer data) {

    float w = widget->allocation.width;
    float h = widget->allocation.height;
    float x = event->x;
    float y = event->y;
    gboolean redraw = FALSE;

    /* Rotation. */
    if (event->state & GDK_BUTTON1_MASK) {
      trackball( view_quat_diff, (2.0 * begin_x - w) / w, (h - 2.0 * begin_y) / h, (2.0 * x - w) / w, (h - 2.0 * y) / h, aInfo->Rcam );
      dx = x - begin_x;
      dy = y - begin_y;
      redraw = TRUE;
    }

    /* Scaling. */
    if (event->state & GDK_BUTTON2_MASK) {
        view_scale = view_scale * (1.0 + (y - begin_y) / h);
        if (view_scale > VIEW_SCALE_MAX) {
	        view_scale = VIEW_SCALE_MAX;
        } else if (view_scale < VIEW_SCALE_MIN) {
	        view_scale = VIEW_SCALE_MIN;
        }
        redraw = TRUE;
    }

    /* Second Rotation -- allows off-center viewing */
    if (event->state & GDK_BUTTON3_MASK) {
          trackball2( view_quat_diff3, &aInfo->Up, &aInfo->Right, (2.0 * begin_x - w) / w, (h - 2.0 * begin_y) / h, (2.0 * x - w) / w, (h - 2.0 * y) / h );
          dx = x - begin_x;
          dy = y - begin_y;
          redraw = TRUE;
    }


    begin_x = x;
    begin_y = y;

    if ( redraw && !AnimateView ) gdk_window_invalidate_rect( widget->window, &widget->allocation, FALSE );

    return TRUE;
}



static gboolean idle( GtkWidget *widget ) {

    static int  AllowSave = TRUE;
    double      oJD;
    char        Str[256];
    Lgm_CTrans *c = Lgm_init_ctrans( 0 );


    if ( AnimateTime != TIME_REALTIMEPLAY ) {
        if      ( (AnimateTime == TIME_PLAY_FOREWARD) && (CurrentJD >= EndJD+1e-8))   AnimateTime  = TIME_STOP;
        else if ( (AnimateTime == TIME_PLAY_BACKWARD) && (CurrentJD <= StartJD-1e-8)) AnimateTime  = TIME_STOP;
        else if ( (AnimateTime == TIME_STEP_FOREWARD) && (CurrentJD >= EndJD+1e-8))   AnimateTime  = TIME_STOP;
        else if ( (AnimateTime == TIME_STEP_BACKWARD) && (CurrentJD <= StartJD-1e-8)) AnimateTime  = TIME_STOP;

        if ( AnimateTime == TIME_STOP ) {
            AdjustIdle( drawing_area );
        }
    }


    if (    (AnimateTime == TIME_STEP_FOREWARD) || (AnimateTime == TIME_STEP_BACKWARD)
         || (AnimateTime == TIME_PLAY_FOREWARD) || (AnimateTime == TIME_PLAY_BACKWARD)
         || (AnimateTime == TIME_REALTIMEPLAY) ) {


        ++cFrame;
        oJD = CurrentJD;
        if ( AnimateTime == TIME_REALTIMEPLAY ) {
            nFramesLeft = 999999;
            CurrentJD = Lgm_GetCurrentJD( c );
        } else {
            --nFramesLeft;
            CurrentJD += TimeInc/86400.0;
        }

        if ( oJD != CurrentJD ){
            Lgm_jd_to_ymdh( CurrentJD, &CurrentDate, &CurrentYear, &CurrentMonth, &CurrentDay, &CurrentUT );
            Lgm_UT_to_hmsms( CurrentUT, &CurrentHour, &CurrentMin, &CurrentSec, &CurrentMilliSec );
            UpdateTimeDepQuants( CurrentDate, CurrentUT );
//MGH MGH MGH MGH UpdateTimeDepQuants( ObjInfo->MagEphemInfo->Date, ObjInfo->MagEphemInfo->UTC  );

            if (ShowStars) ReLoadStars( );
            ReCreateSats();
            ReCreateSatOrbits();
            ReCreateEarth( );
            ReCreateMoon( );
            ReCreateThemisFovs();
        }

        expose_event( drawing_area, NULL, NULL );

        
        if ( DumpFrames ){
            if ( (CurrentSec >  59.5) || (CurrentSec < 0.5) ){
                if ( AllowSave ) {
                    int x, y, width, height, depth;
                    GdkPixbuf *pixbuf;
                    char PngFile[40];
                    gdk_window_get_geometry( widget->window, &x, &y, &width, &height, &depth );
                    pixbuf = gdk_pixbuf_get_from_drawable( NULL, GDK_DRAWABLE(widget->window), NULL, 0, 0, 0, 0, width, height);
                    sprintf( PngFile, "%04ld.png", cFrame );
                    //sprintf( PngFile, "Latest.png" );
printf("Writing: %s\n", PngFile );
                    gdk_pixbuf_save( pixbuf, PngFile, "png", NULL, "compression", "0", NULL);
                    g_object_unref( pixbuf );
                    AllowSave = FALSE;
                }
            } else {
                AllowSave = TRUE;
            }
        }


        if (oJD != CurrentJD){
            sprintf(Str, "Current Time: %4d/%02d/%02d  %02d:%02d:%02d.%03d\n  UT = %.8lf  JD = %.8lf", CurrentYear, CurrentMonth, CurrentDay,
                                CurrentHour, CurrentMin, CurrentSec, CurrentMilliSec, CurrentUT, CurrentJD );
            gtk_misc_set_alignment(GTK_MISC(CurrentTimeLabel), 0, 0.5);
            gtk_label_set_text( GTK_LABEL(CurrentTimeLabel), Str );

            sprintf(Str, "Frames Done: %ld    Frames Remaining: %ld    (Total: %ld)", cFrame, nFramesLeft, nFrames);
            gtk_label_set_text( GTK_LABEL(cFramesLabel), Str );
        }

        if ( (AnimateTime == TIME_STEP_FOREWARD) || (AnimateTime == TIME_STEP_BACKWARD) ) {
            AnimateTime = TIME_STOP;
            AdjustIdle( drawing_area );
        }

    } else if ( AnimateView ) {
        expose_event( drawing_area, NULL, NULL );
    }

    /* Update synchronously. */
// check this one out...
//    gdk_window_process_updates( widget->window, FALSE );

    Lgm_free_ctrans( c );

    return TRUE;

}


static guint idle_id = 0;

static void idle_add (GtkWidget *widget) {
    if (idle_id == 0) {
        idle_id = g_idle_add_full( GDK_PRIORITY_REDRAW, (GSourceFunc)idle, widget, NULL );
    }
}

static void idle_remove (GtkWidget *widget) {
  if (idle_id != 0) {
      g_source_remove( idle_id );
      idle_id = 0;
    }
}


gboolean quit_app( GtkWidget *widget, gpointer data ) {

    idle_remove( widget );
    gtk_widget_destroy( ViewDriftShellWindow );
    gtk_main_quit();

    return TRUE;
}



static gboolean key_press_event( GtkWidget *widget, GdkEventKey *event, gpointer data) {

    switch (event->keyval) {

        case GDK_Escape:
            quit_app( ViewDriftShellWindow, NULL );
            break;

        default:
            return FALSE;
    }

    return TRUE;

}

static gboolean map_event( GtkWidget *widget, GdkEvent *event, gpointer data ) {
  if (AnimateView) idle_add( widget );
  return TRUE;
}


static gboolean unmap_event( GtkWidget *widget, GdkEvent *event, gpointer data ) {
    idle_remove( widget );
    return TRUE;
}




static gboolean visibility_notify_event( GtkWidget *widget, GdkEventVisibility *event, gpointer data) {

    if ( AnimateView ) {
        if (event->state == GDK_VISIBILITY_FULLY_OBSCURED){
	        idle_remove( widget );
        } else {
	        idle_add( widget );
        }
    }

    return TRUE;
}


/*
 *  Examine states for AnimateView and AnimateTime
 *  and perform actions accordingly.
 */
static void AdjustIdle( GtkWidget *widget ){

    int RunIdle = FALSE;

    // idle should ben running for any of these cases
    if ( AnimateView || (AnimateTime == TIME_PLAY_FOREWARD) || (AnimateTime == TIME_PLAY_BACKWARD)
            || (AnimateTime == TIME_STEP_FOREWARD) || (AnimateTime == TIME_STEP_BACKWARD)
            || (AnimateTime == TIME_REALTIMEPLAY) ) {

        RunIdle = TRUE;
    }

    if ( RunIdle ) {
        //start it (idle_add wont do anythinf if its already running)
        idle_add( widget );
    } else if ( !RunIdle ) {
        //stop it (idle_add wont do anythinf if its already stopped)
        idle_remove( widget );
//        gdk_window_invalidate_rect( widget->window, &widget->allocation, FALSE );
    }

}



/* Toggle animation.*/
static void toggle_animation( GtkWidget *widget ) {


    if ( !AnimateView ) {
        // view animation rotations to zero
        view_quat_diff[0] = 0.0;
        view_quat_diff[1] = 0.0;
        view_quat_diff[2] = 0.0;
        view_quat_diff[3] = 1.0;
    }

    if ( AnimateView || AnimateTime ) {
        // keep idle going if either view or time is being animated
        idle_add( widget );
    } else {
        // remove idle calls if neither view nor time is being animated
        idle_remove( widget );
        printf("FOO\n"); gdk_window_invalidate_rect( widget->window, &widget->allocation, FALSE );
    }

}


void  ChangeMapImage( GtkFileChooser *chooser,  gpointer user_data){

    char *filename;
    filename = gtk_file_chooser_get_filename( chooser );
    if ( filename == NULL ) return;

    strcpy( MapImageFilename, filename );
    printf( "******************************************** You selected file: %s\n", MapImageFilename );
    LoadTextures();

    return;

}




static void ChangeTimeInc( GtkWidget  *widget, gpointer data ) {

    char    Str[256];
    double  t;

    t = gtk_spin_button_get_value( GTK_SPIN_BUTTON(widget) );
    TimeInc = ((long int)( t*1000.0 )) / 1000.0;

    if ( CurrentJD >= EndJD ) {
        nFramesLeft = 0;
    } else {
        nFramesLeft = (long int)((EndJD - CurrentJD) / (TimeInc/86400.0) +1.5);
    }
    nFrames = cFrame + nFramesLeft;

    sprintf(Str, "Frames Done: %ld    Frames Remaining: %ld    (Total: %ld)", cFrame, nFramesLeft, nFrames);
    gtk_label_set_text( GTK_LABEL(cFramesLabel), Str );

}

static void ToggleDumpFrames( GtkWidget  *widget, gpointer data ) {
    DumpFrames = gtk_toggle_button_get_active(  GTK_TOGGLE_BUTTON( DumpFramesCheckbutton ) );
}




void SetStartDate( Lgm_DateTime *dt ) {


    gtk_spin_button_set_value( GTK_SPIN_BUTTON(StartYearSpinbutton),  dt->Year  );
    gtk_spin_button_set_value( GTK_SPIN_BUTTON(StartMonthSpinbutton), dt->Month );
    gtk_spin_button_set_value( GTK_SPIN_BUTTON(StartDaySpinbutton),   dt->Day   );


}



static void ChangeStartOrEndDate( GtkWidget  *widget, gpointer data ) {

    char        Str[256];
    int         k, MaxDay;
    Lgm_CTrans *c = Lgm_init_ctrans( 0 );

    k = GPOINTER_TO_INT( data );
    switch ( k ) {
        case 1:
                StartYear = (int) gtk_spin_button_get_value( GTK_SPIN_BUTTON(widget) );
                // changes needed for leap years
                if ( StartMonth == 2 ) {
                    MaxDay = 28+Lgm_LeapYear( StartYear );
                    gtk_adjustment_set_upper( GTK_ADJUSTMENT( StartDaySpinbutton_adj ), MaxDay );
                    if ( StartDay > MaxDay ) gtk_spin_button_set_value( GTK_SPIN_BUTTON(StartDaySpinbutton), MaxDay );
                }
                break;
        case 2:
                StartMonth = (int) gtk_spin_button_get_value( GTK_SPIN_BUTTON(widget) );
                MaxDay = MonthDays[StartMonth-1] + ( (StartMonth == 2) ? Lgm_LeapYear(StartYear) : 0 );
                gtk_adjustment_set_upper( GTK_ADJUSTMENT( StartDaySpinbutton_adj ), MaxDay );
                if ( StartDay > MaxDay ) gtk_spin_button_set_value( GTK_SPIN_BUTTON(StartDaySpinbutton), MaxDay );
                break;
        case 3:
                StartDay = (int) gtk_spin_button_get_value( GTK_SPIN_BUTTON(widget) );
                break;

        case 11:
                EndYear = (int) gtk_spin_button_get_value( GTK_SPIN_BUTTON(widget) );
                // changes needed for leap years
                if ( EndMonth == 2 ) {
                    MaxDay = 28+Lgm_LeapYear( EndYear );
                    gtk_adjustment_set_upper( GTK_ADJUSTMENT( EndDaySpinbutton_adj ), MaxDay );
                    if ( EndDay > MaxDay ) gtk_spin_button_set_value( GTK_SPIN_BUTTON(EndDaySpinbutton), MaxDay );
                }
                break;
        case 12:
                EndMonth = (int) gtk_spin_button_get_value( GTK_SPIN_BUTTON(widget) );
                MaxDay = MonthDays[EndMonth-1] + ( (EndMonth == 2) ? Lgm_LeapYear(EndYear) : 0 );
                gtk_adjustment_set_upper( GTK_ADJUSTMENT( EndDaySpinbutton_adj ), MonthDays[EndMonth-1] );
                if ( EndDay > MaxDay ) gtk_spin_button_set_value( GTK_SPIN_BUTTON(StartDaySpinbutton), MaxDay );
                break;
        case 13:
                EndDay = (int) gtk_spin_button_get_value( GTK_SPIN_BUTTON(widget) );
                break;

    }

    if ( k < 10 ) {
        StartDate = StartYear*10000 + StartMonth*100 + StartDay;
        StartJD = Lgm_JD( StartYear, StartMonth, StartDay, StartUT, LGM_TIME_SYS_UTC, c );
        printf("StartDate = %ld\n", StartDate);
        sprintf(Str, "Start Time: %4d/%02d/%02d  %02d:%02d:%02d.%03d\n  UT = %.8lf  JD = %.8lf", StartYear, StartMonth, StartDay,
                                    StartHour, StartMin, StartSec, StartMilliSec, StartUT, StartJD );
        gtk_misc_set_alignment(GTK_MISC(StartTimeLabel), 0, 0.5);
        gtk_label_set_text( GTK_LABEL(StartTimeLabel), Str );
    } else {
        EndDate = EndYear*10000 + EndMonth*100 + EndDay;
        EndJD = Lgm_JD( EndYear, EndMonth, EndDay, EndUT, LGM_TIME_SYS_UTC, c );
        printf("EndDate = %ld\n", EndDate);
        sprintf(Str, "End Time: %4d/%02d/%02d  %02d:%02d:%02d.%03d\n  UT = %.8lf  JD = %.8lf", EndYear, EndMonth, EndDay,
                                    EndHour, EndMin, EndSec, EndMilliSec, EndUT, EndJD );
        gtk_misc_set_alignment(GTK_MISC(EndTimeLabel), 0, 0.5);
        gtk_label_set_text( GTK_LABEL(EndTimeLabel), Str );
    }

    if ( CurrentJD >= EndJD ) {
        nFramesLeft = 0;
    } else {
        nFramesLeft = (long int)((EndJD - CurrentJD) / (TimeInc/86400.0) +1.5);
    }
    nFrames = cFrame + nFramesLeft;
    sprintf(Str, "Frames Done: %ld    Frames Remaining: %ld    (Total: %ld)", cFrame, nFramesLeft, nFrames);
    gtk_label_set_text( GTK_LABEL(cFramesLabel), Str );

    Lgm_free_ctrans( c );

}

static void ChangeStartOrEndTime( GtkWidget  *widget, gpointer data ) {

    char        Str[256];
    int         k;
    Lgm_CTrans *c = Lgm_init_ctrans( 0 );

    k = GPOINTER_TO_INT( data );
    switch ( k ) {
        case 1:
                StartHour = (int) gtk_spin_button_get_value( GTK_SPIN_BUTTON(widget) );
                break;
        case 2:
                StartMin = (int) gtk_spin_button_get_value( GTK_SPIN_BUTTON(widget) );
                break;
        case 3:
                StartSec = gtk_spin_button_get_value( GTK_SPIN_BUTTON(widget) );
                break;

        case 11:
                EndHour = (int) gtk_spin_button_get_value( GTK_SPIN_BUTTON(widget) );
                break;
        case 12:
                EndMin = (int) gtk_spin_button_get_value( GTK_SPIN_BUTTON(widget) );
                break;
        case 13:
                EndSec = gtk_spin_button_get_value( GTK_SPIN_BUTTON(widget) );
                break;

    }

    if ( k < 10 ) {
        StartUT = (double)StartHour + (double)StartMin/60.0 + (double)StartSec/3600.0;
        StartJD = Lgm_JD( StartYear, StartMonth, StartDay, StartUT, LGM_TIME_SYS_UTC, c );
        printf("StartTime = %lf\n", StartUT);
        sprintf(Str, "Start Time: %4d/%02d/%02d  %02d:%02d:%02d.%03d\n  UT = %.8lf  JD = %.8lf", StartYear, StartMonth, StartDay,
                                    StartHour, StartMin, StartSec, StartMilliSec, StartUT, StartJD );
printf("%s\n", Str);
        gtk_label_set_text( GTK_LABEL(StartTimeLabel), Str );
    } else {
        EndUT = (double)EndHour + (double)EndMin/60.0 + (double)EndSec/3600.0;
        EndJD = Lgm_JD( EndYear, EndMonth, EndDay, EndUT, LGM_TIME_SYS_UTC, c );
        printf("EndTime = %lf\n", EndUT);
        sprintf(Str, "End Time: %4d/%02d/%02d  %02d:%02d:%02d.%03d\n  UT = %.8lf  JD = %.8lf", EndYear, EndMonth, EndDay,
                                    EndHour, EndMin, EndSec, EndMilliSec, EndUT, EndJD );
        gtk_label_set_text( GTK_LABEL(EndTimeLabel), Str );
    }

    if ( CurrentJD >= EndJD ) {
        nFramesLeft = 0;
    } else {
        nFramesLeft = (long int)((EndJD - CurrentJD) / (TimeInc/86400.0) +1.5);
    }
    nFrames = cFrame + nFramesLeft;
    sprintf(Str, "Frames Done: %ld    Frames Remaining: %ld    (Total: %ld)", cFrame, nFramesLeft, nFrames);
    gtk_label_set_text( GTK_LABEL(cFramesLabel), Str );

    Lgm_free_ctrans( c );


}


void TimeAction( GtkWidget *widget, gpointer data ) {

    char    Str[256];
    int     k;

    k = GPOINTER_TO_INT( data );
    switch ( k ) {
        case TIME_RESET_BACKWARD_TO_START:
        case TIME_RESET_FOREWARD_TO_END:

                if ( k == TIME_RESET_BACKWARD_TO_START ) {
                    printf("Reset Current Time to Start Time...\n");
                    CurrentYear  = StartYear; CurrentMonth = StartMonth; CurrentDay = StartDay;
                    CurrentDate  = StartDate;
                    CurrentHour  = StartHour; CurrentMin = StartMin; CurrentSec = StartSec;
                    CurrentUT    = StartUT;
                    CurrentJD    = StartJD;
                } else {
                    printf("Reset Current Time to End Time...\n");
                    CurrentYear  = EndYear; CurrentMonth = EndMonth; CurrentDay  = EndDay;
                    CurrentDate  = EndDate;
                    CurrentHour  = EndHour; CurrentMin = EndMin; CurrentSec = EndSec;
                    CurrentUT    = EndUT;
                    CurrentJD    = EndJD;
                }

                RunTime     = FALSE;
                AnimateTime = TIME_STOP;
if (k== TIME_RESET_BACKWARD_TO_START)
printf("Back to Start CurrentDate, CurrentUT = %ld %g\n", CurrentDate, CurrentUT);
                UpdateTimeDepQuants( CurrentDate, CurrentUT );

                cFrame = 0;
                if ( CurrentJD >= EndJD ) { nFramesLeft = 0;
                } else {                    nFramesLeft = (long int)((EndJD - CurrentJD) / (TimeInc/86400.0) +1.5); }
                nFrames = cFrame + nFramesLeft;

                ReCreateEarth( );
                ReCreateMoon( );
                if (ShowThemisFovs) ReCreateThemisFovs();
                if (ShowStars) ReLoadStars( );
                ReCreateSats();
                ReCreateSatOrbits();
                printf("A\n"); expose_event( drawing_area, NULL, NULL );

                sprintf(Str, "Current Time: %4d/%02d/%02d  %02d:%02d:%02d.%03d\n  UT = %.8lf  JD = %.8lf", CurrentYear, CurrentMonth, CurrentDay,
                                                    CurrentHour, CurrentMin, CurrentSec, CurrentMilliSec, CurrentUT, CurrentJD );
                gtk_misc_set_alignment(GTK_MISC(CurrentTimeLabel), 0, 0.5);
                gtk_label_set_text( GTK_LABEL(CurrentTimeLabel), Str );

                //sprintf(Str, "Current Frame: %ld / %ld", cFrame, nFrames);
                sprintf(Str, "Frames Done: %ld    Frames Remaining: %ld    (Total: %ld)", cFrame, nFramesLeft, nFrames);
                gtk_label_set_text( GTK_LABEL(cFramesLabel), Str );

                break;

        case TIME_STEP_FOREWARD:
                printf("Step Current Time foreward by TimeInc seconds...\n");
                StepTime    = TRUE;
                RunTime     = FALSE;
                AnimateTime = TIME_STEP_FOREWARD;
                AdjustIdle( drawing_area );
//                if (!AnimateView) toggle_animation( drawing_area );

                break;

        case TIME_STEP_BACKWARD:
                printf("Step Current Time back by TimeInc seconds...\n");

                break;

        case TIME_PLAY_FOREWARD:
                printf("Run to the end in steps of TimeInc seconds...\n");
                RunTime     = TRUE;
                AnimateTime = TIME_PLAY_FOREWARD;
                AdjustIdle( drawing_area );
//                if (!AnimateView) toggle_animation( drawing_area );
                break;

        case TIME_PLAY_BACKWARD:
                printf("Run backward to the start in steps of TimeInc seconds...\n");
                RunTime     = TRUE;
                AnimateTime = TIME_PLAY_BACKWARD;
                AdjustIdle( drawing_area );
//                if (!AnimateView) toggle_animation( drawing_area );
                break;

        case TIME_STOP:
                printf("Stop ...\n");
                RunTime     = FALSE;
                StepTime    = FALSE;
                AnimateTime = TIME_STOP;
                AdjustIdle( drawing_area );
                break;

        case TIME_REALTIMEPLAY:
                printf("Run in real time mode...\n");
                RunTime     = TRUE;
                AnimateTime = TIME_REALTIMEPLAY;
                AdjustIdle( drawing_area );
                break;

        case 5:
                UpdateTimeDepQuants( CurrentDate, CurrentUT );

                ReCreateEarth( );
                ReCreateMoon( );
                if (ShowThemisFovs) ReCreateThemisFovs();
                if (ShowStars) ReLoadStars( );
                ReCreateSats();
                ReCreateSatOrbits();
                printf("3.\n"); expose_event( drawing_area, NULL, NULL );


                break;
    }


}
static void TLE_TraceFieldLines( GtkWidget *widget, gpointer data) {

    printf("HEre I am\n");


}

static void ChangeTLE( GtkWidget *widget, gpointer data) {
    int     i;
    double  val;
    char    Filename[1024];

    i   = GPOINTER_TO_INT( data );
    val = gtk_spin_button_get_value( GTK_SPIN_BUTTON(widget) );

    switch ( i ) {  
        
        case 0:
                printf("Setting Name\n");
                break;

        case 1:
                printf("Setting SatNum\n");
                tle.SatNum = (int)val;
                break;

        case 2:
                printf("Setting IntDesig\n");
                break;

        case 3:
                printf("Setting Epoch\n");
                tle.Epoch = val;
                break;

        case 4:
                printf("Setting d1MeanMotion\n");
                tle.d1MeanMotion = val;
                break;

        case 5:
                printf("Setting d2MeanMotion\n");
                tle.d2MeanMotion = val;
                break;

        case 6:
                printf("Setting BSTAR\n");
                tle.BSTAR = val;
                break;

        case 7:
                printf("Setting ElemNum\n");
                tle.ElemNum = (int)val;
                break;

        case 8:
                printf("Setting Inc\n");
                tle.Inc = val;
                break;

        case 9:
                printf("Setting RAAN\n");
                tle.RAAN = val;
                break;

        case 10:
                printf("Setting Ecc\n");
                tle.Ecc = val;
                break;

        case 11:
                printf("Setting AoP\n");
                tle.AoP = val;
                break;

        case 12:
                printf("Setting MeanAnomaly\n");
                tle.MeanAnomaly = val;
                break;

        case 13:
                printf("Setting MeanMotion\n");
                tle.MeanMotion = val;
                break;

        case 14:
                printf("Setting RevNum\n");
                tle.RevNum = (int)val;
                break;

    }

    char Str[80], Str2[80];

    //strcpy( tle.Line1, "1 23455U 94089A   97320.90946019  .00000140  00000-0  10191-3 0  2621");
    for (i=0; i<69; i++ ) tle.Line1[i] = ' '; tle.Line1[69] = '\0';
    tle.Line1[0] = '1';
    tle.Line1[1] = ' ';

    sprintf( Str, "%5ldU", tle.SatNum );
    for (i=0; i<6; i++ ) tle.Line1[2+i] = Str[i];

    sprintf( Str, "%-8s", tle.IntDesig );
    for (i=0; i<8; i++ ) tle.Line1[9+i] = Str[i];


    int iE = (int)tle.Epoch;
    long int fE = (long int)((tle.Epoch - iE)*1.0e8 +.5);

    sprintf( Str, "%05d", iE );
    for (i=0; i<5; i++ ) tle.Line1[18+i] = Str[i];
    tle.Line1[23] = '.';
    sprintf( Str, "%08d", fE );
    for (i=0; i<8; i++ ) tle.Line1[24+i] = Str[i];



    sprintf( Str, "%.8lf", tle.d1MeanMotion );
    if ( tle.d1MeanMotion >= 0.0 ) {
        Str[0] = ' ';
        strcpy( Str2, Str );
    } else {
        Str[0] = ' ';
        Str[1] = '-';
        strcpy( Str2, Str+1 );
    }
    for (i=0; i<10; i++ ) tle.Line1[33+i] = Str2[i];

    tle.Line1[62] = '0';


    BSTAR_to_STR( tle.d2MeanMotion, Str2 );
    for (i=0; i<8; i++ ) tle.Line1[44+i] = Str2[i];

    BSTAR_to_STR( tle.BSTAR, Str2 );
    for (i=0; i<8; i++ ) tle.Line1[53+i] = Str2[i];



    sprintf( Str, "%4d", tle.ElemNum );
    for (i=0; i<4; i++ ) tle.Line1[64+i] = Str[i];


//strcpy( tle.Line1, "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927");
//strcpy( tle.Line1, "1 23455U 94089A   97320.90946019  .00000140  00000-0  10191-3 0  2621");
//strcpy( tle.Line1, "2 23455  99.0090 272.6745 0008546 223.1686 136.8816 14.11711747148495");
//1 23455U 94089A   97320.90946019  .00000140  00000-0  10191-3 0  2621
//2 23455  99.0090 272.6745 0008546 223.1686 136.8816 14.11711747148495

    int sum = 0;
    for (i=0; i<68; i++ ) {
        if ( isdigit(tle.Line1[i]) ){ Str[0] = tle.Line1[i]; Str[1] = '\0';
            sum += atoi( Str );
        } else if ( tle.Line1[i] == '-' ) {
            ++sum;
        }
        
    }
    printf( "sum = %d sum%10 = %d\n", sum, sum%10 );

    tle.Line1[68] = '0' + sum%10;
    tle.Line1[69] = '\0';




    /*
     *    01      Line Number of Element Data
     *    03-07   Satellite Number
     *    09-16   Inclination [Degrees]
     *    18-25   Right Ascension of the Ascending Node [Degrees]
     *    27-33   Eccentricity (decimal point assumed)
     *    35-42   Argument of Perigee [Degrees]
     *    44-51   Mean Anomaly [Degrees]
     *    53-63   Mean Motion [Revs per day]
     *    64-68   Revolution number at epoch [Revs]
     *    69      Checksum (Modulo 10)
     */
    //strcpy( tle.Line2, "2 23455  99.0090 272.6745 0008546 223.1686 136.8816 14.11711747148495" );
    printf("Sample Line2\n2 23455  99.0090 272.6745 0008546 223.1686 136.8816 14.11711747148495\n" );
    for (i=0; i<69; i++ ) tle.Line2[i] = ' '; tle.Line2[69] = '\0';
    tle.Line2[0] = '2';
    tle.Line2[1] = ' ';
    
    sprintf( Str, "%5ldU", tle.SatNum );
    for (i=0; i<6; i++ ) tle.Line2[2+i] = Str[i];

    sprintf( Str, "%8.4lf", tle.Inc );
    for (i=0; i<8; i++ ) tle.Line2[8+i] = Str[i];

    sprintf( Str, "%8.4lf", tle.RAAN );
    for (i=0; i<8; i++ ) tle.Line2[17+i] = Str[i];

    sprintf( Str, "%8.7lf", tle.Ecc );
    for (i=0; i<7; i++ ) tle.Line2[26+i] = Str[i+2];

    sprintf( Str, "%8.4lf", tle.AoP );
    for (i=0; i<8; i++ ) tle.Line2[34+i] = Str[i];

    sprintf( Str, "%8.4lf", tle.MeanAnomaly );
    for (i=0; i<8; i++ ) tle.Line2[43+i] = Str[i];

    sprintf( Str, "%11.8lf", tle.MeanMotion );
    for (i=0; i<11; i++ ) tle.Line2[52+i] = Str[i];

    sprintf( Str, "%05d", tle.RevNum );
    for (i=0; i<5; i++ ) tle.Line2[63+i] = Str[i];

    sum = 0;
    for (i=0; i<68; i++ ) {
        if ( isdigit(tle.Line2[i]) ){ Str[0] = tle.Line2[i]; Str[1] = '\0';
            sum += atoi( Str );
        } else if ( tle.Line2[i] == '-' ) {
            ++sum;
        }
        
    }
    printf( "sum = %d sum%10 = %d\n", sum, sum%10 );

    tle.Line2[68] = '0' + sum%10;

    strcpy( tle.Line0, "0 TLE" );



    printf("Line0: |%s|\n", tle.Line0 );
    printf("Line1: |%s|\n", tle.Line1 );
    printf("Line2: |%s|\n", tle.Line2 );

    sprintf( Filename, "%s/SAT_GROUPS/TLE.txt", getenv("HOME") );
    FILE *fp_tle = fopen( Filename, "w" );
    fprintf( fp_tle, "%s\n", tle.Line0 );
    fprintf( fp_tle, "%s\n", tle.Line1 );
    fprintf( fp_tle, "%s\n", tle.Line2 );
    fclose( fp_tle );




    int     nTLEs;
    _SgpTLE TLEs[1];
    nTLEs = 0; // must set this
    LgmSgp_ReadTlesFromStrings( tle.Line0, tle.Line1, tle.Line2, &nTLEs, &TLEs[0], 1 );

    _GroupNode      *g;
    _SpaceObjects   *Group;

    g = SatSelectorInfo->SatGroupList;
    while ( g != NULL ) {
        Group = g->Group;

        for (i=0; i<Group->nSat; i++){
            if ( strcmp( Group->Sat[i].TLE.Name , "0 TLE" ) == 0 ) {
                printf("Group->Sat[i].TLE.Name = %s\n", Group->Sat[i].TLE.Name);
                printf("FOUND IT. i = %d\n", i);
                Group->Sat[i].TLE = TLEs[0];;
            }
        }

        g = g->Next;
        
    }

    sprintf( Str, "<b><tt><big><span>%s</span></big></tt></b>", TLEs[0].Line0 ); gtk_label_set_markup( GTK_LABEL(tle.Line0Label), Str ); 
    sprintf( Str, "<b><tt><big><span>%s</span></big></tt></b>", TLEs[0].Line1 ); gtk_label_set_markup( GTK_LABEL(tle.Line1Label), Str ); 
    sprintf( Str, "<b><tt><big><span>%s</span></big></tt></b>", TLEs[0].Line2 ); gtk_label_set_markup( GTK_LABEL(tle.Line2Label), Str ); 




//    ReLoadSpaceObjects();
    ReCreateSats();
    ReCreateSatOrbits();

//    view_quat[i] = val;
    expose_event( drawing_area, NULL, NULL );

}

static void ChangeViewQuat( GtkWidget *widget, gpointer data) {
    int     i;
    double  val;

    i   = GPOINTER_TO_INT( data );
    val = gtk_spin_button_get_value( GTK_SPIN_BUTTON(widget) );

    view_quat[i] = val;
    expose_event( drawing_area, NULL, NULL );
}


static void ChangenFieldPnts( GtkWidget *widget, gpointer data) {
    int n;

    // get new n value
    n = gtk_spin_button_get_value( GTK_SPIN_BUTTON(widget) );
    printf("n = %d\n", n);

    // re-create field lines and drift shell meshes with new n
    MakeFieldLines( n, ObjInfo );
    MakeDriftShellMesh( ObjInfo );


    // re-generate display lists
    //ReGenerateFieldLineLists( ObjInfo );
    ReGenerateDriftShellLists( ObjInfo );

    printf("4.\n"); expose_event( drawing_area, NULL, NULL );
}
static void ChangeStarsMaxMag( GtkWidget *widget, gpointer data) {
    MaxStarMagnitude = gtk_spin_button_get_value( GTK_SPIN_BUTTON(widget) );
    ReLoadStars( );
    printf("5.\n"); expose_event( drawing_area, NULL, NULL );
}

static void ChangeLighting( GtkMenuItem  *menuitem, gpointer data ) {
    int i;
    for (i=0; i< 3; i++){
        if ( gtk_toggle_button_get_active(  GTK_TOGGLE_BUTTON( RadioLightingStyle[i] ) ) ) break;
    }
    LightingStyle = GPOINTER_TO_INT( i );
    printf("6.\n"); expose_event( drawing_area, NULL, NULL );
}

static void ChangeRayleighScattConst( GtkWidget *widget, gpointer data) {
    aInfo->Kr = gtk_spin_button_get_value( GTK_SPIN_BUTTON(widget) );
    aInfo->Kr4PI = aInfo->Kr*4.0*M_PI;
    printf("7.\n"); expose_event( drawing_area, NULL, NULL );
}

static void ChangeRayleighScaleHeight( GtkWidget *widget, gpointer data) {

    free( aInfo->DepthBuf );
    aInfo->rScaleHeight = gtk_spin_button_get_value( GTK_SPIN_BUTTON(widget) );
    MakeOpticalDepthBuffer();
    printf("8.\n"); expose_event( drawing_area, NULL, NULL );
}

static void ChangeMieScattConst( GtkWidget *widget, gpointer data) {
    aInfo->Km = gtk_spin_button_get_value( GTK_SPIN_BUTTON(widget) );
    aInfo->Km4PI = aInfo->Km*4.0*M_PI;
    printf("9.\n"); expose_event( drawing_area, NULL, NULL );
}

static void ChangeMieScaleHeight( GtkWidget *widget, gpointer data) {
    free( aInfo->DepthBuf );
    aInfo->mScaleHeight = gtk_spin_button_get_value( GTK_SPIN_BUTTON(widget) );
    MakeOpticalDepthBuffer();
    printf("10.\n"); expose_event( drawing_area, NULL, NULL );
}

static void ChangeSampleRays( GtkWidget *widget, gpointer data) {
    aInfo->nSampleRays = gtk_spin_button_get_value( GTK_SPIN_BUTTON(widget) );
    printf("11.\n"); expose_event( drawing_area, NULL, NULL );
}

static void ChangeDepthBufSamples( GtkWidget *widget, gpointer data) {
    free( aInfo->DepthBuf );
    aInfo->DepthBuf_nSamples = gtk_spin_button_get_value( GTK_SPIN_BUTTON(widget) );
    MakeOpticalDepthBuffer();
    printf("12.\n"); expose_event( drawing_area, NULL, NULL );
}

static void ChangeSunIntensity( GtkWidget *widget, gpointer data) {
    aInfo->SunIntensity = gtk_spin_button_get_value( GTK_SPIN_BUTTON(widget) );
    printf("13.\n"); expose_event( drawing_area, NULL, NULL );
}

static void ChangeSunLambdaRed( GtkWidget *widget, gpointer data) {
    double  lam2;
    aInfo->Wavelength[0] = 0.001*gtk_spin_button_get_value( GTK_SPIN_BUTTON(widget) );
    lam2 = aInfo->Wavelength[0]*aInfo->Wavelength[0];
    aInfo->Wavelength4[0] = lam2*lam2;
    printf("14.\n"); expose_event( drawing_area, NULL, NULL );
}

static void ChangeSunLambdaGrn( GtkWidget *widget, gpointer data) {
    double  lam2;
    aInfo->Wavelength[1] = 0.001*gtk_spin_button_get_value( GTK_SPIN_BUTTON(widget) );
    lam2 = aInfo->Wavelength[1]*aInfo->Wavelength[1];
    aInfo->Wavelength4[1] = lam2*lam2;
    printf("15.\n"); expose_event( drawing_area, NULL, NULL );
}

static void ChangeSunLambdaBlu( GtkWidget *widget, gpointer data) {
    double  lam2;
    aInfo->Wavelength[2] = 0.001*gtk_spin_button_get_value( GTK_SPIN_BUTTON(widget) );
    lam2 = aInfo->Wavelength[2]*aInfo->Wavelength[2];
    aInfo->Wavelength4[2] = lam2*lam2;
    printf("16.\n"); expose_event( drawing_area, NULL, NULL );
}

static void ChangeMieAsymFactor( GtkWidget *widget, gpointer data) {
    double  u, u2, u3, u4, x;
//    free( aInfo->DepthBuf );
    aInfo->um = gtk_spin_button_get_value( GTK_SPIN_BUTTON(widget) );
    u = aInfo->um;
    u2 = u*u;
    u3 = u2*u;
    u4 = u2*u2;
    x = 5.0/9.0*u + 125.0/729.0*u3
        + sqrt( 64.0/27.0 - 325.0/243.0*u2 + 1250.0/2187.0*u4);
    aInfo->gm = 5.0/9.0*u
        - (4.0/3.0 - 25.0/81.0*u2)*pow(x, -1.0/3.0)
        + pow(x, 1.0/3.0);
    aInfo->gm2 = aInfo->gm*aInfo->gm;
//    MakeOpticalDepthBuffer();
    printf("1.\n"); expose_event( drawing_area, NULL, NULL );
}

static void SelectPitchAngles( GtkMenuItem  *menuitem, gpointer data ) {

    int i, j, AllState, n;

    i = GPOINTER_TO_INT( data );
    if (i<0) {
        /*
         * Select/deselect all of them
         */
        AllState = gtk_check_menu_item_get_active(  GTK_CHECK_MENU_ITEM( PitchAngleCheckMenuItem[0] ) ); // get state of "All" check item
        ShowAllPitchAngles = AllState;
        for (j=0; j<ObjInfo->MagEphemInfo->nAlpha; j++) {
            ShowPitchAngle[j] = AllState;
            g_signal_handler_block( G_OBJECT( PitchAngleCheckMenuItem[j+1] ), PitchAngleHandler[j] );
            gtk_check_menu_item_set_active( GTK_CHECK_MENU_ITEM( PitchAngleCheckMenuItem[j+1] ), AllState );
            g_signal_handler_unblock( G_OBJECT( PitchAngleCheckMenuItem[j+1] ), PitchAngleHandler[j] );
        }

    } else {

        ShowPitchAngle[i] = gtk_check_menu_item_get_active(  GTK_CHECK_MENU_ITEM( PitchAngleCheckMenuItem[i+1] ) );

        /*
         *  Test the individual check boxes. If they are all on,
         *  then set the All item to on as well.
         *
         *  If they are all off, then set the All item to off.
         */
        for (n=0, j=0; j<ObjInfo->MagEphemInfo->nAlpha; j++) n += gtk_check_menu_item_get_active(  GTK_CHECK_MENU_ITEM( PitchAngleCheckMenuItem[j+1] ) );
        if ( n == ObjInfo->MagEphemInfo->nAlpha ){
            ShowAllPitchAngles = TRUE;
            g_signal_handler_block( G_OBJECT( PitchAngleCheckMenuItem[0] ), PitchAngleAllHandler );
            gtk_check_menu_item_set_active( GTK_CHECK_MENU_ITEM( PitchAngleCheckMenuItem[0] ), TRUE );
            g_signal_handler_unblock( G_OBJECT( PitchAngleCheckMenuItem[0] ), PitchAngleAllHandler );

        } else {
            ShowAllPitchAngles = FALSE;
            g_signal_handler_block( G_OBJECT( PitchAngleCheckMenuItem[0] ), PitchAngleAllHandler );
            gtk_check_menu_item_set_active( GTK_CHECK_MENU_ITEM( PitchAngleCheckMenuItem[0] ), FALSE );
            g_signal_handler_unblock( G_OBJECT( PitchAngleCheckMenuItem[0] ), PitchAngleAllHandler );
        }


    }


    printf("17.\n"); expose_event( drawing_area, NULL, NULL );



}

static void ToggleEarth( GtkMenuItem  *menuitem, const GLuint *data ) {
    ShowEarth = !ShowEarth;
    printf("18.\n"); expose_event( drawing_area, NULL, NULL );
}

static void ToggleIridiumLines( GtkWidget  *w, const GLuint *data ) {
    int j, State;
    j = GPOINTER_TO_INT( data );

    State = gtk_toggle_button_get_active(  GTK_TOGGLE_BUTTON( w ) );

    switch ( j ) {

        case 0:
                ShowIridiumSatToGround  =  State;
                break;
        case 1:
                ShowIridiumSatToSun     =  State;
                break;

    }

    expose_event( drawing_area, NULL, NULL );

}

static void ToggleAtmosphere( GtkMenuItem  *menuitem, const GLuint *data ) {
    ShowAtmosphere = !ShowAtmosphere;
    printf("20.\n"); expose_event( drawing_area, NULL, NULL );
}

static void ToggleStars( GtkMenuItem  *menuitem, const GLuint *data ) {
    ShowStars = !ShowStars;
    printf("21.\n"); expose_event( drawing_area, NULL, NULL );
}

static void ToggleFieldLineMode( GtkMenuItem  *menuitem, const GLuint *data ) {
    ShowFullFieldLine = !ShowFullFieldLine;
    printf("22.\n"); expose_event( drawing_area, NULL, NULL );
}

static void ToggleDriftShellSurface( GtkMenuItem  *menuitem, const GLuint *data ) {
    ShowDriftShellSurface = !ShowDriftShellSurface;
    printf("23.\n"); expose_event( drawing_area, NULL, NULL );
}



static void change_material (GtkMenuItem  *menuitem, MaterialProp *mat) {
    mat_current = mat;
    printf("24.\n"); expose_event( drawing_area, NULL, NULL );
}





/* For popup menu. */
static gboolean button_press_event_popup_menu (GtkWidget *widget, GdkEventButton *event, gpointer data) {

    if (event->button == 3) {
        /* Popup menu. */
        gtk_menu_popup (GTK_MENU (widget), NULL, NULL, NULL, NULL, event->button, event->time);
        return TRUE;
    }

    return FALSE;

}






/*
 * Creates the popup menu.
 */
static GtkWidget * create_popup_menu (GtkWidget *drawing_area) {



    GtkWidget *PitchAngle_menu;
    GtkWidget *tearoff;
    GtkWidget *materials_menu;
    GtkWidget *options_menu;
    GtkWidget *lighting_menu;
    GtkWidget *menu;
    GtkWidget *menu_item;
    int       i, n;
    char      Str[256];


    /*
     * Pitch Angle submenu.
     */
    PitchAngle_menu = gtk_menu_new( );

    tearoff = gtk_tearoff_menu_item_new( );
    gtk_menu_shell_append( GTK_MENU_SHELL( PitchAngle_menu ), tearoff );
    gtk_widget_show( tearoff );


    n = 0;
    sprintf(Str, "Show All Pitch Angles" );
    PitchAngleCheckMenuItem[n] = gtk_check_menu_item_new_with_label( Str );
    gtk_menu_shell_append( GTK_MENU_SHELL( PitchAngle_menu ), PitchAngleCheckMenuItem[n] );
    if (ShowAllPitchAngles) gtk_check_menu_item_set_active( GTK_CHECK_MENU_ITEM( PitchAngleCheckMenuItem[n] ), TRUE );
    PitchAngleAllHandler = g_signal_connect( G_OBJECT( PitchAngleCheckMenuItem[n] ), "activate", G_CALLBACK( SelectPitchAngles ), GINT_TO_POINTER(-1) );
    gtk_widget_show( PitchAngleCheckMenuItem[n] );
    ++n;

    for (i=0; i<ObjInfo->MagEphemInfo->nAlpha; i++) {
        sprintf(Str, "Show %g Pitch Angle", ObjInfo->MagEphemInfo->Alpha[i] );
        PitchAngleCheckMenuItem[n] = gtk_check_menu_item_new_with_label( Str );
        gtk_menu_shell_append( GTK_MENU_SHELL( PitchAngle_menu ), PitchAngleCheckMenuItem[n] );
        if (ShowPitchAngle[i]) gtk_check_menu_item_set_active( GTK_CHECK_MENU_ITEM( PitchAngleCheckMenuItem[n] ), TRUE );
        PitchAngleHandler[i] = g_signal_connect( G_OBJECT( PitchAngleCheckMenuItem[n] ), "activate", G_CALLBACK( SelectPitchAngles ), GINT_TO_POINTER(i) );
        gtk_widget_show( PitchAngleCheckMenuItem[n] );
        ++n;
    }






  /*
   * Materials submenu.
   */
  materials_menu = gtk_menu_new ();

  /* Emerald */
  menu_item = gtk_menu_item_new_with_label ("Emerald");
  gtk_menu_shell_append (GTK_MENU_SHELL (materials_menu), menu_item);
  g_signal_connect (G_OBJECT (menu_item), "activate", G_CALLBACK (change_material), &mat_emerald);
  gtk_widget_show (menu_item);

  /* Jade */
  menu_item = gtk_menu_item_new_with_label ("Jade");
  gtk_menu_shell_append (GTK_MENU_SHELL (materials_menu), menu_item);
  g_signal_connect (G_OBJECT (menu_item), "activate", G_CALLBACK (change_material), &mat_jade);
  gtk_widget_show (menu_item);

  /* Obsidian */
  menu_item = gtk_menu_item_new_with_label ("Obsidian");
  gtk_menu_shell_append (GTK_MENU_SHELL (materials_menu), menu_item);
  g_signal_connect (G_OBJECT (menu_item), "activate", G_CALLBACK (change_material), &mat_obsidian);
  gtk_widget_show (menu_item);

  /* Pearl */
  menu_item = gtk_menu_item_new_with_label ("Pearl");
  gtk_menu_shell_append (GTK_MENU_SHELL (materials_menu), menu_item);
  g_signal_connect (G_OBJECT (menu_item), "activate", G_CALLBACK (change_material), &mat_pearl);
  gtk_widget_show (menu_item);

  /* Ruby */
  menu_item = gtk_menu_item_new_with_label ("Ruby");
  gtk_menu_shell_append (GTK_MENU_SHELL (materials_menu), menu_item);
  g_signal_connect (G_OBJECT (menu_item), "activate", G_CALLBACK (change_material), &mat_ruby);
  gtk_widget_show (menu_item);

  /* Turquoise */
  menu_item = gtk_menu_item_new_with_label ("Turquoise");
  gtk_menu_shell_append (GTK_MENU_SHELL (materials_menu), menu_item);
  g_signal_connect (G_OBJECT (menu_item), "activate", G_CALLBACK (change_material), &mat_turquoise);
  gtk_widget_show (menu_item);

  /* Brass */
  menu_item = gtk_menu_item_new_with_label ("Brass");
  gtk_menu_shell_append (GTK_MENU_SHELL (materials_menu), menu_item);
  g_signal_connect (G_OBJECT (menu_item), "activate", G_CALLBACK (change_material), &mat_brass);
  gtk_widget_show (menu_item);

  /* Bronze */
  menu_item = gtk_menu_item_new_with_label ("Bronze");
  gtk_menu_shell_append (GTK_MENU_SHELL (materials_menu), menu_item);
  g_signal_connect (G_OBJECT (menu_item), "activate", G_CALLBACK (change_material), &mat_bronze);
  gtk_widget_show (menu_item);

  /* Chrome */
  menu_item = gtk_menu_item_new_with_label ("Chrome");
  gtk_menu_shell_append (GTK_MENU_SHELL (materials_menu), menu_item);
  g_signal_connect (G_OBJECT (menu_item), "activate", G_CALLBACK (change_material), &mat_chrome);
  gtk_widget_show (menu_item);

  /* Copper */
  menu_item = gtk_menu_item_new_with_label ("Copper");
  gtk_menu_shell_append (GTK_MENU_SHELL (materials_menu), menu_item);
  g_signal_connect (G_OBJECT (menu_item), "activate", G_CALLBACK (change_material), &mat_copper);
  gtk_widget_show (menu_item);

  /* Gold */
  menu_item = gtk_menu_item_new_with_label ("Gold");
  gtk_menu_shell_append (GTK_MENU_SHELL (materials_menu), menu_item);
  g_signal_connect (G_OBJECT (menu_item), "activate", G_CALLBACK (change_material), &mat_gold);
  gtk_widget_show (menu_item);

  /* Silver */
  menu_item = gtk_menu_item_new_with_label ("Silver");
  gtk_menu_shell_append (GTK_MENU_SHELL (materials_menu), menu_item);
  g_signal_connect (G_OBJECT (menu_item), "activate", G_CALLBACK (change_material), &mat_silver);
  gtk_widget_show (menu_item);




  /*
   * Options submenu.
   */
  options_menu = gtk_menu_new( );

  menu_item = gtk_menu_item_new_with_label( "Toggle Field Lines" );
  gtk_menu_shell_append( GTK_MENU_SHELL( options_menu ), menu_item );
  g_signal_connect( G_OBJECT( menu_item ), "activate", G_CALLBACK( ToggleFieldLineMode ), NULL );
  gtk_widget_show( menu_item );

  menu_item = gtk_menu_item_new_with_label( "Toggle Drift Shell" );
  gtk_menu_shell_append( GTK_MENU_SHELL( options_menu ), menu_item );
  g_signal_connect( G_OBJECT( menu_item ), "activate", G_CALLBACK( ToggleDriftShellSurface ), NULL );
  gtk_widget_show( menu_item );

  menu_item = gtk_menu_item_new_with_label( "Toggle Earth" );
  gtk_menu_shell_append( GTK_MENU_SHELL( options_menu ), menu_item );
  g_signal_connect( G_OBJECT( menu_item ), "activate", G_CALLBACK( ToggleEarth ), NULL );
  gtk_widget_show( menu_item );

  menu_item = gtk_menu_item_new_with_label( "Toggle Atmosphere" );
  gtk_menu_shell_append( GTK_MENU_SHELL( options_menu ), menu_item );
  g_signal_connect( G_OBJECT( menu_item ), "activate", G_CALLBACK( ToggleAtmosphere ), NULL );
  gtk_widget_show( menu_item );


  /*
   * Lighting submenu.
   */
  lighting_menu = gtk_menu_new( );

  menu_item = gtk_menu_item_new_with_label( "Solar Illumination" );
  gtk_menu_shell_append( GTK_MENU_SHELL( lighting_menu ), menu_item );
  g_signal_connect( G_OBJECT( menu_item ), "activate", G_CALLBACK( ChangeLighting ), GINT_TO_POINTER( 0 ) );
  gtk_widget_show( menu_item );

  menu_item = gtk_menu_item_new_with_label( "Fixed Illumination" );
  gtk_menu_shell_append( GTK_MENU_SHELL( lighting_menu ), menu_item );
  g_signal_connect( G_OBJECT( menu_item ), "activate", G_CALLBACK( ChangeLighting ), GINT_TO_POINTER( 1 ) );
  gtk_widget_show( menu_item );

  menu_item = gtk_menu_item_new_with_label( "Vertex Shader Test" );
  gtk_menu_shell_append( GTK_MENU_SHELL( lighting_menu ), menu_item );
  g_signal_connect( G_OBJECT( menu_item ), "activate", G_CALLBACK( ChangeLighting ), GINT_TO_POINTER( 2 ) );
  gtk_widget_show( menu_item );







  /*
   * Root popup menu.
   */
  menu = gtk_menu_new ();

  /* Pitch Angle */
  menu_item = gtk_menu_item_new_with_label( "Pitch Angles" );
  gtk_menu_item_set_submenu( GTK_MENU_ITEM( menu_item ), PitchAngle_menu );
  gtk_menu_shell_append( GTK_MENU_SHELL( menu ), menu_item );
  gtk_widget_show( menu_item );

  /* Options */
  menu_item = gtk_menu_item_new_with_label( "Options" );
  gtk_menu_item_set_submenu( GTK_MENU_ITEM( menu_item ), options_menu );
  gtk_menu_shell_append( GTK_MENU_SHELL( menu ), menu_item );
  gtk_widget_show( menu_item );

  /* Lighting */
  menu_item = gtk_menu_item_new_with_label( "Lighting" );
  gtk_menu_item_set_submenu( GTK_MENU_ITEM( menu_item ), lighting_menu );
  gtk_menu_shell_append( GTK_MENU_SHELL( menu ), menu_item );
  gtk_widget_show( menu_item );



  /* Materials */
  menu_item = gtk_menu_item_new_with_label ("Materials");
  gtk_menu_item_set_submenu (GTK_MENU_ITEM (menu_item), materials_menu);
  gtk_menu_shell_append (GTK_MENU_SHELL (menu), menu_item);
  gtk_widget_show (menu_item);

  /* Quit */
  menu_item = gtk_menu_item_new_with_label ("Quit");
  gtk_menu_shell_append (GTK_MENU_SHELL (menu), menu_item);
  g_signal_connect (G_OBJECT (menu_item), "activate", G_CALLBACK (quit_app), NULL);
  gtk_widget_show (menu_item);
	
  return menu;
}




static void print_gl_config_attrib (GdkGLConfig *glconfig, const gchar *attrib_str, int attrib, gboolean is_boolean) {

    int value;

    g_print ("%s = ", attrib_str);
    if (gdk_gl_config_get_attrib (glconfig, attrib, &value)) {
        if (is_boolean){
            g_print ("%s\n", value == TRUE ? "TRUE" : "FALSE");
        } else {
            g_print ("%d\n", value);
        }
    } else {
        g_print ("*** Cannot get %s attribute value\n", attrib_str);
    }

}




static void examine_gl_config_attrib (GdkGLConfig *glconfig) {

    g_print ("\nOpenGL visual configurations :\n\n");
    g_print ("gdk_gl_config_is_rgba (glconfig) = %s\n", gdk_gl_config_is_rgba (glconfig) ? "TRUE" : "FALSE");
    g_print ("gdk_gl_config_is_double_buffered (glconfig) = %s\n", gdk_gl_config_is_double_buffered (glconfig) ? "TRUE" : "FALSE");
    g_print ("gdk_gl_config_is_stereo (glconfig) = %s\n", gdk_gl_config_is_stereo (glconfig) ? "TRUE" : "FALSE");
    g_print ("gdk_gl_config_has_alpha (glconfig) = %s\n", gdk_gl_config_has_alpha (glconfig) ? "TRUE" : "FALSE");
    g_print ("gdk_gl_config_has_depth_buffer (glconfig) = %s\n", gdk_gl_config_has_depth_buffer (glconfig) ? "TRUE" : "FALSE");
    g_print ("gdk_gl_config_has_stencil_buffer (glconfig) = %s\n", gdk_gl_config_has_stencil_buffer (glconfig) ? "TRUE" : "FALSE");
    g_print ("gdk_gl_config_has_accum_buffer (glconfig) = %s\n", gdk_gl_config_has_accum_buffer (glconfig) ? "TRUE" : "FALSE");
    g_print ("\n");

    print_gl_config_attrib (glconfig, "GDK_GL_USE_GL",           GDK_GL_USE_GL,           TRUE);
    print_gl_config_attrib (glconfig, "GDK_GL_BUFFER_SIZE",      GDK_GL_BUFFER_SIZE,      FALSE);
    print_gl_config_attrib (glconfig, "GDK_GL_LEVEL",            GDK_GL_LEVEL,            FALSE);
    print_gl_config_attrib (glconfig, "GDK_GL_RGBA",             GDK_GL_RGBA,             TRUE);
    print_gl_config_attrib (glconfig, "GDK_GL_DOUBLEBUFFER",     GDK_GL_DOUBLEBUFFER,     TRUE);
    print_gl_config_attrib (glconfig, "GDK_GL_STEREO",           GDK_GL_STEREO,           TRUE);
    print_gl_config_attrib (glconfig, "GDK_GL_AUX_BUFFERS",      GDK_GL_AUX_BUFFERS,      FALSE);
    print_gl_config_attrib (glconfig, "GDK_GL_RED_SIZE",         GDK_GL_RED_SIZE,         FALSE);
    print_gl_config_attrib (glconfig, "GDK_GL_GREEN_SIZE",       GDK_GL_GREEN_SIZE,       FALSE);
    print_gl_config_attrib (glconfig, "GDK_GL_BLUE_SIZE",        GDK_GL_BLUE_SIZE,        FALSE);
    print_gl_config_attrib (glconfig, "GDK_GL_ALPHA_SIZE",       GDK_GL_ALPHA_SIZE,       FALSE);
    print_gl_config_attrib (glconfig, "GDK_GL_DEPTH_SIZE",       GDK_GL_DEPTH_SIZE,       FALSE);
    print_gl_config_attrib (glconfig, "GDK_GL_STENCIL_SIZE",     GDK_GL_STENCIL_SIZE,     FALSE);
    print_gl_config_attrib (glconfig, "GDK_GL_ACCUM_RED_SIZE",   GDK_GL_ACCUM_RED_SIZE,   FALSE);
    print_gl_config_attrib (glconfig, "GDK_GL_ACCUM_GREEN_SIZE", GDK_GL_ACCUM_GREEN_SIZE, FALSE);
    print_gl_config_attrib (glconfig, "GDK_GL_ACCUM_BLUE_SIZE",  GDK_GL_ACCUM_BLUE_SIZE,  FALSE);
    print_gl_config_attrib (glconfig, "GDK_GL_ACCUM_ALPHA_SIZE", GDK_GL_ACCUM_ALPHA_SIZE, FALSE);
    g_print ("\n");

}









void create_ViewDriftShell( void *data ) {

    GdkGLConfig         *glconfig;
    gint                major, minor;
    GtkWidget           *vbox;
    GtkWidget           *menu;
    GtkWidget           *Menubar;
    GtkAccelGroup       *AccelGroup;
    GtkItemFactory      *ItemFactory;
    int                 nMenuItems;

    int                 i;


    /*
     * Initialize GtkGLExt.
     */
//    gtk_gl_init( &xInfo->argc, &xInfo->argv );





    /*
     * Query OpenGL extension versioy.
     */
    gdk_gl_query_version (&major, &minor);
    g_print ("\nOpenGL extension version - %d.%d\n", major, minor);




    /*
     * Configure OpenGL-capable visual.
     */
    // Try double-buffered visual
    glconfig = gdk_gl_config_new_by_mode( GDK_GL_MODE_RGB | GDK_GL_MODE_DEPTH | GDK_GL_MODE_DOUBLE | GDK_GL_MODE_ALPHA );
    if ( glconfig == NULL ) {
        g_print( "*** Cannot find the double-buffered visual.\n" );
        g_print( "*** Trying single-buffered visual.\n" );

        // Try single-buffered visual
        glconfig = gdk_gl_config_new_by_mode( GDK_GL_MODE_RGB | GDK_GL_MODE_DEPTH );
        if ( glconfig == NULL ) {
            g_print( "*** No appropriate OpenGL-capable visual found.\n" );
            exit( 1 );
        }
    }
    examine_gl_config_attrib( glconfig );






    /*
     * Top-level window.
     */
    ViewDriftShellWindow = gtk_window_new( GTK_WINDOW_TOPLEVEL );
    gtk_window_set_title( GTK_WINDOW (ViewDriftShellWindow), "Drift Shell" );
    g_signal_connect( G_OBJECT(ViewDriftShellWindow), "delete_event", G_CALLBACK(quit_app), NULL );



    /*
     * Get automatically redrawn if any of their children changed allocation.
     */
    gtk_container_set_reallocate_redraws( GTK_CONTAINER(ViewDriftShellWindow), TRUE );




    /*
     * VBox.
     */
    vbox = gtk_vbox_new( FALSE, 0 );
    gtk_container_add( GTK_CONTAINER(ViewDriftShellWindow), vbox );
    gtk_widget_show( vbox );



   /*
     *  Add a menubar to the vbox
     */
    AccelGroup  = gtk_accel_group_new();
    ItemFactory = gtk_item_factory_new( GTK_TYPE_MENU_BAR, "<main>", AccelGroup);
    nMenuItems  = sizeof( MenuItems ) / sizeof( MenuItems[0] );
    gtk_item_factory_create_items( ItemFactory, nMenuItems, MenuItems, NULL );
    Menubar = gtk_item_factory_get_widget( ItemFactory, "<main>" );
    gtk_box_pack_start( GTK_BOX(vbox) , Menubar, FALSE, FALSE, 0);
    gtk_window_add_accel_group( GTK_WINDOW( ViewDriftShellWindow ), AccelGroup);
    gtk_widget_show( Menubar );




    /*
     * Drawing area for drawing OpenGL scene.
     */
    drawing_area = gtk_drawing_area_new( );
    gtk_widget_set_size_request( drawing_area, 500, 500 );



    /*
     * Set OpenGL-capability to the widget.
     */
    gtk_widget_set_gl_capability( drawing_area, glconfig, NULL, TRUE, GDK_GL_RGBA_TYPE );

    // should not have to include GDK_POINTER_MOTION_MASK here but openSUSE 11.2 (or KDE 4.3) doesnt seem to pass me the
    // GDK_BUTTON1_MOTION_MASK type masks...
    gtk_widget_add_events( drawing_area, GDK_POINTER_MOTION_MASK | GDK_BUTTON1_MOTION_MASK | GDK_BUTTON2_MOTION_MASK | GDK_BUTTON3_MOTION_MASK
                                         | GDK_BUTTON_PRESS_MASK | GDK_BUTTON_RELEASE_MASK | GDK_VISIBILITY_NOTIFY_MASK );

    g_signal_connect_after(     G_OBJECT( drawing_area ),           "realize",                  G_CALLBACK( realize ),                   NULL );
    g_signal_connect(           G_OBJECT( drawing_area ),           "configure_event",          G_CALLBACK( configure_event ),           NULL );
    g_signal_connect(           G_OBJECT( drawing_area ),           "expose_event",             G_CALLBACK( expose_event ),              NULL );
    g_signal_connect(           G_OBJECT( drawing_area ),           "button_press_event",       G_CALLBACK( button_press_event ),        NULL );
    g_signal_connect(           G_OBJECT( drawing_area ),           "button_release_event",     G_CALLBACK( button_release_event ),      NULL );
    g_signal_connect(           G_OBJECT( drawing_area ),           "motion_notify_event",      G_CALLBACK( motion_notify_event ),       NULL );
    g_signal_connect(           G_OBJECT( drawing_area ),           "map_event",                G_CALLBACK( map_event ),                 NULL );
    g_signal_connect(           G_OBJECT( drawing_area ),           "unmap_event",              G_CALLBACK( unmap_event ),               NULL );
    g_signal_connect(           G_OBJECT( drawing_area ),           "visibility_notify_event",  G_CALLBACK( visibility_notify_event ),   NULL );
    g_signal_connect_swapped(   G_OBJECT( ViewDriftShellWindow ),   "key_press_event",          G_CALLBACK( key_press_event ),           drawing_area );

    gtk_box_pack_start( GTK_BOX(vbox), drawing_area, TRUE, TRUE, 0 );
    gtk_widget_show( drawing_area );

printf("1. HERE\n");

    /*
     *  Set Pitch Angles to Show
     */
    ShowAllPitchAngles = FALSE;
    for (i=0; i<=ObjInfo->MagEphemInfo->nAlpha; i++){
        //ShowPitchAngle[i] = TRUE;
        ShowPitchAngle[i] = FALSE;
        ShowPitchAngle2[i] = FALSE;
    }



    /*
     * Popup menu.
     */
    menu = create_popup_menu( drawing_area );

    /*
     * Signal handler
     */
//    g_signal_connect_swapped( G_OBJECT(drawing_area), "button_press_event", G_CALLBACK(button_press_event_popup_menu), menu );



    /*
     * Simple quit button.
     */
//    button = gtk_button_new_with_label( "Quit" );
//    g_signal_connect_swapped( GTK_OBJECT( button ), "clicked", G_CALLBACK( quit_app ), NULL );
//    gtk_box_pack_start( GTK_BOX(vbox), button, FALSE, FALSE, 0 );
//    gtk_widget_show( button );



//    aInfo = New_aInfo();
//    InitAtmosphere();




printf("2. HERE ViewDriftShellWindow = %p\n", ViewDriftShellWindow);
    gtk_widget_show( ViewDriftShellWindow );


    return;

}



static void GetOneMaterialColor( GLfloat *c, GtkWidget *b1 ){

    GdkColor    Color;
    guint16     Alpha;

    /*
     *  Get diffuse color
     */
    gtk_color_button_get_color( GTK_COLOR_BUTTON(b1), &Color );
    Alpha = gtk_color_button_get_alpha( GTK_COLOR_BUTTON(b1) );
    c[0] = Color.red/65535.0;
    c[1] = Color.green/65535.0;
    c[2] = Color.blue/65535.0;
    c[3] = Alpha/65535.0;

}


static void GetAllMaterialColors( MaterialProp *material, GtkWidget *b1, GtkWidget *b2, GtkWidget *b3, GtkWidget *s1 ){

    GdkColor    Color;
    guint16     Alpha;


    /*
     *  Get diffuse color
     */
    gtk_color_button_get_color( GTK_COLOR_BUTTON(b1), &Color );
    Alpha = gtk_color_button_get_alpha( GTK_COLOR_BUTTON(b1) );
    material->diffuse[0] = Color.red/65535.0;
    material->diffuse[1] = Color.green/65535.0;
    material->diffuse[2] = Color.blue/65535.0;
    material->diffuse[3] = Alpha/65535.0;

    /*
     *  Get ambient color
     */
    gtk_color_button_get_color( GTK_COLOR_BUTTON(b2), &Color );
    Alpha = gtk_color_button_get_alpha( GTK_COLOR_BUTTON(b2) );
    material->ambient[0] = Color.red/65535.0;
    material->ambient[1] = Color.green/65535.0;
    material->ambient[2] = Color.blue/65535.0;
    material->ambient[3] = Alpha/65535.0;

    /*
     *  Get specular color
     */
    gtk_color_button_get_color( GTK_COLOR_BUTTON(b3), &Color );
    Alpha = gtk_color_button_get_alpha( GTK_COLOR_BUTTON(b3) );
    material->specular[0] = Color.red/65535.0;
    material->specular[1] = Color.green/65535.0;
    material->specular[2] = Color.blue/65535.0;
    material->specular[3] = Alpha/65535.0;


    /*
     *  Get specular color
     */
    material->shininess = (GLfloat)gtk_spin_button_get_value( GTK_SPIN_BUTTON(s1) );


}

static void ChangeMaterialShininess( GtkMenuItem  *menuitem, gpointer data ) {

    int             k, kk, indx;
    MaterialProp    *material;

    k = GPOINTER_TO_INT( data );

    if ( k< 100 ) {
        gInfo->FieldLineMaterial[k].shininess = (GLfloat)gtk_spin_button_get_value( GTK_SPIN_BUTTON(gInfo->FieldLineShininessButton[k]) );
        indx = gtk_combo_box_get_active( GTK_COMBO_BOX(gInfo->FieldLineMaterialButton[k]) );
        if (indx >= 0) {
            /*
             * If we are here, it means we are adjusting the spin button while there is a named material selected.
             * So we need to set it back to "Custom" and re-read the color buttons too...
             */
            material = &gInfo->FieldLineMaterial[k];
//            GetAllMaterialColors( material, gInfo->FieldLineDiffuseColorButton[k], gInfo->FieldLineAmbientColorButton[k],
//                                             gInfo->FieldLineSpecularColorButton[k], gInfo->FieldLineShininessButton[k]);

            g_signal_handler_block( G_OBJECT( gInfo->FieldLineMaterialButton[k] ), gInfo->FieldLineMaterialButtonHandler[k] );
            gtk_combo_box_set_active( GTK_COMBO_BOX(gInfo->FieldLineMaterialButton[k]), 0);
            g_signal_handler_unblock( G_OBJECT( gInfo->FieldLineMaterialButton[k] ), gInfo->FieldLineMaterialButtonHandler[k] );
        }

    } else {
        kk = k-100;
        gInfo->DriftShellMaterial[kk].shininess = (GLfloat)gtk_spin_button_get_value( GTK_SPIN_BUTTON(gInfo->DriftShellShininessButton[kk]) );
        indx = gtk_combo_box_get_active( GTK_COMBO_BOX(gInfo->DriftShellMaterialButton[kk]) );
        if (indx >= 0) {
            /*
             * If we are here, it means we are adjusting the spin button while there is a named material selected.
             * So we need to set it back to "Custom" and re-read the color buttons too...
             */
            material = &gInfo->DriftShellMaterial[kk];
//            GetAllMaterialColors( material, gInfo->DriftShellDiffuseColorButton[kk], gInfo->DriftShellAmbientColorButton[kk],
//                                             gInfo->DriftShellSpecularColorButton[kk], gInfo->DriftShellShininessButton[kk]);

            g_signal_handler_block( G_OBJECT( gInfo->DriftShellMaterialButton[kk] ), gInfo->DriftShellMaterialButtonHandler[kk] );
            gtk_combo_box_set_active( GTK_COMBO_BOX(gInfo->DriftShellMaterialButton[kk]), 0);
            g_signal_handler_unblock( G_OBJECT( gInfo->DriftShellMaterialButton[kk] ), gInfo->DriftShellMaterialButtonHandler[kk] );
        }
    }




    printf("25.\n"); expose_event( drawing_area, NULL, NULL );

}


static void ChangeMaterial( GtkWidget  *widget, gpointer data ) {
    int             i, k, Flag, indx;
    MaterialProp    *material;
    GtkWidget       *button;
    GdkColor        color;

    k = GPOINTER_TO_INT( data );
    Flag = 0;
    if ( k >= 100 ){
        k -= 100;
        Flag = 1;
    }

    if (Flag){
        button   = gInfo->DriftShellMaterialButton[k];
        material = &gInfo->DriftShellMaterial[k];
    } else {
        button   = gInfo->FieldLineMaterialButton[k];
        material = &gInfo->FieldLineMaterial[k];
    }


    /*
     * Get material
     */
    indx = gtk_combo_box_get_active( GTK_COMBO_BOX(button) );
    --indx;
printf("indx = %d\n", indx);

    if ( indx == 0 ) {
        // If "Custom:"
        if ( Flag ){
//            gtk_color_button_get_color( GTK_COLOR_BUTTON(gInfo->DriftShellDiffuseColorButton[k]), &material->diffuse );
//            gtk_color_button_get_color( GTK_COLOR_BUTTON(gInfo->DriftShellAmbientColorButton[k]), &material->ambient );
//            gtk_color_button_get_color( GTK_COLOR_BUTTON(gInfo->DriftShellSpecularColorButton[k]), &material->specular );
        } else {
//            gtk_color_button_get_color( GTK_COLOR_BUTTON(gInfo->FieldLineDiffuseColorButton[k]), &material->diffuse );
//            gtk_color_button_get_color( GTK_COLOR_BUTTON(gInfo->FieldLineAmbientColorButton[k]), &material->ambient );
//            gtk_color_button_get_color( GTK_COLOR_BUTTON(gInfo->FieldLineSpecularColorButton[k]), &material->specular );
        }
    } else if ( indx < nNamedMaterials ) {

        for (i=0; i<4; i++) material->ambient[i]  = NamedMaterials[indx].Ambient[i];
        for (i=0; i<4; i++) material->diffuse[i]  = NamedMaterials[indx].Diffuse[i];
        for (i=0; i<4; i++) material->specular[i] = NamedMaterials[indx].Specular[i];
        material->shininess = NamedMaterials[indx].Shininess;
        if ( Flag ){

            color.red = material->diffuse[0]*65535; color.green = material->diffuse[1]*65535; color.blue = material->diffuse[2]*65535;
            gtk_color_button_set_color( GTK_COLOR_BUTTON(gInfo->DriftShellDiffuseColorButton[k]), &color );
            gtk_color_button_set_alpha( GTK_COLOR_BUTTON(gInfo->DriftShellDiffuseColorButton[k]), material->diffuse[3]*65535 );

            color.red = material->ambient[0]*65535; color.green = material->ambient[1]*65535; color.blue = material->ambient[2]*65535;
            gtk_color_button_set_color( GTK_COLOR_BUTTON(gInfo->DriftShellAmbientColorButton[k]), &color );
            gtk_color_button_set_alpha( GTK_COLOR_BUTTON(gInfo->DriftShellAmbientColorButton[k]), material->ambient[3]*65535 );

            color.red = material->specular[0]*65535; color.green = material->specular[1]*65535; color.blue = material->specular[2]*65535;
            gtk_color_button_set_color( GTK_COLOR_BUTTON(gInfo->DriftShellSpecularColorButton[k]), &color );
            gtk_color_button_set_alpha( GTK_COLOR_BUTTON(gInfo->DriftShellSpecularColorButton[k]), material->specular[3]*65535 );

            g_signal_handler_block( G_OBJECT( gInfo->DriftShellShininessButton[k] ), gInfo->DriftShellShininessButtonHandler[k] );
            gtk_spin_button_set_value( GTK_SPIN_BUTTON(gInfo->DriftShellShininessButton[k]), material->shininess );
            g_signal_handler_unblock( G_OBJECT( gInfo->DriftShellShininessButton[k] ), gInfo->DriftShellShininessButtonHandler[k] );

        } else {

            color.red = material->diffuse[0]*65535; color.green = material->diffuse[1]*65535; color.blue = material->diffuse[2]*65535;
            gtk_color_button_set_color( GTK_COLOR_BUTTON(gInfo->FieldLineDiffuseColorButton[k]), &color );
            gtk_color_button_set_alpha( GTK_COLOR_BUTTON(gInfo->FieldLineDiffuseColorButton[k]), material->diffuse[3]*65535 );

            color.red = material->ambient[0]*65535; color.green = material->ambient[1]*65535; color.blue = material->ambient[2]*65535;
            gtk_color_button_set_color( GTK_COLOR_BUTTON(gInfo->FieldLineAmbientColorButton[k]), &color );
            gtk_color_button_set_alpha( GTK_COLOR_BUTTON(gInfo->FieldLineAmbientColorButton[k]), material->ambient[3]*65535 );

            color.red = material->specular[0]*65535; color.green = material->specular[1]*65535; color.blue = material->specular[2]*65535;
            gtk_color_button_set_color( GTK_COLOR_BUTTON(gInfo->FieldLineSpecularColorButton[k]), &color );
            gtk_color_button_set_alpha( GTK_COLOR_BUTTON(gInfo->FieldLineSpecularColorButton[k]), material->specular[3]*65535 );

            g_signal_handler_block( G_OBJECT( gInfo->FieldLineShininessButton[k] ), gInfo->FieldLineShininessButtonHandler[k] );
            gtk_spin_button_set_value( GTK_SPIN_BUTTON(gInfo->FieldLineShininessButton[k]), material->shininess );
            g_signal_handler_unblock( G_OBJECT( gInfo->FieldLineShininessButton[k] ), gInfo->FieldLineShininessButtonHandler[k] );

        }

    } else {

        /*
         * Get colors and shininess
         */
        if ( Flag) {
            GetAllMaterialColors( material, gInfo->DriftShellDiffuseColorButton[k], gInfo->DriftShellAmbientColorButton[k],
                                             gInfo->DriftShellSpecularColorButton[k], gInfo->DriftShellShininessButton[k]);
        } else {
            GetAllMaterialColors( material, gInfo->FieldLineDiffuseColorButton[k], gInfo->FieldLineAmbientColorButton[k],
                                             gInfo->FieldLineSpecularColorButton[k], gInfo->FieldLineShininessButton[k]);
        }

    }

    printf("26.\n"); expose_event( drawing_area, NULL, NULL );


}



static void ChangeMaterialColor( GtkMenuItem  *menuitem, gpointer data ) {

    int             i, j, k, Flag, indx;
    GtkWidget       *button;
    MaterialProp    *material;
    guint16         Alpha;

    /*
     * k is decoded to give Flag, i and j.
     * Flag determines if we are adjusting FieldLine or DriftShell buttons.
     * i is the pitch angle row.
     * j is the button type (diffuse, ambient, specular)
     */
    k = GPOINTER_TO_INT( data );
    Flag = 0;
    if ( k >= 100 ){
        k -= 100;
        Flag = 1;
    }
    i = k/3;
    j = k-i*3;




    if (Flag){

        material = &gInfo->DriftShellMaterial[i];

        indx = gtk_combo_box_get_active( GTK_COMBO_BOX(gInfo->DriftShellMaterialButton[i]) );
printf("Flag = %d   indx = %d\n", Flag, indx );
        if (indx >= 0) {
            /*
             *  If we are here, it means we are adjusting a color button while
             *  there is a named material selected.  So we need to set it back
             *  to "Custom".
             */
            g_signal_handler_block( G_OBJECT( gInfo->DriftShellMaterialButton[i] ), gInfo->DriftShellMaterialButtonHandler[i] );
            gtk_combo_box_set_active( GTK_COMBO_BOX(gInfo->DriftShellMaterialButton[i]), 0);
            g_signal_handler_unblock( G_OBJECT( gInfo->DriftShellMaterialButton[i] ), gInfo->DriftShellMaterialButtonHandler[i] );
        }

        if (j==0) {
            button = gInfo->DriftShellDiffuseColorButton[i];
            GetOneMaterialColor( material->diffuse, button );
        } else if (j==1) {
            GetOneMaterialColor( material->ambient, button );
            button = gInfo->DriftShellAmbientColorButton[i];
        } else {
            GetOneMaterialColor( material->specular, button );
            button = gInfo->DriftShellSpecularColorButton[i];
        }
        Alpha = gtk_color_button_get_alpha( GTK_COLOR_BUTTON(button) );

    } else {

        material = &gInfo->FieldLineMaterial[i];
        indx = gtk_combo_box_get_active( GTK_COMBO_BOX(gInfo->FieldLineMaterialButton[i]) );
        if (indx >= 0) {
            /*
             * If we are here, it means we are adjusting a color button while there is a named material selected.
             * So we need to set it back to "Custom".
             */
            g_signal_handler_block( G_OBJECT( gInfo->FieldLineMaterialButton[i] ), gInfo->FieldLineMaterialButtonHandler[i] );
            gtk_combo_box_set_active( GTK_COMBO_BOX(gInfo->FieldLineMaterialButton[i]), 0);
            g_signal_handler_unblock( G_OBJECT( gInfo->FieldLineMaterialButton[i] ), gInfo->FieldLineMaterialButtonHandler[i] );
        }

        if (j==0) {
            button = gInfo->FieldLineDiffuseColorButton[i];
            GetOneMaterialColor( material->diffuse, button );
        } else if (j==1) {
            button = gInfo->FieldLineAmbientColorButton[i];
            GetOneMaterialColor( material->ambient, button );
        } else           {
            button = gInfo->FieldLineSpecularColorButton[i];
            GetOneMaterialColor( material->specular, button );
        }
        Alpha = gtk_color_button_get_alpha( GTK_COLOR_BUTTON(button));
    }

    printf("27.\n"); expose_event( drawing_area, NULL, NULL );

}

static void SelectPitchAngles2( GtkWidget  *widget, gpointer data ) {

    int i, j, AllState, AllState2, n;

    i = GPOINTER_TO_INT( data );
    if (i<100){
        if (i==ObjInfo->MagEphemInfo->nAlpha) {
            /*
             * Select/deselect all of them
             */
            AllState = gtk_toggle_button_get_active(  GTK_TOGGLE_BUTTON( gInfo->FieldLineShowPitchAngleButton[ObjInfo->MagEphemInfo->nAlpha] ) ); // get state of "All" check item
            printf("AllState  = %d\n", AllState);
            ShowAllPitchAngles = AllState;
            for (j=0; j<ObjInfo->MagEphemInfo->nAlpha; j++) {
                ShowPitchAngle[j] = AllState;
                g_signal_handler_block( G_OBJECT( gInfo->FieldLineShowPitchAngleButton[j] ), gInfo->FieldLineShowPitchAngleButtonHandler[j] );
                gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON( gInfo->FieldLineShowPitchAngleButton[j] ), AllState );
                g_signal_handler_unblock( G_OBJECT( gInfo->FieldLineShowPitchAngleButton[j] ), gInfo->FieldLineShowPitchAngleButtonHandler[j] );
            }

        } else {

           ShowPitchAngle[i] = gtk_toggle_button_get_active(  GTK_TOGGLE_BUTTON( gInfo->FieldLineShowPitchAngleButton[i] ) );

            /*
             *  Test the individual check boxes. If they are all on,
             *  then set the All item to on as well.
             *
             *  If they are all off, then set the All item to off.
             */
            for (n=0, j=0; j<ObjInfo->MagEphemInfo->nAlpha; j++) n += gtk_toggle_button_get_active(  GTK_TOGGLE_BUTTON( gInfo->FieldLineShowPitchAngleButton[j] ) );
            if ( n == ObjInfo->MagEphemInfo->nAlpha ){
                ShowAllPitchAngles = TRUE;
                g_signal_handler_block( G_OBJECT( gInfo->FieldLineShowPitchAngleButton[ObjInfo->MagEphemInfo->nAlpha] ), gInfo->FieldLineShowPitchAngleButtonHandler[ObjInfo->MagEphemInfo->nAlpha] );
                gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON( gInfo->FieldLineShowPitchAngleButton[ObjInfo->MagEphemInfo->nAlpha] ), TRUE );
                g_signal_handler_unblock( G_OBJECT( gInfo->FieldLineShowPitchAngleButton[ObjInfo->MagEphemInfo->nAlpha] ), gInfo->FieldLineShowPitchAngleButtonHandler[ObjInfo->MagEphemInfo->nAlpha] );

            } else {
                ShowAllPitchAngles = FALSE;
                g_signal_handler_block( G_OBJECT( gInfo->FieldLineShowPitchAngleButton[ObjInfo->MagEphemInfo->nAlpha] ), gInfo->FieldLineShowPitchAngleButtonHandler[ObjInfo->MagEphemInfo->nAlpha] );
                gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON( gInfo->FieldLineShowPitchAngleButton[ObjInfo->MagEphemInfo->nAlpha] ), FALSE );
                g_signal_handler_unblock( G_OBJECT( gInfo->FieldLineShowPitchAngleButton[ObjInfo->MagEphemInfo->nAlpha] ), gInfo->FieldLineShowPitchAngleButtonHandler[ObjInfo->MagEphemInfo->nAlpha] );
            }


        }
    } else {
        i -= 100;
        if (i==ObjInfo->MagEphemInfo->nAlpha) {
            /*
             * Select/deselect all of them
             */
            AllState2 = gtk_toggle_button_get_active(  GTK_TOGGLE_BUTTON( gInfo->DriftShellShowPitchAngleButton[ObjInfo->MagEphemInfo->nAlpha] ) ); // get state of "All" check item
            printf("AllState2 = %d\n", AllState2);
            ShowAllPitchAngles2 = AllState2;
            for (j=0; j<ObjInfo->MagEphemInfo->nAlpha; j++) {
                ShowPitchAngle2[j] = AllState2;
                g_signal_handler_block( G_OBJECT( gInfo->DriftShellShowPitchAngleButton[j] ), gInfo->DriftShellShowPitchAngleButtonHandler[j] );
                gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON( gInfo->DriftShellShowPitchAngleButton[j] ), AllState2 );
                g_signal_handler_unblock( G_OBJECT( gInfo->DriftShellShowPitchAngleButton[j] ), gInfo->DriftShellShowPitchAngleButtonHandler[j] );
            }

        } else {

           ShowPitchAngle2[i] = gtk_toggle_button_get_active(  GTK_TOGGLE_BUTTON( gInfo->DriftShellShowPitchAngleButton[i] ) );

            /*
             *  Test the individual check boxes. If they are all on,
             *  then set the All item to on as well.
             *
             *  If they are all off, then set the All item to off.
             */
            for (n=0, j=0; j<ObjInfo->MagEphemInfo->nAlpha; j++) n += gtk_toggle_button_get_active(  GTK_TOGGLE_BUTTON( gInfo->DriftShellShowPitchAngleButton[j] ) );
            if ( n == ObjInfo->MagEphemInfo->nAlpha ){
                ShowAllPitchAngles = TRUE;
                g_signal_handler_block( G_OBJECT( gInfo->DriftShellShowPitchAngleButton[ObjInfo->MagEphemInfo->nAlpha] ), gInfo->DriftShellShowPitchAngleButtonHandler[ObjInfo->MagEphemInfo->nAlpha] );
                gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON( gInfo->DriftShellShowPitchAngleButton[ObjInfo->MagEphemInfo->nAlpha] ), TRUE );
                g_signal_handler_unblock( G_OBJECT( gInfo->DriftShellShowPitchAngleButton[ObjInfo->MagEphemInfo->nAlpha] ), gInfo->DriftShellShowPitchAngleButtonHandler[ObjInfo->MagEphemInfo->nAlpha] );

            } else {
                ShowAllPitchAngles = FALSE;
                g_signal_handler_block( G_OBJECT( gInfo->DriftShellShowPitchAngleButton[ObjInfo->MagEphemInfo->nAlpha] ), gInfo->DriftShellShowPitchAngleButtonHandler[ObjInfo->MagEphemInfo->nAlpha] );
                gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON( gInfo->DriftShellShowPitchAngleButton[ObjInfo->MagEphemInfo->nAlpha] ), FALSE );
                g_signal_handler_unblock( G_OBJECT( gInfo->DriftShellShowPitchAngleButton[ObjInfo->MagEphemInfo->nAlpha] ), gInfo->DriftShellShowPitchAngleButtonHandler[ObjInfo->MagEphemInfo->nAlpha] );
            }


        }
    }
Lgm_CTrans      *c = Lgm_init_ctrans(0);
Lgm_DateTime    dt;
int doy;
StartDate = ObjInfo->MagEphemInfo->Date;
Lgm_Doy( StartDate, &StartYear, &StartMonth, &StartDay, &doy);
StartUT = ObjInfo->MagEphemInfo->UTC;
Lgm_UT_to_hmsms( StartUT, &StartHour, &StartMin, &StartSec, &StartMilliSec );
StartJD    = Lgm_JD( StartYear, StartMonth, StartDay, StartUT, LGM_TIME_SYS_UTC, c );
Lgm_Make_UTC( StartDate, StartUT, &dt, c );
SetStartDate( &dt );
TimeAction( (GtkWidget *)NULL, GINT_TO_POINTER( TIME_RESET_BACKWARD_TO_START ) );
TimeAction( (GtkWidget *)NULL, GINT_TO_POINTER( TIME_RESET_BACKWARD_TO_START ) );
Lgm_free_ctrans( c );
    printf("28.\n"); expose_event( drawing_area, NULL, NULL );

}


static void SetCoordSystem( GtkWidget  *widget, gpointer data ) {

    int i;

    for (i=0; i< 5; i++){
        if ( gtk_toggle_button_get_active(  GTK_TOGGLE_BUTTON( RadioCoordSystem[i] ) ) ) break;
    }

    switch (i){
        case 0: // GSM
                ObserverCoords    = GSM_COORDS;
                StarsConvertFlag  = GEI2000_TO_GSM;
                SatsConvertFlag   = GEI2000_TO_GSM;
                AtmosConvertFlag  = GSM_TO_GSM;
                ObsToGeoConvertFlag = GSM_TO_GEO;
                GeoToObsConvertFlag = GEO_TO_GSM;
                break;
        case 1: // SM
                ObserverCoords    = SM_COORDS;
                StarsConvertFlag  = GEI2000_TO_SM;
                SatsConvertFlag   = GEI2000_TO_SM;
                AtmosConvertFlag  = GSM_TO_SM;
                ObsToGeoConvertFlag = SM_TO_GEO;
                GeoToObsConvertFlag = GEO_TO_SM;
                break;
        case 2: // GSE
                ObserverCoords    = GSE_COORDS;
                StarsConvertFlag  = GEI2000_TO_GSE;
                SatsConvertFlag   = GEI2000_TO_GSE;
                AtmosConvertFlag  = GSM_TO_GSE;
                ObsToGeoConvertFlag = GSE_TO_GEO;
                GeoToObsConvertFlag = GEO_TO_GSE;
                break;
        case 3: // GEI
                ObserverCoords    = GEI2000_COORDS;
                StarsConvertFlag  = GEI2000_TO_GEI2000;
                SatsConvertFlag   = GEI2000_TO_GEI2000;
                AtmosConvertFlag  = GSM_TO_GEI2000;
                ObsToGeoConvertFlag = GEI2000_TO_GEO;
                GeoToObsConvertFlag = GEO_TO_GEI2000;
                break;
        case 4: // GEO
                ObserverCoords    = GEO_COORDS;
                StarsConvertFlag  = GEI2000_TO_GEO;
                SatsConvertFlag   = GEI2000_TO_GEO;
                AtmosConvertFlag  = GSM_TO_GEO;
                ObsToGeoConvertFlag = GEO_TO_GEO;
                GeoToObsConvertFlag = GEO_TO_GEO;
                break;
    }

    TimeAction( (GtkWidget *)NULL, GINT_TO_POINTER(5) );

}

enum {
    IRIDIUM_FLARE_COLUMN_DATE,
    IRIDIUM_FLARE_COLUMN_TIME,
    IRIDIUM_FLARE_COLUMN_MAGNITUDE,
    IRIDIUM_FLARE_COLUMN_ALT,
    IRIDIUM_FLARE_COLUMN_AZ,
    IRIDIUM_FLARE_COLUMN_DIST_TO_CENTER,
    IRIDIUM_FLARE_COLUMN_MAGNITUDE_AT_CENTER,
    IRIDIUM_FLARE_COLUMN_SATELLITE,
    IRIDIUM_FLARE_NUM_COLUMNS,
};


void PredictIridiumFlares( GtkWidget *widget, gpointer data) {

    int         Hours = GPOINTER_TO_INT( data );
    int         Flag1, Flag2, Flag3, i, tYear, tMonth, tDay;
    long int    tDate;
    double      tUT;

    double 		th,ph, r, JD, sJD, eJD, JDinc, diff;
    Lgm_Vector 	u, uu, EarthToSun, g1, g2, g3, P;
    Lgm_CTrans  *tc = Lgm_init_ctrans( 0 );


    printf("PredictIridiumFlares: Hours = %d\n", Hours);
    th =  35.888; ph = -106.306; r  = 0.000; Lgm_GEOD_to_WGS84( th, ph, r, &u );

    sJD = CurrentJD;
    eJD = CurrentJD+Hours/24.0;
    JDinc = (eJD > sJD) ? 10.0/86400.0 : -10.0/86400.0;
    for ( JD = sJD; JD <= eJD; JD += JDinc ){
        //printf("JD = %lf\n", JD);
        Lgm_jd_to_ymdh( JD, &tDate, &tYear, &tMonth, &tDay, &tUT );
        Lgm_Set_Coord_Transforms( tDate, tUT, tc );
        Lgm_Convert_Coords( &u, &uu, GEO_TO_MOD, tc );

printf("here\n");

printf("SpaceObjects->nSat = %d\n", SpaceObjects->nSat);
        for (i=0; i<SpaceObjects->nSat; i++){
            if ( strstr( SpaceObjects->Sat[i].TLE.Name, "IRIDIUM" ) ) {
                EarthToSun = tc->Sun; Lgm_ScaleVector( &EarthToSun, tc->earth_sun_dist );
                IridiumFlare( JD, &SpaceObjects->Sat[i].TLE, &EarthToSun, &Flag1, &Flag2, &Flag3, &g1, &g2, &g3, &P);
diff = Re*Lgm_VecDiffMag( &uu, &g1 );
printf("g1 diff: %g  Date = %8ld  UT = %g    Sat = %s  ( ", diff, tDate, tUT, SpaceObjects->Sat[i].TLE.Name );
                if (Flag1){
                    diff = Re*Lgm_VecDiffMag( &uu, &g1 );
                    if ( diff < 40.0 ) {
                        printf("g1 diff: %g  Date = %8ld  UT = %g    Sat = %s  ( ", diff, tDate, tUT, SpaceObjects->Sat[i].TLE.Name );
                        Lgm_Print_HMSd( tUT ); printf(" )\n");
                    }
                }

                if (Flag2){
                    diff = Re*Lgm_VecDiffMag( &uu, &g2 );
                    if ( diff < 40.0 ) {
                        printf("g2 diff: %g  Date = %8ld  UT = %g    Sat = %s  ( ", diff, tDate, tUT, SpaceObjects->Sat[i].TLE.Name );
                        Lgm_Print_HMSd( tUT ); printf(" )\n");
                    }
                }

                if (Flag3){
                    diff = Re*Lgm_VecDiffMag( &uu, &g3 );
                    if ( diff < 40.0 ) {
                        printf("g3 diff: %g  Date = %8ld  UT = %g    Sat = %s  ( ", diff, tDate, tUT, SpaceObjects->Sat[i].TLE.Name );
                        Lgm_Print_HMSd( tUT ); printf(" )\n");
                    }
                }

            }
        }
    }
printf("done\n");

    Lgm_free_ctrans( tc );

}


GtkWidget *PitchAngleDisplayProperties(){

    int         i, col, row, MaxDay;
    char        Str[320];
    int         ii;
    GtkObject   *spinbutton1_adj;
    GdkColor    color;
    GtkWidget   *frame, *treeview, *scrolledwindow, *filechooserbutton;
    GtkWidget   *window1, *vseparator, *hseparator, *alignment, *image;
    GtkWidget   *vbox2, *vbox3, *table1, *table3, *label, *RowHbox, *hbox;
    GtkWidget   *notebook, *checkbutton, *spinbutton, *button, *colorbutton;
    GSList      *RadioCoordSystemGroup = NULL;
    GSList      *RadioLightingGroup = NULL;
    GtkCellRenderer     *renderer;
    GtkTreeViewColumn   *column;
    GtkTreeModel        *model;
//    GtkWidget           *grid1;




    window1 = gtk_window_new( GTK_WINDOW_TOPLEVEL); gtk_widget_show( window1 );
    gtk_window_set_title( GTK_WINDOW(window1), _("Drift Shell Settings"));
//    gtk_window_set_policy( GTK_WINDOW(window1), FALSE, FALSE, TRUE );
    gtk_window_set_resizable( GTK_WINDOW(window1), TRUE );

    notebook = gtk_notebook_new(); gtk_widget_show( notebook );
    gtk_container_add (GTK_CONTAINER (window1), notebook);


    /*
     * ******* Date/Time  Page *****************************************************
     */
    vbox2 = gtk_vbox_new (FALSE, 0); gtk_widget_show (vbox2);
    gtk_container_set_border_width (GTK_CONTAINER (vbox2), 20);

    // Page title
    label = gtk_label_new (_("<b><span size=\"small\">Date/Time</span></b>")); gtk_widget_show (label);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_CENTER);
    gtk_label_set_use_markup( GTK_LABEL(label), TRUE );
    gtk_notebook_append_page( GTK_NOTEBOOK(notebook), vbox2, label );

    // Table for widgets
    table1 = gtk_table_new (8, 4, FALSE); gtk_widget_show (table1);
    gtk_box_pack_start (GTK_BOX (vbox2), table1, TRUE, TRUE, 15);
    gtk_table_set_row_spacings (GTK_TABLE (table1), 30);
    gtk_table_set_col_spacings (GTK_TABLE (table1), 10);
    row = 0;

    // start date
    label = gtk_label_new (_("Start\nDate:"));
    gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, 0, 1, row, row+1, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);
    gtk_label_set_use_markup (GTK_LABEL (label), TRUE);

    table3 = gtk_table_new (2, 3, FALSE);
    gtk_widget_show (table3);
    gtk_table_attach (GTK_TABLE (table1), table3, 1, 2, row, row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (GTK_FILL), 0, 0);
    gtk_table_set_col_spacings (GTK_TABLE (table3), 10);

    label = gtk_label_new (_("year")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table3), label, 0, 1, 1, 2, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);

    label = gtk_label_new (_("month")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table3), label, 1, 2, 1, 2, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);

    label = gtk_label_new (_("day")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table3), label, 2, 3, 1, 2, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);

    StartYearSpinbutton_adj = gtk_adjustment_new (StartYear, 1950, 2050, 1, 5, 0); // year
    StartYearSpinbutton = gtk_spin_button_new (GTK_ADJUSTMENT (StartYearSpinbutton_adj), 1, 0); gtk_widget_show (StartYearSpinbutton);
    gtk_table_attach (GTK_TABLE (table3), StartYearSpinbutton, 0, 1, 0, 1, (GtkAttachOptions) (GTK_EXPAND | GTK_FILL), (GtkAttachOptions) (0), 0, 0);
    gtk_spin_button_set_numeric (GTK_SPIN_BUTTON (StartYearSpinbutton), TRUE);
    g_signal_connect( G_OBJECT( StartYearSpinbutton ), "value_changed", G_CALLBACK( ChangeStartOrEndDate ), GINT_TO_POINTER(1) );

    StartMonthSpinbutton_adj = gtk_adjustment_new (StartMonth, 1, 12, 1, 6, 0); // month
    StartMonthSpinbutton = gtk_spin_button_new (GTK_ADJUSTMENT (StartMonthSpinbutton_adj), 1, 0); gtk_widget_show (StartMonthSpinbutton);
    gtk_table_attach (GTK_TABLE (table3), StartMonthSpinbutton, 1, 2, 0, 1, (GtkAttachOptions) (GTK_EXPAND | GTK_FILL), (GtkAttachOptions) (0), 0, 0);
    gtk_spin_button_set_numeric (GTK_SPIN_BUTTON (StartMonthSpinbutton), TRUE);
    gtk_spin_button_set_wrap (GTK_SPIN_BUTTON (StartMonthSpinbutton), TRUE);
    g_signal_connect( G_OBJECT( StartMonthSpinbutton ), "value_changed", G_CALLBACK( ChangeStartOrEndDate ), GINT_TO_POINTER(2) );

    MaxDay = MonthDays[StartMonth-1] + ( (StartMonth == 2) ? Lgm_LeapYear(StartYear) : 0 );
    StartDaySpinbutton_adj = gtk_adjustment_new (StartDay, 1, MaxDay, 1, 7, 0); // day
    StartDaySpinbutton = gtk_spin_button_new (GTK_ADJUSTMENT (StartDaySpinbutton_adj), 1, 0); gtk_widget_show (StartDaySpinbutton);
    gtk_table_attach (GTK_TABLE (table3), StartDaySpinbutton, 2, 3, 0, 1, (GtkAttachOptions) (GTK_EXPAND | GTK_FILL), (GtkAttachOptions) (0), 0, 0);
    gtk_spin_button_set_numeric (GTK_SPIN_BUTTON (StartDaySpinbutton), TRUE);
    gtk_spin_button_set_wrap (GTK_SPIN_BUTTON (StartDaySpinbutton), TRUE);
    g_signal_connect( G_OBJECT( StartDaySpinbutton ), "value_changed", G_CALLBACK( ChangeStartOrEndDate ), GINT_TO_POINTER(3) );


    // start time
    label = gtk_label_new (_("Start\nTime:")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, 2, 3, row, row+1, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);
    gtk_label_set_use_markup (GTK_LABEL (label), TRUE);
    gtk_misc_set_alignment (GTK_MISC (label), 1, 0.5);

    table3 = gtk_table_new (2, 3, FALSE); gtk_widget_show (table3);
    gtk_table_attach (GTK_TABLE (table1), table3, 3, 4, row, row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (GTK_FILL), 0, 0);
    gtk_table_set_col_spacings (GTK_TABLE (table3), 10);

    spinbutton1_adj = gtk_adjustment_new (StartHour, 0, 23, 1, 6, 0); // hour
    spinbutton = gtk_spin_button_new (GTK_ADJUSTMENT (spinbutton1_adj), 1, 0); gtk_widget_show (spinbutton);
    gtk_table_attach (GTK_TABLE (table3), spinbutton, 0, 1, 0, 1, (GtkAttachOptions) (GTK_EXPAND | GTK_FILL), (GtkAttachOptions) (0), 0, 0);
    gtk_spin_button_set_numeric (GTK_SPIN_BUTTON (spinbutton), TRUE);
    gtk_spin_button_set_wrap (GTK_SPIN_BUTTON (spinbutton), TRUE);
    g_signal_connect( G_OBJECT( spinbutton ), "value_changed", G_CALLBACK( ChangeStartOrEndTime ), GINT_TO_POINTER(1) );

    spinbutton1_adj = gtk_adjustment_new (StartMin, 0, 59, 1, 10, 0); // min
    spinbutton = gtk_spin_button_new (GTK_ADJUSTMENT (spinbutton1_adj), 1, 0); gtk_widget_show (spinbutton);
    gtk_table_attach (GTK_TABLE (table3), spinbutton, 1, 2, 0, 1, (GtkAttachOptions) (GTK_EXPAND | GTK_FILL), (GtkAttachOptions) (0), 0, 0);
    gtk_spin_button_set_numeric (GTK_SPIN_BUTTON (spinbutton), TRUE);
    gtk_spin_button_set_wrap (GTK_SPIN_BUTTON (spinbutton), TRUE);
    g_signal_connect( G_OBJECT( spinbutton ), "value_changed", G_CALLBACK( ChangeStartOrEndTime ), GINT_TO_POINTER(2) );

    spinbutton1_adj = gtk_adjustment_new (StartSec, 0, 59.99999, 1.0, 10.0, 0); // sec
    spinbutton = gtk_spin_button_new (GTK_ADJUSTMENT (spinbutton1_adj), 1, 4); gtk_widget_show (spinbutton);
    gtk_table_attach (GTK_TABLE (table3), spinbutton, 2, 3, 0, 1, (GtkAttachOptions) (GTK_EXPAND | GTK_FILL), (GtkAttachOptions) (0), 0, 0);
    gtk_spin_button_set_numeric (GTK_SPIN_BUTTON (spinbutton), TRUE);
    gtk_spin_button_set_wrap (GTK_SPIN_BUTTON (spinbutton), TRUE);
    g_signal_connect( G_OBJECT( spinbutton ), "value_changed", G_CALLBACK( ChangeStartOrEndTime ), GINT_TO_POINTER(3) );

    label = gtk_label_new (_("hour")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table3), label, 0, 1, 1, 2, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);
    gtk_label_set_use_markup (GTK_LABEL (label), TRUE);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_CENTER);

    label = gtk_label_new (_("min")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table3), label, 1, 2, 1, 2, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);

    label = gtk_label_new (_("sec")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table3), label, 2, 3, 1, 2, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);




    // end date
    ++row;
    label = gtk_label_new (_("End\nDate:"));
    gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, 0, 1, row, row+1, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);
    gtk_label_set_use_markup (GTK_LABEL (label), TRUE);

    table3 = gtk_table_new (2, 3, FALSE);
    gtk_widget_show (table3);
    gtk_table_attach (GTK_TABLE (table1), table3, 1, 2, row, row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (GTK_FILL), 0, 0);
    gtk_table_set_col_spacings (GTK_TABLE (table3), 10);

    label = gtk_label_new (_("year")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table3), label, 0, 1, 1, 2, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);

    label = gtk_label_new (_("month")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table3), label, 1, 2, 1, 2, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);

    label = gtk_label_new (_("day")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table3), label, 2, 3, 1, 2, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);

    spinbutton1_adj = gtk_adjustment_new (EndYear, 1950, 2050, 1, 5, 0); // year
    spinbutton = gtk_spin_button_new (GTK_ADJUSTMENT (spinbutton1_adj), 1, 0); gtk_widget_show (spinbutton);
    gtk_table_attach (GTK_TABLE (table3), spinbutton, 0, 1, 0, 1, (GtkAttachOptions) (GTK_EXPAND | GTK_FILL), (GtkAttachOptions) (0), 0, 0);
    gtk_spin_button_set_numeric (GTK_SPIN_BUTTON (spinbutton), TRUE);
    g_signal_connect( G_OBJECT( spinbutton ), "value_changed", G_CALLBACK( ChangeStartOrEndDate ), GINT_TO_POINTER(11) );

    spinbutton1_adj = gtk_adjustment_new (EndMonth, 1, 12, 1, 6, 0); // month
    spinbutton = gtk_spin_button_new (GTK_ADJUSTMENT (spinbutton1_adj), 1, 0); gtk_widget_show (spinbutton);
    gtk_table_attach (GTK_TABLE (table3), spinbutton, 1, 2, 0, 1, (GtkAttachOptions) (GTK_EXPAND | GTK_FILL), (GtkAttachOptions) (0), 0, 0);
    gtk_spin_button_set_numeric (GTK_SPIN_BUTTON (spinbutton), TRUE);
    gtk_spin_button_set_wrap (GTK_SPIN_BUTTON (spinbutton), TRUE);
    g_signal_connect( G_OBJECT( spinbutton ), "value_changed", G_CALLBACK( ChangeStartOrEndDate ), GINT_TO_POINTER(12) );

    MaxDay = MonthDays[EndMonth-1] + ( (EndMonth == 2) ? Lgm_LeapYear(EndYear) : 0 );
    EndDaySpinbutton_adj = gtk_adjustment_new (EndDay, 1, MaxDay, 1, 7, 0); // day
    EndDaySpinbutton = gtk_spin_button_new (GTK_ADJUSTMENT (EndDaySpinbutton_adj), 1, 0); gtk_widget_show (EndDaySpinbutton);
    gtk_table_attach (GTK_TABLE (table3), EndDaySpinbutton, 2, 3, 0, 1, (GtkAttachOptions) (GTK_EXPAND | GTK_FILL), (GtkAttachOptions) (0), 0, 0);
    gtk_spin_button_set_numeric (GTK_SPIN_BUTTON (EndDaySpinbutton), TRUE);
    gtk_spin_button_set_wrap (GTK_SPIN_BUTTON (EndDaySpinbutton), TRUE);
    g_signal_connect( G_OBJECT( EndDaySpinbutton ), "value_changed", G_CALLBACK( ChangeStartOrEndDate ), GINT_TO_POINTER(13) );


    // end time
    label = gtk_label_new (_("End\nTime:")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, 2, 3, row, row+1, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);
    gtk_label_set_use_markup (GTK_LABEL (label), TRUE);
    gtk_misc_set_alignment (GTK_MISC (label), 1, 0.5);

    table3 = gtk_table_new (2, 3, FALSE); gtk_widget_show (table3);
    gtk_table_attach (GTK_TABLE (table1), table3, 3, 4, row, row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (GTK_FILL), 0, 0);
    gtk_table_set_col_spacings (GTK_TABLE (table3), 10);

    spinbutton1_adj = gtk_adjustment_new (EndHour, 0, 23, 1, 6, 0); // hour
    spinbutton = gtk_spin_button_new (GTK_ADJUSTMENT (spinbutton1_adj), 1, 0); gtk_widget_show (spinbutton);
    gtk_table_attach (GTK_TABLE (table3), spinbutton, 0, 1, 0, 1, (GtkAttachOptions) (GTK_EXPAND | GTK_FILL), (GtkAttachOptions) (0), 0, 0);
    gtk_spin_button_set_numeric (GTK_SPIN_BUTTON (spinbutton), TRUE);
    gtk_spin_button_set_wrap (GTK_SPIN_BUTTON (spinbutton), TRUE);
    g_signal_connect( G_OBJECT( spinbutton ), "value_changed", G_CALLBACK( ChangeStartOrEndTime ), GINT_TO_POINTER(11) );

    spinbutton1_adj = gtk_adjustment_new (EndMin, 0, 59, 1, 10, 0); // min
    spinbutton = gtk_spin_button_new (GTK_ADJUSTMENT (spinbutton1_adj), 1, 0); gtk_widget_show (spinbutton);
    gtk_table_attach (GTK_TABLE (table3), spinbutton, 1, 2, 0, 1, (GtkAttachOptions) (GTK_EXPAND | GTK_FILL), (GtkAttachOptions) (0), 0, 0);
    gtk_spin_button_set_numeric (GTK_SPIN_BUTTON (spinbutton), TRUE);
    gtk_spin_button_set_wrap (GTK_SPIN_BUTTON (spinbutton), TRUE);
    g_signal_connect( G_OBJECT( spinbutton ), "value_changed", G_CALLBACK( ChangeStartOrEndTime ), GINT_TO_POINTER(12) );

    spinbutton1_adj = gtk_adjustment_new (EndSec, 0, 59.99999, 1.0, 10.0, 0); // sec
    spinbutton = gtk_spin_button_new (GTK_ADJUSTMENT (spinbutton1_adj), 1, 4); gtk_widget_show (spinbutton);
    gtk_table_attach (GTK_TABLE (table3), spinbutton, 2, 3, 0, 1, (GtkAttachOptions) (GTK_EXPAND | GTK_FILL), (GtkAttachOptions) (0), 0, 0);
    gtk_spin_button_set_numeric (GTK_SPIN_BUTTON (spinbutton), TRUE);
    gtk_spin_button_set_wrap (GTK_SPIN_BUTTON (spinbutton), TRUE);
    g_signal_connect( G_OBJECT( spinbutton ), "value_changed", G_CALLBACK( ChangeStartOrEndTime ), GINT_TO_POINTER(13) );

    label = gtk_label_new (_("hour")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table3), label, 0, 1, 1, 2, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);
    gtk_label_set_use_markup (GTK_LABEL (label), TRUE);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_CENTER);

    label = gtk_label_new (_("min")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table3), label, 1, 2, 1, 2, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);

    label = gtk_label_new (_("sec")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table3), label, 2, 3, 1, 2, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);



    // time increment + number of frames
    ++row;
    label = gtk_label_new (_("Time Inc.\n(seconds):")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, 0, 1, row, row+1, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);
    gtk_label_set_use_markup (GTK_LABEL (label), TRUE);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_RIGHT);
    gtk_misc_set_alignment (GTK_MISC (label), 0, 0.5);

    spinbutton1_adj = gtk_adjustment_new (TimeInc, 0, 86400.0, 1, 60, 0);
    spinbutton = gtk_spin_button_new (GTK_ADJUSTMENT (spinbutton1_adj), 1, 4);
    gtk_widget_show (spinbutton);
    gtk_table_attach (GTK_TABLE (table1), spinbutton, 1, 2, row, row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
    gtk_spin_button_set_numeric (GTK_SPIN_BUTTON (spinbutton), TRUE);
    g_signal_connect( G_OBJECT( spinbutton ), "value_changed", G_CALLBACK( ChangeTimeInc ), NULL );


    // dump frames
    DumpFramesCheckbutton = gtk_check_button_new_with_mnemonic (_("Dump Frames")); gtk_widget_show (DumpFramesCheckbutton);
    gtk_table_attach (GTK_TABLE (table1), DumpFramesCheckbutton, 2, 4, row, row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
    if (DumpFrames) gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON( DumpFramesCheckbutton ), DumpFrames );
    g_signal_connect( G_OBJECT( DumpFramesCheckbutton ), "toggled", G_CALLBACK( ToggleDumpFrames ), NULL );

    // Current Frame Label
    ++row;
    //sprintf(Str, "Current Frame: %ld / %ld", cFrame, nFrames);
    sprintf(Str, "Frames Done: %ld    Frames Remaining: %ld    (Total: %ld)", cFrame, nFramesLeft, nFrames);
    cFramesLabel = gtk_label_new (Str); gtk_widget_show (cFramesLabel);
    gtk_table_attach (GTK_TABLE (table1), cFramesLabel, 0, 4, row, row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
    gtk_misc_set_alignment (GTK_MISC (cFramesLabel), 0, 0.5);


    // Start Time Label
    ++row;
    sprintf(Str, "Start Time: %4d/%02d/%02d  %02d:%02d:%02d.%03d\n  UT = %.8lf  JD = %.8lf", StartYear, StartMonth, StartDay,
                                    StartHour, StartMin, StartSec, StartMilliSec, StartUT, StartJD );
    StartTimeLabel = gtk_label_new( Str ); gtk_widget_show(StartTimeLabel);
    gtk_table_attach (GTK_TABLE (table1), StartTimeLabel, 0, 4, row, row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
    gtk_misc_set_alignment (GTK_MISC (StartTimeLabel), 0, 0.5);

    // End Time Label
    ++row;
    sprintf(Str, "End Time: %4d/%02d/%02d  %02d:%02d:%02d.%03d\n  UT = %.8lf  JD = %.8lf", EndYear, EndMonth, EndDay,
                                    EndHour, EndMin, EndSec, EndMilliSec, EndUT, EndJD );
    EndTimeLabel = gtk_label_new( Str ); gtk_widget_show(EndTimeLabel);
    gtk_table_attach (GTK_TABLE (table1), EndTimeLabel, 0, 4, row, row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
    gtk_misc_set_alignment (GTK_MISC (EndTimeLabel), 0, 0.5);

    // Current Time Label
    ++row;
    sprintf(Str, "Current Time: %4d/%02d/%02d  %02d:%02d:%02d.%03d\n  UT = %.8lf  JD = %.8lf", CurrentYear, CurrentMonth, CurrentDay,
                                    CurrentHour, CurrentMin, CurrentSec, CurrentMilliSec, CurrentUT, CurrentJD );
    CurrentTimeLabel = gtk_label_new( Str ); gtk_widget_show(CurrentTimeLabel);
    gtk_table_attach (GTK_TABLE (table1), CurrentTimeLabel, 0, 4, row, row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
    gtk_misc_set_alignment (GTK_MISC(CurrentTimeLabel), 0, 0.5);


    // ********* action buttons ************
    ++row;
    hbox = gtk_hbox_new (FALSE, 10); gtk_widget_show (hbox);
    gtk_table_attach (GTK_TABLE (table1), hbox, 0, 4, row, row+1, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);

//    button = gtk_button_new_with_mnemonic (_("Reset Current Time\n To Start Time")); gtk_widget_show (button);

    button = gtk_button_new( ); gtk_widget_show( button );
    gtk_box_pack_start (GTK_BOX (hbox), button, FALSE, FALSE, 0);
    alignment = gtk_alignment_new( 0.5, 0.5, 0, 0 ); gtk_widget_show( alignment );
    gtk_container_add( GTK_CONTAINER( button ), alignment );
    vbox2 = gtk_vbox_new( FALSE, 2 ); gtk_widget_show( vbox2 );
    gtk_container_add( GTK_CONTAINER( alignment ), vbox2 );
//    image = gtk_image_new_from_stock( "gtk-goto-first", GTK_ICON_SIZE_BUTTON); gtk_widget_show( image );
    image = create_pixmap(window1, "go-first.png"); gtk_widget_show (image);
    gtk_box_pack_start( GTK_BOX( vbox2 ), image, FALSE, FALSE, 0 );
    label = gtk_label_new_with_mnemonic( _("Start") ); gtk_widget_show( label );
    gtk_box_pack_start( GTK_BOX( vbox2 ), label, FALSE, FALSE, 0 );
    g_signal_connect( G_OBJECT( button ), "clicked", G_CALLBACK( TimeAction ), GINT_TO_POINTER( TIME_RESET_BACKWARD_TO_START ) );

    button = gtk_button_new( ); gtk_widget_show( button );
    gtk_box_pack_start (GTK_BOX (hbox), button, FALSE, FALSE, 0);
    alignment = gtk_alignment_new( 0.5, 0.5, 0, 0 ); gtk_widget_show( alignment );
    gtk_container_add( GTK_CONTAINER( button ), alignment );
    vbox2 = gtk_vbox_new( FALSE, 2 ); gtk_widget_show( vbox2 );
    gtk_container_add( GTK_CONTAINER( alignment ), vbox2 );
//    image = gtk_image_new_from_stock( "gtk-go-back", GTK_ICON_SIZE_BUTTON); gtk_widget_show( image );
    image = create_pixmap(window1, "go-previous.png"); gtk_widget_show (image);
    gtk_box_pack_start( GTK_BOX( vbox2 ), image, FALSE, FALSE, 0 );
    label = gtk_label_new_with_mnemonic( _("Play") ); gtk_widget_show( label );
    gtk_box_pack_start( GTK_BOX( vbox2 ), label, FALSE, FALSE, 0 );
    g_signal_connect( G_OBJECT( button ), "clicked", G_CALLBACK( TimeAction ), GINT_TO_POINTER( TIME_PLAY_BACKWARD ) );

    button = gtk_button_new( ); gtk_widget_show( button );
    gtk_box_pack_start (GTK_BOX (hbox), button, FALSE, FALSE, 0);
    alignment = gtk_alignment_new( 0.5, 0.5, 0, 0 ); gtk_widget_show( alignment );
    gtk_container_add( GTK_CONTAINER( button ), alignment );
    vbox2 = gtk_vbox_new( FALSE, 2 ); gtk_widget_show( vbox2 );
    gtk_container_add( GTK_CONTAINER( alignment ), vbox2 );
//    image = gtk_image_new_from_stock( "gtk-go-back", GTK_ICON_SIZE_BUTTON); gtk_widget_show( image );
    image = create_pixmap(window1, "go-previous.png"); gtk_widget_show (image);
    gtk_box_pack_start( GTK_BOX( vbox2 ), image, FALSE, FALSE, 0 );
    label = gtk_label_new_with_mnemonic( _("Step") ); gtk_widget_show( label );
    gtk_box_pack_start( GTK_BOX( vbox2 ), label, FALSE, FALSE, 0 );
    g_signal_connect( G_OBJECT( button ), "clicked", G_CALLBACK( TimeAction ), GINT_TO_POINTER( TIME_STEP_BACKWARD ) );

    button = gtk_button_new( ); gtk_widget_show( button );
    gtk_box_pack_start (GTK_BOX (hbox), button, FALSE, FALSE, 0);
    alignment = gtk_alignment_new( 0.5, 0.5, 0, 0 ); gtk_widget_show( alignment );
    gtk_container_add( GTK_CONTAINER( button ), alignment );
    vbox2 = gtk_vbox_new( FALSE, 2 ); gtk_widget_show( vbox2 );
    gtk_container_add( GTK_CONTAINER( alignment ), vbox2 );
//    image = gtk_image_new_from_stock( "gtk-stop", GTK_ICON_SIZE_BUTTON); gtk_widget_show( image );
    image = create_pixmap(window1, "window-close.png"); gtk_widget_show (image);
    gtk_box_pack_start( GTK_BOX( vbox2 ), image, FALSE, FALSE, 0 );
    label = gtk_label_new_with_mnemonic( _("Stop") ); gtk_widget_show( label );
    gtk_box_pack_start( GTK_BOX( vbox2 ), label, FALSE, FALSE, 0 );
    g_signal_connect( G_OBJECT( button ), "clicked", G_CALLBACK( TimeAction ), GINT_TO_POINTER( TIME_STOP ) );


    button = gtk_button_new( ); gtk_widget_show( button );
    gtk_box_pack_start (GTK_BOX (hbox), button, FALSE, FALSE, 0);
    alignment = gtk_alignment_new( 0.5, 0.5, 0, 0 ); gtk_widget_show( alignment );
    gtk_container_add( GTK_CONTAINER( button ), alignment );
    vbox2 = gtk_vbox_new( FALSE, 2 ); gtk_widget_show( vbox2 );
    gtk_container_add( GTK_CONTAINER( alignment ), vbox2 );
//    image = gtk_image_new_from_stock( "gtk-go-forward", GTK_ICON_SIZE_BUTTON); gtk_widget_show( image );
    image = create_pixmap(window1, "go-next.png"); gtk_widget_show (image);
    gtk_box_pack_start( GTK_BOX( vbox2 ), image, FALSE, FALSE, 0 );
    label = gtk_label_new_with_mnemonic( _("Step") ); gtk_widget_show( label );
    gtk_box_pack_start( GTK_BOX( vbox2 ), label, FALSE, FALSE, 0 );
    g_signal_connect( G_OBJECT( button ), "clicked", G_CALLBACK( TimeAction ), GINT_TO_POINTER( TIME_STEP_FOREWARD ) );

    button = gtk_button_new( ); gtk_widget_show( button );
    gtk_box_pack_start (GTK_BOX (hbox), button, FALSE, FALSE, 0);
    alignment = gtk_alignment_new( 0.5, 0.5, 0, 0 ); gtk_widget_show( alignment );
    gtk_container_add( GTK_CONTAINER( button ), alignment );
    vbox2 = gtk_vbox_new( FALSE, 2 ); gtk_widget_show( vbox2 );
    gtk_container_add( GTK_CONTAINER( alignment ), vbox2 );
//    image = gtk_image_new_from_stock( "gtk-go-forward", GTK_ICON_SIZE_BUTTON); gtk_widget_show( image );
    image = create_pixmap(window1, "go-next.png"); gtk_widget_show (image);
    gtk_box_pack_start( GTK_BOX( vbox2 ), image, FALSE, FALSE, 0 );
    label = gtk_label_new_with_mnemonic( _("Play") ); gtk_widget_show( label );
    gtk_box_pack_start( GTK_BOX( vbox2 ), label, FALSE, FALSE, 0 );
    g_signal_connect( G_OBJECT( button ), "clicked", G_CALLBACK( TimeAction ), GINT_TO_POINTER( TIME_PLAY_FOREWARD ) );

    button = gtk_button_new( ); gtk_widget_show( button );
    gtk_box_pack_start (GTK_BOX (hbox), button, FALSE, FALSE, 0);
    alignment = gtk_alignment_new( 0.5, 0.5, 0, 0 ); gtk_widget_show( alignment );
    gtk_container_add( GTK_CONTAINER( button ), alignment );
    vbox2 = gtk_vbox_new( FALSE, 2 ); gtk_widget_show( vbox2 );
    gtk_container_add( GTK_CONTAINER( alignment ), vbox2 );
    //image = gtk_image_new_from_stock( "gtk-go-last.png", GTK_ICON_SIZE_BUTTON); gtk_widget_show( image );
    image = create_pixmap(window1, "go-last.png"); gtk_widget_show (image);
    gtk_box_pack_start( GTK_BOX( vbox2 ), image, FALSE, FALSE, 0 );
    label = gtk_label_new_with_mnemonic( _("End") ); gtk_widget_show( label );
    gtk_box_pack_start( GTK_BOX( vbox2 ), label, FALSE, FALSE, 0 );
    g_signal_connect( G_OBJECT( button ), "clicked", G_CALLBACK( TimeAction ), GINT_TO_POINTER( TIME_RESET_FOREWARD_TO_END ) );


    button = gtk_button_new(); gtk_widget_show(button);
    gtk_box_pack_start (GTK_BOX (hbox), button, FALSE, FALSE, 0);
    alignment = gtk_alignment_new(0.5, 0.5, 0, 0); gtk_widget_show(alignment);
    gtk_container_add (GTK_CONTAINER (button), alignment);
    vbox2 = gtk_vbox_new (FALSE, 2); gtk_widget_show (vbox2);
    gtk_container_add (GTK_CONTAINER (alignment), vbox2);
    image = create_pixmap(window1, "system-run.png"); gtk_widget_show (image);
    gtk_box_pack_start (GTK_BOX (vbox2), image, FALSE, FALSE, 0);
    label = gtk_label_new_with_mnemonic (_("RealTime")); gtk_widget_show (label);
    gtk_box_pack_start (GTK_BOX (vbox2), label, FALSE, FALSE, 0);
    g_signal_connect( G_OBJECT( button ), "clicked", G_CALLBACK( TimeAction ), GINT_TO_POINTER( TIME_REALTIMEPLAY ) );




    /*
     * ******* Coordinate Systems Page *****************************************************
     */
    vbox2 = gtk_vbox_new (FALSE, 0); gtk_widget_show (vbox2);
    gtk_container_set_border_width (GTK_CONTAINER (vbox2), 20);

    // Page title
    label = gtk_label_new (_("<b><span size=\"small\">Coord Systems</span></b>")); gtk_widget_show (label);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_CENTER);
    gtk_label_set_use_markup( GTK_LABEL(label), TRUE );
    gtk_notebook_append_page( GTK_NOTEBOOK(notebook), vbox2, label );

    // Table for widgets
    table1 = gtk_table_new (8, 4, FALSE); gtk_widget_show (table1);
    gtk_box_pack_start (GTK_BOX (vbox2), table1, TRUE, TRUE, 15);
    gtk_container_set_border_width (GTK_CONTAINER (table1), 30);
    gtk_table_set_row_spacings (GTK_TABLE (table1), 10);
    gtk_table_set_col_spacings (GTK_TABLE (table1), 40);
    row = 0;

    /*
     *  Coordinate System
     */
    label = gtk_label_new (_("<b>Coordinate System</b>")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, 0, 1, 0, 1, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);
    gtk_label_set_use_markup (GTK_LABEL (label), TRUE);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_CENTER);

    RadioCoordSystem[0] = gtk_radio_button_new_with_mnemonic (NULL, _("GSM")); gtk_widget_show (RadioCoordSystem[0]);
    gtk_table_attach (GTK_TABLE (table1), RadioCoordSystem[0], 0, 1, 1, 2, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 20, 0);
    gtk_radio_button_set_group (GTK_RADIO_BUTTON (RadioCoordSystem[0]), RadioCoordSystemGroup);
    RadioCoordSystemGroup = gtk_radio_button_get_group (GTK_RADIO_BUTTON (RadioCoordSystem[0]));
    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (RadioCoordSystem[0]), TRUE);
    g_signal_connect( G_OBJECT( RadioCoordSystem[0] ), "toggled", G_CALLBACK( SetCoordSystem ), GINT_TO_POINTER( 1 ) );

    RadioCoordSystem[1] = gtk_radio_button_new_with_mnemonic (NULL, _("SM")); gtk_widget_show (RadioCoordSystem[1]);
    gtk_table_attach (GTK_TABLE (table1), RadioCoordSystem[1], 0, 1, 2, 3, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 20, 0);
    gtk_radio_button_set_group (GTK_RADIO_BUTTON (RadioCoordSystem[1]), RadioCoordSystemGroup);
    RadioCoordSystemGroup = gtk_radio_button_get_group (GTK_RADIO_BUTTON (RadioCoordSystem[1]));
    g_signal_connect( G_OBJECT( RadioCoordSystem[1] ), "toggled", G_CALLBACK( SetCoordSystem ), GINT_TO_POINTER( 2 ) );

    RadioCoordSystem[2] = gtk_radio_button_new_with_mnemonic (NULL, _("GSE")); gtk_widget_show (RadioCoordSystem[2]);
    gtk_table_attach (GTK_TABLE (table1), RadioCoordSystem[2], 0, 1, 3, 4, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 20, 0);
    gtk_radio_button_set_group (GTK_RADIO_BUTTON (RadioCoordSystem[2]), RadioCoordSystemGroup);
    RadioCoordSystemGroup = gtk_radio_button_get_group (GTK_RADIO_BUTTON (RadioCoordSystem[2]));
    g_signal_connect( G_OBJECT( RadioCoordSystem[2] ), "toggled", G_CALLBACK( SetCoordSystem ), GINT_TO_POINTER( 3 ) );

    RadioCoordSystem[3] = gtk_radio_button_new_with_mnemonic (NULL, _("GEI2000")); gtk_widget_show (RadioCoordSystem[3]);
    gtk_table_attach (GTK_TABLE (table1), RadioCoordSystem[3], 0, 1, 4, 5, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 20, 0);
    gtk_radio_button_set_group (GTK_RADIO_BUTTON (RadioCoordSystem[3]), RadioCoordSystemGroup);
    RadioCoordSystemGroup = gtk_radio_button_get_group (GTK_RADIO_BUTTON (RadioCoordSystem[3]));
    g_signal_connect( G_OBJECT( RadioCoordSystem[3] ), "toggled", G_CALLBACK( SetCoordSystem ), GINT_TO_POINTER( 4 ) );

    RadioCoordSystem[4] = gtk_radio_button_new_with_mnemonic (NULL, _("WGS84 (GEO)")); gtk_widget_show (RadioCoordSystem[4]);
    gtk_table_attach (GTK_TABLE (table1), RadioCoordSystem[4], 0, 1, 5, 6, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 20, 0);
    gtk_radio_button_set_group (GTK_RADIO_BUTTON (RadioCoordSystem[4]), RadioCoordSystemGroup);
    RadioCoordSystemGroup = gtk_radio_button_get_group (GTK_RADIO_BUTTON (RadioCoordSystem[4]));
    g_signal_connect( G_OBJECT( RadioCoordSystem[4] ), "toggled", G_CALLBACK( SetCoordSystem ), GINT_TO_POINTER( 5 ) );



    /*
     *  Show Axes
     */
    label = gtk_label_new (_("<b>Show Axes</b>")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, 1, 2, 0, 1, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);
    gtk_label_set_use_markup (GTK_LABEL (label), TRUE);

    checkbutton = gtk_check_button_new_with_mnemonic (_("GSM")); gtk_widget_show (checkbutton);
    gtk_table_attach (GTK_TABLE (table1), checkbutton, 1, 2, 1, 2, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);
    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (checkbutton), TRUE);

    checkbutton = gtk_check_button_new_with_mnemonic (_("SM")); gtk_widget_show (checkbutton);
    gtk_table_attach (GTK_TABLE (table1), checkbutton, 1, 2, 2, 3, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);

    checkbutton = gtk_check_button_new_with_mnemonic (_("GSE")); gtk_widget_show (checkbutton);
    gtk_table_attach (GTK_TABLE (table1), checkbutton, 1, 2, 3, 4, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);

    checkbutton = gtk_check_button_new_with_mnemonic (_("GEI")); gtk_widget_show (checkbutton);
    gtk_table_attach (GTK_TABLE (table1), checkbutton, 1, 2, 4, 5, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);

    checkbutton = gtk_check_button_new_with_mnemonic (_("GEO")); gtk_widget_show (checkbutton);
    gtk_table_attach (GTK_TABLE (table1), checkbutton, 1, 2, 5, 6, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);




    /*
     *  Show Eq Plane grid
     */
    label = gtk_label_new (_("<b>Show Eq. Grid</b>")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, 2, 3, 0, 1, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);
    gtk_label_set_use_markup (GTK_LABEL (label), TRUE);

    checkbutton = gtk_check_button_new_with_mnemonic (_("GSM")); gtk_widget_show (checkbutton);
    gtk_table_attach (GTK_TABLE (table1), checkbutton, 2, 3, 1, 2, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);
    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (checkbutton), TRUE);

    checkbutton = gtk_check_button_new_with_mnemonic (_("SM")); gtk_widget_show (checkbutton);
    gtk_table_attach (GTK_TABLE (table1), checkbutton, 2, 3, 2, 3, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);

    checkbutton = gtk_check_button_new_with_mnemonic (_("GSE")); gtk_widget_show (checkbutton);
    gtk_table_attach (GTK_TABLE (table1), checkbutton, 2, 3, 3, 4, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);

    checkbutton = gtk_check_button_new_with_mnemonic (_("GEI")); gtk_widget_show (checkbutton);
    gtk_table_attach (GTK_TABLE (table1), checkbutton, 2, 3, 4, 5, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);

    checkbutton = gtk_check_button_new_with_mnemonic (_("GEO")); gtk_widget_show (checkbutton);
    gtk_table_attach (GTK_TABLE (table1), checkbutton, 2, 3, 5, 6, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);


    /*
     *  Color of Eq Plane grid
     */
    label = gtk_label_new (_("<b>Grid Color</b>")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, 3, 4, 0, 1, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);
    gtk_label_set_use_markup (GTK_LABEL (label), TRUE);

    colorbutton = gtk_color_button_new (); gtk_widget_show (colorbutton);
    gtk_table_attach (GTK_TABLE (table1), colorbutton, 3, 4, 1, 2, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);

    colorbutton = gtk_color_button_new (); gtk_widget_show (colorbutton);
    gtk_table_attach (GTK_TABLE (table1), colorbutton, 3, 4, 2, 3, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);

    colorbutton = gtk_color_button_new (); gtk_widget_show (colorbutton);
    gtk_table_attach (GTK_TABLE (table1), colorbutton, 3, 4, 3, 4, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);

    colorbutton = gtk_color_button_new (); gtk_widget_show (colorbutton);
    gtk_table_attach (GTK_TABLE (table1), colorbutton, 3, 4, 4, 5, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);

    colorbutton = gtk_color_button_new (); gtk_widget_show (colorbutton);
    gtk_table_attach (GTK_TABLE (table1), colorbutton, 3, 4, 5, 6, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);









    /*
     * ******** Field Lines Page ****************************************************
     */
    vbox2 = gtk_vbox_new (FALSE, 0); gtk_widget_show (vbox2);
    gtk_container_set_border_width (GTK_CONTAINER (vbox2), 20);

    // Page title
    label = gtk_label_new (_("<b><span size=\"small\">Field Lines</span></b>")); gtk_widget_show (label);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_CENTER);
    gtk_label_set_use_markup( GTK_LABEL(label), TRUE );
    gtk_notebook_append_page( GTK_NOTEBOOK(notebook), vbox2, label );

    // Table for widgets
    table1 = gtk_table_new( 40, 10, FALSE ); gtk_widget_show (table1);
    gtk_box_pack_start (GTK_BOX (vbox2), table1, TRUE, TRUE, 15);
    gtk_table_set_row_spacings( GTK_TABLE (table1), 0 );
    gtk_table_set_col_spacings( GTK_TABLE (table1), 10 );


    col = 0;
    // Pitch angle column
    label = gtk_label_new (_("<b>Pitch\nAngle</b>")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, col, col+1, 0, 2, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions)(GTK_FILL), 0, 0);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_CENTER);
    gtk_label_set_use_markup( GTK_LABEL(label), TRUE );
    ++col;

    // vertical separator
    vseparator = gtk_vseparator_new (); gtk_widget_show (vseparator);
    gtk_table_attach (GTK_TABLE (table1), vseparator, col, col+1, 0, 13, (GtkAttachOptions)(GTK_FILL), (GtkAttachOptions)(GTK_FILL), 0, 0);
    ++col;

    // Column Headers for Field Lines
    label = gtk_label_new (_("<b>Show\nField\nLines</b>")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, col, col+1, 0, 2, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions)(GTK_FILL), 0, 0);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_CENTER);
    gtk_label_set_use_markup( GTK_LABEL(label), TRUE );
    ++col;

    label = gtk_label_new (_("<b><span size=\"large\">Field Line Material</span></b>")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, 3, 8, 0, 1, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (GTK_FILL), 0, 0);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_CENTER);
    gtk_label_set_use_markup( GTK_LABEL(label), TRUE );

    label = gtk_label_new (_("<b><span size=\"small\">Material</span></b>")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, col, col+1, 1, 2, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (GTK_FILL), 0, 0);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_CENTER);
    gtk_label_set_use_markup( GTK_LABEL(label), TRUE );
    ++col;

    label = gtk_label_new (_("<b><span size=\"small\">Diffuse</span></b>")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, col, col+1, 1, 2, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (GTK_FILL), 0, 0);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_CENTER);
    gtk_label_set_use_markup( GTK_LABEL(label), TRUE );
    ++col;

    label = gtk_label_new (_("<b><span size=\"small\">Ambient</span></b>")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, col, col+1, 1, 2, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (GTK_FILL), 0, 0);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_CENTER);
    gtk_label_set_use_markup( GTK_LABEL(label), TRUE );
    ++col;

    label = gtk_label_new (_("<b><span size=\"small\">Specular</span></b>")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, col, col+1, 1, 2, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (GTK_FILL), 0, 0);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_CENTER);
    gtk_label_set_use_markup( GTK_LABEL(label), TRUE );
    ++col;

    label = gtk_label_new (_("<b><span size=\"small\">Shininess</span></b>")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, col, col+1, 1, 2, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (GTK_FILL), 0, 0);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_CENTER);
    gtk_label_set_use_markup( GTK_LABEL(label), TRUE );
    ++col;

    // Separator
    hseparator = gtk_hseparator_new (); gtk_widget_show (hseparator);
    gtk_table_attach (GTK_TABLE (table1), hseparator, 0, 8, 2, 3, (GtkAttachOptions) (GTK_EXPAND | GTK_FILL), (GtkAttachOptions) (GTK_FILL), 0, 0);



    /*
     *  Fill in rows
     */
    for (i=0; i<=ObjInfo->MagEphemInfo->nAlpha; i++){

        col = 0;

        if (i==ObjInfo->MagEphemInfo->nAlpha){
            sprintf( Str, "<small><b>All</b></small>" );
        } else {
            sprintf( Str, "<small><b>%g\u00b0</b></small>", ObjInfo->MagEphemInfo->Alpha[i] );
        }
        label = gtk_label_new( Str ); gtk_widget_show( label );
        gtk_table_attach( GTK_TABLE(table1), label, col, col+1, 3+i, 4+i, (GtkAttachOptions)(GTK_FILL), (GtkAttachOptions)(GTK_SHRINK), 0, 0 );
        gtk_label_set_use_markup( GTK_LABEL(label), TRUE );
        gtk_label_set_justify( GTK_LABEL(label), GTK_JUSTIFY_CENTER);
        ++col;
        ++col;


        /*
         *  Field Lines
         */
        gInfo->FieldLineShowPitchAngleButton[i] = gtk_check_button_new(); gtk_widget_show( gInfo->FieldLineShowPitchAngleButton[i] );
        gtk_widget_set_size_request( gInfo->FieldLineShowPitchAngleButton[i], 30, 20);
        if (i==ObjInfo->MagEphemInfo->nAlpha) gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON( gInfo->FieldLineShowPitchAngleButton[i] ), ShowAllPitchAngles );
        else      gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON( gInfo->FieldLineShowPitchAngleButton[i] ), ShowPitchAngle[i] );
        gtk_table_attach( GTK_TABLE(table1), gInfo->FieldLineShowPitchAngleButton[i], col, col+1, 3+i, 4+i, (GtkAttachOptions)(GTK_FILL), (GtkAttachOptions)(GTK_SHRINK), 0, 0 );
        gInfo->FieldLineShowPitchAngleButtonHandler[i] = g_signal_connect( G_OBJECT( gInfo->FieldLineShowPitchAngleButton[i] ), "toggled", G_CALLBACK( SelectPitchAngles2 ), GINT_TO_POINTER(i) );
        ++col;



        GtkTreeIter     iter;
        GtkTreeStore    *ts = gtk_tree_store_new(1, G_TYPE_STRING);
        GtkCellRenderer *cr = gtk_cell_renderer_text_new();
        g_object_set( G_OBJECT(cr), "font", "Arial bold 8", NULL );
        gtk_tree_store_clear( GTK_TREE_STORE(ts) );

        for (ii=nNamedMaterials-1; ii>=0; ii--){
            gtk_tree_store_insert( GTK_TREE_STORE(ts), &iter, NULL, 0 );
            gtk_tree_store_set(GTK_TREE_STORE(ts), &iter, 0, NamedMaterials[ii].Name, -1);
        }
        gtk_tree_store_insert( GTK_TREE_STORE(ts), &iter, NULL, 0 );
        gtk_tree_store_set(GTK_TREE_STORE(ts), &iter, 0, "Custom:", -1);

        gInfo->FieldLineMaterialButton[i] = gtk_combo_box_new_with_model( GTK_TREE_MODEL(ts) );
        gtk_widget_set_size_request( gInfo->FieldLineMaterialButton[i], 30, 20);
        gtk_cell_layout_pack_start( GTK_CELL_LAYOUT(gInfo->FieldLineMaterialButton[i]), cr, FALSE);
        gtk_cell_layout_set_attributes( GTK_CELL_LAYOUT(gInfo->FieldLineMaterialButton[i]), cr, "text", 0, NULL);

        gtk_combo_box_set_active( GTK_COMBO_BOX(gInfo->FieldLineMaterialButton[i]), 0);
        gtk_widget_show( gInfo->FieldLineMaterialButton[i] );
        gtk_table_attach( GTK_TABLE(table1), gInfo->FieldLineMaterialButton[i], col, col+1, 3+i, 4+i, (GtkAttachOptions)(GTK_FILL), (GtkAttachOptions)(GTK_SHRINK), 0, 0 );
        gInfo->FieldLineMaterialButtonHandler[i] = g_signal_connect( G_OBJECT( gInfo->FieldLineMaterialButton[i] ), "changed", G_CALLBACK( ChangeMaterial ), GINT_TO_POINTER(i) );
        ++col;


        color.red = gInfo->FieldLineMaterial[i].diffuse[0]*65535;
        color.green = gInfo->FieldLineMaterial[i].diffuse[1]*65535;
        color.blue = gInfo->FieldLineMaterial[i].diffuse[2]*65535;
        gInfo->FieldLineDiffuseColorButton[i] = gtk_color_button_new_with_color( &color ); gtk_widget_show( gInfo->FieldLineDiffuseColorButton[i] );
        gtk_widget_set_size_request( gInfo->FieldLineDiffuseColorButton[i], 30, 20);
        gtk_table_attach( GTK_TABLE(table1), gInfo->FieldLineDiffuseColorButton[i], col, col+1, 3+i, 4+i, (GtkAttachOptions)(GTK_FILL), (GtkAttachOptions)(GTK_SHRINK), 0, 0 );
        gtk_color_button_set_use_alpha( GTK_COLOR_BUTTON(gInfo->FieldLineDiffuseColorButton[i]), TRUE );
        g_signal_connect( G_OBJECT( gInfo->FieldLineDiffuseColorButton[i] ), "color-set", G_CALLBACK( ChangeMaterialColor ), GINT_TO_POINTER(i*3+0) );
        ++col;

        color.red = gInfo->FieldLineMaterial[i].ambient[0]*65535;
        color.green = gInfo->FieldLineMaterial[i].ambient[1]*65535;
        color.blue = gInfo->FieldLineMaterial[i].ambient[2]*65535;
        gInfo->FieldLineAmbientColorButton[i] = gtk_color_button_new_with_color( &color ); gtk_widget_show( gInfo->FieldLineAmbientColorButton[i] );
        gtk_widget_set_size_request( gInfo->FieldLineAmbientColorButton[i], 30, 20);
        gtk_table_attach( GTK_TABLE(table1), gInfo->FieldLineAmbientColorButton[i], col, col+1, 3+i, 4+i, (GtkAttachOptions)(GTK_FILL), (GtkAttachOptions)(GTK_SHRINK), 0, 0 );
        gtk_color_button_set_use_alpha( GTK_COLOR_BUTTON(gInfo->FieldLineAmbientColorButton[i]), TRUE );
        g_signal_connect( G_OBJECT( gInfo->FieldLineAmbientColorButton[i] ), "color-set", G_CALLBACK( ChangeMaterialColor ), GINT_TO_POINTER(i*3+1) );
        ++col;

        color.red = gInfo->FieldLineMaterial[i].specular[0]*65535;
        color.green = gInfo->FieldLineMaterial[i].specular[1]*65535;
        color.blue = gInfo->FieldLineMaterial[i].specular[2]*65535;
        gInfo->FieldLineSpecularColorButton[i] = gtk_color_button_new_with_color( &color ); gtk_widget_show( gInfo->FieldLineSpecularColorButton[i] );
        gtk_widget_set_size_request( gInfo->FieldLineSpecularColorButton[i], 30, 20);
        gtk_table_attach( GTK_TABLE(table1), gInfo->FieldLineSpecularColorButton[i], col, col+1, 3+i, 4+i, (GtkAttachOptions)(GTK_FILL), (GtkAttachOptions)(GTK_SHRINK), 0, 0 );
        gtk_color_button_set_use_alpha( GTK_COLOR_BUTTON(gInfo->FieldLineSpecularColorButton[i]), TRUE );
        g_signal_connect( G_OBJECT( gInfo->FieldLineSpecularColorButton[i] ), "color-set", G_CALLBACK( ChangeMaterialColor ), GINT_TO_POINTER(i*3+2) );
        ++col;

        spinbutton1_adj = gtk_adjustment_new( 0.15, 0, 1, 0.01, 0.1, 0 );
        gInfo->FieldLineShininessButton[i]  = gtk_spin_button_new( GTK_ADJUSTMENT(spinbutton1_adj), 1, 2); gtk_widget_show( gInfo->FieldLineShininessButton[i]  );
        gtk_widget_set_size_request( gInfo->FieldLineShininessButton[i], 30, 20);
        gtk_table_attach( GTK_TABLE(table1), gInfo->FieldLineShininessButton[i] , col, col+1, 3+i, 4+i, (GtkAttachOptions)(GTK_FILL), (GtkAttachOptions)(GTK_SHRINK), 0, 0 );
        gtk_spin_button_set_numeric( GTK_SPIN_BUTTON(gInfo->FieldLineShininessButton[i] ), TRUE );
        gInfo->FieldLineShininessButtonHandler[i] = g_signal_connect( G_OBJECT( gInfo->FieldLineShininessButton[i] ), "value_changed", G_CALLBACK( ChangeMaterialShininess ), GINT_TO_POINTER(i) );
        ++col;




    }







    /*****************************************************************************************************
     *
     *                                  Drift Shell Notebook Page
     *
     *****************************************************************************************************/
    vbox2 = gtk_vbox_new (FALSE, 0); gtk_widget_show (vbox2);
    gtk_container_set_border_width (GTK_CONTAINER (vbox2), 20);

    // Page title
    label = gtk_label_new (_("<b><span size=\"small\">Drift Shells</span></b>"));
    gtk_widget_show (label);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_CENTER);
    gtk_label_set_use_markup( GTK_LABEL(label), TRUE );
    gtk_notebook_append_page( GTK_NOTEBOOK(notebook), vbox2, label );

    // Table for widgets
    table1 = gtk_table_new (40, 10, FALSE);
    gtk_widget_show (table1);
    gtk_box_pack_start( GTK_BOX (vbox2), table1, TRUE, TRUE, 15 );
    gtk_table_set_row_spacings (GTK_TABLE (table1), 0);
    gtk_table_set_col_spacings (GTK_TABLE (table1), 10);


    col = 0;
    // Pitch angle column
    label = gtk_label_new (_("<b>Pitch\nAngle</b>")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, col, col+1, 0, 2, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions)(GTK_FILL), 0, 0);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_CENTER);
    gtk_label_set_use_markup( GTK_LABEL(label), TRUE );
    ++col;

    // vertical separator
    vseparator = gtk_vseparator_new (); gtk_widget_show (vseparator);
    gtk_table_attach (GTK_TABLE (table1), vseparator, col, col+1, 0, 13, (GtkAttachOptions)(GTK_FILL), (GtkAttachOptions)(GTK_FILL), 0, 0);
    ++col;


    /*
     *  Column Headers for Drift Shells
     */
    label = gtk_label_new (_("<b>Show\nDrift\nShells</b>")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, col, col+1, 0, 2, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_CENTER);
    gtk_label_set_use_markup( GTK_LABEL(label), TRUE );
    ++col;

    label = gtk_label_new (_("<b><span size=\"large\">Drift Shell Material</span></b>")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, 3, 8, 0, 1, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_CENTER);
    gtk_label_set_use_markup( GTK_LABEL(label), TRUE );

    label = gtk_label_new (_("<b><span size=\"small\">Material</span></b>")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, col, col+1, 1, 2, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_CENTER);
    gtk_label_set_use_markup( GTK_LABEL(label), TRUE );
    ++col;

    label = gtk_label_new (_("<b><span size=\"small\">Diffuse</span></b>")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, col, col+1, 1, 2, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_CENTER);
    gtk_label_set_use_markup( GTK_LABEL(label), TRUE );
    ++col;

    label = gtk_label_new (_("<b><span size=\"small\">Ambient</span></b>")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, col, col+1, 1, 2, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_CENTER);
    gtk_label_set_use_markup( GTK_LABEL(label), TRUE );
    ++col;

    label = gtk_label_new (_("<b><span size=\"small\">Specular</span></b>")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, col, col+1, 1, 2, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_CENTER);
    gtk_label_set_use_markup( GTK_LABEL(label), TRUE );
    ++col;

    label = gtk_label_new (_("<b><span size=\"small\">Shininess</span></b>")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, col, col+1, 1, 2, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_CENTER);
    gtk_label_set_use_markup( GTK_LABEL(label), TRUE );
    ++col;

    // Separator
    hseparator = gtk_hseparator_new (); gtk_widget_show (hseparator);
    gtk_table_attach (GTK_TABLE (table1), hseparator, 0, 8, 2, 3, (GtkAttachOptions) (GTK_EXPAND | GTK_FILL), (GtkAttachOptions) (GTK_FILL), 0, 0);


    /*
     *  Fill in rows
     */
    for (i=0; i<=ObjInfo->MagEphemInfo->nAlpha; i++){

        col = 0;

        if (i==ObjInfo->MagEphemInfo->nAlpha){
            sprintf( Str, "<small><b>All</b></small>" );
        } else {
            sprintf( Str, "<small><b>%g\u00b0</b></small>", ObjInfo->MagEphemInfo->Alpha[i] );
        }
        label = gtk_label_new( Str ); gtk_widget_show( label );
        gtk_table_attach( GTK_TABLE(table1), label, col, col+1, 3+i, 4+i, (GtkAttachOptions)(GTK_FILL), (GtkAttachOptions)(0), 0, 0 );
        gtk_label_set_use_markup( GTK_LABEL(label), TRUE );
        gtk_label_set_justify( GTK_LABEL(label), GTK_JUSTIFY_CENTER);
        ++col;
        ++col;



        /*
         *  Drift Shells
         */
        gInfo->DriftShellShowPitchAngleButton[i] = gtk_check_button_new();
        gtk_widget_show( gInfo->DriftShellShowPitchAngleButton[i] );
        gtk_widget_set_size_request( gInfo->DriftShellShowPitchAngleButton[i], 30, 20);
        if (i==0) gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON( gInfo->DriftShellShowPitchAngleButton[i] ), ShowAllPitchAngles2 );
        else      gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON( gInfo->DriftShellShowPitchAngleButton[i] ), ShowPitchAngle2[i] );
        gtk_table_attach( GTK_TABLE(table1), gInfo->DriftShellShowPitchAngleButton[i], col, col+1, 3+i, 4+i, (GtkAttachOptions)(GTK_FILL), (GtkAttachOptions)(GTK_SHRINK), 0, 0 );
        gInfo->DriftShellShowPitchAngleButtonHandler[i] = g_signal_connect( G_OBJECT( gInfo->DriftShellShowPitchAngleButton[i] ), "toggled", G_CALLBACK( SelectPitchAngles2 ), GINT_TO_POINTER(100+i) );
        ++col;


        GtkTreeIter     iter;
        GtkTreeStore    *ts = gtk_tree_store_new(1, G_TYPE_STRING);
        GtkCellRenderer *cr = gtk_cell_renderer_text_new();
        g_object_set( G_OBJECT(cr), "font", "Arial bold 8", NULL );
        gtk_tree_store_clear( GTK_TREE_STORE(ts) );

        for (ii=nNamedMaterials-1; ii>=0; ii--){
            gtk_tree_store_insert( GTK_TREE_STORE(ts), &iter, NULL, 0 );
            gtk_tree_store_set(GTK_TREE_STORE(ts), &iter, 0, NamedMaterials[ii].Name, -1);
        }
        gtk_tree_store_insert( GTK_TREE_STORE(ts), &iter, NULL, 0 );
        gtk_tree_store_set(GTK_TREE_STORE(ts), &iter, 0, "Custom:", -1);

        gInfo->DriftShellMaterialButton[i] = gtk_combo_box_new_with_model( GTK_TREE_MODEL(ts) );
        gtk_widget_set_size_request( gInfo->DriftShellMaterialButton[i], 30, 20);
        gtk_cell_layout_pack_start( GTK_CELL_LAYOUT(gInfo->DriftShellMaterialButton[i]), cr, FALSE);
        gtk_cell_layout_set_attributes( GTK_CELL_LAYOUT(gInfo->DriftShellMaterialButton[i]), cr, "text", 0, NULL);

        gtk_combo_box_set_active( GTK_COMBO_BOX(gInfo->DriftShellMaterialButton[i]), 0);
        gtk_widget_show( gInfo->DriftShellMaterialButton[i] );
        gtk_table_attach( GTK_TABLE(table1), gInfo->DriftShellMaterialButton[i], col, col+1, 3+i, 4+i, (GtkAttachOptions)(GTK_FILL), (GtkAttachOptions)(GTK_SHRINK), 0, 0 );
        gInfo->DriftShellMaterialButtonHandler[i] = g_signal_connect( G_OBJECT( gInfo->DriftShellMaterialButton[i] ), "changed", G_CALLBACK( ChangeMaterial ), GINT_TO_POINTER(100+i) );
        ++col;











        color.red = gInfo->DriftShellMaterial[i].diffuse[0]*65535;
        color.green = gInfo->DriftShellMaterial[i].diffuse[1]*65535;
        color.blue = gInfo->DriftShellMaterial[i].diffuse[2]*65535;
        gInfo->DriftShellDiffuseColorButton[i] = gtk_color_button_new_with_color( &color ); gtk_widget_show( gInfo->DriftShellDiffuseColorButton[i] );
        gtk_widget_set_size_request( gInfo->DriftShellDiffuseColorButton[i], 30, 20);
        gtk_table_attach( GTK_TABLE(table1), gInfo->DriftShellDiffuseColorButton[i], col, col+1, 3+i, 4+i, (GtkAttachOptions)(GTK_FILL), (GtkAttachOptions)(GTK_SHRINK), 0, 0 );
        gtk_color_button_set_use_alpha( GTK_COLOR_BUTTON(gInfo->DriftShellDiffuseColorButton[i]), TRUE );
        gtk_color_button_set_alpha( GTK_COLOR_BUTTON(gInfo->DriftShellDiffuseColorButton[i]), gInfo->DriftShellMaterial[i].diffuse[3]*65535 );
        g_signal_connect( G_OBJECT( gInfo->DriftShellDiffuseColorButton[i] ), "color-set", G_CALLBACK( ChangeMaterialColor ), GINT_TO_POINTER(100+(i)*3+0) );
        ++col;

        color.red = gInfo->DriftShellMaterial[i].ambient[0]*65535;
        color.green = gInfo->DriftShellMaterial[i].ambient[1]*65535;
        color.blue = gInfo->DriftShellMaterial[i].ambient[2]*65535;
        gInfo->DriftShellAmbientColorButton[i] = gtk_color_button_new_with_color( &color ); gtk_widget_show( gInfo->DriftShellAmbientColorButton[i] );
        gtk_widget_set_size_request( gInfo->DriftShellAmbientColorButton[i], 30, 20);
        gtk_table_attach( GTK_TABLE(table1), gInfo->DriftShellAmbientColorButton[i], col, col+1, 3+i, 4+i, (GtkAttachOptions)(GTK_FILL), (GtkAttachOptions)(GTK_SHRINK), 0, 0 );
        gtk_color_button_set_use_alpha( GTK_COLOR_BUTTON(gInfo->DriftShellAmbientColorButton[i]), TRUE );
        gtk_color_button_set_alpha( GTK_COLOR_BUTTON(gInfo->DriftShellAmbientColorButton[i]), gInfo->DriftShellMaterial[i].ambient[3]*65535 );
        g_signal_connect( G_OBJECT( gInfo->DriftShellAmbientColorButton[i] ), "color-set", G_CALLBACK( ChangeMaterialColor ), GINT_TO_POINTER(100+(i)*3+1) );
        ++col;

        color.red = gInfo->DriftShellMaterial[i].specular[0]*65535;
        color.green = gInfo->DriftShellMaterial[i].specular[1]*65535;
        color.blue = gInfo->DriftShellMaterial[i].specular[2]*65535;
        gInfo->DriftShellSpecularColorButton[i] = gtk_color_button_new_with_color( &color ); gtk_widget_show( gInfo->DriftShellSpecularColorButton[i] );
        gtk_widget_set_size_request( gInfo->DriftShellSpecularColorButton[i], 30, 20);
        gtk_table_attach( GTK_TABLE(table1), gInfo->DriftShellSpecularColorButton[i], col, col+1, 3+i, 4+i, (GtkAttachOptions)(GTK_FILL), (GtkAttachOptions)(GTK_SHRINK), 0, 0 );
        gtk_color_button_set_use_alpha( GTK_COLOR_BUTTON(gInfo->DriftShellSpecularColorButton[i]), TRUE );
        gtk_color_button_set_alpha( GTK_COLOR_BUTTON(gInfo->DriftShellSpecularColorButton[i]), gInfo->DriftShellMaterial[i].specular[3]*65535 );
        g_signal_connect( G_OBJECT( gInfo->DriftShellSpecularColorButton[i] ), "color-set", G_CALLBACK( ChangeMaterialColor ), GINT_TO_POINTER(100+(i)*3+2) );
        ++col;

        spinbutton1_adj = gtk_adjustment_new( 0.15, 0, 1, 0.01, 0.1, 0 );
        gInfo->DriftShellShininessButton[i]  = gtk_spin_button_new( GTK_ADJUSTMENT(spinbutton1_adj), 1, 2); gtk_widget_show( gInfo->DriftShellShininessButton[i]  );
        gtk_widget_set_size_request( gInfo->DriftShellShininessButton[i], 30, 20);
        gtk_table_attach( GTK_TABLE(table1), gInfo->DriftShellShininessButton[i] , col, col+1, 3+i, 4+i, (GtkAttachOptions)(GTK_FILL), (GtkAttachOptions)(GTK_SHRINK), 0, 0 );
        gtk_spin_button_set_numeric( GTK_SPIN_BUTTON(gInfo->DriftShellShininessButton[i] ), TRUE );
        gInfo->DriftShellShininessButtonHandler[i] = g_signal_connect( G_OBJECT( gInfo->DriftShellShininessButton[i] ), "value_changed", G_CALLBACK( ChangeMaterialShininess ), GINT_TO_POINTER(100+i) );
        ++col;


    }




    /*******************
     * Satellites Page *
     *******************/
    //SatSelectorInfo = New_SpaceObjects();
    //InitSatSelectorInfo();

    vbox2 = CreateSatSelector( nSatTypes, SatTypes );
PUKE_SATSEL_VBOX = vbox2;
    gtk_container_set_border_width (GTK_CONTAINER (vbox2), 20);

    // Page title
    label = gtk_label_new (_("<b><span size=\"small\">Satellites</span></b>")); gtk_widget_show (label);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_CENTER);
    gtk_label_set_use_markup( GTK_LABEL(label), TRUE );
    gtk_notebook_append_page( GTK_NOTEBOOK(notebook), vbox2, label );






    /**********************
     * Iridium Flare Page *
     **********************/
    vbox2 = gtk_vbox_new (FALSE, 0); gtk_widget_show (vbox2);
    gtk_container_set_border_width (GTK_CONTAINER (vbox2), 20);

    // Page title
    label = gtk_label_new (_("<b><span size=\"small\">Iridium Flares</span></b>")); gtk_widget_show (label);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_CENTER);
    gtk_label_set_use_markup( GTK_LABEL(label), TRUE );
    gtk_notebook_append_page( GTK_NOTEBOOK(notebook), vbox2, label );

    // top frame
    frame = gtk_frame_new (NULL); gtk_widget_show (frame);
    gtk_box_pack_start (GTK_BOX (vbox2), frame, FALSE, FALSE, 0);
    gtk_container_set_border_width (GTK_CONTAINER (frame), 10);
    gtk_frame_set_shadow_type (GTK_FRAME (frame), GTK_SHADOW_IN);

    label = gtk_label_new (_("<b>  Lines Drawn When Iridium Satellites are Visible </b>")); gtk_widget_show (label);
    gtk_frame_set_label_widget (GTK_FRAME (frame), label);
    gtk_label_set_use_markup (GTK_LABEL (label), TRUE);


    alignment = gtk_alignment_new (0.5, 0.5, 1, 1); gtk_widget_show (alignment);
    gtk_container_add (GTK_CONTAINER (frame), alignment);
    gtk_alignment_set_padding (GTK_ALIGNMENT (alignment), 0, 0, 12, 0);

    hbox = gtk_hbox_new (FALSE, 30); gtk_widget_show (hbox);
    gtk_container_add (GTK_CONTAINER (alignment), hbox);
    gtk_container_set_border_width (GTK_CONTAINER (hbox), 10);

    checkbutton = gtk_check_button_new_with_mnemonic (_("Sat to Ground Ray")); gtk_widget_show (checkbutton);
    gtk_box_pack_start( GTK_BOX(hbox), checkbutton, FALSE, FALSE, 0 );
    gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(checkbutton), ShowIridiumSatToGround );
    g_signal_connect( G_OBJECT( checkbutton ), "toggled", G_CALLBACK( ToggleIridiumLines ), GINT_TO_POINTER( 0 ) );

    checkbutton = gtk_check_button_new_with_mnemonic (_("Sun to Sat Ray")); gtk_widget_show (checkbutton);
    gtk_box_pack_start( GTK_BOX(hbox), checkbutton, FALSE, FALSE, 0 );
    gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(checkbutton), ShowIridiumSatToSun );
    g_signal_connect( G_OBJECT( checkbutton ), "toggled", G_CALLBACK( ToggleIridiumLines ), GINT_TO_POINTER( 1 ) );



    // bottom frame
    frame = gtk_frame_new (NULL); gtk_widget_show (frame);
    gtk_box_pack_start (GTK_BOX (vbox2), frame, TRUE, TRUE, 0);
    gtk_container_set_border_width (GTK_CONTAINER (frame), 10);
    gtk_frame_set_shadow_type (GTK_FRAME (frame), GTK_SHADOW_IN);

    label = gtk_label_new (_("<b>  Predict Iridium Flares at Location  </b>")); gtk_widget_show (label);
    gtk_frame_set_label_widget (GTK_FRAME (frame), label);
    gtk_label_set_use_markup (GTK_LABEL (label), TRUE);

    alignment = gtk_alignment_new (0.5, 0.5, 1, 1); gtk_widget_show (alignment);
    gtk_container_add (GTK_CONTAINER (frame), alignment);
    gtk_container_set_border_width (GTK_CONTAINER (alignment), 20);
    gtk_alignment_set_padding (GTK_ALIGNMENT (alignment), 0, 0, 12, 0);

    vbox2 = gtk_vbox_new (FALSE, 22); gtk_widget_show (vbox2);
    gtk_container_add (GTK_CONTAINER (alignment), vbox2);

    table1 = gtk_table_new (3, 4, FALSE); gtk_widget_show (table1);
    gtk_box_pack_start (GTK_BOX (vbox2), table1, FALSE, FALSE, 0);
    gtk_table_set_row_spacings (GTK_TABLE (table1), 10);
    gtk_table_set_col_spacings (GTK_TABLE (table1), 30);



    // Latitude
    hbox = gtk_hbox_new (FALSE, 0); gtk_widget_show (hbox);
    gtk_table_attach (GTK_TABLE (table1), hbox, 0, 1, 0, 1, (GtkAttachOptions) (GTK_EXPAND | GTK_FILL), (GtkAttachOptions) (GTK_EXPAND | GTK_FILL), 0, 0);

    label = gtk_label_new (_("Latitude (\302\260):")); gtk_widget_show (label);
    gtk_box_pack_start (GTK_BOX (hbox), label, FALSE, FALSE, 0);

    spinbutton1_adj = gtk_adjustment_new (35.8880004883, -90, 90, 0.0010000000475, 0.10000000149, 0);
    spinbutton = gtk_spin_button_new (GTK_ADJUSTMENT (spinbutton1_adj), 1, 4); gtk_widget_show (spinbutton);
    gtk_box_pack_start (GTK_BOX (hbox), spinbutton, TRUE, TRUE, 0);
    gtk_spin_button_set_numeric (GTK_SPIN_BUTTON (spinbutton), TRUE);

    // Longitude
    hbox = gtk_hbox_new (FALSE, 0); gtk_widget_show (hbox);
    gtk_table_attach (GTK_TABLE (table1), hbox, 1, 2, 0, 1, (GtkAttachOptions) (GTK_EXPAND | GTK_FILL), (GtkAttachOptions) (GTK_FILL), 0, 0);

    label = gtk_label_new (_("Longitude (\302\260):")); gtk_widget_show (label);
    gtk_box_pack_start (GTK_BOX (hbox), label, FALSE, FALSE, 0);

    spinbutton1_adj = gtk_adjustment_new (106.305999756, 0, 360, 0.0010000000475, 0.10000000149, 0);
    spinbutton = gtk_spin_button_new (GTK_ADJUSTMENT (spinbutton1_adj), 1, 4); gtk_widget_show (spinbutton);
    gtk_box_pack_start (GTK_BOX (hbox), spinbutton, TRUE, TRUE, 0);
    gtk_spin_button_set_numeric (GTK_SPIN_BUTTON (spinbutton), TRUE);

    // Altitude
    hbox = gtk_hbox_new (FALSE, 0); gtk_widget_show (hbox);
    gtk_table_attach (GTK_TABLE (table1), hbox, 2, 3, 0, 1, (GtkAttachOptions) (GTK_EXPAND | GTK_FILL), (GtkAttachOptions) (GTK_FILL), 0, 0);

    label = gtk_label_new (_("Altitude (m):")); gtk_widget_show (label);
    gtk_box_pack_start (GTK_BOX (hbox), label, FALSE, FALSE, 0);

    spinbutton1_adj = gtk_adjustment_new (7500, 0, 30000, 1, 10, 0);
    spinbutton = gtk_spin_button_new (GTK_ADJUSTMENT (spinbutton1_adj), 1, 0); gtk_widget_show (spinbutton);
    gtk_box_pack_start (GTK_BOX (hbox), spinbutton, TRUE, TRUE, 0);
    gtk_spin_button_set_numeric (GTK_SPIN_BUTTON (spinbutton), TRUE);

    // Find City
    button = gtk_button_new_with_mnemonic (_("Find City")); gtk_widget_show (button);
    gtk_table_attach (GTK_TABLE (table1), button, 3, 4, 0, 1, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);



    // Magnitude
    hbox = gtk_hbox_new (FALSE, 10); gtk_widget_show (hbox);
    gtk_table_attach (GTK_TABLE (table1), hbox, 0, 1, 1, 2, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (GTK_EXPAND | GTK_FILL), 0, 0);

    label = gtk_label_new (_("Find Magnitudes\nBrighter Than:")); gtk_widget_show (label);
    gtk_box_pack_start (GTK_BOX (hbox), label, FALSE, FALSE, 0);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_RIGHT);
    gtk_misc_set_alignment (GTK_MISC (label), 1, 0.5);

    spinbutton1_adj = gtk_adjustment_new (-6, -10, 10, 0.10000000149, 0.5, 0);
    spinbutton = gtk_spin_button_new (GTK_ADJUSTMENT (spinbutton1_adj), 1, 1); gtk_widget_show (spinbutton);
    gtk_box_pack_start (GTK_BOX (hbox), spinbutton, TRUE, TRUE, 0);
    gtk_spin_button_set_numeric (GTK_SPIN_BUTTON (spinbutton), TRUE);


    // Mirror Angle
    hbox = gtk_hbox_new (FALSE, 10); gtk_widget_show (hbox);
    gtk_table_attach (GTK_TABLE (table1), hbox, 1, 2, 1, 2, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (GTK_FILL), 0, 0);

    label = gtk_label_new (_("or Mirror Angle\nLess Than (\302\260):")); gtk_widget_show (label);
    gtk_box_pack_start (GTK_BOX (hbox), label, FALSE, FALSE, 0);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_RIGHT);
    gtk_misc_set_alignment (GTK_MISC (label), 1, 0.5);

    spinbutton1_adj = gtk_adjustment_new (0.8, 0, 10, 0.00999999977648, 0.10000000149, 0);
    spinbutton = gtk_spin_button_new (GTK_ADJUSTMENT (spinbutton1_adj), 1, 4); gtk_widget_show (spinbutton);
    gtk_box_pack_start (GTK_BOX (hbox), spinbutton, TRUE, TRUE, 0);
    gtk_spin_button_set_numeric (GTK_SPIN_BUTTON (spinbutton), TRUE);


    // Search
    hbox = gtk_hbox_new (TRUE, 0); gtk_widget_show (hbox);
    gtk_table_attach (GTK_TABLE (table1), hbox, 0, 4, 2, 3, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (GTK_EXPAND | GTK_FILL), 0, 0);

    label = gtk_label_new (_("Search:")); gtk_widget_show (label);
    gtk_box_pack_start (GTK_BOX (hbox), label, FALSE, FALSE, 0);

    button = gtk_button_new_with_mnemonic (_("Prev 7 Days")); gtk_widget_show (button);
    gtk_box_pack_start (GTK_BOX (hbox), button, FALSE, FALSE, 0);
    g_signal_connect( G_OBJECT( button ), "clicked", G_CALLBACK( PredictIridiumFlares ), GINT_TO_POINTER( -168 ) );

    button = gtk_button_new_with_mnemonic (_("Prev 48 Hrs")); gtk_widget_show (button);
    gtk_box_pack_start (GTK_BOX (hbox), button, FALSE, FALSE, 0);
    g_signal_connect( G_OBJECT( button ), "clicked", G_CALLBACK( PredictIridiumFlares ), GINT_TO_POINTER( -48 ) );

    button = gtk_button_new_with_mnemonic (_("Prev 24 Hrs")); gtk_widget_show (button);
    gtk_box_pack_start (GTK_BOX (hbox), button, FALSE, FALSE, 0);
    g_signal_connect( G_OBJECT( button ), "clicked", G_CALLBACK( PredictIridiumFlares ), GINT_TO_POINTER( -24 ) );

    button = gtk_button_new_with_mnemonic (_("Next 24 Hrs")); gtk_widget_show (button);
    gtk_box_pack_start (GTK_BOX (hbox), button, FALSE, FALSE, 0);
    g_signal_connect( G_OBJECT( button ), "clicked", G_CALLBACK( PredictIridiumFlares ), GINT_TO_POINTER( 24 ) );

    button = gtk_button_new_with_mnemonic (_("Next 48 Hrs")); gtk_widget_show (button);
    gtk_box_pack_start (GTK_BOX (hbox), button, FALSE, FALSE, 0);
    g_signal_connect( G_OBJECT( button ), "clicked", G_CALLBACK( PredictIridiumFlares ), GINT_TO_POINTER( 48 ) );

    button = gtk_button_new_with_mnemonic (_("Next 7 Days")); gtk_widget_show (button);
    gtk_box_pack_start (GTK_BOX (hbox), button, FALSE, FALSE, 0);
    g_signal_connect( G_OBJECT( button ), "clicked", G_CALLBACK( PredictIridiumFlares ), GINT_TO_POINTER( 168 ) );


    // List View
    scrolledwindow = gtk_scrolled_window_new (NULL, NULL); gtk_widget_show (scrolledwindow);
    gtk_box_pack_start (GTK_BOX (vbox2), scrolledwindow, TRUE, TRUE, 0);
    gtk_scrolled_window_set_shadow_type (GTK_SCROLLED_WINDOW (scrolledwindow), GTK_SHADOW_IN);

    /* create list store model */
    model = gtk_list_store_new( IRIDIUM_FLARE_NUM_COLUMNS,
                                            G_TYPE_STRING,
                                            G_TYPE_STRING,
                                            G_TYPE_DOUBLE,
                                            G_TYPE_DOUBLE,
                                            G_TYPE_DOUBLE,
                                            G_TYPE_DOUBLE,
                                            G_TYPE_DOUBLE,
                                            G_TYPE_STRING);

    /* create tree view */
    treeview = gtk_tree_view_new_with_model (model);
    gtk_tree_view_set_rules_hint( GTK_TREE_VIEW (treeview), TRUE );
    gtk_tree_view_set_search_column( GTK_TREE_VIEW (treeview), IRIDIUM_FLARE_COLUMN_DATE );
    gtk_tree_view_set_enable_tree_lines(GTK_TREE_VIEW (treeview), TRUE);

    gtk_tree_view_set_enable_search( GTK_TREE_VIEW (treeview), TRUE );
    gtk_tree_view_set_show_expanders( GTK_TREE_VIEW (treeview), TRUE );

    gtk_container_add( GTK_CONTAINER(scrolledwindow), treeview );


    /* column for Date */
    renderer = gtk_cell_renderer_text_new ();
    column = gtk_tree_view_column_new_with_attributes( "Date", renderer, "text", IRIDIUM_FLARE_COLUMN_DATE, NULL );
    gtk_tree_view_column_set_alignment( column, 0.5 );
    gtk_tree_view_column_set_resizable( column, TRUE );
    gtk_tree_view_column_set_sort_column_id( column, IRIDIUM_FLARE_COLUMN_DATE );
    gtk_tree_view_append_column( GTK_TREE_VIEW(treeview), column );

    /* column for Time */
    renderer = gtk_cell_renderer_text_new ();
    column = gtk_tree_view_column_new_with_attributes( "Time", renderer, "text", IRIDIUM_FLARE_COLUMN_TIME, NULL );
    gtk_tree_view_column_set_alignment( column, 0.5 );
    gtk_tree_view_column_set_resizable( column, TRUE );
    gtk_tree_view_column_set_sort_column_id( column, IRIDIUM_FLARE_COLUMN_TIME );
    gtk_tree_view_append_column( GTK_TREE_VIEW(treeview), column );

    /* column for Magnitude */
    renderer = gtk_cell_renderer_text_new ();
    column = gtk_tree_view_column_new_with_attributes( "Mag.", renderer, "text", IRIDIUM_FLARE_COLUMN_MAGNITUDE, NULL );
    gtk_tree_view_column_set_alignment( column, 0.5 );
    gtk_tree_view_column_set_resizable( column, TRUE );
    gtk_tree_view_column_set_sort_column_id( column, IRIDIUM_FLARE_COLUMN_MAGNITUDE );
    gtk_tree_view_append_column( GTK_TREE_VIEW(treeview), column );

    /* column for Alt */
    renderer = gtk_cell_renderer_text_new ();
    column = gtk_tree_view_column_new_with_attributes( "Alt", renderer, "text", IRIDIUM_FLARE_COLUMN_ALT, NULL );
    gtk_tree_view_column_set_alignment( column, 0.5 );
    gtk_tree_view_column_set_resizable( column, TRUE );
    gtk_tree_view_column_set_sort_column_id( column, IRIDIUM_FLARE_COLUMN_ALT );
    gtk_tree_view_append_column( GTK_TREE_VIEW(treeview), column );

    /* column for Az */
    renderer = gtk_cell_renderer_text_new ();
    column = gtk_tree_view_column_new_with_attributes( "Az", renderer, "text", IRIDIUM_FLARE_COLUMN_AZ, NULL );
    gtk_tree_view_column_set_alignment( column, 0.5 );
    gtk_tree_view_column_set_resizable( column, TRUE );
    gtk_tree_view_column_set_sort_column_id( column, IRIDIUM_FLARE_COLUMN_AZ );
    gtk_tree_view_append_column( GTK_TREE_VIEW(treeview), column );

    /* column for Dist to Center */
    renderer = gtk_cell_renderer_text_new ();
    column = gtk_tree_view_column_new_with_attributes( "Center Dist.", renderer, "text", IRIDIUM_FLARE_COLUMN_DIST_TO_CENTER, NULL );
    gtk_tree_view_column_set_alignment( column, 0.5 );
    gtk_tree_view_column_set_resizable( column, TRUE );
    gtk_tree_view_column_set_sort_column_id( column, IRIDIUM_FLARE_COLUMN_DIST_TO_CENTER );
    gtk_tree_view_append_column( GTK_TREE_VIEW(treeview), column );

    /* column for Magnitude at Center */
    renderer = gtk_cell_renderer_text_new ();
    column = gtk_tree_view_column_new_with_attributes( "Center Mag.", renderer, "text", IRIDIUM_FLARE_COLUMN_MAGNITUDE_AT_CENTER, NULL );
    gtk_tree_view_column_set_alignment( column, 0.5 );
    gtk_tree_view_column_set_resizable( column, TRUE );
    gtk_tree_view_column_set_sort_column_id( column, IRIDIUM_FLARE_COLUMN_MAGNITUDE_AT_CENTER );
    gtk_tree_view_append_column( GTK_TREE_VIEW(treeview), column );

    /* column for Dist to Center */
    renderer = gtk_cell_renderer_text_new ();
    column = gtk_tree_view_column_new_with_attributes( "Spacecraft", renderer, "text", IRIDIUM_FLARE_COLUMN_SATELLITE, NULL );
    gtk_tree_view_column_set_alignment( column, 0.5 );
    gtk_tree_view_column_set_resizable( column, TRUE );
    gtk_tree_view_column_set_sort_column_id( column, IRIDIUM_FLARE_COLUMN_SATELLITE );
    gtk_tree_view_append_column( GTK_TREE_VIEW(treeview), column );

    gtk_widget_show_all( scrolledwindow );
    gtk_widget_show( treeview );











    /*******************
     * Atmosphere Page *
     *******************/
    vbox2 = gtk_vbox_new (FALSE, 0); gtk_widget_show (vbox2);
    gtk_container_set_border_width (GTK_CONTAINER (vbox2), 20);

    // Page title
    label = gtk_label_new (_("<b><span size=\"small\">Atmosphere</span></b>")); gtk_widget_show (label);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_CENTER);
    gtk_label_set_use_markup( GTK_LABEL(label), TRUE );
    gtk_notebook_append_page( GTK_NOTEBOOK(notebook), vbox2, label );


    table1 = gtk_table_new (7, 2, FALSE); gtk_widget_show (table1);
    gtk_box_pack_start (GTK_BOX (vbox2), table1, TRUE, TRUE, 15);
    gtk_container_set_border_width (GTK_CONTAINER (table1), 20);
    gtk_table_set_row_spacings (GTK_TABLE (table1), 20);
    gtk_table_set_col_spacings (GTK_TABLE (table1), 10);

    /*
     * Row for Toggle button
     * The rows on this page are two cols: col 1 is label col two holds hbox (which holds othwer hboxes)
     */
    label = gtk_label_new (_("<b>Show\nAtmosphere:</b>")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, 0, 1, 0, 1, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);
    gtk_label_set_use_markup (GTK_LABEL (label), TRUE);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_CENTER);
    gtk_misc_set_alignment (GTK_MISC (label), 1, 0.5);

    RowHbox = gtk_hbox_new (FALSE, 10); gtk_widget_show (RowHbox);
    gtk_table_attach (GTK_TABLE (table1), RowHbox, 1, 2, 0, 1, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (GTK_FILL), 0, 0);

    checkbutton = gtk_check_button_new(); gtk_widget_show( checkbutton );
    gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON( checkbutton ), ShowAtmosphere );
    gtk_box_pack_start (GTK_BOX (RowHbox), checkbutton, FALSE, FALSE, 0);
    g_signal_connect( G_OBJECT( checkbutton ), "toggled", G_CALLBACK( ToggleAtmosphere ), NULL );


    /*
     * Row for Rayleigh Scattering
     * The rows on this page are two cols: col 1 is label col two holds hbox (which holds othwer hboxes)
     */
    label = gtk_label_new (_("<b>Rayleigh\nScattering:</b>")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, 0, 1, 1, 2, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);
    gtk_label_set_use_markup (GTK_LABEL (label), TRUE);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_CENTER);
    gtk_misc_set_alignment (GTK_MISC (label), 1, 0.5);

    RowHbox = gtk_hbox_new (FALSE, 10); gtk_widget_show (RowHbox);
    gtk_table_attach (GTK_TABLE (table1), RowHbox, 1, 2, 1, 2, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (GTK_FILL), 0, 0);


    // scattering constant label + button
    hbox = gtk_hbox_new (FALSE, 5); gtk_widget_show (hbox);
    gtk_box_pack_start (GTK_BOX (RowHbox), hbox, FALSE, TRUE, 0);

    label = gtk_label_new (_("Scatt.\nConst.")); gtk_widget_show (label);
    gtk_box_pack_start (GTK_BOX (hbox), label, FALSE, FALSE, 0);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_RIGHT);
    gtk_misc_set_alignment (GTK_MISC (label), 1, 0.5);

    spinbutton1_adj = gtk_adjustment_new (aInfo->Kr, 0.0001, 0.1000, 0.0001, 0.0010, 0);
    GtkWidget *RayleighScattConstButton = gtk_spin_button_new (GTK_ADJUSTMENT (spinbutton1_adj), 1, 5);
    gtk_widget_show (RayleighScattConstButton);
    gtk_box_pack_start (GTK_BOX (hbox), RayleighScattConstButton, FALSE, TRUE, 0);
    gtk_spin_button_set_numeric (GTK_SPIN_BUTTON (RayleighScattConstButton), TRUE);
    g_signal_connect( G_OBJECT( RayleighScattConstButton ), "value_changed", G_CALLBACK( ChangeRayleighScattConst ), NULL );



    // scale height label + button
    hbox = gtk_hbox_new (FALSE, 5); gtk_widget_show (hbox);
    gtk_box_pack_start (GTK_BOX (RowHbox), hbox, TRUE, TRUE, 0);

    label = gtk_label_new (_("Scale\nHeight")); gtk_widget_show (label);
    gtk_box_pack_start (GTK_BOX (hbox), label, FALSE, FALSE, 0);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_RIGHT);
    gtk_misc_set_alignment (GTK_MISC (label), 1, 0.5);

    spinbutton1_adj = gtk_adjustment_new (aInfo->rScaleHeight, 0.01, 1.0, 0.01, .05, 0);
    GtkWidget *RayleighScaleHeightButton = gtk_spin_button_new (GTK_ADJUSTMENT (spinbutton1_adj), 1, 3);
    gtk_widget_show (RayleighScaleHeightButton);
    gtk_box_pack_start (GTK_BOX (hbox), RayleighScaleHeightButton, FALSE, TRUE, 0);
    gtk_spin_button_set_numeric (GTK_SPIN_BUTTON (RayleighScaleHeightButton), TRUE);
    g_signal_connect( G_OBJECT( RayleighScaleHeightButton ), "value_changed", G_CALLBACK( ChangeRayleighScaleHeight ), NULL );


    /*
     * Row for Mie Scattering
     * The rows on this page are two cols: col 1 is label col two holds hbox (which holds othwer hboxes)
     */
    label = gtk_label_new (_("<b>Mie\nScattering:</b>")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, 0, 1, 2, 3, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);
    gtk_label_set_use_markup (GTK_LABEL (label), TRUE);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_CENTER);
    gtk_misc_set_alignment (GTK_MISC (label), 1, 0.5);

    RowHbox = gtk_hbox_new (FALSE, 10); gtk_widget_show (RowHbox);
    gtk_table_attach (GTK_TABLE (table1), RowHbox, 1, 2, 2, 3, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (GTK_FILL), 0, 0);


    // scattering constant label + button
    hbox = gtk_hbox_new (FALSE, 5); gtk_widget_show (hbox);
    gtk_box_pack_start (GTK_BOX (RowHbox), hbox, FALSE, TRUE, 0);

    label = gtk_label_new (_("Scatt.\nConst.")); gtk_widget_show (label);
    gtk_box_pack_start (GTK_BOX (hbox), label, FALSE, FALSE, 0);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_RIGHT);
    gtk_misc_set_alignment (GTK_MISC (label), 1, 0.5);

    spinbutton1_adj = gtk_adjustment_new (aInfo->Km, 0.0001, 0.1000, 0.0001, 0.0010, 0);
    GtkWidget *MieScattConstButton = gtk_spin_button_new (GTK_ADJUSTMENT (spinbutton1_adj), 1, 5);
    gtk_widget_show (MieScattConstButton);
    gtk_box_pack_start (GTK_BOX (hbox), MieScattConstButton, FALSE, TRUE, 0);
    gtk_spin_button_set_numeric (GTK_SPIN_BUTTON (MieScattConstButton), TRUE);
    g_signal_connect( G_OBJECT( MieScattConstButton ), "value_changed", G_CALLBACK( ChangeMieScattConst ), NULL );


    // scale height label + button
    hbox = gtk_hbox_new (FALSE, 5); gtk_widget_show (hbox);
    gtk_box_pack_start (GTK_BOX (RowHbox), hbox, TRUE, TRUE, 0);

    label = gtk_label_new (_("Scale\nHeight")); gtk_widget_show (label);
    gtk_box_pack_start (GTK_BOX (hbox), label, FALSE, FALSE, 0);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_RIGHT);
    gtk_misc_set_alignment (GTK_MISC (label), 1, 0.5);

    spinbutton1_adj = gtk_adjustment_new (aInfo->mScaleHeight, 0.01, 1.0, 0.01, .05, 0);
    GtkWidget *MieScaleHeightButton = gtk_spin_button_new (GTK_ADJUSTMENT (spinbutton1_adj), 1, 3);
    gtk_widget_show (MieScaleHeightButton);
    gtk_box_pack_start (GTK_BOX (hbox), MieScaleHeightButton, FALSE, TRUE, 0);
    gtk_spin_button_set_numeric (GTK_SPIN_BUTTON (MieScaleHeightButton), TRUE);
    g_signal_connect( G_OBJECT( MieScaleHeightButton ), "value_changed", G_CALLBACK( ChangeMieScaleHeight ), NULL );


    // assym factor label + button
    hbox = gtk_hbox_new (FALSE, 5); gtk_widget_show (hbox);
    gtk_box_pack_start (GTK_BOX (RowHbox), hbox, TRUE, TRUE, 0);

    label = gtk_label_new (_("Asym.\nFactor")); gtk_widget_show (label);
    gtk_box_pack_start (GTK_BOX (hbox), label, FALSE, FALSE, 0);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_RIGHT);
    gtk_misc_set_alignment (GTK_MISC (label), 1, 0.5);

    spinbutton1_adj = gtk_adjustment_new (aInfo->um, 0.01, 1.00, 0.01, 0.05, 0);
    GtkWidget *MieAsymFactorButton = gtk_spin_button_new (GTK_ADJUSTMENT (spinbutton1_adj), 1, 3);
    gtk_widget_show (MieAsymFactorButton);
    gtk_box_pack_start (GTK_BOX (hbox), MieAsymFactorButton, FALSE, TRUE, 0);
    gtk_spin_button_set_numeric (GTK_SPIN_BUTTON (MieAsymFactorButton), TRUE);
    g_signal_connect( G_OBJECT( MieAsymFactorButton ), "value_changed", G_CALLBACK( ChangeMieAsymFactor ), NULL );




    /*
     * Row for Sun
     * The rows on this page are two cols: col 1 is label col two holds hbox (which holds othwer hboxes)
     */
    label = gtk_label_new (_("<b>Solar Intensity:</b>")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, 0, 1, 3, 4, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);
    gtk_label_set_use_markup (GTK_LABEL (label), TRUE);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_CENTER);
    gtk_misc_set_alignment (GTK_MISC (label), 1, 0.5);

    RowHbox = gtk_hbox_new (FALSE, 10); gtk_widget_show (RowHbox);
    gtk_table_attach (GTK_TABLE (table1), RowHbox, 1, 2, 3, 4, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (GTK_FILL), 0, 0);

    // Intensity label+button
    hbox = gtk_hbox_new (FALSE, 5); gtk_widget_show (hbox);
    gtk_box_pack_start (GTK_BOX (RowHbox), hbox, FALSE, FALSE, 10);

//    label = gtk_label_new ("Intensity"); gtk_widget_show (label);
//    gtk_box_pack_start (GTK_BOX (hbox), label, FALSE, FALSE, 0);
//    gtk_label_set_use_markup (GTK_LABEL (label), TRUE);
//    gtk_misc_set_alignment (GTK_MISC (label), 1, 0.5);

    spinbutton1_adj = gtk_adjustment_new (15.0, 0, 100, 1, 10, 0);
    GtkWidget *SunIntensityButton = gtk_spin_button_new (GTK_ADJUSTMENT (spinbutton1_adj), 1, 2);
    gtk_widget_show (SunIntensityButton);
    gtk_box_pack_start (GTK_BOX (hbox), SunIntensityButton, FALSE, TRUE, 0);
    gtk_spin_button_set_numeric (GTK_SPIN_BUTTON (SunIntensityButton), TRUE);
    g_signal_connect( G_OBJECT( SunIntensityButton ), "value_changed", G_CALLBACK( ChangeSunIntensity ), NULL );


    /*
     * Row for Sun
     * The rows on this page are two cols: col 1 is label col two holds hbox (which holds othwer hboxes)
     */
    label = gtk_label_new (_("<b>Solar Color:</b>\n(\u019B's in nm)")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, 0, 1, 4, 5, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);
    gtk_label_set_use_markup (GTK_LABEL (label), TRUE);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_CENTER);
    gtk_misc_set_alignment (GTK_MISC (label), 1, 0.5);

    RowHbox = gtk_hbox_new (FALSE, 10); gtk_widget_show (RowHbox);
    gtk_table_attach (GTK_TABLE (table1), RowHbox, 1, 2, 4, 5, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (GTK_FILL), 0, 0);


    // Lambda_red label + button
    hbox = gtk_hbox_new (FALSE, 5); gtk_widget_show (hbox);
    gtk_box_pack_start (GTK_BOX (RowHbox), hbox, FALSE, TRUE, 0);

    label = gtk_label_new ("\u019B<sub>red</sub>"); gtk_widget_show (label);
    gtk_box_pack_start (GTK_BOX (hbox), label, FALSE, FALSE, 0);
    gtk_label_set_use_markup (GTK_LABEL (label), TRUE);
    gtk_misc_set_alignment (GTK_MISC (label), 1, 0.5);

    spinbutton1_adj = gtk_adjustment_new (650, 400, 800, 1, 10, 0);
    GtkWidget *SunLambdaRedButton = gtk_spin_button_new (GTK_ADJUSTMENT (spinbutton1_adj), 1, 0);
    gtk_widget_show (SunLambdaRedButton);
    gtk_box_pack_start (GTK_BOX (hbox), SunLambdaRedButton, FALSE, TRUE, 0);
    gtk_spin_button_set_numeric (GTK_SPIN_BUTTON (SunLambdaRedButton), TRUE);
    g_signal_connect( G_OBJECT( SunLambdaRedButton ), "value_changed", G_CALLBACK( ChangeSunLambdaRed ), NULL );

    // Lambda_grn label + button
    hbox = gtk_hbox_new (FALSE, 5); gtk_widget_show (hbox);
    gtk_box_pack_start (GTK_BOX (RowHbox), hbox, FALSE, TRUE, 0);

    label = gtk_label_new ("\u019B<sub>grn</sub>"); gtk_widget_show (label);
    gtk_box_pack_start (GTK_BOX (hbox), label, FALSE, FALSE, 0);
    gtk_label_set_use_markup (GTK_LABEL (label), TRUE);
    gtk_misc_set_alignment (GTK_MISC (label), 1, 0.5);

    spinbutton1_adj = gtk_adjustment_new (570, 400, 800, 1, 10, 0);
    GtkWidget *SunLambdaGrnButton = gtk_spin_button_new (GTK_ADJUSTMENT (spinbutton1_adj), 1, 0);
    gtk_widget_show (SunLambdaGrnButton);
    gtk_box_pack_start (GTK_BOX (hbox), SunLambdaGrnButton, FALSE, TRUE, 0);
    gtk_spin_button_set_numeric (GTK_SPIN_BUTTON (SunLambdaGrnButton), TRUE);
    g_signal_connect( G_OBJECT( SunLambdaGrnButton ), "value_changed", G_CALLBACK( ChangeSunLambdaGrn ), NULL );

    // Lambda_red label + button
    hbox = gtk_hbox_new (FALSE, 5); gtk_widget_show (hbox);
    gtk_box_pack_start (GTK_BOX (RowHbox), hbox, FALSE, TRUE, 0);

    label = gtk_label_new ("\u019B<sub>blu</sub>"); gtk_widget_show (label);
    gtk_box_pack_start (GTK_BOX (hbox), label, FALSE, FALSE, 0);
    gtk_label_set_use_markup (GTK_LABEL (label), TRUE);
    gtk_misc_set_alignment (GTK_MISC (label), 1, 0.5);

    spinbutton1_adj = gtk_adjustment_new (475, 400, 800, 1, 10, 0);
    GtkWidget *SunLambdaBluButton = gtk_spin_button_new (GTK_ADJUSTMENT (spinbutton1_adj), 1, 0);
    gtk_widget_show (SunLambdaBluButton);
    gtk_box_pack_start (GTK_BOX (hbox), SunLambdaBluButton, FALSE, TRUE, 0);
    gtk_spin_button_set_numeric (GTK_SPIN_BUTTON (SunLambdaBluButton), TRUE);
    g_signal_connect( G_OBJECT( SunLambdaBluButton ), "value_changed", G_CALLBACK( ChangeSunLambdaBlu ), NULL );



    /*
     * Row for Outer Sphere
     * The rows on this page are two cols: col 1 is label col two holds hbox (which holds othwer hboxes)
     */
    label = gtk_label_new (_("<b>Outer Sphere:</b>")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, 0, 1, 5, 6, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);
    gtk_label_set_use_markup (GTK_LABEL (label), TRUE);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_CENTER);
    gtk_misc_set_alignment (GTK_MISC (label), 1, 0.5);

    RowHbox = gtk_hbox_new (FALSE, 10); gtk_widget_show (RowHbox);
    gtk_table_attach (GTK_TABLE (table1), RowHbox, 1, 2, 5, 6, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (GTK_FILL), 0, 0);


    // Radius label+button
    hbox = gtk_hbox_new (FALSE, 5); gtk_widget_show (hbox);
    gtk_box_pack_start (GTK_BOX (RowHbox), hbox, FALSE, TRUE, 0);

    label = gtk_label_new ("Radius, R<sub>e</sub>"); gtk_widget_show (label);
    gtk_box_pack_start (GTK_BOX (hbox), label, FALSE, FALSE, 0);
    gtk_label_set_use_markup (GTK_LABEL (label), TRUE);
    gtk_misc_set_alignment (GTK_MISC (label), 1, 0.5);

spinbutton1_adj = gtk_adjustment_new(1.03, 1.0, 10.0, 1, 10, 0);
GtkWidget *OuterSphereRButton = gtk_spin_button_new (GTK_ADJUSTMENT (spinbutton1_adj), 1, 3);
gtk_widget_show (OuterSphereRButton);
gtk_box_pack_start (GTK_BOX (hbox), OuterSphereRButton, FALSE, TRUE, 0);
gtk_spin_button_set_numeric (GTK_SPIN_BUTTON (OuterSphereRButton), TRUE);
// HOOK ME UP

    // n label+button
    hbox = gtk_hbox_new (FALSE, 5); gtk_widget_show (hbox);
    gtk_box_pack_start (GTK_BOX (RowHbox), hbox, FALSE, TRUE, 0);

    label = gtk_label_new ("n"); gtk_widget_show (label);
    gtk_box_pack_start (GTK_BOX (hbox), label, FALSE, FALSE, 0);
    gtk_label_set_use_markup (GTK_LABEL (label), TRUE);
    gtk_misc_set_alignment (GTK_MISC (label), 1, 0.5);

spinbutton1_adj = gtk_adjustment_new (200, 10, 500, 1, 10, 0);
GtkWidget *OuterSpherenButton = gtk_spin_button_new (GTK_ADJUSTMENT (spinbutton1_adj), 1, 0);
gtk_widget_show (OuterSpherenButton);
gtk_box_pack_start (GTK_BOX (hbox), OuterSpherenButton, FALSE, TRUE, 0);
gtk_spin_button_set_numeric (GTK_SPIN_BUTTON (OuterSpherenButton), TRUE);
// HOOK ME UP


    /*
     * Row for Inner Sphere
     * The rows on this page are two cols: col 1 is label col two holds hbox (which holds othwer hboxes)
     */
    label = gtk_label_new (_("<b>Inner Sphere:</b>")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, 0, 1, 6, 7, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);
    gtk_label_set_use_markup (GTK_LABEL (label), TRUE);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_CENTER);
    gtk_misc_set_alignment (GTK_MISC (label), 1, 0.5);

    RowHbox = gtk_hbox_new (FALSE, 10); gtk_widget_show (RowHbox);
    gtk_table_attach (GTK_TABLE (table1), RowHbox, 1, 2, 6, 7, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (GTK_FILL), 0, 0);


    // Radius label+button
    hbox = gtk_hbox_new (FALSE, 5); gtk_widget_show (hbox);
    gtk_box_pack_start (GTK_BOX (RowHbox), hbox, FALSE, TRUE, 0);

    label = gtk_label_new ("Radius, R<sub>e</sub>"); gtk_widget_show (label);
    gtk_box_pack_start (GTK_BOX (hbox), label, FALSE, FALSE, 0);
    gtk_label_set_use_markup (GTK_LABEL (label), TRUE);
    gtk_misc_set_alignment (GTK_MISC (label), 1, 0.5);

spinbutton1_adj = gtk_adjustment_new (1.001, 1.0, 10.0, 0.01, 10, 0);
GtkWidget *InnerSphereRButton = gtk_spin_button_new (GTK_ADJUSTMENT (spinbutton1_adj), 1, 3);
gtk_widget_show (InnerSphereRButton);
gtk_box_pack_start (GTK_BOX (hbox), InnerSphereRButton, FALSE, TRUE, 0);
gtk_spin_button_set_numeric (GTK_SPIN_BUTTON (InnerSphereRButton), TRUE);
// HOOK ME UP

    // n label+button
    hbox = gtk_hbox_new (FALSE, 5); gtk_widget_show (hbox);
    gtk_box_pack_start (GTK_BOX (RowHbox), hbox, FALSE, TRUE, 0);

    label = gtk_label_new ("n"); gtk_widget_show (label);
    gtk_box_pack_start (GTK_BOX (hbox), label, FALSE, FALSE, 0);
    gtk_label_set_use_markup (GTK_LABEL (label), TRUE);
    gtk_misc_set_alignment (GTK_MISC (label), 1, 0.5);

spinbutton1_adj = gtk_adjustment_new (100, 10, 500, 1, 10, 0);
GtkWidget *InnerSpherenButton = gtk_spin_button_new (GTK_ADJUSTMENT (spinbutton1_adj), 1, 0);
gtk_widget_show (InnerSpherenButton);
gtk_box_pack_start (GTK_BOX (hbox), InnerSpherenButton, FALSE, TRUE, 0);
gtk_spin_button_set_numeric (GTK_SPIN_BUTTON (InnerSpherenButton), TRUE);
// HOOK ME UP




    /*
     * Row for Samples
     * The rows on this page are two cols: col 1 is label col two holds hbox (which holds othwer hboxes)
     */
    label = gtk_label_new (_("<b># Samples:</b>")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, 0, 1, 7, 8, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);
    gtk_label_set_use_markup (GTK_LABEL (label), TRUE);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_CENTER);
    gtk_misc_set_alignment (GTK_MISC (label), 1, 0.5);

    RowHbox = gtk_hbox_new (FALSE, 10); gtk_widget_show (RowHbox);
    gtk_table_attach (GTK_TABLE (table1), RowHbox, 1, 2, 7, 8, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (GTK_FILL), 0, 0);


    // Sample Rays label+button
    hbox = gtk_hbox_new (FALSE, 5); gtk_widget_show (hbox);
    gtk_box_pack_start (GTK_BOX (RowHbox), hbox, FALSE, TRUE, 0);

    label = gtk_label_new ("Sample Rays"); gtk_widget_show (label);
    gtk_box_pack_start (GTK_BOX (hbox), label, FALSE, FALSE, 0);
    gtk_label_set_use_markup (GTK_LABEL (label), TRUE);
    gtk_misc_set_alignment (GTK_MISC (label), 1, 0.5);

    spinbutton1_adj = gtk_adjustment_new (6, 1, 100, 1, 10, 0);
    GtkWidget *SampleRaysButton = gtk_spin_button_new (GTK_ADJUSTMENT (spinbutton1_adj), 1, 0);
    gtk_widget_show (SampleRaysButton);
    gtk_box_pack_start (GTK_BOX (hbox), SampleRaysButton, FALSE, TRUE, 0);
    gtk_spin_button_set_numeric (GTK_SPIN_BUTTON (SampleRaysButton), TRUE);
    g_signal_connect( G_OBJECT( SampleRaysButton ), "value_changed", G_CALLBACK( ChangeSampleRays ), NULL );

    // Depth Buffer Samples  label+button
    hbox = gtk_hbox_new (FALSE, 5); gtk_widget_show (hbox);
    gtk_box_pack_start (GTK_BOX (RowHbox), hbox, FALSE, TRUE, 0);

    label = gtk_label_new ("Depth Buffer Samples"); gtk_widget_show (label);
    gtk_box_pack_start (GTK_BOX (hbox), label, FALSE, FALSE, 0);
    gtk_label_set_use_markup (GTK_LABEL (label), TRUE);
    gtk_misc_set_alignment (GTK_MISC (label), 1, 0.5);

    spinbutton1_adj = gtk_adjustment_new (20, 2, 100, 1, 10, 0);
    GtkWidget *DepthBufSamplesButton = gtk_spin_button_new (GTK_ADJUSTMENT (spinbutton1_adj), 1, 0);
    gtk_widget_show (DepthBufSamplesButton);
    gtk_box_pack_start (GTK_BOX (hbox), DepthBufSamplesButton, FALSE, TRUE, 0);
    gtk_spin_button_set_numeric (GTK_SPIN_BUTTON (DepthBufSamplesButton), TRUE);
    g_signal_connect( G_OBJECT( DepthBufSamplesButton ), "value_changed", G_CALLBACK( ChangeDepthBufSamples ), NULL );





    /**************
     * Stars Page *
     **************/
    vbox2 = gtk_vbox_new (FALSE, 0); gtk_widget_show (vbox2);
    gtk_container_set_border_width (GTK_CONTAINER (vbox2), 20);

    // Page title
    label = gtk_label_new (_("<b><span size=\"small\">Misc</span></b>")); gtk_widget_show (label);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_CENTER);
    gtk_label_set_use_markup( GTK_LABEL(label), TRUE );
    gtk_notebook_append_page( GTK_NOTEBOOK(notebook), vbox2, label );

    table1 = gtk_table_new (7, 2, FALSE); gtk_widget_show (table1);
    gtk_box_pack_start (GTK_BOX (vbox2), table1, TRUE, TRUE, 15);
    gtk_container_set_border_width (GTK_CONTAINER (table1), 20);
    gtk_table_set_row_spacings (GTK_TABLE (table1), 20);
    gtk_table_set_col_spacings (GTK_TABLE (table1), 40);

    /*
     * Row for Toggle button
     * The rows on this page are two cols: col 1 is label col two holds hbox (which holds othwer hboxes)
     */
    RowHbox = gtk_hbox_new (FALSE, 10); gtk_widget_show (RowHbox);
    gtk_table_attach (GTK_TABLE (table1), RowHbox, 0, 1, 0, 1, (GtkAttachOptions) (0), (GtkAttachOptions) (GTK_FILL), 0, 0);

    checkbutton = gtk_check_button_new(); gtk_widget_show( checkbutton );
    gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON( checkbutton ), ShowStars );
    gtk_box_pack_start (GTK_BOX (RowHbox), checkbutton, FALSE, FALSE, 0);
    g_signal_connect( G_OBJECT( checkbutton ), "toggled", G_CALLBACK( ToggleStars ), NULL );

    label = gtk_label_new (_("<b>Show Stars</b>")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, 1, 2, 0, 1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
    gtk_label_set_use_markup (GTK_LABEL (label), TRUE);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_CENTER);
    gtk_misc_set_alignment (GTK_MISC (label), 1, 0.5);


    /*
     * Row for Max Magnitude
     */
    label = gtk_label_new (_("<b>Star Magnitude:</b>")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, 2, 3, 0, 1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
    gtk_label_set_use_markup (GTK_LABEL (label), TRUE);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_CENTER);
    gtk_misc_set_alignment (GTK_MISC (label), 1, 0.5);

    RowHbox = gtk_hbox_new (FALSE, 10); gtk_widget_show (RowHbox);
    gtk_table_attach (GTK_TABLE (table1), RowHbox, 3, 4, 0, 1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);

    // Max Magnitude button
    hbox = gtk_hbox_new (FALSE, 5); gtk_widget_show (hbox);
    gtk_box_pack_start (GTK_BOX (RowHbox), hbox, FALSE, TRUE, 0);

    label = gtk_label_new ("Maximum"); gtk_widget_show (label);
    gtk_box_pack_start (GTK_BOX (hbox), label, FALSE, FALSE, 0);
    gtk_label_set_use_markup (GTK_LABEL (label), TRUE);
    gtk_misc_set_alignment (GTK_MISC (label), 1, 0.5);

    spinbutton1_adj = gtk_adjustment_new(8.0, -4.0, 10.0, .1, .5, 0);
    GtkWidget *StarsMaxMagButton = gtk_spin_button_new (GTK_ADJUSTMENT (spinbutton1_adj), 1, 3);
    gtk_widget_show (StarsMaxMagButton);
    gtk_box_pack_start (GTK_BOX (hbox), StarsMaxMagButton, FALSE, TRUE, 0);
    gtk_spin_button_set_numeric (GTK_SPIN_BUTTON (StarsMaxMagButton), TRUE);
    g_signal_connect( G_OBJECT( StarsMaxMagButton ), "value_changed", G_CALLBACK( ChangeStarsMaxMag ), NULL );




    /*************
     * View Page *
     *************/
    vbox2 = gtk_vbox_new (FALSE, 0); gtk_widget_show (vbox2);
    gtk_container_set_border_width (GTK_CONTAINER (vbox2), 20);

    // Page title
    label = gtk_label_new (_("<b><span size=\"small\">View</span></b>")); gtk_widget_show (label);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_CENTER);
    gtk_label_set_use_markup( GTK_LABEL(label), TRUE );
    gtk_notebook_append_page( GTK_NOTEBOOK(notebook), vbox2, label );



    table1 = gtk_table_new (3, 3, FALSE);
    gtk_widget_show (table1);
    gtk_container_add (GTK_CONTAINER (vbox2), table1);

    filechooserbutton = gtk_file_chooser_button_new (_("Select A File"), GTK_FILE_CHOOSER_ACTION_OPEN);
    gtk_file_chooser_set_current_folder_file( GTK_FILE_CHOOSER(filechooserbutton), g_file_new_for_path( "/data1/mgh/BlueMarble/5400x2700" ), NULL );
    gtk_file_chooser_set_filename( GTK_FILE_CHOOSER(filechooserbutton), "/data1/mgh/BlueMarble/5400x2700/world.topo.bathy.200406.3x5400x2700.png");
    gtk_file_chooser_select_filename( GTK_FILE_CHOOSER(filechooserbutton), "/data1/mgh/BlueMarble/5400x2700/world.topo.bathy.200406.3x5400x2700.png");
    gtk_file_chooser_set_file( GTK_FILE_CHOOSER(filechooserbutton), g_file_new_for_path("/data1/mgh/BlueMarble/5400x2700/world.topo.bathy.200406.3x5400x2700.png"), NULL );

    PngFilter = gtk_file_filter_new( );
    gtk_file_filter_set_name( PngFilter, "PNG Images" );
    gtk_file_filter_add_pattern( PngFilter, "*.png" );
    gtk_file_chooser_add_filter( GTK_FILE_CHOOSER(filechooserbutton), PngFilter );
    gtk_file_chooser_set_filter( GTK_FILE_CHOOSER(filechooserbutton), PngFilter );

    gtk_widget_show( filechooserbutton );
    gtk_table_attach( GTK_TABLE (table1), filechooserbutton, 1, 2, 0, 1, (GtkAttachOptions) (GTK_EXPAND | GTK_FILL), (GtkAttachOptions) (GTK_FILL), 0, 5);

    g_signal_connect( G_OBJECT( filechooserbutton ), "file-set" , G_CALLBACK( ChangeMapImage ), NULL );



    label = gtk_label_new (_("Map Image: "));
    gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, 0, 1, 0, 1, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 10);
    gtk_misc_set_alignment (GTK_MISC (label), 0, 0.5);



    /*
     * Row for quat_view Quaternion
     */
    label = gtk_label_new (_("<b>View Quat:</b>")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, 0, 1, 2, 3, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);
    gtk_label_set_use_markup (GTK_LABEL (label), TRUE);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_CENTER);
    gtk_misc_set_alignment (GTK_MISC (label), 1, 0.5);

    RowHbox = gtk_hbox_new (FALSE, 10); gtk_widget_show (RowHbox);
    gtk_table_attach (GTK_TABLE (table1), RowHbox, 1, 2, 2, 3, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (GTK_FILL), 0, 0);

    hbox = gtk_hbox_new (FALSE, 5); gtk_widget_show (hbox);
    gtk_box_pack_start (GTK_BOX (RowHbox), hbox, FALSE, TRUE, 10);

    spinbutton1_adj = gtk_adjustment_new ( 0.0, -1.0, 1.0, 0.01, 0.10, 0);
    GtkWidget *ViewQuatButton0 = gtk_spin_button_new (GTK_ADJUSTMENT (spinbutton1_adj), 0.001, 5);
    gtk_widget_show( ViewQuatButton0 );
    gtk_box_pack_start( GTK_BOX (hbox), ViewQuatButton0, FALSE, TRUE, 0);
    gtk_spin_button_set_numeric( GTK_SPIN_BUTTON(ViewQuatButton0), TRUE);
    g_signal_connect( G_OBJECT( ViewQuatButton0 ), "value_changed", G_CALLBACK( ChangeViewQuat ), GINT_TO_POINTER( 0 ) );

    spinbutton1_adj = gtk_adjustment_new ( -sqrt(0.5), -1.0, 1.0, 0.01, 0.10, 0);
    GtkWidget *ViewQuatButton1 = gtk_spin_button_new (GTK_ADJUSTMENT (spinbutton1_adj), 0.001, 5);
    gtk_widget_show( ViewQuatButton1 );
    gtk_box_pack_start( GTK_BOX (hbox), ViewQuatButton1, FALSE, TRUE, 0);
    gtk_spin_button_set_numeric( GTK_SPIN_BUTTON(ViewQuatButton1), TRUE);
    g_signal_connect( G_OBJECT( ViewQuatButton1 ), "value_changed", G_CALLBACK( ChangeViewQuat ), GINT_TO_POINTER( 1 ) );

    spinbutton1_adj = gtk_adjustment_new ( -sqrt(0.5), -1.0, 1.0, 0.01, 0.10, 0);
    GtkWidget *ViewQuatButton2 = gtk_spin_button_new (GTK_ADJUSTMENT (spinbutton1_adj), 0.001, 5);
    gtk_widget_show( ViewQuatButton2 );
    gtk_box_pack_start( GTK_BOX (hbox), ViewQuatButton2, FALSE, TRUE, 0);
    gtk_spin_button_set_numeric( GTK_SPIN_BUTTON(ViewQuatButton2), TRUE);
    g_signal_connect( G_OBJECT( ViewQuatButton2 ), "value_changed", G_CALLBACK( ChangeViewQuat ), GINT_TO_POINTER( 2 ) );

    spinbutton1_adj = gtk_adjustment_new ( 0.0, -1.0, 1.0, 0.01, 0.10, 0);
    GtkWidget *ViewQuatButton3 = gtk_spin_button_new (GTK_ADJUSTMENT (spinbutton1_adj), 0.001, 5);
    gtk_widget_show( ViewQuatButton3 );
    gtk_box_pack_start( GTK_BOX (hbox), ViewQuatButton3, FALSE, TRUE, 0);
    gtk_spin_button_set_numeric( GTK_SPIN_BUTTON(ViewQuatButton3), TRUE);
    g_signal_connect( G_OBJECT( ViewQuatButton3 ), "value_changed", G_CALLBACK( ChangeViewQuat ), GINT_TO_POINTER( 3 ) );




    /**********************
     * TLE Generater Page *
     **********************/
    vbox2 = gtk_vbox_new( FALSE, 40 ); gtk_widget_show( vbox2 );
    gtk_container_set_border_width(GTK_CONTAINER( vbox2 ), 20);

    // Page title
    label = gtk_label_new( _("<b><span size=\"small\">TLE</span></b>")); gtk_widget_show( label);
    gtk_label_set_justify( GTK_LABEL( label), GTK_JUSTIFY_CENTER);
    gtk_label_set_use_markup( GTK_LABEL(label), TRUE );
    gtk_notebook_append_page( GTK_NOTEBOOK(notebook), vbox2, label );


    
    // Vbox with no spacing to hold the 3 line element strings
    vbox3 = gtk_vbox_new( FALSE, 0 ); gtk_widget_show( vbox3 );
    gtk_container_set_border_width( GTK_CONTAINER( vbox3 ), 60 );
    gtk_box_pack_start( GTK_BOX(vbox2), vbox3, FALSE, TRUE, 10 );

    sprintf( Str, "<b><big><tt><span>%s</span></tt></big></b>", tle.Line0 );
    label = gtk_label_new( Str ); gtk_widget_show( label);
    gtk_box_pack_start( GTK_BOX(vbox3), label, FALSE, TRUE, 0);

    gtk_label_set_use_markup( GTK_LABEL( label), TRUE);
    gtk_label_set_justify( GTK_LABEL( label), GTK_JUSTIFY_CENTER);
    gtk_misc_set_alignment( GTK_MISC( label), 0.0, 0.5);
    tle.Line0Label = label;

    sprintf( Str, "<b><big><tt><span>%s</span></tt></big></b>", tle.Line1 );
    label = gtk_label_new( Str ); gtk_widget_show( label);
    gtk_box_pack_start( GTK_BOX(vbox3), label, FALSE, TRUE, 0);
    gtk_label_set_use_markup( GTK_LABEL( label), TRUE);
    gtk_label_set_justify( GTK_LABEL( label), GTK_JUSTIFY_CENTER);
    gtk_misc_set_alignment( GTK_MISC( label), 0.0, 0.5);
    tle.Line1Label = label;

    sprintf( Str, "<b><big><tt><span>%s</span></tt></big></b>\n\n", tle.Line2 );
    label = gtk_label_new( Str ); gtk_widget_show( label);
    gtk_box_pack_start( GTK_BOX(vbox3), label, FALSE, TRUE, 0);
    gtk_label_set_use_markup( GTK_LABEL( label), TRUE);
    gtk_label_set_justify( GTK_LABEL( label), GTK_JUSTIFY_CENTER);
    gtk_misc_set_alignment( GTK_MISC( label), 0.0, 0.5);
    tle.Line2Label = label;



    ChangeTLE( NULL, GINT_TO_POINTER(-1) );

    int r = 2;

    table1 = gtk_table_new( 3, 3, FALSE );
    gtk_widget_show( table1);
    gtk_container_add( GTK_CONTAINER( vbox2), table1);

    /*
     * Section for Line 1 Entries
     *  01  Line Number of Element Data
     *  03-07   Satellite Number
     *  08  Classification (U=Unclassified)
     *
     *  10-11   International Designator (Last two digits of launch year)
     *  12-14   International Designator (Launch number of the year)
     *  15-17   International Designator (Piece of the launch)
     *
     *  19-20   Epoch Year (Last two digits of year)
     *  21-32   Epoch (Day of the year and fractional portion of the day)
     *
     *  34-43   First Time Derivative of the Mean Motion
     *  45-52   Second Time Derivative of Mean Motion (decimal point assumed)
     *  54-61   BSTAR drag term (decimal point assumed)
     *  63  Ephemeris type
     *  65-68   Element number
     *  69  Checksum (Modulo 10)
     *  (Letters, blanks, periods, plus signs = 0; minus signs = 1)
     */
    label = gtk_label_new( _("<b>Line 1 Entries:</b>")); gtk_widget_show( label);
    gtk_table_attach( GTK_TABLE( table1), label, 0, 1, r, r+1,( GtkAttachOptions)( GTK_FILL),( GtkAttachOptions)( 0), 0, 0);
    gtk_label_set_use_markup( GTK_LABEL( label), TRUE);
    gtk_label_set_justify( GTK_LABEL( label), GTK_JUSTIFY_CENTER);
    gtk_misc_set_alignment( GTK_MISC( label), 1, 0.5);
    ++r;

    // Satellite Number
    label = gtk_label_new( _("<b>Satellite Number:</b>")); gtk_widget_show( label);
    gtk_table_attach( GTK_TABLE( table1), label, 1, 2, r, r+1,( GtkAttachOptions)( GTK_FILL),( GtkAttachOptions)( 0), 0, 0);
    gtk_label_set_use_markup( GTK_LABEL( label), TRUE);
    gtk_label_set_justify( GTK_LABEL( label), GTK_JUSTIFY_CENTER);
    gtk_misc_set_alignment( GTK_MISC( label), 1, 0.5);

    RowHbox = gtk_hbox_new( FALSE, 10); gtk_widget_show( RowHbox);
    gtk_table_attach( GTK_TABLE( table1), RowHbox, 2, 3, r, r+1,( GtkAttachOptions)( GTK_FILL),( GtkAttachOptions)( GTK_FILL), 0, 0);

    hbox = gtk_hbox_new( FALSE, 2); gtk_widget_show( hbox);
    gtk_box_pack_start( GTK_BOX( RowHbox), hbox, FALSE, TRUE, 10);

    spinbutton1_adj = gtk_adjustment_new(  tle.SatNum, 0, 99999, 1, 1, 0);
    GtkWidget *TleGenButton0 = gtk_spin_button_new( GTK_ADJUSTMENT( spinbutton1_adj), 0.001, 0);
    gtk_widget_show( TleGenButton0 );
    gtk_box_pack_start( GTK_BOX( hbox), TleGenButton0, FALSE, TRUE, 0);
    gtk_spin_button_set_numeric( GTK_SPIN_BUTTON(TleGenButton0), TRUE);
    g_signal_connect( G_OBJECT( TleGenButton0 ), "value_changed", G_CALLBACK( ChangeTLE ), GINT_TO_POINTER( 1 ) );
    ++r;

    // International Designator
    label = gtk_label_new( _("<b>International Designator:</b>")); gtk_widget_show( label);
    gtk_table_attach( GTK_TABLE( table1), label, 1, 2, r, r+1,( GtkAttachOptions)( GTK_FILL),( GtkAttachOptions)( 0), 0, 0);
    gtk_label_set_use_markup( GTK_LABEL( label), TRUE);
    gtk_label_set_justify( GTK_LABEL( label), GTK_JUSTIFY_CENTER);
    gtk_misc_set_alignment( GTK_MISC( label), 1, 0.5);

    RowHbox = gtk_hbox_new( FALSE, 10); gtk_widget_show( RowHbox);
    gtk_table_attach( GTK_TABLE( table1), RowHbox, 2, 3, r, r+1,( GtkAttachOptions)( GTK_FILL),( GtkAttachOptions)( GTK_FILL), 0, 0);

    hbox = gtk_hbox_new( FALSE, 2); gtk_widget_show( hbox);
    gtk_box_pack_start( GTK_BOX( RowHbox), hbox, FALSE, TRUE, 10);

//    spinbutton1_adj = gtk_adjustment_new(  22.3, -90.0, 90.0, 0.01, 0.10, 0);
//    GtkWidget *TleGenButton1 = gtk_spin_button_new( GTK_ADJUSTMENT( spinbutton1_adj), 0.001, 5);
//    gtk_widget_show( TleGenButton1 );
//    gtk_box_pack_start( GTK_BOX( hbox), TleGenButton1, FALSE, TRUE, 0);
//    gtk_spin_button_set_numeric( GTK_SPIN_BUTTON(TleGenButton1), TRUE);
//    g_signal_connect( G_OBJECT( TleGenButton1 ), "value_changed", G_CALLBACK( ChangeTLE ), GINT_TO_POINTER( 2 ) );
    ++r;

    // Epoch
    label = gtk_label_new( _("<b>Epoch:</b>")); gtk_widget_show( label);
    gtk_table_attach( GTK_TABLE( table1), label, 1, 2, r, r+1,( GtkAttachOptions)( GTK_FILL),( GtkAttachOptions)( 0), 0, 0);
    gtk_label_set_use_markup( GTK_LABEL( label), TRUE);
    gtk_label_set_justify( GTK_LABEL( label), GTK_JUSTIFY_CENTER);
    gtk_misc_set_alignment( GTK_MISC( label), 1, 0.5);

    RowHbox = gtk_hbox_new( FALSE, 10); gtk_widget_show( RowHbox);
    gtk_table_attach( GTK_TABLE( table1), RowHbox, 2, 3, r, r+1,( GtkAttachOptions)( GTK_FILL),( GtkAttachOptions)( GTK_FILL), 0, 0);

    hbox = gtk_hbox_new( FALSE, 2); gtk_widget_show( hbox);
    gtk_box_pack_start( GTK_BOX( RowHbox), hbox, FALSE, TRUE, 10);

    spinbutton1_adj = gtk_adjustment_new(  tle.Epoch, 1000.0, 99999.0, 0.01, 0.10, 0);
    GtkWidget *TleGenButton2 = gtk_spin_button_new( GTK_ADJUSTMENT( spinbutton1_adj), 0.001, 8);
    gtk_widget_show( TleGenButton2 );
    gtk_box_pack_start( GTK_BOX( hbox), TleGenButton2, FALSE, TRUE, 0);
    gtk_spin_button_set_numeric( GTK_SPIN_BUTTON(TleGenButton2), TRUE);
    g_signal_connect( G_OBJECT( TleGenButton2 ), "value_changed", G_CALLBACK( ChangeTLE ), GINT_TO_POINTER( 3 ) );
    ++r;

    // First Time Derivative of the Mean Motion
    label = gtk_label_new( _("<b>d/dt of Mean Motion:</b>")); gtk_widget_show( label);
    gtk_table_attach( GTK_TABLE( table1), label, 1, 2, r, r+1,( GtkAttachOptions)( GTK_FILL),( GtkAttachOptions)( 0), 0, 0);
    gtk_label_set_use_markup( GTK_LABEL( label), TRUE);
    gtk_label_set_justify( GTK_LABEL( label), GTK_JUSTIFY_CENTER);
    gtk_misc_set_alignment( GTK_MISC( label), 1, 0.5);

    RowHbox = gtk_hbox_new( FALSE, 10); gtk_widget_show( RowHbox);
    gtk_table_attach( GTK_TABLE( table1), RowHbox, 2, 3, r, r+1,( GtkAttachOptions)( GTK_FILL),( GtkAttachOptions)( GTK_FILL), 0, 0);

    hbox = gtk_hbox_new( FALSE, 2); gtk_widget_show( hbox);
    gtk_box_pack_start( GTK_BOX( RowHbox), hbox, FALSE, TRUE, 10);

    spinbutton1_adj = gtk_adjustment_new(  tle.d1MeanMotion, -0.5, 0.5, 0.01, 0.10, 0);
    GtkWidget *TleGenButton3 = gtk_spin_button_new( GTK_ADJUSTMENT( spinbutton1_adj), 0.001, 8);
    gtk_widget_show( TleGenButton3 );
    gtk_box_pack_start( GTK_BOX( hbox), TleGenButton3, FALSE, TRUE, 0);
    gtk_spin_button_set_numeric( GTK_SPIN_BUTTON(TleGenButton3), TRUE);
    g_signal_connect( G_OBJECT( TleGenButton3 ), "value_changed", G_CALLBACK( ChangeTLE ), GINT_TO_POINTER( 4 ) );
    ++r;

    // Second Time Derivative of the Mean Motion
    label = gtk_label_new( _("<b>d^2/dt^2 of Mean Motion:</b>")); gtk_widget_show( label);
    gtk_table_attach( GTK_TABLE( table1), label, 1, 2, r, r+1,( GtkAttachOptions)( GTK_FILL),( GtkAttachOptions)( 0), 0, 0);
    gtk_label_set_use_markup( GTK_LABEL( label), TRUE);
    gtk_label_set_justify( GTK_LABEL( label), GTK_JUSTIFY_CENTER);
    gtk_misc_set_alignment( GTK_MISC( label), 1, 0.5);

    RowHbox = gtk_hbox_new( FALSE, 10); gtk_widget_show( RowHbox);
    gtk_table_attach( GTK_TABLE( table1), RowHbox, 2, 3, r, r+1,( GtkAttachOptions)( GTK_FILL),( GtkAttachOptions)( GTK_FILL), 0, 0);

    hbox = gtk_hbox_new( FALSE, 2); gtk_widget_show( hbox);
    gtk_box_pack_start( GTK_BOX( RowHbox), hbox, FALSE, TRUE, 10);

    spinbutton1_adj = gtk_adjustment_new(  tle.d2MeanMotion, -0.1, 0.1, 0.01, 0.10, 0);
    GtkWidget *TleGenButton4 = gtk_spin_button_new( GTK_ADJUSTMENT( spinbutton1_adj), 0.001, 8);
    gtk_widget_show( TleGenButton4 );
    gtk_box_pack_start( GTK_BOX( hbox), TleGenButton4, FALSE, TRUE, 0);
    gtk_spin_button_set_numeric( GTK_SPIN_BUTTON(TleGenButton4), TRUE);
    g_signal_connect( G_OBJECT( TleGenButton4 ), "value_changed", G_CALLBACK( ChangeTLE ), GINT_TO_POINTER( 5 ) );
    ++r;

    // BSTAR
    label = gtk_label_new( _("<b>BSTAR:</b>")); gtk_widget_show( label);
    gtk_table_attach( GTK_TABLE( table1), label, 1, 2, r, r+1,( GtkAttachOptions)( GTK_FILL),( GtkAttachOptions)( 0), 0, 0);
    gtk_label_set_use_markup( GTK_LABEL( label), TRUE);
    gtk_label_set_justify( GTK_LABEL( label), GTK_JUSTIFY_CENTER);
    gtk_misc_set_alignment( GTK_MISC( label), 1, 0.5);

    RowHbox = gtk_hbox_new( FALSE, 10); gtk_widget_show( RowHbox);
    gtk_table_attach( GTK_TABLE( table1), RowHbox, 2, 3, r, r+1,( GtkAttachOptions)( GTK_FILL),( GtkAttachOptions)( GTK_FILL), 0, 0);

    hbox = gtk_hbox_new( FALSE, 2); gtk_widget_show( hbox);
    gtk_box_pack_start( GTK_BOX( RowHbox), hbox, FALSE, TRUE, 10);

    spinbutton1_adj = gtk_adjustment_new(  tle.BSTAR, 0.0, 1.0, 0.001e-3, 0.010e-3, 0);
    GtkWidget *TleGenButton5 = gtk_spin_button_new( GTK_ADJUSTMENT( spinbutton1_adj), 0.0001e-3, 12);
    gtk_widget_show( TleGenButton5 );
    gtk_box_pack_start( GTK_BOX( hbox), TleGenButton5, FALSE, TRUE, 0);
    gtk_spin_button_set_numeric( GTK_SPIN_BUTTON(TleGenButton5), TRUE);
    g_signal_connect( G_OBJECT( TleGenButton5 ), "value_changed", G_CALLBACK( ChangeTLE ), GINT_TO_POINTER( 6 ) );
    ++r;


    // Element Number
    label = gtk_label_new( _("<b>Element Number:</b>")); gtk_widget_show( label);
    gtk_table_attach( GTK_TABLE( table1), label, 1, 2, r, r+1,( GtkAttachOptions)( GTK_FILL),( GtkAttachOptions)( 0), 0, 0);
    gtk_label_set_use_markup( GTK_LABEL( label), TRUE);
    gtk_label_set_justify( GTK_LABEL( label), GTK_JUSTIFY_CENTER);
    gtk_misc_set_alignment( GTK_MISC( label), 1, 0.5);

    RowHbox = gtk_hbox_new( FALSE, 10); gtk_widget_show( RowHbox);
    gtk_table_attach( GTK_TABLE( table1), RowHbox, 2, 3, r, r+1,( GtkAttachOptions)( GTK_FILL),( GtkAttachOptions)( GTK_FILL), 0, 0);

    hbox = gtk_hbox_new( FALSE, 2); gtk_widget_show( hbox);
    gtk_box_pack_start( GTK_BOX( RowHbox), hbox, FALSE, TRUE, 10);

    spinbutton1_adj = gtk_adjustment_new(  tle.ElemNum, 0, 9999, 1, 1, 0);
    GtkWidget *TleGenButton6 = gtk_spin_button_new( GTK_ADJUSTMENT( spinbutton1_adj), 0.001, 0);
    gtk_widget_show( TleGenButton6 );
    gtk_box_pack_start( GTK_BOX( hbox), TleGenButton6, FALSE, TRUE, 0);
    gtk_spin_button_set_numeric( GTK_SPIN_BUTTON(TleGenButton6), TRUE);
    g_signal_connect( G_OBJECT( TleGenButton6 ), "value_changed", G_CALLBACK( ChangeTLE ), GINT_TO_POINTER( 7 ) );
    ++r;




    /*
     * Section for Line 2 Entries
     *    01      Line Number of Element Data
     *    03-07   Satellite Number
     *    09-16   Inclination [Degrees]
     *    18-25   Right Ascension of the Ascending Node [Degrees]
     *    27-33   Eccentricity (decimal point assumed)
     *    35-42   Argument of Perigee [Degrees]
     *    44-51   Mean Anomaly [Degrees]
     *    53-63   Mean Motion [Revs per day]
     *    64-68   Revolution number at epoch [Revs]
     *    69      Checksum (Modulo 10)
     */
    label = gtk_label_new( _("<b>Line 2 Entries:</b>")); gtk_widget_show( label);
    gtk_table_attach( GTK_TABLE( table1), label, 0, 1, r, r+1,( GtkAttachOptions)( GTK_FILL),( GtkAttachOptions)( 0), 0, 0);
    gtk_label_set_use_markup( GTK_LABEL( label), TRUE);
    gtk_label_set_justify( GTK_LABEL( label), GTK_JUSTIFY_CENTER);
    gtk_misc_set_alignment( GTK_MISC( label), 1, 0.5);
    ++r;

    // Satellite Number
//    label = gtk_label_new( _("<b>Satellite Number:</b>")); gtk_widget_show( label);
//    gtk_table_attach( GTK_TABLE( table1), label, 1, 2, r, r+1,( GtkAttachOptions)( GTK_FILL),( GtkAttachOptions)( 0), 0, 0);
//    gtk_label_set_use_markup( GTK_LABEL( label), TRUE);
//    gtk_label_set_justify( GTK_LABEL( label), GTK_JUSTIFY_CENTER);
//    gtk_misc_set_alignment( GTK_MISC( label), 1, 0.5);
//
//    RowHbox = gtk_hbox_new( FALSE, 10); gtk_widget_show( RowHbox);
//    gtk_table_attach( GTK_TABLE( table1), RowHbox, 2, 3, r, r+1,( GtkAttachOptions)( GTK_FILL),( GtkAttachOptions)( GTK_FILL), 0, 0);
//
//    hbox = gtk_hbox_new( FALSE, 2); gtk_widget_show( hbox);
//    gtk_box_pack_start( GTK_BOX( RowHbox), hbox, FALSE, TRUE, 10);
//
//    spinbutton1_adj = gtk_adjustment_new(  99999, 0, 999999, 1, 1, 0);
//    GtkWidget *TleGenButton7 = gtk_spin_button_new( GTK_ADJUSTMENT( spinbutton1_adj), 0.001, 0);
//    gtk_widget_show( TleGenButton7 );
//    gtk_box_pack_start( GTK_BOX( hbox), TleGenButton7, FALSE, TRUE, 0);
//    gtk_spin_button_set_numeric( GTK_SPIN_BUTTON(TleGenButton7), TRUE);
//    //g_signal_connect( G_OBJECT( TleGenButton7 ), "value_changed", G_CALLBACK( ChangeViewQuat ), GINT_TO_POINTER( 0 ) );
//    ++r;

    // Inclination
    label = gtk_label_new( _("<b>Inclination [Degrees]:</b>")); gtk_widget_show( label);
    gtk_table_attach( GTK_TABLE( table1), label, 1, 2, r, r+1,( GtkAttachOptions)( GTK_FILL),( GtkAttachOptions)( 0), 0, 0);
    gtk_label_set_use_markup( GTK_LABEL( label), TRUE);
    gtk_label_set_justify( GTK_LABEL( label), GTK_JUSTIFY_CENTER);
    gtk_misc_set_alignment( GTK_MISC( label), 1, 0.5);

    RowHbox = gtk_hbox_new( FALSE, 10); gtk_widget_show( RowHbox);
    gtk_table_attach( GTK_TABLE( table1), RowHbox, 2, 3, r, r+1,( GtkAttachOptions)( GTK_FILL),( GtkAttachOptions)( GTK_FILL), 0, 0);

    hbox = gtk_hbox_new( FALSE, 2); gtk_widget_show( hbox);
    gtk_box_pack_start( GTK_BOX( RowHbox), hbox, FALSE, TRUE, 10);

    spinbutton1_adj = gtk_adjustment_new(  tle.Inc, 0.0, 180, 1.0, 10.0, 0);
    GtkWidget *TleGenButton8 = gtk_spin_button_new( GTK_ADJUSTMENT( spinbutton1_adj), 0.001, 4);
    gtk_widget_show( TleGenButton8 );
    gtk_box_pack_start( GTK_BOX( hbox), TleGenButton8, FALSE, TRUE, 0);
    gtk_spin_button_set_numeric( GTK_SPIN_BUTTON(TleGenButton8), TRUE);
    g_signal_connect( G_OBJECT( TleGenButton8 ), "value_changed", G_CALLBACK( ChangeTLE ), GINT_TO_POINTER( 8 ) );
    ++r;

    // RAAN
    label = gtk_label_new( _("<b>RAAN [Degrees]:</b>")); gtk_widget_show( label);
    gtk_table_attach( GTK_TABLE( table1), label, 1, 2, r, r+1,( GtkAttachOptions)( GTK_FILL),( GtkAttachOptions)( 0), 0, 0);
    gtk_label_set_use_markup( GTK_LABEL( label), TRUE);
    gtk_label_set_justify( GTK_LABEL( label), GTK_JUSTIFY_CENTER);
    gtk_misc_set_alignment( GTK_MISC( label), 1, 0.5);

    RowHbox = gtk_hbox_new( FALSE, 10); gtk_widget_show( RowHbox);
    gtk_table_attach( GTK_TABLE( table1), RowHbox, 2, 3, r, r+1,( GtkAttachOptions)( GTK_FILL),( GtkAttachOptions)( GTK_FILL), 0, 0);

    hbox = gtk_hbox_new( FALSE, 2); gtk_widget_show( hbox);
    gtk_box_pack_start( GTK_BOX( RowHbox), hbox, FALSE, TRUE, 10);

    spinbutton1_adj = gtk_adjustment_new(  tle.RAAN, 0.0, 360.0, 1.0, 10.0, 0);
    GtkWidget *TleGenButton9 = gtk_spin_button_new( GTK_ADJUSTMENT( spinbutton1_adj), 0.001, 4);
    gtk_widget_show( TleGenButton9 );
    gtk_box_pack_start( GTK_BOX( hbox), TleGenButton9, FALSE, TRUE, 0);
    gtk_spin_button_set_numeric( GTK_SPIN_BUTTON(TleGenButton9), TRUE);
    g_signal_connect( G_OBJECT( TleGenButton9 ), "value_changed", G_CALLBACK( ChangeTLE ), GINT_TO_POINTER( 9 ) );
    ++r;

    // Eccentricity
    label = gtk_label_new( _("<b>Eccentricity:</b>")); gtk_widget_show( label);
    gtk_table_attach( GTK_TABLE( table1), label, 1, 2, r, r+1,( GtkAttachOptions)( GTK_FILL),( GtkAttachOptions)( 0), 0, 0);
    gtk_label_set_use_markup( GTK_LABEL( label), TRUE);
    gtk_label_set_justify( GTK_LABEL( label), GTK_JUSTIFY_CENTER);
    gtk_misc_set_alignment( GTK_MISC( label), 1, 0.5);

    RowHbox = gtk_hbox_new( FALSE, 10); gtk_widget_show( RowHbox);
    gtk_table_attach( GTK_TABLE( table1), RowHbox, 2, 3, r, r+1,( GtkAttachOptions)( GTK_FILL),( GtkAttachOptions)( GTK_FILL), 0, 0);

    hbox = gtk_hbox_new( FALSE, 2); gtk_widget_show( hbox);
    gtk_box_pack_start( GTK_BOX( RowHbox), hbox, FALSE, TRUE, 10);

    spinbutton1_adj = gtk_adjustment_new(  tle.Ecc, 0.0, 0.9999999, 0.01, 0.10, 0);
    GtkWidget *TleGenButton10 = gtk_spin_button_new( GTK_ADJUSTMENT( spinbutton1_adj), 0.001, 7);
    gtk_widget_show( TleGenButton10 );
    gtk_box_pack_start( GTK_BOX( hbox), TleGenButton10, FALSE, TRUE, 0);
    gtk_spin_button_set_numeric( GTK_SPIN_BUTTON(TleGenButton10), TRUE);
    g_signal_connect( G_OBJECT( TleGenButton10 ), "value_changed", G_CALLBACK( ChangeTLE ), GINT_TO_POINTER( 10 ) );
    ++r;

    // Arg. Of Perigee
    label = gtk_label_new( _("<b>Arg. Of Perigee [Degrees]:</b>")); gtk_widget_show( label);
    gtk_table_attach( GTK_TABLE( table1), label, 1, 2, r, r+1,( GtkAttachOptions)( GTK_FILL),( GtkAttachOptions)( 0), 0, 0);
    gtk_label_set_use_markup( GTK_LABEL( label), TRUE);
    gtk_label_set_justify( GTK_LABEL( label), GTK_JUSTIFY_CENTER);
    gtk_misc_set_alignment( GTK_MISC( label), 1, 0.5);

    RowHbox = gtk_hbox_new( FALSE, 10); gtk_widget_show( RowHbox);
    gtk_table_attach( GTK_TABLE( table1), RowHbox, 2, 3, r, r+1,( GtkAttachOptions)( GTK_FILL),( GtkAttachOptions)( GTK_FILL), 0, 0);

    hbox = gtk_hbox_new( FALSE, 2); gtk_widget_show( hbox);
    gtk_box_pack_start( GTK_BOX( RowHbox), hbox, FALSE, TRUE, 10);

    spinbutton1_adj = gtk_adjustment_new(  tle.AoP, 0.0, 360.0, 1.0, 10.0, 0);
    GtkWidget *TleGenButton11 = gtk_spin_button_new( GTK_ADJUSTMENT( spinbutton1_adj), 0.001, 4);
    gtk_widget_show( TleGenButton11 );
    gtk_box_pack_start( GTK_BOX( hbox), TleGenButton11, FALSE, TRUE, 0);
    gtk_spin_button_set_numeric( GTK_SPIN_BUTTON(TleGenButton11), TRUE);
    g_signal_connect( G_OBJECT( TleGenButton11 ), "value_changed", G_CALLBACK( ChangeTLE ), GINT_TO_POINTER( 11 ) );
    ++r;

    // Mean Anomaly
    label = gtk_label_new( _("<b>Mean Anomaly [Degrees]:</b>")); gtk_widget_show( label);
    gtk_table_attach( GTK_TABLE( table1), label, 1, 2, r, r+1,( GtkAttachOptions)( GTK_FILL),( GtkAttachOptions)( 0), 0, 0);
    gtk_label_set_use_markup( GTK_LABEL( label), TRUE);
    gtk_label_set_justify( GTK_LABEL( label), GTK_JUSTIFY_CENTER);
    gtk_misc_set_alignment( GTK_MISC( label), 1, 0.5);

    RowHbox = gtk_hbox_new( FALSE, 10); gtk_widget_show( RowHbox);
    gtk_table_attach( GTK_TABLE( table1), RowHbox, 2, 3, r, r+1,( GtkAttachOptions)( GTK_FILL),( GtkAttachOptions)( GTK_FILL), 0, 0);

    hbox = gtk_hbox_new( FALSE, 2); gtk_widget_show( hbox);
    gtk_box_pack_start( GTK_BOX( RowHbox), hbox, FALSE, TRUE, 10);

    spinbutton1_adj = gtk_adjustment_new(  tle.MeanAnomaly, 0.0, 360.0, 1.0, 10.0, 0);
    GtkWidget *TleGenButton12 = gtk_spin_button_new( GTK_ADJUSTMENT( spinbutton1_adj), 0.001, 4);
    gtk_widget_show( TleGenButton12 );
    gtk_box_pack_start( GTK_BOX( hbox), TleGenButton12, FALSE, TRUE, 0);
    gtk_spin_button_set_numeric( GTK_SPIN_BUTTON(TleGenButton12), TRUE);
    g_signal_connect( G_OBJECT( TleGenButton12 ), "value_changed", G_CALLBACK( ChangeTLE ), GINT_TO_POINTER( 12 ) );
    ++r;

    // Mean Motion
    label = gtk_label_new( _("<b>Mean Motion [Revs/day]:</b>")); gtk_widget_show( label);
    gtk_table_attach( GTK_TABLE( table1), label, 1, 2, r, r+1,( GtkAttachOptions)( GTK_FILL),( GtkAttachOptions)( 0), 0, 0);
    gtk_label_set_use_markup( GTK_LABEL( label), TRUE);
    gtk_label_set_justify( GTK_LABEL( label), GTK_JUSTIFY_CENTER);
    gtk_misc_set_alignment( GTK_MISC( label), 1, 0.5);

    RowHbox = gtk_hbox_new( FALSE, 10); gtk_widget_show( RowHbox);
    gtk_table_attach( GTK_TABLE( table1), RowHbox, 2, 3, r, r+1,( GtkAttachOptions)( GTK_FILL),( GtkAttachOptions)( GTK_FILL), 0, 0);

    hbox = gtk_hbox_new( FALSE, 2); gtk_widget_show( hbox);
    gtk_box_pack_start( GTK_BOX( RowHbox), hbox, FALSE, TRUE, 10);

    spinbutton1_adj = gtk_adjustment_new(  tle.MeanMotion, 0.1, 20.0, 0.0001, 0.10, 0);
    GtkWidget *TleGenButton13 = gtk_spin_button_new( GTK_ADJUSTMENT( spinbutton1_adj), 0.001, 8);
    gtk_widget_show( TleGenButton13 );
    gtk_box_pack_start( GTK_BOX( hbox), TleGenButton13, FALSE, TRUE, 0);
    gtk_spin_button_set_numeric( GTK_SPIN_BUTTON(TleGenButton13), TRUE);
    g_signal_connect( G_OBJECT( TleGenButton13 ), "value_changed", G_CALLBACK( ChangeTLE ), GINT_TO_POINTER( 13 ) );
    ++r;

    // Rev Number at Epoch
    label = gtk_label_new( _("<b>Rev Number at Epoch [Revs]:</b>")); gtk_widget_show( label);
    gtk_table_attach( GTK_TABLE( table1), label, 1, 2, r, r+1,( GtkAttachOptions)( GTK_FILL),( GtkAttachOptions)( 0), 0, 0);
    gtk_label_set_use_markup( GTK_LABEL( label), TRUE);
    gtk_label_set_justify( GTK_LABEL( label), GTK_JUSTIFY_CENTER);
    gtk_misc_set_alignment( GTK_MISC( label), 1, 0.5);

    RowHbox = gtk_hbox_new( FALSE, 10); gtk_widget_show( RowHbox);
    gtk_table_attach( GTK_TABLE( table1), RowHbox, 2, 3, r, r+1,( GtkAttachOptions)( GTK_FILL),( GtkAttachOptions)( GTK_FILL), 0, 0);

    hbox = gtk_hbox_new( FALSE, 2); gtk_widget_show( hbox);
    gtk_box_pack_start( GTK_BOX( RowHbox), hbox, FALSE, TRUE, 10);

    spinbutton1_adj = gtk_adjustment_new(  tle.RevNum, 0, 99999, 1, 1, 0);
    GtkWidget *TleGenButton14 = gtk_spin_button_new( GTK_ADJUSTMENT( spinbutton1_adj), 0.001, 0);
    gtk_widget_show( TleGenButton14 );
    gtk_box_pack_start( GTK_BOX( hbox), TleGenButton14, FALSE, TRUE, 0);
    gtk_spin_button_set_numeric( GTK_SPIN_BUTTON(TleGenButton14), TRUE);
    g_signal_connect( G_OBJECT( TleGenButton14 ), "value_changed", G_CALLBACK( ChangeTLE ), GINT_TO_POINTER( 14 ) );
    ++r;



hbox = gtk_hbox_new( FALSE, 10); gtk_widget_show( hbox );
button = gtk_button_new_with_mnemonic( _("Trace Field Lines") ); gtk_widget_show( button );
gtk_box_pack_start( GTK_BOX(hbox), button, FALSE, TRUE, 10 );
g_signal_connect( G_OBJECT( button ), "clicked", G_CALLBACK( TLE_TraceFieldLines ), GINT_TO_POINTER( -1 ) );
gtk_box_pack_start( GTK_BOX(vbox2), hbox, FALSE, TRUE, 10 );







    /***********
     * General *
     ***********/
    vbox2 = gtk_vbox_new( FALSE, 0); gtk_widget_show( vbox2);
    gtk_container_set_border_width( GTK_CONTAINER( vbox2), 20);

    // Page title
    label = gtk_label_new (_("<b><span size=\"small\">General</span></b>")); gtk_widget_show (label);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_CENTER);
    gtk_label_set_use_markup( GTK_LABEL(label), TRUE );
    gtk_notebook_append_page( GTK_NOTEBOOK(notebook), vbox2, label );


    table1 = gtk_table_new (7, 2, FALSE); gtk_widget_show (table1);
    gtk_box_pack_start (GTK_BOX (vbox2), table1, TRUE, TRUE, 15);
    gtk_container_set_border_width (GTK_CONTAINER (table1), 20);
    gtk_table_set_row_spacings (GTK_TABLE (table1), 20);
    gtk_table_set_col_spacings (GTK_TABLE (table1), 10);


    /*
     * Row for nFieldLine Points
     */
    label = gtk_label_new (_("<b>Fieldline Pnts:</b>")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, 0, 1, 1, 2, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);
    gtk_label_set_use_markup (GTK_LABEL (label), TRUE);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_CENTER);
    gtk_misc_set_alignment (GTK_MISC (label), 1, 0.5);

    RowHbox = gtk_hbox_new (FALSE, 10); gtk_widget_show (RowHbox);
    gtk_table_attach (GTK_TABLE (table1), RowHbox, 1, 2, 1, 2, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (GTK_FILL), 0, 0);

    // nFieldline Points button
    hbox = gtk_hbox_new (FALSE, 5); gtk_widget_show (hbox);
    gtk_box_pack_start (GTK_BOX (RowHbox), hbox, FALSE, TRUE, 0);

    label = gtk_label_new ("nFieldPnts"); gtk_widget_show (label);
    gtk_box_pack_start (GTK_BOX (hbox), label, FALSE, FALSE, 0);
    gtk_label_set_use_markup (GTK_LABEL (label), TRUE);
    gtk_misc_set_alignment (GTK_MISC (label), 1, 0.5);

    spinbutton1_adj = gtk_adjustment_new(80, 2, 180, 1, 5, 0);
    GtkWidget *nFieldPntsButton = gtk_spin_button_new (GTK_ADJUSTMENT (spinbutton1_adj), 1, 0);
    gtk_widget_show (nFieldPntsButton);
    gtk_box_pack_start (GTK_BOX (hbox), nFieldPntsButton, FALSE, TRUE, 0);
    gtk_spin_button_set_numeric (GTK_SPIN_BUTTON (nFieldPntsButton), TRUE);
    g_signal_connect( G_OBJECT( nFieldPntsButton ), "value_changed", G_CALLBACK( ChangenFieldPnts ), NULL );



    /*
     * Row for Lighting
     */
    label = gtk_label_new (_("<b>Lighting:</b>")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, 0, 1, 2, 3, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);
    gtk_label_set_use_markup (GTK_LABEL (label), TRUE);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_CENTER);
    gtk_misc_set_alignment (GTK_MISC (label), 1, 0.5);

    RowHbox = gtk_hbox_new (FALSE, 10); gtk_widget_show (RowHbox);
    gtk_table_attach (GTK_TABLE (table1), RowHbox, 1, 2, 2, 3, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (GTK_FILL), 0, 0);

    // nFieldline Points button
    hbox = gtk_hbox_new (FALSE, 5); gtk_widget_show (hbox);
    gtk_box_pack_start (GTK_BOX (RowHbox), hbox, FALSE, TRUE, 0);

    RadioLightingStyle[0] = gtk_radio_button_new_with_mnemonic (NULL, _("Solar Illumination")); gtk_widget_show(RadioLightingStyle[0]);
    gtk_box_pack_start (GTK_BOX (hbox), RadioLightingStyle[0], FALSE, FALSE, 0);
    gtk_radio_button_set_group (GTK_RADIO_BUTTON(RadioLightingStyle[0]), RadioLightingGroup);
    RadioLightingGroup = gtk_radio_button_get_group (GTK_RADIO_BUTTON(RadioLightingStyle[0]));
    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON(RadioLightingStyle[0]), TRUE);
    g_signal_connect( G_OBJECT( RadioLightingStyle[0]), "toggled", G_CALLBACK( ChangeLighting ), GINT_TO_POINTER( 0 ) );

    RadioLightingStyle[2] = gtk_radio_button_new_with_mnemonic (NULL, _("Solar Illumination 2")); gtk_widget_show(RadioLightingStyle[2]);
    gtk_box_pack_start (GTK_BOX (hbox), RadioLightingStyle[2], FALSE, FALSE, 0);
    gtk_radio_button_set_group (GTK_RADIO_BUTTON(RadioLightingStyle[2]), RadioLightingGroup);
    RadioLightingGroup = gtk_radio_button_get_group (GTK_RADIO_BUTTON(RadioLightingStyle[2]));
    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON(RadioLightingStyle[2]), FALSE);
    g_signal_connect( G_OBJECT( RadioLightingStyle[2] ), "toggled", G_CALLBACK( ChangeLighting ), GINT_TO_POINTER( 2 ) );

    RadioLightingStyle[1] = gtk_radio_button_new_with_mnemonic (NULL, _("Fixed Illumination")); gtk_widget_show(RadioLightingStyle[1]);
    gtk_box_pack_start (GTK_BOX (hbox), RadioLightingStyle[1], FALSE, FALSE, 0);
    gtk_radio_button_set_group (GTK_RADIO_BUTTON(RadioLightingStyle[1]), RadioLightingGroup);
    RadioLightingGroup = gtk_radio_button_get_group (GTK_RADIO_BUTTON(RadioLightingStyle[1]));
    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON(RadioLightingStyle[1]), FALSE);
    g_signal_connect( G_OBJECT( RadioLightingStyle[1] ), "toggled", G_CALLBACK( ChangeLighting ), GINT_TO_POINTER( 1 ) );



    return( window1 );



}





int main( int argc, char *argv[] ) {

    int         i;
    Lgm_CTrans *c = Lgm_init_ctrans( 0 );
    char       Command[2048];

    
    /*
     * Initialize the ObjInfo structure
     */
    ObjInfo = Vds_InitObjectInfo();

    /*
     * Create the FLs and DSs from data in a file.
     */
    CreateFieldLinesAndDriftShells( "test.dat", ObjInfo );






    aInfo = New_aInfo();
    InitAtmosphere();

    InitSatSelectorInfo();



    StartYear  = 2017; StartMonth = 1; StartDay   = 1;
//StartYear  = 2009; StartMonth = 4; StartDay   = 20;
    StartDate  = StartYear*10000 + StartMonth*100 + StartDay;
    StartHour  = 6; StartMin = 0; StartSec   = 0; StartMilliSec = 0;
//StartHour  = 0; StartMin = 2; StartSec   = 30; StartMilliSec = 0;
    StartUT    = (double)StartHour + (double)StartMin/60.0 + (double)(StartSec+StartMilliSec/1000.0)/3600.0;
    StartJD    = Lgm_JD( StartYear, StartMonth, StartDay, StartUT, LGM_TIME_SYS_UTC, c );

    EndYear  = 2025; EndMonth = 12; EndDay   = 31;
    EndDate  = EndYear*10000 + EndMonth*100 + EndDay;
    EndHour  = 23; EndMin = 58; EndSec   = 30; EndMilliSec = 0;
    EndDate  = EndYear*10000 + EndMonth*100 + EndDay;
    EndUT    = (double)EndHour + (double)EndMin/60.0 + (double)(EndSec+EndMilliSec/1000.0)/3600.0;
    EndJD    = Lgm_JD( EndYear, EndMonth, EndDay, EndUT, LGM_TIME_SYS_UTC, c );

    CurrentYear  = StartYear; CurrentMonth = StartMonth; CurrentDay   = StartDay;
    CurrentDate  = CurrentYear*10000 + CurrentMonth*100 + CurrentDay;
    CurrentHour  = StartHour; CurrentMin = StartMin; CurrentSec   = StartSec; CurrentMilliSec = StartMilliSec;
    CurrentUT    = (double)CurrentHour + (double)CurrentMin/60.0 + (double)(CurrentSec+CurrentMilliSec/1000.0)/3600.0;
    CurrentJD    = Lgm_JD( CurrentYear, CurrentMonth, CurrentDay, CurrentUT, LGM_TIME_SYS_UTC, c );

    TimeInc     = 5.0;
    cFrame      = 0;
    nFramesLeft = (long int)((EndJD - CurrentJD) / (TimeInc/86400.0) + 1.5);
    if (nFramesLeft<0) nFramesLeft = 0;
    nFrames     = cFrame + nFramesLeft;
    RunTime     = TRUE;
    DumpFrames  = FALSE;
printf("EndJD, StartJD = %lf %lf\n", EndJD, StartJD);
printf("nFramesLeft, nFrames = %ld %ld\n", nFramesLeft, nFrames);

    mInfo  = Lgm_InitMagInfo( );
    UpdateTimeDepQuants( CurrentDate, CurrentUT );


    SpaceObjects = (_SpaceObjects *)calloc( 1, sizeof(_SpaceObjects) );

    MaxStarMagnitude = 8.0;

    for (i=0; i<9; i++) QuadsLoaded[i] = -1;




    for (i=0; i<10; i++) ShowSatellites[i] = 0;

    gInfo = (GuiInfo *)calloc(1, sizeof(GuiInfo));

    Lgm_free_ctrans( c );

    /*
     *  Initialize GTK
     */
    gtk_set_locale( );
    gtk_init( &argc, &argv );
    //add_pixmap_directory( PACKAGE_DATA_DIR "/" PACKAGE "/pixmaps" );
    add_pixmap_directory( "/home/mgh/Desktop" );
    add_pixmap_directory( "." );
    add_pixmap_directory( "/home/mgh/DREAM/Dream/Dream/Images" );



    for (i=0; i<ObjInfo->nPitchAngles+1; i++){
        ShowPitchAngle[i] = 0;
        ShowPitchAngle2[i] = 0;
    }


    /*
     * Set up materials
     */
    LGM_ARRAY_1D( gInfo->FieldLineMaterial,                     ObjInfo->nPitchAngles+1,    MaterialProp );
    LGM_ARRAY_1D( gInfo->FieldLineShininessButton,              ObjInfo->nPitchAngles+1,    GtkWidget * );
    LGM_ARRAY_1D( gInfo->FieldLineShininessButtonHandler,       ObjInfo->nPitchAngles+1,    gulong );
    LGM_ARRAY_1D( gInfo->FieldLineMaterialButton,               ObjInfo->nPitchAngles+1,    GtkWidget * );
    LGM_ARRAY_1D( gInfo->FieldLineMaterialButtonHandler,        ObjInfo->nPitchAngles+1,    gulong );
    LGM_ARRAY_1D( gInfo->FieldLineShowPitchAngleButton,         ObjInfo->nPitchAngles+1,    GtkWidget * );
    LGM_ARRAY_1D( gInfo->FieldLineShowPitchAngleButtonHandler,  ObjInfo->nPitchAngles+1,    gulong );
    LGM_ARRAY_1D( gInfo->FieldLineDiffuseColorButton,           ObjInfo->nPitchAngles+1,    GtkWidget * );
    LGM_ARRAY_1D( gInfo->FieldLineAmbientColorButton,           ObjInfo->nPitchAngles+1,    GtkWidget * );
    LGM_ARRAY_1D( gInfo->FieldLineSpecularColorButton,          ObjInfo->nPitchAngles+1,    GtkWidget * );

    LGM_ARRAY_1D( gInfo->DriftShellMaterial,                    ObjInfo->nPitchAngles+1,    MaterialProp );
    LGM_ARRAY_1D( gInfo->DriftShellShininessButton,             ObjInfo->nPitchAngles+1,    GtkWidget * );
    LGM_ARRAY_1D( gInfo->DriftShellShininessButtonHandler,      ObjInfo->nPitchAngles+1,    gulong );
    LGM_ARRAY_1D( gInfo->DriftShellMaterialButton,              ObjInfo->nPitchAngles+1,    GtkWidget * );
    LGM_ARRAY_1D( gInfo->DriftShellMaterialButtonHandler,       ObjInfo->nPitchAngles+1,    gulong );
    LGM_ARRAY_1D( gInfo->DriftShellShowPitchAngleButton,        ObjInfo->nPitchAngles+1,    GtkWidget * );
    LGM_ARRAY_1D( gInfo->DriftShellShowPitchAngleButtonHandler, ObjInfo->nPitchAngles+1,    gulong );
    LGM_ARRAY_1D( gInfo->DriftShellDiffuseColorButton,          ObjInfo->nPitchAngles+1,    GtkWidget * );
    LGM_ARRAY_1D( gInfo->DriftShellAmbientColorButton,          ObjInfo->nPitchAngles+1,    GtkWidget * );
    LGM_ARRAY_1D( gInfo->DriftShellSpecularColorButton,         ObjInfo->nPitchAngles+1,    GtkWidget * );

    int ii;
    for (i=0; i<ObjInfo->nPitchAngles; i++){

        ii = i%19;

        gInfo->FieldLineMaterial[i] = mat_silver;
        gInfo->FieldLineMaterial[i].diffuse[0] = colors[ii][0];
        gInfo->FieldLineMaterial[i].diffuse[1] = colors[ii][1];
        gInfo->FieldLineMaterial[i].diffuse[2] = colors[ii][2];
        gInfo->FieldLineMaterial[i].diffuse[3] = 1.0;

        gInfo->DriftShellMaterial[i] = mat_pearl;
        gInfo->DriftShellMaterial[i].diffuse[0] = colors[ii][0];
        gInfo->DriftShellMaterial[i].diffuse[1] = colors[ii][1];
        gInfo->DriftShellMaterial[i].diffuse[2] = colors[ii][2];
        gInfo->DriftShellMaterial[i].diffuse[3] = 0.6;

    }


    /*
     * SeT default TLE params
     */
    strcpy( tle.Name, "AMPTE_CCE" );

    tle.SatNum = 15199;
    strcpy( tle.IntDesig, "84088A" );
    tle.Epoch = 17078.89650569;
    tle.d1MeanMotion = 0.00000110;
    tle.d2MeanMotion = 0.0;
    tle.BSTAR = 0.0;
    tle.ElemNum = 262;

    tle.Inc = 10.0;
    tle.RAAN = 180.0;
    tle.Ecc = 0.495323;
    tle.AoP = 120.0;
    tle.MeanAnomaly = 10.0000;
    tle.MeanMotion = 01.33;
    tle.RevNum = 331;


    sprintf( Command, "mkdir -p %s/SAT_GROUPS", getenv("HOME") );
    system( Command );



    /*
     * Create the main window
     */
    create_ViewDriftShell( NULL );








    PitchAngleDisplayProperties();

    gtk_main( );

    Lgm_FreeMagInfo( mInfo );


    return( 0 );


}
