#include "SatSelector.h"

extern float IllumFL_ka;
extern float IllumFL_kd;
extern float IllumFL_ks;
extern double IllumFL_n;
extern double IllumFL_w;
extern GLuint  Texture_Fd;
extern GLuint  Texture_Fs;

extern GtkWidget *PUKE_SATSEL_VBOX;


void ChangeSatColor( GtkWidget  *button, gpointer data ) {
    int i, j, GroupNumber, ButtonNumber;
     _GroupNode *g;
    _SpaceObjects   *Group;
    GdkColor    Color;
    guint16     Alpha;

    j = GPOINTER_TO_INT( data );
    GroupNumber = j/100;
    ButtonNumber = j - GroupNumber*100;

    // Find out which group called us
    i = 0; g = SatSelectorInfo->SatGroupList;
    while ( (g->Next != NULL) && (GroupNumber > i) ) {
        g = g->Next; ++i;
    }
    Group = g->Group;

    /*
     *  Get color
     */
    gtk_color_button_get_color( GTK_COLOR_BUTTON(button), &Color );
    Alpha = gtk_color_button_get_alpha( GTK_COLOR_BUTTON(button) );
    

    switch ( ButtonNumber ) {

        // For Sat or Single position group of buttons
        case 11:
                Group->SatRed_OverRide = Color.red/65535.0;
                Group->SatGrn_OverRide = Color.green/65535.0;
                Group->SatBlu_OverRide = Color.blue/65535.0;
                Group->SatAlf_OverRide = Alpha/65535.0;
                break;
        case 12:
                Group->ssglRed_OverRide = Color.red/65535.0;
                Group->ssglGrn_OverRide = Color.green/65535.0;
                Group->ssglBlu_OverRide = Color.blue/65535.0;
                Group->ssglAlf_OverRide = Alpha/65535.0;
                break;
        case 13:
                Group->ssflRed_OverRide = Color.red/65535.0;
                Group->ssflGrn_OverRide = Color.green/65535.0;
                Group->ssflBlu_OverRide = Color.blue/65535.0;
                Group->ssflAlf_OverRide = Alpha/65535.0;
                break;
        case 14:
                Group->ssfpRed_OverRide = Color.red/65535.0;
                Group->ssfpGrn_OverRide = Color.green/65535.0;
                Group->ssfpBlu_OverRide = Color.blue/65535.0;
                Group->ssfpAlf_OverRide = Alpha/65535.0;
                break;


        // For Orbit group of buttons
        case 20:
                Group->oRed_OverRide = Color.red/65535.0;
                Group->oGrn_OverRide = Color.green/65535.0;
                Group->oBlu_OverRide = Color.blue/65535.0;
                Group->oAlf_OverRide = Alpha/65535.0;
                break;
        case 21:
                Group->ogpRed_OverRide = Color.red/65535.0;
                Group->ogpGrn_OverRide = Color.green/65535.0;
                Group->ogpBlu_OverRide = Color.blue/65535.0;
                Group->ogpAlf_OverRide = Alpha/65535.0;
                break;
        case 22:
                Group->oglRed_OverRide = Color.red/65535.0;
                Group->oglGrn_OverRide = Color.green/65535.0;
                Group->oglBlu_OverRide = Color.blue/65535.0;
                Group->oglAlf_OverRide = Alpha/65535.0;
                break;
        case 23:
                Group->oflRed_OverRide = Color.red/65535.0;
                Group->oflGrn_OverRide = Color.green/65535.0;
                Group->oflBlu_OverRide = Color.blue/65535.0;
                Group->oflAlf_OverRide = Alpha/65535.0;
                break;
        case 24:
                Group->ofpRed_OverRide = Color.red/65535.0;
                Group->ofpGrn_OverRide = Color.green/65535.0;
                Group->ofpBlu_OverRide = Color.blue/65535.0;
                Group->ofpAlf_OverRide = Alpha/65535.0;
                break;



        // For Streak group of buttons
        case 30:
                Group->sRed_OverRide = Color.red/65535.0;
                Group->sGrn_OverRide = Color.green/65535.0;
                Group->sBlu_OverRide = Color.blue/65535.0;
                Group->sAlf_OverRide = Alpha/65535.0;
                break;
        case 31:
                Group->sglRed_OverRide = Color.red/65535.0;
                Group->sglGrn_OverRide = Color.green/65535.0;
                Group->sglBlu_OverRide = Color.blue/65535.0;
                Group->sglAlf_OverRide = Alpha/65535.0;
                break;
        case 32:
                Group->sflRed_OverRide = Color.red/65535.0;
                Group->sflGrn_OverRide = Color.green/65535.0;
                Group->sflBlu_OverRide = Color.blue/65535.0;
                Group->sflAlf_OverRide = Alpha/65535.0;
                break;
        case 33:
                Group->sfpRed_OverRide = Color.red/65535.0;
                Group->sfpGrn_OverRide = Color.green/65535.0;
                Group->sfpBlu_OverRide = Color.blue/65535.0;
                Group->sfpAlf_OverRide = Alpha/65535.0;
                break;
        case 34:
                Group->sfpRed_OverRide = Color.red/65535.0;
                Group->sfpGrn_OverRide = Color.green/65535.0;
                Group->sfpBlu_OverRide = Color.blue/65535.0;
                Group->sfpAlf_OverRide = Alpha/65535.0;
                break;

    }
            
    ReCreateSats();
    ReCreateSatOrbits();
    expose_event( drawing_area, NULL, NULL );

}


/*
 * Callback for orbit check button options
 */
void ToggleOrbitOptions( GtkWidget  *w, unsigned int *data ) {
    int i, j, State, GroupNumber, CheckButtonNumber;
     _GroupNode *g;
    _SpaceObjects   *Group;

    j = GPOINTER_TO_INT( data );
    GroupNumber = j/100;
    CheckButtonNumber = j - GroupNumber*100;

    State = gtk_toggle_button_get_active(  GTK_TOGGLE_BUTTON( w ) );


    // Find out which group called us
    i = 0; g = SatSelectorInfo->SatGroupList;
    while ( (g->Next != NULL) && (GroupNumber > i) ) {
        g = g->Next; ++i;
    }
    Group = g->Group;

    switch ( CheckButtonNumber ) {

        case 0:
                // Overall Draw Group Button
                Group->DrawGroup = State;
                break;




        case 10:
                // Draw Label Button
                for (i=0; i<Group->nSat; i++) Group->Sat[i].DrawLabel = State;
                break;
        case 11:
                // Draw Current Position Button
                Group->DrawPosition = State;
                break;
        case 12:
                // Draw Sat to Ground Lines Button
                for (i=0; i<Group->nSat; i++) Group->Sat[i].DrawSatToGroundLine = State;
                break;
        case 13:
                // Draw Sat Field Lines Button
                for (i=0; i<Group->nSat; i++) Group->Sat[i].DrawSatFieldLines = State;
                break;
        case 14:
                // Draw Sat Field Line Footpoints Button
                for (i=0; i<Group->nSat; i++) Group->Sat[i].DrawSatFLFootpoints = State;
                break;
        case 15:
                // Draw point sprites Button
                Group->DrawAsPointSprites = State;
                break;
        case 16:
                // Override Group Colors
                Group->OverRideSatColors = State;
                break;



        case 20:
                // Draw Orbit
                for (i=0; i<Group->nSat; i++) Group->Sat[i].DrawOrbit = State;
                break;
        case 21:
                // Draw Orbit Ground Path
                for (i=0; i<Group->nSat; i++) Group->Sat[i].DrawGroundPathOfOrbit = State;
                break;
        case 22:
                // Draw Orbit to Ground lines
                for (i=0; i<Group->nSat; i++) Group->Sat[i].DrawOrbitToGroundLines = State;
                break;
        case 23:
                // Draw Orbit Field Lines
                for (i=0; i<Group->nSat; i++) Group->Sat[i].DrawOrbitFieldLines = State;
                break;
        case 24:
                // Draw Orbit Field Line Footpoints
                for (i=0; i<Group->nSat; i++) Group->Sat[i].DrawOrbitFLFootpoints = State;
                break;
        case 26:
                // Override Group Colors
                Group->OverRideOrbitColors = State;
                break;


        case 30:
                for (i=0; i<Group->nSat; i++) Group->Sat[i].DrawStreak = State;
                break;
        case 31:
                for (i=0; i<Group->nSat; i++) Group->Sat[i].DrawGroundPathOfStreak = State;
                break;
        case 32:
                for (i=0; i<Group->nSat; i++) Group->Sat[i].DrawStreakToGroundLines = State;
                break;
        case 33:
                for (i=0; i<Group->nSat; i++) Group->Sat[i].DrawStreakFieldLines = State;
                break;
        case 34:
                for (i=0; i<Group->nSat; i++) Group->Sat[i].DrawStreakFLFootpoints = State;
                break;
        case 36:
                // Override Group Colors
                Group->OverRideStreakColors = State;
                break;







//        case 3:
//                Group->DrawSatellites = State;
//                // loop through all Sats in this group and set the Draw flags appropriately
//                for (i=0; i<Group->nSat; i++) {
//                    if ( !strstr(Group->Sat[i].TLE.Name, " R/B") && !strstr(Group->Sat[i].TLE.Name, " DEB") ) {
//                        Group->Sat[i].Draw = State;
//                    }
//                }
//                break;
//        case 4:
//                Group->DrawRocketBodies = State;
//                // loop through all Sats in this group and set the Draw flags appropriately
//                for (i=0; i<Group->nSat; i++) {
//                    if ( strstr(Group->Sat[i].TLE.Name, " R/B") ) {
//                        Group->Sat[i].Draw = State;
//                    }
//                }
//                break;
//        case 5:
//                Group->DrawDebris = State;
//                // loop through all Sats in this group and set the Draw flags appropriately
//                for (i=0; i<Group->nSat; i++) {
//                    if ( strstr(Group->Sat[i].TLE.Name, " DEB") ) {
//                        Group->Sat[i].Draw = State;
//                    }
//                }
//                break;








    }

    ReCreateSats();
    ReCreateSatOrbits();
    expose_event( drawing_area, NULL, NULL );
}



/*
 * Callback for changing Illuminatehd Field Line mShader params, ka, kd, ks
 */
void ChangeIllumFLParams( GtkWidget  *w, unsigned int *data ) {

    int     SpinButtonNumber;
    double  Value;
    float   *FdImage;
    float   *FsImage;

    SpinButtonNumber = GPOINTER_TO_INT( data );
    Value = gtk_spin_button_get_value(  GTK_SPIN_BUTTON( w ) );


    switch ( SpinButtonNumber ) {

        case 1:
                /*
                 *  Change ka value
                 */
                IllumFL_ka = Value;
                break;

        case 2:
                /*
                 *  Change kd value
                 */
                IllumFL_kd = Value;
                break;

        case 3:
                /*
                 *  Change ks value
                 */
                IllumFL_ks = Value;
                break;

        case 4:
                /*
                 *  Change n value
                 */
                IllumFL_n = Value;
                if ( GenIllumTextures( 256, IllumFL_n, &FdImage, &FsImage ) ) {

                    glGenTextures( 1, &Texture_Fd );
                    glBindTexture( GL_TEXTURE_2D, Texture_Fd );
                    glTexImage2D( GL_TEXTURE_2D, 0, GL_R32F, 256, 256, 0, GL_RED,  GL_FLOAT, FdImage );
                    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
                    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
                    glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE );
                    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
                    free( FdImage );

                    glGenTextures( 1, &Texture_Fs );
                    glBindTexture( GL_TEXTURE_2D, Texture_Fs );
                    glTexImage2D( GL_TEXTURE_2D, 0, GL_R32F, 256, 256, 0, GL_RED,  GL_FLOAT, FsImage );
                    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
                    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
                    glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE );
                    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
                    free( FsImage );

                }
                break;

        case 5:
                /*
                 *  Change ks value
                 */
                IllumFL_w = Value;
                break;


    }



    ReCreateSats();
    ReCreateSatOrbits();
    DrawScene();
    expose_event( drawing_area, NULL, NULL );

        

}





/*
 * Callback for orbit spin  button options
 */
void ChangeOrbitOptions( GtkWidget  *w, unsigned int *data ) {
    int     i, j, GroupNumber, SpinButtonNumber;
    double  Value;
     _GroupNode *g;
    _SpaceObjects   *Group;

    j = GPOINTER_TO_INT( data );
    GroupNumber = j/100;
    SpinButtonNumber = j - GroupNumber*100;

    Value = gtk_spin_button_get_value(  GTK_SPIN_BUTTON( w ) );


    // Find out which group called us
    i = 0; g = SatSelectorInfo->SatGroupList;
    while ( (g->Next != NULL) && (GroupNumber > i) ) {
        g = g->Next; ++i;
    }
    Group = g->Group;

    switch ( SpinButtonNumber ) {

        case 25:
                /*
                 *  Change oPeriodFrac value
                 */
                for (i=0; i<Group->nSat; i++) Group->Sat[i].oPeriodFrac = Value;
                break;
        case 35:
                /*
                 *  Change sPeriodFrac value
                 */
                for (i=0; i<Group->nSat; i++) Group->Sat[i].sPeriodFrac = Value;
                break;
    }

    ReCreateSats();
    ReCreateSatOrbits();
    expose_event( drawing_area, NULL, NULL );
}


static void ToggleSatellites( GtkMenuItem  *menuitem, const GLuint *data ) {
    int i;
    i = GPOINTER_TO_INT( data );
    ShowSatellites[i] = !ShowSatellites[i];
    ReLoadSpaceObjects();
//    ReLoadTLEs( );
    ReCreateSats();
    ReCreateSatOrbits();
    printf("19.\n"); expose_event( drawing_area, NULL, NULL );
}
