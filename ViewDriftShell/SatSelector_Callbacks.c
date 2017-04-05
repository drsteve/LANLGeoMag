#include "SatSelector.h"

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

        case 0:
                Group->SatRed_OverRide = Color.red/65535.0;
                Group->SatGrn_OverRide = Color.green/65535.0;
                Group->SatBlu_OverRide = Color.blue/65535.0;
                Group->SatAlf_OverRide = Alpha/65535.0;
                break;

    }
            
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
