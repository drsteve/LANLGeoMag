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
                Group->SatRed = Color.red/65535.0;
                Group->SatGrn = Color.green/65535.0;
                Group->SatBlu = Color.blue/65535.0;
                Group->SatAlf = Alpha/65535.0;
                break;
        case 1:
                Group->RbRed = Color.red/65535.0;
                Group->RbGrn = Color.green/65535.0;
                Group->RbBlu = Color.blue/65535.0;
                Group->RbAlf = Alpha/65535.0;
                break;
        case 2:
                Group->DebRed = Color.red/65535.0;
                Group->DebGrn = Color.green/65535.0;
                Group->DebBlu = Color.blue/65535.0;
                Group->DebAlf = Alpha/65535.0;
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
//                if ( State ) ReLoadSpaceObjects();
//gtk_widget_destroy( PUKE_SATSEL_VBOX);


                Group->DrawGroup = State;
                break;
        case 1:
                for (i=0; i<Group->nSat; i++) Group->Sat[i].DrawLabel = State;
                break;
        case 2:
                for (i=0; i<Group->nSat; i++) Group->Sat[i].DrawSatToGroundLine = State;
                break;
        case 3:
                Group->DrawSatellites = State;
                // loop through all Sats in this group and set the Draw flags appropriately
                for (i=0; i<Group->nSat; i++) {
                    if ( !strstr(Group->Sat[i].TLE.Name, " R/B") && !strstr(Group->Sat[i].TLE.Name, " DEB") ) {
                        Group->Sat[i].Draw = State;
                    }
                }
                break;
        case 4:
                Group->DrawRocketBodies = State;
                // loop through all Sats in this group and set the Draw flags appropriately
                for (i=0; i<Group->nSat; i++) {
                    if ( strstr(Group->Sat[i].TLE.Name, " R/B") ) {
                        Group->Sat[i].Draw = State;
                    }
                }
                break;
        case 5:
                Group->DrawDebris = State;
                // loop through all Sats in this group and set the Draw flags appropriately
                for (i=0; i<Group->nSat; i++) {
                    if ( strstr(Group->Sat[i].TLE.Name, " DEB") ) {
                        Group->Sat[i].Draw = State;
                    }
                }
                break;
        case 6:
                for (i=0; i<Group->nSat; i++) Group->Sat[i].DrawOrbit = State;
                break;
        case 7:
                for (i=0; i<Group->nSat; i++) Group->Sat[i].DrawGroundPathOfOrbit = State;
                break;
        case 8:
                for (i=0; i<Group->nSat; i++) Group->Sat[i].DrawOrbitToGroundLines = State;
                break;
        case 9:
                for (i=0; i<Group->nSat; i++) Group->Sat[i].DrawStreak = State;
                break;
        case 10:
                for (i=0; i<Group->nSat; i++) Group->Sat[i].DrawGroundPathOfStreak = State;
                break;
        case 11:
                for (i=0; i<Group->nSat; i++) Group->Sat[i].DrawStreakToGroundLines = State;
                break;
        case 12:
                Group->DrawAsPointSprites = State;
                break;
        case 13:
                Group->DrawPosition = State;
                break;
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

        case 0:
                /*
                 *  Change oPeriodFrac value
                 */
                for (i=0; i<Group->nSat; i++) Group->Sat[i].oPeriodFrac = Value;
                break;
        case 1:
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
