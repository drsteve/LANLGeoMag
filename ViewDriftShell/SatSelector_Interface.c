#include "SatSelector.h"


void ToggleOrbitOptions( GtkWidget  *w, unsigned int *data );
void ChangeOrbitOptions( GtkWidget  *w, unsigned int *data );
void ChangeSatColor( GtkWidget  *button, gpointer data );

int SatGroupFilter( const struct dirent *d ) {

    int     len = strlen( d->d_name );

    if ( len < 4 ) return(0);
    if ( strstr( ".txt", &(d->d_name[len-4]) ) != NULL ) return(1);

    return(0);

}


int LoadSpaceObjects( ){

    int             i, j;
    double          r, g, b, mag;
    _SgpTLE         *TLEs;
    _SpaceObjects   *Group;
    int             nTLEs=0;
    char            Filename[1024];
    struct dirent **NameList;
    int             nFiles;
    char            SAT_GRPS_DIR[2000];

    sprintf( SAT_GRPS_DIR, "%s/SAT_GROUPS", getenv("HOME") );
    
    /*
     * Read in TLEs from the SatelliteGroups directory
     */
    if ( (nFiles = scandir( SAT_GRPS_DIR, &NameList, SatGroupFilter, alphasort )) < 0 ) {

        printf( "Problem reading Sat Groups Directory: %s\n", SAT_GRPS_DIR );
        return(-1);
    
    } else {

        // temporarily allocate enough elements to be able to read in any file...
        TLEs = (_SgpTLE *)calloc( 20000, sizeof(_SgpTLE) );
        for (i=0; i<nFiles; i++){
            nTLEs = 0;
            sprintf( Filename, "%s/%s", SAT_GRPS_DIR, NameList[i]->d_name );
            LgmSgp_ReadTlesFromFile( Filename, &nTLEs, TLEs, 0 );
            printf("File = %s nTLEs found=%d\n", Filename, nTLEs);


            /*
             *  Alloc mem for group
             */
            Group       = (_SpaceObjects *) calloc( 1, sizeof( _SpaceObjects ) );
            Group->nSat = nTLEs;
            Group->Sat  = (_SpaceObjectItem *) calloc( nTLEs, sizeof( _SpaceObjectItem ) );
            Group->DrawGroup          = FALSE;
            Group->DrawSatellites     = TRUE;
            Group->DrawRocketBodies   = FALSE;
            Group->DrawDebris         = FALSE;
            Group->DrawAsPointSprites = TRUE;
            Group->DrawPosition       = TRUE;
            Group->SatRed = 1.0; Group->SatGrn = 1.0; Group->SatBlu = 1.0; Group->SatAlf = 0.5*0.75;
            Group->RbRed  = 1.0; Group->RbGrn  = 1.0; Group->RbBlu  = 1.0; Group->RbAlf  = 0.5*0.75;
            Group->DebRed = 1.0; Group->DebGrn = 1.0; Group->DebBlu = 1.0; Group->DebAlf = 0.5*0.75;

            /*
             *  Fill in data for each item in the group
             */
            for (j=0; j<nTLEs; j++){
                // copy over the TLE
                Group->Sat[j].TLE = TLEs[j];

                if ( !strstr(Group->Sat[j].TLE.Name, " R/B") && !strstr(Group->Sat[j].TLE.Name, " DEB") ) Group->Sat[j].Draw = Group->DrawSatellites;
                if ( strstr(Group->Sat[j].TLE.Name, " R/B") ) Group->Sat[j].Draw = Group->DrawRocketBodies;
                if ( strstr(Group->Sat[j].TLE.Name, " DEB") ) Group->Sat[j].Draw = Group->DrawDebris;


                Group->Sat[j].DrawStreak              = TRUE;
                Group->Sat[j].sPeriodFrac             = 25.0;
                Group->Sat[j].DrawGroundPathOfStreak  = FALSE;
                Group->Sat[j].DrawStreakToGroundLines = FALSE;

                Group->Sat[j].DrawOrbit               = TRUE;
                Group->Sat[j].oPeriodFrac             = 100.0;
                Group->Sat[j].DrawGroundPathOfOrbit   = FALSE;
                Group->Sat[j].DrawOrbitToGroundLines  = FALSE;

                Group->Sat[j].DrawSatToGroundLine     = FALSE;

                Group->Sat[j].DrawLabel               = TRUE;


                // Streak Colors
                r = (double)rand()/(double)RAND_MAX;
                g = (double)rand()/(double)RAND_MAX;
                b = (double)rand()/(double)RAND_MAX;
                mag = sqrt(r*r+g*g+b*b); r /= mag; g /= mag; b /= mag;
                Group->Sat[j].sRed = r;
                Group->Sat[j].sGrn = g;
                Group->Sat[j].sBlu = b;
                Group->Sat[j].sAlf = 0.75;

                Group->Sat[j].sgpRed = 0.75*r;
                Group->Sat[j].sgpGrn = 0.75*g;
                Group->Sat[j].sgpBlu = 0.75*b;
                Group->Sat[j].sgpAlf = 0.75*0.75;

                Group->Sat[j].sglRed = 0.75*r;
                Group->Sat[j].sglGrn = 0.75*g;
                Group->Sat[j].sglBlu = 0.75*b;
                Group->Sat[j].sglAlf = 0.75*0.75;

                // Orbit Colors
                Group->Sat[j].oRed = 0.8;
                Group->Sat[j].oGrn = 0.8;
                Group->Sat[j].oBlu = 0.8;
                Group->Sat[j].oAlf = 0.20;

                Group->Sat[j].ogpRed = 0.2;
                Group->Sat[j].ogpGrn = 0.2;
                Group->Sat[j].ogpBlu = 0.2;
                Group->Sat[j].ogpAlf = 0.8;

                Group->Sat[j].oglRed = 0.7;
                Group->Sat[j].oglGrn = 0.7;
                Group->Sat[j].oglBlu = 0.7;
                Group->Sat[j].oglAlf = 0.1;

                // Single Line
                Group->Sat[j].ssglRed = r;
                Group->Sat[j].ssglGrn = g;
                Group->Sat[j].ssglBlu = b;
                Group->Sat[j].ssglAlf = 0.5*0.75;

            }

            /*
             *  Add this group into the SatGroupList
             */
            SatSelectorInfo->nSatGroups += 1;
            Add_SatGroupNode( &SatSelectorInfo->SatGroupList, NameList[i]->d_name, Group );

            
        }
        free( TLEs );

    }

    return(1);

}


void ReLoadSpaceObjects( ) {
    if ( SpaceObjects->Sat ) free( SpaceObjects->Sat );
    LoadSpaceObjects( );
}




int Delete_SatGroupNode( _GroupNode **List, char *GroupName ) {

    _GroupNode *p;
    int        flag;

    /*
     * Find node to delete
     */
    if ( List == NULL ) {

        // List is already empty
        return(-1);

    } else {

        p = *List;
        flag = 0;
        while ( flag == 0 ) {

            if ( strcmp( GroupName, p->GroupName ) == 0 ) {
                // found a match
                flag = 1;
            } else if ( p->Next == NULL ) {
                // return -- Group doesnt exist
                return(-1);
            } else {
                // Increment pointer to next node in list
                p = p->Next;
            }

        }

        // Remove node from list
        if ( (p->Prev == NULL) && (p->Next == NULL) ) {
            // The node is the only one
            *List = NULL;
        } else if ( p->Prev == NULL ) {
            // The node is the start 
            p->Next->Prev = NULL;
            *List = p;
        } else if ( p->Next == NULL ) {
            // The node is the end 
            p->Prev->Next = NULL;
        } else {
            // The node is interior
            p->Prev->Next = p->Next;
            p->Next->Prev = p->Prev;
        }
        //FreeGroup( p->Group );
        free( p->GroupName );
        free( p );
        return(1);
            
    }

}


void Add_SatGroupNode( _GroupNode **List, char *GroupName, _SpaceObjects *Group ) {

    _GroupNode *NewNode, *p;
    int        flag;

    // Alloc memory for a new node
    NewNode = (_GroupNode *) calloc( 1, sizeof( _GroupNode ) );

    // Alloc memory for GroupName string
    NewNode->GroupName = (char *) calloc( strlen(GroupName)+1, sizeof( char ) );
    strcpy( NewNode->GroupName, GroupName );


    NewNode->Group = Group;

    /*
     * Find place in list to insert it
     */
    if ( *List == NULL ) {

        // Nothing in the list yet -- just insert at start.
        *List = NewNode;
        NewNode->Prev = NULL;
        NewNode->Next = NULL;

    } else {

        // Insert in alphabetical order
        p = *List;
        flag = 0;
        while ( flag == 0 ) {

            if ( strcmp( GroupName, p->GroupName ) <= 0 ) {
                // Bail out and indicate that we need to insert before p
                flag = -1;
            } else if ( p->Next == NULL ) {
                // Bail out and indicate that we need to insert after p
                // -- this will be at end of list
                flag = 1;
            } else {
                // Increment pointer to next node in list
                p = p->Next;
            }

        }

        if ( flag < 0 ) {
            // insert before p
            NewNode->Next = p;
            NewNode->Prev = p->Prev;
            if ( p->Prev == NULL ) {
                //NewNode becomes first on the list
                p->Prev = NewNode;
                *List    = NewNode;
            } else {
                //NewNode is inserted between two existing nodes
                p->Prev->Next = NewNode;
                p->Prev       = NewNode;
            }
        } else {
            // insert at end of list p should be last node in current list
            NewNode->Next = NULL;
            NewNode->Prev = p;
            p->Next = NewNode;
            
        }

    }


}



_SpaceObjects *New_SpaceObjects() {
    _SpaceObjects  *s;
    if ( (s = (_SpaceObjects *) calloc( 1, sizeof( _SpaceObjects ) )) == NULL ) {
        printf("New_SpaceObjects(): Could not allocate memory for _SpaceObjects structure\n");
    }
    s->nSat = 0;
    return( s );
}

_SatSelectorInfo *New_SatSelectorInfo() {
    _SatSelectorInfo  *s;
    if ( (s = (_SatSelectorInfo *) calloc( 1, sizeof( _SatSelectorInfo ) )) == NULL ) {
        printf("New_SatSelectorInfo(): Could not allocate memory for _SatSelectorInfo structure\n");
    }
    s->nSatGroups = 0;
    return( s );
}



void InitSatSelectorInfo() {

    /*
     *  Figure out what Satellite Groups we have. These are determined by the
     *  files that exist in the SatGroups directory
     */
    SatSelectorInfo = New_SatSelectorInfo();
    LoadSpaceObjects();

}




GtkWidget *CreateSatSelector( int nGroups, char *GroupNames[] ) {


    int         i, Row, Col, TotalCols;

//    GtkWidget *window1;
    GtkWidget *vbox1;
    GtkWidget *hbox1;
    GtkWidget *table1;
    GtkWidget *label;
    GtkWidget *checkbutton;
    GtkWidget *button;
    GtkWidget *vseparator;
    GtkWidget *hseparator;
    GtkWidget *colorbutton;
    GtkObject *spinbutton_adj;
    GtkWidget *spinbutton;
    GdkColor   color;




    vbox1 = gtk_vbox_new (FALSE, 0); gtk_widget_show (vbox1);





    table1 = gtk_table_new (7, 26, FALSE); gtk_widget_show (table1);
    gtk_box_pack_start (GTK_BOX (vbox1), table1, TRUE, TRUE, 0);
    gtk_container_set_border_width (GTK_CONTAINER (table1), 10);
    gtk_table_set_row_spacings (GTK_TABLE (table1), 10);
    gtk_table_set_col_spacings (GTK_TABLE (table1), 5);



    /*
     * Set up Column Headers
     */
    Col = 0;
    label = gtk_label_new (_("Show Group")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, Col, Col+1, 0, 1, (GtkAttachOptions) (0), (GtkAttachOptions) (GTK_FILL), 0, 0);
    gtk_misc_set_alignment (GTK_MISC (label), 0.5, 1);
    gtk_label_set_angle (GTK_LABEL (label), 90);
    ++Col;

    label = gtk_label_new (_("Group\nName")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, Col, Col+1, 0, 1, (GtkAttachOptions) (0), (GtkAttachOptions) (GTK_FILL), 0, 0);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_CENTER);
    gtk_misc_set_alignment (GTK_MISC (label), 0, 1);
    ++Col;

    label = gtk_label_new (_("Show Satellites")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, Col, Col+1, 0, 1, (GtkAttachOptions) (0), (GtkAttachOptions) (GTK_FILL), 0, 0);
    gtk_misc_set_alignment (GTK_MISC (label), 0.5, 1);
    gtk_label_set_angle (GTK_LABEL (label), 90);
    ++Col;

    label = gtk_label_new (_("Show Rocket Bodies")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, Col, Col+1, 0, 1, (GtkAttachOptions) (0), (GtkAttachOptions) (GTK_FILL), 0, 0);
    gtk_misc_set_alignment (GTK_MISC (label), 0.5, 1);
    gtk_label_set_angle (GTK_LABEL (label), 90);
    ++Col;

    label = gtk_label_new (_("Show Debris")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, Col, Col+1, 0, 1, (GtkAttachOptions) (0), (GtkAttachOptions) (GTK_FILL), 0, 0);
    gtk_misc_set_alignment (GTK_MISC (label), 0.5, 1);
    gtk_label_set_angle (GTK_LABEL (label), 90);
    ++Col;

    vseparator = gtk_vseparator_new ();
    gtk_widget_show (vseparator);
    gtk_table_attach (GTK_TABLE (table1), vseparator, Col, Col+1, 0, 1, (GtkAttachOptions) (0), (GtkAttachOptions) (GTK_FILL), 0, 0);
    gtk_widget_set_size_request (vseparator, 2, -1);
    ++Col;

    label = gtk_label_new (_("Label")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, Col, Col+1, 0, 1, (GtkAttachOptions) (0), (GtkAttachOptions) (GTK_FILL), 0, 0);
    gtk_misc_set_alignment (GTK_MISC (label), 0.5, 1);
    gtk_label_set_angle (GTK_LABEL (label), 90);
    ++Col;

    label = gtk_label_new (_("Current Position")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, Col, Col+1, 0, 1, (GtkAttachOptions) (0), (GtkAttachOptions) (GTK_FILL), 0, 0);
    gtk_misc_set_alignment (GTK_MISC (label), 0.5, 1);
    gtk_label_set_angle (GTK_LABEL (label), 90);
    ++Col;

    label = gtk_label_new (_("Sat to Ground Lines")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, Col, Col+1, 0, 1, (GtkAttachOptions) (0), (GtkAttachOptions) (GTK_FILL), 0, 0);
    gtk_misc_set_alignment (GTK_MISC (label), 0.5, 1);
    gtk_label_set_angle (GTK_LABEL (label), 90);
    ++Col;

label = gtk_label_new (_("Sat Field Lines")); gtk_widget_show (label);
gtk_table_attach (GTK_TABLE (table1), label, Col, Col+1, 0, 1, (GtkAttachOptions) (0), (GtkAttachOptions) (GTK_FILL), 0, 0);
gtk_misc_set_alignment (GTK_MISC (label), 0.5, 1);
gtk_label_set_angle (GTK_LABEL (label), 90);
++Col;

label = gtk_label_new (_("Sat FL Footpoints")); gtk_widget_show (label);
gtk_table_attach (GTK_TABLE (table1), label, Col, Col+1, 0, 1, (GtkAttachOptions) (0), (GtkAttachOptions) (GTK_FILL), 0, 0);
gtk_misc_set_alignment (GTK_MISC (label), 0.5, 1);
gtk_label_set_angle (GTK_LABEL (label), 90);
++Col;



    vseparator = gtk_vseparator_new (); gtk_widget_show (vseparator);
    gtk_table_attach (GTK_TABLE (table1), vseparator, Col, Col+1, 0, 1, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (GTK_FILL), 0, 0);
    gtk_widget_set_size_request (vseparator, 2, -1);
    ++Col;

    label = gtk_label_new (_("Orbit")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, Col, Col+1, 0, 1, (GtkAttachOptions) (0), (GtkAttachOptions) (GTK_FILL), 0, 0);
    gtk_misc_set_alignment (GTK_MISC (label), 0.5, 1);
    gtk_label_set_angle (GTK_LABEL (label), 90);
    ++Col;

    label = gtk_label_new (_("Ground Path of Orbit")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, Col, Col+1, 0, 1, (GtkAttachOptions) (0), (GtkAttachOptions) (GTK_FILL), 0, 0);
    gtk_misc_set_alignment (GTK_MISC (label), 0.5, 1);
    gtk_label_set_angle (GTK_LABEL (label), 90);
    ++Col;

    label = gtk_label_new (_("Orbit to Ground Lines")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, Col, Col+1, 0, 1, (GtkAttachOptions) (0), (GtkAttachOptions) (GTK_FILL), 0, 0);
    gtk_misc_set_alignment (GTK_MISC (label), 0.5, 1);
    gtk_label_set_angle (GTK_LABEL (label), 90);
    ++Col;

label = gtk_label_new (_("Orbit Field Lines")); gtk_widget_show (label);
gtk_table_attach (GTK_TABLE (table1), label, Col, Col+1, 0, 1, (GtkAttachOptions) (0), (GtkAttachOptions) (GTK_FILL), 0, 0);
gtk_misc_set_alignment (GTK_MISC (label), 0.5, 1);
gtk_label_set_angle (GTK_LABEL (label), 90);
++Col;

label = gtk_label_new (_("Orbit FL Footpoints")); gtk_widget_show (label);
gtk_table_attach (GTK_TABLE (table1), label, Col, Col+1, 0, 1, (GtkAttachOptions) (0), (GtkAttachOptions) (GTK_FILL), 0, 0);
gtk_misc_set_alignment (GTK_MISC (label), 0.5, 1);
gtk_label_set_angle (GTK_LABEL (label), 90);
++Col;


    label = gtk_label_new (_("% of Orbit to Draw")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, Col, Col+1, 0, 1, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (GTK_FILL), 0, 0);
    gtk_misc_set_alignment (GTK_MISC (label), 0.5, 1);
    gtk_label_set_angle (GTK_LABEL (label), 90);
    ++Col;

    vseparator = gtk_vseparator_new (); gtk_widget_show (vseparator);
    gtk_table_attach (GTK_TABLE (table1), vseparator, Col, Col+1, 0, 1, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (GTK_FILL), 0, 0);
    gtk_widget_set_size_request (vseparator, 2, -1);
    ++Col;

    label = gtk_label_new (_("Streak")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, Col, Col+1, 0, 1, (GtkAttachOptions) (0), (GtkAttachOptions) (GTK_FILL), 0, 0);
    gtk_label_set_use_markup (GTK_LABEL (label), TRUE);
    gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_CENTER);
    gtk_misc_set_alignment (GTK_MISC (label), 0.5, 1);
    gtk_label_set_angle (GTK_LABEL (label), 90);
    ++Col;

    label = gtk_label_new (_("Ground Path of Streak")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, Col, Col+1, 0, 1, (GtkAttachOptions) (0), (GtkAttachOptions) (GTK_FILL), 0, 0);
    gtk_misc_set_alignment (GTK_MISC (label), 0.5, 1);
    gtk_label_set_angle (GTK_LABEL (label), 90);
    ++Col;

    label = gtk_label_new (_("Streak to Ground Lines")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, Col, Col+1, 0, 1, (GtkAttachOptions) (0), (GtkAttachOptions) (GTK_FILL), 0, 0);
    gtk_misc_set_alignment (GTK_MISC (label), 0.5, 1);
    gtk_label_set_angle (GTK_LABEL (label), 90);
    ++Col;

label = gtk_label_new (_("Streak Field Lines")); gtk_widget_show (label);
gtk_table_attach (GTK_TABLE (table1), label, Col, Col+1, 0, 1, (GtkAttachOptions) (0), (GtkAttachOptions) (GTK_FILL), 0, 0);
gtk_misc_set_alignment (GTK_MISC (label), 0.5, 1);
gtk_label_set_angle (GTK_LABEL (label), 90);
++Col;

label = gtk_label_new (_("Streak FL Footpoints")); gtk_widget_show (label);
gtk_table_attach (GTK_TABLE (table1), label, Col, Col+1, 0, 1, (GtkAttachOptions) (0), (GtkAttachOptions) (GTK_FILL), 0, 0);
gtk_misc_set_alignment (GTK_MISC (label), 0.5, 1);
gtk_label_set_angle (GTK_LABEL (label), 90);
++Col;


    label = gtk_label_new (_("% of Orbit to Draw")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, Col, Col+1, 0, 1, (GtkAttachOptions) (0), (GtkAttachOptions) (GTK_FILL), 0, 0);
    gtk_misc_set_alignment (GTK_MISC (label), 0.5, 1);
    gtk_label_set_angle (GTK_LABEL (label), 90);
    ++Col;

    vseparator = gtk_vseparator_new (); gtk_widget_show (vseparator);
    gtk_table_attach (GTK_TABLE (table1), vseparator, Col, Col+1, 0, 1, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (GTK_FILL), 0, 0);
    gtk_widget_set_size_request (vseparator, 2, -1);
    ++Col;

    label = gtk_label_new (_("Draw as Point Sprites")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, Col, Col+1, 0, 1, (GtkAttachOptions) (0), (GtkAttachOptions) (GTK_FILL), 0, 0);
    gtk_misc_set_alignment (GTK_MISC (label), 0.5, 1);
    gtk_label_set_angle (GTK_LABEL (label), 90);
    ++Col;

    label = gtk_label_new (_("Sat Color")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, Col, Col+1, 0, 1, (GtkAttachOptions) (0), (GtkAttachOptions) (GTK_FILL), 0, 0);
    gtk_misc_set_alignment (GTK_MISC (label), 0.5, 1);
    gtk_label_set_angle (GTK_LABEL (label), 90);
    ++Col;

    label = gtk_label_new (_("Rocket Body Color")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, Col, Col+1, 0, 1, (GtkAttachOptions) (0), (GtkAttachOptions) (GTK_FILL), 0, 0);
    gtk_misc_set_alignment (GTK_MISC (label), 0.5, 1);
    gtk_label_set_angle (GTK_LABEL (label), 90);
    ++Col;

    label = gtk_label_new (_("Debris Color\n")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, Col, Col+1, 0, 1, (GtkAttachOptions) (0), (GtkAttachOptions) (GTK_FILL), 0, 0);
    gtk_misc_set_alignment (GTK_MISC (label), 0.5, 1);
    gtk_label_set_angle (GTK_LABEL (label), 90);
    ++Col;

    TotalCols = Col;


    /*
     * Horizontal seperator
     */
    hseparator = gtk_hseparator_new (); gtk_widget_show (hseparator);
    gtk_table_attach (GTK_TABLE (table1), hseparator, 0, TotalCols, 1, 2, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);
    gtk_widget_set_size_request (hseparator, -1, 4);



    /*
     * Do Rows of widgets
     */
    _GroupNode   *g;
    _SpaceObjects *Group;
    g   = SatSelectorInfo->SatGroupList;
    Row = 2;
    i   = 0;
    while ( g != NULL ) {

        Group = g->Group;

        Col = 0;
        checkbutton = gtk_check_button_new_with_mnemonic (""); gtk_widget_show (checkbutton);
        gtk_table_attach (GTK_TABLE (table1), checkbutton, Col, Col+1, Row, Row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
        g_signal_connect( G_OBJECT( checkbutton ), "toggled", G_CALLBACK( ToggleOrbitOptions ), GINT_TO_POINTER( i*100+0 ) );
        ++Col;

        // for Group button
        button = gtk_button_new_with_mnemonic( g->GroupName ); gtk_widget_show (button);
        gtk_table_attach (GTK_TABLE (table1), button, Col, Col+1, Row, Row+1, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);
        ++Col;

        // for Show Sats
        checkbutton = gtk_check_button_new_with_mnemonic (""); gtk_widget_show (checkbutton);
        gtk_table_attach (GTK_TABLE (table1), checkbutton, Col, Col+1, Row, Row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
        gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(checkbutton), Group->DrawSatellites );
        g_signal_connect( G_OBJECT( checkbutton ), "toggled", G_CALLBACK( ToggleOrbitOptions ), GINT_TO_POINTER( i*100+3 ) );
        ++Col;

        // for Show R/B
        checkbutton = gtk_check_button_new_with_mnemonic (""); gtk_widget_show (checkbutton);
        gtk_table_attach (GTK_TABLE (table1), checkbutton, Col, Col+1, Row, Row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
        gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(checkbutton), Group->DrawRocketBodies );
        g_signal_connect( G_OBJECT( checkbutton ), "toggled", G_CALLBACK( ToggleOrbitOptions ), GINT_TO_POINTER( i*100+4 ) );
        ++Col;

        // for Show Deb
        checkbutton = gtk_check_button_new_with_mnemonic (""); gtk_widget_show (checkbutton);
        gtk_table_attach (GTK_TABLE (table1), checkbutton, Col, Col+1, Row, Row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
        gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(checkbutton), Group->DrawDebris );
        g_signal_connect( G_OBJECT( checkbutton ), "toggled", G_CALLBACK( ToggleOrbitOptions ), GINT_TO_POINTER( i*100+5 ) );
        ++Col;

        vseparator = gtk_vseparator_new (); gtk_widget_show (vseparator);
        gtk_table_attach (GTK_TABLE (table1), vseparator, Col, Col+1, Row, Row+1, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (GTK_FILL), 0, 0);
        gtk_widget_set_size_request (vseparator, 2, -1);
        ++Col;


        // for Label
        checkbutton = gtk_check_button_new_with_mnemonic (""); gtk_widget_show (checkbutton);
        gtk_table_attach (GTK_TABLE (table1), checkbutton, Col, Col+1, Row, Row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
        gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(checkbutton), Group->Sat[0].DrawLabel );
        g_signal_connect( G_OBJECT( checkbutton ), "toggled", G_CALLBACK( ToggleOrbitOptions ), GINT_TO_POINTER( i*100+1 ) );
        ++Col;

        // for Current Position
        checkbutton = gtk_check_button_new_with_mnemonic (""); gtk_widget_show (checkbutton);
        gtk_table_attach (GTK_TABLE (table1), checkbutton, Col, Col+1, Row, Row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
        gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(checkbutton), Group->DrawPosition );
        g_signal_connect( G_OBJECT( checkbutton ), "toggled", G_CALLBACK( ToggleOrbitOptions ), GINT_TO_POINTER( i*100+13 ) );
        ++Col;

        // for Sat to Ground Lines
        checkbutton = gtk_check_button_new_with_mnemonic (""); gtk_widget_show (checkbutton);
        gtk_table_attach (GTK_TABLE (table1), checkbutton, Col, Col+1, Row, Row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
        g_signal_connect( G_OBJECT( checkbutton ), "toggled", G_CALLBACK( ToggleOrbitOptions ), GINT_TO_POINTER( i*100+2 ) );
        ++Col;

// for Sat Field Lines
checkbutton = gtk_check_button_new_with_mnemonic (""); gtk_widget_show (checkbutton);
gtk_table_attach (GTK_TABLE (table1), checkbutton, Col, Col+1, Row, Row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(checkbutton), Group->Sat[0].DrawSatFieldLines );
g_signal_connect( G_OBJECT( checkbutton ), "toggled", G_CALLBACK( ToggleOrbitOptions ), GINT_TO_POINTER( i*100+14 ) );
++Col;


// for Sat FL Footpoints
checkbutton = gtk_check_button_new_with_mnemonic (""); gtk_widget_show (checkbutton);
gtk_table_attach (GTK_TABLE (table1), checkbutton, Col, Col+1, Row, Row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(checkbutton), Group->Sat[0].DrawSatFLFootpoints );
g_signal_connect( G_OBJECT( checkbutton ), "toggled", G_CALLBACK( ToggleOrbitOptions ), GINT_TO_POINTER( i*100+15 ) );
++Col;

        vseparator = gtk_vseparator_new (); gtk_widget_show (vseparator);
        gtk_table_attach (GTK_TABLE (table1), vseparator, Col, Col+1, Row, Row+1, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (GTK_FILL), 0, 0);
        gtk_widget_set_size_request (vseparator, 2, -1);
        ++Col;


        // for Orbit
        checkbutton = gtk_check_button_new_with_mnemonic (""); gtk_widget_show (checkbutton);
        gtk_table_attach (GTK_TABLE (table1), checkbutton, Col, Col+1, Row, Row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
        gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(checkbutton), Group->Sat[0].DrawOrbit );
        g_signal_connect( G_OBJECT( checkbutton ), "toggled", G_CALLBACK( ToggleOrbitOptions ), GINT_TO_POINTER( i*100+6 ) );
        ++Col;

        // for Ground Path of Orbit
        checkbutton = gtk_check_button_new_with_mnemonic (""); gtk_widget_show (checkbutton);
        gtk_table_attach (GTK_TABLE (table1), checkbutton, Col, Col+1, Row, Row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
        gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(checkbutton), Group->Sat[0].DrawGroundPathOfOrbit );
        g_signal_connect( G_OBJECT( checkbutton ), "toggled", G_CALLBACK( ToggleOrbitOptions ), GINT_TO_POINTER( i*100+7 ) );
        ++Col;

        // for Orbit to Ground Lines
        checkbutton = gtk_check_button_new_with_mnemonic (""); gtk_widget_show (checkbutton);
        gtk_table_attach (GTK_TABLE (table1), checkbutton, Col, Col+1, Row, Row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
        gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(checkbutton), Group->Sat[0].DrawOrbitToGroundLines );
        g_signal_connect( G_OBJECT( checkbutton ), "toggled", G_CALLBACK( ToggleOrbitOptions ), GINT_TO_POINTER( i*100+8 ) );
        ++Col;

// for Orbit Field Lines
checkbutton = gtk_check_button_new_with_mnemonic (""); gtk_widget_show (checkbutton);
gtk_table_attach (GTK_TABLE (table1), checkbutton, Col, Col+1, Row, Row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(checkbutton), Group->Sat[0].DrawOrbitFieldLines );
g_signal_connect( G_OBJECT( checkbutton ), "toggled", G_CALLBACK( ToggleOrbitOptions ), GINT_TO_POINTER( i*100+16 ) );
++Col;

// for Orbit FL Footpoints
checkbutton = gtk_check_button_new_with_mnemonic (""); gtk_widget_show (checkbutton);
gtk_table_attach (GTK_TABLE (table1), checkbutton, Col, Col+1, Row, Row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(checkbutton), Group->Sat[0].DrawOrbitFLFootpoints );
g_signal_connect( G_OBJECT( checkbutton ), "toggled", G_CALLBACK( ToggleOrbitOptions ), GINT_TO_POINTER( i*100+17 ) );
++Col;

        // for % Orbit to Draw
        spinbutton_adj = gtk_adjustment_new (100, 1, 2000, 1, 10, 0);
        spinbutton = gtk_spin_button_new (GTK_ADJUSTMENT (spinbutton_adj), 1, 0); gtk_widget_show (spinbutton);
        gtk_table_attach (GTK_TABLE (table1), spinbutton, Col, Col+1, Row, Row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
        gtk_spin_button_set_numeric (GTK_SPIN_BUTTON (spinbutton), TRUE);
        gtk_spin_button_set_value( GTK_SPIN_BUTTON (spinbutton), Group->Sat[0].oPeriodFrac );
        g_signal_connect( G_OBJECT( spinbutton ), "value-changed", G_CALLBACK( ChangeOrbitOptions ), GINT_TO_POINTER( i*100+0 ) );
        ++Col;

        vseparator = gtk_vseparator_new (); gtk_widget_show (vseparator);
        gtk_table_attach (GTK_TABLE (table1), vseparator, Col, Col+1, Row, Row+1, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (GTK_FILL), 0, 0);
        gtk_widget_set_size_request (vseparator, 2, -1);
        ++Col;


        // for Streak
        checkbutton = gtk_check_button_new_with_mnemonic (""); gtk_widget_show (checkbutton);
        gtk_table_attach (GTK_TABLE (table1), checkbutton, Col, Col+1, Row, Row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
        gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(checkbutton), Group->Sat[0].DrawStreak );
        g_signal_connect( G_OBJECT( checkbutton ), "toggled", G_CALLBACK( ToggleOrbitOptions ), GINT_TO_POINTER( i*100+9 ) );
        ++Col;

        // for Ground Path of Streak
        checkbutton = gtk_check_button_new_with_mnemonic (""); gtk_widget_show (checkbutton);
        gtk_table_attach (GTK_TABLE (table1), checkbutton, Col, Col+1, Row, Row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
        gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(checkbutton), Group->Sat[0].DrawGroundPathOfStreak );
        g_signal_connect( G_OBJECT( checkbutton ), "toggled", G_CALLBACK( ToggleOrbitOptions ), GINT_TO_POINTER( i*100+10 ) );
        ++Col;

        // for Streak to Ground Lines
        checkbutton = gtk_check_button_new_with_mnemonic (""); gtk_widget_show (checkbutton);
        gtk_table_attach (GTK_TABLE (table1), checkbutton, Col, Col+1, Row, Row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
        gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(checkbutton), Group->Sat[0].DrawStreakToGroundLines );
        g_signal_connect( G_OBJECT( checkbutton ), "toggled", G_CALLBACK( ToggleOrbitOptions ), GINT_TO_POINTER( i*100+11 ) );
        ++Col;

// for Streak Field Lines
checkbutton = gtk_check_button_new_with_mnemonic (""); gtk_widget_show (checkbutton);
gtk_table_attach (GTK_TABLE (table1), checkbutton, Col, Col+1, Row, Row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(checkbutton), Group->Sat[0].DrawStreakFieldLines );
g_signal_connect( G_OBJECT( checkbutton ), "toggled", G_CALLBACK( ToggleOrbitOptions ), GINT_TO_POINTER( i*100+18 ) );
++Col;

// for Streak FL Footpoints
checkbutton = gtk_check_button_new_with_mnemonic (""); gtk_widget_show (checkbutton);
gtk_table_attach (GTK_TABLE (table1), checkbutton, Col, Col+1, Row, Row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(checkbutton), Group->Sat[0].DrawStreakFLFootpoints );
g_signal_connect( G_OBJECT( checkbutton ), "toggled", G_CALLBACK( ToggleOrbitOptions ), GINT_TO_POINTER( i*100+19 ) );
++Col;


        // for % Orbit to Draw
        spinbutton_adj = gtk_adjustment_new (25, 1, 2000, 1, 10, 0);
        spinbutton = gtk_spin_button_new (GTK_ADJUSTMENT (spinbutton_adj), 1, 0); gtk_widget_show (spinbutton);
        gtk_table_attach (GTK_TABLE (table1), spinbutton, Col, Col+1, Row, Row+1, (GtkAttachOptions) (GTK_EXPAND), (GtkAttachOptions) (0), 0, 0);
        gtk_spin_button_set_numeric (GTK_SPIN_BUTTON (spinbutton), TRUE);
        gtk_spin_button_set_value( GTK_SPIN_BUTTON (spinbutton), Group->Sat[0].sPeriodFrac );
        g_signal_connect( G_OBJECT( spinbutton ), "value-changed", G_CALLBACK( ChangeOrbitOptions ), GINT_TO_POINTER( i*100+1 ) );
        ++Col;

        vseparator = gtk_vseparator_new (); gtk_widget_show (vseparator);
        gtk_table_attach (GTK_TABLE (table1), vseparator, Col, Col+1, Row, Row+1, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (GTK_FILL), 0, 0);
        gtk_widget_set_size_request (vseparator, 2, -1);
        ++Col;

        // for Draw as Point Sprites
        checkbutton = gtk_check_button_new_with_mnemonic (""); gtk_widget_show (checkbutton);
        gtk_table_attach (GTK_TABLE (table1), checkbutton, Col, Col+1, Row, Row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
        gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(checkbutton), Group->DrawAsPointSprites );
        g_signal_connect( G_OBJECT( checkbutton ), "toggled", G_CALLBACK( ToggleOrbitOptions ), GINT_TO_POINTER( i*100+12 ) );
        ++Col;

        // for Sat Color
        color.red = Group->SatRed*65535; color.green = Group->SatGrn*65535; color.blue = Group->SatBlu*65535; 
        colorbutton = gtk_color_button_new (); gtk_widget_show (colorbutton);
        gtk_color_button_set_color( GTK_COLOR_BUTTON(colorbutton), &color );
        gtk_color_button_set_use_alpha( GTK_COLOR_BUTTON(colorbutton), TRUE );
        gtk_color_button_set_alpha( GTK_COLOR_BUTTON(colorbutton), Group->SatAlf*65535 );
        gtk_table_attach (GTK_TABLE (table1), colorbutton, Col, Col+1, Row, Row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
        gtk_widget_set_size_request (colorbutton, 30, 30);
        gtk_color_button_set_use_alpha (GTK_COLOR_BUTTON (colorbutton), TRUE);
        g_signal_connect( G_OBJECT( colorbutton ), "color-set", G_CALLBACK( ChangeSatColor ), GINT_TO_POINTER(i*100+0) );
        ++Col;


        // for R/B Color
        color.red = Group->RbRed*65535; color.green = Group->RbGrn*65535; color.blue = Group->RbBlu*65535; 
        colorbutton = gtk_color_button_new (); gtk_widget_show (colorbutton);
        gtk_color_button_set_color( GTK_COLOR_BUTTON(colorbutton), &color );
        gtk_color_button_set_use_alpha( GTK_COLOR_BUTTON(colorbutton), TRUE );
        gtk_color_button_set_alpha( GTK_COLOR_BUTTON(colorbutton), Group->RbAlf*65535 );
        gtk_table_attach (GTK_TABLE (table1), colorbutton, Col, Col+1, Row, Row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
        gtk_widget_set_size_request (colorbutton, 30, 30);
        gtk_color_button_set_use_alpha (GTK_COLOR_BUTTON (colorbutton), TRUE);
        g_signal_connect( G_OBJECT( colorbutton ), "color-set", G_CALLBACK( ChangeSatColor ), GINT_TO_POINTER(i*100+1) );
        ++Col;

        // for Deb Color
        color.red = Group->DebRed*65535; color.green = Group->DebGrn*65535; color.blue = Group->DebBlu*65535; 
        colorbutton = gtk_color_button_new (); gtk_widget_show (colorbutton);
        gtk_color_button_set_color( GTK_COLOR_BUTTON(colorbutton), &color );
        gtk_color_button_set_use_alpha( GTK_COLOR_BUTTON(colorbutton), TRUE );
        gtk_color_button_set_alpha( GTK_COLOR_BUTTON(colorbutton), Group->DebAlf*65535 );
        gtk_table_attach (GTK_TABLE (table1), colorbutton, Col, Col+1, Row, Row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
        gtk_widget_set_size_request (colorbutton, 30, 30);
        gtk_color_button_set_use_alpha (GTK_COLOR_BUTTON (colorbutton), TRUE);
        g_signal_connect( G_OBJECT( colorbutton ), "color-set", G_CALLBACK( ChangeSatColor ), GINT_TO_POINTER(i*100+2) );
        ++Col;

        // for More Options button
//        button = gtk_button_new_with_mnemonic (_("More Options")); gtk_widget_show (button);
//        gtk_table_attach (GTK_TABLE (table1), button, 24, 25, Row, Row+1, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);
//        ++Col;

        g = g->Next;
        ++i; 
        ++Row;

    }


    /*
     * Horizontal seperator
     */
    hseparator = gtk_hseparator_new (); gtk_widget_show (hseparator);
    gtk_table_attach (GTK_TABLE (table1), hseparator, 0, 25, Row, Row+1, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);
    gtk_widget_set_size_request (hseparator, -1, 4);


    /*
     * for set of buttons at bottom
     */
    hbox1 = gtk_hbox_new (FALSE, 0); gtk_widget_show (hbox1);
    gtk_box_pack_start (GTK_BOX (vbox1), hbox1, TRUE, TRUE, 0);

    button = gtk_button_new_with_mnemonic (_("Delete Group")); gtk_widget_show (button);
    gtk_box_pack_end (GTK_BOX (hbox1), button, FALSE, FALSE, 0);

    button = gtk_button_new_with_mnemonic (_("Edit Group")); gtk_widget_show (button);
    gtk_box_pack_end (GTK_BOX (hbox1), button, FALSE, FALSE, 0);

    button = gtk_button_new_with_mnemonic (_("Create Group")); gtk_widget_show (button);
    gtk_box_pack_end (GTK_BOX (hbox1), button, FALSE, FALSE, 0);


    return vbox1;




}

