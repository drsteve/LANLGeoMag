#include "SatSelector.h"

extern float IllumFL_ka;
extern float IllumFL_kd;
extern float IllumFL_ks;
extern double IllumFL_n;
extern double IllumFL_w;


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

            Group->oflRed_OverRide = Group->Sat[0].oflRed;
            Group->oflGrn_OverRide = Group->Sat[0].oflGrn;
            Group->oflBlu_OverRide = Group->Sat[0].oflBlu;
            Group->oflAlf_OverRide = Group->Sat[0].oflAlf;

            Group->ofpRed_OverRide = Group->Sat[0].ofpRed;
            Group->ofpGrn_OverRide = Group->Sat[0].ofpGrn;
            Group->ofpBlu_OverRide = Group->Sat[0].ofpBlu;
            Group->ofpAlf_OverRide = Group->Sat[0].ofpAlf;

            Group->OverRideSatColors    = FALSE;
            Group->OverRideOrbitColors  = FALSE;
            Group->OverRideStreakColors = FALSE;

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

                Group->Sat[j].sflRed = 0.6;
                Group->Sat[j].sflGrn = 0.6;
                Group->Sat[j].sflBlu = 0.6;
                Group->Sat[j].sflAlf = 0.7;

                Group->Sat[j].sfpRed = 0.6;
                Group->Sat[j].sfpGrn = 0.6;
                Group->Sat[j].sfpBlu = 0.6;
                Group->Sat[j].sfpAlf = 0.7;

                // Orbit Colors
                Group->Sat[j].oRed = 0.6;
                Group->Sat[j].oGrn = 0.6;
                Group->Sat[j].oBlu = 0.6;
                Group->Sat[j].oAlf = 0.7;

                Group->Sat[j].ogpRed = 0.6;
                Group->Sat[j].ogpGrn = 0.6;
                Group->Sat[j].ogpBlu = 0.6;
                Group->Sat[j].ogpAlf = 0.7;

                Group->Sat[j].oglRed = 0.6;
                Group->Sat[j].oglGrn = 0.6;
                Group->Sat[j].oglBlu = 0.6;
                Group->Sat[j].oglAlf = 0.7;

                Group->Sat[j].oflRed = 0.6;
                Group->Sat[j].oflGrn = 0.6;
                Group->Sat[j].oflBlu = 0.6;
                Group->Sat[j].oflAlf = 0.7;

                Group->Sat[j].ofpRed = 0.6;
                Group->Sat[j].ofpGrn = 0.6;
                Group->Sat[j].ofpBlu = 0.6;
                Group->Sat[j].ofpAlf = 0.7;

                // Single Line
                Group->Sat[j].ssglRed = r;
                Group->Sat[j].ssglGrn = g;
                Group->Sat[j].ssglBlu = b;
                Group->Sat[j].ssglAlf = 0.5*0.75;

            }
            Group->oRed_OverRide = Group->Sat[0].sRed;
            Group->oGrn_OverRide = Group->Sat[0].sGrn;
            Group->oBlu_OverRide = Group->Sat[0].sBlu;
            Group->oAlf_OverRide = Group->Sat[0].sAlf;

            Group->ogpRed_OverRide = Group->Sat[0].sRed;
            Group->ogpGrn_OverRide = Group->Sat[0].sGrn;
            Group->ogpBlu_OverRide = Group->Sat[0].sBlu;
            Group->ogpAlf_OverRide = Group->Sat[0].sAlf;

            Group->oglRed_OverRide = Group->Sat[0].sRed;
            Group->oglGrn_OverRide = Group->Sat[0].sGrn;
            Group->oglBlu_OverRide = Group->Sat[0].sBlu;
            Group->oglAlf_OverRide = Group->Sat[0].sAlf;

            Group->oflRed_OverRide = Group->Sat[0].sRed;
            Group->oflGrn_OverRide = Group->Sat[0].sGrn;
            Group->oflBlu_OverRide = Group->Sat[0].sBlu;
            Group->oflAlf_OverRide = Group->Sat[0].sAlf;

            Group->ofpRed_OverRide = Group->Sat[0].sRed;
            Group->ofpGrn_OverRide = Group->Sat[0].sGrn;
            Group->ofpBlu_OverRide = Group->Sat[0].sBlu;
            Group->ofpAlf_OverRide = Group->Sat[0].sAlf;

            Group->sRed_OverRide = Group->Sat[0].sRed;
            Group->sGrn_OverRide = Group->Sat[0].sGrn;
            Group->sBlu_OverRide = Group->Sat[0].sBlu;
            Group->sAlf_OverRide = Group->Sat[0].sAlf;

            Group->sgpRed_OverRide = Group->Sat[0].sRed;
            Group->sgpGrn_OverRide = Group->Sat[0].sGrn;
            Group->sgpBlu_OverRide = Group->Sat[0].sBlu;
            Group->sgpAlf_OverRide = Group->Sat[0].sAlf;

            Group->sglRed_OverRide = Group->Sat[0].sRed;
            Group->sglGrn_OverRide = Group->Sat[0].sGrn;
            Group->sglBlu_OverRide = Group->Sat[0].sBlu;
            Group->sglAlf_OverRide = Group->Sat[0].sAlf;

            Group->sflRed_OverRide = Group->Sat[0].sRed;
            Group->sflGrn_OverRide = Group->Sat[0].sGrn;
            Group->sflBlu_OverRide = Group->Sat[0].sBlu;
            Group->sflAlf_OverRide = Group->Sat[0].sAlf;

            Group->sfpRed_OverRide = Group->Sat[0].sRed;
            Group->sfpGrn_OverRide = Group->Sat[0].sGrn;
            Group->sfpBlu_OverRide = Group->Sat[0].sBlu;
            Group->sfpAlf_OverRide = Group->Sat[0].sAlf;

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
    //gtk_container_set_border_width (GTK_CONTAINER (table1), 10);
    gtk_container_set_border_width (GTK_CONTAINER (table1), 0);
    gtk_table_set_row_spacings( GTK_TABLE(table1), 0 );
    gtk_table_set_col_spacings( GTK_TABLE(table1), 0 );



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

/*
Never really use these anymore...
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
*/

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

    label = gtk_label_new (_("OverRide Sat Color")); gtk_widget_show (label);
    gtk_table_attach (GTK_TABLE (table1), label, Col, Col+1, 0, 1, (GtkAttachOptions) (0), (GtkAttachOptions) (GTK_FILL), 0, 0);
    gtk_misc_set_alignment (GTK_MISC (label), 0.5, 1);
    gtk_label_set_angle (GTK_LABEL (label), 90);
    ++Col;

    label = gtk_label_new (_("Draw as Point Sprites")); gtk_widget_show (label);
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

    label = gtk_label_new (_("Override Colors")); gtk_widget_show (label);
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

    label = gtk_label_new (_("OverRide Streak Color\n")); gtk_widget_show (label);
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

        /******************************
         *
         *  Overall Group Quantities
         *
         ******************************/
        Col = 0;
        checkbutton = gtk_check_button_new_with_mnemonic (""); gtk_widget_show (checkbutton);
        gtk_table_attach (GTK_TABLE (table1), checkbutton, Col, Col+1, Row, Row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
        g_signal_connect( G_OBJECT( checkbutton ), "toggled", G_CALLBACK( ToggleOrbitOptions ), GINT_TO_POINTER( i*100 + 0 ) );
        ++Col;

        // for Group button
        button = gtk_button_new_with_mnemonic( g->GroupName ); gtk_widget_show (button);
        gtk_table_attach (GTK_TABLE (table1), button, Col, Col+1, Row, Row+1, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);
        ++Col;

        vseparator = gtk_vseparator_new (); gtk_widget_show (vseparator);
        gtk_table_attach (GTK_TABLE (table1), vseparator, Col, Col+1, Row, Row+1, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (GTK_FILL), 0, 0);
        gtk_widget_set_size_request (vseparator, 2, -1);
        ++Col;

        /******************************
         *
         *      Position Quantities
         *
         ******************************/
        // for Label
        checkbutton = gtk_check_button_new_with_mnemonic (""); gtk_widget_show (checkbutton);
        gtk_table_attach (GTK_TABLE (table1), checkbutton, Col, Col+1, Row, Row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
        gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(checkbutton), Group->Sat[0].DrawLabel );
        g_signal_connect( G_OBJECT( checkbutton ), "toggled", G_CALLBACK( ToggleOrbitOptions ), GINT_TO_POINTER( i*100 + 10 ) );
        ++Col;

        // for Current Position
        checkbutton = gtk_check_button_new_with_mnemonic (""); gtk_widget_show (checkbutton);
        gtk_table_attach (GTK_TABLE (table1), checkbutton, Col, Col+1, Row, Row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
        gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(checkbutton), Group->DrawPosition );
        g_signal_connect( G_OBJECT( checkbutton ), "toggled", G_CALLBACK( ToggleOrbitOptions ), GINT_TO_POINTER( i*100 + 11 ) );
        //++Col;
        color.red   = Group->Sat[0].oRed*65535; color.green = Group->Sat[0].oGrn*65535; color.blue  = Group->Sat[0].oBlu*65535; 
        colorbutton = gtk_color_button_new(); gtk_widget_show(colorbutton);
        gtk_color_button_set_color( GTK_COLOR_BUTTON(colorbutton), &color );
        gtk_color_button_set_use_alpha( GTK_COLOR_BUTTON(colorbutton), TRUE );
        gtk_color_button_set_alpha( GTK_COLOR_BUTTON(colorbutton), Group->RbAlf*65535 );
        gtk_table_attach( GTK_TABLE(table1), colorbutton, Col, Col+1, Row+1, Row+2, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0 );
        gtk_widget_set_size_request( colorbutton, 30, 25 );
        gtk_color_button_set_use_alpha( GTK_COLOR_BUTTON(colorbutton), TRUE );
        g_signal_connect( G_OBJECT(colorbutton), "color-set", G_CALLBACK(ChangeSatColor), GINT_TO_POINTER(i*100 + 11) );
        ++Col;


        // for Sat to Ground Lines
        checkbutton = gtk_check_button_new_with_mnemonic (""); gtk_widget_show (checkbutton);
        gtk_table_attach (GTK_TABLE (table1), checkbutton, Col, Col+1, Row, Row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
        g_signal_connect( G_OBJECT( checkbutton ), "toggled", G_CALLBACK( ToggleOrbitOptions ), GINT_TO_POINTER( i*100 + 12 ) );
        //++Col;
        color.red   = Group->Sat[0].oRed*65535; color.green = Group->Sat[0].oGrn*65535; color.blue  = Group->Sat[0].oBlu*65535; 
        colorbutton = gtk_color_button_new(); gtk_widget_show(colorbutton);
        gtk_color_button_set_color( GTK_COLOR_BUTTON(colorbutton), &color );
        gtk_color_button_set_use_alpha( GTK_COLOR_BUTTON(colorbutton), TRUE );
        gtk_color_button_set_alpha( GTK_COLOR_BUTTON(colorbutton), Group->RbAlf*65535 );
        gtk_table_attach( GTK_TABLE(table1), colorbutton, Col, Col+1, Row+1, Row+2, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0 );
        gtk_widget_set_size_request( colorbutton, 30, 25 );
        gtk_color_button_set_use_alpha( GTK_COLOR_BUTTON(colorbutton), TRUE );
        g_signal_connect( G_OBJECT(colorbutton), "color-set", G_CALLBACK(ChangeSatColor), GINT_TO_POINTER(i*100 + 12) );
        ++Col;


        // for Sat Field Lines
        checkbutton = gtk_check_button_new_with_mnemonic (""); gtk_widget_show (checkbutton);
        gtk_table_attach (GTK_TABLE (table1), checkbutton, Col, Col+1, Row, Row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
        gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(checkbutton), Group->Sat[0].DrawSatFieldLines );
        g_signal_connect( G_OBJECT( checkbutton ), "toggled", G_CALLBACK( ToggleOrbitOptions ), GINT_TO_POINTER( i*100 + 13 ) );
        //++Col;
        color.red   = Group->Sat[0].oRed*65535; color.green = Group->Sat[0].oGrn*65535; color.blue  = Group->Sat[0].oBlu*65535; 
        colorbutton = gtk_color_button_new(); gtk_widget_show(colorbutton);
        gtk_color_button_set_color( GTK_COLOR_BUTTON(colorbutton), &color );
        gtk_color_button_set_use_alpha( GTK_COLOR_BUTTON(colorbutton), TRUE );
        gtk_color_button_set_alpha( GTK_COLOR_BUTTON(colorbutton), Group->RbAlf*65535 );
        gtk_table_attach( GTK_TABLE(table1), colorbutton, Col, Col+1, Row+1, Row+2, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0 );
        gtk_widget_set_size_request( colorbutton, 30, 25 );
        gtk_color_button_set_use_alpha( GTK_COLOR_BUTTON(colorbutton), TRUE );
        g_signal_connect( G_OBJECT(colorbutton), "color-set", G_CALLBACK(ChangeSatColor), GINT_TO_POINTER(i*100 + 13) );
        ++Col;



        // for Sat FL Footpoints
        checkbutton = gtk_check_button_new_with_mnemonic (""); gtk_widget_show (checkbutton);
        gtk_table_attach (GTK_TABLE (table1), checkbutton, Col, Col+1, Row, Row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
        gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(checkbutton), Group->Sat[0].DrawSatFLFootpoints );
        g_signal_connect( G_OBJECT( checkbutton ), "toggled", G_CALLBACK( ToggleOrbitOptions ), GINT_TO_POINTER( i*100 + 14 ) );
        //++Col;
        color.red   = Group->Sat[0].oRed*65535; color.green = Group->Sat[0].oGrn*65535; color.blue  = Group->Sat[0].oBlu*65535; 
        colorbutton = gtk_color_button_new(); gtk_widget_show(colorbutton);
        gtk_color_button_set_color( GTK_COLOR_BUTTON(colorbutton), &color );
        gtk_color_button_set_use_alpha( GTK_COLOR_BUTTON(colorbutton), TRUE );
        gtk_color_button_set_alpha( GTK_COLOR_BUTTON(colorbutton), Group->RbAlf*65535 );
        gtk_table_attach( GTK_TABLE(table1), colorbutton, Col, Col+1, Row+1, Row+2, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0 );
        gtk_widget_set_size_request( colorbutton, 30, 25 );
        gtk_color_button_set_use_alpha( GTK_COLOR_BUTTON(colorbutton), TRUE );
        g_signal_connect( G_OBJECT(colorbutton), "color-set", G_CALLBACK(ChangeSatColor), GINT_TO_POINTER(i*100 + 14) );
        ++Col;


        // Override Sat Colors
        checkbutton = gtk_check_button_new_with_mnemonic (""); gtk_widget_show (checkbutton);
        gtk_table_attach (GTK_TABLE (table1), checkbutton, Col, Col+1, Row+1, Row+2, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
        gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(checkbutton),  Group->OverRideSatColors );
        g_signal_connect( G_OBJECT( checkbutton ), "toggled", G_CALLBACK( ToggleOrbitOptions ), GINT_TO_POINTER( i*100 + 16 ) );
        ++Col;


        // for Draw as Point Sprites
        checkbutton = gtk_check_button_new_with_mnemonic (""); gtk_widget_show (checkbutton);
        gtk_table_attach (GTK_TABLE (table1), checkbutton, Col, Col+1, Row, Row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
        gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(checkbutton), Group->DrawAsPointSprites );
        g_signal_connect( G_OBJECT( checkbutton ), "toggled", G_CALLBACK( ToggleOrbitOptions ), GINT_TO_POINTER( i*100 + 15 ) );
        ++Col;




        vseparator = gtk_vseparator_new (); gtk_widget_show (vseparator);
        gtk_table_attach (GTK_TABLE (table1), vseparator, Col, Col+1, Row, Row+1, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (GTK_FILL), 0, 0);
        gtk_widget_set_size_request (vseparator, 2, -1);
        ++Col;



        /******************************
         *
         *      Orbit Quantities
         *
         ******************************/
        // Draw orbit
        checkbutton = gtk_check_button_new_with_mnemonic (""); gtk_widget_show (checkbutton);
        gtk_table_attach (GTK_TABLE (table1), checkbutton, Col, Col+1, Row, Row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
        gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(checkbutton), Group->Sat[0].DrawOrbit );
        g_signal_connect( G_OBJECT( checkbutton ), "toggled", G_CALLBACK( ToggleOrbitOptions ), GINT_TO_POINTER( i*100 + 20 ) );
//        ++Col;
        color.red   = Group->oRed_OverRide*65535; color.green = Group->oGrn_OverRide*65535; color.blue  = Group->oBlu_OverRide*65535; 
        colorbutton = gtk_color_button_new(); gtk_widget_show(colorbutton);
        gtk_color_button_set_color( GTK_COLOR_BUTTON(colorbutton), &color );
        gtk_color_button_set_use_alpha( GTK_COLOR_BUTTON(colorbutton), TRUE );
        gtk_color_button_set_alpha( GTK_COLOR_BUTTON(colorbutton), Group->oAlf_OverRide*65535 );
        gtk_table_attach( GTK_TABLE(table1), colorbutton, Col, Col+1, Row+1, Row+2, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0 );
        gtk_widget_set_size_request( colorbutton, 30, 25 );
        gtk_color_button_set_use_alpha( GTK_COLOR_BUTTON(colorbutton), TRUE );
        g_signal_connect( G_OBJECT(colorbutton), "color-set", G_CALLBACK(ChangeSatColor), GINT_TO_POINTER(i*100 + 20) );
        ++Col;


        // for Ground Path of Orbit
        checkbutton = gtk_check_button_new_with_mnemonic (""); gtk_widget_show (checkbutton);
        gtk_table_attach (GTK_TABLE (table1), checkbutton, Col, Col+1, Row, Row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
        gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(checkbutton), Group->Sat[0].DrawGroundPathOfOrbit );
        g_signal_connect( G_OBJECT( checkbutton ), "toggled", G_CALLBACK( ToggleOrbitOptions ), GINT_TO_POINTER( i*100 + 21 ) );
//        ++Col;
        color.red   = Group->ogpRed_OverRide*65535; color.green = Group->ogpGrn_OverRide*65535; color.blue  = Group->ogpBlu_OverRide*65535; 
        colorbutton = gtk_color_button_new(); gtk_widget_show(colorbutton);
        gtk_color_button_set_color( GTK_COLOR_BUTTON(colorbutton), &color );
        gtk_color_button_set_use_alpha( GTK_COLOR_BUTTON(colorbutton), TRUE );
        gtk_color_button_set_alpha( GTK_COLOR_BUTTON(colorbutton), Group->ogpAlf_OverRide*65535 );
        gtk_table_attach( GTK_TABLE(table1), colorbutton, Col, Col+1, Row+1, Row+2, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0 );
        gtk_widget_set_size_request( colorbutton, 30, 25 );
        gtk_color_button_set_use_alpha( GTK_COLOR_BUTTON(colorbutton), TRUE );
        g_signal_connect( G_OBJECT(colorbutton), "color-set", G_CALLBACK(ChangeSatColor), GINT_TO_POINTER(i*100 + 21) );
        ++Col;



        // for Orbit to Ground Lines
        checkbutton = gtk_check_button_new_with_mnemonic (""); gtk_widget_show (checkbutton);
        gtk_table_attach (GTK_TABLE (table1), checkbutton, Col, Col+1, Row, Row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
        gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(checkbutton), Group->Sat[0].DrawOrbitToGroundLines );
        g_signal_connect( G_OBJECT( checkbutton ), "toggled", G_CALLBACK( ToggleOrbitOptions ), GINT_TO_POINTER( i*100 + 22 ) );
//        ++Col;
        color.red   = Group->oglRed_OverRide*65535; color.green = Group->oglGrn_OverRide*65535; color.blue  = Group->oglBlu_OverRide*65535; 
        colorbutton = gtk_color_button_new(); gtk_widget_show(colorbutton);
        gtk_color_button_set_color( GTK_COLOR_BUTTON(colorbutton), &color );
        gtk_color_button_set_use_alpha( GTK_COLOR_BUTTON(colorbutton), TRUE );
        gtk_color_button_set_alpha( GTK_COLOR_BUTTON(colorbutton), Group->oglAlf_OverRide*65535 );
        gtk_table_attach( GTK_TABLE(table1), colorbutton, Col, Col+1, Row+1, Row+2, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0 );
        gtk_widget_set_size_request( colorbutton, 30, 25 );
        gtk_color_button_set_use_alpha( GTK_COLOR_BUTTON(colorbutton), TRUE );
        g_signal_connect( G_OBJECT(colorbutton), "color-set", G_CALLBACK(ChangeSatColor), GINT_TO_POINTER(i*100 + 22) );
        ++Col;



        // for Orbit Field Lines
        checkbutton = gtk_check_button_new_with_mnemonic (""); gtk_widget_show (checkbutton);
        gtk_table_attach (GTK_TABLE (table1), checkbutton, Col, Col+1, Row, Row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
        gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(checkbutton), Group->Sat[0].DrawOrbitFieldLines );
        g_signal_connect( G_OBJECT( checkbutton ), "toggled", G_CALLBACK( ToggleOrbitOptions ), GINT_TO_POINTER( i*100 + 23 ) );
        //++Col;
        color.red   = Group->oflRed_OverRide*65535; color.green = Group->oflGrn_OverRide*65535; color.blue  = Group->oflBlu_OverRide*65535; 
        colorbutton = gtk_color_button_new(); gtk_widget_show(colorbutton);
        gtk_color_button_set_color( GTK_COLOR_BUTTON(colorbutton), &color );
        gtk_color_button_set_use_alpha( GTK_COLOR_BUTTON(colorbutton), TRUE );
        gtk_color_button_set_alpha( GTK_COLOR_BUTTON(colorbutton), Group->oflAlf_OverRide*65535 );
        gtk_table_attach( GTK_TABLE(table1), colorbutton, Col, Col+1, Row+1, Row+2, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0 );
        gtk_widget_set_size_request( colorbutton, 30, 25 );
        gtk_color_button_set_use_alpha( GTK_COLOR_BUTTON(colorbutton), TRUE );
        g_signal_connect( G_OBJECT(colorbutton), "color-set", G_CALLBACK(ChangeSatColor), GINT_TO_POINTER(i*100 + 23) );
        ++Col;



        // for Orbit FL Footpoints
        checkbutton = gtk_check_button_new_with_mnemonic (""); gtk_widget_show (checkbutton);
        gtk_table_attach (GTK_TABLE (table1), checkbutton, Col, Col+1, Row, Row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
        gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(checkbutton), Group->Sat[0].DrawOrbitFLFootpoints );
        g_signal_connect( G_OBJECT( checkbutton ), "toggled", G_CALLBACK( ToggleOrbitOptions ), GINT_TO_POINTER( i*100 + 24 ) );
        //++Col;
        color.red   = Group->ofpRed_OverRide*65535; color.green = Group->ofpGrn_OverRide*65535; color.blue  = Group->ofpBlu_OverRide*65535; 
        colorbutton = gtk_color_button_new(); gtk_widget_show(colorbutton);
        gtk_color_button_set_color( GTK_COLOR_BUTTON(colorbutton), &color );
        gtk_color_button_set_use_alpha( GTK_COLOR_BUTTON(colorbutton), TRUE );
        gtk_color_button_set_alpha( GTK_COLOR_BUTTON(colorbutton), Group->ofpAlf_OverRide*65535 );
        gtk_table_attach( GTK_TABLE(table1), colorbutton, Col, Col+1, Row+1, Row+2, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0 );
        gtk_widget_set_size_request( colorbutton, 30, 25 );
        gtk_color_button_set_use_alpha( GTK_COLOR_BUTTON(colorbutton), TRUE );
        g_signal_connect( G_OBJECT(colorbutton), "color-set", G_CALLBACK(ChangeSatColor), GINT_TO_POINTER(i*100 + 24) );
        ++Col;


        // Override Orbit Colors
        checkbutton = gtk_check_button_new_with_mnemonic (""); gtk_widget_show (checkbutton);
        gtk_table_attach (GTK_TABLE (table1), checkbutton, Col, Col+1, Row+1, Row+2, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
        gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(checkbutton),  Group->OverRideOrbitColors );
        g_signal_connect( G_OBJECT( checkbutton ), "toggled", G_CALLBACK( ToggleOrbitOptions ), GINT_TO_POINTER( i*100 + 26 ) );
        ++Col;


        // for % Orbit to Draw
        spinbutton_adj = gtk_adjustment_new (100, 1, 2000, 1, 10, 0);
        spinbutton = gtk_spin_button_new (GTK_ADJUSTMENT (spinbutton_adj), 1, 0); gtk_widget_show (spinbutton);
        gtk_table_attach (GTK_TABLE (table1), spinbutton, Col, Col+1, Row, Row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
        gtk_spin_button_set_numeric (GTK_SPIN_BUTTON (spinbutton), TRUE);
        gtk_spin_button_set_value( GTK_SPIN_BUTTON (spinbutton), Group->Sat[0].oPeriodFrac );
        g_signal_connect( G_OBJECT( spinbutton ), "value-changed", G_CALLBACK( ChangeOrbitOptions ), GINT_TO_POINTER( i*100 + 25 ) );
        ++Col;



        vseparator = gtk_vseparator_new (); gtk_widget_show (vseparator);
        gtk_table_attach (GTK_TABLE (table1), vseparator, Col, Col+1, Row, Row+1, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (GTK_FILL), 0, 0);
        gtk_widget_set_size_request (vseparator, 2, -1);
        ++Col;

        /******************************
         *
         *      Streak Quantities
         *
         ******************************/
        // for Streak
        checkbutton = gtk_check_button_new_with_mnemonic (""); gtk_widget_show (checkbutton);
        gtk_table_attach (GTK_TABLE (table1), checkbutton, Col, Col+1, Row, Row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
        gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(checkbutton), Group->Sat[0].DrawStreak );
        g_signal_connect( G_OBJECT( checkbutton ), "toggled", G_CALLBACK( ToggleOrbitOptions ), GINT_TO_POINTER( i*100 + 30 ) );
        //++Col;
        color.red   = Group->sRed_OverRide*65535; color.green = Group->sGrn_OverRide*65535; color.blue  = Group->sBlu_OverRide*65535; 
        colorbutton = gtk_color_button_new(); gtk_widget_show(colorbutton);
        gtk_color_button_set_color( GTK_COLOR_BUTTON(colorbutton), &color );
        gtk_color_button_set_use_alpha( GTK_COLOR_BUTTON(colorbutton), TRUE );
        gtk_color_button_set_alpha( GTK_COLOR_BUTTON(colorbutton), Group->sAlf_OverRide*65535 );
        gtk_table_attach( GTK_TABLE(table1), colorbutton, Col, Col+1, Row+1, Row+2, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0 );
        gtk_widget_set_size_request( colorbutton, 30, 25 );
        gtk_color_button_set_use_alpha( GTK_COLOR_BUTTON(colorbutton), TRUE );
        g_signal_connect( G_OBJECT(colorbutton), "color-set", G_CALLBACK(ChangeSatColor), GINT_TO_POINTER(i*100 + 30) );
        ++Col;

        // for Ground Path of Streak
        checkbutton = gtk_check_button_new_with_mnemonic (""); gtk_widget_show (checkbutton);
        gtk_table_attach (GTK_TABLE (table1), checkbutton, Col, Col+1, Row, Row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
        gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(checkbutton), Group->Sat[0].DrawGroundPathOfStreak );
        g_signal_connect( G_OBJECT( checkbutton ), "toggled", G_CALLBACK( ToggleOrbitOptions ), GINT_TO_POINTER( i*100 + 31 ) );
        //++Col;
        color.red   = Group->sgpRed_OverRide*65535; color.green = Group->sgpGrn_OverRide*65535; color.blue  = Group->sgpBlu_OverRide*65535; 
        colorbutton = gtk_color_button_new(); gtk_widget_show(colorbutton);
        gtk_color_button_set_color( GTK_COLOR_BUTTON(colorbutton), &color );
        gtk_color_button_set_use_alpha( GTK_COLOR_BUTTON(colorbutton), TRUE );
        gtk_color_button_set_alpha( GTK_COLOR_BUTTON(colorbutton), Group->sgpAlf_OverRide*65535 );
        gtk_table_attach( GTK_TABLE(table1), colorbutton, Col, Col+1, Row+1, Row+2, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0 );
        gtk_widget_set_size_request( colorbutton, 30, 25 );
        gtk_color_button_set_use_alpha( GTK_COLOR_BUTTON(colorbutton), TRUE );
        g_signal_connect( G_OBJECT(colorbutton), "color-set", G_CALLBACK(ChangeSatColor), GINT_TO_POINTER(i*100 + 31) );
        ++Col;

        // for Streak to Ground Lines
        checkbutton = gtk_check_button_new_with_mnemonic (""); gtk_widget_show (checkbutton);
        gtk_table_attach (GTK_TABLE (table1), checkbutton, Col, Col+1, Row, Row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
        gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(checkbutton), Group->Sat[0].DrawStreakToGroundLines );
        g_signal_connect( G_OBJECT( checkbutton ), "toggled", G_CALLBACK( ToggleOrbitOptions ), GINT_TO_POINTER( i*100 + 32 ) );
        //++Col;
        color.red   = Group->sglRed_OverRide*65535; color.green = Group->sglGrn_OverRide*65535; color.blue  = Group->sglBlu_OverRide*65535; 
        colorbutton = gtk_color_button_new(); gtk_widget_show(colorbutton);
        gtk_color_button_set_color( GTK_COLOR_BUTTON(colorbutton), &color );
        gtk_color_button_set_use_alpha( GTK_COLOR_BUTTON(colorbutton), TRUE );
        gtk_color_button_set_alpha( GTK_COLOR_BUTTON(colorbutton), Group->sglAlf_OverRide*65535 );
        gtk_table_attach( GTK_TABLE(table1), colorbutton, Col, Col+1, Row+1, Row+2, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0 );
        gtk_widget_set_size_request( colorbutton, 30, 25 );
        gtk_color_button_set_use_alpha( GTK_COLOR_BUTTON(colorbutton), TRUE );
        g_signal_connect( G_OBJECT(colorbutton), "color-set", G_CALLBACK(ChangeSatColor), GINT_TO_POINTER(i*100 + 32) );
        ++Col;

        // for Streak Field Lines
        checkbutton = gtk_check_button_new_with_mnemonic (""); gtk_widget_show (checkbutton);
        gtk_table_attach (GTK_TABLE (table1), checkbutton, Col, Col+1, Row, Row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
        gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(checkbutton), Group->Sat[0].DrawStreakFieldLines );
        g_signal_connect( G_OBJECT( checkbutton ), "toggled", G_CALLBACK( ToggleOrbitOptions ), GINT_TO_POINTER( i*100 + 33 ) );
        //++Col;
        color.red   = Group->sflRed_OverRide*65535; color.green = Group->sflGrn_OverRide*65535; color.blue  = Group->sflBlu_OverRide*65535; 
        colorbutton = gtk_color_button_new(); gtk_widget_show(colorbutton);
        gtk_color_button_set_color( GTK_COLOR_BUTTON(colorbutton), &color );
        gtk_color_button_set_use_alpha( GTK_COLOR_BUTTON(colorbutton), TRUE );
        gtk_color_button_set_alpha( GTK_COLOR_BUTTON(colorbutton), Group->sflAlf_OverRide*65535 );
        gtk_table_attach( GTK_TABLE(table1), colorbutton, Col, Col+1, Row+1, Row+2, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0 );
        gtk_widget_set_size_request( colorbutton, 30, 25 );
        gtk_color_button_set_use_alpha( GTK_COLOR_BUTTON(colorbutton), TRUE );
        g_signal_connect( G_OBJECT(colorbutton), "color-set", G_CALLBACK(ChangeSatColor), GINT_TO_POINTER(i*100 + 33) );
        ++Col;

        // for Streak FL Footpoints
        checkbutton = gtk_check_button_new_with_mnemonic (""); gtk_widget_show (checkbutton);
        gtk_table_attach (GTK_TABLE (table1), checkbutton, Col, Col+1, Row, Row+1, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
        gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(checkbutton), Group->Sat[0].DrawStreakFLFootpoints );
        g_signal_connect( G_OBJECT( checkbutton ), "toggled", G_CALLBACK( ToggleOrbitOptions ), GINT_TO_POINTER( i*100 + 34 ) );
        //++Col;
        color.red   = Group->sfpRed_OverRide*65535; color.green = Group->sfpGrn_OverRide*65535; color.blue  = Group->sfpBlu_OverRide*65535; 
        colorbutton = gtk_color_button_new(); gtk_widget_show(colorbutton);
        gtk_color_button_set_color( GTK_COLOR_BUTTON(colorbutton), &color );
        gtk_color_button_set_use_alpha( GTK_COLOR_BUTTON(colorbutton), TRUE );
        gtk_color_button_set_alpha( GTK_COLOR_BUTTON(colorbutton), Group->sfpAlf_OverRide*65535 );
        gtk_table_attach( GTK_TABLE(table1), colorbutton, Col, Col+1, Row+1, Row+2, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0 );
        gtk_widget_set_size_request( colorbutton, 30, 25 );
        gtk_color_button_set_use_alpha( GTK_COLOR_BUTTON(colorbutton), TRUE );
        g_signal_connect( G_OBJECT(colorbutton), "color-set", G_CALLBACK(ChangeSatColor), GINT_TO_POINTER(i*100 + 34) );
        ++Col;



        // Override Streak Colors
        checkbutton = gtk_check_button_new_with_mnemonic (""); gtk_widget_show (checkbutton);
        gtk_table_attach (GTK_TABLE (table1), checkbutton, Col, Col+1, Row+1, Row+2, (GtkAttachOptions) (0), (GtkAttachOptions) (0), 0, 0);
        gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(checkbutton),  Group->OverRideStreakColors );
        g_signal_connect( G_OBJECT( checkbutton ), "toggled", G_CALLBACK( ToggleOrbitOptions ), GINT_TO_POINTER( i*100 + 36 ) );
        ++Col;



        // for % Orbit to Draw
        spinbutton_adj = gtk_adjustment_new (25, 1, 2000, 1, 10, 0);
        spinbutton = gtk_spin_button_new (GTK_ADJUSTMENT (spinbutton_adj), 1, 0); gtk_widget_show (spinbutton);
        gtk_table_attach (GTK_TABLE (table1), spinbutton, Col, Col+1, Row, Row+1, (GtkAttachOptions) (GTK_EXPAND), (GtkAttachOptions) (0), 0, 0);
        gtk_spin_button_set_numeric (GTK_SPIN_BUTTON (spinbutton), TRUE);
        gtk_spin_button_set_value( GTK_SPIN_BUTTON (spinbutton), Group->Sat[0].sPeriodFrac );
        g_signal_connect( G_OBJECT( spinbutton ), "value-changed", G_CALLBACK( ChangeOrbitOptions ), GINT_TO_POINTER( i*100 + 35 ) );
        ++Col;



        hseparator = gtk_hseparator_new (); gtk_widget_show (hseparator);
        gtk_table_attach( GTK_TABLE(table1), hseparator, 0, TotalCols, Row+2, Row+3, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);
        gtk_widget_set_size_request (hseparator, -1, 4);

        vseparator = gtk_vseparator_new (); gtk_widget_show (vseparator);
        gtk_table_attach (GTK_TABLE (table1), vseparator, Col, Col+1, Row, Row+1, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (GTK_FILL), 0, 0);
        gtk_widget_set_size_request (vseparator, 2, -1);
        ++Col;

        // for More Options button
//        button = gtk_button_new_with_mnemonic (_("More Options")); gtk_widget_show (button);
//        gtk_table_attach (GTK_TABLE (table1), button, 24, 25, Row, Row+1, (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);
//        ++Col;

        g = g->Next;
        ++i; 
        Row += 3;

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

    spinbutton_adj = gtk_adjustment_new( IllumFL_ka, 0.0, 5.0, 0.05, .1, 0);
    spinbutton = gtk_spin_button_new( GTK_ADJUSTMENT(spinbutton_adj), 1, 3); gtk_widget_show(spinbutton);
    gtk_box_pack_start( GTK_BOX(hbox1), spinbutton, FALSE, FALSE, 0 );
    gtk_spin_button_set_numeric( GTK_SPIN_BUTTON(spinbutton), TRUE );
    gtk_spin_button_set_value( GTK_SPIN_BUTTON(spinbutton), IllumFL_ka );
    g_signal_connect( G_OBJECT( spinbutton ), "value-changed", G_CALLBACK( ChangeIllumFLParams ), GINT_TO_POINTER( 1 ) );

    spinbutton_adj = gtk_adjustment_new( IllumFL_kd, 0.0, 5.0, 0.05, .1, 0);
    spinbutton = gtk_spin_button_new( GTK_ADJUSTMENT(spinbutton_adj), 1, 3); gtk_widget_show(spinbutton);
    gtk_box_pack_start( GTK_BOX(hbox1), spinbutton, FALSE, FALSE, 0 );
    gtk_spin_button_set_numeric( GTK_SPIN_BUTTON(spinbutton), TRUE );
    gtk_spin_button_set_value( GTK_SPIN_BUTTON(spinbutton), IllumFL_kd );
    g_signal_connect( G_OBJECT( spinbutton ), "value-changed", G_CALLBACK( ChangeIllumFLParams ), GINT_TO_POINTER( 2 ) );

    spinbutton_adj = gtk_adjustment_new( IllumFL_ks, 0.0, 5.0, 0.05, .1, 0);
    spinbutton = gtk_spin_button_new( GTK_ADJUSTMENT(spinbutton_adj), 1, 3); gtk_widget_show(spinbutton);
    gtk_box_pack_start( GTK_BOX(hbox1), spinbutton, FALSE, FALSE, 0 );
    gtk_spin_button_set_numeric( GTK_SPIN_BUTTON(spinbutton), TRUE );
    gtk_spin_button_set_value( GTK_SPIN_BUTTON(spinbutton), IllumFL_ks );
    g_signal_connect( G_OBJECT( spinbutton ), "value-changed", G_CALLBACK( ChangeIllumFLParams ), GINT_TO_POINTER( 3 ) );

    spinbutton_adj = gtk_adjustment_new( IllumFL_n, 0.0, 256.0, 0.1, 1.0, 0);
    spinbutton = gtk_spin_button_new( GTK_ADJUSTMENT(spinbutton_adj), 1, 5); gtk_widget_show(spinbutton);
    gtk_box_pack_start( GTK_BOX(hbox1), spinbutton, FALSE, FALSE, 0 );
    gtk_spin_button_set_numeric( GTK_SPIN_BUTTON(spinbutton), TRUE );
    gtk_spin_button_set_value( GTK_SPIN_BUTTON(spinbutton), IllumFL_n );
    g_signal_connect( G_OBJECT( spinbutton ), "value-changed", G_CALLBACK( ChangeIllumFLParams ), GINT_TO_POINTER( 4 ) );

    spinbutton_adj = gtk_adjustment_new( IllumFL_w, 0.0, 10.0, 0.05, .1, 0);
    spinbutton = gtk_spin_button_new( GTK_ADJUSTMENT(spinbutton_adj), 1, 3); gtk_widget_show(spinbutton);
    gtk_box_pack_start( GTK_BOX(hbox1), spinbutton, FALSE, FALSE, 0 );
    gtk_spin_button_set_numeric( GTK_SPIN_BUTTON(spinbutton), TRUE );
    gtk_spin_button_set_value( GTK_SPIN_BUTTON(spinbutton), IllumFL_w );
    g_signal_connect( G_OBJECT( spinbutton ), "value-changed", G_CALLBACK( ChangeIllumFLParams ), GINT_TO_POINTER( 5 ) );


    button = gtk_button_new_with_mnemonic (_("Delete Group")); gtk_widget_show (button);
    gtk_box_pack_end (GTK_BOX (hbox1), button, FALSE, FALSE, 0);

    button = gtk_button_new_with_mnemonic (_("Edit Group")); gtk_widget_show (button);
    gtk_box_pack_end (GTK_BOX (hbox1), button, FALSE, FALSE, 0);

    button = gtk_button_new_with_mnemonic (_("Create Group")); gtk_widget_show (button);
    gtk_box_pack_end (GTK_BOX (hbox1), button, FALSE, FALSE, 0);


    return vbox1;




}

