#include "ViewDriftShell.h"

#define GLADE_HOOKUP_OBJECT(component,widget,name ) \
  g_object_set_data_full( G_OBJECT( component ), name, \
    gtk_widget_ref( widget ),( GDestroyNotify ) gtk_widget_unref )

#define GLADE_HOOKUP_OBJECT_NO_REF(component,widget,name ) \
  g_object_set_data( G_OBJECT( component ), name, widget )


GtkWidget *checkbutton1;
GtkWidget *combobox1;

gboolean SavePng( GtkWidget *widget, gpointer data ) {

    GdkPixbuf   *pixbuf;
    GError      *error = NULL;
    gboolean    UseExt;
    gchar       *Filename;
    gchar       *Filename2;
    gchar       *ExtText;
    int         x, y, width, height, depth;

    Filename  = gtk_file_chooser_get_filename( GTK_FILE_CHOOSER( widget ) );
    if (Filename == NULL){
        return( FALSE );
    }
    UseExt    = gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON( checkbutton1 ) );
    ExtText = gtk_combo_box_get_active_text( GTK_COMBO_BOX( combobox1 ) );
    Filename2 = (gchar *) calloc( strlen(Filename)+strlen(ExtText)+4, sizeof( gchar ) );
    
    printf( "Filename = %s\n", Filename);
    printf( "UseExt = %d\n", UseExt);
    printf( "ExtText = %s\n", ExtText);
    
    if ( UseExt ) {
        sprintf(Filename2, "%s.%s", Filename, ExtText );
    } else {
        sprintf(Filename2, "%s", Filename );
    }
    printf( "Filename2 = %s\n", Filename2);


    gdk_window_get_geometry( drawing_area->window, &x, &y, &width, &height, &depth );
    pixbuf = gdk_pixbuf_get_from_drawable( NULL, GDK_DRAWABLE(drawing_area->window), NULL, 0, 0, 0, 0, width, height);
    
    if        ( !strcmp( ExtText, "png" ) ) {
        gdk_pixbuf_save( pixbuf, Filename2, "png", &error, "compression", "0", NULL );
    } else if ( !strcmp( ExtText, "jpg" ) ) {
        gdk_pixbuf_save( pixbuf, Filename2, "jpeg", &error, "quality", "100", NULL );
    } else if ( !strcmp( ExtText, "tiff" ) ) {
        gdk_pixbuf_save( pixbuf, Filename2, "tiff", &error, NULL );
    } else if ( !strcmp( ExtText, "bmp" ) ) {
        gdk_pixbuf_save( pixbuf, Filename2, "bmp", &error, NULL );
    }

    


    g_object_unref( pixbuf );
    g_free( Filename );
    g_free( Filename2 );
    g_free( ExtText );
    //gtk_widget_destroy( widget );

    return( TRUE );

}


GtkWidget *CreateSaveRasterFileDialog( void ) {

    GtkWidget *filechooserdialog1;
    GtkWidget *dialog_vbox1;
    GtkWidget *hbox1;
    GtkWidget *dialog_action_area1;
    GtkWidget *button1;
    GtkWidget *button2;


    //filechooserdialog1 = gtk_file_chooser_dialog_new( "Save Raster Image to File", NULL, GTK_FILE_CHOOSER_ACTION_SAVE, NULL );
    filechooserdialog1 = gtk_file_chooser_dialog_new ( "Save Raster Image to File", NULL,
                      //parent_window,
                      GTK_FILE_CHOOSER_ACTION_OPEN,
                      GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
                      GTK_STOCK_OPEN, GTK_RESPONSE_ACCEPT,
                      NULL);
    //g_object_set( filechooserdialog1, "local-only", FALSE, NULL );
    gtk_window_set_type_hint( GTK_WINDOW(filechooserdialog1), GDK_WINDOW_TYPE_HINT_DIALOG );
    gtk_file_chooser_set_do_overwrite_confirmation( GTK_FILE_CHOOSER(filechooserdialog1), TRUE );

    dialog_vbox1 = GTK_DIALOG( filechooserdialog1 )->vbox;
    gtk_widget_show( dialog_vbox1 );

    hbox1 = gtk_hbox_new( FALSE, 52 );
    gtk_widget_show( hbox1 );
    gtk_box_pack_start( GTK_BOX(dialog_vbox1), hbox1, FALSE, TRUE, 0 );

    checkbutton1 = gtk_check_button_new_with_mnemonic( _("Append Extension Automatically") );
    gtk_widget_show( checkbutton1);
    gtk_box_pack_start( GTK_BOX (hbox1), checkbutton1, FALSE, FALSE, 20 );
    gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON (checkbutton1), TRUE );

    combobox1 = gtk_combo_box_new_text( );
    gtk_widget_show( combobox1 );
    gtk_box_pack_start( GTK_BOX(hbox1), combobox1, TRUE, TRUE, 20 );
    gtk_combo_box_append_text( GTK_COMBO_BOX(combobox1), _("png") );
    gtk_combo_box_append_text( GTK_COMBO_BOX(combobox1), _("jpg") );
    //gtk_combo_box_append_text( GTK_COMBO_BOX(combobox1), _("tiff") );
    gtk_combo_box_append_text( GTK_COMBO_BOX(combobox1), _("bmp") );
    gtk_combo_box_set_active( GTK_COMBO_BOX(combobox1), 0 );


    dialog_action_area1 = GTK_DIALOG(filechooserdialog1)->action_area;
    gtk_widget_show( dialog_action_area1 );
    gtk_button_box_set_layout( GTK_BUTTON_BOX(dialog_action_area1), GTK_BUTTONBOX_END );

    button1 = gtk_button_new_from_stock( "gtk-cancel" );
    gtk_widget_show( button1 );
    gtk_dialog_add_action_widget( GTK_DIALOG(filechooserdialog1), button1, GTK_RESPONSE_CANCEL );
    GTK_WIDGET_SET_FLAGS( button1, GTK_CAN_DEFAULT );
//    g_signal_connect( GTK_OBJECT( button1 ), "clicked", G_CALLBACK( Destroy ), GINT_TO_POINTER( i ) );
    g_signal_connect_swapped( GTK_OBJECT( button1 ), "clicked", G_CALLBACK( gtk_widget_destroy), GTK_OBJECT( filechooserdialog1 ));

    button2 = gtk_button_new_from_stock( "gtk-save" );
    gtk_widget_show( button2 );
    gtk_dialog_add_action_widget( GTK_DIALOG(filechooserdialog1), button2, GTK_RESPONSE_OK );
    GTK_WIDGET_SET_FLAGS( button2, GTK_CAN_DEFAULT );
    g_signal_connect_swapped( GTK_OBJECT( button2 ), "clicked", G_CALLBACK( SavePng ), GTK_OBJECT( filechooserdialog1 ));


    gtk_widget_grab_default( button2 );

    return( filechooserdialog1 );

}

