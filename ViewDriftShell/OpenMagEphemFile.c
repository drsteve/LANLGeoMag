#include "ViewDriftShell.h"

#define GLADE_HOOKUP_OBJECT(component,widget,name ) \
  g_object_set_data_full( G_OBJECT( component ), name, \
    gtk_widget_ref( widget ),( GDestroyNotify ) gtk_widget_unref )

#define GLADE_HOOKUP_OBJECT_NO_REF(component,widget,name ) \
  g_object_set_data( G_OBJECT( component ), name, widget )


gboolean OpenMyFile( GtkWidget *widget, gpointer data ) {

    GError      *error = NULL;
    gchar       *Filename;
    gchar       *ExtText;
    int         x, y, width, height, depth;

    if ( (Filename  = gtk_file_chooser_get_filename( GTK_FILE_CHOOSER( widget ) )) == NULL ) return( FALSE );
    printf( "Open Filename = %s\n", Filename);


    gtk_widget_destroy( widget );

    return( TRUE );

}


GtkWidget *CreateOpenMagEphemFileDialog( void ) {

    GtkWidget *filechooserdialog1;
    GtkWidget *dialog_vbox1;
    GtkWidget *hbox1;
    GtkWidget *dialog_action_area1;
    GtkWidget *button1;
    GtkWidget *button2;


    filechooserdialog1 = gtk_file_chooser_dialog_new( "Open MagEphem File", NULL, GTK_FILE_CHOOSER_ACTION_OPEN, NULL );
    //g_object_set( filechooserdialog1, "local-only", FALSE, NULL );
    gtk_window_set_type_hint( GTK_WINDOW(filechooserdialog1), GDK_WINDOW_TYPE_HINT_DIALOG );
    gtk_file_chooser_set_do_overwrite_confirmation( GTK_FILE_CHOOSER(filechooserdialog1), TRUE );

    dialog_vbox1 = GTK_DIALOG( filechooserdialog1 )->vbox;
    gtk_widget_show( dialog_vbox1 );

    hbox1 = gtk_hbox_new( FALSE, 52 );
    gtk_widget_show( hbox1 );
    gtk_box_pack_start( GTK_BOX(dialog_vbox1), hbox1, FALSE, TRUE, 0 );

    dialog_action_area1 = GTK_DIALOG(filechooserdialog1)->action_area;
    gtk_widget_show( dialog_action_area1 );
    gtk_button_box_set_layout( GTK_BUTTON_BOX(dialog_action_area1), GTK_BUTTONBOX_END );

    button1 = gtk_button_new_from_stock( "gtk-cancel" );
    gtk_widget_show( button1 );
    gtk_dialog_add_action_widget( GTK_DIALOG(filechooserdialog1), button1, GTK_RESPONSE_CANCEL );
    GTK_WIDGET_SET_FLAGS( button1, GTK_CAN_DEFAULT );
    g_signal_connect_swapped( GTK_OBJECT( button1 ), "clicked", G_CALLBACK( gtk_widget_destroy), GTK_OBJECT( filechooserdialog1 ));

    button2 = gtk_button_new_from_stock( "gtk-open" );
    gtk_widget_show( button2 );
    gtk_dialog_add_action_widget( GTK_DIALOG(filechooserdialog1), button2, GTK_RESPONSE_OK );
    GTK_WIDGET_SET_FLAGS( button2, GTK_CAN_DEFAULT );
    g_signal_connect_swapped( GTK_OBJECT( button2 ), "clicked", G_CALLBACK( OpenMyFile ), GTK_OBJECT( filechooserdialog1 ));


    gtk_widget_grab_default( button2 );

    return( filechooserdialog1 );

}

