#ifndef VIEWDRIFTSHELL_H
#define VIEWDRIFTSHELL_H

#define GL_GLEXT_PROTOTYPES 1

#include <gtk/gtkgl.h>

#include <GL/glew.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glext.h>
#include <GL/glx.h>

#include "support.h"
#include <Lgm/Lgm_MagEphemInfo.h>
#include "PsdAssim.h"
#include <Lgm/Lgm_Quat.h>
#include <Lgm/Lgm_Sgp.h>
#include "Atmosphere.h"
#include <Lgm/Lgm_DynamicMemory.h>
#include "NamedMaterials.h"

#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <math.h>

#include <gtk/gtk.h>
#include <gdk/gdkkeysyms.h>


#ifdef G_OS_WIN32
#define WIN32_LEAN_AND_MEAN 1
#include <windows.h>
#endif



#include "SatSelector.h"


#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
//#define GSL_INTERP  gsl_interp_linear
#define GSL_INTERP  gsl_interp_akima
//#define GSL_INTERP  gsl_interp_cspline


/*
 * For compiling/loading/etc. custom
 * vertex and fragment shaders
 */
//#include <Cg/cg.h>
//#include <Cg/cgGL.h>

#define SHADER_PATH "/home/mgh/DREAM/Dream/Dream/Shaders"

#define TIME_RESET_BACKWARD_TO_START    -3
#define TIME_PLAY_BACKWARD              -2
#define TIME_STEP_BACKWARD              -1
#define TIME_STOP                        0
#define TIME_STEP_FOREWARD               1
#define TIME_PLAY_FOREWARD               2
#define TIME_RESET_FOREWARD_TO_END       3
#define TIME_REALTIMEPLAY                4


/*
 * Function Prototypes
 */
void raise_OpenMagEphemFileDialog( gpointer callback_data, guint callback_action, GtkWidget *menu_item ) ;
void raise_SaveRasterFileDialog( gpointer callback_data, guint callback_action, GtkWidget *menu_item ) ;
void GoFullScreen( gpointer callback_data, guint callback_action, GtkWidget *menu_item );
void GoNormalScreen( gpointer callback_data, guint callback_action, GtkWidget *menu_item );
GtkWidget *CreateOpenMagEphemFileDialog( void );
GtkWidget *CreateSaveRasterFileDialog( void );




/*
 * Global defs
 */
#ifdef MAIN
/*
 *  Make aInfo Global. We could get around this,
 *  but it adds complexity.
 */
GtkWidget   *drawing_area;
GtkWidget   *ViewDriftShellWindow;
long int    nFrames, cFrame, nFramesLeft;
#else
extern GtkWidget   *drawing_area;
extern GtkWidget   *ViewDriftShellWindow;
extern long int    nFrames, cFrame, nFramesLeft;
#endif



#endif


