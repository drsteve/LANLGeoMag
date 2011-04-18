#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct CoordInfo {

    double  x0, x1;
    double  px0, px1;
    double  sx;
    double  y0, y1;
    double  py0, py1;
    double  sy;

    int     ImageX;
    int     ImageY;
    int     ImageWidth;
    int     ImageHeight;

    char    ImageFile[4096];

    int     Image_n;
    int     Box_n;
    int     Path_n;
    int     Tspan_n;

} CoordInfo;

CoordInfo *Create_CoordInfo( ){

    CoordInfo *g = (CoordInfo *)calloc( 1, sizeof(*g) );

    g->Image_n = 0;
    g->Box_n   = 0;
    g->Path_n = 0;
    g->Tspan_n = 0;

    return( g );


}

void SetupCoords( double x0, double x1, double px0, double px1, double y0, double y1, double py0, double py1, CoordInfo *g ) {

    g->x0  = x0;
    g->x1  = x1;
    g->px0 = px0;
    g->px1 = px1;
    g->sx  = (g->px1 - g->px0)/(g->x1 - g->x0);

    g->y0  = y0;
    g->y1  = y1;
    g->py0 = py0;
    g->py1 = py1;
    g->sy  = (g->py1 - g->py0)/(g->y1 - g->y0);

}

void PageCoords( double x, double y, double *px, double *py, CoordInfo *g ) {

    *px = g->px0 + g->sx*(x - g->x0);
    *py = (g->py0 + g->sy*(y - g->y0));

}

void SvgImage( FILE *fp, int ImageX, int ImageY, int ImageWidth, int ImageHeight, char *ImageFile, CoordInfo *g ) {

    g->ImageX      = ImageX; 
    g->ImageY      = ImageY;
    g->ImageWidth  = ImageWidth; 
    g->ImageHeight = ImageHeight;
    strcpy( g->ImageFile, ImageFile );

    // Flux versus E and Alpha image
    fprintf( fp, "        <image\n" );
//    fprintf( fp, "           image-rendering=\"optimizeSpeed\"\n" );
    fprintf( fp, "           xlink:href=\"%s\"\n", g->ImageFile );
    fprintf( fp, "           x=\"%d\"\n", g->ImageX );
    fprintf( fp, "           y=\"%d\"\n", g->ImageY );
    fprintf( fp, "           width=\"%d\"\n", g->ImageWidth );
    fprintf( fp, "           height=\"%d\"\n", g->ImageHeight );
    fprintf( fp, "           id=\"Image_%04d\"\n", g->Image_n++ );
    fprintf( fp, "           style=\"fill:none;stroke:#950000;stroke-opacity:1\"\n" );
    fprintf( fp, "        />\n" );


    // White rectangle around FLUX versus E and a image
    fprintf( fp, "        <rect\n"  );
    fprintf( fp, "           x=\"%d\"\n", g->ImageX );
    fprintf( fp, "           y=\"%d\"\n", g->ImageY );
    fprintf( fp, "           width=\"%d\"\n", g->ImageWidth );
    fprintf( fp, "           height=\"%d\"\n", g->ImageHeight );
    fprintf( fp, "           id=\"Box_%04d\"\n", g->Box_n++ );
    fprintf( fp, "           style=\"fill:none;stroke:#ffffff;stroke-width:3;stroke-miterlimit:3.4857142;stroke-opacity:1;stroke-dasharray:none\"\n" );
    fprintf( fp, "        />\n" );
    return;
}

SvgXticks( FILE *fp, double x0, double x1, double xinc, CoordInfo *g ) {

    double  x, y, px, py;

    // X-ticks
    for ( y=g->y0, x=x0; x<= x1; x += xinc ){
        PageCoords( x, y, &px, &py, g );
        fprintf( fp, "       <path\n" );
        fprintf( fp, "       d=\"m %g,%g 0,-10\"\n", px, py );
        fprintf( fp, "       id=\"%04d\"\n", g->Path_n++ );
        fprintf( fp, "       style=\"fill:none;stroke:#ffffff;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\" />\n" );
    }
    for ( y=g->y1, x=x0; x<= x1; x += xinc ){
        PageCoords( x, y, &px, &py, g );
        fprintf( fp, "       <path\n" );
        fprintf( fp, "       d=\"m %g,%g 0,10\"\n", px, py );
        fprintf( fp, "       id=\"%04d\"\n", g->Path_n++ );
        fprintf( fp, "       style=\"fill:none;stroke:#ffffff;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\" />\n" );
    }

}

SvgYticks( FILE *fp, double y0, double y1, double yinc, CoordInfo *g ) {

    double  x, y, px, py;

    // Y-ticks
    for ( x=g->x0, y=y0; y<= y1; y += yinc ){
        PageCoords( x, y, &px, &py, g );
        fprintf( fp, "       <path\n" );
        fprintf( fp, "       d=\"m %g,%g 10,0\"\n", px, py );
        fprintf( fp, "       id=\"%04d\"\n", g->Path_n++ );
        fprintf( fp, "       style=\"fill:none;stroke:#ffffff;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\" />\n" );
    }
    for ( x=g->x1, y=y0; y<= y1; y += yinc ){
        PageCoords( x, y, &px, &py, g );
        fprintf( fp, "       <path\n" );
        fprintf( fp, "       d=\"m %g,%g -10,0\"\n", px, py );
        fprintf( fp, "       id=\"%04d\"\n", g->Path_n++ );
        fprintf( fp, "       style=\"fill:none;stroke:#ffffff;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\" />\n" );
    }

}

SvgTitle( FILE *fp, int FontSize, double px, double py, char *Title, CoordInfo *g ) {


    // Title
    fprintf( fp, "       <text\n" );
    fprintf( fp, "          xml:space=\"preserve\"\n");
        fprintf( fp, "      style=\"font-size:%dpx;font-style:normal;font-variant:normal;font-weight:normal;font-stretch:normal;vertical-align:middle;text-align:left;line-height:100%;letter-spacing:0px;word-spacing:0px;writing-mode:lr-tb;text-anchor:middle;fill:#ffffff;fill-opacity:1;stroke:none;font-family:Myriad Pro;-inkscape-font-specification:Myriad Pro\"\n", FontSize);
    fprintf( fp, "          x=\"%g\"\n", px);
    fprintf( fp, "          y=\"%g\"\n", py);
    fprintf( fp, "          ><tspan \n");
    fprintf( fp, "          id=\"%04d\"\n", g->Tspan_n++);
    fprintf( fp, "          x=\"%g\"\n", px );
    fprintf( fp, "          y=\"%g\">%s</tspan></text>\n", py, Title);


}

SvgXnumbers( FILE *fp, double x0, double x1, double xinc, double *a, char *Title, CoordInfo *g ) {

    double  x, y, px, py;
    int     i;

    // X-Numbers
    for ( y=g->y0, x=x0; x<= x1; x += xinc ){
        i = (int)(x+0.5);
        PageCoords( x, y, &px, &py, g );
        fprintf( fp, "       <text\n" );
        fprintf( fp, "          xml:space=\"preserve\"\n");
        fprintf( fp, "          style=\"font-size:12px;font-style:normal;font-variant:normal;font-weight:normal;font-stretch:normal;vertical-align:middle;text-align:left;baseline-shift:-1em; line-height:100%;letter-spacing:0px;word-spacing:0px;writing-mode:lr-tb;text-anchor:middle;fill:#ffffff;fill-opacity:1;stroke:none;font-family:Myriad Pro;-inkscape-font-specification:Myriad Pro\"\n");
        fprintf( fp, "          x=\"%g\"\n", px);
        fprintf( fp, "          y=\"%g\"\n", py);
        fprintf( fp, "          ><tspan \n");
        fprintf( fp, "          sodipodi:role=\"line\"\n");
        fprintf( fp, "          id=\"%04d\"\n", g->Tspan_n++);
        fprintf( fp, "          x=\"%g\"\n", px );
        fprintf( fp, "          y=\"%g\">%.2g</tspan></text>\n", py, a[i]);
    }

    // Title
    x = 0.5*(g->x0+g->x1); y = g->y0;
    PageCoords( x, y, &px, &py, g );
    fprintf( fp, "       <text\n" );
    fprintf( fp, "          xml:space=\"preserve\"\n");
        fprintf( fp, "      style=\"font-size:18px;font-style:normal;font-variant:normal;font-weight:normal;font-stretch:normal;vertical-align:middle;text-align:left;baseline-shift:-2.0em;line-height:100%;letter-spacing:0px;word-spacing:0px;writing-mode:lr-tb;text-anchor:middle;fill:#ffffff;fill-opacity:1;stroke:none;font-family:Myriad Pro;-inkscape-font-specification:Myriad Pro\"\n");
    fprintf( fp, "          x=\"%g\"\n", px);
    fprintf( fp, "          y=\"%g\"\n", py);
    fprintf( fp, "          ><tspan \n");
    fprintf( fp, "          id=\"%04d\"\n", g->Tspan_n++);
    fprintf( fp, "          x=\"%g\"\n", px );
    fprintf( fp, "          y=\"%g\">%s</tspan>\n", py, Title);
    fprintf( fp, "       </text>\n" );

}


SvgYnumbers( FILE *fp, double y0, double y1, double yinc, double *a, char *Title, CoordInfo *g ) {

    double  x, y, px, py;
    int     i;


    // Y-Numbers
    for ( x=g->x0, y=y0; y<= y1; y += yinc ){
        i = (int)(y+0.5);
        PageCoords( x, y, &px, &py, g );
        fprintf( fp, "       <text\n" );
        fprintf( fp, "          xml:space=\"preserve\"\n");
        fprintf( fp, "          style=\"font-size:12px;font-style:normal;font-variant:normal;font-weight:normal;font-stretch:normal;text-align:center;vertical-align:middle;baseline-shift:-.25em;line-height:100%;letter-spacing:0px;word-spacing:0px;writing-mode:lr-tb;text-anchor:end;fill:#ffffff;fill-opacity:1;stroke:none;font-family:Myriad Pro;-inkscape-font-specification:Myriad Pro\"\n");
        fprintf( fp, "          x=\"%g\"\n", px-5);
        fprintf( fp, "          y=\"%g\"\n", py);
        fprintf( fp, "          ><tspan\n");
        fprintf( fp, "          sodipodi:role=\"line\"\n");
        fprintf( fp, "          id=\"%04d\"\n", g->Tspan_n++);
        fprintf( fp, "          x=\"%g\"\n", px-5 );
        fprintf( fp, "          y=\"%g\">%.4g</tspan></text>\n", py, a[i]);
    }


    y = 0.5*(g->y0+g->y1); x = 0.0;
    PageCoords( x, y, &px, &py, g );
px -= 60;
    fprintf( fp, "<text\n");
    fprintf( fp, "       x=\"%g\"\n", px);
    fprintf( fp, "       y=\"%g\"\n", py);
    fprintf( fp, "       transform = \"rotate(-90 %g %g)\"\n", px, py);
    fprintf( fp, "       id=\"text5213\"\n");
    fprintf( fp, "       xml:space=\"preserve\"\n");
    fprintf( fp, "       style=\"font-size:18px;font-style:normal;font-variant:normal;font-weight:normal;font-stretch:normal;text-align:center;line-height:125%;letter-spacing:0px;word-spacing:0px;writing-mode:lr-tb;text-anchor:middle;fill:#ffffff;fill-opacity:1;stroke:none;font-family:Myriad Pro;-inkscape-font-specification:Myriad Pro\"><tspan\n");
    fprintf( fp, "         x=\"%g\"\n", px);
    fprintf( fp, "         y=\"%g\"\n", py);
    fprintf( fp, "         id=\"%04d\">%s</tspan></text>\n", g->Tspan_n++, Title );


}

SvgBarLabels( FILE *fp, double y0, double y1, char *Title, CoordInfo *g ) {

    double  x, y, px, py, y0i, y1i;

    y0i = ceil(y0);
    y1i = floor(y1);

    // Y-ticks
    for ( x=g->x1, y=y0i; y<= y1i; y += 1.0 ){
        PageCoords( x, y, &px, &py, g );
        fprintf( fp, "       <path\n" );
        fprintf( fp, "       d=\"m %g,%g -5,0\"\n", px, py );
        fprintf( fp, "       id=\"%04d\"\n", g->Path_n++ );
        fprintf( fp, "       style=\"fill:none;stroke:#ffffff;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\" />\n" );

        fprintf( fp, "       <text\n" );
        fprintf( fp, "          xml:space=\"preserve\"\n");
        fprintf( fp, "          style=\"font-size:12px;font-style:normal;font-variant:normal;font-weight:normal;font-stretch:normal;text-align:left;vertical-align:middle;baseline-shift:-.25em;line-height:100%;letter-spacing:0px;word-spacing:0px;writing-mode:lr-tb;text-anchor:end;fill:#ffffff;fill-opacity:1;stroke:none;font-family:Myriad Pro;-inkscape-font-specification:Myriad Pro\"\n");
        fprintf( fp, "          x=\"%g\"\n", px+15);
        fprintf( fp, "          y=\"%g\"\n", py);
        fprintf( fp, "          ><tspan\n");
        fprintf( fp, "          sodipodi:role=\"line\"\n");
        fprintf( fp, "          id=\"%04d\"\n", g->Tspan_n++);
        fprintf( fp, "          x=\"%g\"\n", px+15 );
        fprintf( fp, "          y=\"%g\">%.4g</tspan></text>\n", py, y );
    }

    y = 0.5*(g->y0+g->y1); x = 0.0;
    PageCoords( x, y, &px, &py, g );
px += 40;
    fprintf( fp, "<text\n");
    fprintf( fp, "       x=\"%g\"\n", px);
    fprintf( fp, "       y=\"%g\"\n", py);
    fprintf( fp, "       transform = \"rotate(-90 %g %g)\"\n", px, py);
    fprintf( fp, "       id=\"text5213\"\n");
    fprintf( fp, "       xml:space=\"preserve\"\n");
    fprintf( fp, "       style=\"font-size:14px;font-style:normal;font-variant:normal;font-weight:normal;font-stretch:normal;text-align:center;line-height:125%;letter-spacing:0px;word-spacing:0px;writing-mode:lr-tb;text-anchor:middle;fill:#ffffff;fill-opacity:1;stroke:none;font-family:Myriad Pro;-inkscape-font-specification:Myriad Pro\"><tspan\n");
    fprintf( fp, "         x=\"%g\"\n", px);
    fprintf( fp, "         y=\"%g\"\n", py);
    fprintf( fp, "         id=\"%04d\">%s</tspan></text>\n", g->Tspan_n++, Title );

}


void DumpXtitle( FILE *fp ) {
    fprintf( fp,  "    <g\n" );
    fprintf( fp,  "       style=\"font-style:normal;font-variant:normal;font-weight:normal;font-stretch:normal;letter-spacing:normal;word-spacing:normal;text-anchor:start;fill:#ffffff;fill-opacity:1;fill-rule:evenodd;stroke:#000000;stroke-width:4.48816681;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10.43299961;stroke-opacity:1;stroke-dasharray:none;stroke-dashoffset:0\"\n" );
    fprintf( fp,  "       ns0:preamble=\"\"\n" );
    fprintf( fp,  "       ns0:text=\"$K = \\int_{s_s}^{s_n} ( B_m - B )  ds \\:\\:\\:\\:\\: ( R_e\\;G^{1/2})$\"\n" );
    fprintf( fp,  "       word-spacing=\"normal\"\n" );
    fprintf( fp,  "       letter-spacing=\"normal\"\n" );
    fprintf( fp,  "       font-size-adjust=\"none\"\n" );
    fprintf( fp,  "       font-stretch=\"normal\"\n" );
    fprintf( fp,  "       font-weight=\"normal\"\n" );
    fprintf( fp,  "       font-variant=\"normal\"\n" );
    fprintf( fp,  "       font-style=\"normal\"\n" );
    fprintf( fp,  "       stroke-miterlimit=\"10.433\"\n" );
    fprintf( fp,  "       xml:space=\"preserve\"\n" );
    fprintf( fp,  "       transform=\"matrix(1.8947658,0,0,-2.6200317,-158.60379,2848.0556)\"\n" );
    fprintf( fp,  "       id=\"content\">\n" );
    fprintf( fp,  "<path\n" );
    fprintf( fp,  "   id=\"path5032\"\n" );
    fprintf( fp,  "   d=\"m 228.51,711.15 -0.01,0 0,0 0,0.01 0,0 0,0 0,0.01 0,0 -0.01,0.01 0,0 0,0.01 0,0 0,0.01 -0.01,0.01 0,0.01 -0.01,0.01 0,0.01 0,0 0,0.01 -0.01,0 0,0 0,0.01 0,0 0,0.01 0,0 0,0.01 0,0 0,0 0,0 c 0,0.01 0.17,0.15 0.28,0.23 l 1.75,1.35 c 0.94,0.68 1.33,0.72 1.62,0.75 0.08,0.01 0.18,0.02 0.18,0.2 0,0.04 -0.03,0.11 -0.11,0.11 -0.22,0 -0.46,-0.03 -0.7,-0.03 -0.36,0 -0.75,0.03 -1.11,0.03 -0.07,0 -0.19,0 -0.19,-0.2 0,-0.07 0.05,-0.1 0.12,-0.11 0.22,-0.02 0.31,-0.07 0.31,-0.2 0,-0.18 -0.3,-0.41 -0.36,-0.46 l -3.89,-2.99 0.8,3.2 c 0.09,0.36 0.11,0.45 0.84,0.45 0.25,0 0.34,0 0.34,0.2 0,0.09 -0.08,0.11 -0.14,0.11 -0.28,0 -1,-0.03 -1.28,-0.03 -0.29,0 -1,0.03 -1.29,0.03 -0.06,0 -0.19,0 -0.19,-0.19 0,-0.12 0.08,-0.12 0.29,-0.12 0.12,0 0.3,-0.01 0.42,-0.02 0.16,-0.02 0.22,-0.05 0.22,-0.16 0,-0.03 -0.01,-0.07 -0.04,-0.18 L 225,707.9 c -0.1,-0.39 -0.12,-0.47 -0.91,-0.47 -0.17,0 -0.28,0 -0.28,-0.19 0,-0.12 0.12,-0.12 0.15,-0.12 0.28,0 0.99,0.04 1.27,0.04 0.21,0 0.42,-0.01 0.63,-0.01 0.22,0 0.44,-0.03 0.65,-0.03 0.07,0 0.2,0 0.2,0.2 0,0.11 -0.09,0.11 -0.28,0.11 -0.37,0 -0.64,0 -0.64,0.18 0,0.07 0.05,0.29 0.08,0.44 0.14,0.52 0.27,1.05 0.4,1.56 l 1.49,1.16 1.15,-2.68 c 0.12,-0.27 0.12,-0.29 0.12,-0.35 0,-0.3 -0.42,-0.31 -0.51,-0.31 -0.11,0 -0.22,0 -0.22,-0.2 0,-0.11 0.12,-0.11 0.14,-0.11 0.4,0 0.81,0.04 1.21,0.04 0.22,0 0.76,-0.04 0.98,-0.04 0.05,0 0.18,0 0.18,0.2 0,0.11 -0.11,0.11 -0.2,0.11 -0.41,0.01 -0.54,0.1 -0.69,0.45 z\"\n" );
    fprintf( fp,  "   inkscape:connector-curvature=\"0\"\n" );
    fprintf( fp,  "   style=\"fill:#ffffff;fill-opacity:1;stroke-width:0\" />\n" );
    fprintf( fp,  "<path\n" );
    fprintf( fp,  "   id=\"path5034\"\n" );
    fprintf( fp,  "   d=\"m 242.22,710.38 0.03,0 0.03,0 0.03,0 0.01,0 0.02,0 0.01,0 0.02,0.01 0.01,0 0.01,0 0.02,0 0.01,0.01 0.01,0 0.02,0 0.01,0.01 0.01,0 0.01,0.01 0.01,0.01 0.01,0.01 0.01,0 0.01,0.01 0.01,0.01 0,0.01 0.01,0 0,0.01 0,0.01 0,0 0.01,0.01 0,0.01 0,0 0,0.01 0.01,0.01 0,0 0,0.01 0,0.01 0,0.01 0,0.01 0,0.01 c 0,0.2 -0.19,0.2 -0.33,0.2 h -5.97 c -0.14,0 -0.33,0 -0.33,-0.2 0,-0.2 0.19,-0.2 0.34,-0.2 z\"\n" );
    fprintf( fp,  "   inkscape:connector-curvature=\"0\"\n" );
    fprintf( fp,  "   style=\"fill:#ffffff;fill-opacity:1;stroke-width:0\" />\n" );
    fprintf( fp,  "<path\n" );
    fprintf( fp,  "   id=\"path5036\"\n" );
    fprintf( fp,  "   d=\"m 242.23,708.45 0.03,0 0.03,0 0.02,0 0.02,0 0.01,0 0.02,0 0.01,0.01 0.01,0 0.02,0 0.01,0 0.01,0.01 0.02,0 0.01,0 0.01,0.01 0.01,0 0.01,0.01 0.01,0.01 0.01,0.01 0.01,0 0.01,0.01 0,0.01 0.01,0 0,0.01 0.01,0 0,0.01 0,0.01 0,0 0.01,0.01 0,0.01 0,0 0,0.01 0.01,0.01 0,0 0,0.01 0,0.01 0,0.01 0,0.01 0,0.01 c 0,0.2 -0.19,0.2 -0.34,0.2 h -5.95 c -0.15,0 -0.34,0 -0.34,-0.2 0,-0.2 0.19,-0.2 0.33,-0.2 z\"\n" );
    fprintf( fp,  "   inkscape:connector-curvature=\"0\"\n" );
    fprintf( fp,  "   style=\"fill:#ffffff;fill-opacity:1;stroke-width:0\" />\n" );
    fprintf( fp,  "<path\n" );
    fprintf( fp,  "   id=\"path5038\"\n" );
    fprintf( fp,  "   d=\"m 248.6,706.39 -0.01,-0.07 0,-0.08 -0.02,-0.14 -0.01,-0.14 -0.02,-0.12 -0.02,-0.13 -0.02,-0.11 -0.02,-0.12 -0.02,-0.1 -0.02,-0.1 -0.03,-0.09 -0.02,-0.09 -0.03,-0.08 -0.03,-0.08 -0.03,-0.07 -0.03,-0.07 -0.03,-0.06 -0.04,-0.06 -0.03,-0.05 -0.04,-0.05 -0.04,-0.04 -0.04,-0.04 -0.04,-0.03 -0.04,-0.03 -0.04,-0.03 -0.04,-0.02 -0.05,-0.02 -0.05,-0.02 -0.04,-0.01 -0.05,-0.01 -0.05,-0.01 -0.05,-0.01 -0.06,0 c -0.13,0 -0.42,0.03 -0.63,0.21 0.29,0.04 0.36,0.27 0.36,0.41 0,0.29 -0.22,0.42 -0.4,0.42 -0.2,0 -0.42,-0.13 -0.42,-0.43 0,-0.48 0.5,-0.83 1.09,-0.83 0.95,0 1.43,0.87 1.65,1.76 0.13,0.52 0.49,3.42 0.57,4.52 l 0.19,2.48 c 0.14,1.83 0.48,2.08 0.92,2.08 0.1,0 0.41,-0.02 0.63,-0.21 -0.29,-0.04 -0.37,-0.27 -0.37,-0.4 0,-0.29 0.22,-0.42 0.41,-0.42 0.2,0 0.42,0.13 0.42,0.43 0,0.47 -0.5,0.82 -1.1,0.82 -0.94,0 -1.33,-0.96 -1.5,-1.72 -0.12,-0.55 -0.48,-3.35 -0.57,-4.56 z\"\n" );
    fprintf( fp,  "   inkscape:connector-curvature=\"0\"\n" );
    fprintf( fp,  "   style=\"fill:#ffffff;fill-opacity:1;stroke-width:0\" />\n" );
    fprintf( fp,  "<path\n" );
    fprintf( fp,  "   id=\"path5040\"\n" );
    fprintf( fp,  "   d=\"m 255.53,715.3 -0.01,0 -0.02,-0.01 -0.02,0 -0.01,-0.01 -0.02,0 -0.01,-0.01 -0.02,-0.01 -0.01,-0.01 -0.01,-0.01 -0.02,0 -0.01,-0.01 -0.01,-0.01 -0.01,-0.01 -0.01,-0.01 -0.01,-0.01 -0.01,-0.01 -0.01,-0.02 -0.01,-0.01 -0.01,-0.01 -0.01,-0.01 0,-0.01 -0.01,-0.02 0,-0.01 -0.01,-0.01 0,-0.01 -0.01,-0.02 0,-0.01 0,-0.01 -0.01,-0.01 0,-0.02 0,-0.01 0,-0.01 c 0,-0.16 0.14,-0.23 0.23,-0.23 0.07,0 0.33,0.04 0.33,0.4 0,0.45 -0.5,0.6 -0.92,0.6 -1.09,0 -1.28,-0.79 -1.28,-1.01 0,-0.26 0.14,-0.42 0.24,-0.51 0.18,-0.14 0.31,-0.16 0.79,-0.25 0.14,-0.03 0.59,-0.11 0.59,-0.46 0,-0.12 -0.08,-0.38 -0.37,-0.55 -0.27,-0.16 -0.61,-0.16 -0.7,-0.16 -0.27,0 -0.67,0.07 -0.83,0.29 0.23,0.03 0.38,0.21 0.38,0.4 0,0.17 -0.12,0.26 -0.27,0.26 -0.2,0 -0.4,-0.16 -0.4,-0.47 0,-0.41 0.44,-0.67 1.12,-0.67 1.28,0 1.51,0.87 1.51,1.15 0,0.64 -0.7,0.76 -0.96,0.81 -0.06,0.01 -0.23,0.04 -0.27,0.06 -0.26,0.04 -0.39,0.19 -0.39,0.34 0,0.16 0.13,0.35 0.28,0.45 0.19,0.12 0.43,0.13 0.55,0.13 0.15,0 0.51,-0.02 0.66,-0.26 z\"\n" );
    fprintf( fp,  "   inkscape:connector-curvature=\"0\"\n" );
    fprintf( fp,  "   style=\"fill:#ffffff;fill-opacity:1;stroke-width:0\" />\n" );
    fprintf( fp,  "<path\n" );
    fprintf( fp,  "   id=\"path5042\"\n" );
    fprintf( fp,  "   d=\"m 257.26,712.18 0,-0.01 0,-0.02 -0.01,-0.01 0,-0.01 0,-0.02 -0.01,-0.01 0,-0.03 -0.01,-0.03 -0.01,-0.03 -0.01,-0.03 0,-0.01 0,-0.01 -0.01,-0.02 0,-0.01 0,-0.01 -0.01,-0.02 0,-0.01 0,-0.01 0,-0.01 0,-0.01 -0.01,-0.01 0,-0.01 0,-0.01 0,0 0,-0.01 0,0 0,-0.01 0,0 0,0 c 0,-0.11 0.09,-0.17 0.18,-0.17 0.07,0 0.19,0.04 0.24,0.16 0,0.03 0.06,0.26 0.09,0.39 0.17,0.67 0.17,0.68 0.18,0.7 0.07,0.13 0.43,0.83 1.03,0.83 0.25,0 0.35,-0.14 0.35,-0.36 0,-0.28 -0.21,-0.79 -0.32,-1.08 -0.03,-0.06 -0.05,-0.12 -0.05,-0.21 0,-0.25 0.25,-0.43 0.54,-0.43 0.53,0 0.81,0.64 0.81,0.77 0,0.07 -0.09,0.07 -0.11,0.07 -0.08,0 -0.09,-0.02 -0.12,-0.12 -0.09,-0.3 -0.32,-0.56 -0.56,-0.56 -0.11,0 -0.15,0.09 -0.15,0.19 0,0.11 0.02,0.16 0.07,0.27 0.08,0.2 0.31,0.76 0.31,1.03 0,0.37 -0.27,0.6 -0.74,0.6 -0.45,0 -0.76,-0.28 -0.94,-0.49 -0.02,0.4 -0.42,0.49 -0.6,0.49 -0.5,0 -0.65,-0.76 -0.65,-0.77 0,-0.07 0.08,-0.07 0.11,-0.07 0.08,0 0.09,0.03 0.11,0.09 0.08,0.29 0.19,0.58 0.4,0.58 0.17,0 0.2,-0.16 0.2,-0.27 0,-0.06 -0.04,-0.25 -0.08,-0.38 -0.03,-0.14 -0.08,-0.34 -0.1,-0.45 z\"\n" );
    fprintf( fp,  "   inkscape:connector-curvature=\"0\"\n" );
    fprintf( fp,  "   style=\"fill:#ffffff;fill-opacity:1;stroke-width:0\" />\n" );
    fprintf( fp,  "<path\n" );
    fprintf( fp,  "   id=\"path5044\"\n" );
    fprintf( fp,  "   d=\"m 253.6,706.2 -0.02,-0.01 -0.02,0 -0.01,-0.01 -0.02,0 -0.01,-0.01 -0.02,-0.01 -0.01,0 -0.02,-0.01 -0.01,-0.01 -0.01,-0.01 -0.01,-0.01 -0.02,-0.01 -0.01,-0.01 -0.01,-0.01 -0.01,-0.01 -0.01,-0.01 -0.01,-0.01 -0.01,-0.01 0,-0.01 -0.01,-0.02 -0.01,-0.01 -0.01,-0.01 0,-0.01 -0.01,-0.02 0,-0.01 -0.01,-0.01 0,-0.01 0,-0.02 0,-0.01 -0.01,-0.01 0,-0.02 0,-0.01 c 0,-0.16 0.14,-0.22 0.23,-0.22 0.07,0 0.34,0.04 0.34,0.39 0,0.46 -0.5,0.61 -0.93,0.61 -1.09,0 -1.28,-0.8 -1.28,-1.02 0,-0.25 0.15,-0.42 0.24,-0.5 0.18,-0.14 0.31,-0.17 0.79,-0.25 0.15,-0.03 0.59,-0.11 0.59,-0.46 0,-0.12 -0.07,-0.39 -0.37,-0.56 -0.27,-0.15 -0.61,-0.15 -0.69,-0.15 -0.28,0 -0.68,0.06 -0.84,0.29 0.23,0.03 0.38,0.2 0.38,0.4 0,0.17 -0.12,0.25 -0.27,0.25 -0.2,0 -0.4,-0.16 -0.4,-0.46 0,-0.42 0.44,-0.68 1.12,-0.68 1.28,0 1.52,0.88 1.52,1.15 0,0.64 -0.71,0.77 -0.96,0.81 -0.07,0.02 -0.24,0.05 -0.28,0.06 -0.26,0.05 -0.39,0.19 -0.39,0.35 0,0.16 0.13,0.34 0.28,0.44 0.19,0.12 0.44,0.14 0.55,0.14 0.15,0 0.51,-0.03 0.67,-0.26 z\"\n" );
    fprintf( fp,  "   inkscape:connector-curvature=\"0\"\n" );
    fprintf( fp,  "   style=\"fill:#ffffff;fill-opacity:1;stroke-width:0\" />\n" );
    fprintf( fp,  "<path\n" );
    fprintf( fp,  "   id=\"path5046\"\n" );
    fprintf( fp,  "   d=\"m 256.88,704.46 -0.01,0 -0.01,0 -0.01,0 0,-0.01 -0.01,0 -0.01,0 -0.01,-0.01 -0.02,-0.01 -0.01,0 -0.01,-0.01 -0.01,-0.01 -0.01,-0.01 -0.01,-0.01 -0.01,0 0,-0.01 -0.01,-0.01 -0.01,-0.01 0,-0.01 -0.01,-0.01 0,-0.01 -0.01,-0.01 0,-0.01 -0.01,-0.01 0,-0.01 0,0 0,-0.01 -0.01,-0.01 0,-0.01 0,-0.01 0,0 0,-0.01 0,-0.01 0,0 0,-0.01 0,0 c 0,-0.13 0.1,-0.18 0.19,-0.18 0.06,0 0.26,0.04 0.26,0.3 0,0.34 -0.39,0.44 -0.77,0.44 -0.84,0 -1,-0.51 -1,-0.69 0,-0.18 0.11,-0.28 0.2,-0.35 0.14,-0.1 0.29,-0.12 0.55,-0.16 0.31,-0.05 0.66,-0.11 0.66,-0.36 0,-0.01 0,-0.52 -0.91,-0.52 -0.19,0 -0.5,0.02 -0.66,0.17 0.18,0.04 0.26,0.18 0.26,0.3 0,0.14 -0.11,0.21 -0.22,0.21 -0.17,0 -0.31,-0.14 -0.31,-0.36 0,-0.49 0.77,-0.49 0.92,-0.49 1.02,0 1.19,0.59 1.19,0.79 0,0.45 -0.55,0.53 -0.81,0.58 -0.35,0.05 -0.59,0.09 -0.59,0.3 0,0 0,0.41 0.71,0.41 0.13,0 0.4,-0.01 0.53,-0.15 z\"\n" );
    fprintf( fp,  "   inkscape:connector-curvature=\"0\"\n" );
    fprintf( fp,  "   style=\"fill:#ffffff;fill-opacity:1;stroke-width:0\" />\n" );
    fprintf( fp,  "<path\n" );
    fprintf( fp,  "   id=\"path5048\"\n" );
    fprintf( fp,  "   d=\"m 264.98,704.74 -0.01,0 0,0.01 0,0 0,0 0,0.01 0,0 0,0 0,0.01 0,0 0,0 -0.01,0.01 0,0 0,0.01 -0.01,0 0,0.01 -0.01,0.01 0,0 -0.01,0.01 0,0 0,0.01 -0.01,0 0,0.01 0,0 -0.01,0 0,0.01 -0.01,0 0,0.01 -0.01,0 0,0.01 -0.01,0 0,0.01 -0.01,0.01 0,0 -0.01,0.01 -0.01,0 0,0.01 -0.01,0.01 -0.01,0.01 0,0 -0.01,0.01 -0.01,0.01 -0.01,0.01 c -1.24,1.25 -1.56,3.13 -1.56,4.65 0,1.73 0.38,3.47 1.61,4.71 0.13,0.12 0.13,0.14 0.13,0.17 0,0.07 -0.05,0.1 -0.1,0.1 -0.11,0 -1,-0.68 -1.59,-1.94 -0.51,-1.1 -0.63,-2.2 -0.63,-3.04 0,-0.77 0.11,-1.98 0.66,-3.1 0.6,-1.23 1.45,-1.87 1.56,-1.87 0.05,0 0.1,0.03 0.1,0.1 z\"\n" );
    fprintf( fp,  "   inkscape:connector-curvature=\"0\"\n" );
    fprintf( fp,  "   style=\"fill:#ffffff;fill-opacity:1;stroke-width:0\" />\n" );
    fprintf( fp,  "<path\n" );
    fprintf( fp,  "   id=\"path5050\"\n" );
    fprintf( fp,  "   d=\"m 267.14,707.9 0,-0.03 -0.01,-0.04 -0.01,-0.03 -0.01,-0.03 -0.01,-0.03 -0.01,-0.03 -0.01,-0.02 -0.01,-0.03 -0.01,-0.02 -0.01,-0.02 -0.01,-0.02 -0.02,-0.02 -0.01,-0.01 -0.02,-0.02 -0.02,-0.01 -0.02,-0.02 -0.02,-0.01 -0.03,-0.01 -0.03,-0.01 -0.03,-0.01 -0.03,-0.01 -0.02,0 -0.01,-0.01 -0.02,0 -0.02,0 -0.02,0 -0.02,-0.01 -0.02,0 -0.03,0 -0.02,0 -0.02,-0.01 -0.03,0 -0.03,0 -0.02,0 -0.03,0 -0.03,0 -0.03,0 -0.03,0 -0.03,0 -0.04,-0.01 -0.03,0 -0.04,0 -0.03,0 c -0.17,0 -0.27,0 -0.27,-0.2 0,-0.11 0.09,-0.11 0.27,-0.11 h 3.56 c 1.57,0 2.75,1.18 2.75,2.16 0,0.71 -0.58,1.29 -1.55,1.4 1.04,0.19 2.09,0.92 2.09,1.87 0,0.74 -0.66,1.37 -1.86,1.37 l -0.14,-0.31 c 0.88,0 1.09,-0.58 1.09,-1.02 0,-0.88 -0.86,-1.82 -2.07,-1.82 h -1.46 l -0.06,-0.22 h 1.88 c 0.96,0 1.15,-0.74 1.15,-1.17 0,-0.98 -0.89,-1.95 -2.06,-1.95 h -1.36 c -0.14,0 -0.16,0 -0.22,0.01 -0.1,0.01 -0.13,0.02 -0.13,0.1 0,0.03 0,0.05 0.05,0.23 l 0.69,2.78 0.06,0.22 0.62,2.47 c 0.09,0.34 0.11,0.37 0.54,0.37 h 1.28 l 0.14,0.31 h -3.35 c -0.19,0 -0.29,0 -0.29,-0.2 0,-0.11 0.09,-0.11 0.28,-0.11 0.02,0 0.21,0 0.38,-0.02 0.18,-0.02 0.27,-0.03 0.27,-0.16 0,-0.03 -0.01,-0.07 -0.04,-0.18 z\"\n" );
    fprintf( fp,  "   inkscape:connector-curvature=\"0\"\n" );
    fprintf( fp,  "   style=\"fill:#ffffff;fill-opacity:1;stroke-width:0\" />\n" );
    fprintf( fp,  "<path\n" );
    fprintf( fp,  "   id=\"path5052\"\n" );
    fprintf( fp,  "   d=\"m 276.13,706.19 0,-0.01 -0.01,-0.01 0,-0.02 -0.01,-0.01 0,-0.02 0,-0.01 -0.01,-0.03 -0.01,-0.04 -0.01,-0.03 -0.01,-0.03 0,-0.01 0,-0.02 -0.01,-0.01 0,-0.02 0,-0.01 -0.01,-0.02 0,-0.01 0,-0.01 -0.01,-0.02 0,-0.01 0,-0.01 0,-0.01 0,-0.01 -0.01,0 0,-0.01 0,-0.01 0,0 0,-0.01 c 0,-0.15 0.13,-0.22 0.24,-0.22 0.12,0 0.24,0.09 0.27,0.15 0.03,0.07 0.09,0.29 0.13,0.44 0.03,0.13 0.11,0.45 0.15,0.62 0.04,0.15 0.08,0.31 0.12,0.47 0.07,0.3 0.07,0.31 0.21,0.53 0.23,0.34 0.57,0.73 1.12,0.73 0.39,0 0.41,-0.32 0.41,-0.48 0,-0.42 -0.3,-1.19 -0.41,-1.49 -0.08,-0.19 -0.11,-0.25 -0.11,-0.37 0,-0.37 0.31,-0.6 0.66,-0.6 0.7,0 1.01,0.96 1.01,1.06 0,0.1 -0.09,0.1 -0.12,0.1 -0.09,0 -0.1,-0.05 -0.13,-0.12 -0.16,-0.56 -0.46,-0.84 -0.73,-0.84 -0.15,0 -0.18,0.09 -0.18,0.24 0,0.16 0.03,0.25 0.16,0.56 0.08,0.22 0.37,0.96 0.37,1.35 0,0.11 0,0.4 -0.26,0.6 -0.12,0.09 -0.32,0.19 -0.65,0.19 -0.62,0 -1,-0.41 -1.23,-0.71 -0.05,0.6 -0.55,0.71 -0.9,0.71 -0.58,0 -0.97,-0.36 -1.18,-0.63 -0.05,0.47 -0.46,0.63 -0.75,0.63 -0.3,0 -0.46,-0.22 -0.55,-0.38 -0.15,-0.25 -0.25,-0.65 -0.25,-0.69 0,-0.09 0.1,-0.09 0.12,-0.09 0.1,0 0.11,0.02 0.15,0.21 0.11,0.41 0.24,0.75 0.51,0.75 0.18,0 0.23,-0.15 0.23,-0.34 0,-0.13 -0.06,-0.39 -0.11,-0.57 -0.05,-0.19 -0.12,-0.48 -0.15,-0.63 l -0.22,-0.89 c -0.03,-0.09 -0.07,-0.27 -0.07,-0.29 0,-0.15 0.12,-0.22 0.23,-0.22 0.13,0 0.24,0.09 0.27,0.15 0.04,0.07 0.09,0.29 0.13,0.44 0.03,0.13 0.11,0.45 0.15,0.62 0.04,0.15 0.09,0.31 0.12,0.47 0.08,0.28 0.09,0.34 0.29,0.62 0.2,0.28 0.53,0.64 1.05,0.64 0.4,0 0.41,-0.35 0.41,-0.48 0,-0.18 -0.02,-0.27 -0.12,-0.66 z\"\n" );
    fprintf( fp,  "   inkscape:connector-curvature=\"0\"\n" );
    fprintf( fp,  "   style=\"fill:#ffffff;fill-opacity:1;stroke-width:0\" />\n" );
    fprintf( fp,  "<path\n" );
    fprintf( fp,  "   id=\"path5054\"\n" );
    fprintf( fp,  "   d=\"m 289.46,709.41 0.03,0 0.03,0 0.03,0.01 0.02,0 0.01,0 0.02,0 0.01,0 0.02,0 0.01,0 0.02,0.01 0.01,0 0.01,0 0.02,0.01 0.01,0 0.01,0.01 0.01,0.01 0.01,0 0.01,0.01 0.01,0.01 0.01,0.01 0.01,0.01 0,0 0,0.01 0.01,0 0,0.01 0,0.01 0,0 0.01,0.01 0,0.01 0,0 0,0.01 0,0.01 0,0.01 0,0.01 0.01,0 0,0.01 0,0.01 c 0,0.2 -0.18,0.2 -0.35,0.2 h -5.4 c -0.17,0 -0.35,0 -0.35,-0.2 0,-0.2 0.18,-0.2 0.35,-0.2 z\"\n" );
    fprintf( fp,  "   inkscape:connector-curvature=\"0\"\n" );
    fprintf( fp,  "   style=\"fill:#ffffff;fill-opacity:1;stroke-width:0\" />\n" );
    fprintf( fp,  "<path\n" );
    fprintf( fp,  "   id=\"path5056\"\n" );
    fprintf( fp,  "   d=\"m 294.45,707.9 -0.01,-0.03 -0.01,-0.04 -0.01,-0.03 -0.01,-0.03 -0.01,-0.03 -0.01,-0.03 -0.01,-0.02 -0.01,-0.03 -0.01,-0.02 -0.01,-0.02 -0.01,-0.02 -0.02,-0.02 -0.01,-0.01 -0.02,-0.02 -0.02,-0.01 -0.02,-0.02 -0.02,-0.01 -0.03,-0.01 -0.02,-0.01 -0.03,-0.01 -0.04,-0.01 -0.01,0 -0.02,-0.01 -0.02,0 -0.02,0 -0.02,0 -0.02,-0.01 -0.02,0 -0.03,0 -0.02,0 -0.02,-0.01 -0.03,0 -0.03,0 -0.02,0 -0.03,0 -0.03,0 -0.03,0 -0.03,0 -0.03,0 -0.04,-0.01 -0.03,0 -0.03,0 -0.04,0 c -0.17,0 -0.27,0 -0.27,-0.2 0,-0.11 0.09,-0.11 0.27,-0.11 h 3.56 c 1.57,0 2.75,1.18 2.75,2.16 0,0.71 -0.58,1.29 -1.55,1.4 1.04,0.19 2.09,0.92 2.09,1.87 0,0.74 -0.66,1.37 -1.86,1.37 l -0.14,-0.31 c 0.88,0 1.09,-0.58 1.09,-1.02 0,-0.88 -0.86,-1.82 -2.07,-1.82 h -1.46 l -0.06,-0.22 h 1.88 c 0.96,0 1.15,-0.74 1.15,-1.17 0,-0.98 -0.88,-1.95 -2.06,-1.95 h -1.36 c -0.14,0 -0.16,0 -0.22,0.01 -0.1,0.01 -0.13,0.02 -0.13,0.1 0,0.03 0,0.05 0.05,0.23 l 0.69,2.78 0.06,0.22 0.62,2.47 c 0.09,0.34 0.11,0.37 0.54,0.37 h 1.28 l 0.14,0.31 h -3.35 c -0.19,0 -0.29,0 -0.29,-0.2 0,-0.11 0.09,-0.11 0.28,-0.11 0.02,0 0.21,0 0.38,-0.02 0.18,-0.02 0.27,-0.03 0.27,-0.16 0,-0.03 -0.01,-0.07 -0.04,-0.18 z\"\n" );
    fprintf( fp,  "   inkscape:connector-curvature=\"0\"\n" );
    fprintf( fp,  "   style=\"fill:#ffffff;fill-opacity:1;stroke-width:0\" />\n" );
    fprintf( fp,  "<path\n" );
    fprintf( fp,  "   id=\"path5058\"\n" );
    fprintf( fp,  "   d=\"m 303.79,709.61 0,0.08 0,0.07 -0.01,0.08 0,0.08 -0.01,0.17 -0.01,0.18 -0.01,0.18 -0.03,0.19 -0.02,0.2 -0.04,0.19 -0.03,0.21 -0.05,0.21 -0.05,0.21 -0.07,0.21 -0.07,0.21 -0.07,0.21 -0.05,0.11 -0.04,0.11 -0.05,0.1 -0.05,0.11 c -0.6,1.22 -1.46,1.87 -1.55,1.87 -0.06,0 -0.1,-0.04 -0.1,-0.1 0,-0.03 0,-0.05 0.19,-0.23 0.97,-0.98 1.54,-2.57 1.54,-4.65 0,-1.7 -0.37,-3.45 -1.6,-4.7 -0.13,-0.12 -0.13,-0.14 -0.13,-0.17 0,-0.06 0.04,-0.1 0.1,-0.1 0.09,0 0.99,0.67 1.58,1.94 0.51,1.09 0.63,2.2 0.63,3.03 z\"\n" );
    fprintf( fp,  "   inkscape:connector-curvature=\"0\"\n" );
    fprintf( fp,  "   style=\"fill:#ffffff;fill-opacity:1;stroke-width:0\" />\n" );
    fprintf( fp,  "<path\n" );
    fprintf( fp,  "   id=\"path5060\"\n" );
    fprintf( fp,  "   d=\"m 309.92,713.92 0,0 0,0 0,0.01 0,0 0,0 0,0 0,0.01 0,0 0,0 0,0.01 0,0 0,0.01 0,0 0,0.01 -0.01,0 0,0.01 0,0 0,0.01 -0.01,0 0,0.01 -0.01,0 0,0 0,0 0,0.01 -0.01,0 0,0 0,0 -0.01,0 0,0.01 0,0 -0.01,0 0,0 0,0 -0.01,0 0,0.01 -0.01,0 0,0 -0.01,0 0,0 -0.01,0 0,0 -0.01,0 0,0 -0.01,0 0,0 c -0.15,0 -1.1,-0.09 -1.27,-0.11 -0.08,-0.01 -0.14,-0.06 -0.14,-0.19 0,-0.12 0.09,-0.12 0.24,-0.12 0.48,0 0.5,-0.07 0.5,-0.16 l -0.03,-0.2 -0.6,-2.36 c -0.18,0.37 -0.47,0.63 -0.92,0.63 -1.16,0 -2.4,-1.46 -2.4,-2.91 0,-0.94 0.55,-1.59 1.33,-1.59 l 0.01,0.21 c -0.49,0 -0.63,0.55 -0.63,0.94 0,0.5 0.32,1.72 0.55,2.18 0.3,0.59 0.75,0.95 1.15,0.95 0.65,0 0.79,-0.81 0.79,-0.87 0,-0.06 -0.02,-0.12 -0.03,-0.17 l -0.5,-1.95 c -0.05,-0.18 -0.05,-0.2 -0.2,-0.37 -0.44,-0.55 -0.85,-0.71 -1.13,-0.71 l -0.01,-0.21 c 0.19,0 0.69,0.03 1.29,0.74 0.08,-0.42 0.43,-0.74 0.91,-0.74 0.35,0 0.58,0.22 0.74,0.54 0.17,0.36 0.29,0.97 0.29,0.99 0,0.1 -0.08,0.1 -0.11,0.1 -0.1,0 -0.11,-0.04 -0.14,-0.18 -0.17,-0.65 -0.35,-1.24 -0.76,-1.24 -0.27,0 -0.3,0.26 -0.3,0.46 0,0.24 0.02,0.31 0.06,0.48 z\"\n" );
    fprintf( fp,  "   inkscape:connector-curvature=\"0\"\n" );
    fprintf( fp,  "   style=\"fill:#ffffff;fill-opacity:1;stroke-width:0\" />\n" );
    fprintf( fp,  "<path\n" );
    fprintf( fp,  "   id=\"path5062\"\n" );
    fprintf( fp,  "   d=\"m 313.86,710.85 -0.02,0 -0.03,-0.01 -0.02,0 -0.03,-0.01 -0.02,0 -0.02,-0.01 -0.03,-0.01 -0.02,-0.01 -0.02,-0.01 -0.02,-0.01 -0.02,-0.01 -0.02,-0.01 -0.02,-0.01 -0.02,-0.02 -0.01,-0.01 -0.02,-0.02 -0.02,-0.01 -0.01,-0.02 -0.01,-0.02 -0.02,-0.02 -0.01,-0.01 -0.01,-0.02 -0.01,-0.02 -0.01,-0.02 -0.01,-0.02 0,-0.02 -0.01,-0.02 0,-0.02 -0.01,-0.02 0,-0.02 0,-0.02 0,-0.02 c 0,-0.14 0.09,-0.29 0.31,-0.29 0.21,0 0.45,0.17 0.45,0.56 0,0.45 -0.42,0.85 -1.18,0.85 -1.32,0 -1.69,-1.01 -1.69,-1.45 0,-0.78 0.74,-0.92 1.03,-0.98 0.52,-0.1 1.04,-0.21 1.04,-0.76 0,-0.26 -0.23,-1.1 -1.43,-1.1 -0.14,0 -0.91,0 -1.14,0.53 0.38,-0.05 0.63,0.25 0.63,0.53 0,0.23 -0.16,0.35 -0.37,0.35 -0.26,0 -0.56,-0.21 -0.56,-0.66 0,-0.57 0.57,-0.96 1.43,-0.96 1.61,0 2,1.2 2,1.65 0,0.36 -0.18,0.61 -0.3,0.72 -0.27,0.28 -0.56,0.33 -1,0.42 -0.36,0.08 -0.76,0.15 -0.76,0.6 0,0.29 0.24,0.89 1.12,0.89 0.25,0 0.74,-0.07 0.89,-0.45 z\"\n" );
    fprintf( fp,  "   inkscape:connector-curvature=\"0\"\n" );
    fprintf( fp,  "   style=\"fill:#ffffff;fill-opacity:1;stroke-width:0\" />\n" );
    fprintf( fp,  "<path\n" );
    fprintf( fp,  "   id=\"path5064\"\n" );
    fprintf( fp,  "   d=\"m 329.01,704.74 0,0 0,0.01 0,0 -0.01,0 0,0.01 0,0 0,0 0,0.01 0,0 0,0 -0.01,0.01 0,0 0,0.01 -0.01,0 0,0.01 0,0.01 -0.01,0 -0.01,0.01 0,0 0,0.01 -0.01,0 0,0.01 0,0 -0.01,0 0,0.01 -0.01,0 0,0.01 -0.01,0 0,0.01 -0.01,0 0,0.01 -0.01,0.01 0,0 -0.01,0.01 0,0 -0.01,0.01 -0.01,0.01 0,0.01 -0.01,0 -0.01,0.01 -0.01,0.01 0,0.01 c -1.25,1.25 -1.57,3.13 -1.57,4.65 0,1.73 0.38,3.47 1.61,4.71 0.13,0.12 0.13,0.14 0.13,0.17 0,0.07 -0.04,0.1 -0.1,0.1 -0.1,0 -1,-0.68 -1.59,-1.94 -0.51,-1.1 -0.63,-2.2 -0.63,-3.04 0,-0.77 0.11,-1.98 0.66,-3.1 0.6,-1.23 1.46,-1.87 1.56,-1.87 0.06,0 0.1,0.03 0.1,0.1 z\"\n" );
    fprintf( fp,  "   inkscape:connector-curvature=\"0\"\n" );
    fprintf( fp,  "   style=\"fill:#ffffff;fill-opacity:1;stroke-width:0\" />\n" );
    fprintf( fp,  "<path\n" );
    fprintf( fp,  "   id=\"path5066\"\n" );
    fprintf( fp,  "   d=\"m 334.76,710.51 0.09,0.02 0.1,0.03 0.09,0.02 0.1,0.03 0.09,0.04 0.1,0.03 0.1,0.04 0.09,0.04 0.1,0.04 0.09,0.05 0.09,0.05 0.09,0.05 0.09,0.05 0.08,0.06 0.09,0.05 0.08,0.06 0.08,0.07 0.07,0.06 0.07,0.07 0.07,0.07 0.06,0.07 0.06,0.07 0.06,0.08 0.05,0.08 0.04,0.08 0.04,0.08 0.04,0.08 0.03,0.09 0.02,0.09 0.01,0.09 0.01,0.09 0.01,0.09 c 0,0.86 -0.9,1.49 -2.21,1.49 h -2.84 c -0.2,0 -0.29,0 -0.29,-0.2 0,-0.11 0.09,-0.11 0.28,-0.11 0.02,0 0.21,0 0.38,-0.02 0.18,-0.02 0.27,-0.03 0.27,-0.16 0,-0.03 -0.01,-0.07 -0.04,-0.18 l -1.33,-5.35 c -0.1,-0.39 -0.12,-0.47 -0.91,-0.47 -0.18,0 -0.27,0 -0.27,-0.2 0,-0.11 0.12,-0.11 0.14,-0.11 0.28,0 0.98,0.04 1.26,0.04 0.27,0 0.98,-0.04 1.26,-0.04 0.08,0 0.2,0 0.2,0.2 0,0.11 -0.09,0.11 -0.28,0.11 -0.37,0 -0.65,0 -0.65,0.18 0,0.06 0.02,0.11 0.03,0.17 l 0.66,2.64 h 1.19 l -0.01,0.22 h -1.12 l 0.65,2.6 c 0.06,0.23 0.09,0.33 0.28,0.36 0.09,0.01 0.41,0.01 0.61,0.01 0.7,0 1.81,0 1.81,-0.98 0,-0.34 -0.16,-1.03 -0.55,-1.41 -0.26,-0.26 -0.79,-0.58 -1.68,-0.58 l 0.01,-0.22 c 0.9,0 1.08,-0.56 1.08,-0.91 0,-0.15 -0.08,-0.45 -0.14,-0.68 -0.07,-0.28 -0.16,-0.65 -0.16,-0.85 0,-1.07 1.2,-1.07 1.33,-1.07 0.85,0 1.2,1 1.2,1.14 0,0.12 -0.11,0.12 -0.12,0.12 -0.09,0 -0.11,-0.07 -0.13,-0.14 -0.25,-0.74 -0.68,-0.91 -0.91,-0.91 -0.33,0 -0.4,0.22 -0.4,0.61 0,0.31 0.06,0.82 0.1,1.14 0.02,0.14 0.04,0.33 0.04,0.46 0,0.77 -0.67,1.08 -0.93,1.18 z\"\n" );
    fprintf( fp,  "   inkscape:connector-curvature=\"0\"\n" );
    fprintf( fp,  "   style=\"fill:#ffffff;fill-opacity:1;stroke-width:0\" />\n" );
    fprintf( fp,  "<path\n" );
    fprintf( fp,  "   id=\"path5068\"\n" );
    fprintf( fp,  "   d=\"m 338.7,707.22 0.02,0 0.03,0 0.02,0 0.02,0 0.03,0 0.03,0 0.03,0 0.03,0 0.03,0 0.04,0 0.03,0.01 0.04,0 0.03,0 0.04,0 0.08,0.01 0.08,0.01 0.08,0.01 0.08,0.01 0.08,0.01 0.08,0.02 0.04,0.01 0.04,0.01 0.04,0.01 0.03,0.01 0.04,0.01 0.04,0.01 c 0.51,0.18 0.61,0.52 0.61,0.73 0,0.38 -0.38,0.61 -0.86,0.61 v -0.2 c 0.3,0 0.56,-0.15 0.56,-0.41 0,-0.68 -1.18,-0.68 -1.48,-0.68 h -0.34 c 0.29,1.01 1.06,1.09 1.26,1.09 v 0.2 c -0.84,0 -1.97,-0.64 -1.97,-1.83 0,-0.7 0.44,-1.31 1.25,-1.31 1.18,0 1.73,0.69 1.73,0.79 0,0.04 -0.06,0.12 -0.12,0.12 -0.04,0 -0.05,-0.02 -0.11,-0.07 -0.54,-0.64 -1.35,-0.64 -1.48,-0.64 -0.42,0 -0.7,0.27 -0.7,0.84 0,0.1 0,0.23 0.09,0.62 z\"\n" );
    fprintf( fp,  "   inkscape:connector-curvature=\"0\"\n" );
    fprintf( fp,  "   style=\"fill:#ffffff;fill-opacity:1;stroke-width:0\" />\n" );
    fprintf( fp,  "<path\n" );
    fprintf( fp,  "   id=\"path5070\"\n" );
    fprintf( fp,  "   d=\"m 351.77,714.04 0,0 0,0.01 0,0 0,0 0,0.01 -0.01,0 0,0.01 0,0 0,0 0,0.01 0,0 0,0.01 -0.01,0 0,0.01 0,0 0,0 -0.01,0.01 0,0 0,0 -0.01,0.01 0,0 -0.01,0 0,0.01 -0.01,0 0,0 -0.01,0 0,0.01 0,0 -0.01,0 0,0 0,0 -0.01,0 0,0 0,0 -0.01,0 0,0 0,0 -0.01,0 0,0 c -0.03,0 -0.04,-0.01 -0.15,-0.12 l -0.7,-0.76 c -0.09,0.14 -0.55,0.88 -1.66,0.88 -2.22,0 -4.46,-2.2 -4.46,-4.51 0,-1.58 1.11,-2.72 2.72,-2.72 0.44,0 0.89,0.09 1.25,0.24 0.49,0.19 0.69,0.4 0.86,0.6 0.09,-0.25 0.35,-0.61 0.45,-0.61 0.05,0 0.07,0.02 0.07,0.03 0.02,0.03 0.12,0.4 0.17,0.61 l 0.19,0.77 c 0.04,0.17 0.09,0.34 0.13,0.51 0.11,0.44 0.12,0.46 0.69,0.47 0.05,0 0.16,0.01 0.16,0.2 0,0.07 -0.05,0.11 -0.13,0.11 -0.23,0 -0.82,-0.03 -1.05,-0.03 -0.31,0 -1.09,0.03 -1.39,0.03 -0.09,0 -0.21,0 -0.21,-0.2 0,-0.11 0.08,-0.11 0.3,-0.11 0.01,0 0.3,0 0.52,-0.02 0.26,-0.03 0.31,-0.06 0.31,-0.18 0,-0.1 -0.11,-0.54 -0.21,-0.91 -0.28,-1.1 -1.57,-1.2 -1.92,-1.2 -0.96,0 -2,0.56 -2,2.08 0,0.31 0.1,1.96 1.14,3.25 0.54,0.68 1.51,1.28 2.49,1.28 1.02,0 1.61,-0.76 1.61,-1.92 0,-0.4 -0.03,-0.41 -0.03,-0.5 0,-0.11 0.11,-0.11 0.15,-0.11 0.13,0 0.13,0.03 0.18,0.2 z\"\n" );
    fprintf( fp,  "   inkscape:connector-curvature=\"0\"\n" );
    fprintf( fp,  "   style=\"fill:#ffffff;fill-opacity:1;stroke-width:0\" />\n" );
    fprintf( fp,  "<path\n" );
    fprintf( fp,  "   id=\"path5072\"\n" );
    fprintf( fp,  "   d=\"m 354.36,715.16 0,0.01 0,0.01 0,0.01 0,0.01 0,0 0,0.01 0,0.01 0,0.01 0,0 0,0.01 -0.01,0 0,0.01 0,0.01 0,0 0,0.01 0,0 0,0.01 0,0 0,0.01 0,0.01 -0.01,0.01 0,0 -0.01,0.01 0,0.01 0,0 -0.01,0.01 -0.01,0 0,0 -0.01,0.01 -0.01,0 -0.01,0 -0.01,0 0,0 -0.01,0 0,0 -0.01,0.01 0,0 -0.01,0 -0.01,0 0,0 -0.01,0 -0.01,0 0,0 -0.01,0 -0.01,0 -0.01,0 -0.01,0 0,0 -0.01,0 -0.01,0 c -0.45,-0.44 -1.08,-0.45 -1.37,-0.45 v -0.25 c 0.17,0 0.63,0 1.01,0.2 v -3.55 c 0,-0.23 0,-0.32 -0.69,-0.32 h -0.27 v -0.25 c 0.13,0.01 0.98,0.03 1.24,0.03 0.22,0 1.1,-0.02 1.25,-0.03 v 0.25 h -0.27 c -0.69,0 -0.69,0.09 -0.69,0.32 z\"\n" );
    fprintf( fp,  "   inkscape:connector-curvature=\"0\"\n" );
    fprintf( fp,  "   style=\"fill:#ffffff;fill-opacity:1;stroke-width:0\" />\n" );
    fprintf( fp,  "<path\n" );
    fprintf( fp,  "   id=\"path5074\"\n" );
    fprintf( fp,  "   d=\"m 359.48,715.65 0.01,0 0,0.01 0,0 0,0.01 0,0 0,0 0.01,0.01 0,0 0,0 0,0.01 0,0 0,0 0.01,0.01 0,0.01 0,0 0,0.01 0,0 0.01,0 0,0.01 0,0 0,0.01 0,0 0,0 0,0.01 0,0 0,0 0,0 0,0.01 0.01,0 0,0 0,0.01 0,0 0,0.01 0,0 c 0,0.11 -0.1,0.18 -0.17,0.18 -0.12,0 -0.15,-0.09 -0.19,-0.18 l -2.6,-6.47 c -0.04,-0.09 -0.04,-0.11 -0.04,-0.13 0,-0.11 0.09,-0.18 0.17,-0.18 0.12,0 0.15,0.09 0.19,0.18 z\"\n" );
    fprintf( fp,  "   inkscape:connector-curvature=\"0\"\n" );
    fprintf( fp,  "   style=\"fill:#ffffff;fill-opacity:1;stroke-width:0\" />\n" );
    fprintf( fp,  "<path\n" );
    fprintf( fp,  "   id=\"path5076\"\n" );
    fprintf( fp,  "   d=\"m 363.59,712.01 h -0.23 c -0.02,-0.16 -0.09,-0.57 -0.18,-0.64 -0.06,-0.04 -0.6,-0.04 -0.69,-0.04 h -1.29 c 0.74,0.65 0.98,0.84 1.4,1.17 0.51,0.41 0.99,0.84 0.99,1.5 0,0.84 -0.73,1.36 -1.63,1.36 -0.86,0 -1.44,-0.61 -1.44,-1.25 0,-0.35 0.3,-0.39 0.37,-0.39 0.16,0 0.36,0.12 0.36,0.37 0,0.13 -0.05,0.37 -0.41,0.37 0.22,0.5 0.69,0.65 1.02,0.65 0.7,0 1.06,-0.54 1.06,-1.11 0,-0.6 -0.43,-1.08 -0.66,-1.33 l -1.68,-1.66 c -0.06,-0.06 -0.06,-0.07 -0.06,-0.27 h 2.87 z\"\n" );
    fprintf( fp,  "   inkscape:connector-curvature=\"0\"\n" );
    fprintf( fp,  "   style=\"fill:#ffffff;fill-opacity:1;stroke-width:0\" />\n" );
    fprintf( fp,  "<path\n" );
    fprintf( fp,  "   id=\"path5078\"\n" );
    fprintf( fp,  "   d=\"m 367.43,709.61 0,0.08 0,0.07 -0.01,0.08 0,0.08 -0.01,0.17 -0.01,0.18 -0.02,0.18 -0.02,0.19 -0.02,0.2 -0.04,0.19 -0.04,0.21 -0.04,0.21 -0.06,0.21 -0.06,0.21 -0.07,0.21 -0.08,0.21 -0.04,0.11 -0.04,0.11 -0.05,0.1 -0.05,0.11 c -0.6,1.22 -1.46,1.87 -1.56,1.87 -0.06,0 -0.1,-0.04 -0.1,-0.1 0,-0.03 0,-0.05 0.19,-0.23 0.98,-0.98 1.55,-2.57 1.55,-4.65 0,-1.7 -0.37,-3.45 -1.61,-4.7 -0.13,-0.12 -0.13,-0.14 -0.13,-0.17 0,-0.06 0.04,-0.1 0.1,-0.1 0.1,0 1,0.67 1.59,1.94 0.51,1.09 0.63,2.2 0.63,3.03 z\"\n" );
    fprintf( fp,  "   inkscape:connector-curvature=\"0\"\n" );
    fprintf( fp,  "   style=\"fill:#ffffff;fill-opacity:1;stroke-width:0\" />\n" );
    fprintf( fp,  "</g>\n" );
}
void DumpSvg( int n, char *Time ){

    double      a, b, r, S[36], K[36], Mu[36];
    double      x, y, px, py, Min, Max;
    int         i, j, nK, nMu;
    int         tspan_n = 0;
    int         path_n  = 0;
    char        Line[512], Filename[1024], FilenameBase[1024];
    FILE        *fp, *fp_info;
    CoordInfo   *g = Create_CoordInfo();

    nK = 36;
    a = 0.01; b = 10.0;
    r = pow( b/a, 1.0/((double)(nK-1)));
printf("r = %g\n", r);
    S[0] = a;
    K[0] = a;
//    for (j=1; j<nK; j++) S[j] = S[j-1]*r;
//    for (j=0; j<nK; j++) K[j] = S[nK-1-j];
    for (j=1; j<nK; j++) K[j] = K[j-1]*r;

    nMu = 36;
    a = 1.0; b = 2000.0;
    r = pow( b/a, 1.0/((double)(nMu-1)));
    S[0] = a;
    for (j=1; j<nMu; j++) S[j] = S[j-1]*r;
    for (j=0; j<nMu; j++) Mu[j] = S[j];


    double E[10], A[18];
    E[0] = 0.050000; E[1] = 0.075000; E[2] = 0.10500; E[3] = 0.15000; E[4] = 0.22500; E[5] = 0.31500; E[6] = 0.50000; E[7] = 0.75000; E[8] = 1.1000; E[9] = 1.5000;
    for (i=0; i<18; i++) A[i] = 5.0+i*5.0;




    fp = fopen("mike.svg", "w");

    fprintf( fp, "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n"
                 "<!-- Created with Inkscape (http://www.inkscape.org/) -->\n"
                 "\n"
                 "<svg\n"
                 "    xmlns:ns0=\"http://www.iki.fi/pav/software/textext/\""
                 "    xmlns:dc=\"http://purl.org/dc/elements/1.1/\"\n"
                 "    xmlns:cc=\"http://creativecommons.org/ns#\"\n"
                 "    xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\"\n"
                 "    xmlns:svg=\"http://www.w3.org/2000/svg\"\n"
                 "    xmlns=\"http://www.w3.org/2000/svg\"\n"
                 "    xmlns:sodipodi=\"http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd\""
                 "    xmlns:inkscape=\"http://www.inkscape.org/namespaces/inkscape\""
                 "    xmlns:xlink=\"http://www.w3.org/1999/xlink\"\n"
                 "    id=\"svg2\">\n"
                 "    width=\"744.09448\"\n"
                 "    height=\"1052.3622\"\n"
                 "    version=\"1.1\"\n"
                 "    inkscape:version=\"0.48.1 r9760\"\n"
                 "    sodipodi:docname=\"New document 1\">\n"
                 "\n"
                 "    <defs id=\"defs4\" />\n"
                 "\n"
                 "    <metadata id=\"metadata7\">\n"
                 "        <rdf:RDF>\n"
                 "            <cc:Work rdf:about=\"\">\n"
                 "                <dc:format>image/svg+xml</dc:format>\n"
                 "                <dc:type rdf:resource=\"http://purl.org/dc/dcmitype/StillImage\" />\n"
                 "                <dc:title></dc:title>\n"
                 "            </cc:Work>\n"
                 "        </rdf:RDF>\n"
                 "    </metadata>\n"
                 "\n");



    // start graphics (layer1)
    fprintf( fp, "    <g id=\"layer1\">\n" );

    // Page Background (black filled rectangle)
    fprintf( fp, "        <rect\n"
                 "           width=\"744.09448\"\n"
                 "           height=\"1052.1366\"\n"
                 "           x=\"0\"\n"
                 "           y=\"0\"\n"
                 "           id=\"Background\"\n"
                 "           style=\"fill:#000000;fill-opacity:1;stroke:none\"\n"
                 "        />\n" );


    /*
     *  Title
     */
    SvgTitle( fp, 32, 744.09448/2.0, 50.0, "Flux to Phase Space Density Conversion", g ); // Plot title
    sprintf( Line, "LANL-97A:  %s", Time );
    SvgTitle( fp, 28, 744.09448/2.0, 100.0, Line, g ); // Plot title








    /********************
     * Original FLUX 
     ********************/
    strcpy( FilenameBase, "/home/mgh/git/LanlGeoMag/Examples/FluxToPSD/Lgm_FluxToPsd_FLUX_EA" );
    sprintf( Filename, "file://%s_%03d.gif", FilenameBase, n );
    SvgImage( fp, 90, 200, 200, 200, Filename, g );
    SetupCoords( -0.5, 17.5, g->ImageX, g->ImageX+g->ImageWidth, -0.5, 9.5, g->ImageY+g->ImageHeight, g->ImageY, g );
    SvgXticks( fp, 0.0, 17.0, 2.0, g ); // X-ticks
    SvgYticks( fp, 0.0, 9.0, 1.0, g ); // Y-ticks
    SvgXnumbers( fp, 0.0, 17.0, 2.0, A, "\u03b1, Degrees", g ); // X-Numbers
    SvgYnumbers( fp, 0.0, 9.0, 1.0, E, "Energy, MeV", g );  // Y-numbers

    // Color Bar
    sprintf( Filename, "%s_%03d.info", FilenameBase, n );
    fp_info = fopen( Filename, "r" );
    fscanf( fp_info, "%*[^:]:%lf", &Min );
    fscanf( fp_info, "%*[^:]:%lf", &Max );
    fclose( fp_info );
    sprintf( Filename, "file://%s_Bar.gif", FilenameBase );
    SvgImage( fp, 90+200+10, 200, 10, 200, Filename, g );
    SetupCoords( -0.5, 9.5, g->ImageX, g->ImageX+g->ImageWidth, Min, Max, g->ImageY+g->ImageHeight, g->ImageY, g );
    SvgBarLabels( fp, Min, Max, "#/cm\u00b2/s/sr/MeV", g );




    /********************
     * Computed PSD 
     ********************/
    strcpy( FilenameBase, "/home/mgh/git/LanlGeoMag/Examples/FluxToPSD/Lgm_FluxToPsd_PSD_MK" );
    sprintf( Filename, "file://%s_%03d.gif", FilenameBase, n );
    SvgImage( fp, 460, 200, 200, 200, Filename, g );
    SetupCoords( -0.5, 35.5, g->ImageX, g->ImageX+g->ImageWidth, -0.5, 35.5, g->ImageY+g->ImageHeight, g->ImageY, g );
    SvgXticks( fp, 0.0, 35.0, 5.0, g ); // X-ticks
    SvgYticks( fp, 0.0, 35.0, 5.0, g ); // Y-ticks
    SvgXnumbers( fp, 0.0, 35.0, 5.0, K, "K, Re \u221AG", g );  // X-numbers
    SvgYnumbers( fp, 0.0, 35.0, 5.0, Mu, "\u00B5, MeV/G", g );  // Y-numbers

    
    // Color Bar
    sprintf( Filename, "%s_%03d.info", FilenameBase, n );
printf("Filename = %s\n", Filename);
    fp_info = fopen( Filename, "r" );
    fscanf( fp_info, "%*[^:]:%lf", &Min );
    fscanf( fp_info, "%*[^:]:%lf", &Max );
    fclose( fp_info );
    sprintf( Filename, "file://%s_Bar.gif", FilenameBase );
    SvgImage( fp, 460+200+10, 200, 10, 200, Filename, g );
    SetupCoords( -0.5, 9.5, g->ImageX, g->ImageX+g->ImageWidth, Min, Max, g->ImageY+g->ImageHeight, g->ImageY, g );
    SvgBarLabels( fp, Min, Max, "(c/cm/MeV)\u00b3", g );



    /********************
     * Re-derived  FLUX 
     ********************/
    strcpy( FilenameBase, "/home/mgh/git/LanlGeoMag/Examples/FluxToPSD/Lgm_PsdToFlux_FLUX_EA" );
    sprintf( Filename, "file://%s_%03d.gif", FilenameBase, n );
    SvgImage( fp, 90, 550, 200, 200, Filename, g );
    SetupCoords( -0.5, 17.5, g->ImageX, g->ImageX+g->ImageWidth, -0.5, 9.5, g->ImageY+g->ImageHeight, g->ImageY, g );
    SvgXticks( fp, 0.0, 17.0, 2.0, g ); // X-ticks
    SvgYticks( fp, 0.0, 9.0, 1.0, g ); // Y-ticks
    SvgXnumbers( fp, 0.0, 17.0, 2.0, A, "\u03b1, Degrees", g ); // X-Numbers
    SvgYnumbers( fp, 0.0, 9.0, 1.0, E, "Energy, MeV", g );  // Y-numbers

    // Color Bar
    sprintf( Filename, "%s_%03d.info", FilenameBase, n );
    fp_info = fopen( Filename, "r" );
    fscanf( fp_info, "%*[^:]:%lf", &Min );
    fscanf( fp_info, "%*[^:]:%lf", &Max );
    fclose( fp_info );
    sprintf( Filename, "file://%s_Bar.gif", FilenameBase );
    SvgImage( fp, 90+200+10, 550, 10, 200, Filename, g );
    SetupCoords( -0.5, 9.5, g->ImageX, g->ImageX+g->ImageWidth, Min, Max, g->ImageY+g->ImageHeight, g->ImageY, g );
    SvgBarLabels( fp, Min, Max, "#/cm\u00b2/s/sr/MeV", g );



    /********************
     * Difference FLUX 
     ********************/
    strcpy( FilenameBase, "/home/mgh/git/LanlGeoMag/Examples/FluxToPSD/J_DIFF" );
    sprintf( Filename, "file://%s_%03d.gif", FilenameBase, n );
    SvgImage( fp, 460, 550, 200, 200, Filename, g );
    SetupCoords( -0.5, 17.5, g->ImageX, g->ImageX+g->ImageWidth, -0.5, 9.5, g->ImageY+g->ImageHeight, g->ImageY, g );
    SvgXticks( fp, 0.0, 17.0, 2.0, g ); // X-ticks
    SvgYticks( fp, 0.0, 9.0, 1.0, g ); // Y-ticks
    SvgXnumbers( fp, 0.0, 17.0, 2.0, A, "\u03b1, Degrees", g ); // X-Numbers
    SvgYnumbers( fp, 0.0, 9.0, 1.0, E, "Energy, MeV", g );  // Y-numbers

    // Color Bar
    sprintf( Filename, "%s_%03d.info", FilenameBase, n );
    fp_info = fopen( Filename, "r" );
    fscanf( fp_info, "%*[^:]:%lf", &Min );
    fscanf( fp_info, "%*[^:]:%lf", &Max );
    fclose( fp_info );
    sprintf( Filename, "file://%s_Bar.gif", FilenameBase );
    SvgImage( fp, 460+200+10, 550, 10, 200, Filename, g );
    SetupCoords( -0.5, 9.5, g->ImageX, g->ImageX+g->ImageWidth, Min, Max, g->ImageY+g->ImageHeight, g->ImageY, g );
    SvgBarLabels( fp, Min, Max, "#/cm\u00b2/s/sr/MeV", g );








    // end graphics (layer 1)
    fprintf( fp, "    </g>\n" );

    fprintf( fp, "</svg>\n" );




    fclose(fp);

}

/*
main(){
    char    Command[4096];
    DumpSvg( 0 );
    sprintf( Command, "inkscape --export-ps=mike.ps mike.svg" ); system( Command );
    
}
*/
