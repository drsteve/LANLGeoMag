
import ctypes




class Lgm_Vector(ctypes.Structure):
    _fields_ = [ ( "x", ctypes.c_double ),
        ("y", ctypes.c_double),
        ("z", ctypes.c_double) ]

#set up types
LgmLong = ctypes.c_long
LgmBoolean = ctypes.c_int
ConstLgmBoolean = LgmBoolean
LgmFALSE = 0
LgmTRUE = 1
LgmChar = ctypes.c_char
ConstLgmChar = LgmChar
LgmCharP = ctypes.c_char_p
ConstLgmCharP = ctypes.c_char_p
LgmDouble = ctypes.c_double
ConstLgmDouble = LgmDouble
LgmDoubleP = ctypes.POINTER(LgmDouble)
c_types = [ctypes.c_byte, ctypes.c_short, ctypes.c_int,
           ctypes.c_long, ctypes.c_longlong]
c_sizes = [ctypes.sizeof(i) for i in c_types]
idx = c_sizes.index(ctypes.sizeof(LgmDouble) / 2)
LgmInt = c_types[idx]
ConstLgmInt = LgmInt
LgmIntP = ctypes.POINTER(LgmInt)
Lgm_Vector_p = ctypes.POINTER(Lgm_Vector)



#  TODO still working on this part
#
#class Lgm_MagModelInfo(ctypes.Structure):
#    _fields_ = [ \
#                ('c', Lgm_CTrans), # This contains all time info and a bunch more stuff
#                ('nFunc', LgmLong),
#                ('int                 (*Bfield)();
#    int                 SavePoints;
#    double              Hmax;
#    FILE            *fp;
#    double      W[6];
#    double      G1, G2;
#    int                 Kp;
#    double      Dst;
#    double      P;
#    double      Bx, By, Bz;
#    double      T96MOD_V[11];       /* free params for T96_MOD */
#
#    double      Trace_s;
#
#    double      B0, B1, B2;             // params for simplified mead field
#
#
#    /*
#     * Variable to control which internal model to use
#     */
#    int         InternalModel;          // Can be LGM_CDIP, LGM_EDIP or LGM_IGRF
#
#
#
#    /*
#     * Stuff that may be neeeded for field integrals, etc.
#     */
#    double              KineticEnergy;  // Particle kinetic energy
#    double              Mass;           // Particle mass
#    double              PitchAngle;     // Particle Pitch Angle
#
#    Lgm_Vector          Pm_South;
#    Lgm_Vector          Pm_North;
#    double              Bm, Sm_South, Sm_North, Blocal;
#    int                 FirstCall;
#//    double              epsabs, epsrel;
#
#    /*
#     * Arrays containing FL vals
#     */
#    double              s[LGM_MAX_INTERP_PNTS];    // distance along FL
#    double              Px[LGM_MAX_INTERP_PNTS];   // Px along FL  (in GSM)
#    double              Py[LGM_MAX_INTERP_PNTS];   // Py along FL  (in GSM)
#    double              Pz[LGM_MAX_INTERP_PNTS];   // Pz along FL  (in GSM)
#    Lgm_Vector          Bvec[LGM_MAX_INTERP_PNTS]; // 3D B-field vector   (in GSM)
#    double              Bmag[LGM_MAX_INTERP_PNTS]; // magnitude of B
#    double              BminusBcdip[LGM_MAX_INTERP_PNTS]; // magnitude of B minus magnitude of Cent. Dipole
#    int                 nPnts;      // actual number of points defined
#    double              ds;         // spacing in s (dist. along FL)
#                                    // this will help in seacrhing the
#                                    // arrays (e.g. for interpolation).
#
#    Lgm_Vector      P_gsm;          //< S/C position in GSM
#    double          S;              //< Distance along FL from southern footpoint to S/C location in Re.
#    double          B;              //< Local (model) B-field magnitude (i.e. at S/C position)
#
#    Lgm_Vector      Pmin;           //< position of minimum |B| in GSM
#    Lgm_Vector      Bvecmin;        //< value of Bvecmin
#    double          Bmin;           //< Value of |Bmin|
#    double          Smin;           //< Distance from southern footpoint to Pmin along FL.
#
#    Lgm_Vector      Spherical_Footprint_Pn;   //< position of northern footpoint (at 120km)
#    double          Spherical_Footprint_Sn;   //< Distance along FL from southern foorpoint in Re
#    double          Spherical_Footprint_Bn;   //< Value of |B| at Footprint_Pn
#
#    Lgm_Vector      Spherical_Footprint_Ps;   //< position of southern footpoint (at 120km)
#    double          Spherical_Footprint_Ss;   //< Distance along FL from southern foorpoint in Re (i.e. this one is zero by definition)
#    double          Spherical_Footprint_Bs;   //< Value of |B| at Footprint_Ps
#
#    Lgm_Vector      Ellipsoid_Footprint_Pn;   //< position of northern footpoint (at 120km above surface of WGS84 ellipsoid)
#    double          Ellipsoid_Footprint_Sn;   //< Distance along FL from southern footpoint in Re
#    double          Ellipsoid_Footprint_Bn;   //< Value of |B| at Footprint_Pn
#
#    Lgm_Vector      Ellipsoid_Footprint_Ps;   //< position of southern footpoint (at 120km)
#    double          Ellipsoid_Footprint_Ss;   //< Distance along FL from southern foorpoint in Re (i.e. this one is zero by definition)
#    double          Ellipsoid_Footprint_Bs;   //< Value of |B| at Footprint_Ps
#
#    int             FieldLineType;  //< Field line type. (I.e., LGM_OPEN_IMF, LGM_CLOSED, LGM_OPEN_N_LOBE, LGM_OPEN_S_LOBE, LGM_INSIDE_EARTH, LGM_TARGET_HEIGHT_UNREACHABLE)
#
#
#
#
#    double              d2B_ds2;    // second derivative of B wrt to s at smin.
#    double              Sb0;        // value of Sb integral for eq. mirroring particles.
#    int                 imin1;      // imin1 and imin2 are the indices in the
#    int                 imin2;      //   array between which smin is located.
#
#
#
#    /*
#     *  GSL defs for spline interpolation
#     */
#    gsl_interp_accel    *acc;       // accelerator
#    gsl_interp_accel    *accPx;     // accelerator
#    gsl_interp_accel    *accPy;     // accelerator
#    gsl_interp_accel    *accPz;     // accelerator
#    gsl_spline          *spline;    // spline object
#    gsl_spline          *splinePx;  // spline object
#    gsl_spline          *splinePy;  // spline object
#    gsl_spline          *splinePz;  // spline object
#
#
#
#    /*
#     *  Other stuff
#     */
#    int                 VerbosityLevel;
#    int                 UseInterpRoutines; // whether to use fast I and Sb routines.
#
#
#    /*
#     *  These variables are needed to make Lgm_MagStep() reentrant/thread-safe.
#     *  They basically used to be static declarations within Lgm_MagStep()
#     */
#    double  Lgm_MagStep_eps_old;
#    int     Lgm_MagStep_FirstTimeThrough;
#    int     Lgm_MagStep_kmax;
#    int     Lgm_MagStep_kopt;
#    double  Lgm_MagStep_snew;
#    double  Lgm_MagStep_A[LGM_MAGSTEP_JMAX+1];
#    double  Lgm_MagStep_alpha[LGM_MAGSTEP_IMAX+1][LGM_MAGSTEP_IMAX+1];
#    double  Lgm_MagStep_d[LGM_MAGSTEP_JMAX][LGM_MAGSTEP_JMAX];
#    double  Lgm_MagStep_x[LGM_MAGSTEP_JMAX];
#
#
#    /*
#     *  These variables are needed to make I_integrand() reentrant/thread-safe.
#     *  They basically used to be static declarations.
#     */
#    int         Lgm_I_integrand_FirstCall;
#    int         Lgm_I_integrand_JumpMethod;
#    double      Lgm_I_integrand_S;
#    Lgm_Vector  Lgm_I_integrand_P;
#    Lgm_Vector  Lgm_I_integrand_u_scale;
#    int         Lgm_n_I_integrand_Calls;
#    int         Lgm_I_Integrator;
#    double      Lgm_I_Integrator_epsrel;
#    double      Lgm_I_Integrator_epsabs;
#
#    /*
#     *  These variables are needed to make Sb_integrand() reentrant/thread-safe.
#     *  They basically used to be static declarations.
#     */
#    int         Lgm_Sb_integrand_FirstCall;
#    double      Lgm_Sb_integrand_S;
#    Lgm_Vector  Lgm_Sb_integrand_P;
#    Lgm_Vector  Lgm_Sb_integrand_u_scale;
#    int         Lgm_n_Sb_integrand_Calls;
#    int         Lgm_Sb_Integrator;              // This variable is not currently used (for now I force
#                                                // Sb to use DQAGP. )
#    double      Lgm_Sb_Integrator_epsrel;
#    double      Lgm_Sb_Integrator_epsabs;
#
#
#    /*
#     * Variables to control MagFlux integration
#     */
#    int         Lgm_MagFlux_Integrator;         // This variable is not currently used (for now I force
#                                                // MagFlux to use DQAGS. )
#    double      Lgm_MagFlux_Integrator_epsrel;
#    double      Lgm_MagFlux_Integrator_epsabs;
#
#
#    /*
#     * Variables to control LambdaIntegral integration
#     */
#    int         Lgm_LambdaIntegral_Integrator;  // This variable is not currently used (for now I force
#                                                // LambdaIntegral to use DQAGS. )
#    double      Lgm_LambdaIntegral_Integrator_epsrel;
#    double      Lgm_LambdaIntegral_Integrator_epsabs;
#
#
#    /*
#     * Some other tolerances
#     */
#    double      Lgm_FindBmRadius_Tol;
#    double      Lgm_FindShellLine_I_Tol;
#    double      Lgm_TraceToMirrorPoint_Tol;
#
#
#
#    /*
#     * Variables for defining Octree stuff
#     */
#    Lgm_OctreeCell  *OctreeRoot;
#    int             Octree_kNN_k;
#    int             Octree_kNN_InterpMethod;  // 0 = linear Div Feee; 1 = Quadratic Div Free; 2 = Newton 4-point
#    double          Octree_kNN_MaxDist;       // in physical units
#    double          OctreeScaleMin;
#    double          OctreeScaleMax;
#    double          OctreeScaleDiff;
#
#
#    /*
#     * limits for tracing. If you go outside this box, we consider it open
#     */
#    double OpenLimit_xMin;
#    double OpenLimit_xMax;
#    double OpenLimit_yMin;
#    double OpenLimit_yMax;
#    double OpenLimit_zMin;
#    double OpenLimit_zMax;
#
#    double  Lgm_LossConeHeight;
#
#
#    /*
#     *  Globals for OP77 Model
#     */
#    double OP77_TILTL;
#    double OP77_A[65], OP77_B[65], OP77_C[45], OP77_D[45], OP77_E[65], OP77_F[65], OP77_TT[5];
#
#
#
#
#} Lgm_MagModelInfo;


