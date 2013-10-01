/*
 * See 2010 CODATA reccomendations. See; Mohr P.J., et al.,
 * "CODATA recommended values of the fundamental physical constants: 2010",
 * Rev. Mod. Phys. 84, 2012, doi: 10.1103/RevModPhys.84.1527
 *
 */
#ifndef LGM_CONSTANTS_H
#define LGM_CONSTANTS_H
#include <math.h>

#ifndef LGM_c
#define LGM_c               (2.99792458e8)          // Speed of light  m/s
#endif

#ifndef LGM_e
#define LGM_e               (1.602176565e-19)       // Elementary charge, Coulombs (or s·A)
#endif


#ifndef LGM_mu0
#define LGM_mu0             (12.5663706144e-7)      // N · A^-2 (defined to be exactly 4*pi e-7)
#endif


#ifndef LGM_Eps0
#define LGM_Eps0            (8.854187817e-12)       // Permittivity of free space, A^2·s^4·kg^−1·m^−3
#endif

#ifndef LGM_U
#define LGM_U               (1.660538921e-27)       // Atomic mass unit, kg
#endif

#ifndef LGM_ELECTRON_MASS
#define LGM_ELECTRON_MASS   (9.10938291e-31)        // Rest mass of electron, kg
#endif

#ifndef LGM_PROTON_MASS
#define LGM_PROTON_MASS     (1.672621777e-27)       // Proton rest mass,  kg
#endif

#ifndef LGM_ALPHA
#define LGM_ALPHA           (7.2973525698e-3)       // Fine-structure constant
#endif

#ifndef LGM_ALPHA_INV
#define LGM_ALPHA_INV       (137.035999074)         // Inverse of Fine-structure constant
#endif

#ifndef LGM_ELECTRON_RADIUS
#define LGM_ELECTRON_RADIUS (2.8179403267e-15)      // Classical electron radius, m
#endif

#ifndef LGM_BOLTZMANN_CONSTANT
#define LGM_BOLTZMANN_CONSTANT (1.3806488e-23)      // Boltzmann's constant, J/K
#endif


#ifndef LGM_PLANCK_CONSTANT_H
#define LGM_PLANCK_CONSTANT_H (6.62606957e-34)      // Planck's constant h, J·s
#endif

#ifndef LGM_PLANCK_CONSTANT_HBAR
#define LGM_PLANCK_CONSTANT_HBAR (1.054571726e-34)  // Planck's constant h-bar = h/(2*pi), J·s
#endif







#ifndef LGM_EPS
#define LGM_EPS            (0.000544617021954)          // me/mp
#endif

#ifndef LGM_EPS2
#define LGM_EPS2           (2.966077006016041e-7)       // (me/mp)^2
#endif

#ifndef LGM_EPS3
#define LGM_EPS3           (1.615376025901495e-10)      // (me/mp)^3
#endif

#ifndef LGM_EPS4
#define LGM_EPS4           (8.797612805617079e-14)      // (me/mp)^3
#endif

#ifndef LGM_EPS5
#define LGM_EPS5           (4.791329686495997e-17)      // (me/mp)^3
#endif

#ifndef LGM_OXYGEN_MASS
#define LGM_OXYGEN_MASS     (15.9994*LGM_U)             // Oxygen rest mass,  kg
#endif

#ifndef LGM_Ee0
#define LGM_Ee0            (0.510998910)                // Electron rest energy in MeV
#endif

#ifndef LGM_Ep0
#define LGM_Ep0             (938.27201323)              // Proton rest energy in MeV
#endif

#ifndef nT_Per_Gauss
#define nT_Per_Gauss        (1e5)                       // 1e5 nT == 1 Gauss
#endif

#ifndef Joules_Per_eV
#define Joules_Per_eV        (1.60217648740e-19)        // 1 eV = 1.60217648740e-19 Joules
#endif




#endif
