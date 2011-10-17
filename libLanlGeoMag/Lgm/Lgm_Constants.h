#ifndef LGM_CONSTANTS_H
#define LGM_CONSTANTS_H
#include <math.h>

#ifndef LGM_c
#define LGM_c               2.99792458e8                // Speed of light  m/s
#endif

#ifndef LGM_e
#define LGM_e               1.60217656535e-19           // Elementary charge, Coulombs (or s·A)
#endif


#ifndef LGM_Eps0
#define LGM_Eps0            8.854187817620e-12          // Permittivity of free space, A^2·s^4·kg^−1·m^−3
#endif

#ifndef LGM_U
#define LGM_U               (1.66053878283e-27)         // Atomic mass unit, kg
#endif

#ifndef LGM_ELECTRON_MASS
#define LGM_ELECTRON_MASS   9.1093821545e-31            // Rest mass of electron, kg
#endif

#ifndef LGM_PROTON_MASS
#define LGM_PROTON_MASS     (1.0072764667710*LGM_U)     // Proton rest mass,  kg
#endif

#ifndef LGM_EPS
#define LGM_EPS             0.000544617021983          // me/mp
#endif

#ifndef LGM_EPS2
#define LGM_EPS2            2.966077006336315e-7       // (me/mp)^2
#endif

#ifndef LGM_EPS3
#define LGM_EPS3            1.615376026163136e-10      // (me/mp)^3
#endif

#ifndef LGM_OXYGEN_MASS
#define LGM_OXYGEN_MASS     (15.9994*LGM_U)             // Oxygen rest mass,  kg
#endif

#ifndef LGM_Ee0
#define LGM_Ee0             0.510998910                 // Electron rest energy in MeV
#endif

#ifndef LGM_Ep0
#define LGM_Ep0             938.27201323                // Proton rest energy in MeV
#endif

#ifndef nT_Per_Gauss
#define nT_Per_Gauss        (1e5)                       // 1e5 nT == 1 Gauss
#endif

#ifndef Joules_Per_eV
#define Joules_Per_eV        (1.60217648740e-19)         // 1 eV = 1.60217648740e-19 Joules
#endif


#endif
