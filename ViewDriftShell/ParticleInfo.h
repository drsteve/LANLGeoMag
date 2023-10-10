#ifndef PARTICLE_INFO_H
#define PARTICLE_INFO_H


#define ELECTRON 0
#define PROTON   1

#define GCA_ONLY          0
#define FULL_ORBIT_ONLY   1
#define ADAPTIVE          2

#define GCA         0
#define FULL_ORBIT  1

typedef struct ParticleInfo {

    /*
     * Some housekeeping information to identify the job and where it will run
     * etc..
     */
    int     ID;
    int     Master_wrank; // World rank of the master process
    int     Master_nrank; // Node rank of the master process
    int     Worker_wrank; // world rank of the wroker process (i.e. where the tracing will be done)
    int     Worker_nrank; // Node rank of the wroker process (i.e. where the tracing will be done)


    int     Species;            // ELECTRON, PROTON, etc...
    double  KineticEnergy;      // keep track of KE
    double  Charge;             // In units of fundametal charge (i.e. electron has charge -1.0, proton has +1.0)
    

    double  dT;                 // Time to integrate for. I.e., t_out = t_in + dT
    int     SwitchingMethod;    // E.g., GCA_ONLY, FULL_ORBIT_ONLY, ADAPTIVE
    int     SaveTrajectory;     // TRUE or FALSE. 
    double  dt;                 // interval to save trajectory.
    
    /*
     *  Starting conditions
     */
    int     TracingType_in;             // TracingType E.g. GCA, FULL_ORBIT
    double  t_in;                       // starting time
    double  x_in[3], v_in[3];           // Full Orbit coordinates
    double  X_in[3], v_par_in, mu_in;   // Guiding Center coordinates

    /*
     *  Ending conditions 
     */
    int     TracingType_out;            // TracingType E.g. GCA, FULL_ORBIT (might have change if SwitchingMethod==DYNAMIC)
    double  t_out;                      // endingm time
    double  x_out[3], v_out[3];         // Full Orbit coordinates
    double  X_out[3], v_par_out, mu_out;// Guiding Center coordinates



} ParticleInfo;


int MpiDataType( ParticleInfo *Job, MPI_Datatype *mpi_type );


#endif
