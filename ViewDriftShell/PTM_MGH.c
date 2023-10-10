#include <mpi.h>
#include <limits.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include "ParticleInfo.h"

#define TRUE 1
#define FALSE 0

#define BUSY 1
#define IDLE 0

#define NOT_STARTED 0
#define SENT        1
#define DONE        2


/*
 * Allocate memory that will be shared across the ranks in a node
 *
 *      bufsize: This is the number of bytes to allocate
 *          buf: This is a (void *) pointer the the chunck of memory
 */
int AllocSharedMem(  size_t bufsize, void **mem, MPI_Win *winbuf, int wrank, int nrank, MPI_Comm NodeComm ) {

    int     localbufsize;
    void    *buf;
    void    *localbuf;

    int flag;
    int *model;
    MPI_Aint winsize;
    int windisp;
    int *winptr;

    /*
     * This is a collective routine -- all of the members of NodeComm will call
     * it.  When calling MPI_Win_allocate_shared(), each member rank can
     * contribute to the overall size of the buffer, or just one rank can
     * allocate all the memory while the others contribute zero bytes.  In
     * either case, the end result is a shared chunk of memory.  (It is
     * possible that having all ask for some could be more efficient? Might
     * depend on how memory is attached to the hardware the proc is running
     * on.)
     *
     * For now, we will just have nrank=0 request the full size and all others
     * will request 0 size. This is what localbufsize is for (it'll be full
     * size for nrank=0, and zero for the other ranks.)
     */
    localbufsize = ( nrank == 0 ) ? bufsize : 0; 
    MPI_Win_allocate_shared( localbufsize, 1, MPI_INFO_NULL, NodeComm, &localbuf, winbuf );
    MPI_Win_get_attr( *winbuf, MPI_WIN_MODEL, &model, &flag );

    if ( flag != 1 ) { 
        printf("Attribute MPI_WIN_MODEL not defined\n");
    } else { 
        if (MPI_WIN_UNIFIED == *model) {
            if ( wrank == 0 ) printf("Memory model is MPI_WIN_UNIFIED\n");
        } else {
            if ( wrank == 0 ) printf("Memory model is *not* MPI_WIN_UNIFIED\n");
            MPI_Finalize();
            return 1;
        }
    }
    // on nrank=0, buf is just the local buf
    buf = localbuf; 
    // on other nranks, need to query to get the pointer we should use
    if ( nrank != 0 ) { MPI_Win_shared_query( *winbuf, 0, &winsize, &windisp, &buf ); }

    *mem = buf;

    // All buf pointers should now point to copy on noderank 0
    MPI_Win_fence( 0, *winbuf ); // syncronize

}

int main(int argc, char** argv) {

    ParticleInfo JobList[1000];
    int wsize, wrank;
    int nsize, nrank;
    int mycolor;

    MPI_Status status;
    MPI_Comm AllComm, NodeComm;

    double start, end;

    AllComm = MPI_COMM_WORLD;

    // Initialize the MPI environment
    MPI_Init( &argc, &argv );

    // Get the total number of processes and the world rank
    MPI_Comm_size( AllComm, &wsize );
    MPI_Comm_rank( AllComm, &wrank );


    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);


    // Split up the world communicator by shared memory. (I.e. by node)
    MPI_Comm_split_type( AllComm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &NodeComm );
    MPI_Comm_size( NodeComm, &nsize ); // Number of processors in this shared mem group
    MPI_Comm_rank( NodeComm, &nrank ); // My rank in this shared mem group

    if (nrank == 0 ) mycolor = wrank;  // make the "color" be the wrank of the nrank=0 process.
    //if ( nrank == 0 ){ printf("I am nrank = 0. I am setting mycolor to be my world rank (wrank = %d). I am broadcasting this to others in my NodeComm\n", wrank); }
    MPI_Bcast( &mycolor, 1, MPI_INT, 0, NodeComm ); // Broadcast mycolor to others in this communicator
    //if ( nrank != 0 ) { printf( "I am nrank = %d. I am receiving mycolor (mycolor = %d) (my wrank = %d).\n", nrank, mycolor, wrank ); }
    MPI_Barrier( NodeComm );


    // at this point, every proc in my NodeComm will have mycolor set to the
    // world rank of the local nrank=0 proc
    // Let's Gather all of the colors into the root proc
    int *ProcColors;
    ProcColors = (int *)calloc( wsize, sizeof(int) );
    MPI_Gather( &mycolor, 1, MPI_INT,   ProcColors, wsize, MPI_INT,    0, AllComm );

    int i, j;
    if ( wrank == 0 ){
        for ( i=0; i<wsize; ++i ) {
            printf("ProcColors[%d] = %d\n", i, ProcColors[i] );
        }
    }
//this machine has only 1 node. Lets kludge things up to make it looks as though we have 6 nodes
if (0==1){
ProcColors[0] = 0;
ProcColors[1] = 12;
ProcColors[2] = 12;
ProcColors[3] = 0;

ProcColors[4] = 4;
ProcColors[5] = 4;
ProcColors[6] = 4;
ProcColors[7] = 4;

ProcColors[8] = 8;
ProcColors[9] = 8;
ProcColors[10] = 8;
ProcColors[11] = 8;

ProcColors[12] = 12;
ProcColors[13] = 0;
ProcColors[14] = 12;
ProcColors[15] = 12;

ProcColors[16] = 16;
ProcColors[17] = 16;
ProcColors[18] = 16;
ProcColors[19] = 16;

ProcColors[20] = 20;
ProcColors[21] = 20;
ProcColors[22] = 20;
ProcColors[23] = 20;
}

    int done, c, LastMin, Min, n, Proc[100][100], nProc[100];
    if ( wrank == 0 ){
        // The array ProcColors[] has wsize members. each index corresponds to a world rank.
        // the value stored is the mycolor values (it will be the world rank of the nrank=0 proc in each node.)
        // What we want to do now is sift through this to create a 2D array
        // Proc[i][j] where i denotes a node and j is a proc within the node.
        n = 0; // assigned node number 
        done = FALSE; LastMin = -1;
        while( !done ) {

            // Find minimum value in the ProcColors[] array
            Min = INT_MAX;
            for ( i=0; i<wsize; ++i ) { 
                c = ProcColors[i];
                if ( (c>LastMin) && (c < Min ) ) Min = c; 
                //printf("c = %d, n = %d,  LastMin = %d,  Min = %d\n", c, n, LastMin, Min);
            }
            LastMin = Min;
            //printf("\n");

            if ( Min == INT_MAX ) {
                // no values left to process -- we're done
                done = TRUE;
                break;
            } else {
                nProc[n] = 0;
                for ( i=0; i<wsize; ++i ) { 
                    c = ProcColors[i];
                    if ( c == Min ) {
                        Proc[n][ nProc[n] ] = i; // add another processor (i is the world_rank)
                        ++nProc[n];
                    }
                }
            }
            ++n; // increment the node number;

        } // end while loop

    }

    for ( j=0; j<n; ++j ) {
        printf("Node %d has the following %d shared mem processes in it:\n", j, nProc[j] );
        for ( i=0; i<nProc[j]; ++i ) {
            printf("\tProc[%d][%d] = %d (this number is the process's world_rank)\n", j, i, Proc[j][i] );
        }
    }


    // At this point the array Proc[j][i] holds the world_ranks of all the processes
    // The first index specifies a shared memory node, the second specifies process number within that node.
    // This now allows us to partition up the nodes by shared memory.....

    



    /*
     *  We need to allocate memory for the E(x,y,z,t) and B(x,y,z,t) files
     *  On each node we can only fit a limited number of time slices.
     * 
     *  Each array is 4D. To determine how many nodes we need, we need to know;
     *    1) How large is a single timeslice?
     *    2) How many timeslices do we have?
     *    3) How much memory do we have on a node?
     *
     *  E.g., if there are 90e6 cells, we have ~10*90e6*8bytes = 7.2GB per time slice
     *  on grizzly, each node has 128GB/node ~ 17 time slices/node.
     *  However, there are other variables that are needed.
     *  For AMR, we could also make use of the BATL tree library, but it just
     *  takes a structure as input (I would guess).
     * 
     *  If we haveb 4 hours at 1s resolution, thats 3600*4 = 14400 slices. 
     *  We'd need 14400/17 nodes = 847 nodes on grizzy. Grizzly has 1490 nodes.
     * 
     *  On Chicoma, memory is 512G/node so 70 slices/node.
     *  would need 14400/70 = 205 nodes. Chicoma has 560 nodes.
     * 
     * 
     *  Although we proibably eventually want to use the full AMR files with
     *  the BATL library, lets just add code to allocate simple shared 4D
     *  arrays for now
     * 
     * 
     * 
     * 
     */


    size_t bufsize;
    int *buf;
    void   *mem;
    MPI_Win winbuf;

    bufsize = 5*sizeof(int);
    AllocSharedMem( bufsize, &mem, &winbuf, wrank, nrank, NodeComm ); // call to get shared mem across nranks on a node
    buf = (int *)mem;



    /*
     *  nrank=0 fills the array
     */
    if ( nrank == 0 ) { for ( i=0; i < bufsize; i++ ) { buf[i] = i*bufsize + i; } }
    MPI_Win_fence( 0, winbuf );



    /*
     * Check that it worked. Results should be same on all ranks
     */ 
//    for (i=0; i < bufsize; i++) { printf("wrank %d, nrank %d, buf[%d] = %d\n", wrank, nrank, i, buf[i]); }


    /*
     * Set up a list of jobs to do
     */
    int JobStatus[1000];
    if ( wrank == 0 ) {
        for ( i=0; i<1000; ++i ) {
            JobList[i].ID            = i;
            JobList[i].Master_wrank  = wrank;
            JobList[i].Master_nrank  = nrank;
            JobList[i].Worker_wrank  = -1; // not assigned yet
            JobList[i].Worker_nrank  = -1; // not assigned yet

            JobList[i].Species       = PROTON;
            JobList[i].KineticEnergy = 0.5 + (double)i/1000.0; //MeV
            JobList[i].Charge        = 1.0;

            JobList[i].dT              = 3600.0;   // seconds
            JobList[i].SwitchingMethod = ADAPTIVE; // allow switchingm between GCA nd Full orbit equations
            JobList[i].SaveTrajectory  = FALSE;    // dont save trajectory
            JobList[i].dt              = 0.0;      // not used
            
            JobList[i].TracingType_in  = GCA;      // start by using GC eqns
            JobList[i].t_in            = 0.0;       // start time

            JobList[i].x_in[0] = JobList[i].x_in[1] = JobList[i].x_in[2] = 0.0; // not used
            JobList[i].v_in[0] = JobList[i].v_in[1] = JobList[i].v_in[2] = 0.0; // not used

            JobList[i].X_in[0]  = -6.6;   // Initial GC position in Re
            JobList[i].X_in[1]  =  0.0;   // Initial GC position in Re
            JobList[i].X_in[2]  =  0.0;   // Initial GC position in Re
            JobList[i].v_par_in = 1234.5; // need proper units...
            JobList[i].mu_in    = 67.89;  // need proper units...
            
        }

        JobStatus[i] = NOT_STARTED; // Flag job as "not started"

    }

    ParticleInfo Job;
    MPI_Datatype mpi_type;
    MpiDataType( &Job, &mpi_type ); // This create a data type from the ParticleInfo structure for MPI to send/receice using its rotuines.
    

    /*
     *   Enter Master/Worker loops
     */
    MPI_Request request;
    int received;
    int MessageWaiting;
    if ( wrank == 0 ) {

        /*
         * I am the master
         */

        // for e.g.
        int nWorkers;
        int nWorking;
        int WorkerToUse[36];
        int WorkerStatus[36];
        nWorkers = 0;
        for ( i=0; i<nProc[0]; i++ ) {
            if ( Proc[0][i] != wrank ) {
                WorkerToUse[nWorkers]  = Proc[0][i];
                WorkerStatus[nWorkers] = IDLE;
                ++nWorkers;
            }
        }


        done = 0;
        int k = 0;
        nWorking = 0;
        while ( !done ) {


            /*
             * loop over target workers and send jobs if they are idle.
             */
            for ( j=0; j<nWorkers; j++ ) {
                if ( WorkerStatus[j] == IDLE ) {

                    if ( k >= 2 ) {
                        /*
                         *  No more jobs to send, but we couls still be weaiting for results.
                         */
                        if ( nWorking == 0 ) {
                            // End worker procs
                            Job.ID = -1;
                            for ( i=0; i<nWorkers; i++ )  {
                                MPI_Send( &Job, 1, mpi_type, WorkerToUse[i], 0, AllComm );
                            }
                            done = TRUE;
                        }
                    } else {
                        /*
                         *  Send it the next job
                         */
                        Job = JobList[k]; // has to be copied into the mem space of Job to transmmit across MPI properly
                        Job.Worker_wrank = WorkerToUse[j];
                        printf("Master sending worker %d (wrank=%d), job %d\n", j, WorkerToUse[j], k);
                        MPI_Send( &Job, 1, mpi_type, WorkerToUse[j], 0, AllComm );
                        JobStatus[k] = SENT; // mark it as "SENT"
                        WorkerStatus[j] = BUSY;
                        ++nWorking;
                        ++k;
                    }

                }
            }


            /*
             *  Check for a message. FIC  -- We need to get the result back not just a message syaing its done...
             */
            MessageWaiting = FALSE;
            MPI_Iprobe( MPI_ANY_SOURCE, 0, AllComm, &MessageWaiting, &status );
            if ( MessageWaiting ) {
                MPI_Recv( &received, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status );
                 
                printf("Master recieved message for job #%02d from wrank=%d.\n", received, status.MPI_SOURCE );
                for (j=0; j<nWorkers; j++) { 
                    // dumb way to do this. We really need more accounting of tasks and who we sent it to etc...
                    if ( WorkerToUse[j] == status.MPI_SOURCE ){
                        WorkerStatus[j] = IDLE;
                        JobStatus[ received ] = DONE;
                        --nWorking;
                    }
    
                }
                
            }
            

           



        }





    } else {

        /*
         * I am a worker
         */

        start = MPI_Wtime(); // get wall time start
        done = 0;
        while ( !done ) {
sleep(1);

            // accept new job
            //printf("MPI process %d here.\n", wrank);
            MPI_Recv( &Job, 1, mpi_type, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
            
            // Determine if we should quit.
            if (Job.ID < 0 ) {

                done = TRUE;

            } else {

                printf("MPI process %d working on job #%02d.  Energy to do: %g\n", wrank, Job.ID, Job.KineticEnergy );


                // do work
                double ddd;
                int ii, jj;
                for (jj=0; jj<1000000; ++jj){
                for (ii=0; ii<1000; ++ii){
                    ddd += tan( cos((double)ii/10.0) );
                }
                }
                //sleep(3);

                // send back results
                printf("MPI process %d done working on job #%02d.\n", wrank, Job.ID );
                MPI_Send( &Job.ID, 1, MPI_INT, 0, 0, AllComm );

            }


        }
        end = MPI_Wtime();
        printf("Work took %lf seconds to run.\n", end-start );



    }




    // Finalize the MPI environment.
    MPI_Type_free( &mpi_type );
    MPI_Finalize();
}

