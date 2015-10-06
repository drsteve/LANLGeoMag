# README for LANLGeoMag/Examples/KdTree

##Example programs:
 - DstNNs.c. Finds similar magnetic storms by kNN search on Dst time series.
 - TLE_NNs.c. Finds physical conjunctions between satellites using TLEs to 
   generate orbits (Collision candidates)
 - magConj.c. Builds kd-tree from set of HDF5 MagEphem files, uses kNN to find
   magnetic conjunctions with input query points
 - test.c. Tests implementation of kd tree

## Notes for compiling (Oct 2015)
On recent Debian-based linux systems (and possibly others), the install 
locations for HDF5 have changed.
HDF5 have also (finally) adopted the pkg-config system, so once LGM is properly
compiled a sample command line for compiling an example in this folder would be

> gcc test.c `pkg-config --cflags --libs lgm hdf5` -o test

