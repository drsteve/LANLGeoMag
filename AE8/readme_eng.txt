Input file should contains lines of following format:

L    B/B0  MODEL  FLUXTYPE  E1    E2

Example:

2.15  2.5    2       1      0.5   0.7  

Where: 

L - McIlwein parameter, 

B/B0 - fraction of magnetic field strength in the position of observer
and in the position of field line top.

MODEL - one of following integer numbers, flux model:
    1   AP8MAX, protons, solar maximum
    2   AP8MIN, protons, solar minimum
    3   AE4MAX, electrons, solar maximum
    4   AE4MIN, electrons, solar minimum
    5   AEI7HI, electrons, high estimate
    6   AEI7LO, electrons, low estimate
    7   AE8MAX, electrons, solar maximum
    8   AE8MIN, electrons, solar minimum

FLUXTYPE - one of following integer numbers, selection of flux type:
    1 - integral
    2 - differential

E1 - for integral flux - threshold of energy (MeV), to calculate Flux(>E1), 
for differential - lowest threshold of energy (MeV) to calculate Flux(E1<E<E2)

E2 - for integral flux program doesn't take this number into account
for differential - highest threshold of energy (MeV) to calculate Flux(E1<E<E2)

As the result of program work file of follwoing format will be generated:

L         B/B0      MODEL     FLUXTYPE  E1        E2         FLUX1 
2.1500001 2.5000000 2.0000000 1.0000000 0.5000000 0.7000000  1.2473444E+006

Where L, B/B0, MODEL, FLUXTYPE, E1, E2 are the same as in input file,
but integer numbers obtain decimal part.

The last column FLUX1 means flux calculated using model of selected name
Integral flux units are particles/(cm*cm*sec), 
differential flux units are particles/(cm*cm*sec*MeV)
