#!/usr/bin/python
import sys
import os
os.putenv("CDF_LIB", "/usr/local/cdf/lib")
import spacepy.pycdf as cdf

#cdf_file        = cdf.CDF('~/Download/QinDenton_1min_merged_20101229-v2.cdf')
cdf_file        = cdf.CDF('~/Download/QinDenton_hour_merged_20120206-v5.cdf')
Epoch           = cdf_file['Epoch']
Year         = cdf_file['Year']
DOY          = cdf_file['DOY']
hour         = cdf_file['hour']
min          = cdf_file['min']
ByIMF        = cdf_file['ByIMF']
BzIMF        = cdf_file['BzIMF']
V_SW         = cdf_file['V_SW']
Den_P        = cdf_file['Den_P']
Pdyn         = cdf_file['Pdyn']
G1           = cdf_file['G1']
G2           = cdf_file['G2']
G3           = cdf_file['G3']
ByIMF_status = cdf_file['ByIMF_status']
BzIMF_status = cdf_file['BzIMF_status']
V_SW_status  = cdf_file['V_SW_status']
Den_P_status = cdf_file['Den_P_status']
Pdyn_status  = cdf_file['Pdyn_status']
G1_status    = cdf_file['G1_status']
G2_status    = cdf_file['G2_status']
G3_status    = cdf_file['G3_status']
Kp           = cdf_file['Kp']
akp3         = cdf_file['akp3']
Dst          = cdf_file['Dst']
Bz1          = cdf_file['Bz1']
Bz2          = cdf_file['Bz2']
Bz3          = cdf_file['Bz3']
Bz4          = cdf_file['Bz4']
Bz5          = cdf_file['Bz5']
Bz6          = cdf_file['Bz6']
W1           = cdf_file['W1']
W2           = cdf_file['W2']
W3           = cdf_file['W3']
W4           = cdf_file['W4']
W5           = cdf_file['W5']
W6           = cdf_file['W6']
W1_status    = cdf_file['W1_status']
W2_status    = cdf_file['W2_status']
W3_status    = cdf_file['W3_status']
W4_status    = cdf_file['W4_status']
W5_status    = cdf_file['W5_status']
W6_status    = cdf_file['W6_status']

saveout = sys.stdout


old_doy  = -999.9
old_Year = -999
for i, doy in enumerate(DOY):


    e      = Epoch[i]
    Year   = e.year
    Month  = e.month
    Day    = e.day
    Hour   = e.hour
    Minute = e.minute
    Second = e.second

    if (Year != old_Year):
        dirname = "{0:4d}".format( Year )
        old_Year = Year
        if not os.path.exists(dirname):
            print "Making Directory:",  dirname
            os.makedirs(dirname)


   

    if (doy != old_doy):
        # Close Old file and OIpen a new one
        if (old_doy != -999.9):
            f.close()
        # Open New file
        Filename = "{0:4d}".format( Year ) + "/QinDenton_{0:4d}{1:02d}{2:02d}_1hr.txt".format( Year, Month, Day )
        sys.stdout = saveout
        print "Writing File:",  Filename
        f = open(Filename, 'w')
        sys.stdout = f
        old_doy = doy





    print e.isoformat(), Year, Month, Day, Hour, Minute, Second, ByIMF[i], BzIMF[i], V_SW[i], Den_P[i], Pdyn[i], G1[i], G2[i], G3[i], int(ByIMF_status[i]), int(BzIMF_status[i]), int(V_SW_status[i]), int(Den_P_status[i]), int(Pdyn_status[i]), int(G1_status[i]), int(G2_status[i]), int(G3_status[i]), Kp[i], akp3[i], Dst[i], Bz1[i], Bz2[i], Bz3[i], Bz4[i], Bz5[i], Bz6[i], W1[i], W2[i], W3[i], W4[i], W5[i], W6[i], int(W1_status[i]), int(W2_status[i]), int(W3_status[i]), int(W4_status[i]), int(W5_status[i]), int(W6_status[i])

    

