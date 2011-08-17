#!/usr/bin/python
import sys
import os
os.putenv("CDF_LIB", "/usr/local/cdf/lib")
import spacepy.pycdf as cdf

cdf_file        = cdf.CDF('~/Download/QinDenton_1min_merged_20101229-v2.cdf')
Epoch           = cdf_file['Epoch']
Year_v2         = cdf_file['Year_v2']
DOY_v2          = cdf_file['DOY_v2']
hour_v2         = cdf_file['hour_v2']
min_v2          = cdf_file['min_v2']
ByIMF_v2        = cdf_file['ByIMF_v2']
BzIMF_v2        = cdf_file['BzIMF_v2']
V_SW_v2         = cdf_file['V_SW_v2']
Den_P_v2        = cdf_file['Den_P_v2']
Pdyn_v2         = cdf_file['Pdyn_v2']
G1_v2           = cdf_file['G1_v2']
G2_v2           = cdf_file['G2_v2']
G3_v2           = cdf_file['G3_v2']
ByIMF_status_v2 = cdf_file['ByIMF_status_v2']
BzIMF_status_v2 = cdf_file['BzIMF_status_v2']
V_SW_status_v2  = cdf_file['V_SW_status_v2']
Den_P_status_v2 = cdf_file['Den_P_status_v2']
Pdyn_status_v2  = cdf_file['Pdyn_status_v2']
G1_status_v2    = cdf_file['G1_status_v2']
G2_status_v2    = cdf_file['G2_status_v2']
G3_status_v2    = cdf_file['G3_status_v2']
Kp_v2           = cdf_file['Kp_v2']
akp3_v2         = cdf_file['akp3_v2']
Dst_v2          = cdf_file['Dst_v2']
Bz1_v2          = cdf_file['Bz1_v2']
Bz2_v2          = cdf_file['Bz2_v2']
Bz3_v2          = cdf_file['Bz3_v2']
Bz4_v2          = cdf_file['Bz4_v2']
Bz5_v2          = cdf_file['Bz5_v2']
Bz6_v2          = cdf_file['Bz6_v2']
W1_v2           = cdf_file['W1_v2']
W2_v2           = cdf_file['W2_v2']
W3_v2           = cdf_file['W3_v2']
W4_v2           = cdf_file['W4_v2']
W5_v2           = cdf_file['W5_v2']
W6_v2           = cdf_file['W6_v2']
W1_status_v2    = cdf_file['W1_status_v2']
W2_status_v2    = cdf_file['W2_status_v2']
W3_status_v2    = cdf_file['W3_status_v2']
W4_status_v2    = cdf_file['W4_status_v2']
W5_status_v2    = cdf_file['W5_status_v2']
W6_status_v2    = cdf_file['W6_status_v2']

saveout = sys.stdout


old_doy = -999
for i, doy in enumerate(DOY_v2):

    e      = Epoch[i]
    Year   = e.year
    Month  = e.month
    Day    = e.day
    Hour   = e.hour
    Minute = e.minute
    Second = e.second

    if (doy != old_doy):
        # Close Old file and OIpen a new one
        if (old_doy != -999):
            f.close()
        # Open New file
        Filename = "QinDenton_{0:4d}{1:02d}{2:02d}_1min.txt".format( Year, Month, Day )
        sys.stdout = saveout
        print "Writing File:",  Filename
        f = open(Filename, 'w')
        sys.stdout = f
        old_doy = doy


    print e.isoformat(), Year, Month, Day, Hour, Minute, Second, ByIMF_v2[i], BzIMF_v2[i], V_SW_v2[i], Den_P_v2[i], Pdyn_v2[i], G1_v2[i], G2_v2[i], G3_v2[i], int(ByIMF_status_v2[i]), int(BzIMF_status_v2[i]), int(V_SW_status_v2[i]), int(Den_P_status_v2[i]), int(Pdyn_status_v2[i]), int(G1_status_v2[i]), int(G2_status_v2[i]), int(G3_status_v2[i]), Kp_v2[i], akp3_v2[i], Dst_v2[i], Bz1_v2[i], Bz2_v2[i], Bz3_v2[i], Bz4_v2[i], Bz5_v2[i], Bz6_v2[i], W1_v2[i], W2_v2[i], W3_v2[i], W4_v2[i], W5_v2[i], W6_v2[i], int(W1_status_v2[i]), int(W2_status_v2[i]), int(W3_status_v2[i]), int(W4_status_v2[i]), int(W5_status_v2[i]), int(W6_status_v2[i])

    

