# get data from xyzf file and plot the graphics of
# diameter and number of vessels
# in function of time
# this program needs the f2py library


from get_data_m import *
from math import pi, sqrt
#from matplotlib import *
#from numpy import *
import os

def rc(r, Lsize):
    i = 0
    rn_private =  [0,0,0]
    while(i<=2):
        rn_private[i] = round(r[i]-Lsize[i])
        i+=1
    return rn_private

def rc_inv(r, Lsize):
    i = 0
    rn_private =  [0,0,0]
    while(i<=2):
        rn_private[i] = round(r[i]+Lsize[i])
        i+=1
    return rn_private


dir_name = '001'
time_interval = 2000
first_file_id = 10000
last_file_id = 34000
ifile = first_file_id
x,y,z,f = [],[],[],[]


Lsize =  [100,100,50]
#np = 8*Lsize[0]*Lsize[1]*Lsize[2]
#lxyz = [ [ 0 for j in range(4) ] for i in range(np) ]
#lxyz_inv = [[ [ 0 for k in range(60) ] for j in range(100) ] for i in range(100)]

branches = []

lx = Lsize[0]
ly = Lsize[1]
lz = Lsize[2]


volume_i = float(len(open(dir_name+'/phi 12000.xyz', 'r' ).readlines(  )))
length_i = len(open(dir_name+'/phis 12000.xyz', 'r' ).readlines()) 

file_output = open(dir_name+'/branches','w')
                   
while( ifile<=last_file_id):
    
    file_import = open(dir_name+'/phis '+str(ifile)+'.xyz', 'r' )
    volume = float(len(open(dir_name+'/phi '+str(ifile)+'.xyz', 'r' ).readlines(  )))
    length = len(open(dir_name+'/phis '+str(ifile)+'.xyz', 'r' ).readlines()) 
    x = []
    y = []
    z = []
    f = []
    for line in file_import:

        row = line.split()
        
        x.append(row[0])
        y.append(row[1])
        z.append(row[2])
        f.append(row[3])


    file_import.close() 
    nbranches = get_data_m.run_get_data(x,y,z,f,lx,ly,lz)
    
    #branches.append(nbranches)
    #print volume
    
    if(length>length_i):
        volume = volume - volume_i
        length = length - length_i
    else:
        volume = 0
    print >> file_output, ifile, nbranches, volume, 2.0*(  sqrt( volume/(2*pi*length) ) )


    
    print ifile
    
    ifile += time_interval

file_output.close()


fgnu=os.popen('gnuplot' ,'w')
PLOTFILE = 'branches.pdf'
DATAFILE = dir_name+'/branches'
print >>fgnu, "set xlabel 'time step'"
print >>fgnu, "set terminal pdf"
print >>fgnu, "set grid"
print >>fgnu, "set key right bottom"
print >>fgnu, "set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 pi -1 ps 1"
print >>fgnu, "set pointintervalbox 2"
print >>fgnu, "set output '%s'" %PLOTFILE
print >>fgnu, "plot '%s' using 1:2 title 'branches' w linespoints ls 1  , '%s' using 1:4 title 'diameter' w linespoints ls 2  " %(DATAFILE, DATAFILE)
fgnu.flush()
print "End program"
    
"""
r = [89,50,20]

rn = crng(r,Lsize)

print rn

"""
