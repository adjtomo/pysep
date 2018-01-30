import obspy
import pandas
import os

# Give source parameters
depth = 23    # in km
mag = 5.0
strike = 280
dip = 80
rake = 1
dur = 1       # source duration (in seconds)
event_name = 'test'

velocity_model = 'tactmod'
greens_fun_path = '/store/wf/FK_synthetics'

station_file = '~/REPOSITORIES/pyseis/test_data/20180123093142000_station_list_test.dat'

#--------------------------------------------------------
# Use syn to generate synthetics
'''
Usage: syn -Mmag([[/Strike/Dip]/Rake]|/Mxx/Mxy/Mxz/Myy/Myz/Mzz) -Aazimuth ([-SsrcFunctionName | -Ddura[/rise]] [-Ff1/f2[/n]] [-I | -J] -OoutName.z -GFirstCompOfGreen | -P)
   Compute displacement in cm produced by difference seismic sources
   -M Specify source magnitude and orientation or moment-tensor
      For double-couple, mag is Mw, strike/dip/rake are in A&R convention
      For explosion; mag in in dyne-cm, no strike, dip, and rake needed
      For single-force source; mag is in dyne, only strike and dip are needed
      For moment-tensor; mag in dyne-cm, x=N,y=E,z=Down
   -A Set station azimuth in degree measured from the North
   -D Specify the source time function as a trapezoid,
      give the total duration and rise-time (0-0.5, default 0.5=triangle)
   -F apply n-th order Butterworth band-pass filter, SAC lib required (off, n=4, must be < 10)
   -G Give the name of the first component of the FK Green function
   -I Integration once
   -J Differentiate the synthetics
   -O Output SAC file name
   -P Compute static displacement, input Green functions from stdin in the form
	distance Z45 R45 T45 ZDD RDD TDD ZSS RSS TSS [distance ZEX REX TEX]
      The displacements will be output to stdout in the form of
	distance azimuth z r t
   -Q Convolve a Futterman Q operator of tstar (no)
   -S Specify the SAC file name of the source time function (its sum. must be 1)
   Examples:
   * To compute three-component velocity at N33.5E azimuth from a Mw 4.5
earthquake (strike 355, dip 80, rake -70), use:
	syn -M4.5/355/80/-70 -D1 -A33.5 -OPAS.z -Ghk_15/50.grn.0
   * To compute the static displacements from the same earthquake, use:
	nawk '$1==50' st.out | syn -M4.5/355/80/-70 -A33.5 -P
   * To compute displacement from an explosion, use:
   	syn -M3.3e20 -D1 -A33.5 -OPAS.z -Ghk_15/50.grn.a
      or
        syn -M3.3e20/1/0/0/1/0/1 -D1 -A33.5 -OPAS.z -Ghk_15/50.grn.0
'''

#------------------------------------------------

def read_station_file(filename):
    '''
    Read station file generated while extracting waveforms using pyseis
    '''
    data = pandas.io.parsers.read_table(
        filename, sep=r"\s+", header=None,
        names=["station", "network", "latitude", "longitude", "distance",
               "azimuth"])
    return data

#-------------------------------------------
# create event directory
if not os.path.exists(event_name):
    os.makedirs(event_name)

# greens function path
green = os.path.join(greens_fun_path,velocity_model) + '/' + velocity_model + '_' + str(depth)

# source mechanism
m = str(mag)+'/'+str(strike)+'/'+str(dip)+'/'+str(rake)

# Get station info
df = read_station_file(station_file)
#print(df)

# Loop over stations
for index, row in df.iterrows():
    zcomp = row['station'] + '.z'
    az = row['azimuth']
    dist = round(row['distance']) 
    

    # set syn flags
    syn_command = 'syn' + ' -M'+m + ' -D'+str(dur) + ' -O'+event_name+'/'+zcomp + ' -A'+str(az) + ' -G'+green+'/'+str(dist)+ '.grn.0'

    print(syn_command)

    # call syn
    os.system(syn_command)
