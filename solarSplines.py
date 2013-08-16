"""
solarSplines.py
Authors: David John Gagne II, Andrew MacKenzie, Darrel Kingfield, Ben Herzog, Tim Humphrey, Thibault Lucidarme
Description: Reads Mesonet and GEFS data and applies Catmull-Rom splines to interpolate GEFS data to Mesonet sites.
"""
from netCDF4 import Dataset
import numpy as np

def loadMesonetData(filename,stationFilename="station_info.csv"):
    """
    loadMesonetData(filename,stationFilename)
    Description: loads Mesonet data and station data
    Parameters:
    filename (str) - Name of Mesonet csv file being read.
    stationData (str) - Name of file containing station information. Default station_info.csv.
    Returns:
    data - numpy array containing total daily solar radiation for each date
    dates - numpy array of dates as integers in YYYYMMDD format
    stationData - numpy structured array containing the station information, including lat-lon and elevation.
    """
    data = np.genfromtxt(filename,delimiter=',',skiprows=1,dtype=float)
    dates = np.array(data[:,0].T,dtype=int)
    data = data[:,1:]
    stationData = np.genfromtxt(stationFilename,delimiter=',',dtype=[("stid","S4"),("nlat",float),("elon",float),("elev",float)],skiprows=1)
    return data,dates,stationData

def loadData(filename):
    """
    loadData()
    Description: Creates a netCDF4 file object for the specified file.
    Parameters:
    filename (str) - name of the GEFS netCDF4 file.
    Returns:
    data - Dataset object that allows access of GEFS data.
    """
    data = Dataset(filename)
    return data

def getGrid(data,date,fHour,eMember):
    """
    getGrid()
    Description: Load GEFS data from a specified date, forecast hour, and ensemble member.
    Parameters:
    data (Dataset) - Dataset object from loadData
    date (int) - date of model run in YYYYMMDD format
    fHour (int) - forecast hour
    eMember (int) - ensemble member id.
    Returns: numpy 2d array from the specified output
    """
    dateIdx = np.where(data.variables['intTime'][:] == date)[0]
    fIdx = np.where(data.variables['fhour'][:] == fHour)[0]
    eIdx = np.where(data.variables['ens'][:] == eMember)[0]
    return data.variables.values()[-1][dateIdx,eIdx,fIdx][0]

def getDailyMeanSumGrid(data,date):
    """
    getDailyMeanSumGrid()
    Description: For a particular date, sums over all forecast hours for each ensemble member then takes the 
    mean of the summed data and scales it by the GEFS time step.
    Parameters:
    data (Dataset) - netCDF4 object from loadData
    date (int) - date of model run in YYYYMMDD format
    Returns - numpy 2d array from the specified output
    """
    dateIdx = np.where(data.variables['intTime'][:] == date)[0]
    fIdx = np.where(data.variables['fhour'][:] <= 24)[0]
    return data.variables.values()[-1][dateIdx,:,fIdx,:,:].sum(axis=2).mean(axis=1)[0] * 3600 * 3

def buildSplines(data,grid,stationdata):
    """
    buildSplines()
    Description: For each station in stationdata, a set of Catmull-Rom splines are calculated to interpolate from the
    nearest grid points to the station location. A set of horizontal splines are created at each latitude and 
    interpolated at the station longitude. Then another spline is built from the output of those splines to get the 
    value at the station location.
    Paramters:
    data (Dataset) - netCDF4 object with the GEFS data 
    grid (numpy array) - the grid being interpolated
    stationdata (numpy structured array) - array containing station names, lats, and lons.
    Returns: array with the interpolated values.
    """
    outdata=np.zeros(stationdata.shape[0])
    print stationdata.shape
    for i in xrange(stationdata.shape[0]):
        slat,slon=stationdata['nlat'][i],stationdata['elon'][i]
        nearlat=np.where(np.abs(data.variables['lat'][:]-slat)<2)[0]
        nearlon=np.where(np.abs(data.variables['lon'][:]-slon-360)<2)[0]
        Spline1=np.zeros(nearlon.shape)
        for l,lat in enumerate(nearlat):
            Spline1[l]=Spline(grid[nearlat[l],nearlon],(slon-np.floor(slon))/1)
        outdata[i]=Spline(Spline1,(slat-np.floor(slat))/1)
    return outdata
       
def Spline(y,xi):
    """
    Spline
    Description: Given 4 y values and a xi point, calculate the value at the xi point.
    Parameters:
    y - numpy array with 4 values from the 4 nearest grid points
    xi - index at which the interpolation is occurring.
    Returns: yi - the interpolated value
    """
    return 0.5*((2*y[1])+(y[2]-y[0])*xi+(-y[3]+4*y[2]-5*y[1]+2*y[0])*xi**2+(y[3]-3*y[2]+3*y[1]-y[0])*xi**3)

def main():
    dataPath = "../"
    MesonetData,dates,stationdata = loadMesonetData(dataPath + 'sampleSubmission.csv',dataPath + "station_info.csv")
    data = loadData(dataPath + "dswrf_sfc_latlon_subset_20080101_20121130.nc")
    lons = data.variables['lon'][:]
    lats = data.variables['lat'][:]
    f = open('spline_submission.csv', 'w')
    header = ["Date"]
    header.extend(stationdata['stid'].tolist())
    f.write(",".join(header) + "\n")
    for date in dates:
        print date
        grid = getDailyMeanSumGrid(data,date*100)
        outdata=buildSplines(data,grid,stationdata)
        outdata = outdata.reshape(outdata.shape[0],1)
        f.write("%d" % date + ",")
        np.savetxt(f, outdata.T, delimiter=',',fmt='%7.0f')
    f.close()
    data.close()

if __name__ == "__main__":
    main()
