# -*- coding: utf-8 -*-
import netCDF4
import numpy as np
import sys
sys.path.append('/Users/imerino/Documents/These/python/data')
from data import *


factor= 1e12/365/24/3600 #Convert from GT/a to Kg/s

#---Clean the bathymetry applying a threshold. 0 in land points, value in ocean points.
def cleanBathy(var,value):
    var = var*(var==0.0) + value*(var>0.0)
    return var



# getIceShelfSegmentPoints
#---Get the coastal points belonging to the segments of an ice shelf
#   Inputs:
#       IceShelfLimits:Limits of an ice shelf at its segment [lon1,lat1,lon2,lat2,...,lon2n,lat2n]
#       xP: List of X grid coordinates of the Coastal points previously extracted
#       yP: List of Y grid coordinates of the Coastal points previously extracted
#       varLon: Longitude field over the same grid where XP and YP are defined
#       varLat: Latitude field over the same grid where XP and YP are defined
#
#   Output:
#       XIceShelf:list of X grid coordinates of the ice shelf extracted from xP
#       YIceShelf:list of Y grid coordinates of the ice shelf extracted from yP
#       PointsPerSegment: list with NumOfSegments values of the number of points of each segment composing the ice shelf
#       NumOfSegments: Number of segments composing the ice shelf
  
def getIceShelfSegmentPoints(IceShelfLimits,xP,yP,varLon,varLat):
    
        #Needed data to fill up per Ice Shelf 
        XIceShelf=[] #Container of X grid position of all the segments of an ice shelf
        YIceShelf=[] #Container of Y grid position of all the segments of an ice shelf
        NumOfSegments=0 #Number of segments
        PointsPerSegment=[] #Grid points per segment (dim=NumOfSegments)
    
        if np.mod(np.size(IceShelfLimits),4)!=0:
            print ('ERROR in limits dimension, hsould be x4') #Every segment is deined by 4 values (lonInit,latInit,lonFin,latFin)
            print (limits)
            
        NumOfSegments = np.size(IceShelfLimits)/4 #number of segments
        for i in (np.arange(NumOfSegments)-1): #Loop over segments
            lonInit=IceShelfLimits[int(i)*4]
            latInit=IceShelfLimits[int(i)*4+1]
            lonFin=IceShelfLimits[int(i)*4+2]
            latFin=IceShelfLimits[int(i)*4+3]
            if lonInit > 80:
               lonInit = lonInit - 360
            if lonFin > 80:
               lonFin = lonFin - 360
            #Detect the nearest coastal point belonging to global list xP,yP
            indexInit = getNearestPoint(xP,yP,varLon,varLat,lonInit,latInit) 
            indexFin = getNearestPoint(xP,yP,varLon,varLat,lonFin,latFin)
            #listX and listY constains X,Y coastal grip positions of the current segment
            listX,listY = getSegment(xP,yP,indexInit,indexFin)
            #Filling up
            XIceShelf=XIceShelf+listX
            YIceShelf=YIceShelf+listY
            PointsPerSegment.append(np.size(listX))
            #Just for plotting purpose
            #im=plt.scatter(listX, listY,s=10, marker='s',edgecolor='none',c=colors[indexColor],label=i)  
            
        return XIceShelf,YIceShelf,PointsPerSegment,NumOfSegments


#createIceShelfFluxFile       
#---Filling up a grid field with freshwaterfluxes rates (kg/m2/s) from a list of ice shelves or sectors
#   Inputs:
#       varRunoff:Variable to fill up with the calving rates
#       listIceShelves: list of ice shelves or sectors. Must be defined with the same name in data.py
#       xP: List of X grid coordinates of the Coastal points previously extracted
#       yP: List of Y grid coordinates of the Coastal points previously extracted
#       varLon: Longitude field over the same grid where XP and YP are defined
#       varLat: Latitude field over the same grid where XP and YP are defined
#       area: Grid file in the same dimensions than varRunoff with the area in meters
#
#   Output:
#       CalvingVar: the grid field filled up
def createIceShelfFluxFile(varRunoff,listIceShelves,xP,yP,varLon,varLat,area):
    for shelf in listIceShelves:
        #Needed data to fill up per Ice Shelf 
        XIceShelf=[] #Container of X grid position of all the segments of an ice shelf
        YIceShelf=[] #Container of Y grid position of all the segments of an ice shelf
        NumOfSegments=0 #Number of segments
        PointsPerSegment=[] #Grid points per segment (dim=NumOfSegments)
        
        limits=globals()[shelf]#Get geospatial limits of the iceshelf i from data.py. It may be made by segments
        
        if shelf in valuesIS:
            FWF=-1*valuesIS[shelf]
            print(shelf,FWF)
        elif shelf in  upscaling:
            FWF=-1*upscaling[shelf]
        else:
            print('ERROR shelf or sector not found')
        
        XIceShelf,YIceShelf,PointsPerSegment,NumOfSegments = getIceShelfSegmentPoints(limits,xP,yP,varLon,varLat)
        
        #Avoid points represented twice list=list(set(points))
        points=zip(XIceShelf,YIceShelf)
        listPoints=list(set(points))
        xPoints=list(zip(*listPoints))[0]
        yPoints=list(zip(*listPoints))[1]

        print(shelf,limits,xPoints)
        
        numPoints=np.size(xPoints)
        freshVal=FWF*factor/numPoints #Share the Gt/yr between all the points which belong to the IS
        varRunoff[:,yPoints,xPoints]=varRunoff[:,yPoints,xPoints]+(freshVal)/area[yPoints,xPoints] #units are KG/m2/s, we divide by the area
            
    return varRunoff


#createCalvingFile       
#---Filling up a grid field with calving rates from a list of ice shelves or sectors
#   Inputs:
#       CalvingVar:Variable to fill up with the calving rates
#       listIceShelves: list of ice shelves or sectors. Must be defined with the same name in data.py
#       MaxPointsPerCalvingP: Max number of grid points represented by a single calving point.
#       xP: List of X grid coordinates of the Coastal points previously extracted
#       yP: List of Y grid coordinates of the Coastal points previously extracted
#       varLon: Longitude field over the same grid where XP and YP are defined
#       varLat: Latitude field over the same grid where XP and YP are defined
#
#   Output:
#       CalvingVar: the grid field filled up
def createCalvingFile(CalvingVar,listIceShelves,MaxPointsPerCalvingP,xP,yP,varLon,varLat):
    for shelf in listIceShelves:
        #Needed data to fill up per Ice Shelf 
        XIceShelf=[] #Container of X grid position of all the segments of an ice shelf
        YIceShelf=[] #Container of Y grid position of all the segments of an ice shelf
        NumOfSegments=0 #Number of segments
        PointsPerSegment=[] #Grid points per segment (dim=NumOfSegments)
        
        limits=globals()[shelf]#Get geospatial limits of the iceshelf i from data.py. It may be made by segments
        
        XIceShelf,YIceShelf,PointsPerSegment,NumOfSegments = getIceShelfSegmentPoints(limits,xP,yP,varLon,varLat)
        XC,YC,WC=CalvPoints(XIceShelf,YIceShelf,PointsPerSegment,NumOfSegments,MaxPointsPerCalvingP)
        if shelf in valuesIB:
            CF=valuesIB[shelf]/0.85 #from GT/yr to km3/yr (density 850)
        elif shelf in  upscalingIB:
            CF=upscalingIB[shelf]/0.85 #from GT/yr to km3/yr (density 850)
        else:
            print('ERROR shelf or sector not found')
            
        CalvingVar[YC,XC]=CalvingVar[YC,XC]+CF*np.array(WC)

    return CalvingVar


#writeDepthFwFluxes       
#---Filling up a grid field with a typical grounding line depth and calving front depth of a list of ice shelves
#   Inputs:
#       varGL:Variable to fill up with the grounding line depth definded in data.py
#       varFront:Variable to fill up with the calving front depth defined in data.py
#       varBathy:Variable with the bathymetry depth
#       listIceShelves: list of ice shelves or sectors. Must be defined with the same name in data.py
#       xP: List of X grid coordinates of the Coastal points previously extracted
#       yP: List of Y grid coordinates of the Coastal points previously extracted
#       varLon: Longitude field over the same grid where XP and YP are defined
#       varLat: Latitude field over the same grid where XP and YP are defined
#
#   Output:
#       varGL,varFront
def writeDepthFwFluxes(varGL,varFront,varBathy,listIceShelves,xP,yP,varLon,varLat):
    for shelf in listIceShelves:
        #Needed data to fill up per Ice Shelf 
        XIceShelf=[] #Container of X grid position of all the segments of an ice shelf
        YIceShelf=[] #Container of Y grid position of all the segments of an ice shelf
        NumOfSegments=0 #Number of segments
        PointsPerSegment=[] #Grid points per segment (dim=NumOfSegments)
        
        limits=globals()[shelf]#Get geospatial limits of the iceshelf i from data.py. It may be made by segments
        
        if shelf in valuesIS:
            Gldepth=valuesGL[shelf]
            Frdepth=valuesFront[shelf]
        elif shelf in  upscaling:
            Gldepth=0.0
            Frdepth=0.0
        else:
            print('ERROR shelf or sector not found')
        
        XIceShelf,YIceShelf,PointsPerSegment,NumOfSegments = getIceShelfSegmentPoints(limits,xP,yP,varLon,varLat)

        # the value to write is the shallower depth between the proposed value and the grid bathymetry
        varGL[:,YIceShelf,XIceShelf]=np.minimum(varBathy[YIceShelf,XIceShelf],(varBathy[YIceShelf,XIceShelf]-varBathy[YIceShelf,XIceShelf]+Gldepth))
        varFront[:,YIceShelf,XIceShelf]=np.minimum(varBathy[YIceShelf,XIceShelf],(varBathy[YIceShelf,XIceShelf]-varBathy[YIceShelf,XIceShelf]+Frdepth))

    return varGL,varFront

#getNearestPoint
#---Return the index corresponding to the list of X and Y values (listPx and lisPy respectively) which is closer to the lon lat values
def getNearestPoint(listPx,listPy,varLon,varLat,lonPoint,latPoint):
    
    testSq=(varLat[listPy,listPx]-latPoint)*(varLat[listPy,listPx]-latPoint)+(varLon[listPy,listPx]-lonPoint)*(varLon[listPy,listPx]-lonPoint)
    indexFin=np.argmin(testSq)
    return indexFin

#getSegment
#---Return a sublist of ranged points between Init and Fin indices of a list of coastal points
def getSegment(listPx,listPy,indexInit,indexFin):
    if indexFin>indexInit:
        listX=listPx[indexInit:indexFin]
        listY=listPy[indexInit:indexFin]
    elif indexFin==indexInit:
        listX=listPx[indexInit:indexFin+1]
        listY=listPy[indexInit:indexFin+1]
    else:
        listX=listPx[indexInit:np.size(listPx)]+listPx[0:indexFin]
        listY=listPy[indexInit:np.size(listPy)]+listPy[0:indexFin]
        
    return list(listX),list(listY)
    
#SegmentPoints
#Extract a list of the ocean points touching land between two coastal points. It follows the shore starting from a given point
def SegmentPoints(varIn,xinit,yinit,xfin,yfin,prevInit):
    [ydim,xdim]=varIn.shape
    
    var = np.concatenate((varIn[:,xinit:xdim],varIn[:,0:xinit+1]),axis=1)

    
    listP=[]
    counter = 1
    
    # prev:Land position relative to the last coastal point known
    #         7  6   5 
    #          \ | /
    #        8 – O – 4
    #          / | \
    #         1  2  3
    #
    
    
    prev=prevInit #By definition xinit and yinit must be an ocean point with land point at yinit-1
    
    x=0
    y=yinit
    
    if xfin<xinit:
        xfin=(xdim-xinit)+xfin-1
    else:
        xfin=xfin-xinit
    yfin=yfin

            
    while (x!=xfin or y!=yfin):
        xfin=xdim-1
        yfin=yinit
        
        if (x+xinit)>=xdim:
            xwrite=x+xinit-xdim
        else:
            xwrite=x+xinit
        
        listP.append((xwrite,y))
        
        
        
        #Case block teating the eigth different cases     
        if prev == 2:
            if var[y,x+1]==0:
                if var[y+1,x]>0:            
                    x=x
                    y=y+1
                    prev=3
                else:
                    x=x
                    y=y
                    prev=6            
            elif var[y,x+1]>0:
                x=x+1
                y=y
                prev=1
                      
        elif prev ==4:
            if var[y+1,x,]==0:
                if var[y,x-1]>0:            
                    x=x-1
                    y=y
                    prev=5
                else:
                    x=x
                    y=y
                    prev=8            
            elif var[y+1,x]>0:
                x=x
                y=y+1
                prev=3
    
        elif prev ==6:
            if var[y,x-1]==0:
                if var[y-1,x]>0:            
                    x=x
                    y=y-1
                    prev=7
                else:
                    x=x
                    y=y
                    prev=2            
            elif var[y,x-1]>0:
                x=x-1
                y=y
                prev=5
    
        elif prev ==8:
            if var[y-1,x]==0:
                if var[y,x+1]>0:            
                    x=x+1
                    y=y
                    prev=1
                else:
                    x=x
                    y=y
                    prev=4            
            elif var[y-1,x]>0:
                x=x
                y=y-1
                prev=7
            
        elif prev ==1:
            if var[y-1,x]>0:
                x=x
                y=y-1
                prev=8          
            elif var[y-1,x]==0:
                if var[y,x+1]>0:
                    x=x+1
                    y=y
                    prev=1
                else:
                    x=x
                    y=y
                    prev=4
            
        elif prev ==3:
            if var[y,x+1]>0:
                x=x+1
                y=y
                prev=2          
            elif var[y,x+1]==0:
                if var[y+1,x]>0:
                    x=x
                    y=y+1
                    prev=3
                else:
                    x=x
                    y=y
                    prev=6
            
        elif prev ==5:
            if var[y+1,x]>0:
                x=x
                y=y+1
                prev=4          
            elif var[y+1,x]==0:
                if var[y,x-1]>0:
                    x=x-1
                    y=y
                    prev=5
                else:
                    x=x
                    y=y
                    prev=8
            
        elif prev ==7:
            if var[y,x-1]>0:
                x=x-1
                y=y
                prev=6          
            elif var[y,x-1]==0:
                if var[y-1,x]>0:
                    x=x
                    y=y-1
                    prev=7
                else:
                    x=x
                    y=y
                    prev=2
                    
        counter=counter+1
        if counter==2000000:
            print('extractCoast felt in a loop')
            break

    
    return listP


#extractCoast
#Extract a list of the ocean points touching land. It follows the shore starting from a given point
def extractCoast(varIn,xinit,yinit):
    xfin=xinit-1 #dummy definition in order to get in into de while loop
    yfin=yinit-1
    [ydim,xdim]=varIn.shape
    
    var = np.concatenate((varIn[:,xinit:xdim],varIn[:,0:xinit+1]),axis=1)

    
    listP=[]
    counter = 1
    
    # prev:Land position relative to the last coastal point known
    #         7  6   5 
    #          \ | /
    #        8 – O – 4
    #          / | \
    #         1  2  3
    #
    
    
    prev=2 #By definition xinit and yinit must be an ocean point with land point at yinit-1
    
    x=0
    y=yinit
    xfin=xdim-1
    yfin=ydim

            
    #while (x!=xfin or y!=yfin):
    while (x<xfin and y<yfin):
        print(x,y)
        #xfin=xdim-1
        #yfin=yinit
        
        if (x+xinit)>=xdim:
            xwrite=x+xinit-xdim
        else:
            xwrite=x+xinit
        
        listP.append((xwrite,y))
        
        
        
        #Case block teating the eigth different cases     
        if prev == 2:
            if var[y,x+1]==0:
                if var[y+1,x]>0:            
                    x=x
                    y=y+1
                    prev=3
                else:
                    x=x
                    y=y
                    prev=6            
            elif var[y,x+1]>0:
                x=x+1
                y=y
                prev=1
                      
        elif prev ==4:
            if var[y+1,x,]==0:
                if var[y,x-1]>0:            
                    x=x-1
                    y=y
                    prev=5
                else:
                    x=x
                    y=y
                    prev=8            
            elif var[y+1,x]>0:
                x=x
                y=y+1
                prev=3
    
        elif prev ==6:
            if var[y,x-1]==0:
                if var[y-1,x]>0:            
                    x=x
                    y=y-1
                    prev=7
                else:
                    x=x
                    y=y
                    prev=2            
            elif var[y,x-1]>0:
                x=x-1
                y=y
                prev=5
    
        elif prev ==8:
            if var[y-1,x]==0:
                if var[y,x+1]>0:            
                    x=x+1
                    y=y
                    prev=1
                else:
                    x=x
                    y=y
                    prev=4            
            elif var[y-1,x]>0:
                x=x
                y=y-1
                prev=7
            
        elif prev ==1:
            if var[y-1,x]>0:
                x=x
                y=y-1
                prev=8          
            elif var[y-1,x]==0:
                if var[y,x+1]>0:
                    x=x+1
                    y=y
                    prev=1
                else:
                    x=x
                    y=y
                    prev=4
            
        elif prev ==3:
            if var[y,x+1]>0:
                x=x+1
                y=y
                prev=2          
            elif var[y,x+1]==0:
                if var[y+1,x]>0:
                    x=x
                    y=y+1
                    prev=3
                else:
                    x=x
                    y=y
                    prev=6
            
        elif prev ==5:
            if var[y+1,x]>0:
                x=x
                y=y+1
                prev=4          
            elif var[y+1,x]==0:
                if var[y,x-1]>0:
                    x=x-1
                    y=y
                    prev=5
                else:
                    x=x
                    y=y
                    prev=8
            
        elif prev ==7:
            if var[y,x-1]>0:
                x=x-1
                y=y
                prev=6          
            elif var[y,x-1]==0:
                if var[y-1,x]>0:
                    x=x
                    y=y-1
                    prev=7
                else:
                    x=x
                    y=y
                    prev=2
                    
        counter=counter+1
        if counter==200000:
            print('extractCoast felt in a loop')
            break

    
    return listP


#CalvPoints

def CalvPoints(XPoints,YPoints,PointsPerSegment,NumOfSegments,MaxPointsPerCalvPoint):
    xCalvPoint=[]
    yCalvPoint=[]
    wCalvPoint=[]
    IceShelfPoints=np.size(XPoints)
    first=0
    for i in np.arange(NumOfSegments):
        SegPoints=PointsPerSegment[i]
        nPoints=1+int(np.floor(SegPoints/MaxPointsPerCalvPoint)) #How much calv points I need in the current segment
        distP=int(np.floor(SegPoints/(nPoints+1))) #Estimate of the distance between calving points at the current segment
        #float value of the relative weight of the calving points relative to the global ice shelf 
        #It considers the relative importance of the given segment into the global ice shelf and the number of calving points contained in the segment
        weight=1.0*SegPoints/IceShelfPoints/nPoints 
        for point in np.arange(nPoints):
            xCalvPoint.append(XPoints[(first+distP)+(point*distP)])
            yCalvPoint.append(YPoints[(first+distP)+(point*distP)])
            wCalvPoint.append(weight)
        first=first+PointsPerSegment[i]
         
    
    return xCalvPoint,yCalvPoint,wCalvPoint


