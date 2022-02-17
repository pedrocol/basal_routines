#Ice shelf meltwater in Gt/yr from Depoorter 2013
valuesIS= {'ar':-39, 'ne':-36, 'ais': -39, 'w': -26, 'sha':-76, 'van':-5, 'tot':-64,
'mu':-28, 'por':-18,'ade':-13,'mer':-5,'nin':0,'coo':-3,'ren':-7, 'dry':-5,'ris':-34,
'sul':-28,'lan':-6,'getz':-136,'cd':-78,'thw':-69,'pig':-95,'cosg':-11, 'abb':-86,'ven':-15,
'geo':-144,'wor':-10,'lbc':-18,'fris':-50,'brl':-26,'jf':-24}

#Calving fluxes in Gt/yr from Depoorter 2013
valuesIB= {'ar':44, 'ne':30, 'ais': 50, 'w': 31, 'sha':34, 'van':9, 'tot':28,
'mu':23, 'por':36,'ade':7,'mer':20,'nin':23,'coo':32,'ren':1, 'dry':3,'ris':147,
'sul':3,'lan':11,'getz':56,'cd':18,'thw':62,'pig':50,'cosg':2, 'abb':4,'ven':8,
'geo':17,'wor':2,'lbc':25,'fris':250,'brl':46,'jf':34}

#Typical Grounding Line depth in meters
valuesGL= {'ar':600, 'ne':450, 'ais': 200, 'w': 1000, 'sha':1200, 'van':1000, 'tot':1800,
'mu':2000, 'por':600,'ade':1000,'mer':1000,'nin':1000,'coo':600,'ren':1000, 'dry':900,'ris':1000,
'sul':1000,'lan':600,'getz':700,'cd':1000,'thw':800,'pig':900,'cosg':500, 'abb':500,'ven':400,
'geo':1000,'wor':200,'lbc':700,'fris':1500,'brl':850,'jf':800}

#Typical ice shelf front depth in meters
valuesFront= {'ar':200, 'ne':-200, 'ais': 250, 'w': 200, 'sha':200, 'van':200, 'tot':350,
'mu':350, 'por':300,'ade':200,'mer':250,'nin':300,'coo':250,'ren':200, 'dry':200,'ris':200,
'sul':200,'lan':230,'getz':330,'cd':250,'thw':250,'pig':300,'cosg':280, 'abb':250,'ven':300,
'geo':200,'wor':150,'lbc':200,'fris':200,'brl':20,'jf':200}

#Upscalling values for the ice shelf meltwater in Gt/yr from Depoorter2013
upscaling={'westIndian':-40,'eastIndian':-82,'rossSea':-12,'amundsen':-89,'bellingshausen':-25,'weddell':-13}

#Upscalling values for the calving flux in Gt/yr from Depoorter2013
upscalingIB={'westIndian':49,'eastIndian':93,'rossSea':14,'amundsen':34,'bellingshausen':10,'weddell':16}

#Typical GL depth in meters. Defined arbitrary from BEDMAP2
valuesGL= {'ar':600, 'ne':450, 'ais': 200, 'w': 1000, 'sha':1200, 'van':1000, 'tot':1800,
'mu':2000, 'por':600,'ade':1000,'mer':1000,'nin':1000,'coo':600,'ren':1000, 'dry':900,'ris':1000,
'sul':1000,'lan':600,'getz':700,'cd':1000,'thw':800,'pig':900,'cosg':500, 'abb':500,'ven':400,
'geo':1000,'wor':200,'lbc':700,'fris':1500,'brl':850,'jf':800}

#Typical Calving front depth in meters. Defined arbitrary from BEDMAP2
valuesFront= {'ar':200, 'ne':200, 'ais': 250, 'w': 200, 'sha':200, 'van':200, 'tot':350,
'mu':350, 'por':300,'ade':200,'mer':250,'nin':300,'coo':250,'ren':200, 'dry':200,'ris':200,
'sul':200,'lan':230,'getz':330,'cd':250,'thw':250,'pig':300,'cosg':280, 'abb':250,'ven':300,
'geo':200,'wor':150,'lbc':200,'fris':200,'brl':200,'jf':200}



#Iceshelves of segments,
#!!!! Longitude values defined as -180/180
#every two values is a lon,lat pair definiing a segment
#(ex:8 values means (lon1,lat1,lon2,lat2,lon3,lat3,lon4,lat4) 
#This means that the ice shelf is composed by 2 segments between points 1 and 2, and between 3 and 4)a
# ar, ne, ais Abbot ice shelf, w, sha shackleton,  

lbc=[-61.93,-65.61,-61.74,-69.50]
fris=[-61.56,-74.80,-47.97,-77.90,-44.07,-78.26,-36.11,-78.30]
brl=[-27.15,-76.01,-12.90,-71.98]
jf=[-11.63,-71.24,1.00,-70.00]
ar=[1.5,-70.0,32.17,-69.11]
ne=[34.46,-68.77,38.88,-69.97,47.45,-67.43,48.88,-67.24,49.90,-67.10,50.53,-67.00,56.79,-66.72,59.67,-67.39]
ais=[70.15,-68.49,73.89,-69.73]
w=[81.41,-67.76,88.75,-66.80]
sha=[96.07,-66.48,102.68,-65.90]
van=[108.01,-66.61,110.8,-66.41]
tot=[115.84,-66.53,118.37,-67.00]
mu=[118.37,-67.00,122.60,-66.94]
por=[126.80, -66.56,129.40,-66.80]
ade=[134.18,-66.17,135.39,-66.09]
mer=[144.45,-67.16,145.42,-67.47]
nin=[147.12,-68.05,148.72,-68.39]
coo=[151.74,-68.61,153.74,-68.27]
ren=[161.05,-70.27,162.09,-70.34]
dry=[163.67,-75.10,163.06,-75.49]
ris=[169.24,-77.43,-158.69,-77.88]
sul=[-154.98,-77.10,-149.57,-76.35]
lan=[-142.03,-75.51,-141.22,-75.47]
getz=[-134.82,-74.56,-127.25,-73.72,-125.21,-73.64,-123.44,-73.81,-120.53,-73.83,
     -119.96,-73.83,-119.05,-73.94,-117.31,-74.04,-116.26,-73.94,-115.23,-74.09]
cd=[-113.49,-74.16,-111.90,-74.24,-110.20,-74.66,-109.46,-75.24]
thw=[-108.23,-75.26,-103.76,-75.14]
pig=[-101.90,-75.09,-101.47,-74.30]
cosg=[-101.53,-73.68,-101.68,-73.32]
abb=[-102.82,-72.72,-102.22,-72.12,-100.22,-71.87,-91.76,-72.61,-90.80,-72.69,-89.49,-72.65]
ven=[-88.20,-72.83,-85.76,-73.16]
geo=[-78.77,-73.36,-78.35,-73.01,-77.40,-72.63,-75.75,-72.89,-74.31,-73.05,-72.48,-72.69,-72.86,-72.29,
    -73.58,-72.03,-74.90,-71.55,-71.90,-69.67,-69.16,-70.16,-68.55,-70.10]
wor=[-67.68,-69.38,-67.03,-69.02]

###Sectors for upscalling purpose
westIndian=[1.5,-70.0,96.07,-66.48]
eastIndian=[96.07,-66.48,162.09,-70.34]
rossSea=[162.09,-70.34,-142.03,-75.51]
amundsen=[-142.03,-75.51,-101.68,-73.32]
bellingshausen=[-101.68,-73.32,-67.03,-69.02]
weddell=[-61.93,-65.61,1.5,-70.0]
