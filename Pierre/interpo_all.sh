#!/bin/bash
#OAR -n build_bathy
#OAR -l /nodes=1/core=3,walltime=01:00:00
#OAR --project tipaccs

. ~/bin/load_intelmodule.bash

echo 'build step 1 (make file NEMOBAT ready)'
cp BedMachineAntarctica-2020-07-15.nc BedMachineAntarctica-2020-07-15_NEMOBAT.nc
python BedMachine2NEMOBAT.py -f BedMachineAntarctica-2020-07-15_NEMOBAT.nc -epsg 3031

echo 'build step 2 (interpo)'
./batinterp.exe -f namelist_bathy_BedMachineAntarctica
./batinterp.exe -f namelist_geoid_BedMachineAntarctica
./batinterp.exe -f namelist_isf_BedMachineAntarctica
./batinterp.exe -f namelist_msk_BedMachineAntarctica

echo 'build step 3 (smooth)'
#for nam in `ls namelist_[b-i]*_BedMachineAntarctica`; do
#   ln -sf $nam namelist
#   ./batsmooth
#done

./batsmooth MOM_bathy_step2.nc 2 MOM_bathy_step3.nc Bathymetry 4
./batsmooth MOM_geoid_step2.nc 2 MOM_geoid_step3.nc geoid 4
ncrename -v Bathymetry,geoid MOM_geoid_step3.nc
./batsmooth MOM_isf_step2.nc 2 MOM_isf_step3.nc isf_draft 4
ncrename -v Bathymetry,isf_draft MOM_isf_step3.nc
./batsmooth MOM_msk_step2.nc 2 MOM_msk_step3.nc msk 4
ncrename -v Bathymetry,msk MOM_msk_step3.nc

echo 'build step 4 (gather data in same file)'
cp MOM_bathy_step3.nc MOM_bathy_step4.0.nc
ncks -A -v isf_draft MOM_isf_step3.nc MOM_bathy_step4.0.nc
ncks -A -v msk MOM_msk_step3.nc MOM_bathy_step4.0.nc

echo 'build step 4.1 (set isfd=bathy=0 where grounded and reverse the sign as it is need for BedMachine)'
ncap2 -s 'where ( isf_draft > 0.0 ) msk = 0 ; where ( Bathymetry > 0.0 ) msk = 0 ;' MOM_bathy_step4.0.nc MOM_bathy_step4.1.nc
ncap2 -s 'isf_draft = isf_draft * msk ; Bathymetry = Bathymetry * msk ; isf_draft = - isf_draft ; Bathymetry = - Bathymetry' MOM_bathy_step4.1.nc MOM_bathy_step4.2.nc
ncrename -v Bathymetry,Bathymetry_isf MOM_bathy_step4.2.nc

echo 'interpolation of GEBCO bathymetry into the model grid'
qsub gebco_interp.ksh
mv MOM_bathymetry.nc MOM_bathymetry_0.nc
ncap2 -Os 'where ( Bathymetry > 0.0 ) Bathymetry = 0; Bathymetry = - Bathymetry; where ( Bathymetry > 5808.66 ) Bathymetry = 5808.66' MOM_bathymetry_0.nc MOM_bathymetry_1.nc

echo 'combine bathys'
./combine MOM_bathymetry_1.nc Bathymetry MOM_bathy_step4.2.nc Bathymetry_isf









echo 'build transition mask'
#mask eORCA025 
./cdfmkmask -f eORCA025_bathymetry_b0.1.nc -zoomvar Bathymetry_isf 2000 10000 -fill 625 344 -r -o eORCA025_bathymetry_combine_msk_tmp1.nc
./cdfmkmask -f eORCA025_bathymetry_combine_msk_tmp1.nc -zoomvar tmask 0.0001 1 -fill 448 116 -r -o eORCA025_bathymetry_b0.1_combine_msk_ant2000m.nc
# rm useless time dimension
ncwa -a dummy eORCA025_bathymetry_b0.1_combine_msk_ant2000m.nc -o eORCA025_bathymetry_b0.1_combine_msk_ant2000m_nodummy.nc
cp eORCA025_bathymetry_b0.0.nc eORCA025_bathymetry_b0.1_combine.nc
ncks -A -v tmask eORCA025_bathymetry_b0.1_combine_msk_ant2000m_nodummy.nc eORCA025_bathymetry_b0.1_combine.nc

#mask BedMachine
# get mask 3000m (open ocean)
./cdfmkmask -f eORCA025_bathymetry_b0.1.nc -zoomvar Bathymetry_isf 3000 10000 -fill 625 344 -r -o eORCA025_bathymetry_combine_msk_tmp2.nc
# get Antarctic mask (0-3000m)
./cdfmkmask -f eORCA025_bathymetry_combine_msk_tmp2.nc -zoomvar tmask 0.0001 1 -fill 448 116 -o eORCA025_bathymetry_b0.1_combine_msk_tmp3.nc
# start building mask Antarctic shelf (limited to 1000m isobath on shelf break)
cp eORCA025_bathymetry_b0.1_combine_msk_tmp3.nc eORCA025_bathy_tmp1.nc
# add step4.2 bathy
ncks -A -v Bathymetry_isf eORCA025_bathy_step4.2.nc eORCA025_bathy_tmp1.nc
# mask bathy between where bathy > 3000m
ncap2 -s 'Bathymetry_isf = Bathymetry_isf * tmask' eORCA025_bathy_tmp1.nc eORCA025_bathy_tmp2.nc
# extract 1000m 3000m band on continental slope
./cdfmkmask -f eORCA025_bathy_tmp2.nc -zoomvar Bathymetry_isf 1000 3000 -fill 966 258 -r -o eORCA025_bathy_tmp3.nc
# extract antarctica
./cdfmkmask -f eORCA025_bathy_tmp3.nc -zoomvar tmask 0.5 1.5 -fill 5 5 -o eORCA025_bathy_step4.2_msk_ant1000m.nc
# add mask in eORCA025_bathy_step4.2.nc file => eORCA025_bathy_step4.3.nc
cp eORCA025_bathy_step4.2.nc eORCA025_bathy_step4.3.nc
# rm useless time dimension
ncwa -a dummy eORCA025_bathy_step4.2_msk_ant1000m.nc -o eORCA025_bathy_step4.2_msk_ant1000m_nodummy.nc
ncks -A -v tmask eORCA025_bathy_step4.2_msk_ant1000m_nodummy.nc eORCA025_bathy_step4.3.nc

echo 'combine both files'
ln -sf namelist_bathy_BedMachineAntarctica namelist
./combine
ncks -A -v isf_draft eORCA025_bathy_step4.3.nc eORCA025_bathy_step5.0.nc

echo 'add ice shelf draft'
ncatted -a valid_min,,m,f,0.0 eORCA025_bathy_step5.0.nc
ncap2 -O -s 'Bathymetry = Bathymetry_isf ; where( isf_draft > 0 ) Bathymetry = 0.0 ;' eORCA025_bathy_step5.0.nc eORCA025_bathy_step5.1.nc

echo 'rm weird point in Bathymetry in all the data set (single ocean point close to the claving front)'
rm msk_holeintheisf_nodummy.nc
./cdfisf_2dpoolchk -d eORCA025_bathy_step5.1.nc -v isf_draft -nc4 -o msk_holeintheisf.nc -b eORCA025_bathy_step5.1.nc
ncwa -a dummy msk_holeintheisf.nc msk_holeintheisf_nodummy.nc
ncks -A -v tmask_pool2d msk_holeintheisf_nodummy.nc eORCA025_bathy_step5.1.nc
ncap2 -O -s 'Bathymetry = Bathymetry * (1-tmask_pool2d) ; Bathymetry_isf = Bathymetry_isf * (1-tmask_pool2d) ; isf_draft = isf_draft * (1-tmask_pool2d) ;' eORCA025_bathy_step5.1.nc eORCA025_bathy_step5.2.nc

echo 'build final file'
nccopy -4 -d 1 -c x/200,y/200 -v nav_lat,nav_lon,isf_draft,Bathymetry_isf,Bathymetry eORCA025_bathy_step5.2.nc eORCA025_bathy_step5.3.nc
ncks -x -v tmask_pool2d eORCA025_bathy_step5.3.nc eORCA025_bathymetry_b0.2.nc

