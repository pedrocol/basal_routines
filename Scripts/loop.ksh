init_year=4
end_year=4
var="salt"
conf1=accessom2-GPC016
conf=`echo "${conf1: -2}"`
conf2=accessom2-GPC009
confm=`echo "${conf2: -2}"`

for var in "salt" "temp"; do
for year in $(seq $init_year $end_year ) ; do
    case "$year" in (0|1|2|3|4) digit=0;esac
    case "$year" in (5|6|7|8|9) digit=1;esac
    case "$year" in (10|11|12|13|14) digit=2;esac
    yearp=$(printf "%02d" $year)
    ncdiff "$conf1"/output00"$digit"/ocean/ocean-3d-"$var"-1-monthly-mean-ym_19"$yearp"_01.nc "$conf2"/output00"$digit"/ocean/ocean-3d-"$var"-1-monthly-mean-ym_19"$yearp"_01.nc diff_"$conf"m"$confm"_"$var"_y"$yearp".nc
    ncra -O "$conf1"/output00"$digit"/ocean/ocean-3d-"$var"-1-monthly-mean-ym_19"$yearp"_01.nc "$conf1"_"$var"_y"$yearp"_avet.nc
    ncra -O -d st_ocean,20,28 "$conf1"_"$var"_y"$yearp"_avet.nc "$conf1"_"$var"_y"$yearp"_avet_avek20-28.nc
done
#Personal diags
#ncrcat diff_"$conf"m"$confm"_"$var"_y* diff_"$conf"m"$confm"_"$var"_y00-"$yearp".nc
#Diags for Mathiot 2017
ncra -O diff_"$conf"m"$confm"_"$var"_y"$yearp".nc diff_"$conf"m"$confm"_"$var"_y"$yearp"_avet.nc
ncra -O -d st_ocean,20,28 diff_"$conf"m"$confm"_"$var"_y"$yearp"_avet.nc diff_"$conf"m"$confm"_"$var"_y"$yearp"_avet_avek20-28.nc
done



#ncks -v "$var" -d st_ocean,1 diff_"$conf"m"$confm"_"$var"_y00-"$yearp".nc diff_"$conf"m"$confm"_"$var"_y00-"$yearp"_k1.nc
#ncks -v "$var" -d st_ocean,19 diff_"$conf"m"$confm"_"$var"_y00-"$yearp".nc diff_"$conf"m"$confm"_"$var"_y00-"$yearp"_k19.nc



#3d values of temp and salt - Daily
#conf="acces-om2-01-GPC001"
conf="acces-om2-01-GAM001"
outini=999
outfin=999
year=2150
#root_data=/home/552/pc5520/access-om2/control/01deg_jra55v13_ryf9091_rerun_for_easterlies/archive 
root_data=/g/data/ik11/outputs/access-om2-01/01deg_jra55v13_ryf9091_rerun_for_easterlies/
var='salt'
root_save=/scratch/e14/pc5520/OUTPUT/$conf/extract/y$year/vert_levels/$var
mkdir -p $root_save/2334/
mkdir -p $root_save/3446/
#t0=$(expr 30 + 59 + 91 + 92 + 1)
#t0=0
t0=273
for digit in $(seq $outini $outfin); do
    number=`ncdump -h $root_data/output$digit/ocean/rregionocean_daily_3d_salt.nc | grep UNLIMITED | sed 's/[^0-9]*//g'`
    #init_t=0
    init_t=46
    end_t=$(expr $number - 1)
    for tim in $(seq $init_t $end_t ); do
        tf=$(expr $t0 + $tim)
        tf=$(seq -f "%03g" $tf $tf)
        echo $tf
        ncks -O -v "$var" -d st_ocean,23,34 -d time,$tim,$tim $root_data/output"$digit"/ocean/rregionocean_daily_3d_"$var".nc $root_save/2334/rregionocean_daily_3d_k2334_"$var"_$tf.nc
        ncks -O -v "$var" -d st_ocean,34,46 -d time,$tim,$tim $root_data/output"$digit"/ocean/rregionocean_daily_3d_"$var".nc $root_save/3446/rregionocean_daily_3d_k3446_"$var"_$tf.nc
    done
    t0=$(expr $t0 + $end_t + 1)
done

#Ice Volume - Monthly
ncrcat -v vicen_m $root_data/output*/ice/OUTPUT/iceh.2150-*.nc $root_save/extract/y2150/ice/ice_vicen_y2150.nc

#Differences
conf1="acces-om2-01-GPC001"
conf2="acces-om2-01-GAM001"
year=2150
var='temp'
for levels in 2334 3446; do
    #Perform differences
    root_data1=/scratch/e14/pc5520/OUTPUT/$conf1/extract/y$year/vert_levels/$var/$levels/
    root_data2=/scratch/e14/pc5520/OUTPUT/$conf2/extract/y$year/vert_levels/$var/$levels/
    root_save=/scratch/e14/pc5520/OUTPUT/$conf1/extract/y$year/diffs/$var/$levels/
    mkdir -p $root_save
    cd $root_data1
    for f in *; do
        echo $f
        ncdiff -O $root_data1/$f $root_data2/$f $root_save/ncdiff_$f
    done
    #Perform vertical average
    root_data=/scratch/e14/pc5520/OUTPUT/$conf1/extract/y$year/diffs/$var/$levels/
    root_save=/scratch/e14/pc5520/OUTPUT/$conf1/extract/y$year/diffs_vertave/$var/$levels/
    mkdir -p $root_save
    cd $root_data
    for f in *; do
	echo $f
        ncwa -O -a st_ocean $f $root_save/vertave_$f
    done
    cd $root_save
    ncrcat -O vertave_* avek_cat_ncdiff_rregionocean_daily_3d_k"$levels"_"$var"_y"$year".nc
    ncea -O  vertave_* avet_avek_ncdiff_rregionocean_daily_3d_k"$levels"_"$var"_y"$year".nc
    rm vertave_ncdiff_*
done


