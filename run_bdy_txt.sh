#!/bin/bash
#
# by AC Goglio (CMCC)
# annachiara.goglio@cmcc.it
#
# Written: 30/03/2021
#
#set -u
set -e
#set -x 
########################
echo "*********** Tidal bdy 4 EAS System *********"

####################### EXTRACTION FROM TPXO ############################################

echo "Extraction from TPXO model.."

# Load the environment
module load intel19.5/19.5.281 impi19.5/19.5.281 impi19.5/netcdf/C_4.7.2-F_4.5.2_CXX_4.3.1

# Executable to extract single fieds from TPXO along the BDY
EXE_NAME='extract_tidalbdy_txt'

# Tidal components
TIDAL_U_COMPONENTS=('M2' 'S2' 'N2' 'K2' 'K1' 'O1' 'P1' 'Q1')
# WARNING: the order of tidal components MUST be the same given in the f90 code..
#          Remember to check it also in the py script!  

# Bdy infiles template (NEMO bdy files) 
INBDY_TEMPLATE='mfs_bdy%grid%_u2d_*.nc'
OUTBDY_TEMPLATE_TAB_TEMP='bdytide_tab_grid%grid%.txt' 
OUTBDY_TEMPLATE_PRE='mfs_bdytide'
echo "Input bdy file templat: ${INBDY_TEMPLATE}"
echo "Intermediate txt table: ${OUTBDY_TEMPLATE_TAB}"

# Outfields to grids connection
OUTFIELD=("z" "u" "v")
INGRIDS=("T" "U" "V") 
UDM=("cm" "cm/s" "cm/s")

# Loop on out fields/grids (z/T; u/U; v/V )
for ID_FIELDGRID in 0 1 2; do

   echo "I am building tidal bdy field (Re and Imm): ${OUTFIELD[${ID_FIELDGRID}]} [${UDM[${ID_FIELDGRID}]}] on bdy grid: ${INGRIDS[$ID_FIELDGRID]}"

   # Check bdy input file existence
   INBDY=$( ls $( echo $INBDY_TEMPLATE | sed -e "s/%grid%/${INGRIDS[${ID_FIELDGRID}]}/g") )
   if [[ -e $INBDY ]]; then
      echo "I am reading coordinates from ${INBDY}!"

   else
      echo "ERROR: Input file ${INBDY} NOT found!"
   fi

   # Build the outfiles
   echo "I am building the outfile templates.."
   for ID_TCOMP in ${TIDAL_U_COMPONENTS[@]}; do
       OUTBDY=${OUTBDY_TEMPLATE_PRE}_${ID_TCOMP}_grid_${INGRIDS[$ID_FIELDGRID]}.nc
       cp $INBDY ${OUTBDY}
       echo "${OUTBDY} ready"
   done

   # Set the arg
   ARG_ZUV=${OUTFIELD[${ID_FIELDGRID}]}
   ARG_INFILE=${INBDY}
   ARG_OUTFILE_TAB=$( echo $OUTBDY_TEMPLATE_TAB_TEMP | sed -e "s/%grid%/${INGRIDS[${ID_FIELDGRID}]}/g") 

   echo "I am going to run the extractor with the following args: $ARG_ZUV $ARG_INFILE tmp_$ARG_OUTFILE_TAB"

   # Run the extraction script
   echo "I am going to execute the following command: ./${EXE_NAME} $ARG_ZUV $ARG_INFILE tmp_$ARG_OUTFILE_TAB"
   ./${EXE_NAME} $ARG_ZUV $ARG_INFILE tmp_$ARG_OUTFILE_TAB

   # Define and write the new tab header
   LINE_STR_1=$( head -n 3 tmp_$ARG_OUTFILE_TAB | tail -n 1 )
   LINE_STR_2=$( head -n 4 tmp_$ARG_OUTFILE_TAB | tail -n 1 )
   NEW_HEADER=$( echo "${LINE_STR_1^^}${LINE_STR_2^^}" | sed -e "s/ /;/g" | sed -e "s/;;;/;/g" | sed -e "s/;;/;/g" | sed -e "s/;;/;/g" | sed -e "s/;_/_/g" ) 
   echo $NEW_HEADER > $ARG_OUTFILE_TAB

   # Cut the old header from the out tab
   tail -n +5 tmp_$ARG_OUTFILE_TAB | sed -e "s/ /;/g" | sed -e "s/;;;/;/g" | sed -e "s/;;/;/g" | sed -e "s/;;/;/g" >> $ARG_OUTFILE_TAB
   echo "....Done!!"

   # Read the field into the tab and write it into the nc file
   
done 

# Clean the directory
rm tmp_*.txt

####################### TIDAL BDY FIELDS BUILDING  ############################################

echo "Tidal bdy files building.."

# Load the environment
#module purge
#module load anaconda
#source activate mappyenv

B_EXE_NAME='tbdy_writer.py'
python ${B_EXE_NAME}

# Set title and author in the outfiles and remove extra fields
module purge
module load anaconda/3.7 curl/7.70.0 cmake/3.17.3 gams/28.2.0 gcc_9.1.0/9.1.0 gcc_9.1.0/gempack/12.885 gcc_9.1.0/OpenBLAS/0.3.9 gcc_9.1.0/papi/6.0.0 gcc_9.1.0/R/3.6.1 modules mysql/5.7.28 ncl/6.6.2 sqlite/3.32.2 subversion/1.14.0 wgrib/1.8.1.0b impi20.1/19.7.217 impi20.1/esmf/8.0.1-intelmpi-64-g impi20.1/hdf5/1.12.0 impi20.1/hdf5-threadsafe/1.12.0 impi20.1/netcdf/C_4.7.4-F_4.5.3_CXX_4.3.1 impi20.1/netcdf-threadsafe/C_4.7.4-F_4.5.3_CXX_4.3.1 impi20.1/papi/6.0.0 impi20.1/parallel-netcdf/1.12.1 impi20.1/petsc/3.13.2 impi20.1/zoltan/3.8 intel20.1/20.1.217 intel20.1/advisor intel20.1/boost/1.73.0 intel20.1/cdo/1.9.8 intel20.1/cnvgrib/3.1.1 intel20.1/eccodes/2.17.0 intel20.1/esmf/8.0.1-mpiuni-64-g intel20.1/esmf/8.0.1-mpiuni-64-O intel20.1/exactextract/545f0d6 intel20.1/g2lib/3.1.0 intel20.1/gdal/3.1.0 intel20.1/hdf5/1.12.0 intel20.1/hdf5-threadsafe/1.12.0 intel20.1/inspector intel20.1/itac intel20.1/libemos/4.5.9 intel20.1/libemos/4.5.9 intel20.1/magics/3.3.1 intel20.1/nco/4.9.3 intel20.1/ncview/2.1.8 intel20.1/netcdf/C_4.7.4-F_4.5.3_CXX_4.3.1 intel20.1/netcdf-threadsafe/C_4.7.4-F_4.5.3_CXX_4.3.1 intel20.1/proj/7.0.1 intel20.1/R/4.0.2 intel20.1/szip/2.1.1 intel20.1/udunits/2.2.26 intel20.1/valgrind/3.16.0 intel20.1/vtune intel20.1/w3lib/2.0.6 intel20.1/wgrib2/2.0.8

for TOBEMOD in $( ls ${OUTBDY_TEMPLATE_PRE}*.nc ); do
    mv $TOBEMOD ${TOBEMOD}_tmp
    ncatted -O -h -a title,global,m,c,"tidal bdy data from TPXO9-atlas model" ${TOBEMOD}_tmp
    ncatted -O -h -a author,global,m,c,"Anna Chiara Goglio - CMCC" ${TOBEMOD}_tmp
    cdo delete,name=sossheig,vobtcrtx,vobtcrty ${TOBEMOD}_tmp $TOBEMOD
    rm ${TOBEMOD}_tmp
done
