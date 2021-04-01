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
:
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
       echo "prova $ID_TCOMP"
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
   NEW_HEADER=$( echo "${LINE_STR_1^^}${LINE_STR_2^^}" ) | sed -e "s/ /;/g" | sed -e "s/;;;/;/g" | sed -e "s/;;/;/g" > $ARG_OUTFILE_TAB

   # Cut the old header from the out tab
   tail -n +5 tmp_$ARG_OUTFILE_TAB | sed -e "s/ /;/g" | sed -e "s/;;;/;/g" | sed -e "s/;;/;/g" >> $ARG_OUTFILE_TAB
   echo "....Done!!"

   # Read the field into the tab and write it into the nc file
   
   # Clean the directory
done 
