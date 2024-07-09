#!/bin/ksh
###################             USAGE                  ###################
# S. Rast, MPI-Met, 07/20/2009
# multiplot_latex.sh <output texfile name> <list of figures (ps or eps)>
# the texfiles name has to be given without extension .tex
# example: multiplot_latex.sh stations stn*.ps
##########################################################################
# path of location of vgr_pt file (if it is in the directory from which
# you start the script, it can be empty) 
VPATH=
#number of figures per page (1-9 allowed)
NFIGURES=6
#width of each individual figure in cm
WIDTH=8.0
#angle in degrees of the figure (for some postscript files 270, 
#for eps files: 0, 90, 270 are the most commonly needed angles 
ANGLE=0
#instruction section
TEXFILE=$1
GGFILES="$*"
GFILES=${GGFILES#$TEXFILE}
###over###if [ -f ${TEXFILE}.aux -o -f ${TEXFILE}.dvi -o -f ${TEXFILE}.log -o -f ${TEXFILE}.pdf -o -f ${TEXFILE}.ps -o -f ${TEXFILE}.tex ]; then
###over###echo 'one of the files '${TEXFILE}'.{aux,dvi,log,pdf,ps,tex} exists. ###over###Overwrite (y/n)?'
###over###read AW
###over###if [ $AW != 'y' ]; then
###over###exit
###over###fi
###over###fi
#echo $TEXFILE
#echo $GFILES
set -A PS ""
for GF in $GFILES; do
set -A PS ${PS[*]} $GF
done
#VPATH=${VPATH:=.}
VPATH=/pool/data/ECHAM5/post_processing/scripts/gmt

cat > ${TEXFILE}.tex <<EOF
\input{${VPATH}/vgr_pt.tex}
\begin{document}

EOF
ii=0
NPS=$(( ${#PS[*]} / ${NFIGURES} ))
NREST=$(( ${#PS[*]} % ${NFIGURES} ))
#echo ${#PS[*]} $NREST $NPS
id=0
ips=1
while [ $ips -le $NPS ]; do
#echo $ips $NPS
# choice of number of graphs on one page
if [ $NFIGURES -eq 1 ]; then
   PUT='\pctone'
fi 
if [ $NFIGURES -eq 2 ]; then
   PUT='\pcttwo'
fi 
if [ $NFIGURES -eq 3 ]; then
   PUT='\pctthreebot'
fi 
if [ $NFIGURES -eq 4 ]; then
   PUT='\pctfour'
fi 
if [ $NFIGURES -eq 5 ]; then
   PUT='\pctfivebot'
fi 
if [ $NFIGURES -eq 6 ]; then
   PUT='\pctsix'
fi 
if [ $NFIGURES -eq 7 ]; then
   PUT='\pctsevenbot'
fi 
if [ $NFIGURES -eq 8 ]; then
   PUT='\pcteight'
fi 
if [ $NFIGURES -eq 9 ]; then
   PUT='\pctninet'
fi 
#echo $NFIGURES $PUT
iid=1
while [ $iid -le ${NFIGURES} ]; do
PUT=$PUT"{\pfaw{${PS[$id]}}{${ANGLE}}{${WIDTH}}}"
id=$(( id + 1 ))
iid=$(( iid + 1 ))
done
cat >> ${TEXFILE}.tex <<EOF

$PUT
\newpage
EOF
ips=$(( ips + 1 ))
done
if [ $NREST -gt 0 ]; then
if [ $NREST -eq 1 ]; then
   PUT='\pctone'
fi 
if [ $NREST -eq 2 ]; then
   PUT='\pcttwo'
fi 
if [ $NREST -eq 3 ]; then
   PUT='\pctthreebot'
fi 
if [ $NREST -eq 4 ]; then
   PUT='\pctfour'
fi 
if [ $NREST -eq 5 ]; then
   PUT='\pctfivebot'
fi 
if [ $NREST -eq 6 ]; then
   PUT='\pctsix'
fi 
if [ $NREST -eq 7 ]; then
   PUT='\pctsevenbot'
fi 
if [ $NREST -eq 8 ]; then
   PUT='\pcteight'
fi 
if [ $NREST -eq 9 ]; then
   PUT='\pctninet'
fi 
iid=1
while [ $iid -le $NREST ]; do
PUT=$PUT"{\pfaw{${PS[$id]}}{${ANGLE}}{${WIDTH}}}"
id=$(( id + 1 ))
iid=$(( iid + 1 ))
done

cat >> ${TEXFILE}.tex <<EOF

$PUT
\newpage
EOF
fi
cat >> ${TEXFILE}.tex <<EOF
\end{document}
EOF
latex ${TEXFILE}.tex
dvips -Ppdf ${TEXFILE}.dvi

#make pdf-file
#ps2pdf ${TEXFILE}.ps
#show plot
#gv ${TEXFILE}.pdf &
exit
