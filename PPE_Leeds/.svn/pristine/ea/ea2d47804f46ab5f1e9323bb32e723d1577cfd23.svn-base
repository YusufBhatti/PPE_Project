#!/bin/ksh
#-----------------------------------------------------------------------------
FILE=$(ls echam6_userguide.tex)
#-----------------------------------------------------------------------------
ECHAMVERSION=$(grep '\\newcommand{\\echamversion}' ${FILE})
GUIDEVERSION=$(grep '\\newcommand{\\guideversion}' ${FILE})
#-----------------------------------------------------------------------------
ECHAMVERSION=${ECHAMVERSION#'\newcommand{\echamversion}{'}
ECHAMVERSION=${ECHAMVERSION%'}'}
echo 'ECHAMVERSION='$ECHAMVERSION
echo 'Is echam version number correct? (y,j/[n])'
read aw
case "$aw" in
    [yYjJ]*)
        NECHAMVERSION=$ECHAMVERSION ;;
    *)
        echo 'New echam version number:'
        read NECHAMVERSION ;;
esac
#-----------------------------------------------------------------------------
GUIDEVERSION=${GUIDEVERSION#'\newcommand{\guideversion}{'}
GUIDEVERSION=${GUIDEVERSION%'}'}
echo 'GUIDEVERSION='$GUIDEVERSION
echo 'Is guide version number correct? (y,j/[n])'
read aw
case "$aw" in
    [yYjJ]*)
        NGUIDEVERSION=$GUIDEVERSION ;;
    *)
        echo 'New guide version number:'
        read NGUIDEVERSION ;;
esac
#-----------------------------------------------------------------------------
DATE=$(grep '\\newcommand{\\creationdate}' ${FILE})
DATE=${DATE#'\newcommand{\creationdate}{'}
DATE=${DATE%'}'}
echo 'Date of original creation: '$DATE
echo 'Is the date of creation of original creation correct? (y,j/[n])'
read aw
case "$aw" in
    [yYjJ]*)
        NDATE=$DATE ;;
    *)
        echo 'New date of original creation (yyyy-mm-dd):'
        read NDATE ;;
esac
#-----------------------------------------------------------------------------
sed -i -e "s/\\\\newcommand{\\\\echamversion}{${ECHAMVERSION}}/\\\\newcommand{\\\\echamversion}{${NECHAMVERSION}}/" $FILE
sed -i -e "s/\\\\newcommand{\\\\guideversion}{${GUIDEVERSION}}/\\\\newcommand{\\\\guideversion}{${NGUIDEVERSION}}/" $FILE
sed -i -e "s/\\\\newcommand{\\\\creationdate}{${DATE}}/\\\\newcommand{\\\\creationdate}{${NDATE}}/" $FILE
#-----------------------------------------------------------------------------
exit
#-----------------------------------------------------------------------------
