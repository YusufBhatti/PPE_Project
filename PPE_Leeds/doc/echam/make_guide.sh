#!/bin/ksh
#-----------------------------------------------------------------------------
FILE=$(ls echam6_userguide.tex)
#-----------------------------------------------------------------------------
./update_version.sh
#-----------------------------------------------------------------------------
pdflatex ${FILE%.tex}
makeindex ${FILE%.tex}
bibtex ${FILE%.tex}
pdflatex ${FILE%.tex}
pdflatex ${FILE%.tex}
acroread ${FILE%.tex}.pdf &
echo 'print this file? (y,j/[n])'
read aw
case "$aw" in
    [yYjJ]*)
        echo 'Printer name for pdf file:' 
        read aww
        lpr -P "$aww" ${FILE%.tex}.pdf ;;
esac
exit
#-----------------------------------------------------------------------------
