#! /bin/bash
for f in $(find . ! -path ./include/cfortran.h -type f \( -name *.c -o -name *.h \) -print | sort) ; do
  declare -i problem_found
  problem_found=0
  if ! grep '#include *["<]config.h[>"]' "$f" >/dev/null ; then
    echo $f does not contain include of config.h! >&2
    problem_found=1
  fi
  if ! grep '\* *[@\\]file '"$(basename $f)" $f >/dev/null ; then
    echo $f contains no/malformed license/documentation header! >&2
    problem_found=1
  fi
  if ((problem_found)); then
    ${VISUAL-${EDITOR-vi}} "$f"
  fi
done
