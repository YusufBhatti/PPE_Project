#define _XOPEN_SOURCE 700
#define _GNU_SOURCE // needed for getline(3) on some systems it seems
#include <ctype.h>
#include <errno.h>
#include <regex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

#ifndef HAVE_GETLINE
#include "getline.c"
#endif

//#include "config.h"
#define VERSION "1.7.0"
typedef struct
{
  size_t naline;
  char *fname;
  char *aline[99];
  char *text;
}
Docu;

// Example: ./make_fint -d ../doc/pio/ cdipio.h

static const char default_doc_root[] = "../doc";
static struct {
  const char *name;
  size_t len;
} fname_list[] = {
  { "c_quick_ref.txt", 0 },
  { "f_quick_ref.txt", 0 },
  { "tex/c_quick_ref.tex", 0 },
  { "tex/f_quick_ref.tex", 0 },
};
enum {
  NAME_C_QUICK_REF,
  NAME_F_QUICK_REF,
  NAME_C_QUICK_REF_TEX,
  NAME_F_QUICK_REF_TEX,
  fname_list_size = sizeof(fname_list)/sizeof(fname_list[0]),
};


static Docu cdoc[9999], fdoc[9999];
static size_t ncdoc = 0, nfdoc = 0;
static int debug = 0, verbose = 0;

static int doccmp(const void *s1, const void *s2)
{
  Docu *x = (Docu *) s1;
  Docu *y = (Docu *) s2;

  return (strcmp(x->fname, y->fname));
}

static void doctotex(FILE *fp, Docu *doc, size_t ndoc)
{
  size_t i, k;

  for ( i = 0; i < ndoc; i++ )
    {
      fprintf(fp, "\\section*{\\tt \\htmlref{%s}{%s}}\n\n", doc[i].fname, doc[i].fname);
      fprintf(fp, "\\begin{verbatim}\n");
      for ( k = 0; k < doc[i].naline; k++ )
	fprintf(fp, "    %s\n", doc[i].aline[k]);
      fprintf(fp, "\\end{verbatim}\n");
      fprintf(fp, "\n%s.\n\n\n", doc[i].text);
    }
}

static void doctotxt(FILE *fp, Docu *doc, size_t ndoc)
{
  size_t i, k;

  for ( i = 0; i < ndoc; i++ )
    {
      fprintf(fp, "%s\n\n", doc[i].fname);
      for ( k = 0; k < doc[i].naline; k++ )
	fprintf(fp, "    %s\n", doc[i].aline[k]);
      fprintf(fp, "\n  %s.\n\n", doc[i].text);
    }
}

enum cftype {ISVOID, ISCONSTSTRING, ISINT, ISREAL, ISDOUBLE, ISMPI_COMM,
             ISXT_IDXLIST, ISCHOICE, ISINTP, ISFLOATV, ISFLOATVV,
             ISDOUBLEV, ISDOUBLEVV, ISINTV, ISINTVV, ISREALP,
             ISDOUBLEP, ISCBUF, ISUUID, ISUCHAR, ISSTRING, ISSTRINGP,
             VOIDFUNCVOID,
             NUM_KNOWN_ARG_TYPES};

static inline int
isArrayArgType(int argType);

enum conversionType { CONV_ARG, CONV_RET };


typedef int (*cfConversionEmitter)(FILE *outfp, const char *argName,
                                   size_t argNameLen, enum conversionType part);
typedef int (*cfPrologueEmitter)(FILE *outfp, size_t argNum);


static int cfMPICommConvert(FILE *outfp, const char *argName,
                            size_t argNameLen, enum conversionType part);

static int cfXtIdxlistConvert(FILE *outfp, const char *argName,
                            size_t argNameLen, enum conversionType part);

static int cfVoidFuncPrologue(FILE *outfp, size_t argNum);

struct symbol {
  const char *f77name, *cfint, *cfmt, *parseRE;
  /* pair of parentheses which matches the argument name */
  size_t nameMatch;
  int needsExtraWrapper, needsPrologue;
  cfConversionEmitter convert;
  const char *convcfmt;
  cfPrologueEmitter prologue;
  regex_t preg;
};

/* C symbol names */
#define SYMRE "([A-Za-z_][A-Za-z_0-9]*)"
/* white-space */
#define WS "[[:blank:]\n]"
#define NWS "[^[:blank:]\n]"
/* Note: size of this table must match the cftype enum */
static struct symbol funArgSym[]
  = { { "",                "",        "void",
        "^"WS"*void"WS"*)", 0, 0, 0 },
      { "CHARACTER(80)",    "STRING",  "char *%.*s",
        "^"WS"*const"WS"+char"WS"+\\*"SYMRE WS"*\\(", 1, 0, 0 },
      { "INTEGER",         "INT",     "int %.*s",
        "^"WS"*(const"WS"+)?int("WS"+"SYMRE")?"WS"*[,\\)]", 3, 0, 0 },
      { "REAL",            "FLOAT",   "float %.*s",
        "^"WS"*(const"WS"+)?float"WS"+"SYMRE"?"WS"*[,\\)]", 2, 0, 0 },
      { "DOUBLEPRECISION", "DOUBLE",  "double %.*s",
        "^"WS"*(const"WS"+)?double"WS"+"SYMRE"?"WS"*[,\\)]", 2, 0, 0 },
      { "INTEGER",         "INT", "MPI_Comm %.*s",
        "^"WS"*MPI_Comm"WS"+"SYMRE"?"WS"*[,\\)]", 1, 1, 0,
        cfMPICommConvert, "int %.*s" },
      { "TYPE(XT_IDXLIST)", "PVOID", "Xt_idxlist %.*s",
        "^"WS"*Xt_idxlist"WS"+"SYMRE"?"WS"*[,\\)]", 1, 1, 0,
        cfXtIdxlistConvert, "void *%.*s" },
      { "CHOICE", "PVOID", "const void *%.*s",
        "^"WS"*const"WS"+void"WS"*\\*"WS"*"SYMRE"?"WS"*[,\\)]", 1, 0, 0 },
      { "INTEGER",         "PINT",    "int *%.*s",
        "^"WS"*(const"WS"+)?int"WS"+\\*"SYMRE"?"WS"*[,\\)]", 2, 0, 0 },
      { "REAL",            "FLOATV",  "float %.*s[]",
        "^"WS"*(const"WS"+)?float("WS"+"SYMRE")?"WS"*\\[[^]]*\\]"
        WS"*[,\\)]", 3, 0, 0 },
      { "REAL",            "FLOATVV", "float %.*s[][]",
        "^"WS"*(const"WS"+)?float("WS"+"SYMRE")?"WS"*\\[[^]]*\\]"
        WS"*\\[[^]]*\\]"WS"*[,\\)]", 3, 0, 0 },
      { "DOUBLEPRECISION", "DOUBLEV",  "double %.*s[]",
        "^"WS"*(const"WS"+)?double("WS"+"SYMRE")?"WS"*\\[[^]]*\\]"
        WS"*[,\\)]", 3, 0, 0 },
      { "DOUBLEPRECISION", "DOUBLEVV", "double %.*s[][]",
        "^"WS"*(const"WS"+)?double("WS"+"SYMRE")?"WS"*\\[[^]]*\\]"
        WS"*\\[[^]]*\\]"WS"*[,\\)]", 3, 0, 0 },
      { "INTEGER",         "INTV",    "int  %.*s[]",
        "^"WS"*(const"WS"+)?int("WS"+"SYMRE")?"WS"*\\[[^]]*\\]"
        WS"*[,\\)]", 3, 0, 0 },
      { "INTEGER",         "INTVV",    "int %.*s[][]",
        "^"WS"*(const"WS"+)?int("WS"+"SYMRE")?"WS"*\\[[^]]*\\]"
        WS"*\\[[^]]*\\]"WS"*[,\\)]", 3, 0, 0 },
      { "REAL",            "PFLOAT",  "float *%.*s",
        "^"WS"*(const"WS"+)?float"WS"+\\*"SYMRE"?"WS"*[,\\)]", 2, 0, 0 },
      { "DOUBLEPRECISION", "PDOUBLE", "double *%.*s",
        "^"WS"*(const"WS"+)?double"WS"+\\*"SYMRE"?"WS"*[,\\)]", 2, 0, 0 },
      { "CHARACTER*(*)",   "PPSTRING",    "char *%.*s",
        "^"WS"*(const"WS"+)?char"WS"+\\*""([A-Za-z_][A-Za-z_0-9]*_cbuf)"
        WS"*[,\\)]", 2, 0, 0 },
      { "INTEGER*1(16)",   "PVOID",    "unsigned char %.*s[16]",
        "^"WS"*(const"WS"+)?unsigned"WS"+char"WS"+"SYMRE"?\\[(16|CDI_UUID_SIZE)\\]"WS"*[,\\)]", 2, 0, 0 },
      { "INTEGER*1(*)",   "PVOID",    "unsigned char *%.*s",
        "^"WS"*(const"WS"+)?unsigned"WS"+char"WS"+\\*"SYMRE"?"WS"*[,\\)]", 2, 0, 0 },
      { "CHARACTER*(*)",   "STRING",  "char *%.*s",
        "^"WS"*const"WS"+char"WS"+\\*"WS"*"SYMRE"?"WS"*[,\\)]", 1, 0, 0 },
      { "CHARACTER*(*)",   "PSTRING", "char *%.*s",
        "^"WS"*char"WS"+\\*"SYMRE"?"WS"*[,\\)]", 1, 0, 0 },
      { "PROCEDURE", "ROUTINE", "void (*%.*s)(void)",
        "^"WS"*void"WS"*\\("WS"*\\*"WS"*"SYMRE"?"WS"*\\)"
        WS"*\\("WS"*void"WS"*\\)"WS"*[,\\)]", 1, 0, 1,
        NULL, NULL, cfVoidFuncPrologue },
};

static struct symbol funRet[] = {
  { "",                "",        "void %.*s",
    "void"WS"+"SYMRE WS"*\\(", 1, 0, 0 },
  { "CHARACTER",       "STRING",  "char *%.*s",
    "char"WS"+\\*"WS"*"SYMRE WS"*\\(", 1, 0, 0 },
  { "INTEGER",         "INT",     "int %.*s",
    "(const"WS"+)?int"WS"+"SYMRE WS"*\\(", 2, 0, 0 },
  { "REAL",            "FLOAT",   "float %.*s",
    "(const"WS"+)?float"WS"+"SYMRE WS"*\\(", 2, 0, 0 },
  { "DOUBLEPRECISION", "DOUBLE",  "double %.*s",
    "(const"WS"+)?double"WS"+"SYMRE WS"*\\(", 2, 0, 0 },
  { "INTEGER",         "INT",     "MPI_Comm %.*s",
    "MPI_Comm"WS"+"SYMRE WS"*\\(", 1, 0, 0, cfMPICommConvert, "int %.*s" },
};

enum { NUM_RET_TYPES = sizeof (funRet) / sizeof (funRet[0]) };
enum decl { UNKNOWN_DECL, FUNC_DECL, PARAM_DECL };

enum {
  MAX_FUNC_ARGS = 200,
  MAX_FUNC_NAME_LEN = 127,
};

static inline size_t
compress_whitespace(size_t len, char str[]);

static int
reCompile(regex_t *restrict RE, const char *restrict REstring,
          char * restrict *restrict lineBuf, size_t * restrict lineBufSize);
static size_t
symRegexCompile(size_t numSyms, struct symbol symList[],
                char **line, size_t *lineBufSize);

static void
build_header_name(const char *fname, char *cppMacro);

static int detectComment(char **line_, ssize_t *lineLen, size_t *lineBufSize,
                         size_t maxMatch, regmatch_t reMatch[],
                         char *xname, size_t *xnameLen,
                         char *xdes,
                         FILE *fpin, FILE *fpinc, FILE *fpint);

static regex_t commentStartRE, commentEndRE, commentRE, docCommentRE;

static void fortran_interface(char *fname, char *fnameinc, char *fnameint,
                              const char *doc_root)
{
  FILE *fpin, *fpinc, *fpint;
  FILE *fp;
  char *line = NULL, *pline;
  size_t lineBufSize = 0;
  char sname[128], *parname;
  char xname[128], xdes[128];
  xname[0] = 0;
  size_t xnameLen = 0;
  int parvalue;
  enum cftype functype;
  int lineno = 0;

  char funcname[MAX_FUNC_NAME_LEN];
  regmatch_t funcargfull[MAX_FUNC_ARGS];
  regmatch_t funcargname[MAX_FUNC_ARGS];
  int  funcargtype[MAX_FUNC_ARGS];
  /* char *strsort[99999]; */
  char timestr[30];
  time_t date_and_time_in_sec;
  struct tm *date_and_time;
  regmatch_t *reMatch = NULL;
  size_t maxMatch = 0;

  date_and_time_in_sec = time(NULL);
  timestr[0] = 0;

  if ( date_and_time_in_sec != -1 )
    {
      date_and_time = localtime(&date_and_time_in_sec);
      (void) strftime(timestr, sizeof(timestr), "%B %Y", date_and_time);
    }

  fpin = fopen(fname, "r");
  if ( fpin == NULL ) { perror(fname); return; }

  fpinc = fopen(fnameinc, "w");
  if ( fpinc == NULL ) { perror(fnameinc); return; }

  fpint = fopen(fnameint, "w");
  if ( fpint == NULL ) { perror(fnameint); return; }

  /* complete symbol table data */
  {
    maxMatch = symRegexCompile(NUM_KNOWN_ARG_TYPES, funArgSym,
                             &line, &lineBufSize);
    size_t maxFunMatch = symRegexCompile(NUM_RET_TYPES, funRet,
                                         &line, &lineBufSize);
    if (maxFunMatch > maxMatch)
      maxMatch = maxFunMatch;
  }
  ++maxMatch;
  reMatch = (regmatch_t *)malloc((size_t)maxMatch * sizeof (reMatch[0]));
  /* compile comment start regular expression */
  {
    static const char commentStartREString[] = "^"WS"*/\\*"WS"*(.*"NWS")"WS"*";
    if (reCompile(&commentStartRE, commentStartREString, &line, &lineBufSize))
      exit(EXIT_FAILURE);
  }
  /* compile comment end regular expression */
  {
    static const char commentEndREString[] = "\\*/";
    if (reCompile(&commentEndRE, commentEndREString, &line, &lineBufSize))
      exit(EXIT_FAILURE);
  }
  /* compile complete comment regular expression */
  {
    static const char commentREString[] = "^"WS"*/\\*"WS"*(.*"NWS")"WS"*\\*/";
    if (reCompile(&commentRE, commentREString, &line, &lineBufSize))
      exit(EXIT_FAILURE);
  }
  /* compile documentation comment regular expression */
  {
    static const char docCommentREString[] = "^"WS"*/\\*"WS"*"SYMRE":"
      WS"*("NWS".*"NWS")"WS"*\\*/";
    if (reCompile(&docCommentRE, docCommentREString, &line, &lineBufSize))
      exit(EXIT_FAILURE);
  }
  regex_t cppCondRE;
  {
    static const char cppCondREString[]
      = "^"WS"*#"WS"*((ifn?def)"WS"+\\(?"SYMRE"\\)?|endif)"WS"*(/\\*[^*]*\\*/|//.*)?";
    if (reCompile(&cppCondRE, cppCondREString, &line, &lineBufSize))
      exit(EXIT_FAILURE);
  }
  regex_t cppElseRE;
  {
    static const char cppElseREString[]
      = "^"WS"*#"WS"*else"WS"*(/\\*[^*]*\\*/|//.*)?";
    if (reCompile(&cppElseRE, cppElseREString, &line, &lineBufSize))
      exit(EXIT_FAILURE);
  }
  regex_t emptyStringRE;
  {
    static const char emptyStringREString[] = "^"WS"*";
    if (reCompile(&emptyStringRE, emptyStringREString, &line, &lineBufSize))
      exit(EXIT_FAILURE);
  }
  /* fortran include */

  fprintf(fpinc, "! This file was automatically generated, don't edit!\n");
  fprintf(fpinc, "!\n");
  fprintf(fpinc, "! Fortran interface for CDI library version %s\n", VERSION);
  fprintf(fpinc, "!\n");
  fprintf(fpinc, "! Author:\n");
  fprintf(fpinc, "! -------\n");
  fprintf(fpinc, "! Uwe Schulzweida, MPI-MET, Hamburg,   %s\n", timestr);
  fprintf(fpinc, "!\n\n");

  /* fortran interface */

  fprintf(fpint, "/* Automatically generated by make_fint.c, don't edit! */\n");
  fprintf(fpint, "\n");
  fprintf(fpint, "#if defined (HAVE_CONFIG_H)\n");
  fprintf(fpint, "#  include \"config.h\"\n");
  fprintf(fpint, "#endif\n");
  fprintf(fpint, "\n");
  {
    char *cppMacro = (char*) malloc(strlen(fname) + 2);
    build_header_name(fname, cppMacro);
    fprintf(fpint, "#if ! defined (%s)\n"
            "#  include \"%s\"\n"
            "#endif\n"
            "\n", cppMacro, fname);
  }
  fprintf(fpint, "#if defined (HAVE_CF_INTERFACE)\n");
  fprintf(fpint, "\n");
  fprintf(fpint, "#if ! defined (__CFORTRAN_LOADED)\n");
  fprintf(fpint, "#  include \"cfortran.h\"\n");
  fprintf(fpint, "#endif\n");
  fprintf(fpint, "\n");

  ssize_t lineLen;
  while ((lineLen = getline(&line, &lineBufSize, fpin)) >= 0)
    {
      static const char cplusplus_macro[] = "__cplusplus";
      lineno++;
      functype = ISVOID;
      size_t funcargc = 0;
      pline = line;
      int needsExtraWrapper = 0, needsPrologue = 0;
      size_t funcnameLen;
      enum decl declType = UNKNOWN_DECL;
      do {
        for (int retType = 0; retType < NUM_RET_TYPES; ++retType)
          if (!regexec(&funRet[retType].preg, pline, maxMatch, reMatch, 0))
            {
              functype = retType;
              declType = FUNC_DECL;
              needsExtraWrapper
                = needsExtraWrapper || funRet[retType].needsExtraWrapper;
              break;
            }
        if (declType == UNKNOWN_DECL)
          break;
        regmatch_t *nameMatch = reMatch + funRet[functype].nameMatch;
        if (debug)
          printf("Found: %.*s\n",
                 (int) (nameMatch->rm_eo - nameMatch->rm_so),
                 pline + nameMatch->rm_so);
        ssize_t funNameLast = reMatch[0].rm_eo - 1;
        ssize_t nameLen = nameMatch->rm_eo - nameMatch->rm_so;
        funcnameLen = (size_t)nameLen;
        if ( pline[funNameLast] != '(' )
          {
            printf("%s\n>(< not found!", line);
            return;
          }
        memcpy(funcname, pline + nameMatch->rm_so, funcnameLen);
        funcname[funcnameLen] = 0;
        pline += reMatch[0].rm_eo;
      } while (0);
      int cppSwitchLen, cppSymLen;

      if (declType == FUNC_DECL)
        {
	  funcargname[funcargc].rm_so = (regoff_t)(pline - line);
          {
            ssize_t i = 0;
            size_t innerParens = 0;
            do {
              ssize_t restLen = lineLen - (ssize_t)(pline - line);
              for (; i < restLen; i++ )
                {
                  switch (pline[i])
                    {
                    case ',':
                      if (!innerParens)
                        {
                          funcargc++;
                          funcargname[funcargc].rm_so
                            = (regoff_t)(pline - line + i + 1);
                        }
                      break;
                    case '(':
                      ++innerParens;
                      break;
                    case ')':
                      if (!innerParens)
                        {
                          funcargc++;
                          funcargname[funcargc].rm_so
                            = (regoff_t)(pline - line + i + 1);
                          goto endOfArgSearch;
                        }
                      else
                        --innerParens;
                      break;
                    }
                }
              endOfArgSearch:
              if (i < restLen)
                break;
              char *lineExtension = NULL;
              size_t extSize = 0, plineOff = (size_t)(pline - line);
              ssize_t extLen;
              if ((extLen = getline(&lineExtension, &extSize, fpin)) <= 0)
                break;
              if ((size_t)(lineLen + extLen) >= lineBufSize)
                if (!(line = (char*) realloc(line, (size_t)(lineLen + extLen + 1))))
                  exit(EXIT_FAILURE);
              memcpy(line + lineLen, lineExtension, (size_t)extLen + 1);
              lineLen += extLen;
              pline = line + plineOff;
            } while (1);
          }

	  /*  printf("funcargc = %d\n", funcargc);*/
            /* test if argument list is actually empty */
          if (funcargc == 1
              && !regexec(&emptyStringRE, line + funcargname[0].rm_so, 1,
                          reMatch, 0)
              && (funcargname[0].rm_so + reMatch[0].rm_eo
                  == funcargname[funcargc].rm_so - 1))
            {
              funcargc = 0;
            }
          {
            size_t i;
            for (i = 0; i < funcargc; ++i )
              {
                pline = line + funcargname[i].rm_so;
                int argtype;
                regoff_t argStart = (regoff_t)(pline - line);
                for (argtype = ISVOID;
                     argtype < NUM_KNOWN_ARG_TYPES;
                     ++argtype)
                  if (!regexec(&funArgSym[argtype].preg, pline, maxMatch,
                               reMatch, 0))
                    {
                      funcargtype[i] = argtype;
                      funcargfull[i].rm_so = reMatch[0].rm_so + argStart;
                      funcargfull[i].rm_eo = reMatch[0].rm_eo + argStart;
                      regmatch_t *nameMatch = reMatch + funArgSym[argtype].nameMatch;
                      funcargname[i].rm_so = nameMatch->rm_so + argStart;
                      funcargname[i].rm_eo = nameMatch->rm_eo + argStart;
                      needsExtraWrapper
                        = needsExtraWrapper || funArgSym[argtype].needsExtraWrapper;
                      needsPrologue = needsPrologue
                        || funArgSym[argtype].needsPrologue;
                      break;
                    }
                if (argtype == NUM_KNOWN_ARG_TYPES)
                  {
                    printf("%s not implemented\n", line + funcargname[i].rm_so);
                    break;
                  }
              }
            if ( i != funcargc )
              {
                printf("problem parsing line: %s\n", line);
                continue;
              }
          }

	  strcpy(sname, funcname);

	  /* fortran include */

	  if ( functype == ISVOID )
	    fprintf(fpinc, "!     %-16s", "");
	  else
	    fprintf(fpinc, "      %-16s", funArgSym[functype].f77name);

          fprintf(fpinc, "%s", sname);
	  fprintf(fpinc, "\n");
	  if ( (funcargc == 1 && funcargtype[0] == ISVOID) ) funcargc = 0;
	  for (size_t i = 0; i < funcargc; i++ )
	    {
	      if ( i == 0 )
		fprintf(fpinc, "!%36s(", "");
	      else
		fprintf(fpinc, ",\n!%36s ", "");
              int argType = funcargtype[i];
              int isArray = isArrayArgType(argType);
	      fprintf(fpinc, "%-16s%.*s%s", funArgSym[argType].f77name,
                      (int)(funcargname[i].rm_eo - funcargname[i].rm_so),
                      line + funcargname[i].rm_so,
                      isArray ? "(*)" : "");
	    }
	  if ( funcargc )
	    fprintf(fpinc, ")\n");
	  fprintf(fpinc, "      %-16s%s\n\n", "EXTERNAL", sname);

	  /* fortran interface */
          const char *delegateName;
          char delegateNameBuf[MAX_FUNC_NAME_LEN + 7];
          size_t delegateNameLen = funcnameLen;
          /* emit prologue if needed */
          if (needsPrologue)
            {
              if (funRet[functype].needsPrologue)
                funRet[functype].prologue(fpint, (size_t)-1);
              for (size_t i = 0; i < funcargc; i++ )
                {
                  if (funArgSym[funcargtype[i]].needsPrologue)
                    funArgSym[funcargtype[i]].prologue(fpint, i);
                }
            }
          /* emit wrapper for type conversions if needed */
          if (needsExtraWrapper)
            {
              strcpy(delegateNameBuf, funcname);
              strcat(delegateNameBuf, "_fwrap");
              delegateNameLen += 6;
              delegateName = delegateNameBuf;
              fputs("static ", fpint);
              fprintf(fpint, (funRet[functype].convert
                              ?funRet[functype].convcfmt:funRet[functype].cfmt),
                      (int)delegateNameLen, delegateName);
              fputs("(", fpint);
              for (size_t i = 0; i < funcargc; i++ )
                {
                  if (i > 0)
                    fputs(", ", fpint);
                  fprintf(fpint, (funArgSym[funcargtype[i]].convert
                                  ?funArgSym[funcargtype[i]].convcfmt
                                  :funArgSym[funcargtype[i]].cfmt),
                          (int)(funcargname[i].rm_eo - funcargname[i].rm_so),
                          line + funcargname[i].rm_so);
                }
              fputs(")\n{\n", fpint);
              if (functype != ISVOID)
                {
                  fputs("  ", fpint);
                  fprintf(fpint, funRet[functype].cfmt, 1, "v");
                  fprintf(fpint, ";\n"
                          "  v = %s(", funcname);
                }
              else
                fprintf(fpint, "  %s(", funcname);
              for (size_t i = 0; i < funcargc; i++ )
                {
                  if (i > 0)
                    fputs(", ", fpint);
                  if (funArgSym[funcargtype[i]].convert)
                    {
                      funArgSym[funcargtype[i]]
                        .convert(fpint,
                                 line + funcargname[i].rm_so,
                                 (size_t)(funcargname[i].rm_eo
                                          - funcargname[i].rm_so), CONV_ARG);
                    }
                  else
                    fprintf(fpint, "%.*s",
                            (int)(funcargname[i].rm_eo - funcargname[i].rm_so),
                            line + funcargname[i].rm_so);
                }
              fputs(");\n", fpint);
              if (functype != ISVOID)
                {
                  fputs("  return ", fpint);
                  if (funRet[functype].convert)
                    funRet[functype].convert(fpint, "v", 1, CONV_RET);
                  else
                    fputc('v', fpint);
                  fputs(";\n", fpint);
                }
              fputs("}\n", fpint);
            }
          else
            delegateName = funcname;
	  if ( functype == ISVOID )
	    fprintf(fpint, "FCALLSCSUB");
	  else
	    fprintf(fpint, "FCALLSCFUN");
	  fprintf(fpint, "%zd ", funcargc);
	  fprintf(fpint, "(");
	  if ( functype != ISVOID )
	    fprintf(fpint, "%s, ", funRet[functype].cfint);
	  fprintf(fpint, "%s, ", delegateName);
	  for (size_t i = 0; i < funcnameLen; ++i)
            sname[i] = (char)toupper((int) sname[i]);
	  fprintf(fpint, "%s, ", sname);
	  for (size_t i = 0; i < funcnameLen; ++i)
            sname[i] = (char)tolower((int) sname[i]);
	  fprintf(fpint, "%s", sname);
	  for (size_t i = 0; i < funcargc; i++ )
	    {
	      fprintf(fpint, ", %s", funArgSym[funcargtype[i]].cfint);
	    }
	  fprintf(fpint, ")\n");


	  if ( funcnameLen == xnameLen
               && memcmp(funcname, xname, funcnameLen) == 0 )
	    {
	      char xline[128];
              size_t xlineLen = 0;
	      int nch;

	      /* C Quick Guide */

	      cdoc[ncdoc].naline = 0;
	      cdoc[ncdoc].text   = NULL;
	      cdoc[ncdoc].fname  = strdup(funcname);

	      nch = sprintf(xline, funRet[functype].cfmt,
                            (int)funcnameLen, funcname);
              xline[nch++] = ' ';
              xline[nch++] = '(';
              xline[nch] = '\0';
              xlineLen = (size_t)nch;

	      if ( (funcargc == 1 && funcargtype[0] == ISVOID) ) funcargc = 0;

	      for (size_t i = 0; i < funcargc; i++ )
		{
		  if (i)
                    {
                      strcat(xline, ", ");
                      xlineLen += 2;
                    }

                  /* extract full argument text from match */
                  char farg[128];
                  /* - 1 to omit closing paren ) or comma , */
                  int nchn = snprintf(farg, sizeof (farg), "%.*s",
                                      (int)(funcargfull[i].rm_eo
                                            - funcargfull[i].rm_so - 1),
                                      line + funcargfull[i].rm_so);
                  if (nchn < 0)
                    abort();
                  /* compress white-space */
                  nchn = (int)compress_whitespace((size_t)nchn, farg);
		  if ( (xlineLen + (size_t)nchn) > (size_t)80 )
		    {
                      if (i) xline[--xlineLen] = 0;
		      cdoc[ncdoc].aline[cdoc[ncdoc].naline++] = strdup(xline);
		      sprintf(xline, "%*s", nch, "");
                      xlineLen = (size_t)nch;
		    }
		  strcat(xline, farg);
                  xlineLen += (size_t)nchn;
		}
	      strcat(xline, ");");
	      cdoc[ncdoc].aline[cdoc[ncdoc].naline++] = strdup(xline);
	      cdoc[ncdoc].text  = strdup(xdes);

	      ncdoc++;

	      /* Fortran Quick Guide */

	      fdoc[nfdoc].naline = 0;
	      fdoc[nfdoc].text   = NULL;
	      fdoc[nfdoc].fname  = strdup(funcname);

	      if ( functype == ISVOID )
		nch = sprintf(xline, "SUBROUTINE %s", xname);
	      else
		nch = sprintf(xline, "%s FUNCTION %s", funArgSym[functype].f77name, xname);

	      if ( (funcargc == 1 && funcargtype[0] == ISVOID) ) funcargc = 0;
              if (funcargc) strcat(xline, " ("), nch += 2;

              xlineLen = (size_t)nch;

	      for (size_t i = 0; i < funcargc; i++ )
		{
		  if (i)
                    {
                      strcat(xline, ", ");
                      xlineLen += 2U;
                    }

                  char farg[128];
                  /* FIXME: optional empty argument name unhandled */
                  int argType = funcargtype[i];
                  int isArray = isArrayArgType(argType);
		  int nchn
                    = snprintf(farg, sizeof (farg), "%s %.*s%s",
                               funArgSym[argType].f77name,
                               (int)(funcargname[i].rm_eo
                                     - funcargname[i].rm_so),
                               line + funcargname[i].rm_so,
                               isArray ? "(*)" : "");
                  if (nchn < 0)
                    abort();
		  if ( (xlineLen + (size_t)nchn) > 80 )
		    {
                      if (i) xline[--xlineLen] = 0;
		      fdoc[nfdoc].aline[fdoc[nfdoc].naline++] = strdup(xline);
		      sprintf(xline, "%*s", nch, "");
                      xlineLen = (size_t)nch;
		    }
		  strcat(xline, farg);
                  xlineLen += (size_t)nchn;
		}
	      if ( funcargc ) strcat(xline, ")");
	      fdoc[nfdoc].aline[fdoc[nfdoc].naline++] = strdup(xline);
	      fdoc[nfdoc].text  = strdup(xdes);

	      nfdoc++;
	    }
	}
      else if ( memcmp(line, "#define", 7) == 0 )
	{
	  pline = line;
	  pline += 7;
	  while ( isspace((int) *pline) ) pline++;
	  parname = pline;
	  size_t len = strlen(pline);
          size_t i = 0;
	  for (; i < len; i++ )
	    {
	      if ( isspace((int) pline[i]) ) break;
	    }
	  if ( i == len ) continue;
	  parname[i] = 0;
	  pline += i+1;
	  while ( isspace((int) *pline) ) pline++;
	  if ( ! (isdigit((int) *pline) || *pline == '-') ) continue;
	  parvalue = atoi(pline);

	  /* fortran include */
	  fprintf(fpinc, "      INTEGER    %-22s\n"
                  "      PARAMETER (%-22s = %2d)\n", parname, parname,
                  parvalue);
	}
      else if (!regexec(&cppCondRE, line, maxMatch, reMatch, 0)
               && ((cppSwitchLen = reMatch[2].rm_eo - reMatch[2].rm_so) == 5)
               && ((size_t)(cppSymLen = reMatch[3].rm_eo - reMatch[3].rm_so)
                   == sizeof (cplusplus_macro) - 1)
               && !memcmp(line + reMatch[3].rm_so, cplusplus_macro,
                          sizeof (cplusplus_macro) - 1))
	{
          fprintf(stderr, "Found conditional C++ block, skipping to #else\n");
          while ((lineLen = getline(&line, &lineBufSize, fpin)) >= 0)
            {
              ++lineno;
              static const char endif_str[] = "endif";
              if (!regexec(&cppElseRE, line, maxMatch, reMatch, 0)
                  || (!regexec(&cppCondRE, line, maxMatch, reMatch, 0)
                      && ((size_t)(cppSymLen = reMatch[1].rm_eo - reMatch[1].rm_so)
                          == sizeof (endif_str) - 1)
                      && !memcmp(line + reMatch[1].rm_so, endif_str,
                                 sizeof (endif_str) - 1)))
                break;
            }
        }
      else if (detectComment(&line, &lineLen, &lineBufSize,
                             maxMatch, reMatch,
                             xname, &xnameLen, xdes,
                             fpin, fpinc, fpint))
        ;
      else
	{
	  if ( lineLen > 1 )
	    printf("skip: line %3d  size %3zd  %s%s", lineno, lineLen, line,
                   line[lineLen-1]=='\n'?"":"missing new-line\n");
	}
    }

  fputs("\n"
        "#endif\n", fpint);

  fclose(fpin);
  fclose(fpinc);
  fclose(fpint);

  qsort(cdoc, ncdoc, sizeof(Docu), doccmp);
  qsort(fdoc, nfdoc, sizeof(Docu), doccmp);


  char *filename;
  size_t doc_root_len = strlen(doc_root);
  {
    size_t max_len = 0;
    for (size_t i = 0; i < (size_t)fname_list_size; ++i)
      {
        size_t len = strlen(fname_list[i].name);
        fname_list[i].len = len;
        if (len > max_len)
          max_len = len;
      }
    /* path to document root, separating /, max of path within root,
     * terminating '\0'  */
    filename = (char*) malloc(doc_root_len + 1 + max_len + 1);
  }

  memcpy(filename, doc_root, doc_root_len);
  filename[doc_root_len] = '/';
  memcpy(filename + doc_root_len + 1,
         fname_list[NAME_C_QUICK_REF].name,
         fname_list[NAME_C_QUICK_REF].len + 1);
  fp = fopen(filename, "w");
  if ( fp )
    {
      fputs("C Quick Reference\n"
            "-----------------\n\n", fp);

      doctotxt(fp, cdoc, ncdoc);
      fclose(fp);
    }
  else
    {
      fprintf(stderr, "warning: cannot open documentation output file %s, %s",
              filename, strerror(errno));
    }

  memcpy(filename + doc_root_len + 1,
         fname_list[NAME_F_QUICK_REF].name,
         fname_list[NAME_F_QUICK_REF].len + 1);
  fp = fopen(filename, "w");
  if ( fp )
    {
      fputs("Fortran Quick Reference\n"
            "-----------------------\n\n", fp);

      doctotxt(fp, fdoc, nfdoc);
      fclose(fp);
    }
  else
    {
      fprintf(stderr, "warning: cannot open documentation output file %s, %s",
              filename, strerror(errno));
    }

  memcpy(filename + doc_root_len + 1,
         fname_list[NAME_C_QUICK_REF_TEX].name,
         fname_list[NAME_C_QUICK_REF_TEX].len + 1);
  fp = fopen(filename, "w");
  if ( fp )
    {
      fputs("\\chapter*{Quick Reference}\n"
            "\\addcontentsline{toc}{chapter}{Quick Reference}\n"
            "\n"
            "This appendix provide a brief listing of the C language bindings of the\n"
            "CDI library routines:\n"
            "\n", fp);

      doctotex(fp, cdoc, ncdoc);
      fclose(fp);
    }
  else
    {
      fprintf(stderr, "warning: cannot open documentation output file %s, %s",
              filename, strerror(errno));
    }

  memcpy(filename + doc_root_len + 1,
         fname_list[NAME_F_QUICK_REF_TEX].name,
         fname_list[NAME_F_QUICK_REF_TEX].len + 1);
  fp = fopen(filename, "w");
  if ( fp )
    {
      fputs("\\chapter*{Quick Reference}\n"
            "\\addcontentsline{toc}{chapter}{Quick Reference}\n"
            "\n"
            "This appendix provide a brief listing of the Fortran language bindings of the\n"
            "CDI library routines:\n"
            "\n", fp);

      doctotex(fp, fdoc, nfdoc);
      fclose(fp);
    }
  else
    {
      fprintf(stderr, "warning: cannot open documentation output file %s, %s",
              filename, strerror(errno));
    }
  free(filename);
}

static void
build_header_name(const char *fname, char *cppMacro)
{
  size_t len = strlen(fname);
  for (size_t i = 0; i < len; ++i)
    switch (fname[i])
      {
      case '.':
        cppMacro[i] = '_';
        break;
      default:
        cppMacro[i] = (char)toupper((int)fname[i]);
      }
  cppMacro[len] = '_';
  cppMacro[len + 1] = '\0';
}

int main(int argc, char *argv[])
{
  char *fname;
  char *cp;
  const char *doc_root = default_doc_root;
  char fnameinc[128], fnameint[128];
  size_t len;

  int optargCount = 0;
  {
    int opt;
    while ((opt = getopt(argc, argv, "d:")) != -1)
      switch (opt) {
      case 'd':
        doc_root = optarg;
        optargCount = 2;
        break;
      default: /* '?' */
        fprintf(stderr, "Usage: %s [-d DOCROOT] includefile\n", argv[0]);
        return(EXIT_FAILURE);
      }
  }


  if ( argc != 2 + optargCount)
    {
      printf("Usage:  %s [-d DOCROOT] includefile\n", argv[0]);
      return (1);
    }

  fname = argv[1 + optargCount];

  cp = strrchr(fname, '.');
  if ( cp == NULL ) len = strlen(fname);
  else              len = (size_t)(cp - fname);

  memcpy(fnameinc, fname, len);
  memcpy(fnameint, fname, len);

  strcpy(fnameinc+len, ".inc");
  strcpy(fnameint+len, "Fortran.c");

  fortran_interface(fname, fnameinc, fnameint, doc_root);

  return (0);
}

static inline size_t
compress_whitespace(size_t len, char str[])
{
  size_t wpos = 0;
  size_t i = 0;
  /* skip leading white-space */
  while (i < len && (isblank(str[i]) || str[i] == '\n'))
    ++i;
  /* after the leading white-space the following is
   * an alternation of white- and non-white-space
   * characters, where sequences of the latter will
   * be compressed to a single space */
  while (i < len)
    {
      /* handle white-space */
      while (i < len && !(isblank(str[i]) || str[i] == '\n'))
        str[wpos++] = str[i++];
      /* skip non-white-space */
      size_t wscount = 0;
      while (i < len && (isblank(str[i]) || str[i] == '\n'))
        ++i, ++wscount;
      if (wscount)
        str[wpos++] = ' ';
    }
  str[wpos] = '\0';
  return wpos;
}

enum {
  REGEX_MAX_ERRSTRLEN = 1024,
};


static size_t
symRegexCompile(size_t numSyms, struct symbol symList[],
                char **line, size_t *lineBufSize)
{
  size_t maxMatch = 0;
  for (size_t sym = 0; sym < numSyms; ++sym)
    {
      if (reCompile(&symList[sym].preg, symList[sym].parseRE,
                    line, lineBufSize))
        exit(EXIT_FAILURE);
      if (symList[sym].nameMatch > maxMatch)
        maxMatch = symList[sym].nameMatch;
    }
  return maxMatch;
}

static int
reCompile(regex_t *restrict RE, const char *restrict REstring,
          char * restrict *restrict lineBuf, size_t * restrict lineBufSize)
{
  int errcode;
  if ((errcode = regcomp(RE, REstring, REG_EXTENDED)))
    {
      char *restrict line;
      size_t resize;
      if (*lineBufSize < REGEX_MAX_ERRSTRLEN
          && (line = (char*) realloc(*lineBuf, resize = REGEX_MAX_ERRSTRLEN)))
        {
          *lineBuf = line;
          *lineBufSize = resize;
          regerror(errcode, RE, line, *lineBufSize);
          fprintf(stderr, "Error compiling regular expression: %s: %s\n",
                  REstring, *lineBuf);
        }
    }
  return errcode;
}

/* emit conversion code for MPI_Comm argument */
static int cfMPICommConvert(FILE *outfp, const char *argName,
                            size_t argNameLen, enum conversionType part)
{
  int retval = 0;
  switch (part)
    {
    case CONV_ARG:
      retval
        = fprintf(outfp, "MPI_Comm_f2c(%.*s)", (int)argNameLen, argName);
      break;
    case CONV_RET:
      retval = fprintf(outfp, "MPI_Comm_c2f(%.*s)", (int)argNameLen, argName);
      break;
    }
  return retval;
}

/* emit conversion code for Xt_idxlist argument */
static int cfXtIdxlistConvert(FILE *outfp, const char *argName,
                              size_t argNameLen, enum conversionType part)
{
  int retval = 0;
  switch (part)
    {
    case CONV_ARG:
      retval
        = fprintf(outfp, "(*(Xt_idxlist *)%.*s)", (int)argNameLen, argName);
      break;
    case CONV_RET:
      abort();
      break;
    }
  return retval;
}

static int cfVoidFuncPrologue(FILE *outfp, size_t argNum)
{
  int retval
    = fprintf(outfp, "\n#undef ROUTINE_%zu\n#define ROUTINE_%zu %s\n",
              argNum+1, argNum+1, "(void (*)(void))");
  return retval;
}

enum {
  FOUND_NOTHING,
  FOUND_COMMENT,
  FOUND_DOCCOMMENT,
};

static int detectComment(char **line_, ssize_t *lineLen, size_t *lineBufSize,
                         size_t maxMatch, regmatch_t reMatch[],
                         char *xname, size_t *xnameLen,
                         char *xdes,
                         FILE *fpin, FILE *fpinc, FILE *fpint)
{
  char *restrict line = *line_;
  int matchType;
  do {
    if (!regexec(&docCommentRE, line, maxMatch, reMatch, 0))
      {
        /* found documentation comment */
        size_t nameMatchLen = (size_t)(reMatch[1].rm_eo - reMatch[1].rm_so),
          docMatchLen = (size_t)(reMatch[2].rm_eo - reMatch[2].rm_so);
        memcpy(xname, line + reMatch[1].rm_so, nameMatchLen);
        xname[nameMatchLen] = 0;
        *xnameLen = nameMatchLen;
        memcpy(xdes, line + reMatch[2].rm_so, docMatchLen);
        {
          char *eol = xdes;
          while ((eol = strchr(eol, '\n')))
            {
              ++eol;
              /* delete whitespace following newline */
              size_t squeezeLen = strspn(eol, " \t*");
              char *startoftext = eol + squeezeLen;
              memmove(eol, startoftext, docMatchLen - (size_t)(eol - xdes));
              docMatchLen -= squeezeLen;
            }
        }
        xdes[docMatchLen] = 0;
        if (verbose || debug)
          printf("Found documentation for \"%s\": \"%s\"\n", xname,
                 xdes);
        matchType = FOUND_DOCCOMMENT;
        break;
      }
    else if (!regexec(&commentRE, line, maxMatch, reMatch, 0))
      {
        size_t commentLen = (size_t)(reMatch[1].rm_eo - reMatch[1].rm_so);
        char *comment = line + reMatch[1].rm_so;
        {
          char savedCommentEnd = comment[commentLen];
          comment[commentLen] = '\0';
          /* fortran include */
          fputs("!\n", fpinc);
          char *cline = comment;
          do {
            cline += strspn(cline, " \t*");
            char *eol = strchr(cline, '\n');
            if (!eol)
              eol = comment + commentLen;
            size_t lineLen = (size_t)(eol - cline);
            fprintf(fpinc, "!  %.*s\n", (int)lineLen, cline);
            cline = (eol != comment + commentLen) ? eol + 1: NULL;
          } while (cline);
          fputs("!\n", fpinc);
          comment[commentLen] = savedCommentEnd;
        }
        /* fortran interface */
        fprintf(fpint, "\n/*  %.*s  */\n\n", (int)commentLen, comment);
        matchType = FOUND_COMMENT;
        break;
      }
    /* found comment start, read further lines and retry */
    else if (!regexec(&commentStartRE, line, maxMatch, reMatch, 0))
      {
        int foundCommentEnd = 0;
        char *lineExtension = NULL;
        size_t extSize = 0;
        do {
          ssize_t extLen;
          if ((extLen = getline(&lineExtension, &extSize, fpin)) <= 0)
            break;
          if ((size_t)(*lineLen + extLen) >= *lineBufSize)
            if (!(line = realloc(line, (size_t)(*lineLen + extLen + 1))))
              exit(EXIT_FAILURE);
          memcpy(line + *lineLen, lineExtension, (size_t)extLen + 1);
          *lineLen += extLen;
          foundCommentEnd
            = !regexec(&commentEndRE, lineExtension, maxMatch, reMatch, 0);
        } while (!foundCommentEnd);
      }
    else
      /* found no comment at all */
      break;
  } while (1);
  *line_ = line;
  return matchType;
}

static inline int
isArrayArgType(int argType)
{
  return argType == ISFLOATV
    || argType == ISFLOATVV
    || argType == ISDOUBLEV
    || argType == ISDOUBLEVV
    || argType == ISINTV
    || argType == ISINTVV;
}


/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */

