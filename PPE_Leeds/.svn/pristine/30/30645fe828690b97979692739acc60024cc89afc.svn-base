RM = rm -f
RM_R = $(RM) -r
MV = mv

MKEXP = ../util/mkexp/mkexp
MKEXP_FLAGS = EXP_ID=$${EXP_TYPE} EXP_TYPE=$$EXP_TYPE ENVIRONMENT=$$ENVIRONMENT MODEL_ROOT=. WORK_ROOT=. DATA_ROOT=. MODEL_SUBDIR= EXPERIMENTS_SUBDIR=$(CHECK_TARGET)/$$ENVIRONMENT

DIFF = diff -I '^ *\# *\$$Id: '
DIFF_R = $(DIFF) $(DIFF_FLAGS) -r -x .svn

ALL_TARGET = reference_experiments
CHECK_TARGET = experiments

ENVIRONMENTS = blizzard thunder DEFAULT
KINDS = amip sstClim
QUALITIES = LR MR

TEST = examples/amiptest.config

all: all.make

clean: all.clean check.clean

check: check.make
	$(DIFF_R) $(ALL_TARGET) $(CHECK_TARGET)
	@echo Info: passed all tests

install doc:


all.make: all.clean check.make
	$(MV) $(CHECK_TARGET) $(ALL_TARGET)

all.clean:
	$(RM_R) $(ALL_TARGET)


check.make: check.clean
	@ \
	set -e; \
	for ENVIRONMENT in $(ENVIRONMENTS); do \
	  for EXP_KIND in $(KINDS); do \
	    for EXP_QUALITY in $(QUALITIES); do \
	      EXP_TYPE="$$EXP_KIND-$$EXP_QUALITY"; \
	      echo $(MKEXP) $(TEST) \
	        EXP_ID=$$EXP_TYPE ENVIRONMENT=$$ENVIRONMENT; \
	      $(MKEXP) $(TEST) $(MKEXP_FLAGS) > /dev/null; \
	    done; \
	  done; \
	done

check.clean:
	$(RM_R) $(CHECK_TARGET)
