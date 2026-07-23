This directory is meant to contain various emission specification examples.

#------------------
It must contain the 3 compulsory (for the jobscript templates sanity) examples:

* HAM-only:

  * emi_spec_hammoz_defaults+psoa_ham.txt: hammoz defaults emissions for HAM, with prescribed SOA (with perpetual year, actual year value set to 2000)

* MOZ-only:

  * emi_spec_hammoz_defaults_moz.txt: hammoz defaults emissions for MOZ

* HAMMOZ:

  * emi_spec_hammoz_default_hammoz.txt: hammoz defaults emissions for HAMMOZ

#------------------
It also contain the additional example:

  * emi_spec_hammoz_defaults+psoa_transient_ham.txt: hammoz defaults emissions for HAM, with prescribed SOA, transient emissions.

  * emi_spec_maccity_biomassburning+psoa_transient_ham.txt: same as the above, but with maccity_biomassburning inventory for the FFIRE and GFIRE sectors

  * emi_spec_bs_biomassurning+psoa_transient_ham.txt: same as the above, but with GFAS inventory for the FFIRE and GFIRE sectors. Note that GFAS is used in its 'daily' flavour here. It could also be used as monthly data (simply replace 'daily' by 'monthly' in the relevant sector definition lines.)

#------------------
You may store your own emission specification examples here.

