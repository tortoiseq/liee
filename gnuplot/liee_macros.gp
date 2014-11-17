set macros

# TERMINAL
#---------
TERMINAL='set term png enhanced size width_plt, height_plt'


# SCAN_DATA: extract some parameters and scan through the data-file to find global extrema
#-----------
cmd_extract="grep '##' '%s' | awk '{print $2}' > '%s.param'"
cmd_stats="liee_tools.py stats '%s' '%s.param'"

SCAN_DATA='system( sprintf( cmd_extract, data, data ) );\
system( sprintf( cmd_stats, data, data ) );\
load data.".param";'


# Constants
#----------
CONVm = 5.2917720859e-11
CONVnm = 5.2917720859e-2
CONVeV = 27.211
CONVs = 2.418884326505e-17
CONVfs = 2.418884326505e-2
CONVVoM = 5.1421e11
CONVV = CONVVoM * CONVm

hbar = 1.05457148e-34
me = 9.10938188e-31
Ce = 1.60217646e-19

