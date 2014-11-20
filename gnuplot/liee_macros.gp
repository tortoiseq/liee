set macros

cmd_stats="liee_tools.py stats '%s' '%s.param'"
cmd_extract="grep '##' '%s' | awk '{print $2}' >> '%s.param'"
cmd_pot="liee_tools.py potential %f %f %d %f %f %d '%s'"
cmd_transpose="liee_tools.py transpose '%s' '%s.t'"
cmd_cutrow="tail -n +2 '%s' > '%s'.tmp"
#cmd_cutcol="grep '##' '%s' > '%s'.tmp && grep -v '##' '%s' | cut --fields=2-" >> '%s'.tmp
#cmd_cutswap="rm '%s' && mv '%s'.tmp '%s'"

TERMINAL='set term png enhanced size width_plt, height_plt;'

# @PARAMETER: extract and load ##parameters from filename in variable "maIN"
PARAMETER='system( sprintf( cmd_extract, maIN, maIN ) ); load maIN.".param";'

# @EXTREMA: find global extrema and extract ##parameters from the data-file given by the variable "marco_in" 
EXTREMA='@PARAMETER; system( sprintf( cmd_stats, maIN, maIN ) ); load maIN.".param";'

# @POTENTIAL: write tabulated potential values V(r,t) to the file named by variable "maIN". 
#            the set of variables: [r0, r1, Nr, t0, t1, Nt] defines a regular grid of requested potential values.
#            the maximum and minimum potential values are stored in [V_min, V_max] and the respective time-indices in [V_min_t, V_max_t]
POTENTIAL='system( sprintf( cmd_pot, ma_r0, ma_r1, ma_Nr, ma_t0, ma_t1, ma_Nt, maIN ) ); @PARAMETER'

# @TRANSPOSE: 
TRANSPOSE='system( sprintf( cmd_transpose, maIN, maIN ) );'


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

