#!/bin/bash



## THIS MACRO SHOULD RUN PRINTCANVASES
#root -l -q 'printCanvases.C++(0)' #PY8_UNFDIR_DATA    ## CONSIDER HANDLING THIS OUTPUT DIFF FROM OTHER PLOTS (i.e. use makeFinalPlots for this)
#root -l -q 'printCanvases.C++(1)' #PY8_UNFDIR_CLOSURE ## CONSIDER HANDLING THIS OUTPUT DIFF FROM OTHER PLOTS (i.e. use makeFinalPlots for this)
#root -l -q 'printCanvases.C++(2)' #NLO_UNFDIR_DATA    ## CONSIDER HANDLING THIS OUTPUT DIFF FROM OTHER PLOTS (i.e. use makeFinalPlots for this)
#root -l -q 'printCanvases.C++(3)'  #NLO_UNFDIR_CLOSURE## CONSIDER HANDLING THIS OUTPUT DIFF FROM OTHER PLOTS (i.e. use makeFinalPlots for this)
#root -l -q 'printCanvases.C++(4)' #JERDIR_MC  
#root -l -q 'printCanvases.C++(5)' #VTXWDIR_MC 
#root -l -q 'printCanvases.C++(6)' #JETQA      
#root -l -q 'printCanvases.C++(7)' #JETIDQA    
#root -l -q 'printCanvases.C++(8)' #JETTRIGQA  
root -l -q 'printCanvases.C++(9)' #JETTRIGQA  

## THIS MACRO SHOULD EXECUTE [./] A BASH SCRIPT THAT SETS UP THE NOTE/AN ENVIRONMENT, THEN REMAKES PLOT LISTS, THEN BUILDS THE AN [do_ANStuff.sh]
#./doANStuff.sh