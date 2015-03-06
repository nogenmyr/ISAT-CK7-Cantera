# ISAT-CK7-Cantera

The original library can be found here: https://tcg.mae.cornell.edu/ISATCK7/

That library uses ChemkinII as a backend for thermochemistry. The library provided here uses Cantera 2.1 or 2.2 for thermochemistry. The code in ice-pic/g_jacadf.f was however very dependent on Chemkins internal work arrays. This is yet to port to Cantera. I currently do not see any easy way to do that unfortunately. (one way can be to create those work arrays using Cantera. That is however not straight-forward.) This means that the ice-pic library is not functional. Ths ISAT code does however work without problems.
