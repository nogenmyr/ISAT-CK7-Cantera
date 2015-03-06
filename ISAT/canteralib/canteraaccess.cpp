#include "thermo/ThermoFactory.h"
#include "base/config.h"
#include "clib/Cabinet.h"

using namespace Cantera;

typedef Cabinet<ThermoPhase> ThermoCabinet;
typedef integer status_t;

/*
std::string f2string_loc(const char* s, ftnlen n)
{
    int k;
    std::string ss = "";
    for (k = 0; k < n; k++) {
        if (s[k] == '\0') {
            break;
        }
        ss += s[k];
    }
    return ss;
}
*/

inline ThermoPhase* _fph(const integer* n)
{
    return &ThermoCabinet::item(*n);
}

extern "C" {
    integer nasadata_(
        const integer* n, 
        double* thermo_ns
        )
    {
        try {
            int type, k, id, tid;
            doublereal c[15];
            doublereal minTemp, maxTemp, refPressure;
            int nsp = _fph(n)->nSpecies();
            SpeciesThermo& sp = _fph(n)->speciesThermo();

            for (k = 0; k < nsp; k++) {
                // get the NASA coefficients in array c
                sp.reportParams(k, type, c, minTemp, maxTemp, refPressure);
//                if (type != NASA) stop!
                for (id=0; id<15; id++){
                    tid = k+id*nsp;
                    thermo_ns[tid]=c[id];
                    }
            }
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }
}

