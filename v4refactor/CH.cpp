#include "CH.h"
#include "HeMS.h"




// JR FIX THIS: check this - rewrite <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
STELLAR_TYPE CH::EvolveToNextPhase() {

    STELLAR_TYPE stellarType = STELLAR_TYPE::MS_GT_07;

    if (Age() < m_InterimState.HurleyTimescales(TS::MS)) {       // evolving off because of age?
        stellarType = STELLAR_TYPE::MS_GT_07;                                   // no - must have spun down - evolve as MS star now
        m_InterimState.SetCHE(false);                                           // evolved CH->MS
        
        // if BRCEK core mass calculations enabled, initialise the core mass based on current mass and central helium fraction
        if ((OPTIONS->MainSequenceCoreMassPrescription() == MS_CORE_MASS_PRESCRIPTION::BRCEK) && (utils::Compare(MZAMS(), BRCEK_LOWER_MASS_LIMIT) >= 0))
            m_InterimState.SetMSCoreMass(MainSequence::CalculateCNOprocessedCoreMass_Brcek2025(Mass(), HeAbundanceCore()));
    }
    else {                                                                      // yes
        stellarType = STELLAR_TYPE::NAKED_HELIUM_STAR_MS;                       // evolve as HeMS star now
        m_InterimState.SetAge(0.0);                                             // can't use Hurley et al. 2000, eq 76 here - timescales(tHe) not calculated yet
        m_InterimState.SetTau(0.0);
        m_InterimState.SetCHE(true);                                            // stayed on MS as CH
    }

    return stellarType;
}

