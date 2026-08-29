#include "ONeWD.h"
#include "Star.h"








///////// constituent functions <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<





STELLAR_TYPE ONeWD_Constituent::EvolveToNextPhaseBinaryOutcome() {
    if (!m_Star) return STELLAR_TYPE::OXYGEN_NEON_WHITE_DWARF;
    return m_Star->IsMassAboveChandrasekhar() ? m_Star->ResolveAIC() : STELLAR_TYPE::OXYGEN_NEON_WHITE_DWARF;
}


bool ONeWD_Constituent::IsSupernovaBinaryOutcome() const {
    return m_Star ? m_Star->IsMassAboveChandrasekhar() : false;
}


STELLAR_TYPE ONeWD_Constituent::ResolveSupernovaBinaryOutcome() {
    return m_Star ? m_Star->ResolveAIC() : STELLAR_TYPE::OXYGEN_NEON_WHITE_DWARF;
}
