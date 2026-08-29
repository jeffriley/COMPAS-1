#include "ZDependent.h"
#include "astroUtils.h"
#include "utils.h"

#include <algorithm>
#include <cmath>


using std::max;
using std::min;


/*
 * Resize
 *
 * @brief
 * Resize the metallicity cache.
 * 
 * The cache can be resized up (to a maximum size of MAX_Z_CACHE_SIZE), or
 * down to a minimum size of 1.  The p_CacheCapacity will be bounds checked
 * against these limits.
 * 
 * If the cache is resized down, the tail beyond the new size will be destroyed
 * and memory freed.
 *
 * If the cache is resized up, the existing cache entries are preserved, and the
 * new tail is value-initialised.
 * 
 * void Resize(const SizeT p_CacheCapacity)
 * 
 * @param       p_CacheCapacity                 The new cache capacity (number of slots)
 */
void ZDependent::Resize(const SizeT p_CacheCapacity) {

    const SizeT capacity = std::clamp(p_CacheCapacity, SizeT{1}, MAX_Z_CACHE_SIZE); // Just in case...

    if (capacity == m_CacheCapacity) return;                                        // Nothing to do

    // Resize the vectors - shrink or grow.
    //
    // On shrink: destroy tail [n,old) (runs ~ZDependentT, frees inner vectors)
    // On grow  : preserve [0,old) + value-init new tail
    m_Keys.resize(capacity);
    m_Cache.resize(capacity);

    // Adjust bookkeeping - fix on shrink; no-op on grow
    if (m_CacheCount > capacity)     m_CacheCount = capacity;                       // Keep at most m_CacheCapacity surviving entries
    if (m_LastFound >= m_CacheCount) m_LastFound = 0;                               // Reset last found pointer if necessary

    m_NextEvict    %= capacity;                                                     // Next eviction candidate 
    m_CacheCapacity = capacity;                                                     // Set new cache capacity
}


/*
 * ResolveSlot
 *
 * @brief
 * Returns the cache slot number for the requested metallicity.
 * 
 * Looks p_Z up in the cache and, if found (cache hit), returns the slot number.
 * If p_Z is not found (cache miss), calculates the values, inserts them into the
 * cache (with FIFO eviction if the cache is full), with key p_Z, and returns the
 * slot number.
 * 
 * The lookup of the cache key (p_Z) is performed with an absolute tolerance of
 * Z_CACHE_TOL - see constants.h.
 * 
 * 
 * SizeT ResolveSlot(const double p_Z) const
 *
 * @param       p_Z                             Metallicity to look for in cache
 * @return                                      Cache slot number for p_Z
 */
SizeT ZDependent::ResolveSlot(const double p_Z) const {

    // Note that we have multiple exit points (return statements) in this function.  Ordinarily
    // we'd avoid that, but it saves CPU cycles - important here.

    // Find existing slot for requested metallicity.
    // Return slot number here if cache hit.
    std::optional<SizeT> slot = FindSlot(p_Z);
    if (slot.has_value()) {
        m_LastFound = *slot;
        return *slot;
    }

    // Cache miss.
    // We only fall through to here if the required metallicity was not found.  We need to add
    // the metallicity to the keys vector, and calculate and add the Z-dependent values to the
    // cache vector.  If the cache is not full we just append the key and values to the vectors.
    // If the cache is full we overwrite (evict) an existing entry in both vectors.
    if (m_CacheCount < m_CacheCapacity) {                   // Cache full?
        slot = m_CacheCount++;                              // No - append
    }
    else {                                                  // Yes - cache full
        slot = m_NextEvict;                                 // Element to evict
        m_NextEvict = (m_NextEvict + 1) % m_CacheCapacity;  // Next eviction candidate
    }

    m_Keys[*slot] = p_Z;                                    // Populate key slot with metallicity
    CalculateForZ(p_Z, m_Cache[*slot]);                     // Calculate z-dependent values and populate cache slot
    m_LastFound = *slot;                                    // Last found slot

    return *slot;                                           // Return slot number
}


/*
 * CalculateForZ
 *
 * @brief
 * Calculates the ZDependentT struct values for the given metallicity and writes it into
 * the supplied slot. Called by ResolveSlot() on a cache miss.
 * 
 * We do not populate the optional structs, Hurley and Shikauchi, by default - these
 * are populated on-demand by the first access request.
 *
 *
 * void CalculateForZ(const double p_Z, ZDependentT& p_Slot) const
 *
 * @param       p_Z                             Metallicity to use to calculate the values
 * @param       p_Slot                          Cache slot to write into (overwrites the calculated values)
 */
void ZDependent::CalculateForZ(const double p_Z, ZDependentT& p_Slot) const {

    // Scalar Z-derived values
    p_Slot.log10Z         = std::log10(p_Z);

    p_Slot.zetaHurley     = p_Slot.log10Z - LOG10_ZSOL_HURLEY;
    p_Slot.zScaledHurley  = p_Z / ZSOL_HURLEY;

    p_Slot.zetaAnders     = p_Slot.log10Z - LOG10_ZSOL_ANDERS;
    p_Slot.zetaAsplund    = p_Slot.log10Z - LOG10_ZSOL_ASPLUND;

    p_Slot.zScaledAnders  = p_Z / ZSOL_ANDERS;
    p_Slot.zScaledAsplund = p_Z / ZSOL_ASPLUND;

    // ZAMS abundances per Pols 1998
    p_Slot.HAbundanceAtZAMS   = astro::CalculateZAMSHAbundance_Pols1998(p_Z);
    p_Slot.HeAbundanceAtZAMS  = astro::CalculateZAMSHeAbundance_Pols1998(p_Z);

    // Tout 1996 ZAMS coefficients
    p_Slot.ToutZAMSLCoeffs = astro::CalculateZAMSLuminosityCoefficients_Tout1996(p_Slot.zetaHurley);
    p_Slot.ToutZAMSRCoeffs = astro::CalculateZAMSRadiusCoefficients_Tout1996(p_Slot.zetaHurley);
}


/*
 * CalculateHurleyForZ
 *
 * @brief
 * Calculates the Hurley struct values for the given metallicity and writes it into
 * the supplied slot. Called by getters on first access.
 *
 * The order of operations matters here because later quantities depend on earlier ones.
 *
 *
 * void CalculateHurleyForZ(const double p_Z, ZDependentT& p_Slot) const
 *
 * @param       p_Z                             Metallicity to calculate for.
 * @param       p_Slot                          Cache slot to write into (overwritten in full).
 */
void ZDependent::CalculateHurleyForZ(const double p_Z, ZDependentT& p_Slot) const {

    p_Slot.Hurley.emplace();                 // Construct struct
    HurleyZDependentT& h = *p_Slot.Hurley;

    // Populate it
    h.rho               = p_Slot.zetaHurley + 1.0;
    h.gbRadiusXExponent = CalculateGBRadiusXExponent_Hurley2000(p_Slot.zetaHurley);

    h.massCutoffs       = CalculateMassCutoffs_Hurley2000(p_Z, p_Slot.zetaHurley);
    h.aCoeffs           = CalculateACoeffs_Hurley2000(p_Z, p_Slot.log10Z, p_Slot.zetaHurley);
    h.bCoeffs           = CalculateBCoeffs_Hurley2000(p_Z, p_Slot.log10Z, p_Slot.zetaHurley, h.rho, h.massCutoffs);

    h.alphas            = CalculateAlphaConstants_Hurley2000(h.massCutoffs[HURLEY_MCO::HeF], h.aCoeffs, h.bCoeffs);
    h.gammas            = CalculateGammaConstants_Hurley2000(h.aCoeffs);
    h.LConstants        = CalculateLuminosityConstants_Hurley2000(h.aCoeffs);
    h.RConstants        = CalculateRadiusConstants_Hurley2000(h.aCoeffs);

    h.lMinCHeB          = CHeB::CalculateMinLuminosity_Hurley2000(h.massCutoffs[HURLEY_MCO::HeF], h.massCutoffs[HURLEY_MCO::HeF], h.massCutoffs[HURLEY_MCO::FGB], h.bCoeffs);
}


/*
 * CalculateShikauchiForZ
 *
 * @brief
 * Calculates the Shikauchi struct values for the given metallicity and writes it into
 * the supplied slot. Called by getters on first access.
 *
 *
 * void CalculateShikauchiForZ(const double p_Z, ZDependentT& p_Slot) const
 *
 * @param       p_Z                             Metallicity to calculate for.
 * @param       p_Slot                          Cache slot to write into (overwritten in full).
 */
void ZDependent::CalculateShikauchiForZ(const double p_Z, ZDependentT& p_Slot) const {
    p_Slot.Shikauchi = std::optional<ShikauchiZDependentT>(CalculateCoeffs_Shikauchi2024(p_Z, p_Slot.log10Z));
}


////////////////////////////////////////////////////////////////////////////////////////
//                                                                                    //
//  These functions are here rather than utils or astroUtils because this is the only //
//  place they are called - if that changes in the future they could be moved.        //
//                                                                                    // 
////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////
//                                                                                    //
//                  HURLEY 2000 COEFFICIENTS AND CONSTANTS                            //
//                                                                                    //
////////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateACoeffs_Hurley2000
 *
 * @brief
 * Calculate Hurley et al. 2000 'a' coefficients. Values are from the table in
 * Appendix A of Hurley et al. 2000.
 *
 *
 * HurleyACoeffsT CalculateACoeffs_Hurley2000(const double p_Z, const double p_Sigma, const double p_Zeta) const
 *
 * @param       p_Z                             Metallicity
 * @param       p_Sigma                         Hurley sigma = log10(Z)
 * @param       p_Zeta                          Hurley zeta  = log10(Z / Zsol_Hurley)
 * @return                                      a coefficients array
 */
HurleyACoeffsT ZDependent::CalculateACoeffs_Hurley2000(const double p_Z, const double p_Sigma, const double p_Zeta) const {

    HurleyACoeffsT a;

    const double zeta2 = p_Zeta * p_Zeta;
    const double zeta3 = p_Zeta * zeta2;
    const double zeta4 = p_Zeta * zeta3;

    // First entry (index 0) is a dummy so our indices match Hurley et al. 2000.
    for (const auto& [idx, vals] : A_COEFF) {
        using enum AB_TCoeff;
        auto coeff = [&](AB_TCoeff p) { return vals.at(p); };
        a[idx] = coeff(ALPHA) + coeff(BETA) * p_Zeta + coeff(GAMMA) * zeta2 + coeff(ETA) * zeta3 + coeff(MU) * zeta4;
    }

    // Special cases - see Hurley et al. 2000.
    a[11] *= a[14];
    a[12] *= a[14];
    a[17]  = PPOW(10.0, max((0.097 - (0.1072 * (p_Sigma + 3.0))), max(0.097, min(0.1461, (0.1461 + (0.1237 * (p_Sigma + 2.0)))))));
    a[18] *= a[20];
    a[19] *= a[20];
    a[29]  = PPOW(a[29], (a[32]));
    a[33]  = max(0.6355 - (0.4192 * p_Zeta), max(1.25, min(1.4, 1.5135 + (0.3769 * p_Zeta))));
    a[42]  = min(1.25, max(1.1, a[42]));
    a[44]  = min(1.3, max(0.45, a[44]));
    a[49]  = max(a[49], 0.145);
    a[50]  = min(a[50], (0.306 + (0.053 * p_Zeta)));
    a[51]  = min(a[51], (0.3625 + (0.062 * p_Zeta)));
    a[52]  = max(a[52], 0.9);
    a[52]  = p_Z > 0.01 ? min(a[52], 1.0) : a[52];
    a[53]  = max(a[53], 1.0);
    a[53]  = p_Z > 0.01 ? min(a[53], 1.1) : a[53];
    a[57]  = max((0.6355 - (0.4192 * p_Zeta)), max(1.25, min(1.4, a[57])));
    a[62]  = max(0.065, a[62]);
    a[63]  = p_Z < 0.004 ? min(0.055, a[63]) : a[63];
    a[66]  = max(0.8, min(0.8 - (2.0 * p_Zeta), max(a[66], min(1.6, -0.308 - (1.046 * p_Zeta)))));
    a[68]  = max(0.9, min(a[68], 1.0));
    a[72]  = p_Z > 0.01 ? max(a[72], 0.95) : a[72];
    a[74]  = max(1.4, min(a[74], 1.6));
    a[75]  = max(1.0, min(a[75], 1.27));
    a[75]  = max(a[75], 0.6355 - (0.4192 * p_Zeta));
    a[76]  = max(a[76], -0.1015564 - (0.2161264 * p_Zeta) - (0.05182516 * zeta2));
    a[77]  = max((-0.3868776 - (0.5457078 * p_Zeta) - (0.1463472 * zeta2)), min(0.0, a[77]));
    a[78]  = max(0.0, min(a[78], 7.454 + (9.046 * p_Zeta)));
    a[79]  = min(a[79], max(2.0, -13.3 - (18.6 * p_Zeta)));
    a[80]  = max(0.0585542, a[80]);
    a[81]  = min(1.5, max(0.4, a[81]));

    // (Re)calculate special cases that depend on radius constants.
    // Uses alphaR = (a[58] * PPOW(a[66], a[60])) / (a[59] + PPOW(a[66], a[61])) as given in
    // Hurley et al. 2000, eq 21a (wrong in the arxiv version - says = a59*M**(a61)).
    a[64] = a[68] > a[66] ? (a[58] * PPOW(a[66], a[60])) / (a[59] + PPOW(a[66], a[61])) : max(0.091, min(0.121, a[64]));
    a[68] = min(a[68], a[66]);  // must calculate after calculating a[64]

    return a;
}


/*
 * CalculateBCoeffs_Hurley2000
 *
 * @brief
 * Calculate Hurley et al. 2000 'b' coefficients. Values are from the table in
 * Appendix B of Hurley et al. 2000.
 *
 *
 * HurleyBCoeffsT CalculateBCoeffs_Hurley2000(const double p_Z, const double p_Sigma, const double p_Zeta, const double p_Rho, const HurleyMassCutoffsT& p_MassCutoffs) const
 *
 * @param       p_Z                             Metallicity
 * @param       p_Sigma                         Hurley sigma = log10(Z)
 * @param       p_Zeta                          Hurley zeta  = log10(Z / Zsol_Hurley)
 * @param       p_Rho                           Hurley rho   = zeta + 1.0
 * @param       p_MassCutoffs                   Hurley mass cutoffs (already calculated for this Z)
 * @return                                      b coefficients array
 */
HurleyBCoeffsT ZDependent::CalculateBCoeffs_Hurley2000(const double p_Z, const double p_Sigma, const double p_Zeta, const double p_Rho, const HurleyMassCutoffsT& p_MassCutoffs) const {

    HurleyBCoeffsT b;

    const double zeta2 = p_Zeta * p_Zeta;
    const double zeta3 = p_Zeta * zeta2;
    const double zeta4 = p_Zeta * zeta3;
    const double zeta5 = p_Zeta * zeta4;
    const double rho2  = p_Rho * p_Rho;
    const double rho3  = p_Rho * rho2;

    for (const auto& [idx, vals] : B_COEFF) {
        using enum AB_TCoeff;
        auto coeff = [&](AB_TCoeff p) { return vals.at(p); };
        b[idx] = coeff(ALPHA) + coeff(BETA) * p_Zeta + coeff(GAMMA) * zeta2 + coeff(ETA) * zeta3 + coeff(MU) * zeta4;
    }

    // Special cases - see Hurley et al. 2000.
    b[ 1]  = min(0.54, b[1]);
    b[ 2]  = min(max(PPOW(10.0, (-4.6739 - (0.9394 * p_Sigma))), (-0.04167 + (55.67 * p_Z))), (0.4771 - (9329.21 * PPOW(p_Z, 2.94))));
    b[ 3]  = PPOW(10.0, max(-0.1451, (-2.2794 - (1.5175 * p_Sigma) - (0.254 * p_Sigma * p_Sigma))));
    b[ 3]  = p_Z > 0.004 ? max(b[3], 0.7307 + (14265.1 * PPOW(p_Z, 3.395))) : b[3];
    b[ 4] += 0.1231572 * zeta5;
    b[ 6] += 0.01640687 * zeta5;
    b[11] *= b[11];
    b[13] *= b[13];
    b[14]  = PPOW(b[14], b[15]);
    b[16]  = PPOW(b[16], b[15]);
    // rho = zeta + 1.0, so PPOW(p_Rho, ...) here is identical to the original
    // COMPAS form PPOW((xi + 1.0), ...), and the comparison p_Zeta > -1.0 is
    // identical to the original xi > -1.0. Confirmed equivalent.
    b[17]  = p_Zeta > -1.0 ? 1.0 - (0.3880523 * PPOW(p_Rho, 2.862149)) : 1.0;
    b[24]  = PPOW(b[24], b[28]);
    b[26]  = 5.0 - (0.09138012 * PPOW(p_Z, -0.3671407));
    b[27]  = PPOW(b[27], (2.0 * b[28]));
    b[31]  = PPOW(b[31], b[33]);
    b[34]  = PPOW(b[34], b[33]);
    b[36] *= b[36] * b[36] * b[36];
    b[37] *= 4.0;
    b[38] *= b[38] * b[38] * b[38];
    b[40]  = max(b[40], 1.0);
    b[41]  = PPOW(b[41], b[42]);
    b[44] *= b[44] * b[44] * b[44] * b[44];
    b[45]  = p_Rho <= 0.0 ? 1.0 : 1.0 - ((2.47162 * p_Rho) - (5.401682 * rho2) + (3.247361 * rho3));
    // b[46] needs the HeF and FGB mass cutoffs. These are now passed in via
    // p_MassCutoffs (calculated before this function is called - see CalculateForZ),
    // rather than read from singleton state. This is what keeps the function pure.
    b[46]  = -1.0 * b[46] * std::log10(p_MassCutoffs[utils::to_underlying(HURLEY_MCO::HeF)] /
                                       p_MassCutoffs[utils::to_underlying(HURLEY_MCO::FGB)]);
    b[47]  = (1.127733 * p_Rho) + (0.2344416 * rho2) - (0.3793726 * rho3);
    b[51] -= 0.1343798 * zeta5;
    b[53] += 0.4426929 * zeta5;
    b[55]  = min((0.99164 - (743.123 * PPOW(p_Z, 2.83))), b[55]);
    b[56] += 0.1140142 * zeta5;
    b[57] -= 0.01308728 * zeta5;

    return b;
}


/*
 * CalculateAlphaConstants_Hurley2000
 *
 * @brief
 * Calculate Hurley et al. 2000 alpha constants.
 *
 *
 * HurleyAlphasT CalculateAlphaConstants_Hurley2000(const double p_MHeF, const HurleyACoeffsT& p_aCoeffs, const HurleyBCoeffsT& p_bCoeffs) const
 *
 * @param       p_MHeF                          Maximum initial mass for which helium ignites degenerately in a Helium Flash (HeF)
 * @param       p_aCoeffs                       Hurley a coefficients (must be fully populated)
 * @param       p_bCoeffs                       Hurley b coefficients (must be fully populated)
 * @return                                      alphas array
 */
HurleyAlphasT ZDependent::CalculateAlphaConstants_Hurley2000(const double p_MHeF, const HurleyACoeffsT& p_aCoeffs, const HurleyBCoeffsT& p_bCoeffs) const {

    HurleyAlphasT alphas;

    const auto& b = p_bCoeffs;

    // alpha1
    const double lHeI_MHeF = (b[11] + (b[12] * PPOW(p_MHeF, 3.8))) / (b[13] + (p_MHeF * p_MHeF));
    alphas[0] = ((b[9] * PPOW(p_MHeF, b[10])) - lHeI_MHeF) / lHeI_MHeF;

    // alpha3
    const double lBAGB = (b[31] + (b[32] * PPOW(p_MHeF, (b[33] + 1.8)))) / (b[34] + PPOW(p_MHeF, b[33]));
    alphas[1] = ((b[29] * PPOW(p_MHeF, b[30])) - lBAGB) / lBAGB;

    // alpha4
    const double tBGB_MHeF = astro::CalculateLifetimeToBGB_Hurley2000(p_MHeF, p_aCoeffs);
    const double mHeF5     = utils::IntPow(p_MHeF, 5);
    const double tHe_MHeF  = tBGB_MHeF * (b[41] * PPOW(p_MHeF, b[42]) + b[43] * mHeF5) / (b[44] + mHeF5);
    alphas[2] = ((tHe_MHeF - b[39]) / b[39]);

    return alphas;
}


/*
 * CalculateGammaConstants_Hurley2000
 *
 * @brief
 * Calculate Hurley et al. 2000 gamma constants (eq 23 and discussion following).
 *
 *
 * HurleyGConstsT CalculateGammaConstants_Hurley2000(const HurleyACoeffsT& p_aCoeffs) const
 *
 * @param       p_aCoeffs                       Hurley a coefficients (must be fully populated)
 * @return                                      Gamma constants array
 */
HurleyGConstsT ZDependent::CalculateGammaConstants_Hurley2000(const HurleyACoeffsT& p_aCoeffs) const {


    HurleyGConstsT gammas;

    const auto& a = p_aCoeffs;
    const double bGamma = max(0.0, a[76] + (a[77] * PPOW((1.0 - a[78]), a[79])));

    using enum HURLEY_GAMMA_CONSTANTS;

    gammas[utils::to_underlying(B_GAMMA)] = bGamma;
    gammas[utils::to_underlying(C_GAMMA)] = a[75] <= 1.0 ? bGamma : a[80];

    return gammas;
}


/*
 * CalculateLuminosityConstants_Hurley2000
 *
 * @brief
 * Calculate Hurley et al. 2000 luminosity constants (eq 16, 19, 20).
 *
 *
 * HurleyLConstsT CalculateLuminosityConstants_Hurley2000(const HurleyACoeffsT& p_aCoeffs) const
 *
 * @param       p_aCoeffs                       Hurley a coefficients (must be fully populated)
 * @return                                      Luminosity constants array
 */
HurleyLConstsT ZDependent::CalculateLuminosityConstants_Hurley2000(const HurleyACoeffsT& p_aCoeffs) const {

    HurleyLConstsT lums;

    const auto& a = p_aCoeffs;

    using enum HURLEY_L_CONSTANTS;

    lums[utils::to_underlying(B_ALPHA_L)] = (a[45] + (a[46] * PPOW(2.0, a[48]))) / (PPOW(2.0, 0.4) + (a[47] * PPOW(2.0, 1.9))); // eq 19
    lums[utils::to_underlying(B_BETA_L)]  = max(0.0, (a[54] - (a[55] * PPOW(a[57], a[56]))));                                   // eq 20
    lums[utils::to_underlying(B_DELTA_L)] = min((a[34] / PPOW(a[33], a[35])), (a[36] / PPOW(a[33], a[37])));                    // eq 16

    return lums;
}


/*
 * CalculateMassCutoffs_Hurley2000
 *
 * @brief
 * Calculate Hurley et al. 2000 mass cutoffs:
 *   Hook: mass above which a hook appears on the MS (eq 1)
 *   HeF : maximum initial mass for which He ignites degenerately in the He Flash (eq 2)
 *   FGB : maximum mass at which He ignites degenerately on the FGB (eq 3)
 *
 *
 * HurleyMassCutoffsT CalculateMassCutoffs_Hurley2000(const double p_Z, const double p_Zeta) const
 *
 * @param       p_Z                             Metallicity
 * @param       p_Zeta                          Hurley zeta
 * @return                                      Mass cutoffs array
 */
HurleyMassCutoffsT ZDependent::CalculateMassCutoffs_Hurley2000(const double p_Z, const double p_Zeta) const {

    using enum HURLEY_MASS_CUTOFFS;

    HurleyMassCutoffsT massCutoffs;

    const double zeta2 = p_Zeta * p_Zeta;

    massCutoffs[utils::to_underlying(Hook)] = 1.0185 + (0.16015 * p_Zeta) + (0.0892 * zeta2);
    massCutoffs[utils::to_underlying(HeF)]  = 1.995 + (0.25 * p_Zeta) + (0.087 * zeta2);
    massCutoffs[utils::to_underlying(FGB)]  = 13.048 * PPOW((p_Z / ZSOL_HURLEY), 0.06) / (1.0 + (0.0012 * PPOW((ZSOL_HURLEY / p_Z), 1.27)));

    return massCutoffs;
}


/*
 * CalculateRadiusConstants_Hurley2000
 *
 * @brief
 * Calculate Hurley et al. 2000 radius constants (eq 17, 21, 22).
 *
 *
 * HurleyRConstsT CalculateRadiusConstants_Hurley2000(const HurleyACoeffsT& p_aCoeffs) const
 *
 * @param       p_aCoeffs                       Hurley a coefficients (must be fully populated)
 * @return                                      Radius constants array
 */
HurleyRConstsT ZDependent::CalculateRadiusConstants_Hurley2000(const HurleyACoeffsT& p_aCoeffs) const {

    using enum HURLEY_R_CONSTANTS;

    HurleyRConstsT rads;

    const auto& a = p_aCoeffs;
    rads[utils::to_underlying(B_ALPHA_R)] = (a[58] * PPOW(a[66], a[60])) / (a[59] + PPOW(a[66], a[61]));                // eq 21 (arxiv version wrong)
    rads[utils::to_underlying(C_ALPHA_R)] = (a[58] * PPOW(a[67], a[60])) / (a[59] + PPOW(a[67], a[61]));                // eq 21 (arxiv version wrong)
    rads[utils::to_underlying(B_BETA_R)]  = (a[69] * 8.0 * M_SQRT2) / (a[70] + PPOW(2.0, a[71]));                       // eq 22
    rads[utils::to_underlying(C_BETA_R)]  = (a[69] * 16384.0) / (a[70] + PPOW(16.0, a[71]));                            // eq 22
    rads[utils::to_underlying(B_DELTA_R)] = (a[38] + a[39] * 8.0 * M_SQRT2) / (a[40] * 8.0 + PPOW(2.0, a[41])) - 1.0;   // eq 17

    return rads;
}


/*
 * CalculateGBRadiusXExponent_Hurley2000
 *
 * @brief
 * Calculate the Giant Branch radius 'x' exponent (Hurley et al. 2000 eq 47),
 * a hybrid of b5 and b7. Depends on metallicity only.
 *
 *
 * double CalculateGBRadiusXExponent_Hurley2000(const double p_Zeta) const
 *
 * @param       p_Zeta                          Hurley zeta
 * @return                                      Giant Branch radius 'x' exponent
 */
double ZDependent::CalculateGBRadiusXExponent_Hurley2000(const double p_Zeta) const {

    const double zeta2 = p_Zeta * p_Zeta;
    const double zeta3 = p_Zeta * zeta2;
    const double zeta4 = p_Zeta * zeta3;

    return 0.30406 + (0.0805 * p_Zeta) + (0.0897 * zeta2) + (0.0878 * zeta3) + (0.0222 * zeta4);
}


////////////////////////////////////////////////////////////////////////////////////////
//                                                                                    //
//                            SHIKAUCHI 2024 COEFFICIENTS                             //
//                                                                                    //
////////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateCoeffs_Shikauchi2024
 *
 * @brief
 * Calculate metallicity-dependent coefficients per Shikauchi et al. 2024.
 *
 * The coefficient sets are tabulated at three bracket points (0.1 Zsol,
 * 1/3 Zsol, Zsol). For Z between brackets the coefficients are interpolated
 * linearly in log Z; outside the brackets they are clamped to the nearest one.
 *
 *
 * ShikauchiZDependentT CalculateCoeffs_Shikauchi2024(const double p_Z, const double p_logZ) const
 *
 * @param       p_Z                             Metallicity
 * @param       p_logZ                          log10(p_Z)
 * @return                                      ShikauchiZDependentT (alpha, fMix, luminosity arrays)
 */
ShikauchiZDependentT ZDependent::CalculateCoeffs_Shikauchi2024(const double p_Z, const double p_logZ) const {

    // The SHIKAUCHI_*_COEFF tables are std::vector<DblVectorT>.
    // const access is via .at(). These are references into the (const global) tables - no copies.
    const DblVectorT& alphaLower = SHIKAUCHI_ALPHA_COEFFICIENTS.at(static_cast<SizeT>(SHIKAUCHI_Coeff::ONE_TENTH_Z_SOL));
    const DblVectorT& alphaMid   = SHIKAUCHI_ALPHA_COEFFICIENTS.at(static_cast<SizeT>(SHIKAUCHI_Coeff::ONE_THIRD_Z_SOL));
    const DblVectorT& alphaUpper = SHIKAUCHI_ALPHA_COEFFICIENTS.at(static_cast<SizeT>(SHIKAUCHI_Coeff::Z_SOL));

    const DblVectorT& fMixLower  = SHIKAUCHI_FMIX_COEFFICIENTS.at(static_cast<SizeT>(SHIKAUCHI_Coeff::ONE_TENTH_Z_SOL));
    const DblVectorT& fMixMid    = SHIKAUCHI_FMIX_COEFFICIENTS.at(static_cast<SizeT>(SHIKAUCHI_Coeff::ONE_THIRD_Z_SOL));
    const DblVectorT& fMixUpper  = SHIKAUCHI_FMIX_COEFFICIENTS.at(static_cast<SizeT>(SHIKAUCHI_Coeff::Z_SOL));

    const DblVectorT& lumLower   = SHIKAUCHI_L_COEFFICIENTS.at(static_cast<SizeT>(SHIKAUCHI_Coeff::ONE_TENTH_Z_SOL));
    const DblVectorT& lumMid     = SHIKAUCHI_L_COEFFICIENTS.at(static_cast<SizeT>(SHIKAUCHI_Coeff::ONE_THIRD_Z_SOL));
    const DblVectorT& lumUpper   = SHIKAUCHI_L_COEFFICIENTS.at(static_cast<SizeT>(SHIKAUCHI_Coeff::Z_SOL));

    // Return value vectors, sized from the bracket tables.
    ShikauchiACoeffsT alphaCoeff{};
    ShikauchiFCoeffsT fMixCoeff{};
    ShikauchiLCoeffsT luminosityCoeff{};

    // Bracket points in log Z.
    const double logLowerZ = std::log10(0.1 * ZSOL_HURLEY);          // SHIKAUCHI_Coeff::ONE_TENTH_Z_SOL
    const double logMidZ   = std::log10((1.0 / 3.0) * ZSOL_HURLEY);  // SHIKAUCHI_Coeff::ONE_THIRD_Z_SOL
    const double logUpperZ = LOG10_ZSOL_HURLEY;                      // SHIKAUCHI_Coeff::Z_SOL

    const double logZ_logLowerZ    = p_logZ    - logLowerZ;
    const double logZ_logMidZ      = p_logZ    - logMidZ;
    const double logMidZ_logZ      = logMidZ   - p_logZ;
    const double logMidZ_logLowerZ = logMidZ   - logLowerZ;
    const double logUpperZ_logZ    = logUpperZ - p_logZ;
    const double logUpperZ_logMidZ = logUpperZ - logMidZ;

    if (p_logZ <= logLowerZ) {                                       // At or below lower bound -> clamp to lower
        std::copy(alphaLower.begin(), alphaLower.end(), alphaCoeff.begin());
        std::copy(fMixLower.begin(), fMixLower.end(), fMixCoeff.begin());
        std::copy(lumLower.begin(), lumLower.end(), luminosityCoeff.begin());
    }
    else if (p_logZ <= logMidZ) {                                    // Lower-to-mid band -> interpolate
        for (SizeT i = 0; i < alphaCoeff.size(); i++)
            alphaCoeff[i] = (alphaLower[i] * logMidZ_logZ + alphaMid[i] * logZ_logLowerZ) / logMidZ_logLowerZ;
        for (SizeT i = 0; i < fMixCoeff.size(); i++)
            fMixCoeff[i] = (fMixLower[i]  * logMidZ_logZ + fMixMid[i]  * logZ_logLowerZ) / logMidZ_logLowerZ;
        for (SizeT i = 0; i < luminosityCoeff.size(); i++)
            luminosityCoeff[i] = (lumLower[i]   * logMidZ_logZ + lumMid[i]   * logZ_logLowerZ) / logMidZ_logLowerZ;
    }
    else if (p_logZ < logUpperZ) {                                   // Mid-to-upper band -> interpolate
        for (SizeT i = 0; i < alphaCoeff.size(); i++)
            alphaCoeff[i] = (alphaMid[i] * logUpperZ_logZ + alphaUpper[i] * logZ_logMidZ) / logUpperZ_logMidZ;
        for (SizeT i = 0; i < fMixCoeff.size(); i++)
            fMixCoeff[i] = (fMixMid[i]  * logUpperZ_logZ + fMixUpper[i]  * logZ_logMidZ) / logUpperZ_logMidZ;
        for (SizeT i = 0; i < luminosityCoeff.size(); i++)
            luminosityCoeff[i] = (lumMid[i]   * logUpperZ_logZ + lumUpper[i]   * logZ_logMidZ) / logUpperZ_logMidZ;
    }
    else {                                                           // At or above upper bound -> clamp to upper
        std::copy(alphaUpper.begin(), alphaUpper.end(), alphaCoeff.begin());
        std::copy(fMixUpper.begin(), fMixUpper.end(), fMixCoeff.begin());
        std::copy(lumUpper.begin(), lumUpper.end(), luminosityCoeff.begin());
    }

    return ShikauchiZDependentT{ alphaCoeff, fMixCoeff, luminosityCoeff };
}
