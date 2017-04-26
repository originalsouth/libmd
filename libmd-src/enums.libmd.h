#define __libmd_h_file__

#ifndef enums_libmd_h
#define enums_libmd_h
struct INTEGRATOR {enum integrator:uc {SEULER,VVERLET,FO,FO_OVERDAMPED};};      ///< Integration options
struct MP_INTEGRATOR {enum mp_integrator:uc {VZ,VZ_P,VZ_WFI,SEULER,VVERLET};};  ///< Monge patch integration options
struct BCOND {enum bcond:uc {NONE,PERIODIC,HARD,BOXSHEAR};};                    ///< Boundary condition options
struct INDEX {enum index:uc {CELL,BRUTE_FORCE,KD_TREE};};                       ///< Indexing options
struct POT {enum pot:ui                                                         ///< Potential options
{
    COULOMB,
    YUKAWA,
    HOOKEAN,
    LJ,
    MORSE,
    FORCEDIPOLE,
    HOOKEANFORCEDIPOLE,
    ANHARMONICSPRING
};};
/// External force options
struct EXTFORCE {enum extforce:ui
{
    DAMPING,
    DISSIPATION,
    LANGEVIN
};};
/// Monge patch options
struct MP {enum mp:ui
{
    FLATSPACE,
    GAUSSIANBUMP,
    EGGCARTON,
    MOLLIFIER
};};
#endif
