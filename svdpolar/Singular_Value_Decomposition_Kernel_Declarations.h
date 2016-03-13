//###########################################################
// Local variable declarations
//###########################################################
const float Four_Gamma_Squared=sqrt(8.)+3.;
const float Sine_Pi_Over_Eight=.5*sqrt(2.-sqrt(2.));
const float Cosine_Pi_Over_Eight=.5*sqrt(2.+sqrt(2.));

ENABLE_SCALAR_IMPLEMENTATION(union {float f;unsigned int ui;} Sfour_gamma_squared;)
ENABLE_SCALAR_IMPLEMENTATION(union {float f;unsigned int ui;} Ssine_pi_over_eight;)
ENABLE_SCALAR_IMPLEMENTATION(union {float f;unsigned int ui;} Scosine_pi_over_eight;)
ENABLE_SCALAR_IMPLEMENTATION(union {float f;unsigned int ui;} Sone_half;)
ENABLE_SCALAR_IMPLEMENTATION(union {float f;unsigned int ui;} Sone;)
ENABLE_SCALAR_IMPLEMENTATION(union {float f;unsigned int ui;} Stiny_number;)
ENABLE_SCALAR_IMPLEMENTATION(union {float f;unsigned int ui;} Ssmall_number;)

ENABLE_SCALAR_IMPLEMENTATION(Sfour_gamma_squared.f=Four_Gamma_Squared;)
ENABLE_SCALAR_IMPLEMENTATION(Ssine_pi_over_eight.f=Sine_Pi_Over_Eight;)
ENABLE_SCALAR_IMPLEMENTATION(Scosine_pi_over_eight.f=Cosine_Pi_Over_Eight;)
ENABLE_SCALAR_IMPLEMENTATION(Sone_half.f=.5;)
ENABLE_SCALAR_IMPLEMENTATION(Sone.f=1.;)
ENABLE_SCALAR_IMPLEMENTATION(Stiny_number.f=1.e-20;)
ENABLE_SCALAR_IMPLEMENTATION(Ssmall_number.f=1.e-12;)

ENABLE_SCALAR_IMPLEMENTATION(union {float f;unsigned int ui;} Sa11;)
ENABLE_SCALAR_IMPLEMENTATION(union {float f;unsigned int ui;} Sa21;)
ENABLE_SCALAR_IMPLEMENTATION(union {float f;unsigned int ui;} Sa31;)
ENABLE_SCALAR_IMPLEMENTATION(union {float f;unsigned int ui;} Sa12;)
ENABLE_SCALAR_IMPLEMENTATION(union {float f;unsigned int ui;} Sa22;)
ENABLE_SCALAR_IMPLEMENTATION(union {float f;unsigned int ui;} Sa32;)
ENABLE_SCALAR_IMPLEMENTATION(union {float f;unsigned int ui;} Sa13;)
ENABLE_SCALAR_IMPLEMENTATION(union {float f;unsigned int ui;} Sa23;)
ENABLE_SCALAR_IMPLEMENTATION(union {float f;unsigned int ui;} Sa33;)

#ifdef COMPUTE_V_AS_MATRIX
ENABLE_SCALAR_IMPLEMENTATION(union {float f;unsigned int ui;} Sv11;)
ENABLE_SCALAR_IMPLEMENTATION(union {float f;unsigned int ui;} Sv21;)
ENABLE_SCALAR_IMPLEMENTATION(union {float f;unsigned int ui;} Sv31;)
ENABLE_SCALAR_IMPLEMENTATION(union {float f;unsigned int ui;} Sv12;)
ENABLE_SCALAR_IMPLEMENTATION(union {float f;unsigned int ui;} Sv22;)
ENABLE_SCALAR_IMPLEMENTATION(union {float f;unsigned int ui;} Sv32;)
ENABLE_SCALAR_IMPLEMENTATION(union {float f;unsigned int ui;} Sv13;)
ENABLE_SCALAR_IMPLEMENTATION(union {float f;unsigned int ui;} Sv23;)
ENABLE_SCALAR_IMPLEMENTATION(union {float f;unsigned int ui;} Sv33;)
#endif

#ifdef COMPUTE_U_AS_MATRIX
ENABLE_SCALAR_IMPLEMENTATION(union {float f;unsigned int ui;} Su11;)
ENABLE_SCALAR_IMPLEMENTATION(union {float f;unsigned int ui;} Su21;)
ENABLE_SCALAR_IMPLEMENTATION(union {float f;unsigned int ui;} Su31;)
ENABLE_SCALAR_IMPLEMENTATION(union {float f;unsigned int ui;} Su12;)
ENABLE_SCALAR_IMPLEMENTATION(union {float f;unsigned int ui;} Su22;)
ENABLE_SCALAR_IMPLEMENTATION(union {float f;unsigned int ui;} Su32;)
ENABLE_SCALAR_IMPLEMENTATION(union {float f;unsigned int ui;} Su13;)
ENABLE_SCALAR_IMPLEMENTATION(union {float f;unsigned int ui;} Su23;)
ENABLE_SCALAR_IMPLEMENTATION(union {float f;unsigned int ui;} Su33;)
#endif

ENABLE_SCALAR_IMPLEMENTATION(union {float f;unsigned int ui;} Sc;)
ENABLE_SCALAR_IMPLEMENTATION(union {float f;unsigned int ui;} Ss;)
ENABLE_SCALAR_IMPLEMENTATION(union {float f;unsigned int ui;} Sch;)
ENABLE_SCALAR_IMPLEMENTATION(union {float f;unsigned int ui;} Ssh;)
ENABLE_SCALAR_IMPLEMENTATION(union {float f;unsigned int ui;} Stmp1;)
ENABLE_SCALAR_IMPLEMENTATION(union {float f;unsigned int ui;} Stmp2;)
ENABLE_SCALAR_IMPLEMENTATION(union {float f;unsigned int ui;} Stmp3;)
ENABLE_SCALAR_IMPLEMENTATION(union {float f;unsigned int ui;} Stmp4;)
ENABLE_SCALAR_IMPLEMENTATION(union {float f;unsigned int ui;} Stmp5;)
