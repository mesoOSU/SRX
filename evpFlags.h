#ifndef EVPFLAGS_H
#define EVPFLAGS_H

#define NPHMAX 2	//maximum number of phases allowed
#define NSYSMX 12	// maximum number of active slip+twin systems allowed in any phase

#define DD_BASED_FLAG	// dislocation density based constitutive laws
//#define DD_CU		// used to check power law for Copper case which use T-dependent elastic constants

// In dislocation density based model, one can consider
// both forward and backward jump or only forward
// jump (at low temperature)
#ifdef DD_BASED_FLAG
//#define DD_BCC
//#define DD_AL		// for Aluminium case which use T-dependent elastic constants
#define DD_CU		// for Copper case which use T-dependent elastic constants
//#define FORWARD_ONLY
#define FORWARD_PLUS_BACKWARD

/* Using modified Newton-Raphson method */
//#define NR_MODIFIED

/* use power law instead of exp() */
//#define DD_POWER_LAW
//
/* initially soften some artibrary grains */
//#define COMPOSITE_CHECK

//#define DD_GND
#ifdef DD_GND
#define DD_GND_FIRST_ORDER
#endif

#define PF_DRX	// simulate DRX
#ifdef PF_DRX
#define GB_BULGING_ONLY	// only allow GB bulging?
//#define NUCL_RHO_DIFF	// use disl. density difference to check nucleation
//#define LOCAL_RELAX_FIRST
//#define LOCAL_RELAX_SECOND
//#define DRX_ELASTIC_TREAT
#endif


#endif

#endif
