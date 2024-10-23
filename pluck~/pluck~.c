/**
	@file
	pluck~:
	@ingroup examples
*/

#include "ext.h"			// standard Max include, always required (except in Jitter)
#include "ext_obex.h"		// required for "new" style objects
#include "z_dsp.h"			// required for MSP objects

// calculate a slope for control rate interpolation to audio rate.
#define CALCSLOPE(next, prev) ((next - prev) / sampleframes)
#define LOG001 -6.907755278982137 /* log(0.001) */
#define MIN_DELAY_SAMPLES 1.0  /* Minimum delay in samples - adjust as needed */

// struct to represent the object's state
typedef struct _pluck {
	t_pxobject		ob;			// the object itself (t_pxobject in MSP instead of t_object)
    float* m_dlybuf;
    long m_mask;
    long m_iwrphase;
    float m_lastsamp;
    float m_prevtrig;
    float m_delaytime;
    float m_decaytime;
    float m_dsamp;
    float m_feedbk;
    float m_coef;
    unsigned long m_inputsamps;
    float m_maxdelaytime;
    long m_sr;
} t_pluck;


// method prototypes
void *pluck_new(t_symbol *s, long argc, t_atom *argv);
void pluck_free(t_pluck *x);
void pluck_assist(t_pluck *x, void *b, long m, long a, char *s);
//void pluck_float(t_pluck *x, double f);
void pluck_dsp64(t_pluck *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags);
void pluck_perform64(t_pluck *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam);
/* Forward declarations of all functions */
static inline double clip(double x, double lo, double hi);
static inline double calc_delay(double delaytime, double fdelaylen);
static inline double calc_feedback(double delaytime, double decaytime);



// global class pointer variable
static t_class *pluck_class = NULL;

//***********************************************************************************************

// Utility functions
static float zapgremlins(float x) {
    float absx = fabsf(x);
    return (absx > 1e-15f && absx < 1e15f) ? x : 0.f;
}

/* Implementation of utility functions */
static inline double clip(double x, double lo, double hi) {
    return x < lo ? lo : (x > hi ? hi : x);
}

static inline double calc_delay(double delaytime, double fdelaylen) {
    double next_dsamp = delaytime * sys_getsr();  /* sys_getsr() gets Max's sample rate */
    return clip(next_dsamp, MIN_DELAY_SAMPLES, fdelaylen);
}

//static float calc_feedback(float delaytime, float decaytime) {
//    return decaytime == 0.f ? 0.f : powf(0.001f, delaytime/decaytime);
//}
static inline double calc_feedback(double delaytime, double decaytime) {
    if (delaytime == 0.0 || decaytime == 0.0)
        return 0.0;
    double absret = exp(LOG001 * delaytime / fabs(decaytime));
    /* Use signbit to handle the sign - C equivalent of copysign */
    return signbit(decaytime) ? -absret : absret;
}

static float cubicinterp(float frac, float d0, float d1, float d2, float d3) {
    float a = 0.5f * (d3 - d0) + 1.5f * (d1 - d2);
    float b = d0 - 2.5f * d1 + 2.f * d2 - 0.5f * d3;
    float c = 0.5f * (d2 - d0);
    float d = d1;
    return ((a * frac + b) * frac + c) * frac + d;
}

//***********************************************************************************************

void ext_main(void *r)
{
	// object initialization, note the use of dsp_free for the freemethod, which is required
	// unless you need to free allocated memory, in which case you should call dsp_free from
	// your custom free function.

	t_class *c = class_new("pluck~", (method)pluck_new, (method)dsp_free, (long)sizeof(t_pluck), 0L, A_GIMME, 0);

//	class_addmethod(c, (method)pluck_float,		"float",	A_FLOAT, 0);
	class_addmethod(c, (method)pluck_dsp64,		"dsp64",	A_CANT, 0);
	class_addmethod(c, (method)pluck_assist,	"assist",	A_CANT, 0);

	class_dspinit(c);
	class_register(CLASS_BOX, c);
	pluck_class = c;
}


void *pluck_new(t_symbol *s, long argc, t_atom *argv)
{
	t_pluck *x = (t_pluck *)object_alloc(pluck_class);

	if (x) {
		dsp_setup((t_pxobject *)x, 6);	// MSP inlets: arg is # of inlets and is REQUIRED!
		// use 0 if you don't need inlets
		outlet_new(x, "signal"); 		// signal outlet (note "signal" rather than NULL)
//		x->offset = 0.0;
        // Initialize parameters
        x->m_maxdelaytime = 20.0f; // default max delay time
        x->m_delaytime = 1.0f;     // default delay time
        x->m_decaytime = 1.0f;     // default decay time
        x->m_coef = 0.5f;          // default coefficient
        
        // Allocate delay buffer
        long size = (long)(x->m_maxdelaytime * sys_getsr() + 2.f);
        long realsize = 1;
        while (realsize < size) realsize <<= 1;
        x->m_mask = realsize - 1;
        x->m_dlybuf = (float*)sysmem_newptrclear(realsize * sizeof(float));
        
        if (!x->m_dlybuf) {
            object_free(x);
            return NULL;
        }
        
        // Initialize state variables
        x->m_iwrphase = 0;
        x->m_lastsamp = 0.f;
        x->m_prevtrig = 0.f;
        x->m_inputsamps = 0;
        x->m_sr = sys_getsr();
	}
	return (x);
}


// NOT CALLED!, we use dsp_free for a generic free function
void pluck_free(t_pluck *x)
{
    dsp_free((t_pxobject*)x);
    if (x->m_dlybuf) sysmem_freeptr(x->m_dlybuf);
}


void pluck_assist(t_pluck *x, void *b, long m, long a, char *s)
{
	if (m == ASSIST_INLET) { //inlet
        switch (a) {
            case 0: sprintf(s, "(signal) Input"); break;
            case 1: sprintf(s, "(signal) Trigger"); break;
            case 2: sprintf(s, "(signal) Max Delay Time"); break;
            case 3: sprintf(s, "(signal) Delay Time"); break;
            case 4: sprintf(s, "(signal) Decay Time"); break;
            case 5: sprintf(s, "(signal) Coefficient"); break;
        }
	}
	else {	// outlet
		sprintf(s, "I am outlet %ld", a);
	}
}


//void pluck_float(t_pluck *x, double f)
//{
////	x->offset = f;
//}


// registers a function for the signal chain in Max
void pluck_dsp64(t_pluck *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags)
{
//	post("my sample rate is: %f", samplerate);
    x->m_sr = samplerate;
	// instead of calling dsp_add(), we send the "dsp_add64" message to the object representing the dsp chain
	// the arguments passed are:
	// 1: the dsp64 object passed-in by the calling function
	// 2: the symbol of the "dsp_add64" message we are sending
	// 3: a pointer to your object
	// 4: a pointer to your 64-bit perform method
	// 5: flags to alter how the signal chain handles your object -- just pass 0
	// 6: a generic pointer that you can use to pass any additional data to your perform method

	object_method(dsp64, gensym("dsp_add64"), x, pluck_perform64, 0, NULL);
}


// this is the 64-bit perform method audio vectors
void pluck_perform64(t_pluck *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam)
{
	t_double *in = ins[0];		// we get audio for each inlet of the object from the **ins argument
	t_double *out = outs[0];	// we get audio for each outlet of the object from the **outs argument
    t_double* trig = ins[1];
    t_double delaytime = *ins[3];
    t_double decaytime = *ins[4];
    t_double* coef = ins[5];
	int n = sampleframes;
    
    // Local state
    float* dlybuf = x->m_dlybuf;
    long iwrphase = x->m_iwrphase;
    float lastsamp = x->m_lastsamp;
    float prevtrig = x->m_prevtrig;
    unsigned long inputsamps = x->m_inputsamps;
    float dsamp = x->m_dsamp;
    float feedbk = x->m_feedbk;
    long mask = x->m_mask;
    if (delaytime == x->m_delaytime && decaytime == x->m_decaytime) {
        // Static delay time case
        long idsamp = (long)dsamp;
        float frac = dsamp - idsamp;
        for (int i = 0; i < sampleframes; i++) {
           float curtrig = trig[i];
           if (prevtrig <= 0.f && curtrig > 0.f) {
               inputsamps = (long)(delaytime * x->m_sr + 0.5f);
           }
           prevtrig = curtrig;
           
           long irdphase1 = iwrphase - idsamp;
           long irdphase2 = irdphase1 - 1;
           long irdphase3 = irdphase1 - 2;
           long irdphase0 = irdphase1 + 1;
           
           float thisin = (inputsamps > 0) ? in[i] : 0.f;
           if (inputsamps > 0) --inputsamps;
           
           float d0 = dlybuf[irdphase0 & mask];
           float d1 = dlybuf[irdphase1 & mask];
           float d2 = dlybuf[irdphase2 & mask];
           float d3 = dlybuf[irdphase3 & mask];
           
           float value = cubicinterp(frac, d0, d1, d2, d3);
           float thiscoef = coef[i];
           float onepole = ((1.f - fabsf(thiscoef)) * value) + (thiscoef * lastsamp);
           
           dlybuf[iwrphase & mask] = thisin + feedbk * onepole;
           out[i] = lastsamp = onepole;
           
           iwrphase++;
       }
    } else {
        // Interpolating delay time case
//        calc_delay(delaytime, x->m_maxdelaytime);
        float next_dsamp = delaytime * x->m_sr;
        float dsamp_slope = CALCSLOPE(next_dsamp, dsamp);
        
        float next_feedbk = calc_feedback(delaytime, decaytime);
        float feedbk_slope = CALCSLOPE(next_feedbk, feedbk);
        
        for (int i = 0; i < sampleframes; i++) {
            float curtrig = trig[i];
            if (prevtrig <= 0.f && curtrig > 0.f) {
                inputsamps = (long)(delaytime * x->m_sr + 0.5f);
            }
            prevtrig = curtrig;
            
            dsamp += dsamp_slope;
            long idsamp = (long)dsamp;
            float frac = dsamp - idsamp;
            
            long irdphase1 = iwrphase - idsamp;
            long irdphase2 = irdphase1 - 1;
            long irdphase3 = irdphase1 - 2;
            long irdphase0 = irdphase1 + 1;
            
            float thisin = (inputsamps > 0) ? in[i] : 0.f;
            if (inputsamps > 0) --inputsamps;
            
            float d0 = dlybuf[irdphase0 & mask];
            float d1 = dlybuf[irdphase1 & mask];
            float d2 = dlybuf[irdphase2 & mask];
            float d3 = dlybuf[irdphase3 & mask];
            
            float value = cubicinterp(frac, d0, d1, d2, d3);
            float thiscoef = coef[i];
            float onepole = ((1.f - fabsf(thiscoef)) * value) + (thiscoef * lastsamp);
            
            dlybuf[iwrphase & mask] = thisin + feedbk * onepole;
            out[i] = lastsamp = onepole;
            
            feedbk += feedbk_slope;
            iwrphase++;
        }
        
        // Store updated parameters
        x->m_delaytime = delaytime;
        x->m_decaytime = decaytime;
        x->m_feedbk = feedbk;
        x->m_dsamp = dsamp;
    }
    x->m_iwrphase = iwrphase;
    x->m_lastsamp = zapgremlins(lastsamp);
    x->m_prevtrig = prevtrig;
    x->m_inputsamps = inputsamps;
            // this perform method simply copies the input to the output, offsetting the value
//            while (n--)
//                *out++ = *in++;
}

