/**
	@file
	gverb~:  A two-channel reverb, based on the "GVerb" LADSPA effect by Juhana Sadeharju (kouhia at nic.funet.fi).
	@ingroup examples
*/


#include "ext.h"			// standard Max include, always required (except in Jitter)
#include "ext_obex.h"		// required for "new" style objects
#include "z_dsp.h"			// required for MSP objects
// calculate a slope for control rate interpolation to audio rate.
//#define CALCSLOPE(next, prev) ((next - prev) / sampleframes)
#define LOG001 -6.907755278982137 /* log(0.001) */
#define MIN_DELAY_SAMPLES 1.0  /* Minimum delay in samples - adjust as needed */
#define sc_max(x,y) (((x) >= (y)) ? (x) : (y))
#define sc_min(x,y) (((x) <= (y)) ? (x) : (y))
#define INT_MAX 2147483647

// Memory failure handling
#define ClearUnitOnMemFailed(x)                                                                                      \
    object_error((t_object *)x, "Memory allocation failed, increase Max's buffer size or system memory!");            \
    clear_outputs(x, 0);                                                                                              \
    x->done = true;                                                                                                   \
    return;

#define ClearUnitIfMemFailed(condition, x)                                                                            \
    if (!(condition)) {                                                                                               \
        ClearUnitOnMemFailed(x);                                                                                      \
    }

/*  GVerb work */
#define FDNORDER 4

//***********************************************************************************************
// Struct definition

typedef union {
    float f;
    int32_t i;
} ls_pcast32;

typedef struct {
    int size;
    int idx;
    float* buf;
} g_fixeddelay;

typedef struct {
    int size;
    float coef;
    int idx;
    float* buf;
} g_diffuser;

typedef struct {
    float damping;
    float delay;
} g_damper;


//***********************************************************************************************
// struct to represent the object's state

typedef struct _gverb {
	t_pxobject		ob;			// the object itself (t_pxobject in MSP instead of t_object)
    float *outputs;     // Buffer for outputs
    int num_outputs;    // Number of output channels
    double alpha;
    float roomsize, revtime, damping, spread, inputbandwidth, drylevel, earlylevel, taillevel;
    float maxroomsize;
    float maxdelay, largestdelay;
    //slope
    float earlylevelslope, taillevelslope, drylevelslope;
    float fdngainslopes[FDNORDER], tapgainslopes[FDNORDER];
    //FDN
    g_fixeddelay* fdndels[FDNORDER];
    g_damper* fdndamps[FDNORDER];
    float fdngains[FDNORDER];
    int fdnlens[FDNORDER];
    //diffuser variables
    g_diffuser* ldifs[FDNORDER];
    g_diffuser* rdifs[FDNORDER];
    //tap variables
    g_fixeddelay* tapdelay;
    int taps[FDNORDER];
    float tapgains[FDNORDER];
    //Processing buffers
    float u[FDNORDER], f[FDNORDER], d[FDNORDER];
    
    g_damper* inputdamper;
    
    bool done;          // Indicates if the object is done
} t_gverb;

//***********************************************************************************************
// method prototypes

void *gverb_new(t_symbol *s, long argc, t_atom *argv);
void gverb_free(t_gverb *x);
void gverb_assist(t_gverb *x, void *b, long m, long a, char *s);
//void gverb_float(t_gverb *x, double f);
void gverb_dsp64(t_gverb *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags);
void gverb_perform64(t_gverb *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam);
/* Forward declarations of all functions */
static g_damper* make_damper(t_gverb *x, float damping);
static g_diffuser* make_diffuser(t_gverb *x, int size, float coef);
static g_fixeddelay* make_fixeddelay(t_gverb *x, int size, int maxsize);
// global class pointer variable
static t_class *gverb_class = NULL;

//***********************************************************************************************

// Function to clear all outputs (similar to Unit_ZeroOutputs)
void clear_outputs(t_gverb *x, int num_samples) {
    for (int i = 0; i < x->num_outputs; i++) {
        if (x->outputs) {
            // Set output buffer to zero
            memset(x->outputs, 0, num_samples * sizeof(float));
        }
    }
}

// Utility functions
static float zapgremlins(float x) {
    float absx = fabsf(x);
    return (absx > 1e-15f && absx < 1e15f) ? x : 0.f;
}

// Helper functions from original code, converted to static inline
static inline float flush_to_zero(float f) {
    ls_pcast32 v;
    v.f = f;
    return (v.i & 0x7f800000) < 0x08000000 ? 0.0f : f;
}

static inline int isprime(int n) {
    unsigned int i;
    const unsigned int lim = (int)sqrtf((float)n);

    if (n == 2)
        return TRUE;
    if ((n & 1) == 0)
        return FALSE;
    for (i = 3; i <= lim; i += 2)
        if ((n % i) == 0)
            return FALSE;
    return TRUE;
}

static inline int nearestprime(int n, float rerror) {
    int bound, k;

    if (isprime(n))
        return n;
    bound = n * (int)rerror;
    for (k = 1; k <= bound; k++) {
        if (isprime(n + k))
            return (n + k);
        if (isprime(n - k))
            return (n - k);
    }
    return -1;
}
//changed here from sys_getsr() to sys_getblksize()
//    sys_getblksize()
static inline float calc_slope(float next, float prev) {
    return ((next - prev) /sys_getblksize());
}

static inline int f_round(float f) {
    ls_pcast32 p;
    p.f = f;
    p.f += (3 << 22);
    return p.i - 0x4b400000;
}

//***********************************************************************************************

static inline float damper_do(g_damper* p, float x) {
    float y;
    y = x * (1.0f - p->damping) + p->delay * p->damping;
    p->delay = zapgremlins(y);
    return (y);
}


static inline void fixeddelay_write(g_fixeddelay* p, float x) {
    p->buf[p->idx] = zapgremlins(x);
    p->idx = (p->idx + 1) % p->size;
}

static inline void diffuser_write(g_diffuser* p, float x) {
    p->buf[p->idx] = (isnan(x)) ? 0.0f : x;
    p->idx = (p->idx + 1) % p->size;
}

//static inline float diffuser_do(g_diffuser* p, float x) {
//    float y, w;
//    
//    // Add coefficient scaling based on delay size
//    float size_scale = fmin(1.0f, 1000.0f / (float)p->size);
//    float adjusted_coef = p->coef * size_scale;
//    
//    // Add soft clipping to feedback path
//    float fb = p->buf[p->idx];
//    fb = (fb > 1.0f) ? (2.0f - 1.0f/fb) : ((fb < -1.0f) ? (-2.0f - 1.0f/fb) : fb);
//    
//    w = x - fb * adjusted_coef;
//    w = flush_to_zero(w);
//    
//    // Soft clip the output
//    y = p->buf[p->idx] + w * adjusted_coef;
//    y = (y > 1.0f) ? (2.0f - 1.0f/y) : ((y < -1.0f) ? (-2.0f - 1.0f/y) : y);
//    
//    p->buf[p->idx] = zapgremlins(w);
//    p->idx = (p->idx + 1) % p->size;
//    
//    return y;
//}
static inline float diffuser_do(g_diffuser* p, float x) {
    float y, w;
    w = x - p->buf[p->idx] * p->coef; // Feedback path
    w = flush_to_zero(w);
    y = p->buf[p->idx] + w * p->coef; // Feed forward path
    p->buf[p->idx] = zapgremlins(w);
    p->idx = (p->idx + 1) % p->size;
    return (y);
}



// Function to clear the buffer (similar to the Clear function in SuperCollider)
void clear_buffer(int size, float* buf) {
    // Zero out the buffer for the given size
    memset(buf, 0, size * sizeof(float));
}

g_damper* make_damper(t_gverb *x, float damping) {
    g_damper* p;
    p = (g_damper*)sysmem_newptr(sizeof(g_damper));
    if (p) {
        p->damping = damping;
        p->delay = 0.f;
    }
    return (p);
}

g_diffuser* make_diffuser(t_gverb *x, int size, float coef) {
    g_diffuser* p;
    p = (g_diffuser*)sysmem_newptr(sizeof(g_diffuser));
//    p = (g_diffuser*)RTAlloc(x->mWorld, sizeof(g_diffuser));
    if (!p) return NULL;
    p->size = size;
    p->coef = coef;
    p->idx = 0;
    p->buf = (float*)sysmem_newptr(size * sizeof(float));
    if (p->buf == NULL) {
        sysmem_freeptr(p);
        p= NULL;
        return NULL;
    }
//    Clear(size, p->buf);
//    for (int i = 0; i < size; i++) p->buf[i] = 0.0f;
    clear_buffer(size, p->buf);

    return (p);
}

g_fixeddelay* make_fixeddelay(t_gverb *x, int size, int maxsize) {
    g_fixeddelay* p;
    
    p = (g_fixeddelay*)sysmem_newptr(sizeof(g_fixeddelay));
    if (p == NULL)
        return NULL;
    p->size = size;
    p->idx = 0;
    // Add overflow check for buffer allocation
    if (maxsize > (INT_MAX / sizeof(float))) {
        sysmem_freeptr(p);
        return NULL;
    }
    
    p->buf = (float*)sysmem_newptr(maxsize * sizeof(float));
    if (p->buf == NULL) {
        sysmem_freeptr(p);
        p = NULL;
//        RTFree(x->mWorld, p);
        return NULL;
    }
//    for (int i = 0; i < size; i++) p->buf[i] = 0.0f;
    clear_buffer(maxsize, p->buf);
    return (p);
}

void free_fixeddelay(t_gverb *x, g_fixeddelay* p) {
    if (p)
        sysmem_freeptr(p->buf);
//        RTFree(x->mWorld, p->buf);
    sysmem_freeptr(p);
//    RTFree(x->mWorld, p);
}

void free_diffuser(t_gverb *x, g_diffuser* p) {
    if (p)
        sysmem_freeptr(p->buf);
//        RTFree(x->mWorld, p->buf);
    sysmem_freeptr(p);
//    if (p)
//        RTFree(x->mWorld, p->buf);
//    RTFree(x->mWorld, p);
}

void free_damper(t_gverb* x, g_damper* p) { sysmem_freeptr(p); };

static inline float fixeddelay_read(g_fixeddelay* p, int n) {
    int i;
    i = (p->idx - n + p->size) % p->size;
    return (p->buf[i]);
}

//bad access qui per damping, che è chiamato da gverb_set_damping
static inline void damper_set(g_damper* p, t_double damping) {
    if (p != NULL && isfinite(damping)) {
        p->damping = damping;
    }
}
//
//static inline float limit_in(float x) {
//    return (x > 1.0f) ? (2.0f - 1.0f/x) : ((x < -1.0f) ? (-2.0f - 1.0f/x) : x);
//}
//
//static inline void gverb_fdnmatrix(float* a, float* b) {
//    const float dl0 = a[0], dl1 = a[1], dl2 = a[2], dl3 = a[3];
//    
//    // Calculate matrix operations with scaling
//    const float scale = 0.45f;  // Slightly reduce the matrix gain
//    b[0] = scale * limit_in(+dl0 + dl1 - dl2 - dl3);
//    b[1] = scale * limit_in(+dl0 - dl1 - dl2 + dl3);
//    b[2] = scale * limit_in(-dl0 + dl1 - dl2 + dl3);
//    b[3] = scale * limit_in(+dl0 + dl1 + dl2 + dl3);
//}

static inline void gverb_fdnmatrix(float* a, float* b) {
    const float dl0 = a[0], dl1 = a[1], dl2 = a[2], dl3 = a[3];
    b[0] = 0.5f * (+dl0 + dl1 - dl2 - dl3);
    b[1] = 0.5f * (+dl0 - dl1 - dl2 + dl3);
    b[2] = 0.5f * (-dl0 + dl1 - dl2 + dl3);
    b[3] = 0.5f * (+dl0 + dl1 + dl2 + dl3);
}
//
//static inline void damper_clear(g_damper* p) {
//    // Reset any internal state of your damper
//    p->damping = 0.0f;  // Assuming your damper has a 'value' field
//    p->delay = 0.0f;
//    // Add any other state resets needed for your damper implementation
//}

//valori che sono ok
//largestdelay 1167.352905
//fdndels[0] 62489328
//fdndels[1] 62490000
//fdndels[2] 62489392
//fdndels[3] 62489872
//fdngains[i] -1.000052
//fdngainslopes[i] 0.001024
//fdngains[i] -0.951459
//fdngainslopes[i] 0.000082
//fdngains[i] -0.957839
//fdngainslopes[i] 0.000072
//fdngains[i] -0.962200
//fdngainslopes[i] 0.000064
//taps[0] 483
//taps[1] 355
//taps[2] 185
//taps[3] 5


//valori che fanno glitch
//largestdelay 1426.764648
//fdndels[0] 62489328
//fdndels[1] 62490000
//fdndels[2] 62489392
//fdndels[3] 62489872
//fdngains[0] -1.000052
//fdngainslopes[0] 0.000000
//fdngains[1] -0.941032
//fdngainslopes[1] -0.000163
//fdngains[2] -0.948728
//fdngainslopes[2] -0.000142
//fdngains[3] -0.953994
//fdngainslopes[4] -0.000128
//taps[0] 589
//taps[1] 433
//taps[2] 226
//taps[3] 5


static inline void gverb_set_roomsize(t_gverb *x, const float a) {
    unsigned int i;

    if (a <= 1.0 || isnan(a)) {
        x->roomsize = 1.0;
    } else {
        if (a >= x->maxroomsize)
            x->roomsize = x->maxroomsize - 1.f;
        else
            x->roomsize = a;
    };
    // a roomsize = 11 e largestdelay = 1426.764648 va in tilt

    x->largestdelay = sys_getsr() * x->roomsize / 340.0; // * 0.00294f;
//    post("largestdelay[0] %f", (float)x->largestdelay);

    // the line below causes everything to blow up.... why?????
//    x->fdnlens[0] = nearestprime((int)(x->largestdelay), 0.5);
//    post("fdnlens[0] %d", (int)x->fdnlens[0]);
    x->fdnlens[1] = (int)(0.816490 * x->largestdelay);
//    post("fdnlens[1] %d", (int)x->fdnlens[1]);
    x->fdnlens[2] = (int)(0.707100 * x->largestdelay);
//    post("fdnlens[2] %d", (int)x->fdnlens[2]);
    x->fdnlens[3] = (int)(0.632450 * x->largestdelay);
//    post("fdnlens[3] %d", (int)x->fdnlens[3]);

//working
//    largestdelay[0] 1167.352905
//    fdnlens[0] 1297
//    fdnlens[1] 953
//    fdnlens[2] 825
//    fdnlens[3] 738
//    fdngains[0] -0.934522
//    fdngains[1] -0.951459
//    fdngains[2] -0.957839
//    fdngains[3] -0.962200
//    fdngainslopes[0] 0.000000
//    fdngainslopes[1] 0.000082
//    fdngainslopes[2] 0.000072
//    fdngainslopes[3] 0.000064
//    taps[0] 483
//    taps[1] 355
//    taps[2] 185
//    taps[3] 5
//    tapgains[0] 0.975097
//    tapgains[1] 0.981635
//    tapgains[2] 0.990387
//    tapgains[3] 0.999739
//    tapgainslopes[0] -0.000042
//    tapgainslopes[1] -0.000031
//    tapgainslopes[2] -0.000017
//    tapgainslopes[3] 0.000000

    
//    crash
//    largestdelay[0] 1426.764648
//    fdnlens[0] 1297
//    fdnlens[1] 1164
//    fdnlens[2] 1008
//    fdnlens[3] 902
//    fdngains[0] -0.934519
//    fdngains[1] -0.941032
//    fdngains[2] -0.948728
//    fdngains[3] -0.953994
//    fdngainslopes[0] -0.000000
//    fdngainslopes[1] -0.000163
//    fdngainslopes[2] -0.000142
//    fdngainslopes[3] -0.000128
//    taps[0] 589
//    taps[1] 433
//    taps[2] 226
//    taps[3] 5
//    tapgains[0] 0.969713
//    tapgains[1] 0.977645
//    tapgains[2] 0.988269
//    tapgains[3] 0.999739
//    tapgainslopes[0] 0.000042
//    tapgainslopes[1] 0.000031
//    tapgainslopes[2] 0.000016
//    tapgainslopes[3] 0.000000


    for (i = 0; i < FDNORDER; i++) {
        float oldfdngain = x->fdngains[i];
        x->fdngains[i] = (float)(-pow(x->alpha, (double)(x->fdnlens[i])));
        x->fdngainslopes[i] = calc_slope(x->fdngains[i], oldfdngain);
    }
//    post("fdngains[0] %f", (float)x->fdngains[0]);
//    post("fdngains[1] %f", (float)x->fdngains[1]);
//    post("fdngains[2] %f", (float)x->fdngains[2]);
//    post("fdngains[3] %f", (float)x->fdngains[3]);
//
//    post("fdngainslopes[0] %f", (float)x->fdngainslopes[0]);
//    post("fdngainslopes[1] %f", (float)x->fdngainslopes[1]);
//    post("fdngainslopes[2] %f", (float)x->fdngainslopes[2]);
//    post("fdngainslopes[3] %f", (float)x->fdngainslopes[3]);
//

    x->taps[0] = 5 + (int)(0.410 * x->largestdelay);
    x->taps[1] = 5 + (int)(0.300 * x->largestdelay);
    x->taps[2] = 5 + (int)(0.155 * x->largestdelay);
    x->taps[3] = 5; //+ f_round(0.000 * largestdelay);
//    post("taps[0] %d", (int)x->taps[0]);
//    post("taps[1] %d", (int)x->taps[1]);
//    post("taps[2] %d", (int)x->taps[2]);
//    post("taps[3] %d", (int)x->taps[3]);
    
    for (i = 0; i < FDNORDER; i++) {
        float oldtapgain = x->tapgains[i];
        x->tapgains[i] = (float)(pow(x->alpha, (double)(x->taps[i])));
//        float tap_tmp_slope =sc_max(calc_slope(x->tapgains[i], oldtapgain), 0.f);
        x->tapgainslopes[i] = calc_slope(x->tapgains[i], oldtapgain);
//        post("tapgains[i] %f",(float)x->tapgains[i]);
//        post("tapgainslopes[i] %f",(float)x->tapgainslopes[i]);
    }
//    post("tapgains[0] %f", (float)x->tapgains[0]);
//    post("tapgains[1] %f", (float)x->tapgains[1]);
//    post("tapgains[2] %f", (float)x->tapgains[2]);
//    post("tapgains[3] %f", (float)x->tapgains[3]);
//
//    post("tapgainslopes[0] %f", (float)x->tapgainslopes[0]);
//    post("tapgainslopes[1] %f", (float)x->tapgainslopes[1]);
//    post("tapgainslopes[2] %f", (float)x->tapgainslopes[2]);
//    post("tapgainslopes[3] %f", (float)x->tapgainslopes[3]);
}

static inline void gverb_set_revtime(t_gverb *x, float a) {
    float ga;
    double n;
    unsigned int i;

    x->revtime = a;
    
    ga = 0.001;
    n = sys_getsr() * a;
    x->alpha = (double)powf(ga, (float)(1.f / n));

    for (i = 0; i < FDNORDER; i++) {
        float oldfdngain = x->fdngains[i];
        x->fdngains[i] = (float)(-pow(x->alpha, (double)(x->fdnlens[i])));
        x->fdngainslopes[i] = calc_slope(x->fdngains[i], oldfdngain);
    }
}

static inline void gverb_set_damping(t_gverb *x, float a) {
    unsigned int i;
    x->damping = a;
    for (i = 0; i < FDNORDER; i++) {
        damper_set(x->fdndamps[i], x->damping);
    }
}

static inline void gverb_set_inputbandwidth(t_gverb *x, float a) {
    x->inputbandwidth = a;
    damper_set(x->inputdamper, 1.0 - x->inputbandwidth);
}

static inline float gverb_set_earlylevel(t_gverb *x, float a) {
    float oldearly = x->earlylevel;
    x->earlylevel = a;
    x->earlylevelslope = calc_slope(a, oldearly);
    return (oldearly);
}

static inline float gverb_set_taillevel(t_gverb *x, float a) {
    float oldtail = x->taillevel;
    x->taillevel = a;
    x->taillevelslope = calc_slope(a, oldtail);
    return (oldtail);
}

static inline float gverb_set_drylevel(t_gverb *x, float a) {
    float olddry = x->drylevel;
    x->drylevel = a;
    x->drylevelslope = calc_slope(a, olddry);
    return (olddry);
}


//***********************************************************************************************

void ext_main(void *r)
{
	// object initialization, note the use of dsp_free for the freemethod, which is required
	// unless you need to free allocated memory, in which case you should call dsp_free from
	// your custom free function.

	t_class *c = class_new("scgverb~", (method)gverb_new, (method)dsp_free, (long)sizeof(t_gverb), 0L, A_GIMME, 0);

//	class_addmethod(c, (method)gverb_float,		"float",	A_FLOAT, 0);
	class_addmethod(c, (method)gverb_dsp64,		"dsp64",	A_CANT, 0);
	class_addmethod(c, (method)gverb_assist,	"assist",	A_CANT, 0);
    
//    
//    CLASS_ATTR_FLOAT(c, "maxroomsize", 0, t_gverb, maxroomsize);
//    CLASS_ATTR_BASIC(c, "maxroomsize", 0);
//    CLASS_ATTR_LABEL(c, "maxroomsize", 0, "Maximum Room Size");
//    CLASS_ATTR_ALIAS(c, "maxroomsize", "roomsize");
//    
//    CLASS_ATTR_FLOAT(c, "spread", 0, t_gverb, spread);
//    CLASS_ATTR_BASIC(c, "spread", 0);
//    CLASS_ATTR_LABEL(c, "spread", 0, "CSpread Amount");
//    CLASS_ATTR_ALIAS(c, "spread", "spread");
    

    
	class_dspinit(c);
	class_register(CLASS_BOX, c);
	gverb_class = c;
}


void *gverb_new(t_symbol *s, long argc, t_atom *argv)
{
	t_gverb *x = (t_gverb *)object_alloc(gverb_class);
    
	if (x) {
		dsp_setup((t_pxobject *)x, 8);	// MSP inlets: arg is # of inlets and is REQUIRED!
		// use 0 if you don't need inlets
		outlet_new(x, "signal"); 		// signal outlet (note "signal" rather than NULL)
        outlet_new(x, "signal"); // Right output
        if (argc < 2) {
                    post("scgverb~: Warning - expecting 2 arguments:");
                    post("1. maxroomsize (float > 1.0001)");
                    post("2. spread (float > 10.0)");
                }

        //x->maxroomsize = 300.0f;  // Default max room size in meters
        float revtime = x->revtime = 3.0f;        // Default reverb time in seconds
        float damping =  x->damping = 0.5f;        // Default damping
        float inputbandwidth =x->inputbandwidth = 0.5f; // Default input bandwidth
        float maxroomsize = x->maxroomsize = sc_max(1.0001f, atom_getfloatarg(0, argc, argv));
        float roomsize = x->roomsize = 10.0f;
        float spread = x->spread = sc_max(3.0f, atom_getfloatarg(1, argc, argv));;        // Default stereo spread
        float drylevel = x->drylevel = 0.0f;       // Default dry level
        float earlylevel = x->earlylevel = 0.5f;     // Default early reflections level
        float taillevel = x->taillevel = 0.5f;      // Default tail level

        
        // when roomsize is greater than maxroomsize, it is set to maxroomsize - 1
        // when roomsize is less than 0, it is set to 1, therefore maxroomsize must be at least 1.
//        float maxroomsize = x->maxroomsize = sc_max(1.0001f, IN0(9));
        
//        x->outputs = (float *)sysmem_newptrclear(x->num_outputs * sizeof(float));  // Allocate output buffers
        
//        ClearUnitIfMemFailed(x->outputs, x);  // Check if memory allocation failed
        
        float maxdelay = x->maxdelay = sys_getsr() * maxroomsize / 340.f;
        float largestdelay = x->largestdelay = sys_getsr() * roomsize / 340.f;

        for (int i = 0; i < FDNORDER; ++i) {
            x->fdndels[i] = NULL;
            x->fdndamps[i] = NULL;
            x->ldifs[i] = NULL;
            x->rdifs[i] = NULL;
        }
        x->tapdelay = NULL;
        // make the inputdamper
        x->inputdamper = make_damper(x, 1. - inputbandwidth);
        ClearUnitIfMemFailed(x->inputdamper, x);

        // float ga = powf(10.f, -60.f/20.f);
        float ga = 0.001f;
        float n = sys_getsr() * revtime;
        double alpha = x->alpha = pow((double)ga, 1. / (double)n);
        float gbmul[4] = { 1.000, 0.816490, 0.707100, 0.632450 };
        for (int i = 0; i < FDNORDER; ++i) {
            float gb = gbmul[i] * largestdelay;
            if (i == 0) {
                x->fdnlens[i] = nearestprime((int)gb, 0.5);
            } else {
                x->fdnlens[i] = f_round(gb);
            }
            x->fdngains[i] = (float)(-pow(alpha, (double)(x->fdnlens[i])));
        }
        // make the fixeddelay lines and dampers
        
        //questa parte fa crashare tutto!
        for (int i = 0; i < FDNORDER; i++) {
//            post("Allocating fixed delay for FDNORDER index: %d", i);
            x->fdndels[i] = make_fixeddelay(x, (int)x->fdnlens[i], (int)maxdelay + 1000);
            if (!x->fdndels[i]) {
                object_error((t_object*)x, "Failed to allocate fdndels for index %d", i);
                return;
            }
            
//            post("Allocating damper for FDNORDER index: %d", i);
            x->fdndamps[i] = make_damper(x, damping);
            if (!x->fdndamps[i]) {
                object_error((t_object*)x, "Failed to allocate fdndamps for index %d", i);
                return;
            }
            
            // If both allocations succeeded
            ClearUnitIfMemFailed((x->fdndels[i] && x->fdndamps[i]), x);
        }
         
        // diffuser section
        float diffscale = (float)(x->fdnlens[3] / (210. + 159. + 562. + 410.));
//        float spread1 = spread;
//        float spread2 = 3.0 * spread;
        // Clamp spread and adjust coefficients for stability
        float spread1 = fmin(fmax(spread, 0.0f), 11.0f);
        float spread2 = fmin(3.0f * spread1, 33.0f);  // Clamp the 3x multiplier too
        
        // Adjust diffusion coefficients inversely with spread
        float baseCoef1 = 0.75f;
        float baseCoef2 = 0.625f;
        float coef1 = baseCoef1 * (11.0f - spread1) / 11.0f;  // Reduce coefficient as spread increases
        float coef2 = baseCoef2 * (11.0f - spread1) / 11.0f;

        int b = 210;
        float r = 0.125541;
        int a = (int)(spread1 * r);
        int c = 210 + 159 + a;
        int cc = c - b;
        r = 0.854046;
        a = (int)(spread2 * r);
        int d = 210 + 159 + 562 + a;
        int dd = d - c;
        int e = 1341 - d;

        x->ldifs[0] = make_diffuser(x, f_round(diffscale * b),coef1);// 0.75);
        x->ldifs[1] = make_diffuser(x, f_round(diffscale * cc),coef1);// 0.75);
        x->ldifs[2] = make_diffuser(x, f_round(diffscale * dd),coef2);// 0.625);
        x->ldifs[3] = make_diffuser(x, f_round(diffscale * e),coef2);// 0.625);
        ClearUnitIfMemFailed(x->ldifs[0] && x->ldifs[1] && x->ldifs[2] && x->ldifs[3], x);

        b = 210;
        r = -0.568366;
        a = (int)(spread1 * r);
        c = 210 + 159 + a;
        cc = c - b;
        r = -0.126815;
        a = (int)(spread2 * r);
        d = 210 + 159 + 562 + a;
        dd = d - c;
        e = 1341 - d;

        x->rdifs[0] = make_diffuser(x, f_round(diffscale * b),coef1);// 0.75);
        x->rdifs[1] = make_diffuser(x, f_round(diffscale * cc),coef1);// 0.75);
        x->rdifs[2] = make_diffuser(x, f_round(diffscale * dd),coef2);// 0.625);
        x->rdifs[3] = make_diffuser(x, f_round(diffscale * e),coef2);// 0.625);
        ClearUnitIfMemFailed(x->rdifs[0] && x->rdifs[1] && x->rdifs[2] && x->rdifs[3], x);

        x->taps[0] = 5 + (int)(0.410 * largestdelay);
        x->taps[1] = 5 + (int)(0.300 * largestdelay);
        x->taps[2] = 5 + (int)(0.155 * largestdelay);
        x->taps[3] = 5; //+ f_round(0.000 * largestdelay);

        for (int i = 0; i < FDNORDER; i++) {
            x->tapgains[i] = (float)(pow(alpha, (double)(x->taps[i])));
        }

        x->tapdelay = make_fixeddelay(x, 44000, 44000);
        ClearUnitIfMemFailed(x->tapdelay, x);

        // init the slope values
        x->earlylevelslope = x->drylevelslope = x->taillevelslope = 0.f;
        clear_outputs(x, 1);
//        ClearUnitOutputs(x, 1);
	}
	return (x);
}


// NOT CALLED!, we use dsp_free for a generic free function
void gverb_free(t_gverb *x)
{
    dsp_free((t_pxobject*)x);
    // Free all allocated memory
    if (x->inputdamper) sysmem_freeptr(x->inputdamper);
    if (x->tapdelay) {
        if (x->tapdelay->buf) sysmem_freeptr(x->tapdelay->buf);
        sysmem_freeptr(x->tapdelay);
    }
    
    for (int i = 0; i < FDNORDER; i++) {
        if (x->fdndels[i]) {
            if (x->fdndels[i]->buf) sysmem_freeptr(x->fdndels[i]->buf);
            sysmem_freeptr(x->fdndels[i]);
        }
        if (x->fdndamps[i]) sysmem_freeptr(x->fdndamps[i]);
        if (x->ldifs[i]) {
            if (x->ldifs[i]->buf) sysmem_freeptr(x->ldifs[i]->buf);
            sysmem_freeptr(x->ldifs[i]);
        }
        if (x->rdifs[i]) {
            if (x->rdifs[i]->buf) sysmem_freeptr(x->rdifs[i]->buf);
            sysmem_freeptr(x->rdifs[i]);
        }
    }
//    if (x->m_dlybuf) sysmem_freeptr(x->m_dlybuf);
}


void gverb_assist(t_gverb *x, void *b, long m, long a, char *s)
{
	if (m == ASSIST_INLET) { //inlet
        switch (a) {
            case 0: sprintf(s, "(signal) Input"); break;
            case 1: sprintf(s, "(signal) Room Size"); break;
            case 2: sprintf(s, "(signal) Reverb Time"); break;
            case 3: sprintf(s, "(signal) Damping"); break;
            case 4: sprintf(s, "(signal) Input Bandwidth"); break;
            case 5: sprintf(s, "(signal) Dry Level"); break;
            case 6: sprintf(s, "(signal) Early Level"); break;
            case 7: sprintf(s, "(signal) Tail Level"); break;
        }
	}
	else {	// outlet
		sprintf(s, "(signal) Outlet %ld", a);
	}
}



// registers a function for the signal chain in Max
void gverb_dsp64(t_gverb *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags)
{
//	post("my sample rate is: %f", samplerate);
//    x->m_sr = samplerate;
	// instead of calling dsp_add(), we send the "dsp_add64" message to the object representing the dsp chain
	// the arguments passed are:
	// 1: the dsp64 object passed-in by the calling function
	// 2: the symbol of the "dsp_add64" message we are sending
	// 3: a pointer to your object
	// 4: a pointer to your 64-bit perform method
	// 5: flags to alter how the signal chain handles your object -- just pass 0
	// 6: a generic pointer that you can use to pass any additional data to your perform method

	object_method(dsp64, gensym("dsp_add64"), x, gverb_perform64, 0, NULL);
}


// this is the 64-bit perform method audio vectors
void gverb_perform64(t_gverb *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam)
{
	t_double *in = ins[0];		// we get audio for each inlet of the object from the **ins argument
//    int n = sampleframes;
    float roomsize = sc_min(*ins[1], x->maxroomsize);
    float revtime = *ins[2];
    float damping = *ins[3];
    float inputbandwidth = *ins[4];
    // float spread = IN0(5); // spread can only be set at inittime
    float drylevel = *ins[5];
    float earlylevel = *ins[6];
    float taillevel = *ins[7];

    t_double *outl = outs[0];    // we get audio for each outlet of the object from the **outs argument
    t_double *outr = outs[1];
    
    float earlylevelslope, taillevelslope, drylevelslope;
    float* fdngainslopes;
    float* tapgainslopes;
    g_diffuser** ldifs = x->ldifs;
    g_diffuser** rdifs = x->rdifs;
    float* u = x->u;
    float* f = x->f;
    float* d = x->d;
    g_damper* inputdamper = x->inputdamper;
    float* tapgains = x->tapgains;
    g_fixeddelay* tapdelay = x->tapdelay;
    int* taps = x->taps;
    g_damper** fdndamps = x->fdndamps;
    g_fixeddelay** fdndels = x->fdndels;
    float* fdngains = x->fdngains;
    int* fdnlens = x->fdnlens;
    
//    // At the start of gverb_perform64, add:
//    float energy = 0.0f;
//    const float MAX_ENERGY = 100.0f;  // Adjust this threshold as needed
//    
    if ((roomsize != x->roomsize) || (revtime != x->revtime) || (damping != x->damping)
        || (inputbandwidth != x->inputbandwidth) || (drylevel != x->drylevel) || (earlylevel != x->earlylevel)
        || (taillevel != x->taillevel)) {
        // these should calc slopes for k-rate interpolation
        gverb_set_roomsize(x, roomsize);
        gverb_set_revtime(x, revtime);
        gverb_set_damping(x, damping);
        gverb_set_inputbandwidth(x, inputbandwidth);
        drylevel = gverb_set_drylevel(x, drylevel);
        earlylevel = gverb_set_earlylevel(x, earlylevel);
        taillevel = gverb_set_taillevel(x, taillevel);
    }

    earlylevelslope = x->earlylevelslope;
    taillevelslope = x->taillevelslope;
    drylevelslope = x->drylevelslope;
    fdngainslopes = x->fdngainslopes;
    tapgainslopes = x->tapgainslopes;


    for (int i = 0; i < sampleframes; i++) {
        float sign, sum, lsum, rsum, x;
        if (isnan(in[i]))
            x = 0.f;
        else
            x = in[i];
        sum = 0.f;
        sign = 1.f;

        float z = damper_do(inputdamper, x);
        z = diffuser_do(ldifs[0], z);

        for (int j = 0; j < FDNORDER; j++) {
            u[j] = tapgains[j] * fixeddelay_read(tapdelay, taps[j]);
        }

        fixeddelay_write(tapdelay, z);

        for (int j = 0; j < FDNORDER; j++) {
            d[j] = damper_do(fdndamps[j], fdngains[j] * fixeddelay_read(fdndels[j], fdnlens[j]));
        }

        for (int j = 0; j < FDNORDER; j++) {
            sum += sign * (taillevel * d[j] + earlylevel * u[j]);
            sign = -sign;
        }

        sum += x * earlylevel;
        lsum = sum;
        rsum = sum;

        gverb_fdnmatrix(d, f);
        
        for (int j = 0; j < FDNORDER; j++) {
            fixeddelay_write(fdndels[j], u[j] + f[j]);
        }

        lsum = diffuser_do(ldifs[1], lsum);
        lsum = diffuser_do(ldifs[2], lsum);
        lsum = diffuser_do(ldifs[3], lsum);
        rsum = diffuser_do(rdifs[1], rsum);
        rsum = diffuser_do(rdifs[2], rsum);
        rsum = diffuser_do(rdifs[3], rsum);

        x = x * drylevel;
        outl[i] = lsum + x;
        outr[i] = rsum + x;

        drylevel += drylevelslope;
        taillevel += taillevelslope;
        earlylevel += earlylevelslope;
        for (int j = 0; j < FDNORDER; j++) {
            fdngains[j] += fdngainslopes[j];
            tapgains[j] += tapgainslopes[j];
        }
//        
//        // Inside the sample loop, after FDN processing:
//        float energy = 0.0f;
//        const float MAX_ENERGY = 100.0f;  // Adjust this threshold as needed
//
//        for (int j = 0; j < FDNORDER; j++) {
//            energy = energy * 0.999f + d[j] * d[j] * 0.001f;  // Smooth energy measurement
//        }
//
//        if (energy > MAX_ENERGY) {
//            // Reset the delays and dampers
//            for (int j = 0; j < FDNORDER; j++) {
//                clear_buffer(fdndels[j]->size, fdndels[j]->buf);
//                // Reset damper state - assuming you have a function for this
//                damper_clear(fdndamps[j]);
//            }
//            for (int j = 0; j < 4; j++) {
//                clear_buffer(ldifs[j]->size, ldifs[j]->buf);
//                clear_buffer(rdifs[j]->size, rdifs[j]->buf);
//            }
//            energy = 0.0f;
//        }
    }

    // store vals back to the struct
    for (int i = 0; i < FDNORDER; i++) {
        x->ldifs[i] = ldifs[i];
        x->rdifs[i] = rdifs[i];
        x->u[i] = u[i];
        x->f[i] = f[i];
        x->d[i] = d[i];
        x->tapgains[i] = tapgains[i];
        x->taps[i] = taps[i];
        x->fdndamps[i] = fdndamps[i];
        x->fdndels[i] = fdndels[i];
        x->fdngains[i] = fdngains[i];
        x->fdnlens[i] = fdnlens[i];
        x->fdngainslopes[i] = 0.f;
        x->tapgainslopes[i] = 0.f;
    }
    x->inputdamper = inputdamper;
    x->tapdelay = tapdelay;
    // clear the slopes
    x->earlylevelslope = x->taillevelslope = x->drylevelslope = 0.f;

//             this perform method simply copies the input to the output, offsetting the value
//    while (n--)
//        *outl++ = *in++;
}
