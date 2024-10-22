/**
    @file
    blowshelf~: lowshelf
    The B equalization suite is based on the Second Order Section (SOS) biquad UGen.
    @ingroup examples
*/

#include "ext.h"            // standard Max include, always required (except in Jitter)
#include "ext_obex.h"        // required for "new" style objects
#include "z_dsp.h"            // required for MSP objects


// struct to represent the object's state
typedef struct _blowshelf {
    t_pxobject        ob;            // the object itself (t_pxobject in MSP instead of t_object)
    double samplerate;
    double m_y1, m_y2, m_a0, m_a1, m_a2, m_b1, m_b2;
    float m_freq, m_rs, m_db;
    //double            offset;     // the value of a property of our object
} t_blowshelf;


// method prototypes
void *blowshelf_new(t_symbol *s, long argc, t_atom *argv);
void blowshelf_free(t_blowshelf *x);
void blowshelf_assist(t_blowshelf *x, void *b, long m, long a, char *s);
//void blowshelf_float(t_blowshelf *x, double f);
void blowshelf_dsp64(t_blowshelf *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags);
void blowshelf_perform64(t_blowshelf *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam);


// global class pointer variable
static t_class *blowshelf_class = NULL;


//***********************************************************************************************

void ext_main(void *r)
{
    // object initialization, note the use of dsp_free for the freemethod, which is required
    // unless you need to free allocated memory, in which case you should call dsp_free from
    // your custom free function.

    t_class *c = class_new("blowshelf~", (method)blowshelf_new, (method)dsp_free, (long)sizeof(t_blowshelf), 0L, A_GIMME, 0);

//    class_addmethod(c, (method)blowshelf_float,        "float",    A_FLOAT, 0);
    class_addmethod(c, (method)blowshelf_dsp64,        "dsp64",    A_CANT, 0);
    class_addmethod(c, (method)blowshelf_assist,    "assist",    A_CANT, 0);

    class_dspinit(c);
    class_register(CLASS_BOX, c);
    blowshelf_class = c;
}

//BLowShelf.ar(in, freq: 1200.0, rs: 1.0, db: 0.0, mul: 1.0, add: 0.0)
//Arguments:
//in
//input signal to be processed.
//
//freq
//center frequency. WARNING: due to the nature of its implementation frequency values close to 0 may cause glitches and/or extremely loud audio artifacts!
//
//rs
//the reciprocal of S. Shell boost/cut slope. When S = 1, the shelf slope is as steep as it can be and remain monotonically increasing or decreasing gain with frequency. The shelf slope, in dB/octave, remains proportional to S for all other values for a fixed freq/SampleRate.ir and db.
//
//db
//gain. boost/cut the center frequency in dBs.
//
//mul
//add

void *blowshelf_new(t_symbol *s, long argc, t_atom *argv)
{
    t_blowshelf *x = (t_blowshelf *)object_alloc(blowshelf_class);

    if (x) {
        dsp_setup((t_pxobject *)x, 4);    // MSP inlets: arg is # of inlets and is REQUIRED!
        // use 0 if you don't need inlets
        outlet_new(x, "signal");         // signal outlet (note "signal" rather than NULL)
        x->m_freq = 440.0;
        x->m_rs = 1.0;
        x->m_db = 0.0;
        x->m_y1 = 0.0;
        x->m_y2 = 0.0;
        //atom_arg_getfloat(&x->m_freq, 0, argc, argv);
        //freq = x->m_freq;
    }
    return (x);
}


// NOT CALLED!, we use dsp_free for a generic free function
void blowshelf_free(t_blowshelf *x)
{
    ;
}


void blowshelf_assist(t_blowshelf *x, void *b, long m, long a, char *s)
{
    if (m == ASSIST_INLET) { //inlet
        switch (a) {
            case 0: sprintf(s, "(signal) Input"); break;
            case 1: sprintf(s, "(signal) Frequency"); break;
            case 2: sprintf(s, "(signal) RS"); break;
            case 3: sprintf(s, "(signal) dB"); break;
        }
//        sprintf(s, "I am inlet %ld", a);
    }
    else {    // outlet
        sprintf(s, "I am outlet %ld", a);
    }
}



Float32 zapgremlins(Float32 x) {
    Float32 absx = fabsf(x);
    // very small numbers fail the first test, eliminating denormalized numbers
    //    (zero also fails the first test, but that is OK since it returns zero.)
    // very large numbers fail the second test, eliminating infinities
    // Not-a-Numbers fail both tests and are eliminated.
    return (absx > (Float32)1e-15 && absx < (Float32)1e15) ? x : (Float32)0.;
}
//void blowshelf_float(t_blowshelf *x, double f)
//{
//    x->offset = f;
//}


// registers a function for the signal chain in Max
void blowshelf_dsp64(t_blowshelf *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags)
{
    //post("my sample rate is: %f", samplerate);
//    x->samplerate = samplerate;
//    lowshelf_updateCoefficients(x, x->m_freq, x->m_rs, x->m_db);
    // instead of calling dsp_add(), we send the "dsp_add64" message to the object representing the dsp chain
    // the arguments passed are:
    // 1: the dsp64 object passed-in by the calling function
    // 2: the symbol of the "dsp_add64" message we are sending
    // 3: a pointer to your object
    // 4: a pointer to your 64-bit perform method
    // 5: flags to alter how the signal chain handles your object -- just pass 0
    // 6: a generic pointer that you can use to pass any additional data to your perform method

    object_method(dsp64, gensym("dsp_add64"), x, blowshelf_perform64, 0, NULL);
}


// this is the 64-bit perform method audio vectors
void blowshelf_perform64(t_blowshelf *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam)
{
    
//    t_double *inL = ins[0];        // we get audio for each inlet of the object from the **ins argument
//    t_double *outL = outs[0];    // we get audio for each outlet of the object from the **outs argument
    t_double *in = ins[0];
    t_double *freq = ins[1];
    t_double *rs = ins[2];
    t_double *db = ins[3];
    t_double *out = outs[0];
    int n = sampleframes;
    
    t_double y1 = x->m_y1;
    t_double y2 = x->m_y2;
    t_double a0 = x->m_a0;
    t_double a1 = x->m_a1;
    t_double a2 = x->m_a2;
    t_double b1 = x->m_b1;
    t_double b2 = x->m_b2;
    // this perform method simply copies the input to the output, offsetting the value
    while (n--) {
        t_double nextfreq = *freq++;
        t_double nextrs = *rs++;
        t_double nextdb = *db++;
        if (x->m_freq != nextfreq || x->m_rs != nextrs || x->m_db != nextdb) {
            //            lowshelf_updateCoefficients(x, nextfreq, nextrs, nextdb);
            double a = pow(10., nextdb * 0.025);
            double w0 = TWOPI * nextfreq / sys_getsr();
            double cosw0 = cos(w0);
            double sinw0 = sin(w0);
            double alpha = sinw0 * 0.5 * sqrt((a + (1. / a)) * (nextrs - 1.) + 2.);
            double i = (a + 1.) * cosw0;
            double j = (a - 1.) * cosw0;
            double k = 2. * sqrt(a) * alpha;
            double b0rz = 1. / ((a + 1.) + j + k);
            
            a0 = a * ((a + 1.) - j + k) * b0rz;
            a1 = 2. * a * ((a - 1.) - i) * b0rz;
            a2 = a * ((a + 1.) - j - k) * b0rz;
            b1 = 2. * ((a - 1.) + i) * b0rz;
            b2 = ((a + 1.) + j - k) * -b0rz;
            
            x->m_freq = nextfreq;
            x->m_rs = nextrs;
            x->m_db = nextdb;
            x->m_a0 = a0;
            x->m_a1 = a1;
            x->m_a2 = a2;
            x->m_b1 = b1;
            x->m_b2 = b2;
//            *out++ = *in++ + x->m_freq;
        } else {
            a0 = x->m_a0;
            a1 = x->m_a1;
            a2 = x->m_a2;
            b1 = x->m_b1;
            b2 = x->m_b2;
        }
        // Apply the filter
        t_double y0 = *in++ + b1 * y1 + b2 * y2;
        *out++ = a0 * y0 + a1 * y1 + a2 * y2;

        y2 = y1;
        y1 = y0;
    }
    // Store state
    x->m_y1 = zapgremlins(y1);
    x->m_y2 = zapgremlins(y2);
}
