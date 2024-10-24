/**
    @file
    texprand~: a simple audio object for Max
    original by: jeremy bernstein, jeremy@bootsquad.com
    @ingroup examples
*/

#include "ext.h"            // standard Max include, always required (except in Jitter)
#include "ext_obex.h"        // required for "new" style objects
#include "z_dsp.h"            // required for MSP objects
#include "time.h"


// struct to represent the object's state
typedef struct _texprand {
    t_pxobject        ob;            // the object itself (t_pxobject in MSP instead of t_object)
    double lo;
    double hi;
    double last_trig;
    double last_output;
} t_texprand;


// method prototypes
void *texprand_new(t_symbol *s, long argc, t_atom *argv);
void texprand_free(t_texprand *x);
void texprand_assist(t_texprand *x, void *b, long m, long a, char *s);
void texprand_float(t_texprand *x, double f);
void texprand_dsp64(t_texprand *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags);
void texprand_perform64(t_texprand *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam);


// global class pointer variable
static t_class *texprand_class = NULL;


//***********************************************************************************************

void ext_main(void *r)
{
    // object initialization, note the use of dsp_free for the freemethod, which is required
    // unless you need to free allocated memory, in which case you should call dsp_free from
    // your custom free function.

    t_class *c = class_new("texprand~", (method)texprand_new, (method)dsp_free, (long)sizeof(t_texprand), 0L, A_GIMME, 0);

    class_addmethod(c, (method)texprand_float,        "float",    A_FLOAT, 0);
    class_addmethod(c, (method)texprand_dsp64,        "dsp64",    A_CANT, 0);
    class_addmethod(c, (method)texprand_assist,    "assist",    A_CANT, 0);

    class_dspinit(c);
    class_register(CLASS_BOX, c);
    texprand_class = c;
}


void *texprand_new(t_symbol *s, long argc, t_atom *argv)
{
    t_texprand *x = (t_texprand *)object_alloc(texprand_class);

    if (x) {
        dsp_setup((t_pxobject *)x, 3);    // MSP inlets: arg is # of inlets and is REQUIRED!
        // use 0 if you don't need inlets
        outlet_new(x, "signal");         // signal outlet (note "signal" rather than NULL)
        
        x->lo = 1.0;  // Default low value
        x->hi = 10.0;  // Default high value
        x->last_trig = 0.0;
        x->last_output = x->lo;

        // Parse arguments if provided
        if (argc >= 2) {
            x->lo = atom_getfloat(argv);
            x->hi = atom_getfloat(argv + 1);
        }

        // Ensure lo and hi have the same sign and are non-zero
        if (x->lo * x->hi <= 0 || x->lo == 0 || x->hi == 0) {
            object_error((t_object *)x, "Invalid lo and hi values. They must have the same sign and be non-zero.");
            x->lo = 1.0;
            x->hi = 10.0;
        }

        srand(time(NULL));  // Initialize random seed

    }
    return (x);
}


// NOT CALLED!, we use dsp_free for a generic free function
void texprand_free(t_texprand *x)
{
    dsp_free((t_pxobject *)x);
}


void texprand_assist(t_texprand *x, void *b, long m, long a, char *s)
{
    if (m == ASSIST_INLET) { //inlet
        switch (a) {
            case 0: sprintf(s, "Trigger signal (0 to 1)"); break;
            case 1: sprintf(s, "Low value (signal/float)"); break;
            case 2: sprintf(s, "High value (signal/float)"); break;
        }
    }   else {    // outlet
        sprintf(s, "Random value (signal)");
    }
}


void texprand_float(t_texprand *x, double f)
{
    x->lo = f;
    if (x->lo * x->hi <= 0 || x->lo == 0) {
        object_error((t_object *)x, "Invalid lo value. It must have the same sign as hi and be non-zero.");
        x->lo = (x->hi > 0) ? 1.0 : -10.0;
    }
}


// registers a function for the signal chain in Max
void texprand_dsp64(t_texprand *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags)
{
    //post("my sample rate is: %f", samplerate);

    // instead of calling dsp_add(), we send the "dsp_add64" message to the object representing the dsp chain
    // the arguments passed are:
    // 1: the dsp64 object passed-in by the calling function
    // 2: the symbol of the "dsp_add64" message we are sending
    // 3: a pointer to your object
    // 4: a pointer to your 64-bit perform method
    // 5: flags to alter how the signal chain handles your object -- just pass 0
    // 6: a generic pointer that you can use to pass any additional data to your perform method

    object_method(dsp64, gensym("dsp_add64"), x, texprand_perform64, 0, NULL);
}


// this is the 64-bit perform method audio vectors
void texprand_perform64(t_texprand *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam)
{
    t_double *trig = ins[0];        // we get audio for each inlet of the object from the **ins argument
    t_double *lo_in = ins[1];
    t_double *hi_in = ins[2];
    t_double *out = outs[0];    // we get audio for each outlet of the object from the **outs argument
    t_double last_trig = x->last_trig;
    t_double last_output = x->last_output;
    //t_double hi = x->hi;
    
    int n = sampleframes;

    // this perform method simply copies the input to the output, offsetting the value
    while (n--) {
        t_double current_trig = *trig++;
        t_double current_lo = *lo_in++;
        t_double currenth_hi = *hi_in++;
        
        if (current_trig > 0 && last_trig <= 0) {
            // Check if lo has changed and validate it
            if (current_lo != x->lo) {
                if (current_lo * currenth_hi > 0 && current_lo != 0) {
                    x->lo = current_lo;
                } else {
                    // If invalid, use the last valid lo value
                    current_lo = x->lo;
                }
            }
            if (currenth_hi != x->hi) {
                if (current_lo * currenth_hi > 0 && currenth_hi != 0) {
                    x->hi = currenth_hi;
                } else {
                    // If invalid, use the last valid lo value
                    currenth_hi = x->hi;
                }
            }

            // Generate new random value
            t_double r = (double)rand() / RAND_MAX;
            t_double lambda = log(currenth_hi / current_lo);
            last_output = current_lo * exp(r * lambda);
        }

        
        *out++ = last_output;
        last_trig = current_trig;
    }
    x->last_trig = last_trig;
    x->last_output = last_output;
}

