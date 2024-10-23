#include "ext.h"
#include "ext_obex.h"
#include "z_dsp.h"
#include <math.h>
#include <stdlib.h>
#include <time.h>

typedef struct _exprand {
    t_pxobject obj;
    double lo;
    double hi;
    double last_trig;
    double last_output;
} t_exprand;

void *exprand_new(t_symbol *s, long argc, t_atom *argv);
void exprand_free(t_exprand *x);
void exprand_assist(t_exprand *x, void *b, long m, long a, char *s);
void exprand_float(t_exprand *x, double f);
void exprand_dsp64(t_exprand *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags);
void exprand_perform64(t_exprand *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam);

static t_class *exprand_class;

void ext_main(void *r)
{
    t_class *c = class_new("exprand~", (method)exprand_new, (method)exprand_free, sizeof(t_exprand), 0, A_GIMME, 0);

    class_addmethod(c, (method)exprand_float, "float", A_FLOAT, 0);
    class_addmethod(c, (method)exprand_dsp64, "dsp64", A_CANT, 0);
    class_addmethod(c, (method)exprand_assist, "assist", A_CANT, 0);

    class_dspinit(c);
    class_register(CLASS_BOX, c);
    exprand_class = c;
}

void *exprand_new(t_symbol *s, long argc, t_atom *argv)
{
    t_exprand *x = (t_exprand *)object_alloc(exprand_class);
    if (x) {
        dsp_setup((t_pxobject *)x, 2);  // Two inlets: trigger and lo
        outlet_new((t_object *)x, "signal");  // One outlet for the random value

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
    return x;
}

void exprand_free(t_exprand *x)
{
    dsp_free((t_pxobject *)x);
}

void exprand_assist(t_exprand *x, void *b, long m, long a, char *s)
{
    if (m == ASSIST_INLET) {
        switch (a) {
            case 0: sprintf(s, "Trigger signal (0 to 1)"); break;
            case 1: sprintf(s, "Low value (signal/float)"); break;
        }
    } else {
        sprintf(s, "Random value (signal)");
    }
}

void exprand_float(t_exprand *x, double f)
{
    x->lo = f;
    if (x->lo * x->hi <= 0 || x->lo == 0) {
        object_error((t_object *)x, "Invalid lo value. It must have the same sign as hi and be non-zero.");
        x->lo = (x->hi > 0) ? 1.0 : -10.0;
    }
}

void exprand_dsp64(t_exprand *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags)
{
    object_method(dsp64, gensym("dsp_add64"), x, exprand_perform64, 0, NULL);
}

void exprand_perform64(t_exprand *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam)
{
    double *trig = ins[0];
    double *lo_in = ins[1];
    double *out = outs[0];
    double last_trig = x->last_trig;
    double last_output = x->last_output;
    double hi = x->hi;

    while (sampleframes--) {
        double current_trig = *trig++;
        double current_lo = *lo_in++;

        if (current_trig > 0 && last_trig <= 0) {
            // Check if lo has changed and validate it
            if (current_lo != x->lo) {
                if (current_lo * hi > 0 && current_lo != 0) {
                    x->lo = current_lo;
                } else {
                    // If invalid, use the last valid lo value
                    current_lo = x->lo;
                }
            }

            // Generate new random value
            double r = (double)rand() / RAND_MAX;
            double lambda = log(hi / current_lo);
            last_output = current_lo * exp(r * lambda);
        }

        *out++ = last_output;
        last_trig = current_trig;
    }

    x->last_trig = last_trig;
    x->last_output = last_output;
}
