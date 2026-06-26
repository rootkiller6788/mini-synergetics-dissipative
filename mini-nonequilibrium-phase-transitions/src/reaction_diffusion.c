#include <stdlib.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <stdio.h>
#include <string.h>

/* Reaction-Diffusion Systems */
/* Oregonator (BZ), Brusselator, Schlogl models */

typedef struct { double k1,k2,k3,k4,k5; double A,B,f; } OregonatorParams;
typedef struct { double X,Y,Z,t; } OregonatorState;
typedef struct { double A,B; } BrusselatorParams;
typedef struct { double X,Y,t; } BrusselatorState;
typedef struct { double k1,k2,k3,k4; } SchloglParams;

static void oregonator_rhs(const OregonatorState *s, const OregonatorParams *p,
                            double *dX, double *dY, double *dZ)
{
    double X=s->X, Y=s->Y, Z=s->Z;
    *dX = p->k1*p->A*Y - p->k2*X*Y + p->k3*p->A*X - 2.0*p->k4*X*X;
    *dY = -p->k1*p->A*Y - p->k2*X*Y + p->f*p->k5*p->B*Z;
    *dZ = 2.0*p->k3*p->A*X - p->k5*p->B*Z;
}

double oregonator_rk4_step(OregonatorState *s, const OregonatorParams *p, double dt)
{
    double X0=s->X, Y0=s->Y, Z0=s->Z;
    double dX1,dY1,dZ1;
    oregonator_rhs(s,p,&dX1,&dY1,&dZ1);
    OregonatorState s2={X0+0.5*dt*dX1,Y0+0.5*dt*dY1,Z0+0.5*dt*dZ1,0};
    double dX2,dY2,dZ2;
    oregonator_rhs(&s2,p,&dX2,&dY2,&dZ2);
    OregonatorState s3={X0+0.5*dt*dX2,Y0+0.5*dt*dY2,Z0+0.5*dt*dZ2,0};
    double dX3,dY3,dZ3;
    oregonator_rhs(&s3,p,&dX3,&dY3,&dZ3);
    OregonatorState s4={X0+dt*dX3,Y0+dt*dY3,Z0+dt*dZ3,0};
    double dX4,dY4,dZ4;
    oregonator_rhs(&s4,p,&dX4,&dY4,&dZ4);
    s->X = X0 + (dt/6.0)*(dX1+2.0*dX2+2.0*dX3+dX4);
    s->Y = Y0 + (dt/6.0)*(dY1+2.0*dY2+2.0*dY3+dY4);
    s->Z = Z0 + (dt/6.0)*(dZ1+2.0*dZ2+2.0*dZ3+dZ4);
    s->t += dt;
    if (s->X<0.0) s->X=0.0;
    if (s->Y<0.0) s->Y=0.0;
    if (s->Z<0.0) s->Z=0.0;
    return sqrt(s->X*s->X + s->Y*s->Y + s->Z*s->Z);
}

void oregonator_set_defaults(OregonatorParams *p)
{
    if (!p) return;
    p->k1=1.28; p->k2=2.4e6; p->k3=33.6; p->k4=3.0e3; p->k5=1.0;
    p->A=0.06; p->B=0.02; p->f=0.9;
}

void oregonator_init_state(OregonatorState *s)
{
    if (!s) return;
    s->X=1.0e-7; s->Y=1.0e-7; s->Z=1.0e-7; s->t=0.0;
}

/* Brusselator: dx/dt = A - (B+1)*X + X^2*Y, dy/dt = B*X - X^2*Y */
static void brusselator_rhs(const BrusselatorState *s, const BrusselatorParams *p,
                             double *dX, double *dY)
{
    double X2Y = s->X*s->X*s->Y;
    *dX = p->A - (p->B+1.0)*s->X + X2Y;
    *dY = p->B*s->X - X2Y;
}

double brusselator_rk4_step(BrusselatorState *s, const BrusselatorParams *p, double dt)
{
    double X0=s->X, Y0=s->Y;
    double dX1,dY1;
    brusselator_rhs(s,p,&dX1,&dY1);
    BrusselatorState s2={X0+0.5*dt*dX1,Y0+0.5*dt*dY1,0};
    double dX2,dY2;
    brusselator_rhs(&s2,p,&dX2,&dY2);
    BrusselatorState s3={X0+0.5*dt*dX2,Y0+0.5*dt*dY2,0};
    double dX3,dY3;
    brusselator_rhs(&s3,p,&dX3,&dY3);
    BrusselatorState s4={X0+dt*dX3,Y0+dt*dY3,0};
    double dX4,dY4;
    brusselator_rhs(&s4,p,&dX4,&dY4);
    s->X = X0 + (dt/6.0)*(dX1+2.0*dX2+2.0*dX3+dX4);
    s->Y = Y0 + (dt/6.0)*(dY1+2.0*dY2+2.0*dY3+dY4);
    s->t += dt;
    if (s->X<0.0) s->X=0.0;
    if (s->Y<0.0) s->Y=0.0;
    return sqrt(s->X*s->X + s->Y*s->Y);
}

void brusselator_fixed_point(const BrusselatorParams *p, double *X0, double *Y0)
{
    if (X0) *X0 = p->A;
    if (Y0) *Y0 = p->B / p->A;
}

int brusselator_is_oscillatory(const BrusselatorParams *p)
{
    double Bc = 1.0 + p->A*p->A;
    return (p->B > Bc) ? 1 : 0;
}

double brusselator_max_eigenvalue(const BrusselatorParams *p)
{
    double A2 = p->A*p->A;
    double trace = p->B - 1.0 - A2;
    double disc = trace*trace - 4.0*A2;
    if (disc >= 0.0) return 0.5*(trace + sqrt(disc));
    return 0.5*trace;
}

void brusselator_init_state(BrusselatorState *s, double A, double B)
{
    if (!s) return;
    s->X = A + 0.01;
    s->Y = B/A + 0.01;
    s->t = 0.0;
}

void brusselator_set_params(BrusselatorParams *p, double A, double B)
{
    if (!p) return;
    p->A = A; p->B = B;
}


/* Schlogl model: dx/dt = -k1*x^3 + k2*x^2 - k3*x + k4
 * First-order non-equilibrium phase transition.
 * Bistability when 3 real roots exist with 2 stable.
 * V(x) = (k1/4)*x^4 - (k2/3)*x^3 + (k3/2)*x^2 - k4*x */

static double schlogl_potential(double x, const SchloglParams *p)
{
    return 0.25*p->k1*x*x*x*x - (p->k2/3.0)*x*x*x + 0.5*p->k3*x*x - p->k4*x;
}

int schlogl_steady_states(const SchloglParams *p,
                           double *x1, double *x2, double *x3)
{
    double a=p->k1, b=p->k2, c=p->k3, d=p->k4;
    if (fabs(a)<1e-15) return 0;
    double bm=-b/a, cm=c/a, dm=-d/a;
    double pc=cm - bm*bm/3.0;
    double qc=2.0*bm*bm*bm/27.0 - bm*cm/3.0 + dm;
    double disc=qc*qc/4.0 + pc*pc*pc/27.0;
    if (disc > 1e-15) {
        double u=cbrt(-qc/2.0+sqrt(disc));
        double v=cbrt(-qc/2.0-sqrt(disc));
        if (x1) *x1=u+v-bm/3.0;
        return 1;
    } else if (disc < -1e-15) {
        double phi=acos(3.0*qc/(2.0*pc)*sqrt(-3.0/pc));
        double r=2.0*sqrt(-pc/3.0);
        if (x1) *x1=r*cos(phi/3.0)-bm/3.0;
        if (x2) *x2=r*cos((phi+2.0*M_PI)/3.0)-bm/3.0;
        if (x3) *x3=r*cos((phi+4.0*M_PI)/3.0)-bm/3.0;
        return 3;
    }
    double t1=2.0*cbrt(-qc/2.0);
    double t2=-cbrt(-qc/2.0);
    if (x1) *x1=t1-bm/3.0;
    if (x2) *x2=t2-bm/3.0;
    return 2;
}

int schlogl_classify_stability(const SchloglParams *p,
                                const double *xss, int n, int *stab)
{
    if (!p || !xss || !stab) return -1;
    for (int i=0; i<n; i++) {
        double d2V = 3.0*p->k1*xss[i]*xss[i] - 2.0*p->k2*xss[i] + p->k3;
        stab[i] = (d2V > 1e-12) ? 1 : 0;
    }
    return n;
}

int schlogl_is_bistable(const SchloglParams *p)
{
    double r1,r2,r3;
    int n=schlogl_steady_states(p,&r1,&r2,&r3);
    if (n<3) return 0;
    double rs[3]={r1,r2,r3}; int st[3];
    schlogl_classify_stability(p,rs,n,st);
    int sc=0;
    for (int i=0;i<n;i++) if (st[i]) sc++;
    return (sc>=2)?1:0;
}

void schlogl_set_bistable(SchloglParams *p)
{
    if (!p) return;
    p->k1=1.0; p->k2=3.0; p->k3=3.0; p->k4=0.5;
}

double schlogl_barrier(const SchloglParams *p)
{
    double r1,r2,r3;
    int n=schlogl_steady_states(p,&r1,&r2,&r3);
    if (n<3) return 0.0;
    double xs[3]={r1,r2,r3};
    for (int i=0;i<2;i++)
        for (int j=i+1;j<3;j++)
            if (xs[i]>xs[j]) {double t=xs[i];xs[i]=xs[j];xs[j]=t;}
    int st[3];
    schlogl_classify_stability(p,xs,3,st);
    double V[3];
    for (int i=0;i<3;i++) V[i]=schlogl_potential(xs[i],p);
    if (st[0]&&st[2]&&!st[1]) {
        double b1=V[1]-V[0], b2=V[1]-V[2];
        return (b1>b2)?b1:b2;
    }
    double mx=V[0], mn=V[0];
    for (int i=1;i<3;i++) {if(V[i]>mx)mx=V[i]; if(V[i]<mn)mn=V[i];}
    return mx-mn;
}
