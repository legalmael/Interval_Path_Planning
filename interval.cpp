// Compatiblity layer with ibex for the simple interval library from Luc JAULIN, with minor modifications from Fabrice LE BARS and Jeremy NICOLA.

#include "interval.h"

#if defined(_MSC_VER) || defined(__BORLANDC__) 
// Enable the use of isnan().
#include <float.h>
#ifndef isnan
#define isnan _isnan
#endif // isnan
// Used to define NAN (Not A Number).
#ifdef NAN_CONSTS_NEEDED
const unsigned long nan[2] = {0xffffffff, 0xffffffff};
const double nan_double = -*(double*)nan;
#endif // NAN_CONSTS_NEEDED
#endif // defined(_MSC_VER) || defined(__BORLANDC__) 

#ifdef NAI_CONST_NEEDED
// Used to define NAI (Not An Interval).
const interval nai = interval();
#endif // NAI_CONST_NEEDED

#include <cmath>

using namespace std;
using namespace ibex;

//----------------------------------------------------------------------
// Operators
//----------------------------------------------------------------------
std::ostream& operator<<(std::ostream& os, const interval_empty_default& a)
{
	if (a.is_empty()) os << "EmptyInterval";
	else if (a.lb() != a.ub())
	{ 
		os << "[" << setprecision(4) << a.lb() << ", " << setprecision(4) << a.ub() << "] "; 
	}
	else os << a.lb();
	return os;
}
//----------------------------------------------------------------------
// Member functions
//----------------------------------------------------------------------
interval_empty_default& interval_empty_default::Intersect(const interval_empty_default& Y) 
{ 
	interval_empty_default X = *this; interval_empty_default Z = Inter(X, Y); *this = Z; return *this; 
}
//----------------------------------------------------------------------
double Min(vector<double>& x)
{	double d = oo;
	for (unsigned int i = 0; i < x.size(); i++)
		d = min(x[i], d);
	return d;
}
//----------------------------------------------------------------------
double Max(vector<double>& x)
{
	double d = -oo;
	for (unsigned int i = 0; i < x.size(); i++)
		d = max(x[i], d);
	return d;
}
//----------------------------------------------------------------------
double Sign(const double x)
{
	if (x > 0) return (1);
	else return (0);
}
//----------------------------------------------------------------------
double Chi(const double a, const double b, const double c)
{
	if (a < 0) return (b);
	else return (c);
}
//----------------------------------------------------------------------
double Det(double ux, double uy, double vx, double vy)
{
	return (ux*vy - vx*uy);
}

//----------------------------------------------------------------------
double Rand(interval X)
{	double a=1.0*(rand()%10000)/10000;
        return (1-a)*X.lb()+a*X.ub() ;}
//----------------------------------------------------------------------
double DistanceDirSegment(double mx, double my, double theta, double ax, double ay, double bx, double by)
{      
	// Distance directionnelle du point m au segment [a,b]. La direction est donnee par theta
	double ma_x = ax - mx;
	double ma_y = ay - my;
	double mb_x = bx - mx;
	double mb_y = by - my;
	double ab_x = bx - ax;
	double ab_y = by - ay;
	double ux = cos(theta);
	double uy = sin(theta);
	double z1 = Det(ma_x, ma_y, ux, uy);
	double z2 = Det(ux, uy, mb_x, mb_y);
	double z3 = Det(ma_x, ma_y, ab_x, ab_y);
	double z4 = Det(ux, uy, ab_x, ab_y);
	double z5 = min(z1, min(z2, z3));
	double d1 = 0;
	if (z4 == 0) d1 = oo;
	else d1 = z3 / z4;
	return (Chi(z5, oo, d1));
}
//----------------------------------------------------------------------
void DistanceDirSegment(double& d,double& phi, double mx, double my, double theta, double ax, double ay, double bx, double by)
{      
	// Distance directionnelle du point m au segment [a,b].
	double ma_x=ax-mx;
	double ma_y=ay-my;
	double mb_x=bx-mx;
	double mb_y=by-my;
	double ab_x=bx-ax;
	double ab_y=by-ay;
	double ux=cos(theta);
	double uy=sin(theta);
	double z1=Det(ma_x,ma_y,ux,uy);
	double z2=Det(ux,uy,mb_x,mb_y);
	double z3=Det(ma_x,ma_y,ab_x,ab_y);
	double z4=Det(ux,uy,ab_x,ab_y);
	double z5=min(z1,min(z2,z3));
	double d1=z3/z4;
	d=Chi(z5,oo,d1);
	phi=atan2(-ab_x,ab_y); //phi is the angle of the normal vector of [a,b]
}
//----------------------------------------------------------------------
double DistanceDirSegments(double mx, double my, double theta, vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by)
{      
	// Distance directionnelle relativement a un polygone
	vector<double> dist(ax.size());
	for (unsigned int j = 0; j < ax.size(); j++)
		dist[j] = DistanceDirSegment(mx, my, theta, ax[j], ay[j], bx[j], by[j]);
	double distmin = Min(dist);
	return (distmin);
}
//----------------------------------------------------------------------
void DistanceDirSegments(double& d, double& phi, double mx, double my, double theta, vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by)
{     
	// Distance directionnelle relativement a un polygone (le triedre m-a-b doit etre direct, sinon cette distance est infinie)
	d = oo;
	for (unsigned int j = 0; j < ax.size(); j++)
	{
		double dj,phij;
		DistanceDirSegment(dj,phij,mx,my,theta,ax[j],ay[j],bx[j],by[j]);
		if (dj < d) { d=dj;phi=phij; };
	}
}
//----------------------------------------------------------------------
double DistanceDirCircle(double mx, double my, double theta, double cx, double cy, double r)
{      
	// Distance directionnelle du point m au cercle de centre c et de rayon r.
	double ux, uy, alpha, beta, a, b, c, delta, px1, py1, px2, py2, d1, d2;
	ux = cos(theta);
	uy = sin(theta);
	if (fabs(uy) > 0.00)  //pour eviter la division par zero. Il conviendrait de traiter le cas uy=0
	{
		alpha = ux / uy;
		beta = mx - my*alpha;
		a = alpha*alpha + 1;
		b = 2 * alpha*(beta - cx) - 2 * cy;
		c = (beta - cx)*(beta - cx) + cy*cy - r*r;
		delta = b*b - 4 * a*c;
		if (delta < 0) return(oo);
		py1 = (-b - sqrt(delta)) / (2 * a);
		px1 = alpha*py1 + beta;
		py2 = (-b + sqrt(delta)) / (2 * a);
		px2 = alpha*py2 + beta;
		d1 = Chi((px1 - mx)*ux + (py1 - my)*uy, oo, sqrt((px1 - mx)*(px1 - mx) + (py1 - my)*(py1 - my)));
		d2 = Chi((px2 - mx)*ux + (py2 - my)*uy, oo, sqrt((px2 - mx)*(px2 - mx) + (py2 - my)*(py2 - my)));
		return min(d1, d2);
	}
	return oo;
}
//----------------------------------------------------------------------
void DistanceDirCircle(double& d, double& phi, double mx, double my, double theta, double cx, double cy, double r)
{      
	//  d is the directional distance from m to the circle of center c and radius r.
	//  phi is the angle of the normal vector of of the impact point
	double ux,uy,alpha,beta,a,b,c,delta,px1,py1,px2,py2,d1,d2,phi1,phi2;
	ux=cos(theta);
	uy=sin(theta);
	if (fabs(uy)<0.01)    // pour eviter une division par zero, on permutte x et y
	{ double aux; aux=mx; mx=my; my=aux;  aux=cx; cx=cy; cy=aux;    aux=ux; ux=uy; uy=aux; }
	alpha=ux/uy;
	beta=mx-my*alpha;
	a=alpha*alpha+1;
	b=2*alpha*(beta-cx)-2*cy;
	c=(beta-cx)*(beta-cx)+cy*cy-r*r;
	delta=b*b-4*a*c;
	if (delta<0) {d=oo;phi=-999;return;};
	py1=(-b-sqrt(delta))/(2*a);
	px1=alpha*py1+beta;
	phi1=atan2(cy-py1,cx-px1);
	py2=(-b+sqrt(delta))/(2*a);
	px2=alpha*py2+beta;
	phi2=atan2(cy-py2,cx-px2);

	d1=oo;d2=oo;
	if ((px1-mx)*ux+(py1-my)*uy>0)   d1=sqrt((px1-mx)*(px1-mx)+(py1-my)*(py1-my));
	if ((px2-mx)*ux+(py2-my)*uy>0)   d2=sqrt((px2-mx)*(px2-mx)+(py2-my)*(py2-my));
	if (d1<d2) { d=d1; phi=phi1; } else { d=d2; phi=phi2; }
}
//----------------------------------------------------------------------
double DistanceDirCircles(double mx, double my, double theta, vector<double> cx, vector<double> cy, vector<double> r)
{      
	// Distance directionnelle relativement a plusieurs cercles.
	vector<double> dist(cx.size());
	for (unsigned int j = 0; j < cx.size(); j++)
                dist[j] = DistanceDirCircle(mx, my, theta, cx[j], cy[j], r[j]);
	double distmin = Min(dist);
	return (distmin);
}
//----------------------------------------------------------------------
void DistanceDirCircles(double& d,double& phi, double mx, double my, double theta, vector<double> cx, vector<double> cy, vector<double> r)
{       
	d = oo;
	for (unsigned int j = 0; j < cx.size(); j++)
	{
		double dj, phij;
		DistanceDirCircle(dj,phij,mx,my,theta,cx[j],cy[j],r[j]);
		if (dj < d) { d=dj; phi=phij; };
	}
}
//----------------------------------------------------------------------
double DistanceDirSegmentsOrCircles(double mx, double my, double theta,
									vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by,
									vector<double> cx, vector<double> cy, vector<double> r)
{     
	double d1a, d1b;
	d1a = DistanceDirSegments(mx, my, theta, ax, ay, bx, by);
	d1b = DistanceDirCircles(mx, my, theta, cx, cy, r);
	return min(d1a, d1b);
}
//----------------------------------------------------------------------
void DistanceDirSegmentsOrCircles(double& d, double& phi, double mx, double my, double theta,
								  vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by,
								  vector<double> cx, vector<double> cy, vector<double> r)
{     
	// returns the distance and orientation collected by a laser rangefinder in a room made with segments and circles
	double phi1a, phi1b, d1a, d1b;
	DistanceDirSegments(d1a,phi1a,mx,my,theta,ax,ay,bx,by);
	DistanceDirCircles(d1b,phi1b,mx,my,theta,cx,cy,r);
	if (d1a < d1b) { d = d1a; phi = phi1a-theta; } else { d = d1b; phi = phi1b-theta; }
}
//----------------------------------------------------------------------
// Interval-valued functions
//----------------------------------------------------------------------
interval Min(const interval& x, const interval& y)
{
	return ibex::min(x, y);
}
//----------------------------------------------------------------------
interval Min(const interval& x, const interval& y, const interval& z)
{
	return ibex::min(ibex::min(x, y), z);
}
//----------------------------------------------------------------------
interval Max(const interval& x, const interval& y)
{
	return ibex::max(x, y);
}
//----------------------------------------------------------------------
interval Max(const interval& x, const interval& y, const interval& z)
{
	return ibex::max(ibex::max(x, y), z);
}
//----------------------------------------------------------------------
interval InvSqrt(interval& X)
{
	interval Y(-oo,oo);
	bwd_sqrt(Y, X);
	return Y;
}
//----------------------------------------------------------------------------
interval Det(interval& ux, interval& uy, interval& vx, interval& vy)
{
	return (ux*vy - vx*uy);
}
//----------------------------------------------------------------------
interval Det(interval& ux, interval& uy, double& vx, double& vy)
{
	return (vy*ux - vx*uy);
}
//----------------------------------------------------------------------
interval Step(const interval& X)
{ 
	if (X.is_empty()) return interval();
	if (X.lb() > 0) return (interval(1));
	if (X.ub() < 0) return (interval(0));
	return (interval(0,1));
}
//----------------------------------------------------------------------
interval Inter(const interval& a, const interval& b)
{       return a&b; }
//----------------------------------------------------------------------
interval Inter(vector<interval> x)
{	interval r = interval::EMPTY_SET; // r is empty
	for (unsigned int i = 0 ; i < x.size(); i++)
		r = r&x[i];
	return r;
}
//----------------------------------------------------------------------
interval Union(const interval& a, const interval& b)
{   
	return a|b;
}
//----------------------------------------------------------------------
interval Union(vector<interval> x)
{   
	interval r = interval::EMPTY_SET; // r is empty
	for (unsigned int i = 0 ; i < x.size(); i++)
		r = r|x[i];
	return r;
}
//----------------------------------------------------------------------
interval Inflate(const interval& a, double eps)
{
	interval r(a.lb() - eps, a.ub() + eps);
	if (a.is_empty()) r = interval::EMPTY_SET; // r is empty
	return interval(r);
}
//----------------------------------------------------------------------
// Other functions
//----------------------------------------------------------------------
double Inf(const interval& a)
{
	return a.lb();
}
//-----------------------------------------------------------------------
double Sup(const interval& a)
{
	return a.ub();
}
//-----------------------------------------------------------------------
double Center(const interval& a)
{
	return a.mid();
}
//-----------------------------------------------------------------------
double Width(const interval& a)
{       return Sup(a)-Inf(a);
}
//----------------------------------------------------------------------
double Volume(const interval& a)
{
	return a.diam();
}
//----------------------------------------------------------------------
double Rad(const interval& a)
{
	return a.rad();
}
//----------------------------------------------------------------------
bool Disjoint(const interval& a, const interval& b)
{
	if (a.is_empty()||b.is_empty()) return true;
	return ((a.ub() < b.lb())||(b.ub() < a.lb()));
}
//----------------------------------------------------------------------
bool Subset(const interval& a, const interval& b)
{
	if (a.is_empty()) return true;
	if (b.is_empty()) return false;
	return ((a.lb() >= b.lb())&&(a.ub() <= b.ub()));
}
//----------------------------------------------------------------------
bool In(double a, const interval& b)
{
        if (b.is_empty()) return false;
        return ((b.lb() <= a)&&(a <= b.ub()));
}
//----------------------------------------------------------------------
// Contractors
//----------------------------------------------------------------------
void Cadd(interval& Z, interval& X, interval& Y)
{
        Z &= X+Y; bwd_add(Z,X,Y);
}
//----------------------------------------------------------------------
void Cadd(interval& Z, double x, interval& Y)
{
        Z &= x+Y;
        interval X(x);
        bwd_add(Z,X,Y);
}
//----------------------------------------------------------------------
void Cadd(interval& Z, interval& X, double y)
{
        Z &= X+y;
        interval Y(y);
        bwd_add(Z,X,Y);
}
//----------------------------------------------------------------------
void Csub(interval& Z, interval& X, interval& Y)
{
        Z &= X-Y;
        bwd_sub(Z,X,Y);
}
//----------------------------------------------------------------------
void Csub(interval& Z, double x, interval& Y)
{
        Z &= x-Y;
        interval X(x);
        bwd_sub(Z,X,Y);
}
//----------------------------------------------------------------------
void Csub(interval& Z, interval& X, double y)
{
        Z &= X-y;
        interval Y(y);
        bwd_sub(Z,X,Y);
}
//----------------------------------------------------------------------
void Cmul(interval& Z, interval& X, interval& Y)
{
        Z &= X*Y;
        bwd_mul(Z,X,Y);
}
//----------------------------------------------------------------------
void Cmul(interval& Z, double x, interval& Y)
{
        Z &= x*Y;
        interval X(x);
        bwd_mul(Z,X,Y);

}
//----------------------------------------------------------------------
void Cmul(interval& Z, interval& X, double y)
{
        Z &= X*y;
        interval Y(y);
        bwd_mul(Z,X,Y);
}
//----------------------------------------------------------------------
void Cdiv(interval& Z, interval& X, interval& Y)
{
        Z &= X/Y;
        bwd_div(Z,X,Y);
}
//----------------------------------------------------------------------
void Cequal(interval& Y, interval& X)
{
	Y = Inter(Y, X); X = Y;
}
//----------------------------------------------------------------------
void Cmin(interval& a, interval& b, interval& c)
{   
        a = Inter(a, Min(b, c));
        if (Disjoint(a, b)) c = Inter(c, a);
        else { if (Disjoint(a, c)) b = Inter(b, a); }
        interval temp(a.lb(), oo);
        b = Inter(b, temp); c = Inter(c, temp);
}
//----------------------------------------------------------------------
void Cmin(interval& a, interval& b, interval& c, interval& d)
{
	// contrainte quaternaire   a=min(b,c,d)
	interval z1 = Min(b, c);
        Cmin(a, z1, d);
        Cmin(a, z1, d);
        Cmin(z1, b, c);
}
//----------------------------------------------------------------------
void Cmin(interval& a, interval& b, interval& c, interval& d, interval& e)
{   
	interval z1 = Min(b,c,d);
        Cmin(a, z1, e);
        Cmin(a, z1, e);
        Cmin(z1, b, c, d);
}
//----------------------------------------------------------------------
int Cmin(interval& a, vector<interval>& x)
{
	vector<interval> z(x.size());
	z[0] = x[0];
	for (unsigned int i = 1; i < x.size(); i++)
        {	z[i] = interval(-oo, oo);
                Cmin(z[i], x[i], z[i - 1]);
	}
	Cequal(a, z[x.size() - 1]);
	for (int i = (int)(x.size() - 1); i >= 1; i--)
        {       Cmin(z[i], x[i], z[i - 1]);;
		if (z[i].is_empty()) return -1;
	}
	return 1;
}
//----------------------------------------------------------------------
void Cmax(interval& a, interval& b, interval& c)
{
        a = Inter(a, Max(b, c));
        if (Disjoint(a, b)) c = Inter(c, a);
        else { if (Disjoint(a, c)) b = Inter(b, a); }
        interval temp(-oo, a.ub());
        b = Inter(b, temp); c = Inter(c, temp);

}
//----------------------------------------------------------------------
void Cabs(interval& Y, interval& X)
{
	interval Xd = Inter(X, interval(0, oo)), Xg = Inter(X, interval(-oo, 0));
        interval Yd = Inter(Y, Xd), Yg = Inter(Y, -Xg); Y = Union(Yd, Yg);
        Xd = Inter(Xd, Y); Xg = Inter(Xg, -Y); X = Union(Xd, Xg);
}
//----------------------------------------------------------------------
void Cchi(interval& F, interval& A, interval& B, interval& C)
{
	if (A.ub() < 0) { Cequal(B, F); }
	else if (A.lb() > 0) { Cequal(C, F); };
	if (Disjoint(F, B)) { A = Inter(A, interval(0, oo)); };
	if (Disjoint(F, C)) { A = Inter(A, interval(-oo, 0)); };
	F = Union(Inter(F, B), Inter(F, C));
}
//----------------------------------------------------------------------
void Cgeq(interval& Y, interval& X)
{
	if (Y.lb() >= X.ub()) return;
	interval Z = Y - X;
	Z = Inter(Z, interval(0, oo));
        Csub(Z, Y, X);
}
//----------------------------------------------------------------------
void Cinteger(interval& X)
{
	interval A1(ceil(X.lb()), floor(X.ub()));
	X = Inter(X, A1);
}
//----------------------------------------------------------------------
void Cboolean(interval& X)
{
	interval nul(0), one(1);
	interval A1 = Inter(nul, X);
	interval A2 = Inter(one, X);
	interval A = Union(A1, A2);
	X = Inter(X, A);
}
//----------------------------------------------------------------------
void Csqr(interval &Y, interval &X)
{   
        Y &= sqr(X);
        bwd_sqr(Y,X);
}
//-----------------------------------------------------------------------
void Cexp(interval& Y, interval& X)
{   
        Y &= exp(X);
        bwd_exp(Y,X);
}
//----------------------------------------------------------------------
void Clog(interval& Y, interval& X)
{   
        Y &= log(X);
        bwd_log(Y,X);
}
//----------------------------------------------------------------------
void Cpow(interval& Y, interval& X, int n) 
{
        Y &= pow(X,n);
        bwd_pow(Y,n ,X);
}
//----------------------------------------------------------------------
void Ccos(interval &Y, interval &X)
{   
        Y &= cos(X);
        bwd_cos(Y,X);
}
//-----------------------------------------------------------------------
void Csin(interval &Y, interval &X)
{   
        Y &= sin(X);
        bwd_sin(Y,X);
}
//-----------------------------------------------------------------------
void Ctan(interval& Y, interval& X)
{
        Y &= tan(X);
        bwd_tan(Y,X);
        Y &= tan(X); //On est oblige de le rappeler car sinon par optimal
}
//----------------------------------------------------------------------
void Cnorm(interval& N, interval& X, interval& Y)
{
	interval SqrX, SqrY, SqrN;
	SqrX = Sqr(X);	SqrY = Sqr(Y);	SqrN = Sqr(N);
        Cadd(SqrN, SqrX, SqrY);   Cadd(SqrN, SqrX, SqrY);
        Csqr(SqrY, Y);        Csqr(SqrX, X);        Csqr(SqrN, N);
}
//----------------------------------------------------------------------
void Cnorm(interval& N, interval& X, interval& Y, interval& Z)
{
	interval SqrX = Sqr(X), SqrY = Sqr(Y), SqrZ = Sqr(Z), SqrN = Sqr(N);
        SqrN = Inter(SqrN, SqrX+SqrY+SqrZ);
        Csqr(SqrN, N);
        SqrZ = Inter(SqrZ, SqrN-SqrX-SqrY);
        SqrY = Inter(SqrY, SqrN-SqrX-SqrZ);
        SqrX = Inter(SqrX, SqrN-SqrZ-SqrY);
        Csqr(SqrZ, Z);
        Csqr(SqrY, Y);
        Csqr(SqrX, X);
}
//----------------------------------------------------------------------
void Cdist(interval& R, interval& X1, interval& Y1, interval& X2, interval& Y2)
{
	static Variable x1, y1, x2, y2, r;
	static NumConstraint C(x1,y1,x2,y2,r, pow(x1-x2,2)+pow(y1-y2,2)=pow(r,2));
	static CtcHC4 ctc(C);
	IntervalVector P(5);
	P[0]=X1; P[1]=Y1; P[2]=X2; P[3]=Y2; P[4]=R;
	try {ctc.contract(P);} catch(EmptyBoxException) {X1=interval::EMPTY_SET;Y1=X1; X2=X1; Y2=X1; R=X1; return;}
	X1=P[0]; Y1=P[1]; X2=P[2]; Y2=P[3]; R=P[4];
	return;
}
//----------------------------------------------------------------------
void Cscal(interval& s, interval& ux, interval& uy, interval& vx, interval& vy)
{	interval z1 = ux*vx;
	interval z2 = uy*vy;
	Cadd(s, z1, z2);
        Cmul(z2, uy, vy);
        Cmul(z1, ux, vx);
}
//----------------------------------------------------------------------
void Cscal(interval& s, double& ux, double& uy, interval& vx, interval& vy)
{	interval z1 = ux*vx;
	interval z2 = uy*vy;
	Cadd(s, z1, z2);
        Cmul(z2, uy, vy);
        Cmul(z1, ux, vx);
}
//----------------------------------------------------------------------
void Cdet(interval& det, interval& ux, interval& uy, interval& vx, interval& vy)
{	interval z1 = ux*vy;
	interval z2 = vx*uy;
        Csub(det, z1, z2);
        Csub(det, z1, z2);
        Cmul(z2, vx, uy);
        Cmul(z1, ux, vy);
}
//----------------------------------------------------------------------
void Cdet(interval& det, double& ux, double& uy, interval& vx, interval& vy)
{
	interval z1 = ux*vy;
	interval z2 = uy*vx;
        Csub(det, z1, z2);
        Csub(det, z1, z2);
        Cmul(z2, uy, vx);
        Cmul(z1, ux, vy);
}
//----------------------------------------------------------------------
void Cdet(interval& det, interval& ux, interval& uy, double& vx, double& vy)
{
	interval z1 = vy*ux;
	interval z2 = vx*uy;
        Csub(det, z1, z2);
        Csub(det, z1, z2);
        Cmul(z2, vx, uy);
        Cmul(z1, vy, ux);
}
//----------------------------------------------------------------------
void CDistanceDirLine(interval& dist, interval& mx, interval& my, interval& theta, double& ax, double& ay, double& bx, double& by)
{     
	// la distance dist entre le point m=(mx,my) a la droite [a,b] suivant le vecteur u
	if ((dist.is_empty())||(mx.is_empty())||(my.is_empty())||(theta.is_empty()))
	{  
		dist = interval(); mx = interval(); my = interval(); theta = interval(); return;
	}
	interval ma_x=ax-mx; interval ma_y=ay-my;
	interval mb_x=bx-mx; interval mb_y=by-my;
	double ab_x=bx-ax; double ab_y=by-ay;
	interval ux=Cos(theta); interval uy=Sin(theta);
	interval z3=Det(ma_x,ma_y,ab_x,ab_y);
	interval z4=Det(ux,uy,ab_x,ab_y);
	dist=Inter(dist,z3/z4);

        Cdiv(dist,z3,z4);
        Cdet(z4,ux,uy,ab_x,ab_y);
        Cdet(z3,ma_x,ma_y,ab_x,ab_y);
        Csin(uy,theta); Ccos(ux,theta);
        Csub(mb_y,by,my); Csub(mb_x,bx,mx);
        Csub(ma_y,ay,my); Csub(ma_x,ax,mx);
}
//----------------------------------------------------------------------
int CDistanceDirSegment(interval& dist, interval& mx, interval& my, interval& theta, double ax, double ay, double bx, double by)
{
	// la distance dist entre le point m=(mx,my) au segment [a,b] suivant le vecteur u
	if ((dist.is_empty())||(mx.is_empty())||(my.is_empty())||(theta.is_empty()))
	{
		dist = interval(); mx = interval(); my = interval(); theta = interval(); return -1;
		//dist = interval(); mx = interval(); my = interval(); theta = interval(); return;
	}
	interval ma_x = ax - mx; interval ma_y = ay - my;
	interval mb_x = bx - mx; interval mb_y = by - my;
	double ab_x = bx - ax; double ab_y = by - ay;
	interval ux = Cos(theta); interval uy = Sin(theta);
	interval z1 = Det(ma_x, ma_y, ux, uy);
	interval z2 = Det(ux, uy, mb_x, mb_y);
	interval z3 = Det(ma_x, ma_y, ab_x, ab_y);
	interval z4 = Det(ux, uy, ab_x, ab_y);
	interval z5 = Min(z1, z2, z3);
	interval d1 = z3 / z4;
	interval infty = interval(oo, oo);
	Cchi(dist, z5, infty, d1);
	//REDONDANT
	Cchi(dist, z3, infty, dist);
	Cchi(dist, z4, infty, dist);

        Cdiv(d1, z3, z4);
        Cmin(z5, z1, z2, z3);
        Cdet(z4, ux, uy, ab_x, ab_y);
        Cdet(z3, ma_x, ma_y, ab_x, ab_y);
        Cdet(z2, ux, uy, mb_x, mb_y);
        Cdet(z1, ma_x, ma_y, ux, uy);
        Csin(uy, theta); Ccos(ux, theta);
        Csub(mb_y, by, my); Csub(mb_x, bx, mx);
        Csub(ma_y, ay, my); Csub(ma_x, ax, mx);
	return 1;
}
//----------------------------------------------------------------------
void CDistanceDirSegments(interval& distmin, interval& mx, interval& my, interval& theta, vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by)
{      
	// Distance directionnelle relativement a un polygone
	vector<interval> dist(ax.size());
	for (unsigned int j = 0; j < ax.size(); j++) dist[j] = interval(0, oo);
	for (unsigned int j = 0; j < ax.size(); j++)
                CDistanceDirSegment(dist[j], mx, my, theta, ax[j], ay[j], bx[j], by[j]);
        Cmin(distmin, dist);
	for (int j = (int)(ax.size() - 1); j >= 0; j--)
                CDistanceDirSegment(dist[j], mx, my, theta, ax[j], ay[j], bx[j], by[j]);
}
//----------------------------------------------------------------------
void CinRing(interval& X, interval& Y, double cx, double cy, interval R)
{
	static Variable x, y, r, vcx, vcy;
	static NumConstraint C(x,y,r,vcx,vcy, pow(x-vcx,2)+pow(y-vcy,2)=pow(r,2));
	static CtcHC4 ctc(C);
	IntervalVector P(5); P[0]=X; P[1]=Y; P[2]=R; P[3]=cx; P[4]=cy;
	try {ctc.contract(P);} catch(EmptyBoxException) {X=interval::EMPTY_SET;Y=X; return;}
	X=P[0]; Y=P[1];
	return;
}
//----------------------------------------------------------------------
void CinLine(interval& mx, interval& my, double& ax, double& ay, double& bx, double& by)
{ 	// contracte relativement a la contrainte : "m appartient a la droite (a,b)"
	interval ma_x=ax-mx;
	interval ma_y=ay-my;
	double ab_x=bx-ax;
	double ab_y=by-ay;
	interval z1=interval(0,0);
        Cdet(z1,ab_x,ab_y,ma_x,ma_y);
        Csub(ma_y,ay,my);
        Csub(ma_x,ax,mx);
	if ((mx.is_empty())||(my.is_empty())) { mx = interval(); my = interval(); }
}
//----------------------------------------------------------------------
void CinSegment(interval& mx, interval& my, double ax, double ay, double bx, double by)
{	// contracte relativement a la contrainte : "m appartient au segment [a,b]"
	mx = Inter(mx, Union(interval(ax, ax), interval(bx, bx)));
	my = Inter(my, Union(interval(ay, ay), interval(by, by)));
	interval ma_x = ax - mx;
	interval ma_y = ay - my;
	double ab_x = bx - ax;
	double ab_y = by - ay;
	interval z1 = interval(0, 0);
        Cdet(z1, ab_x, ab_y, ma_x, ma_y);
        Csub(ma_y, ay, my);
        Csub(ma_x, ax, mx);
	if ((mx.is_empty())||(my.is_empty())) { mx = interval::EMPTY_SET; my = interval::EMPTY_SET; }
}
//----------------------------------------------------------------------
void CinSegments(interval& mx, interval& my, vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by)
{      
	// contracte relativement a la contrainte : "m appartient au polygone dont les segments sont les [ai,bi]"
	if ((ax.size() == 0)||(ay.size() == 0)||(bx.size() == 0)||(by.size() == 0)) return;
	vector<interval> Mx(ax.size());
	vector<interval> My(ax.size());
	for (unsigned int j = 0; j < ax.size(); j++)
	{
		interval mx0 = mx;
		interval my0 = my;
                CinSegment(mx0, my0, ax[j], ay[j], bx[j], by[j]);
		Mx[j] = mx0;
		My[j] = my0;
	}
	mx = Inter(mx, Union(Mx));
	my = Inter(my, Union(My));
}
//----------------------------------------------------------------------
void CinCircle(interval& X, interval& Y, double cx, double cy, double r)
{
        static Variable x, y, vr, vcx, vcy;
        static NumConstraint C(x,y,vr,vcx,vcy, pow(x-vcx,2)+pow(y-vcy,2)=pow(vr,2));
        static CtcHC4 ctc(C);
        IntervalVector P(5); P[0]=X; P[1]=Y; P[2]=r; P[3]=cx; P[4]=cy;
        try {ctc.contract(P);} catch(EmptyBoxException) {X=interval::EMPTY_SET;Y=X; return;}
        X=P[0]; Y=P[1];
        return;
}

//----------------------------------------------------------------------
void SinCircles(interval& mx, interval& my, vector<double> cx, vector<double> cy, vector<double> r, bool truth)
{      
	// contracte relativement a la contrainte : "m appartient a un des cercles de centre ci et de rayon ri"
	if ((cx.size() == 0)||(cy.size() == 0)||(r.size() == 0)) return;
	vector<interval> Mx(cx.size());
	vector<interval> My(cx.size());
	for (unsigned int j = 0; j < cx.size(); j++)
        {       interval mx0 = mx;
		interval my0 = my;
                CinCircle(mx0, my0, cx[j], cy[j], r[j]);
		Mx[j] = mx0;
		My[j] = my0;
	}
	if (truth)
	{
		mx = Inter(mx, Union(Mx));
		my = Inter(my, Union(My));
	} 
	else 
	{
		mx = Inter(mx, Inter(Mx));
		my = Inter(my, Inter(My));
	}
}
//----------------------------------------------------------------------
void CinSegmentsOrCircles(interval& mx, interval& my, vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by, vector<double> cx, vector<double> cy, vector<double> r)
{      
	// contracte relativement a la contrainte : "m appartient soit au polygone soit a un des cercles de centre ci et de rayon ri"
	interval mx1 = mx;
	interval my1 = my;
        CinSegments(mx1, my1, ax, ay, bx, by);
	interval mx2 = mx;
	interval my2 = my;
        SinCircles(mx2, my2, cx, cy, r);
	if ((ax.size() > 0)&&(cx.size() > 0))            //cercles et polygone
	{
		mx = Inter(mx, Union(mx1, mx2)); 
		my = Inter(my, Union(my1, my2));
	}
	if ((ax.size() > 0)&&(cx.size() == 0))           // pas de cercles
	{
		mx = mx1; my = my1;
	}
	if ((ax.size() == 0)&&(cx.size() >= 0))          // pas de segments
	{
		mx = mx2; my = my2;
	}
}
//----------------------------------------------------------------------
void CPoseInSegment(interval& mx, interval& my, interval& phi, double& ax, double& ay, double& bx, double& by)
{     
	// contracte relativement a "la pose (m,phi) appartient au segment [a,b]"
        CinSegment(mx,my,ax,ay,bx,by);
	double ab_x=bx-ax; //(bx-ax)*cos(phi)+(by-ay)*sin(phi)=0
	double ab_y=by-ay;
	interval cphi=Cos(phi);
	interval sphi=Sin(phi);
	interval scal=interval(-0.0,0.0);   CScal(scal,ab_x,ab_y,cphi,sphi); // phi orthogonal to the line
	interval det(0,oo);
	Cdet(det,cphi,sphi,ab_x,ab_y); // select the right direction
        Ccos(cphi,phi);
        Csin(sphi,phi);
	if ((mx.is_empty())||(my.is_empty())||(phi.is_empty())) { mx = interval(); my = interval(); phi = interval(); }
}
//----------------------------------------------------------------------
void CPoseInSegments(interval& mx, interval& my, interval& phi,vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by)
{
	if ((ax.size() == 0)||(ay.size() == 0)||(bx.size() == 0)||(by.size() == 0)) return;
	vector<interval> Mx(ax.size());
	vector<interval> My(ax.size());
	vector<interval> Mphi(ax.size());
	for (unsigned int j=0;j<ax.size();j++)
	{ 
		interval mx0=mx;
		interval my0=my;
		interval phi0=phi;
		CPoseInSegment(mx0,my0,phi0,ax[j],ay[j],bx[j],by[j]);
		Mx[j]=mx0;
		My[j]=my0;
		Mphi[j]=phi0;
	}
	mx=Inter(mx,Union(Mx));
	my=Inter(my,Union(My));
	phi=Inter(phi,Union(Mphi));
}
//----------------------------------------------------------------------
void CPoseInCircle(interval& mx, interval& my, interval& phi, double& cx, double& cy, double& r)
{ 
        CinCircle(mx,my,cx,cy,r);
	interval mc_x=cx-mx;
	interval mc_y=cy-my;
	interval cphi=Cos(phi);
	interval sphi=Sin(phi);
	interval scal=interval(0,oo);   CScal(scal,mc_x,mc_y,cphi,sphi); //scal(cphi,sphi,cx_mx,cy-my)>0
	interval det(0,0);
	Cdet(det,cphi,sphi,mc_x,mc_y); //det(cphi,sphi,cx_mx,cy-my)=0
        Ccos(cphi,phi);
        Csin(sphi,phi);
	if ((mx.is_empty())||(my.is_empty())||(phi.is_empty())) { mx = interval(); my = interval(); phi = interval(); }
}
//----------------------------------------------------------------------
void CPoseInCircles(interval& mx, interval& my, interval& phi,vector<double> cx, vector<double> cy, vector<double> r)
{ 
	if ((cx.size() == 0)||(cy.size() == 0)||(r.size() == 0)) return;
	vector<interval> Mx(cx.size());
	vector<interval> My(cx.size());
	vector<interval> Mphi(cx.size());
	for (unsigned int j=0;j<cx.size();j++)
        {	interval mx0=mx;
		interval my0=my;
		interval phi0=phi;
		CPoseInCircle(mx0,my0,phi0,cx[j],cy[j],r[j]);
		Mx[j]=mx0;
		My[j]=my0;
		Mphi[j]=phi0;
	}
	mx=Inter(mx,Union(Mx));
	my=Inter(my,Union(My));
	phi=Inter(phi,Union(Mphi));
}
//----------------------------------------------------------------------
void CPoseInSegmentsOrCircles(interval& mx,interval& my,interval& malpha, vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by,
							  vector<double> cx, vector<double> cy, vector<double> r)
{      	// "la pose (m,alpha) appartient soit au polygone soit a un des cercles de centre ci et de rayon ri"
	vector<interval> Mx(2);
	vector<interval> My(2);
	vector<interval> Malpha(2);
	interval mx0=mx;
	interval my0=my;
	interval malpha0=malpha;
	CPoseInSegments(mx0,my0,malpha0,ax,ay,bx,by);
	Mx[0]=mx0;   My[0]=my0;  Malpha[0]=malpha0;
	mx0=mx;  my0=my;  malpha0=malpha;
	CPoseInCircles(mx0,my0,malpha0,cx,cy,r);
	Mx[1]=mx0;
	My[1]=my0;
	Malpha[1]=malpha0;
	mx=Inter(mx,Union(Mx));
	my=Inter(my,Union(My));
	malpha=Inter(malpha,Union(Malpha));
        // Bug if no segments or no circles, see CinSegmentsOrCircles()...?
}
//----------------------------------------------------------------------
void CPoseTrans(interval& qx,interval& qy,interval& dist, interval& px, interval& py, interval& alpha)
{ 
	interval ux=Cos(alpha);
	interval uy=Sin(alpha);
	interval dx=ux*dist;
	interval dy=uy*dist;
	interval pxdx = px+dx;
	interval pydy = py+dy;
	qx=Inter(qx,pxdx);
	qy=Inter(qy,pydy);
        Cadd(qy,py,dy);
        Cadd(qx,px,dx);
        Cmul(dy,uy,dist);
        Cmul(dx,ux,dist); // Plante ici car dx est negatif et d > 0 et ux > 0?
        Csin(uy,alpha);
        Ccos(ux,alpha);
	if ((qx.is_empty())||(qy.is_empty())||(alpha.is_empty())) { qx = interval(); qy = interval(); alpha = interval(); }
}
//----------------------------------------------------------------------
void CPoseTransRot(interval& qx,interval& qy,interval& beta,interval& d, interval& psi, interval& px,interval& py,interval& alpha)
{  
	CPoseTrans(qx,qy,d, px, py,alpha);
	Cadd(beta,psi,alpha);
}
//----------------------------------------------------------------------
void CPoseRotTrans(interval& qx,interval& qy,interval& beta,interval& phi,interval& d, interval& px,interval& py,interval& alpha)
{   
	Cadd(beta,phi,alpha);
	CPoseTrans(qx,qy,d, px, py,beta);
	Cadd(beta,phi,alpha);
}
//----------------------------------------------------------------------
void CPoseRotTransRot(interval& qx,interval& qy,interval& beta,interval& phi,interval& d, interval& psi, interval& px,interval& py,interval& alpha)
{   
	interval delta=alpha+phi;
	CPoseTransRot(qx,qy,beta,d,psi,px,py,delta);
        Cadd(delta,alpha,phi);
}
//----------------------------------------------------------------------
void CPoseTransInWallsOrCircles(interval& px, interval& py, interval& alpha, interval& d,
								vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by,
								vector<double> cx, vector<double> cy, vector<double> r)
{  	interval qx=interval(-oo,oo);
	interval qy=interval(-oo,oo);
	CPoseTrans(qx,qy,d,px, py,alpha);
	CPoseInSegmentsOrCircles(qx, qy, alpha, ax, ay, bx, by,cx,cy,r);
	CPoseTrans(qx,qy,d,px,py,alpha);
}
//----------------------------------------------------------------------
void CPoseTransRotInWallsOrCircles(interval& px, interval& py, interval& alpha, interval& d, interval& psi, vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by, vector<double> cx, vector<double> cy, vector<double> r)
{ 	// All feet should be supported by one of the circles or one of the segments
	interval qx=interval(-oo,oo);
	interval qy=interval(-oo,oo);
	interval beta=interval(-oo,oo);
	CPoseTransRot(qx,qy,beta,d, psi,px, py,alpha);
	CPoseInSegmentsOrCircles(qx, qy, beta, ax, ay, bx, by,cx,cy,r);
	CPoseTransRot(qx,qy,beta,d, psi,px, py,alpha);
}
//----------------------------------------------------------------------
void CPoseRotTransRotInWallsOrCircles(interval& px,interval& py,interval& alpha,interval& phi,interval& d,interval& psi, vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by, vector<double> cx, vector<double> cy, vector<double> r)
{	interval qx=interval(-oo,oo);
	interval qy=interval(-oo,oo);
	interval beta=interval(-oo,oo);
	CPoseRotTransRot(qx,qy,beta,phi,d,psi,px,py,alpha);
	CPoseInSegmentsOrCircles(qx, qy, beta, ax, ay, bx, by,cx,cy,r);
	CPoseRotTransRot(qx,qy,beta,phi,d,psi,px,py,alpha);
}
//----------------------------------------------------------------------
void CPoseRotTransPointInWallsOrCircles(interval& px,interval& py,interval& alpha,interval& phi,interval& d, vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by, vector<double> cx, vector<double> cy, vector<double> r)
{	interval qx=interval(-oo,oo);
	interval qy=interval(-oo,oo);
	interval beta=interval(-oo,oo);
	CPoseRotTrans(qx,qy,beta,phi,d,px,py,alpha);
        CinSegmentsOrCircles(qx, qy, ax, ay, bx, by,cx,cy,r);
	CPoseRotTrans(qx,qy,beta,phi,d,px,py,alpha);
}
//----------------------------------------------------------------------
void CPoseTransPointInLine(interval& px,interval& py,interval& alpha, interval& d, double ax, double ay, double bx, double by)
{	interval qx=interval(-oo,oo);
	interval qy=interval(-oo,oo);
	CPoseTrans(qx,qy,d, px, py,alpha);
	CDistanceDirLine(d,px, py, alpha,ax, ay,bx,by);
        CinLine(qx,qy,ax,ay,bx,by);
	CPoseTrans(qx,qy,d, px, py,alpha);
}
//----------------------------------------------------------------------
void SPoseTransPointInCircles(interval& px,interval& py,interval& alpha, interval& d, vector<double>& cx, vector<double>& cy, vector<double>& r, bool truth)
{	interval qx=interval(-oo,oo);
	interval qy=interval(-oo,oo);
	CPoseTrans(qx,qy,d, px, py,alpha);
        SinCircles(qx,qy,cx,cy,r,truth);
	CPoseTrans(qx,qy,d, px, py,alpha);
}
//----------------------------------------------------------------------
void SPoseTransPointInWall(interval& px,interval& py,interval& alpha, interval& d0, double ax, double ay, double bx, double by, bool truth)
{    
	interval qx=interval(-oo,oo);
	interval qy=interval(-oo,oo);
	interval d=interval(-oo,oo);
	interval px1(px);
	interval py1(py);
	interval alpha1(alpha);
	interval px2(px);
	interval py2(py);
	interval alpha2(alpha);
	CPoseTransPointInLine(px1,py1,alpha1,d,ax,ay,bx,by);
	if (truth) d=Inter(d,d0); else Cnotin(d,d0);
	CPoseTransPointInLine(px1,py1,alpha1,d,ax,ay,bx,by);
        SPoseTowardSegment(px2, py2, alpha2, ax, ay,bx,by,truth);
	if (truth) {px=Inter(px1,px2);py=Inter(py1,py2);alpha=Inter(alpha1,alpha2);}  // Intersection des deux contracteurs
	else {px=Union(px1,px2);py=Union(py1,py2);alpha=Union(alpha1,alpha2);}        // Union des deux contracteurs

}
//----------------------------------------------------------------------
void SPoseTransPointInWalls(interval& px,interval& py,interval& alpha, interval& d0, vector<double>& ax, vector<double>& ay, vector<double>& bx, vector<double>& by, bool truth)
{	vector<interval> Lpx;
	vector<interval> Lpy;
	for(unsigned int i = 0; i < ax.size(); i++){
		interval px1(px); interval py1(py); interval d(d0); interval theta(alpha);
		double eps = 0.25;
		double dx = (bx[i] - ax[i]);
		double dy = (by[i] - ay[i]);
		double dxn = (dx/hypot(dx,dy))*eps;
		double dyn = (dy/hypot(dx,dy))*eps;
                SPoseTransPointInWall(px1,py1,theta,d,ax[i]-dxn,ay[i]-dyn,bx[i]+dxn,by[i]+dyn,truth);
		Lpx.push_back(px1);
		Lpy.push_back(py1);
	}
	if (truth)
	{
		px = Inter(px, Union(Lpx));
		py = Inter(py, Union(Lpy));
	} 
	else 
	{
		px = Inter(px, Inter(Lpx));
		py = Inter(py, Inter(Lpy));
	}
}
//----------------------------------------------------------------------
void SPoseTowardSegment(interval& mx, interval& my, interval& theta, double& ax, double& ay, double& bx,double& by, bool truth)
{	// La pose m=(mx,my,theta) pointe sur le segment [a,b]
	if ((mx.is_empty())||(my.is_empty())||(theta.is_empty()))
	{ 
		mx = interval(); my = interval(); theta = interval(); return;
	}
	interval ma_x=ax-mx;       interval ma_y=ay-my;
	interval mb_x=bx-mx;       interval mb_y=by-my;
	double ab_x=bx-ax;         double ab_y=by-ay;
	interval ux=Cos(theta);    interval uy=Sin(theta);
	interval z1=Det(ma_x,ma_y,ux,uy);
	interval z2=Det(ux,uy,mb_x,mb_y);
	interval z3=Det(ma_x,ma_y,ab_x,ab_y);
	interval z4=Det(ux,uy,ab_x,ab_y);

	interval z5=interval(0,oo);      // si tous les zi >0 alors on satisfait la contrainte
	if (!truth) z5=interval(-oo,0);
        Cmin(z5,z1,z2,z3,z4);
        Cdet(z4,ux,uy,ab_x,ab_y);
        Cdet(z3,ma_x,ma_y,ab_x,ab_y);
        Cdet(z2,ux,uy,mb_x,mb_y);
        Cdet(z1,ma_x,ma_y,ux,uy);
        Csin(uy,theta);           Ccos(ux,theta);
        Csub(mb_y,by,my);       Csub(mb_x,bx,mx);
        Csub(ma_y,ay,my);       Csub(ma_x,ax,mx);
}
//----------------------------------------------------------------------
void Cnocross(interval& px, interval& py, interval& mx, interval& my, double& ax, double& ay, double& bx, double& by)
{	interval ma_x=ax-mx;	interval ma_y=ay-my;
	interval mb_x=bx-mx;	interval mb_y=by-my;
	interval ap_x=px-ax;	interval ap_y=py-ay;
	interval pb_x=bx-px;	interval pb_y=by-py;
	interval pm_x=mx-px;	interval pm_y=my-py;
	double ab_x=bx-ax;	double ab_y=by-ay;
	interval z1=interval(-oo,oo);        Cdet(z1,ab_x,ab_y,ma_x,ma_y);
	interval z2=interval(-oo,oo);        Cdet(z2,ab_x,ab_y,ap_x,ap_y);
        interval z3=z1*z2;
        interval z4=interval(-oo,oo);        Cdet(z4,pm_x,pm_y,ap_x,ap_y);
	interval z5=interval(-oo,oo);        Cdet(z5,pm_x,pm_y,pb_x,pb_y);
	interval z6=z4*z5;
	interval z7(0,oo);
        Cmax(z7,z3,z6);                   Cmul(z6,z4,z5);
        Cdet(z5,pm_x,pm_y,pb_x,pb_y);     Cdet(z4,pm_x,pm_y,ap_x,ap_y);
        Cmul(z3,z1,z2);
        Cdet(z2,ab_x,ab_y,ap_x,ap_y);     Cdet(z1,ab_x,ab_y,ma_x,ma_y);
        Csub(pm_y,my,py);   Csub(pm_x,mx,px);
        Csub(pb_y,by,py);   Csub(pb_x,bx,px);
        Csub(ap_y,py,ay);   Csub(ap_x,px,ax);
        Csub(mb_y,by,my);   Csub(mb_x,bx,mx);
        Csub(ma_y,ay,my);   Csub(ma_x,ax,mx);
	if ((mx.is_empty())||(my.is_empty())||(px.is_empty())||(py.is_empty()))
        {  mx = interval(); my = interval(); px = interval(); py = interval();	}
}
//----------------------------------------------------------------------
void CLegCrossNoSegment(interval& dist, interval& px, interval& py, interval& theta, vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by)
{      	// Aucun segment ne doit être croise
	if ((ax.size() == 0)||(ay.size() == 0)||(bx.size() == 0)||(by.size() == 0)) return;
	interval ux=Cos(theta);
	interval uy=Sin(theta);
	interval dx=ux*dist;
	interval dy=uy*dist;
	interval leg_x=px+dx;
	interval leg_y=py+dy;
	for (unsigned int j = 0; j < ax.size(); j++)
            Cnocross(px,py,leg_x,leg_y, ax[j], ay[j], bx[j], by[j]);
        Cadd(leg_y,py,dy);
        Cadd(leg_x,px,dx);
        Cmul(dy,uy,dist);
        Cmul(dx,ux,dist);
        Csin(uy,theta);
        Ccos(ux,theta);
}
//----------------------------------------------------------------------
void CLegOnWalls(interval& dist, interval& px, interval& py, interval& theta, vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by)
{      	// Toutes les pattes doivent être sur le mur
	if ((ax.size() == 0)||(ay.size() == 0)||(bx.size() == 0)||(by.size() == 0)) return;
	interval ux = Cos(theta);
	interval uy = Sin(theta);
	interval dx = ux*dist;
	interval dy = uy*dist;
	interval leg_x = px + dx;
	interval leg_y = py + dy;
        CinSegments(leg_x, leg_y, ax, ay, bx, by);
        Cadd(leg_y, py, dy);
        Cadd(leg_x, px, dx);
        Cmul(dy, uy, dist);
        Cmul(dx, ux, dist);
        Csin(uy, theta);
        Ccos(ux, theta);
}
//----------------------------------------------------------------------
void CLegOnWallsOrCircles(interval& dist, interval& px, interval& py, interval& theta, vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by, vector<double> cx, vector<double> cy, vector<double> r)
{      	// Toutes les pattes doivent être sur le mur ou sur un des cercles
	interval ux = Cos(theta);
	interval uy = Sin(theta);
	interval dx = ux*dist;
	interval dy = uy*dist;
	interval leg_x = px + dx;
	interval leg_y = py + dy;
        CinSegmentsOrCircles(leg_x, leg_y, ax, ay, bx, by, cx, cy, r);
        Cadd(leg_y, py, dy);
        Cadd(leg_x, px, dx);
        Cmul(dy, uy, dist);
        Cmul(dx, ux, dist);
        Csin(uy, theta);
        Ccos(ux, theta);
}
//----------------------------------------------------------------------
void Cnotin(interval& X, const interval& Y)
{ 	X=(X&interval(-oo,Y.lb()))|(X&interval(Y.ub(),oo));
}
//----------------------------------------------------------------------
void C_q_in(interval& x, int q, vector<interval>& y)
{	Array<IntervalVector> V((int)y.size());
	for (unsigned int i = 0 ;i < y.size(); i++)
	{
		IntervalVector* Vi = new IntervalVector(1, y[i]);
		V.set_ref(i, *Vi);
	}
	x = (qinter(V, q))[0];
	for (unsigned int i = 0 ;i < y.size(); i++)	
		delete &(V[i]);	
}
//----------------------------------------------------------------------
void SinRing(interval& X, interval& Y, double cx, double cy, interval R, bool outer)
{ 	if (outer==true)
        { CinRing(X,Y,cx,cy,R); return;	}
	interval Xa(X),Ya(Y),Xb(X),Yb(Y);
        CinRing(Xa,Ya,cx,cy,interval(-1,R.lb()));
	CinRing(Xb,Yb,cx,cy,interval(R.ub(),POS_INFINITY));
	X=Xa|Xb;
	Y=Ya|Yb;
}
//----------------------------------------------------------------------
