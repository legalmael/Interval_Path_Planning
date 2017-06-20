// Compatiblity layer with ibex for the simple interval library from Luc JAULIN, with minor modifications from Fabrice LE BARS and Jeremy NICOLA.

#include "box.h"

using namespace std;
using namespace ibex;

//----------------------------------------------------------------------
// Constructors/destructors
//----------------------------------------------------------------------
box::box()
{
	dim = 0;
	data = new interval[1];
}
//----------------------------------------------------------------------
box::box(int a)
{
	dim = a;
	data = new interval[dim];
}
//----------------------------------------------------------------------
box::box(interval x)
{
	dim = 1;
	data = new interval[dim];
	(*this)[1] = x;
}
//----------------------------------------------------------------------
box::box(interval x, interval y)
{
	dim = 2;
	data = new interval[dim];
	(*this)[1] = x;
	(*this)[2] = y;
}
//----------------------------------------------------------------------
box::box(interval x, interval y, interval z)
{
	dim = 3;
	data = new interval[dim];
	(*this)[1] = x;
	(*this)[2] = y;
	(*this)[3] = z;
}
//----------------------------------------------------------------------
box::box(interval x, int n)
{
	dim = n;
	data = new interval[dim];
	for (int k = 1; k <= dim; k++) (*this)[k] = x;
}
//----------------------------------------------------------------------
box::box(const box& X)
{
	if (&X == this) return;
	dim = Size(X);
	data = new interval[dim];
	for (int k = 1; k <= dim; k++) (*this)[k] = X[k];
}
//----------------------------------------------------------------------
box::~box() { delete[] data; }
//----------------------------------------------------------------------
// Operators
//----------------------------------------------------------------------
box& box::operator=(const box& X)
{
	delete[] data;
	dim = Size(X);
	data = new interval[dim];
	for (int k = 1; k <= dim; k++) (*this)[k] = X[k];
	return *this;
}
//----------------------------------------------------------------------
box operator+(const box& X, const box& Y)
{
	box Ans(Size(X));
	for (int k = 1; k <= Size(Ans); k++) Ans[k] = X[k] + Y[k];
	return Ans;
}
//----------------------------------------------------------------------
box operator-(const box& X)
{
	box Ans(Size(X));
	for (int k = 1; k <= Size(Ans); k++) Ans[k] = -X[k];
	return Ans;
}
//----------------------------------------------------------------------
box operator-(const box& X, const box& Y)
{
	box Ans(Size(X));
	for (int k = 1; k <= Size(Ans); k++) Ans[k] = X[k] - Y[k];
	return Ans;
}
//----------------------------------------------------------------------
box operator*(const interval& a, const box& X)
{
	box Ans(Size(X));
	for (int k = 1; k <= Size(Ans); k++) Ans[k] = a*X[k];
	return Ans;
}
//----------------------------------------------------------------------
box operator*(const box& X, const interval& a)
{
	box Ans(Size(X));
	for (int k = 1; k <= Size(Ans); k++) Ans[k] = a*X[k];
	return Ans;
}
//----------------------------------------------------------------------
box operator*(const double a, const box& X)
{
	box Ans(Size(X));
	for (int k = 1; k <= Size(Ans); k++) Ans[k] = a*X[k];
	return Ans;
}
//----------------------------------------------------------------------
box operator*(const box& X, const double a)
{
	box Ans(Size(X));
	for (int k = 1; k <= Size(Ans); k++) Ans[k] = a*X[k];
	return Ans;
}
//----------------------------------------------------------------------
box operator&(const box& X, const box& Y)
{
	return Inter(X, Y);
}
//----------------------------------------------------------------------
box operator|(const box& X, const box& Y)
{  
	return Union(X, Y);
}
//----------------------------------------------------------------------
bool operator==(const box& A, const box& B)
{
	if (A.dim != B.dim) return false;
	for (int k = 1; k <= A.dim; k++)
	{
		if (!(A[k] == B[k]))
			return false;
	}
	return true;
}
//----------------------------------------------------------------------
interval& box::operator[](int i) const
{
	return data[i - 1];
}
//----------------------------------------------------------------------
std::ostream& operator<<(std::ostream& os, const box& X)
{
	cout << "box :" << "\t dim=" << Size(X) << "\n";
	if (X.IsEmpty()) os << "EmptyBox";
	for (int i = 1; i <= Size(X); i++)
		os << "  " << i << ": " << X[i] << "\n";
	return os;
}
//----------------------------------------------------------------------
// Member functions
//----------------------------------------------------------------------
box& box::Intersect(const box& Y)
{
	box X = *this;
	box Z = Inter(X, Y);
	*this = Z;
	return *this;
}
//----------------------------------------------------------------------
double box::Width(void)
{
	box X = *this;
	int i = AxePrincipal(X);
	interval Xi = X[i];
	double w = Xi.ub() - Xi.lb();
	return (w);
}
//----------------------------------------------------------------------
double box::SumWidth(void)
{  
	box X = *this;
	if (X.IsEmpty()) return -oo;
	double w = X[1].ub()-X[1].lb();
	for (int i = 2; i <= Size(X); i++)
		w = w+X[i].ub()-X[i].lb();
	return w;
}
//----------------------------------------------------------------------
bool box::IsEmpty(void) const
{
	if (dim == 0) return true;
	for (int i = 1; i <= dim; i++)
	{
		if ((*this)[i].is_empty()) return true;
	}
	return false;
}
//----------------------------------------------------------------------
void box::Resize(int dim1)
{
	box X(dim1);
	for (int k = 1; k <= min(dim, dim1); k++)
	{
		X[k] = (*this)[k];
	}
	for (int k = dim+1; k <= dim1; k++)
	{
		X[k] = interval(-oo,oo);
	}
	delete[] data; (*this) = X;
}
//----------------------------------------------------------------------
void box::DrawBox(const char* param)
{vibes::drawBox((*this)[1].inf,(*this)[1].sup,(*this)[2].inf,(*this)[2].sup,param);
}
//----------------------------------------------------------------------


void DrawArrow(double x1,double y1,double dx,double dy, double r, const char* param)
{       vector<double> X,Y;
        double x2=x1+dx;
        double y2=y1+dy;
        X.push_back(x1);    Y.push_back(y1);
        X.push_back(x2); Y.push_back(y2);
        double a=3.14-(3.14/4);
        double px=x2+r*(cos(a)*dx-sin(a)*dy);    // cot� de pointe
        double py=y2+r*(sin(a)*dx+cos(a)*dy);
        double qx=x2+r*(cos(-a)*dx-sin(-a)*dy);  // autre cot� de pointe
        double qy=y2+r*(sin(-a)*dx+cos(-a)*dy);
        X.push_back(px); Y.push_back(py);
        X.push_back(x2); Y.push_back(y2);
        X.push_back(qx); Y.push_back(qy);
        vibes::drawLine(X,Y,param);
 }

void DrawPolygon(double x,double y,double theta,vector<double> X, vector<double> Y,const char* param)
{       vector<double> X1,Y1;
        for (int k=0;k<X.size();k++)
        {   double x1=x+cos(theta)*X[k]-sin(theta)*Y[k];
            double y1=y+sin(theta)*X[k]+cos(theta)*Y[k];
            X1.push_back(x1); Y1.push_back(y1);
        }
        vibes::drawLine(X1,Y1,param);
}
//--------------------------------------------------------------------------------------------------
void DrawRobot(double x,double y,double theta, double s,const char* param)
{   vector<double> X,Y;
    X.push_back(0*s); Y.push_back(-1*s);
    X.push_back(3*s); Y.push_back(0*s);
    X.push_back(0*s); Y.push_back(1*s);
    X.push_back(0*s); Y.push_back(-1*s);
    DrawPolygon(x,y,theta,X,Y,param);
}
void DrawSquare(double x,double y,double theta, double s,const char* param)
{   vector<double> X,Y;
    X.push_back(-s); Y.push_back(s);
    X.push_back(s); Y.push_back(s);
    X.push_back(s); Y.push_back(-s);
    X.push_back(-s); Y.push_back(-s);
    X.push_back(-s); Y.push_back(s);
    DrawPolygon(x,y,theta,X,Y,param);
}



//----------------------------------------------------------------------
// Box-valued functions
//----------------------------------------------------------------------
box Inf(const box& X)
{
	box Ans(Size(X));
	for (int k = 1; k <= Size(Ans); k++) Ans[k] = X[k].lb();
	return Ans;
}
//----------------------------------------------------------------------
box Sup(const box& X)
{
	box Ans(Size(X));
	for (int k = 1; k <= Size(Ans); k++) Ans[k] = X[k].ub();
	return Ans;
}
//----------------------------------------------------------------------
box Center(const box& X)
{
	int sizeX = Size(X); box Ans(sizeX);
	//if (X.IsEmpty()) Ans = Empty(sizeX);
	if (X.IsEmpty()) Ans = EmptyBox(X);
	else { for (int k = 1; k <= sizeX; k++) Ans[k] = Center(X[k]); }
	return Ans;
}
//-----------------------------------------------------------------------
box Center(const box& X, vector<int>& v)
{
	int sizev = (int)(v.size() - 1); box Ans = X;
	for (int k = 1; k <= sizev; k++) Ans[v[k]] = Center(X[v[k]]);
	return Ans;
}
//-----------------------------------------------------------------------
box Zeros(int d)
{
	box Ans(d);
	for (int k = 1; k <= d; k++) Ans[k] = 0;
	return Ans;
}
//----------------------------------------------------------------------
box EmptyBox(int n)
    {return (box(n)); }
//----------------------------------------------------------------------
box EmptyBox(const box& X)
    { return (box(Size(X)));}
//----------------------------------------------------------------------
box Rn(int a)
{	box Ans(a);
	for (int k = 1; k <= a; k++) Ans[k] = interval(-oo, oo);
	return Ans;
}
//----------------------------------------------------------------
box Rotate(const box&X, double t, int a)
{return box(cos(t)*X[1]-a*sin(t)*X[2],a*sin(t)*X[1]+cos(t)*X[2]);}

//----------------------------------------------------------------------
// Produit Cartesien ou Concatenation de deuX paves X et y :
// Ans=[X,Y]     =>     Ans=Concat(X,Y); 
box Concat(const box& X, const box& Y)
{
	double dim = X.dim + Y.dim; box Ans((int)dim);
	if ((!X.IsEmpty()) && (!Y.IsEmpty()))
	{
		for (int i = 1; i <= dim; i++)
		{
			if (i <= Size(X)) Ans[i] = X[i]; else Ans[i] = Y[i - X.dim];
		}
	}
	return Ans;
}
//----------------------------------------------------------------------
// Projection du pave X dans un espace de dimension dim=(j-i)+1;
// X=[[X1],[X2],..,[Xi],..,[Xj],..[Xn]]
// =>  Proj(X,i,j)=[[Xi],..,[Xj]] et Proj(X,i,i)=[Xi]
box Proj(const box& X, int i, int j)
{
	int dim = abs(j - i) + 1; box Ans(dim);
	if (!X.IsEmpty())
	{
		int lb = min(i, j);
		for (int k = 1; k <= dim; k++) Ans[k] = X[k + lb - 1];
	}
	return Ans;
}
//----------------------------------------------------------------------
box Inter(const box& X, const box& Y)
{	// Should update with ibex IntervalVector...
	box Ans(Size(X));
	if ((X.IsEmpty()) || (Y.IsEmpty())) { return Ans; }
	for (int k = 1; k <= Size(Ans); k++)
        {       Ans[k] = Inter(X[k], Y[k]);
		if (Ans[k].is_empty()) { Update(Ans); return Ans; }
	}
	return Ans;
}
//----------------------------------------------------------------------
box Inter(vector<box>& x)
{
        box E = EmptyBox(0);
	if (x.size() == 0) return E;
	box r = x[0];
	for (unsigned int i = 1; i < x.size(); i++)
		r = Inter(r, x[i]);
	return r;
}
//----------------------------------------------------------------------
box Union(vector<box>& x)
{
	box E = EmptyBox(0);
	if (x.size() == 0) return E;
	box r = x[0];
	for (unsigned int i = 1; i < x.size(); i++)
		r = Union(r, x[i]);
	return r;
}
//----------------------------------------------------------------------
box Union(const box& X, const box& Y)
{
	box Ans(max(Size(X), Size(Y)));
	if (X.IsEmpty()) return (Y);
	if (Y.IsEmpty()) return (X);
	for (int k = 1; k <= Size(Ans); k++)  Ans[k] = Union(X[k], Y[k]);
	return Ans;
}
//----------------------------------------------------------------------
box Inflate(const box& X, double eps)
{
	Update(X);
	box Ans(Size(X));
	for (int k = 1; k <= Size(Ans); k++)  { Ans[k] = Inflate(X[k], eps); }
	return Ans;
}
//----------------------------------------------------------------------
// Other functions
//----------------------------------------------------------------------
int Size(const box& X) { return (X.dim); }
//----------------------------------------------------------------------
double Width(const box& X) 
{
	return Width(X[AxePrincipal(X)]);
}
//----------------------------------------------------------------------
double Volume(const box& X)
{
	double vol = 1;
	for (int i = 1; i <= Size(X); i++) vol = vol*Width(X[i]);
	return vol;
}
//----------------------------------------------------------------------
void Update(const box& X)
{
	for (int i = 1; i <= Size(X); i++)
	{
		if (X[i].is_empty())
		{
			for (int j = 1; j <= Size(X); j++)
				X[j] = interval(); return;
		}
	}
}
//----------------------------------------------------------------------
interval Norm(const box& X)
{
	if (X.IsEmpty()) return interval();
	interval r = 0;
	for (int i = 1; i <= Size(X); i++) r = r + Sqr(X[i]);
	return (Sqrt(r));
}
//----------------------------------------------------------------------
interval NormEuclid(const box& X, const box& Y)
{
	if (Size(X) != Size(Y)) return interval();
	if (X.IsEmpty()||Y.IsEmpty()) return interval();
	interval r = 0;
	for (int i = 1; i <= Size(X); i++) r = r + Sqr(Y[i] - X[i]);
	return (Sqrt(r));
}
//----------------------------------------------------------------------
interval NormInf(const box& X, const box& Y)
{
	if (Size(X) != Size(Y)) return interval();
	if (X.IsEmpty()||Y.IsEmpty()) return interval();
	interval ans = Abs(Y[1] - X[1]);
	for (int i = 1; i <= Size(X); i++) ans = Max(ans, Abs(Y[i] - X[i]));
	return ans;
}
//----------------------------------------------------------------------
interval Scal(const box& U, const box& V)
{
	interval sum = 0;
	for (int i = 1; i <= Size(U); i++)  sum = sum + U[i] * V[i];
	return (sum);
}
//----------------------------------------------------------------------
interval Det(const box& U, const box& V)
{
	interval u1 = U[1];
	interval v2 = V[2];
	interval v1 = V[1];
	interval u2 = U[2];
	interval r = u1*v2 - v1*u2;
	return u1*v2 - v1*u2;
}
//----------------------------------------------------------------------
bool Disjoint(const box& X, const box& Y)
{
	if (X.IsEmpty() || Y.IsEmpty()) return true;
	for (int i = 1; i <= Size(X); i++)
		if (Disjoint(X[i], Y[i])) return true;
	return false;
}
//----------------------------------------------------------------------
bool Subset(const box& X, const box& Y)
{
	if (X.IsEmpty()) return true;
	if (Y.IsEmpty()) return false;
	bool b = true;
	for (int k = 1; k <= Size(X); k++) b = b && Subset(X[k], Y[k]);
	return (b);
}
//----------------------------------------------------------------------
bool SubsetStrict(const box& X, const box& Y)
{
	if (Y.IsEmpty()) return false;
	if (X.IsEmpty()) return true;
	bool b = true;
	for (int k = 1; k <= Size(X); k++) b = b && SubsetStrict(X[k], Y[k]);
	return (b);
}
//----------------------------------------------------------------------
int AxePrincipal(const box& X)
{
	int kmax = 1;
	double widthmax = Width(X[kmax]);
	for (int k = 2; k <= Size(X); k++)
	{
		if (Width(X[k]) > widthmax)
		{
			kmax = k; widthmax = Width(X[k]);
		}
	}
	return kmax;
}
//----------------------------------------------------------------------
double decrease(const box& X, const box& Y)
{
	//if (X.IsEmpty()) return (-oo);
	double e = 0;
	for (int k = 1; k <= X.dim; k++)
	{
		if ((X[k].is_empty())||(Y[k].is_empty())) return (-1);
		double e1 = 0;
		double Xinf = X[k].lb(), Xsup = X[k].ub();
		double Yinf = Y[k].lb(), Ysup = Y[k].ub();
		if (Xsup >= Ysup) e1 = max(e1, Xsup - Ysup);
		if (Xinf <= Yinf) e1 = max(e1, Yinf - Xinf);
		if (Xsup < Ysup) e1 = max(e1, Ysup - Xsup);
		if (Xinf > Yinf) e1 = max(e1, Xinf - Yinf);
		e = max(e, e1);
	}
	return e;
}
//----------------------------------------------------------------------
// Contractors
//----------------------------------------------------------------------
void Cadd(box& Z, box& X, box& Y)
{
	for (int k = 1; k <= Size(X); k++)
                Cadd(Z[k], X[k], Y[k]);
}
//----------------------------------------------------------------------
void Csub(box& Z, box& X, box& Y)
{
	for (int k = 1; k <= Size(X); k++)
                Csub(Z[k], X[k], Y[k]);
}
//----------------------------------------------------------------------
void Cmul(box& Y, interval& a, box& X)
{
	for (int k = 1; k <= Size(X); k++)
                Cmul(Y[k], a, X[k]);
}
//----------------------------------------------------------------------
void Cnorm(interval& R, box& X)
{
	int imax = X.dim;
	box Sum2(imax);
	box X2(imax);
	for (int i = 1; i <= imax; i++)
	{
		X2[i] = interval(-oo, oo);
                Csqr(X2[i], X[i]);
	}
	Sum2[1] = X2[1];
	for (int i = 2; i <= imax; i++)
	{
		Sum2[i] = interval(-oo, oo);
                Cadd(Sum2[i], Sum2[i - 1], X2[i]);
	}
        Csqr(Sum2[imax], R);
        Csqr(Sum2[imax], R);
	for (int i = imax; i >= 2; i--)
	{
                Cadd(Sum2[i], Sum2[i - 1], X2[i]);
	}
	X2[1] = Sum2[1];
	for (int i = imax; i >= 1; i--)
                Csqr(X2[i], X[i]);
}
//----------------------------------------------------------------------
void Cdist(interval& R, box& X, box& Y)
{	box Z = Y - X;
	Cnorm(R, Z);
        Csub(Z, Y, X);
}
//----------------------------------------------------------------------
void Cscal(interval& R, box& X, box& Y)
{	int imax = X.dim;
	box SumXiYi(imax);
	box XiYi(imax);
	for (int i = 1; i <= imax; i++)
	{
		XiYi[i] = X[i] * Y[i];
	}
	SumXiYi[1] = XiYi[1];
	for (int i = 2; i <= imax; i++)
	{
		SumXiYi[i] = interval(-oo, oo);
                Cadd(SumXiYi[i], SumXiYi[i - 1], XiYi[i]);
	}
	R = Inter(R, SumXiYi[imax]);
	SumXiYi[imax] = R;
	for (int i = imax; i >= 2; i--)
	{
                Cadd(SumXiYi[i], SumXiYi[i - 1], XiYi[i]);
	}
	XiYi[1] = SumXiYi[1];
	for (int i = imax; i >= 1; i--)
	{
                Cmul(XiYi[i], X[i], Y[i]);
	}
}
//----------------------------------------------------------------------
void Cortho(box& X, box& Y)
{	interval S(0, 0);
	Cscal(S, X, Y);
}
//----------------------------------------------------------------------
void SinBox(box& X, const box& Y, bool outer)  // constraint : x belongs to the box Y
    {   if (X.IsEmpty()) {return;}
        if (outer)  {X=X&Y; return; };
        if (Y.IsEmpty()) return;
        box U=EmptyBox(X);
        for (int j=1;j<X.dim+1;j++)
        {   box A(X);
            Cnotin(A[j],Y[j]);
            U=U|A;
        }
        X=U;
    }
//----------------------------------------------------------------------
void S_q_in(box& x, int q, const vector<box>& yj, bool outer)
{   // q : number of constraints to be satisfied
    int m=(int)yj.size();

    if (!outer)  { S_q_in(x,m-q+1    ,yj,1);    return;  }

    Array<IntervalVector> V(m);
        for (unsigned int j = 0 ; j < m; j++)
        {       IntervalVector* Vj = new IntervalVector(yj[j].dim);
                for (int i = 1; i <= yj[j].dim; i++) (*Vj)[i-1] = yj[j][i]&x[i];
                V.set_ref(j, *Vj);
	}
	IntervalVector res = qinter(V, q);
        for (int i = 1; i <= x.dim; i++)
                x[i] = res[i-1];
        for (unsigned int j = 0 ;j < m; j++)
                delete &(V[j]);

}
//----------------------------------------------------------------------
void Bisect(box& X, box& X1, box& X2)
{
	BisectAlong(X, X1, X2, AxePrincipal(X));
}
//----------------------------------------------------------------------
void BisectAlong(box& X, box& X1, box& X2, int i)
{
	X1 = X2 = X;
	double m = ((X[i].lb())+(1.01*X[i].ub()))/2.01;
	X1[i] = interval(X1[i].lb(),m);
	X2[i] = interval(m,X2[i].ub());
}

