// Compatiblity layer with ibex for the simple interval library from Luc JAULIN, with minor modifications from Fabrice LE BARS and Jeremy NICOLA.

#ifndef __BOX__
#define __BOX__

#include "vibes.h";
#include "interval.h"

class box
{
public:
	interval* data;
	int dim;

public:
	//----------------------------------------------------------------------
	// Constructors/destructors
	//----------------------------------------------------------------------
	box();
	box(int);
	box(interval x);
	box(interval x, interval y);
	box(interval x, interval y, interval z);
	box(interval x, int n);
	box(const box&);
	~box();
	//----------------------------------------------------------------------
	// Operators
	//----------------------------------------------------------------------
	box& operator=(const box&);
	friend box operator+(const box&, const box&);
	friend box operator-(const box&);
	friend box operator-(const box&, const box&);
	friend box operator*(const interval&, const box&);
	friend box operator*(const box&, const interval&);
	friend box operator*(const double, const box&);
	friend box operator*(const box&, const double);
	friend box operator&(const box&, const box&);
	friend box operator|(const box&, const box&);
	friend bool operator==(const box&, const box&);
	interval& operator[](int) const;
	friend std::ostream& operator<<(std::ostream& os, const box& X);
#ifdef QT_VERSION 
	friend QDebug operator<<(QDebug os, const box&X)
	{   
		os.nospace() << "box :" << "\t dim=" << X.dim << "\n";
		if (X.IsEmpty()) os.nospace() << "EmptyBox";
		for (int i = 1; i <= X.dim; i++)
			os.nospace() << "  " << i << ": "<< X[i] << "\n";
		return (os.space());
	}
#endif // QT_VERSION 
	//----------------------------------------------------------------------
	// Member functions
	//----------------------------------------------------------------------
	box& Intersect(const box& Y);
	double Width(void);
	double SumWidth(void);
	bool IsEmpty(void) const;
    void Resize(int);           //RqM peut être un pbs
        void DrawBox(const char* param);


};

//----------------------------------------------------------------------
// Box-valued functions
//----------------------------------------------------------------------
box Inf(const box&);
box Sup(const box&);
box Center(const box&);
box Center(const box&, std::vector<int>&);
box Zeros(int);
box Empty(int);
box EmptyBox(int);
box EmptyBox(const box&);
box Rn(int);
box Rotate(const box&X, double t, int a);  // transformation
box Concat(const box&, const box&);
box Proj(const box&, int, int);
box Inter(const box&, const box&);
box Inter(std::vector<box>&);
box Union(const box&, const box&);
box Union(std::vector<box>&);
box Inflate(const box&, double);

//----------------------------------------------------------------------
// Other functions
//----------------------------------------------------------------------
int Size(const box&);
double Width(const box&);
double Volume(const box&);
void Update(const box&);
interval Norm(const box&);
interval Scal(const box&, const box&);
interval Det(const box&, const box&);
bool Disjoint(const box&, const box&);
bool Subset(const box&, const box&);
int AxePrincipal(const box&);
double decrease(const box&, const box&);
void DrawArrow(double x1,double y1,double dx,double dy, double r, const char* param);
void DrawPolygon(double x,double y,double theta,vector<double> X, vector<double> Y,const char* param);
void DrawRobot(double x,double y,double theta, double s,const char* param);
void DrawSquare(double x,double y,double theta, double s,const char* param);



//----------------------------------------------------------------------
// Contractors
//----------------------------------------------------------------------
void Cadd(box& Z, box& X, box& Y);
void Csub(box& Z, box& X, box& Y);
#define CProd Cmul
void Cmul(box& Y, interval& a, box& X);
void Cnorm(interval& R, box& X);
#define Cdistance Cdist
void Cdist(interval& R, box& X, box& Y);
#define CProdScalaire Cscal
void Cscal(interval& R, box& X, box& Y);
#define COrtho Cortho
void Cortho(box& X, box& Y);
void SinBox(box& X, const box& Y, bool outer);
void S_q_in(box& x, int q, const std::vector<box>& yj, bool outer);
//----------------------------------------------------------------------
// Other
//----------------------------------------------------------------------
void Bisect(box& X, box& X1, box& X2);
void BisectAlong(box& X, box& X1, box& X2, int i);

#endif // __BOX__
