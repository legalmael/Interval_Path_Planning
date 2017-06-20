#include <iostream>
#include "vibes.h"
#include "interval.h"
#include "box.h"
#include "node.h"
#include <list>

#include <time.h>
#include <chrono>
#include <thread>

#define PI 3.14159265


using namespace std;


// Global Variables

vector<double> sx,sy,ax,ay,bx,by;
int imax; // number of segments for the object
int jmax;  // number of segments

list<Node*> paving;
//--------------------------------------------------------

void initGlobalVariables(){

    sx.push_back(0);        sy.push_back(0);
    sx.push_back(0);        sy.push_back(14);
    sx.push_back(14);       sy.push_back(14);
    sx.push_back(14);       sy.push_back(6);
    sx.push_back(10);       sy.push_back(6);
    sx.push_back(10);       sy.push_back(8);
    sx.push_back(12);       sy.push_back(8);
    sx.push_back(12);       sy.push_back(12);
    sx.push_back(2);        sy.push_back(12);
    sx.push_back(2);        sy.push_back(2);
    sx.push_back(18);       sy.push_back(2);
    sx.push_back(18);       sy.push_back(18);
    sx.push_back(20);       sy.push_back(18);
    sx.push_back(20);       sy.push_back(0);
    sx.push_back(0);        sy.push_back(0);

    ax.push_back(8);        ay.push_back(10);
    ax.push_back(25);       ay.push_back(10);
    bx.push_back(11);       by.push_back(10);
    bx.push_back(28);       by.push_back(10);

    imax=sx.size()-1;
    jmax=ax.size();
}

//---------------------------------------------------------------------
double ToReel (const interval& a){
    return a.ub() - a.lb();
}

// Calcule l'angle entre 2 vecteur de dimension 2
// Attention, il faut des vecteurs et non des box et les vecteurs doivent etre
// de dimension 2
double Angle (const box& V, const box& W){
    if ((Size(V)!=2)||(Size(W)!=2)) cout<<"error";
    double n2,costeta,sinteta;
    interval nv=Norm(V);
    interval nw=Norm(W);
    interval nvw=nv*nw;
    n2=ToReel(nvw);
    costeta=ToReel(Scal(V,W)/n2);
    sinteta=ToReel(Det(V,W))/n2;
    if (sinteta>0) return(acos(costeta));
    else return(-acos(costeta));
}

//Test si le pave A (de dimension 2) touche ou non les obstacles
iboolean InRobot(box& A){
   box s0(2);box s1(2);
   int imax=sx.size()-1;
   for (int i=0;i<=imax-1;i++)

   {   s0[1]=sx[i];   s0[2]=sy[i]; s1[1]=sx[i+1]; s1[2]=sy[i+1];
       box s1_0=s0-s1;
       box s0_A=s0-A;
      // interval d=Det(s0-s1,s0-A);
       interval d=Det(s1_0,s0_A);
       if ((In(0,d)) && (!Disjoint(A,Union(s0,s1)))) return iperhaps;
   }
   box a=Center(A);
   double total=0;
   double alpha;
   for (int i=0; i<=imax-1; i++)
   {  s0[1]=sx[i];   s0[2]=sy[i]; s1[1]=sx[i+1]; s1[2]=sy[i+1];
       box s0_a=box(2);
       s0_a= s0-a;
       box s1_a=box(2);
       s1_a=s1-a;
      alpha=Angle(s0-a,s1-a);
      // alpha=Angle(s0_a,s1_a);
      total=total + alpha;
   }
   if ((total>1)||(total<-1))  return itrue;
   return ifalse;
}

// Test si le pave P (de dimension 3) est ou n'est pas dans l'ensemble des configs admissibles
iboolean Inside(box &P)
{   interval P1=P[1]*0.1;    interval P2=P[2]*PI/180; // convertion in rad
    vector<interval> Ax,Ay,Bx,By;
    double p2c=Center(P2);
    double p1c=Center(P1);
    interval CosP2=Cos(P2);
    double cosp2=cos(p2c);
    interval SinP2=Sin(P2);
    double sinp2=sin(p2c);
    iboolean result=itrue;
    for (int j=0;j<=jmax-1;j++) //pour chaque primitive
    {   //on calcule les coordonn�es de A1, B1, A2 et B2 (obstacles) dans le repere de l'objet,
        box A=box(2);
        box B=box(2);
        A[1]= CosP2*(ax[j]-P1) + SinP2*(ay[j]);
        double fc=cosp2*(ax[j]-p1c) + sinp2*(ay[j]);      //Les 5 lignes qui suivent diminue le pessimisme
        interval df1=-CosP2;                            //pour A[1] en utilisant la forme centr�e
        interval df2= -(ax[j]-P1)*SinP2+CosP2*(ay[j]);  // Idem pour A[2], B[1] et B[2]
        interval A1c=fc+df1*(P1-p1c)+df2*(P2-p2c);
        A[1]=Inter(A[1],A1c);
        //
        A[2]=-SinP2*(ax[j]-P1) + CosP2*(ay[j]);
        fc=-sinp2*(ax[j]-p1c) + cosp2*(ay[j]);
        df1=SinP2;
        df2= -(ax[j]-P1)*CosP2-SinP2*(ay[j]);
        interval A2c=fc+df1*(P1-p1c)+df2*(P2-p2c);
        A[2]=Inter(A[2],A2c);
        //
        B[1]= CosP2*(bx[j]-P1) + SinP2*(by[j]);
        fc=cosp2*(bx[j]-p1c) + sinp2*(by[j]);
        df1=-CosP2;
        df2= -(bx[j]-P1)*SinP2+CosP2*(by[j]);
        interval B1c=fc+df1*(P1-p1c)+df2*(P2-p2c);
        B[1]=Inter(B[1],B1c);
        //
        B[2]=-SinP2*(bx[j]-P1) + CosP2*(by[j]);
        fc=-sinp2*(bx[j]-p1c) + cosp2*(by[j]);
        df1=SinP2;
        df2= -(bx[j]-P1)*CosP2-SinP2*(by[j]);
        interval B2c=fc+df1*(P1-p1c)+df2*(P2-p2c);
        B[2]=Inter(B[2],B2c);

        if (InRobot(A)==true) return ifalse;
        if (InRobot(B)==true) return ifalse;

        for (int i=0;i<=imax-1;i++)
        {   //a[j]b[j] coupe le segment du robot, on rejete la configuration
            box s0 =box(2);    s0[1]=sx[i];   s0[2]=sy[i];
            box s1=box(2);     s1[1]=sx[i+1]; s1[2]=sy[i+1];
            box ds=box(2); ds=s1-s0;            box B_s1=box(2);B_s1=B-s1;
            interval d1=Det(B-s1,A-s1)*Det(B-s0,A-s0);
            interval d2=Det(ds,s1-A)*Det(ds,s1-B);
            if ((Sup(d1)<=0)&&(Sup(d2)<=0)) return ifalse; // FullCoupure
            bool fullDisjoint=Disjoint(Union(A,B),Union(s0,s1));
            if ((Inf(d1)<=0)&&(Inf(d2)<=0)&&(!fullDisjoint)) result=iperhaps;
        }
    }
    return result;
}

//----------------------------------------------------------------
bool ShortestPath(bool outer, Node* pva, Node* pvb, list<Node*>& path){


    if ((outer==false) && (Inside(pva->space)!=itrue || Inside(pvb->space)!=itrue)){
        return false;
    }
    if(pva==pvb){
        path.push_back(pva);
        return true;
    }

    // Init
    for(list<Node*>::iterator it = paving.begin(); it != paving.end(); it++){
        (*it)->distance = Infinity;
    }
    pva->distance = 0;

    list<Node*> L(0);
    L.push_back(pva);
    Node* pV(0);
    while((pvb->distance == Infinity)){
        if (L.empty()){
            //cout << "failed to find a path" << endl;
            return false;
        }
        pV = L.front();    L.pop_front();
        Node* pN(0);
        for(int it = 0; it != pV->Neighbors.size(); it++){
            pN = pV->Neighbors.at(it);
            if (outer){
                if (((pN->distance) == Infinity)&&(pN->category != ifalse)){
                    pN->distance = pV->distance+1;
                    L.push_back(pN);
                }
            }
            else{
                if (((pN->distance) == Infinity)&&(pN->category == itrue)){
                    pN->distance = pV->distance+1;
                    L.push_back(pN);
                }
            }
        }
    }
    path.push_front(pvb);
        Node* pvl = pvb;
        Node* pN(0);
        for(int i=pvb->distance-1; i>0; i--){
            for(int it = 0; it != pvl->Neighbors.size(); it++){
                pN = pvl->Neighbors.at(it);
                if(pN->distance == i){
                    path.push_front(pN);
                    pvl = pN;
                    break;
                }
            }
        }
        path.push_front(pva);
    return true;
}

void drawPath(list<Node*>& path){

    vector<double> vectx;   vector<double> vecty;

    Node* Vx = path.front();   path.pop_front();
    box boxa = Vx->space;
    double vax = Center(boxa[1]);    double vay = Center(boxa[2]);
    vectx.push_back(vax);    vecty.push_back(vay);

    box boxb;
    for(list<Node*>::iterator it = path.begin(); it != path.end(); it++){
        boxb = (*(it))->space;
        double vbx = Center(boxb[1]);    double vby = Center(boxb[2]);
        vectx.push_back(vbx);    vecty.push_back(vby);

        vibes::drawLine(vectx,vecty,"black");

        vectx.clear();  vecty.clear();
        vectx.push_back(vbx);    vecty.push_back(vby);
    }
    path.push_front(Vx);
}

//----------------------------------------------------------------
list<Node*> FeasiblePath1 (box pointA, box pointB, box startSpace, double epsilon){

    time_t t = clock();

    Node* pva(0);    Node* pvb(0);
    Node* pX = new Node();
    pX->space = startSpace;

    if (!(Subset(pointA,pX->space)) || !(Subset(pointB,pX->space))){
        cout << "Error: pointA and pointB should belong to P0" << endl;
        return list<Node*>();
    }
    if (Inside(pointA)==ifalse || Inside(pointB)==ifalse){
        cout << "Error: pointA and pointB should be feasible" << endl;
        return list<Node*>();
    }

    // SIVIA
    int k=1;
    list<Node*> L;
    L.push_back (pX);

    while ( !L.empty() ){
        pX=L.front();   L.pop_front();

        if (Subset(pointA,pX->space)) pva=pX;
        if (Subset(pointB, pX->space)) pvb=pX;

        iboolean test=Inside(pX->space);
        if (test==ifalse){
            //Box Out
            pX->category = ifalse;
            paving.push_back(pX);
            vibes::drawBox(pX->space[1].inf,pX->space[1].sup,pX->space[2].inf,pX->space[2].sup,"red[red]");
        }
        else if	(test==itrue){
            //Box IN
            pX->category = itrue;
            paving.push_back(pX);
            vibes::drawBox(pX->space[1].inf,pX->space[1].sup,pX->space[2].inf,pX->space[2].sup,"green[green]");
        }
        else if (pX->space.Width()< epsilon){
            //Box Perhaps
            pX->category = iperhaps;
            paving.push_back(pX);
            vibes::drawBox(pX->space[1].inf,pX->space[1].sup,pX->space[2].inf,pX->space[2].sup,"blue[blue]");
        }
        else{
            vibes::drawBox(pX->space[1].inf,pX->space[1].sup,pX->space[2].inf,pX->space[2].sup,"blue");
            Node* pX1 = new Node();    Node* pX2 = new Node();
            BisectNode(pX, pX1, pX2);
            delete pX;  pX=0;
            L.push_back(pX1);  L.push_back(pX2);
            k++;    // Incrementation de l'indice k qui copmte le nombre des boxes totales
        }
    }

    list<Node*> path;
    bool outerPath = ShortestPath(true, pva, pvb, path);
    if(outerPath == false){
        cout << "No path" << endl;
        return list<Node*>();
    }
    path.clear();
    bool innerPath = ShortestPath(false, pva, pvb, path);
    if(innerPath == false){
        cout << "Failure to find a path with epsilon = " << epsilon << endl;
        return list<Node*>();
    }
    t = clock() - t;
    drawPath(path);

    cout << "Sivia To found a feasible path is finished ..."<<endl;
    cout <<"Number of boxes with Sivia : " <<k<<endl;
    cout <<"path.size() with Sivia : " <<path.size()<<endl;
    printf ("time : %f seconds.\n\n",((float)t)/CLOCKS_PER_SEC);
    return path;
}






list<Node*> FeasiblePath2 (box pointA, box pointB, box startSpace){

    time_t t = clock();

    Node* pva(0);    Node* pvb(0);  Node* pvp(0);
    Node* pX = new Node();
    pX->space = startSpace;

    if (!(Subset(pointA,pX->space)) || !(Subset(pointB,pX->space))){
        cout << "Error: pointA and pointB should belong to P0" << endl;
        return list<Node*>();
    }
    if (Inside(pointA)==ifalse || Inside(pointB)==ifalse){
        cout << "Error: pointA and pointB should be feasible" << endl;
        return list<Node*>();
    }

    list<Node*> outpath;
    list<Node*> path;

    int k=1;
    list<Node*> L;
    L.push_back(pX);

    while(k<50000){

        while ( !L.empty() ){
            pX=L.front();   L.pop_front();

            if (Subset(pointA,pX->space)) pva=pX;
            if (Subset(pointB, pX->space)) pvb=pX;

            iboolean test=Inside(pX->space);
            if (test==ifalse){
                //Box Out
                pX->category = ifalse;
                paving.push_back(pX);
                vibes::drawBox(pX->space[1].inf,pX->space[1].sup,pX->space[2].inf,pX->space[2].sup,"red[red]");
            }
            else if	(test==itrue){
                //Box IN
                pX->category = itrue;
                paving.push_back(pX);
                vibes::drawBox(pX->space[1].inf,pX->space[1].sup,pX->space[2].inf,pX->space[2].sup,"green[green]");
            }
            else{
                //Box Perhaps
                pX->category = iperhaps;
                paving.push_back(pX);
                vibes::drawBox(pX->space[1].inf,pX->space[1].sup,pX->space[2].inf,pX->space[2].sup,"blue");
            }
        }

        bool outerPath = ShortestPath(true, pva, pvb, outpath);
        if(outerPath == false){
            cout << "No path" << endl;
            return list<Node*>();
        }
        bool innerPath = ShortestPath(false, pva, pvb, path);
        if(innerPath == true){
            t = clock() - t;
            drawPath(path);
            cout << "Cameleon To found a feasible path is finished ..."<<endl;
            cout <<"Number of boxes with Camelon : " <<k<<endl;
            cout <<"path.size() with Cameleon : " <<path.size()<<endl;
            printf ("time : %f seconds.\n\n",((float)t)/CLOCKS_PER_SEC);
            return path;
        }
        path.clear();

        for(list<Node*>::iterator it = outpath.begin(); it != outpath.end(); it++){
            pvp = (*it);
            if (pvp->category==iperhaps){
                Node* pX1 = new Node();    Node* pX2 = new Node();
                BisectNode(pvp, pX1, pX2);
                for(list<Node*>::iterator j = paving.begin(); j != paving.end(); j++){
                    if((*j)==pvp){
                        paving.erase(j);
                        break;
                    }
                }
                delete pvp;  pvp=0;
                L.push_back(pX1);  L.push_back(pX2);
                k++;    // Incrementation de l'indice k qui compte le nombre des boxes totales
            }
        }
        outpath.clear();
    }
}

void drawMotion(list<Node*>& path){

    for(list<Node*>::iterator it = path.begin(); it != path.end(); it++){
        box bo = (*(it))->space;
        double x = 0.1*Center(bo[1]);    double theta = Center(bo[2])*PI/180;

        vibes::clearFigure("Motion");

        for (int j=0;j<=ax.size()-1;j++){
            vector<double> x, y;
            x.push_back(ax[j]);  x.push_back(bx[j]);
            y.push_back(ay[j]); y.push_back(by[j]);
            vibes::drawLine(x, y,"black");
        }

        vector<double> linex, liney;
        linex.push_back(-20);   linex.push_back(57);
        liney.push_back(0);  liney.push_back(0);
        vibes::drawLine(linex,liney,"red");

        DrawPolygon(x,0,theta,sx,sy,"blue");
        vibes::drawCircle(x,0,0.5,"red[red]");

        std::this_thread::sleep_for(std::chrono::milliseconds(400));
    }
}

//----------------------------------------------------------------

int main(){
    vibes::beginDrawing();

    initGlobalVariables();

    box pointA = box(interval(0),interval(0));
    box pointB = box(interval(17)*10,interval(0));
    box P0 = box(interval(-20,57)*10,interval(-1.4,2.7)*180/PI); //Domaine initiale (x-coordinate (*10), heading angle (in deg))

    //cout << P0 << endl;

    vibes::newFigure("FeasiblePath1");
    list<Node*> path1 = FeasiblePath1 (pointA, pointB, P0, 0.1*10);
    paving.clear();
    vibes::newFigure("FeasiblePath2");
    vibes::drawBox(P0[1].inf,P0[1].sup,P0[2].inf,P0[2].sup,"black");
    list<Node*> path2 = FeasiblePath2 (pointA, pointB, P0);

    vibes::newFigure("Motion");
    drawMotion(path2);


    vibes::endDrawing();
}
