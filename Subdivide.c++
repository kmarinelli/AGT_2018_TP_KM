#include <iostream>
using namespace std;

#include <iomanip> 
#include <iostream> 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>  

#define BOUNDEXPANSION 0.05

#define max(a,b) (a>b?a:b)
#define min(a,b) (a<b?a:b)

#define L 4

#define X 0
#define Y 1
#define Z 2

class Point
{
  public:
   long x;
   long y;
   long z;

   Point(long x, long y, long z) : x(x), y(y),z(z) {};
   Point();
   void Set(long x, long y, long z);
   void Print();

   Point operator*(long v);
   Point operator/(long v);
   Point operator+(Point p);
   Point operator-(Point p);
};

Point::Point()
{
   this->x=0;
   this->y=0;
   this->z=0;
}

void Point::Print()
{
   cout<<x<<", "<<y<<","<<z;
}

void Point::Set(long x, long y, long z)
{
   this->x=x;
   this->y=y;
   this->z=z;
}

Point Point::operator*(long v)
{
   Point result;
 
   result.Set(x*v,y*v,z*v);

   return result;

}

Point Point::operator/(long v)
{
   Point result;

   result.Set(x/v,y/v,z/v);

   return result;

}


Point Point::operator+(Point p)
{
   Point result;

   result.Set(x+p.x, y+p.y, z+p.z);

   return result;

}

Point Point::operator-(Point p)
{
   Point result;

   result.Set(x-p.x, y-p.y, z-p.z);

   return result;

}

class RTree
{
   RTree *Node;
   RTree *Left;
   RTree *Right;
   long xmin; // X ordinate min bound.
   long xmax; // X ordinate max bound.
   long ymin; // Y ordinate min bound.
   long ymax; // Y ordinate max bound.
   long zmin; // Z ordinate min bound.
   long zmax; // Z ordinate max bound.

   public:

   RTree();
};

RTree:: RTree()
{
   Node=NULL;
   Left = NULL;
   Right = NULL;
   xmin = 0;
   xmax = 0;
   ymin = 0;
   ymax = 0;
   zmin = 0;
   zmax = 0;
};

class Bezier
{
   public:
   int degree;
   Point *P;

   long xmin;
   long ymin;
   long zmin;
   long xmax;
   long ymax;
   long zmax;

   double totalCurvature;

   Bezier(int degree, Point *p);
   Bezier();
   Point Evaluate(double t);
   void PrintControlPoints();
   Bezier *DeCasteljau();
   void CreateBounds();
   void ComputeTotalCurvature();
   double Derivative(long degree, double t);
   double Curvature(double t); 
   void PrintBounds();
   string Monotonicity();
   bool isConvex();
};

Bezier *Bezier::DeCasteljau()
{
   Bezier *B;
   Point **T;
   int i;
   int j;

   T = (Point **) new Point*[degree+1];

   T[0] = new Point[degree+1];
   for(i=0;i<degree+1;i++)
   {
      T[0][i]=P[i];
   }

   for(j=1;j<=degree;j++)
   {
      T[j]=new Point[degree-j+1];
      for(i=0;i<=degree-j;i++)
      {
         T[j][i]=(T[j-1][i+1]+T[j-1][i])/2;
      }
   }

   B=new Bezier[2];
   B[0].degree=degree;
   B[1].degree=degree;

   B[0].P = new Point[degree+1];

   B[1].P = new Point[degree+1];

   for(j=0;j<=degree;j++)
   {
      B[0].P[j]=T[j][0];
      B[1].P[j]=T[degree-j][j];
   }
   B[0].CreateBounds();
   B[0].ComputeTotalCurvature();
   B[1].CreateBounds();
   B[1].ComputeTotalCurvature();

   for(j=0;j<=degree;j++)
   {
      delete [] T[j];
   }
  
   delete [] T;

   return B;
}

Bezier::Bezier()
{
   this->degree=0;
   this->P = NULL;
   CreateBounds();
   ComputeTotalCurvature();
}

Bezier::Bezier(int degree, Point *p)
{
   int i;

   this->degree=degree;
   this->P=p;

   CreateBounds();
   ComputeTotalCurvature();
}

void Bezier::PrintBounds()
{
   cout << "Bounds:"<<endl;
   cout << "X : "<<setw(10)<<setprecision( 16 )<<xmin<<", "<<setw(10)<<setprecision( 10 )<<xmax<<endl;
   cout << "Y : "<<setw(10)<<setprecision( 16 )<<ymin<<", "<<setprecision( 10 )<<setw(10)<<ymax<<endl;
   cout << "Z : "<<setw(10)<<setprecision( 16 )<<zmin<<", "<<setw(10)<<setprecision( 10 )<<zmax<<endl;
}

void Bezier::CreateBounds()
{
   int i;

   if( P == NULL)
   {
      return;
   }

   xmin=P[0].x;
   ymin=P[0].y;
   zmin=P[0].z;
   xmax=P[0].x;
   ymax=P[0].y;
   zmax=P[0].z;

   for(i=1;i<=degree;i++)
   {
      xmin=min(xmin,P[i].x);
      ymin=min(ymin,P[i].y);
      zmin=min(zmin,P[i].z);

      xmax=max(xmax,P[i].x);
      ymax=max(ymax,P[i].y);
      zmax=max(zmax,P[i].z);
   }
}

int sign(long x)
{
   if( x < 0) return -1;
   if( x == 0) return 0;
   return 1;
}

double rsign(double x)
{
   if( x < 0) return -1;
   if( x == 0) return 0;
   return 1;
}

string Bezier::Monotonicity()
{
   long x,y,z;
   long sumx, sumy, sumz;
   int i;
   string s="";

   sumx=0;
   sumy=0;
   sumz=0;
   
   for(i=0;i<degree;i++)
   {
      sumx+=sign((P[i].x-P[i+1].x));
      sumy+=sign((P[i].y-P[i+1].y));
      sumz+=sign((P[i].z-P[i+1].z));
   }
   if( ((int)abs(sumx))==degree) s+="X";
   if( ((int)abs(sumy))==degree) s+="Y";
   if( ((int)abs(sumz))==degree) s+="Z";
   if( s=="") s="NONE";
   return s;
}

double dotproduct(double v1[3], double v2[3])
{
   return v1[X]*v2[X] + v1[Y]*v2[Y]+v1[Z]*v2[Z];
}

double zcrossproduct(double v1[3], double v2[3])
{
   return v1[X]*v2[Y]-v1[Y]*v2[X] ;
}

bool::Bezier::isConvex()
{
   int i;
   int n;
   Point a;
   Point b;
   Point c;
   double angle;
   double lastAngle;
   bool convex;

   double v1[3];
   double v2[3];
   double d;

   n=degree+1; // Closed loop is polygon degree + 1 .
   if( n < 4) return true;

   //cout << "nangles = "<<n<<endl;
   lastAngle=0;
   convex=true;
   for(i=0;i<n;i++)
   {
      a=P[i%n];
      b=P[(i+1)%n];
      c=P[(i+2)%n];

      // Vector 1
      v1[X]=b.x-a.x;
      v1[Y]=b.y-a.y;
      v1[Z]=b.z-a.z;

      // Vector length.
      d=sqrt(v1[X]*v1[X]+v1[Y]*v1[Y]+v1[Z]*v1[Z]);

      // scale by 1/d.
      v1[X]=v1[X]/d;
      v1[Y]=v1[Y]/d;
      v1[Z]=v1[Z]/d;
      //cout << "length v1 = "<<v1[X]*v1[X]+v1[Y]*v1[Y]+v1[Z]*v1[Z]<<endl;

      //vector 2
      v2[X]=c.x-b.x;
      v2[Y]=c.y-b.y;
      v2[Z]=c.z-b.z;

      // vector length
      d=sqrt(v2[X]*v2[X]+v2[Y]*v2[Y]+v2[Z]*v2[Z]);

      // scale by 1/d.
      v2[X]=v2[X]/d;
      v2[Y]=v2[Y]/d;
      v2[Z]=v2[Z]/d;
      //cout << "length v2 = "<<v2[X]*v2[X]+v2[Y]*v2[Y]+v2[Z]*v2[Z]<<endl;

      angle=dotproduct(v1,v2);
      angle=zcrossproduct(v1,v2);
      //cout << "Angle = "<<angle<<" Last Angle = "<<lastAngle<<endl;
      //cout << " v2="<<v2[X]<<" "<<v2[Y]<<" "<<v2[Z]<<" "<<endl;
      if( i!=0 && rsign(angle) != rsign(lastAngle))
      {
         //cout << "Not Convex!"<<endl;
         convex=false;
         //return false;
      }
      lastAngle=angle;
   }
   return convex;
}

double Bezier::Derivative(long degree, double t)
{
   long d;
   d=0;

   return d;
}

/*
float dBezier(int n, float t, float P[],int d)
{
   float x;
   float *dP;
   int i;

   if( n <= 0) return 1;

   if(d==0) return Bezier(n,t,P);

   dP=(float *)malloc(sizeof(float)*(n));

   // Computer derivative 
   for(i=0;i<n;i++)
   {
      dP[i]=P[i+1]-P[i];
      //printf("dP[%d] = %f\n",i,dP[i]);
   }

   x=n*dBezier(n-1,t,dP,d-1);

   //printf("dBezier(%d %f, %d)=%f\n",n,t,d,x);

   return x;
}
*/

double Bezier::Curvature(double t)
{
   return Derivative(2,t);
}

void Bezier::ComputeTotalCurvature()
{
   double t;
   double deltaT=0.001;
   totalCurvature=0;
   if( P==NULL) return;

   for(t=0;t<1.0;t+=deltaT)
   {
      totalCurvature=deltaT*(Curvature(t));
   }

}

void Bezier::PrintControlPoints()
{
  int i;

   cout << "{"<<endl;
   for(i=0;i<=degree;i++)
   {
      cout<<"{ "; 
      P[i].Print();
      cout <<" }";
      if( i<degree) cout <<",";
      cout << endl;
   }
   cout <<"}"<<endl;

}

int Factorial(int f)
{
   int F;
   int i;
   
   F = 1;
   for(i=1;i<=f;i++)
   {
     F*=i; 
   }

   return F;
}

int Binomial(int x, int y)
{
   return Factorial(x)/(Factorial(y)*Factorial(x-y));
}

Point Bezier::Evaluate(double t)
{
   int i;
   Point sum;
   Point Q;
   double tk;
   double s;
   int b;

   tk=1;

   sum.Set(0,0,0);

   for(i=0;i<=degree;i++)
   {
      b=Binomial(degree,i);
      tk=pow(t,(double)i) * pow( (1-t), (double)(degree-i));
      s=b*tk;

      Q=P[i]*s;

      sum=sum+Q;
   }

   return sum;
}

int main()
{
   Point *p;

   Bezier *C;

   Bezier ***Q;

   int degree;
   int i;
   int j;

   string prefix;
   string setName;

   double scale=pow(2,29);
   //double scale=1.0;
   //double scale=1024;
   

   Q= new Bezier**[L+1];
   for(i=0;i<(L+1);i++)
   {
      Q[i]=new Bezier*[(int)pow(2,i)];
   }

   degree=6;

   p=new Point[degree+1];

   // Unknot
   setName="Unknot";
   prefix="VU";
   p[0].Set(0, 9, 20);
   p[1].Set(-15, -95, -50); 
   p[2].Set(40, 80, -20);
   p[3].Set(-10, -60, 58);
   p[4].Set(-60, 30, 20);
   p[5].Set(40, -60, -60);
   p[6].Set(0, 9, 20);
/*
   // Trefoil
   setName="Trefoil";
   prefix="VT";
   p[0].Set(0, 9, 20);
   p[1].Set(-15, -95, -50);
   p[2].Set(40, 80, -20);
   p[3].Set(10, -60, 58);
   p[4].Set(-60, 30, 20);
   p[5].Set(40, -60, -60);
   p[6].Set(0, 9, 20);

*/
   p[0]=p[0]*scale;
   p[1]=p[1]*scale;
   p[2]=p[2]*scale;
   p[3]=p[3]*scale;
   p[4]=p[4]*scale;
   p[5]=p[5]*scale;
   p[6]=p[6]*scale;

   Q[0][0] = new Bezier(degree,p);
   cout <<"(*Original "<<setName<<" Control Points*)"<<endl;
   cout<<prefix<<"0s0 = "<<endl;
   Q[0][0]->PrintControlPoints();
   cout <<endl<<endl;;

   for(i=1;i<(L+1);i++)
   {
      cout << "(*Subdivision "<<i<<" *)"<<endl;
      for(j=0;j<(int)pow(2,(i-1));j++)
      {
         
         C=Q[i-1][j]->DeCasteljau();

         cout <<prefix<<i<<"s"<<j*2<<" = "<<endl;
         Q[i][j*2]=&(C[0]);
         Q[i][j*2]->PrintControlPoints();
         cout << "(* Monotonicity "<<prefix<<" l:"<<i<<" h:"<<j*2<<"  = "<<Q[i][j*2]->Monotonicity()<<" *)"<<endl;
         cout << "(* 2D Convexity "<<prefix<<" l:"<<i<<" h:"<<j*2<<"  = "<<Q[i][j*2]->isConvex()<<" *)"<<endl;
         cout << endl;

         cout <<prefix<<i<<"s"<<j*2+1<<" = "<<endl;
         Q[i][j*2+1]=&(C[1]);
         Q[i][j*2+1]->PrintControlPoints();
         cout << "(* Monotonicity "<<prefix<<" l:"<<i<<" h:"<<j*2+1<<"  = "<<Q[i][j*2+1]->Monotonicity()<<" *)"<<endl;
         cout << "(* 2D Convexity "<<prefix<<" l:"<<i<<" h:"<<j*2+1<<"  = "<<Q[i][j*2+1]->isConvex()<<" *)"<<endl;

         cout << endl;

      }
   }

   for(i=0;i<(L+1);i++)
   {
         delete Q[i];
   }


   delete Q;
   delete p;


   exit(0);
}
