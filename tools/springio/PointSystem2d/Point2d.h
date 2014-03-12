#ifndef _POINT_H_
#define _POINT_H_
#include <stdio.h>
#include <math.h>
#include <set>
#include <map>

using namespace std;
typedef set<int> PtSet;


class Point2d {
   private:
      double x_;
      double y_;
      double lx_; // length of domain in x direction for periodic bc
      double ly_; // length of domain in y direction for periodic bc
      double shearx_; // shift for lees-edwards BC in space (not velocity)

   public:
      Point2d(double x=0, double y=0, double lx=0, double ly=0, double shearx=0) { x_ = x; y_ = y; lx_ = lx; ly_ = ly; shearx_ = shearx; }
      void set (double a, double b, double c, double d, double e) { x_ = a; y_ = b; lx_ = c; ly_ = d; shearx_ = e; }
      void set (double a, double b) { x_ = a; y_ = b; }
      void setshear (double shear) {shearx_ = shear;}
      void set (Point2d q) { x_ = q.x(); y_ = q.y(); lx_ = q.lx_; ly_ = q.ly_; shearx_ = q.shearx_; }
      

      double x() const { return x_; }
      double y() const { return y_; }
      //~ double l() const { return side_; }
      
      // vector addition
      Point2d operator+(Point2d pt2) const {
          return Point2d(x_+pt2.x(),y_+pt2.y(),lx_, ly_, shearx_);
      }
      
      // vector subtraction
      Point2d operator-(Point2d pt2) const {
         double dx = x_ - pt2.x(); 
         double dy = y_ - pt2.y(); 
         
         // periodic in x direction
         if(dx > 0.5*lx_) dx -= lx_; 
         else if(dx < -0.5*lx_) dx += lx_; 
         
         // periodic in y direction; include lees-edwards shift
         if(dy > 0.5*ly_) { dy -= ly_; dx -= shearx_; }
         else if(dy < -0.5*ly_) { dy += ly_; dx += shearx_; }
         
         return Point2d(dx,dy,lx_,ly_,shearx_);
      } 
      
      Point2d operator*(double scalar) const {return Point2d(x_*scalar,y_*scalar,lx_,ly_,shearx_);} // scalar multiplication
      Point2d operator/(double scalar) const {return Point2d(x_/scalar,y_/scalar,lx_,ly_,shearx_);}
      double operator*(Point2d pt2) const { return x_*pt2.x()+y_*pt2.y(); } // dot product
      //~ Point2d operator%(Point2d pt2) const { return Point2d(y_*pt2.l()-side_*pt2.y(), side_*pt2.x()-x_*pt2.l(), x_*pt2.y()-y_*pt2.x()) ; } // cross product
      
      void print() const { printf("%2.8f %2.8f %2.8f %2.8f\n", x_, y_, lx_, ly_); }
      
      double distanceTo(Point2d &q) const {
         Point2d dr = Point2d(x_, y_, lx_, ly_,shearx_) - q;
         return sqrt(dr*dr);
      }
};


class PointNeighbours {
    private:
        PtSet nbrs;
        bool exterior;
    
    public:
        PointNeighbours() { nbrs = PtSet(); exterior = false; } 
    
        // HEX direction-agnostic neighbor functions
        int degree() const { return nbrs.size(); }
        bool is_conn() const { return (degree() >= 6 && exterior == false); } 
        bool is_conn_strict() const { return (degree() == 6 && exterior == false); } 

        void add_nbr(int nbr) { nbrs.insert(nbr); }
        void remove_nbr(int nbr) { nbrs.erase(nbr); }
        bool is_nbr(int nbr) const { return (nbrs.count(nbr) == 1); }

        const PtSet* neighbours() const { return &nbrs; }

        bool isExterior() const { return exterior; }
        void set_exterior() { exterior = true; }
        
};



class Bond {
    private:
        int p1_;
        int p2_;
        double k_;
        double l0_;
        
    public:
        Bond(int p1 = -1, int p2 = -1, double k = 0., double l0 = -1) {p1_ = p1; p2_ = p2; k_ = k; l0_ = l0;} // initialize with -1 
        void setlength(double len0) { l0_ = len0; }
        void setk(double kk) { k_ = kk; }
        
        int p1() const { return p1_; }
        int p2() const { return p2_; }
        double k() const { return k_; }
        double l0() const { return l0_; }
        
        //~ double energy();

     };
     
class Facet {
    private:
        int p1_;
        int p2_;
        int p3_;
        
    public:
        Facet(int p1 = -1, int p2 = -1, int p3 = -1) {p1_ = p1; p2_ = p2; p3_=p3;} // initialize with -1 
        
        int p1() const { return p1_; }
        int p2() const { return p2_; }
        int p3() const { return p3_; }

     };
     
class Dihedral {
    private:
        int p1_;
        int p2_;
        int p3_;
        int p4_;
        double k_;
        double l0_;
        
    public:
        Dihedral(int p1 = -1, int p2 = -1, int p3 = -1, int p4 = -1, double k = 0., double l0 = -1) {p1_ = p1; p2_ = p2; p3_=p3; p4_ = p4;k_ = k; l0_ = l0;} // initialize with -1 
        
        void setlength(double len0) { l0_ = len0; }
        void setk(double kk) { k_ = kk; }
        
        int p1() const { return p1_; }
        int p2() const { return p2_; }
        int p3() const { return p3_; }
        int p4() const { return p4_; }
        double k() const { return k_; }
        double l0() const { return l0_; }

     };
#endif
