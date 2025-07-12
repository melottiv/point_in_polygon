#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <random>
#include <chrono>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
using namespace std;



// ANTICLOCKWISE ORIENTATION

struct Point {
    double x, y;
};

/**
 * @brief Contains values referring to a line
 *
 * ax + by + c = 0
 */
struct Line{
    double a,b,c;
};

/**
 * @brief Turns a Point struct into an array
 *
 * @param p Point to convert
 * @param a array to fill
 * 
 * @return void
 */

void pointToArray(Point p, double *a){
    a[0]=p.x;
    a[1]=p.y;
    a[2]=1;
    return;
}

/**
 * @brief Turns a Line struct into an array
 * 
 * @param l Line to convert
 * @param a array to fill
 * 
 * @return void
 */
void lineToArray(Line l, double *a){
    a[0]=l.a;
    a[1]=l.b;
    a[2]=l.c;
    return;
}

/**
 * @brief Turns an array with 3 elements into a Point, normalizing by the last element
 *
 * @return Point
 */
Point arrayToPoint(double *a){
    if(a[2]==0){
        cout<<"ERROR\nPoint goes to infinity"<<endl;
        return {-9999999,-9999999};
    }
    return {a[0]/a[2],a[1]/a[2]};
}

/**
 * @brief Turns an array with 3 elements into a Line
 *
 * @return Point
 */
Line arrayToLine(double *a){
    return{a[0],a[1],a[2]};
}

/**
 * @brief computes distance between two points
 * 
 * @return double
 */
double getDistance(Point p1,Point p2){
    double x=p1.x-p2.x;
    double y=p1.y-p2.y;
    return sqrt(x*x+y*y);
}

/**
 * @brief generates a set of random points
 * @param c: centroid
 * @param R radius or generation
 * 
 * @return vector<Point>
 */
vector<Point> generatePoints(Point c, double R, int n_points){
    double x1=c.x+R;
    double x2=c.x-R;
    double y1=c.y+R;
    double y2=c.y-R;
    uniform_real_distribution<>disx(x1,x2);
    uniform_real_distribution<>disy(y1,y2);
    random_device rd; 
    mt19937 gen(rd());

    vector<Point> v(n_points);

    for(int i=0;i<n_points;i++){
        v[i].x=disx(gen);
        v[i].y=disy(gen);
    
       while(getDistance(v[i],c)>R){
            v[i].x=disx(gen);
            v[i].y=disy(gen);
        }
       
    }
    return v;
}


/**
 * @brief computes the centroid given vector of Points
 * 
 * @return void
 */
Point getCentroid(vector<Point> v){
    double x=0;
    double y=0;
    if(v.size()<=0){
        cout<<"ERRORE\nil vettore non ha lunghezza positiva"<<endl;
        return{-999999,-999999};
    }
    for(unsigned int i=0;i<v.size();i++){
        x+=v[i].x;
        y+=v[i].y;
    }
    return{x/v.size(),y/v.size()};
}

/**
 * @brief generates a regular polygon
 * @param m number of vertices
 * @param c centroid
 * @param r radius of the polygon
 * 
 * @return void
 */
vector<Point> generateRegularPolygon(int m,Point c,double r){
    vector<Point> points(m);
    vector<double> angles(m);
    for (int i = 0; i < m; ++i) {
        angles[i]=i*2*M_PI/m;
    }
    for (int i=0;i<m;i++){
        points[i].x=c.x+r*cos(angles[i]);
        points[i].y=c.y+r*sin(angles[i]);
    }

    return points;
}

/**
 * @brief generates a random polygon
 * @param m number of vertices
 * @param c centroid
 * @param r radius of the polygon
 * 
 * @return void
 */
vector<Point> generatePolygon(int m,Point c,double r){
    vector<Point> points(m);
    vector<double> angles(m);
    random_device rd; 
    mt19937 gen(rd());
    uniform_real_distribution<>dis(0,2*M_PI);
    for (int i = 0; i < m; ++i) {
        angles[i]=dis(gen);
    }

    sort(angles.begin(),angles.end());

    for (int i=0;i<m;i++){
        cout<<angles[i]<<endl;
        points[i].x=c.x+r*cos(angles[i]);
        points[i].y=c.y+r*sin(angles[i]);
    }

    return points;
}


/**
 * @brief computes the vectoral product between two arrays, storing results in the last one
 * 
 * This function solves the linear system Ax=0, works both for finding a line through points and for a point given two lines
 * 
 * @return void
 */
void vecProduct(double a1[3], double a2[3], double ans[3]) {

    ans[0] = a1[1] * a2[2] - a1[2] * a2[1]; // a
    ans[1] = a1[2] * a2[0] - a1[0] * a2[2]; // b
    ans[2] = a1[0] * a2[1] - a1[1] * a2[0]; // c
}

/**
 * @brief which side of the line does a point fall on
 * 
 * @return 1 if internal, -1 if external, 0 if on the line
 */
int whichSide(Point p, Line l){
    double val=p.x*l.a+p.y*l.b+l.c;
    if(val>0) return 1;
    if(val<0) return -1;
    return 0;
}

/**
 * @brief Finds a line given two points
 *
 * @return Line
 */
Line lineThroughPoints(Point p1, Point p2){
    double homop1[3], homop2[3], homol[3];
    pointToArray(p1,homop1);
    pointToArray(p2,homop2);

    vecProduct(homop1,homop2,homol);

    return arrayToLine(homol);
}

/**
 * @brief Finds an intersection Point given two Lines
 *
 * @return Point
 */
Point intersectionPoint(Line l1, Line l2){
    double homol1[3], homol2[3], homop[3];
    lineToArray(l1,homol1);
    lineToArray(l2,homol2);

    vecProduct(homol1,homol2,homop);

    return arrayToPoint(homop);
}

/**
 * @brief Tests if a point is inside a polygon with O(N) complexity
 *
 * @return bool
 */
bool testPolygonON(Point point,vector<Point> polygon){
    polygon.push_back(polygon[0]);
        for(unsigned int i=0;i<polygon.size()-1;i++){
            Point p1=polygon[i];
            Point p2=polygon[i+1];
            // draws the line passing through consecutive points
            Line l=lineThroughPoints(p1,p2);
            // checks which side of the line the point falls on
            if(whichSide(point,l)<0){
                return false;
            }
    }
    return true;
}

bool logSearch(vector<Point> polygon, Point point, int i,int j){
    Point C=getCentroid({polygon[0],polygon[floor(polygon.size()/2)]});
    Line q=lineThroughPoints(C,point);
    int k;
    while((j-i)>1){
        k=(j+i)>>1;     
        if(i==k||k==j) {
            cout<<"ERRORE"<<endl;
            return false;
            }
        if(whichSide(polygon[k],q)>=0) j=k;
        else i=k;        
    }
    Line a= lineThroughPoints(polygon[i],polygon[i+1]);
    if(whichSide(point,a)>0) return true;
    else return false;
}   

bool testPolygonLogN(Point point, vector<Point> polygon){
    int k=polygon.size()>>1;           
    Line p=lineThroughPoints(polygon[0],polygon[k]);
    if(whichSide(point,p)>=0) {
        polygon.push_back(polygon[0]);
        return logSearch(polygon,point,k,polygon.size()-1);}
    else 
    return logSearch(polygon,point,0,k);
}


string formatNumber(double number) {
    ostringstream oss;
    oss << std::fixed << std::setprecision(2) << number;
    string result = oss.str();
    replace(result.begin(), result.end(), '.', ','); 
    return result;
}

void testFunction(ofstream &file, Point c,double r,int m, double R, int n_points, bool (*func)(Point, vector<Point>)){
    vector<Point> test=generatePoints(c,R,n_points);
    vector<Point> polygon=generateRegularPolygon(m,c,r);

    //bool results[n_points];
    int val=0;

    for(int i=0;i<test.size();i++){
            val+=func(test[i],polygon);
    }

    file<<val<<";"<<n_points<<endl;

    return;
}

int main() {

    ofstream file("newfile.csv");

   Point c={1,1};       //centroid
    double r=3;         //radius of polygon
    int m=10;            //number of vertices of the polygon
    double R=5;         //radius of generated points
    int n_points=1000;  //number of generated points
    int n_tests=1000;

    file<<"O(N) algorithm"<<endl<<r<<";"<<R<<endl;
    for(int i=1;i<n_tests+1;i++){
        testFunction(file,c,r,m,R,10*i,&testPolygonLogN);
    }

    file.close();
    return 0;
}