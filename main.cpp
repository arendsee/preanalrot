// main.cpp
// Zebulun Arendsee
// BCB569 HW1
// September 2014

#include <stdio.h>
#include <math.h>

double dist(struct Point, struct Point);
double angle(struct Point, struct Point, struct Point);
struct Point makePoint(double, double, double);

struct Point {
    double x;
    double y;
    double z;
};

struct AA {
    Point N;
    Point CA;
    Point C;
    Point O;
    Point H;
};

struct edge {
    Point a;
    Point b;
};

int main(int argc, char* argv[]){
    struct Point a = makePoint(1,0,0);
    struct Point b = makePoint(0,0,1);
    struct Point c = makePoint(0,1,0);
    printf("%f\n", dist(a, b));
    return 0;
}

struct Point makePoint(double x, double y, double z){
    struct Point p;
    p.x = x;
    p.y = y;
    p.z = z;
    return p;
}

double dist(struct Point a, struct Point b){
    return sqrt((a.x - b.x) * (a.x - b.x) + 
                (a.y - b.y) * (a.y - b.y) +
                (a.z - b.z) * (a.z - b.z));
}

double angle(struct Point a, struct Point b, struct Point c){
    return 0;
}


/*
 * Functions
 *
 * Atom section columns:
 * 1.  ATOM
 * 2.  serial-id
 * 3.  atom name (N, CA, C, O, etc) + ? alternative location indicator
 * 4.  residue name (MET, THR)
 * 5.  residue sequence id ??
 * 6.  code for residue insertion ??
 * 7.  x
 * 8.  y
 * 9.  z
 * 10. occupancy
 * 11. tempFactor
 * 12. element + charge
 *
 * 1. load the data
 *      in - filename
 *      out - 
 */
