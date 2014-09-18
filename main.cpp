// main.cpp
// Zebulun Arendsee
// BCB569 HW1
// September 2014

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <string>
#include <string.h>
#include <vector>
#include <stdio.h>
using namespace std;

struct Point {
    float x;
    float y;
    float z;
};

struct Atom {
    Point pos;
    char atom_name[5];
    int serial_id;
    int aa_id;
};

float dist(struct Point, struct Point);
float angle(struct Point, struct Point, struct Point);
void print_atom(struct Atom);
vector<struct Atom> load_pdb_file(char*);

int main(int argc, char* argv[]){
    // open a file, print its contents
    if(argc == 2){
        vector<struct Atom> atoms = load_pdb_file(argv[1]); 
    } else {
        cerr << "Please provide a filename" << endl;
        return 1;
    }
    return 0;
}

vector<struct Atom> load_pdb_file(char* filename){
    vector<Atom> atoms;
    ifstream infile(filename, ios::in);
    string line;
    while(getline(infile, line)){
        if(line.substr(0,4) != "ATOM")
            continue;
        struct Atom atom;
        sscanf(line.c_str(),
               "%*s %d %5s %*s %*s %d %f %f %f",
               &atom.serial_id, 
               atom.atom_name,
               &atom.aa_id,
               &atom.pos.x, &atom.pos.y, &atom.pos.z
              );
        atoms.push_back(atom);
    }
    return atoms;
}

float dist(struct Point a, struct Point b){
    return sqrt((a.x - b.x) * (a.x - b.x) + 
                (a.y - b.y) * (a.y - b.y) +
                (a.z - b.z) * (a.z - b.z));
}

float angle(struct Point a, struct Point b, struct Point c){
    return 0;
}

void print_atom(struct Atom a){
    cout << a.atom_name;
    cout << "\t(" << a.pos.x << "," << a.pos.y << "," << a.pos.z << ")"; 
    cout << "\t(" << a.serial_id << "," << a.aa_id << ")" << endl;
}


/*
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
 * Rotater function
 * 1. find a rotation matrix
 * 2. select all points to be rotated (i.e. those downstream of rotating bond), express as a matrix
 * 3. multiply rotation matrix by point matrix
 */
// void rotate(int bondid, Matrix protein, rotation angle)

// struct AA {
//     Point N;
//     Point CA;
//     Point C;
//     Point O;
//     Point H;
// };
// 
// struct edge {
//     Point a;
//     Point b;
// };
