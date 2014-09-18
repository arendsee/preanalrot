// main.cpp
// Zebulun Arendsee
// BCB569 HW1
// September 2014

#include <math.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <vector>
#include <stdio.h>
using namespace std;

const float PI = acos(-1);


struct Point
{
    float x;
    float y;
    float z;
};


struct Atom
{
    Point pos;
    char atom_name[5];
    int serial_id;
    int aa_id;
};


struct Peptide {
    Point N;
    Point CA;
    Point C;
    Point O;
};


float dist(struct Point, struct Point);
float angle(struct Point, struct Point, struct Point);
void print_atom(struct Atom);
vector<struct Atom> load_pdb_file(char*);
vector<struct Peptide> get_backbone(vector<struct Atom>);
int print_backbone_statistics(vector<struct Peptide>);
Point copy_point(Point);


int main(int argc, char* argv[])
{
    vector<struct Atom> atoms;
    // open a file, print its contents
    if(argc == 2){
        atoms = load_pdb_file(argv[1]); 
    } else {
        cerr << "Please provide a filename" << endl;
        return 1;
    }
    vector<struct Peptide> backbone = get_backbone(atoms);
    print_backbone_statistics(backbone);
    return 0;
}


vector<struct Atom> load_pdb_file(char* filename)
{
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


float dist(struct Point a, struct Point b)
{
    return sqrt((a.x - b.x) * (a.x - b.x) + 
                (a.y - b.y) * (a.y - b.y) +
                (a.z - b.z) * (a.z - b.z));
}


/*
 * Solve for the angle between points a, b, and c
 */
float angle(struct Point a, struct Point b, struct Point c)
{
    float ab = dist(a, b);
    float bc = dist(b, c);
    float ac = dist(a, c);
    float theta = acos((ab * ab + bc * bc - ac * ac ) / (2 * bc * ac));
    return theta * (180 / PI);
}

/*
 * Setting points a, b, c in a plane, calculates torsion angle of d about the bc axis
 */
float torsion_angle(struct Point a, struct Point b, struct Point c, struct Point d)
{
    float ab = dist(a, b);
    float bc = dist(b, c);
    float ac = dist(a, c);
    float theta = acos((ab * ab + bc * bc - ac * ac ) / (2 * bc * ac));
    return theta * (180 / PI);
}


void print_atom(struct Atom a)
{
    cout << a.atom_name;
    cout << "\t(" << a.pos.x << "," << a.pos.y << "," << a.pos.z << ")"; 
    cout << "\t(" << a.serial_id << "," << a.aa_id << ")" << endl;
}

vector<struct Peptide> get_backbone(vector<struct Atom> atoms){
    vector<struct Peptide> backbone;
    int j = 0;
    for(int i = 0; i < atoms.size(); i++){
        if(strcmp(atoms[i].atom_name, "N") == 0){
            struct Peptide pep;
            backbone.push_back(pep);
            backbone[j].N = atoms[i].pos;
        }
        else if(strcmp(atoms[i].atom_name, "CA") == 0){
            backbone[j].CA = atoms[i].pos;
        }
        else if(strcmp(atoms[i].atom_name, "C") == 0){
            backbone[j].C = atoms[i].pos;
        }
        else if(strcmp(atoms[i].atom_name, "O") == 0){
            backbone[j].O = atoms[i].pos;
            j++;
        }
    }
    return backbone;
}

Point copy_point(Point p)
{
    Point copy;
    copy.x = p.x;
    copy.y = p.y;
    copy.z = p.z;
    return copy;
}

int print_backbone_statistics(vector<struct Peptide> b){
    for(int i = 0; i < b.size(); i++){
        printf(
            "%d\t%f\t%f\t%f\t%f\t%f\n",
            // Peptide number
            i,
            // N -> CA bond length
            dist(b[i].N, b[i].CA),
            // CA -> C bond length
            dist(b[i].CA, b[i].C),
            // C -> N bond length
            (i+1) < b.size() ? dist(b[i].C, b[i+1].N) : 0,
            // N-CA-C angle
            angle(b[i].N, b[i].CA, b[i].C),
            // CA-C-O angle
            angle(b[i].CA, b[i].C, b[i].O)
        );
    }
    return 1;
}
