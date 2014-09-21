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
    char residue[4];
    char atom_name[5];
    int serial_id;
    int aa_id;
};


struct Peptide {
    char residue[4];
    Point N;
    Point CA;
    Point C;
    Point O;
};

void print_atom(struct Atom);
void print_coor(float*);

float dist(struct Point, struct Point);
float angle(struct Point, struct Point, struct Point);
float torsion_angle(struct Point, struct Point, struct Point, struct Point);
vector<struct Atom> load_pdb_file(char*);
vector<struct Peptide> get_backbone(vector<struct Atom>);
int print_backbone_statistics(vector<struct Peptide>);
void tovector(Point, Point, float*);
void vadd(float*, float*, float*);
void vsub(float*, float*, float*);
void cross(float*, float*, float*);
float dot(float*, float*);
float edge_length(float* a);

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
               "%*s %d %5s %3s %*s %d %f %f %f",
               &atom.serial_id, 
               atom.atom_name,
               atom.residue,
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

void print_coor(float* x){
    printf("%f,%f,%f\n", x[0], x[1], x[2]);
}

void tovector(Point a, Point b, float * v){
    v[0] = b.x - a.x;
    v[1] = b.y - a.y;
    v[2] = b.z - a.z;
}

void vadd(float* a, float* b, float * out){
    out[0] = a[0] + b[0];
    out[1] = a[1] + b[1];
    out[2] = a[2] + b[2];
}

void vsub(float* a, float* b, float* out){
    out[0] = a[0] - b[0];
    out[1] = a[1] - b[1];
    out[2] = a[2] - b[2];
}

void cross(float* a, float* b, float* out){
    out[0] = a[1] * b[2] - a[2] * b[1];
    out[1] = a[2] * b[0] - a[0] * b[2];
    out[2] = a[0] * b[1] - a[1] * b[0];
}

float edge_length(float* a){
    return float(pow(a[0] * a[0] + a[1] * a[1] + a[2] * a[2], 0.5));
}

float dot(float* a, float* b){
    return a[0] * b[0] +
           a[1] * b[1] +
           a[2] * b[2];
}

vector<struct Peptide> get_backbone(vector<struct Atom> atoms){
    vector<struct Peptide> backbone;
    int j = 0;
    for(int i = 0; i < atoms.size(); i++){
        if(strcmp(atoms[i].atom_name, "N") == 0){
            struct Peptide pep;
            backbone.push_back(pep);
            backbone[j].N = atoms[i].pos;
            // TODO I'll be damned to hell for doing this ...
            memcpy(backbone[j].residue, atoms[i].residue, 4);
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


/*
 * Setting points a, b, c in a plane, calculates torsion angle of d about the bc axis
 */
float torsion_angle(struct Point a, struct Point b, struct Point c, struct Point d)
{
    float ab[3], ac[3], Ua[3];
    tovector(a, b, ab);
    tovector(a, c, ac);
    cross(ab, ac, Ua);

    float dc[3], db[3], Ub[3];
    tovector(d, c, dc);
    tovector(d, b, db);
    cross(dc, db, Ub);

    float t_angle = acos( dot(Ua, Ub) / (edge_length(Ua) * edge_length(Ub)) );

    return t_angle * (180 / PI);
}


int print_backbone_statistics(vector<struct Peptide> b){
    printf(
        "%s %s %s %s %s %s %s %s %s %s\n",
        "ind", "N-CA", "CA-C", "C-N", "tau", "som", "phi", "psi", "omega", "aa"
    );
    float psi, phi;
    for(int i = 0; i < b.size(); i++){
        printf(
            "%d %f %f %f %f %f %f %f %f %s\n",
            // Peptide number
            i + 1,
            // N -> CA bond length
            dist(b[i].N, b[i].CA),
            // CA -> C bond length
            dist(b[i].CA, b[i].C),
            // C -> N bond length
            (i+1) < b.size() ? dist(b[i].C, b[i+1].N) : 0,
            // N-CA-C angle
            angle(b[i].N, b[i].CA, b[i].C),
            // CA-C-O angle
            angle(b[i].CA, b[i].C, b[i].O),
            // phi torsion angle - C- _ N _ CA _ C
            i > 0 ? torsion_angle(b[i-1].C, b[i].N, b[i].CA, b[i].C) : (float)999,
            // psi torsion angle - N _ CA _ C _ N+
            (i + 1) < b.size() ? torsion_angle(b[i].N, b[i].CA, b[i].C, b[i+1].N) : (float)999,
            // omega torsion angle - CA _ C _ N+ _ CA+
            (i + 1) < b.size() ? torsion_angle(b[i].CA, b[i].C, b[i+1].N, b[i+1].CA) : (float)999,
            // residue three letter name
            b[i].residue
        );
    }
    return 1;
}


void print_atom(struct Atom a)
{
    cout << a.atom_name;
    cout << "\t(" << a.pos.x << "," << a.pos.y << "," << a.pos.z << ")"; 
    cout << "\t(" << a.serial_id << "," << a.aa_id << ")" << endl;
}


