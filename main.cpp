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
typedef struct Point pnt;


struct Atom
{
    pnt pos;
    char residue[4];
    char atom_name[5];
    int serial_id;
    int aa_id;
};


struct Peptide {
    char residue[4];
    pnt N;
    pnt CA;
    pnt C;
};


/*
 * TODO delete this debugging function
 */
void ppnt(pnt p){
    printf("%f, %f, %f\n", p.x, p.y, p.x); 
}

/*
 * Calculate distance between two points
 */
float dist(pnt a, pnt b)
{
    return sqrt((a.x - b.x) * (a.x - b.x) + 
                (a.y - b.y) * (a.y - b.y) +
                (a.z - b.z) * (a.z - b.z));
}

/*
 * Solve for the angle between points a, b, and c
 */
float angle(pnt a, pnt b, pnt c)
{
    float ab = dist(a, b);
    float bc = dist(b, c);
    float ac = dist(a, c);
    float theta = acos((ab * ab + bc * bc - ac * ac ) / (2 * bc * ac));
    return theta * (180 / PI);
}


/*
 * Adds two vectors
 */
pnt vadd(pnt a, pnt b){
    pnt out;
    out.x = a.x + b.x;
    out.y = a.y + b.y;
    out.z = a.z + b.z;
    return out;
}

/*
 * Subtracts two vectors
 */
pnt vsub(pnt a, pnt b){
    pnt out;
    out.x = a.x - b.x;
    out.y = a.y - b.y;
    out.z = a.z - b.z;
    return out;
}

/*
 * Multiply a point by a scalar
 */
pnt vmult(pnt a, float x){
    pnt out;
    out.x = a.x * x;
    out.x = a.y * x;
    out.x = a.z * x;
    return out;
}

/*
 * Calculate cross-product of two vectors
 */
pnt cross(pnt a, pnt b){
    pnt out;
    out.x = a.y * b.z - a.z * b.y;
    out.y = a.z * b.x - a.x * b.z;
    out.z = a.x * b.y - a.y * b.x;
    return out;
}

/*
 * Calculate ||v||
 */
float magnitude(pnt a){
    return float(pow(a.x * a.x + a.y * a.y + a.z * a.z, 0.5));
}

/*
 * Calculates the dot product of two vectors
 */
float dot(pnt a, pnt b){
    return a.x * b.x +
           a.y * b.y +
           a.z * b.z;
}

/*
 * Calculate the angle between planes abc and bcd
 */
float torsion_angle(pnt a, pnt b, pnt c, pnt d)
{
    // Convert to internal coordinates
    pnt ab = vsub(b, a);
    pnt ac = vsub(c, a);
    pnt dc = vsub(c, d);
    pnt db = vsub(b, d);

    // Get vecors normal to each plane
    pnt n1 = cross(ab, ac);
    pnt n2 = cross(dc, db);

    // Calculate torsion angle
    float t_angle = acos( dot(n1, n2) / (magnitude(n1) * magnitude(n2)) );

    // Calculate sign
    int sign = sin(dot(n1, dc)) > 0 ? 1 : -1;

    return t_angle * (180 / PI) * sign - sign * 180;
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
            j++;
        }
    }
    return backbone;
}


int print_backbone_statistics(vector<struct Peptide> b){
    printf(
        "%s %s %s %s %s %s %s %s %s %s %s\n",
        "ind", "aa", "N-CA", "CA-C", "C-N", "N-CA-C", "CA-C-N+", "C-N+-CA+", "phi", "psi", "omega"
    );
    float psi, phi;
    for(int i = 0; i < b.size(); i++){
        bool not_last=(i + 1) < b.size();
        printf(
            "%d %s   %f %f %f   %f %f %f   %f %f %f\n",
            // Peptide number
            i + 1,
            // residue three letter name
            b[i].residue,

            // N -> CA bond length
            dist(b[i].N, b[i].CA),
            // CA -> C bond length
            dist(b[i].CA, b[i].C),
            // C -> N bond length
            (i+1) < b.size() ? dist(b[i].C, b[i+1].N) : 999,

            // N-CA-C angle
            angle(b[i].N, b[i].CA, b[i].C),
            // CA-C-N+ angle
            not_last ? angle(b[i].CA, b[i].C, b[i+1].N) : (float)999,
            // C-N+-CA+ angle
            not_last ? angle(b[i].C, b[i+1].N, b[i+1].CA) : (float)999,
            
            // phi torsion angle - C- _ N _ CA _ C
            i > 0 ? torsion_angle(b[i-1].C, b[i].N, b[i].CA, b[i].C) : (float)999,
            // psi torsion angle - N _ CA _ C _ N+
            not_last ? torsion_angle(b[i].N, b[i].CA, b[i].C, b[i+1].N) : (float)999,
            // omega torsion angle - CA _ C _ N+ _ CA+
            not_last ? torsion_angle(b[i].CA, b[i].C, b[i+1].N, b[i+1].CA) : (float)999
        );
    }
    return 1;
}


vector<struct Atom> load_pdb_file()
{
    vector<Atom> atoms;
    string line;
    while(cin){
        getline(cin, line);
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


int main(int argc, char* argv[])
{
    vector<struct Atom> atoms;
    // open a file, print its contents
    if(argc == 1){
        atoms = load_pdb_file(); 
    } else {
        fprintf(stderr, "Please provide input (from STDIN)\n");
        return 1;
    }
    vector<struct Peptide> backbone = get_backbone(atoms);
    print_backbone_statistics(backbone);
    return 0;
}
