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
#include <cstdlib>
using namespace std;

const double PI = acos(-1);


struct Point
{
    double x;
    double y;
    double z;
};
typedef struct Point pnt;


struct Atom
{
    pnt pos;
    char residue[4];
    char element[3];
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
    printf("%lf, %lf, %lf\n", p.x, p.y, p.z); 
}

/*
 * Calculate distance between two points
 */
double dist(pnt a, pnt b)
{
    return sqrt((a.x - b.x) * (a.x - b.x) + 
                (a.y - b.y) * (a.y - b.y) +
                (a.z - b.z) * (a.z - b.z));
}

/*
 * Solve for the angle between points a, b, and c
 */
double angle(pnt a, pnt b, pnt c)
{
    double ab = dist(a, b);
    double bc = dist(b, c);
    double ac = dist(a, c);
    double theta = acos((ab * ab + bc * bc - ac * ac ) / (2 * bc * ac));
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
pnt vmult(pnt a, double x){
    pnt out;
    out.x = a.x * x;
    out.y = a.y * x;
    out.z = a.z * x;
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
double magnitude(pnt a){
    return double(pow(a.x * a.x + a.y * a.y + a.z * a.z, 0.5));
}

/*
 * Calculates the dot product of two vectors
 */
double dot(pnt a, pnt b){
    return a.x * b.x +
           a.y * b.y +
           a.z * b.z;
}

/*
 * Calculate the angle between planes abc and bcd
 */
double torsion_angle(pnt a, pnt b, pnt c, pnt d)
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
    double t_angle = acos( dot(n1, n2) / (magnitude(n1) * magnitude(n2)) );

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
    double psi, phi;
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
            not_last ? angle(b[i].CA, b[i].C, b[i+1].N) : (double)999,
            // C-N+-CA+ angle
            not_last ? angle(b[i].C, b[i+1].N, b[i+1].CA) : (double)999,
            
            // phi torsion angle - C- _ N _ CA _ C
            i > 0 ? torsion_angle(b[i-1].C, b[i].N, b[i].CA, b[i].C) : (double)999,
            // psi torsion angle - N _ CA _ C _ N+
            not_last ? torsion_angle(b[i].N, b[i].CA, b[i].C, b[i+1].N) : (double)999,
            // omega torsion angle - CA _ C _ N+ _ CA+
            not_last ? torsion_angle(b[i].CA, b[i].C, b[i+1].N, b[i+1].CA) : (double)999
        );
    }
    return 1;
}

/*
 * Get van der Waals radius of given atom (in angstroms)
 *
 * TODO I should extend this to handle other elements they may be in ligands,
 * e.g. Fe
 */
double get_radius(struct Atom a){
    double radius;
    if(strcmp(a.element, "H") == 0){
        radius = 1.2; 
    }
    else if(strcmp(a.element, "C") == 0){
        radius = 1.7; 
    }
    else if(strcmp(a.element, "O") == 0){
        radius = 1.52; 
    }
    else if(strcmp(a.element, "N") == 0){
        radius = 1.55; 
    }
    else if(strcmp(a.element, "S") == 0){
        radius = 1.8; 
    }
    else{
        radius = 1.6;
    }
    return radius;
}


/*
 * Atoms are considered adjacent if
 * 1) They are the same atom
 * 2) They are the C and N joining two peptides
 * 3) They are in the same peptide (exception: the carbonyl oxygen is only
 *    considered adjacent to the carbonyl carbon)
 */
bool are_adjacent(struct Atom a, struct Atom b){
    // If they are identical, they are considered adjacent
    if(a.serial_id == b.serial_id){
        return true;
    }


    // Order inputs be serial id
    if(a.serial_id > b.serial_id){
        struct Atom temp = a;
        a = b;
        b = temp;
    }

    // If atoms are the C and N that join peptides, they are adjacent
    if((strcmp(a.atom_name, "C") && strcmp(b.atom_name, "N") && (b.serial_id == a.serial_id + 1))){
        return true;
    }

    // If they are in the same peptide, they are (naively) considered adjacent
    if(a.aa_id == b.aa_id){
        // Unless one of the atoms is the carbonyl oxygen (and the other isn't
        // the carbonyl carbon)
        if(strcmp(a.atom_name, "O") && ! strcmp(b.atom_name, "C")){
            false;
        }
        return true;
    }

    return false;
}


/*
 * Finds the nearest non-adjacent neighbor (see are_adjacent documentation)
 * Prints the following columns:
 * 1) Atom-1 name (e.g. CA)
 * 2) Residue-1 (e.g. MET-1)
 * 3) Atom-2 name - the name of the nearest non-adjacent atom
 * 4) Residue-2
 * 5) distance between the two in angstroms
 * 6) ratio of the summed van der Waals radii to the distance
 *    (radius(a) + radius(b)) / distance
 */
void print_mindist(vector<struct Atom> a){
    double d;
    double vradius;
    for(int i = 0; i < a.size(); i++){
        struct Atom ma;
        double min = 9999;
        for(int j = 0; j < a.size(); j++){
            if (are_adjacent(a[i], a[j])){
                continue;
            }
            d = dist(a[i].pos, a[j].pos);
            if(d < min){
                min = d;
                ma = a[j];
                vradius = get_radius(a[i]) + get_radius(ma);
            }
        }
        printf("%s\t%s-%d\t%s\t%s-%d\t%f\t%f\n",
                a[i].atom_name,
                a[i].residue,
                a[i].aa_id,
                ma.atom_name,
                ma.residue,
                ma.aa_id,
                min,
                min / (get_radius(a[i]) + get_radius(ma))
              );
    }
}


pnt rotate(pnt inpnt, pnt axis_1, pnt axis_2, double angle){
    pnt u = vsub(axis_2, axis_1);
    u = vmult(u, 1 / magnitude(u));
    pnt p = vsub(inpnt, axis_1);

    double t = angle * ( PI / 180 );
    double ct = cos(angle);
    double st = sin(angle);

    pnt r;

    r.x = p.x * (ct + u.x * u.x * (1 - ct)) +
          p.y * (u.x * u.y * (1 - ct) - u.z * st) +
          p.z * (u.x * u.z * (1 - ct) + u.y * st);

    r.y = p.x * (u.y * u.x * (1 - ct) + u.z * st) +
          p.y * (ct + u.y * u.y * (1 - ct)) +
          p.z * (u.y * u.z * (1 - ct) - u.x * st);

    r.z = p.x * (u.z * u.x * (1 - ct) - u.y * st) +
          p.y * (u.z * u.y * (1 - ct) + u.x * st) +
          p.z * (ct + u.z * u.z * (1 - ct));
   
    pnt rotated = vadd(r, axis_1);

    return rotated;
}

void print_rotated(int id1, int id2, double angle){
    string line;
    int sid = 0;
    pnt axis_1, axis_2, c;
    while(cin){
        getline(cin, line);
        if(line.substr(0,4) == "ATOM"){
            pnt p;
            sscanf(line.c_str(),
                   "%*s %d %*s %*s %*s %*d %lf %lf %lf",
                   &sid, &p.x, &p.y, &p.z);
            if(sid == id1){
                axis_1 = p;
            }
            else if(sid == id2){
                axis_2 = p;
            }
            else if(sid > id2){
                pnt r = rotate(p, axis_1, axis_2, angle);

                printf("%s%8.3f%8.3f%8.3f %s\n",
                       line.substr(0,30).c_str(),
                       r.x, r.y, r.z,
                       line.substr(55,26).c_str()
                      );
                continue;
            }
        }
        printf("%s\n", line.c_str());
    }
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
               "%*s %d %5s %3s %*s %d %lf %lf %lf %*lf %*lf %3s",
               &atom.serial_id, 
               atom.atom_name,
               atom.residue,
               &atom.aa_id,
               &atom.pos.x, &atom.pos.y, &atom.pos.z,
               atom.element
              );
        atoms.push_back(atom);
    }
    return atoms;
}


int main(int argc, char* argv[])
{
    vector<struct Atom> atoms;
    // open a file, print its contents
    if(argc > 1 && strcmp(argv[1], "bstat") == 0){
        atoms = load_pdb_file(); 
        vector<struct Peptide> backbone = get_backbone(atoms);
        print_backbone_statistics(backbone);
    } 
    else if(argc > 1 && strcmp(argv[1], "rotate") == 0){ 
        if(! argc == 5){
            fprintf(stderr, "'rotate' requires 3 arguments: serial_id-1 serial_id-2 angle\n");
            return 1;
        }
        print_rotated(atoi(argv[2]), atoi(argv[3]), atof(argv[4]));
    }
    else if(argc > 1 && strcmp(argv[1], "chi") == 0){ 
        fprintf(stderr, "not implemented\n");
        return 1;
    }
    else if(argc > 1 && strcmp(argv[1], "mindist") == 0){ 
        atoms = load_pdb_file(); 
        print_mindist(atoms);
    }
    else {
        fprintf(stderr, "Please provide an option [bstat, rotate, chi, mindist]\n");
        return 1;
    }
    return 0;
}
