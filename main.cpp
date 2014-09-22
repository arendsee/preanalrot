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
        "%s %s %s %s %s %s %s %s %s %s %s %s\n",
        "ind", "aa", "N-CA", "CA-C", "C-N", "CA-CA", "N-CA-C", "CA-C-N+", "C-N+-CA+", "phi", "psi", "omega"
    );
    double psi, phi;
    for(int i = 0; i < b.size(); i++){
        bool not_last=(i + 1) < b.size();
        printf(
            "%d %s %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n",
            // Peptide number
            i + 1,
            // residue three letter name
            b[i].residue,

            // N -> CA bond length
            dist(b[i].N, b[i].CA),
            // CA -> C bond length
            dist(b[i].CA, b[i].C),
            // C -> N bond length
            not_last ? dist(b[i].C, b[i+1].N) : 999,
            // CA - CA bond length
            not_last ? dist(b[i].CA, b[i+1].CA) : 999,

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


/*
 * Rotate point p theta degrees about the axis formed by points a and c
 */
pnt rotate(pnt p, pnt a, pnt b, double theta){
    double t = theta * ( PI / 180 );
    pnt u = vsub(b, a);
    u = vmult(u, 1 / magnitude(u));
    
    p = vsub(p, a);
    
    double ct = cos(t);
    double st = sin(t);
    
    pnt q;
    q.x = p.x * (ct + (1 - ct) * u.x * u.x) +
          p.y * ((1 - ct) * u.x * u.y - u.z * st) +
          p.z * ((1 - ct) * u.x * u.z + u.y * st);
    
    q.y = p.x * ((1 - ct) * u.x * u.y + u.z * st) +
          p.y * (ct + (1 - ct) * u.y * u.y) +
          p.z * ((1 - ct) * u.y * u.z - u.x * st);
    
    q.z = p.x * ((1 - ct) * u.x * u.z - u.y * st) +
          p.y * ((1 - ct) * u.y * u.z + u.x * st) +
          p.z * (ct + (1 - ct) * u.z * u.z);
    
    q = vadd(q, a);

    return(q);
}

/*
 * Handle rotation logic and print output
 */
bool print_rotated(int residue, char * bond, double theta){

    bool is_psi = strcmp(bond, "psi") == 0 ? true : false;
    bool is_phi = strcmp(bond, "phi") == 0 ? true : false;
    if(! is_psi && ! is_phi){
        fprintf(stderr, "I can only rotate the psi and phi angles\n"); 
        return false;
    }

    string line;
    char atom_buffer[5];
    int rid = 0;
    pnt axis_1, axis_2, c;
    bool set = false; 

    string first = is_psi ? "CA" : "N";
    string second = is_psi ? "C" : "CA";

    while(cin){
        getline(cin, line);
        if(line.substr(0,4) == "ATOM"){
            pnt p;
            sscanf(line.c_str(),
                   "%*s %*d %s %*s %*s %d %lf %lf %lf",
                   atom_buffer, &rid, &p.x, &p.y, &p.z);
            string atom = string(atom_buffer);

            if(rid == residue && atom == first){
                axis_1 = p;
            }
            else if(rid == residue && atom == second){
                axis_2 = p;
                set = true;
            }
            // If the bond or rotation has been found, rotate
            // UNLESS, it is a psi bond AND we are in the bond's peptide
            else if(set){
                if(is_phi || (is_psi && (rid != residue || atom == "O"))){
                    pnt r = rotate(p, axis_1, axis_2, theta);
                    printf("%s%8.3f%8.3f%8.3f %s\n",
                           line.substr(0,30).c_str(),
                           r.x, r.y, r.z,
                           line.substr(55,26).c_str()
                          );
                    continue;
                }
            }
        }
        printf("%s\n", line.c_str());
    }
    return true;
}

    // pnt pos;
    // char residue[4];
    // char element[3];
    // char atom_name[5];
    // int serial_id;
    // int aa_id;

bool print_chi(vector<Atom> atom){
    vector<Atom> four;
    string name;
    int nchi = 0;
    for(int i = 0; i < atom.size(); i++){
        name = atom[i].atom_name; 
        bool poop = true;
        if(name == "N"){
            four.clear();
            four.push_back(atom[i]);
            poop = false;
        } 
        else if(name == "CA" || name == "CB"){
            four.push_back(atom[i]);
            poop = false;
        }
        // Chi1
        else if(name == "CG" || name == "SG" || name == "CG1" || name == "OG" || name == "OG1"){
            four.push_back(atom[i]);
            nchi = 1;
        }
        // Chi2
        else if(name == "CD" || name == "OD1" || name == "ND1" || name == "SD"){
            four.push_back(atom[i]);
            nchi = 2;
        }
        // Chi3
        else if(name == "NE" || name == "OE1" || name == "CE"){
            four.push_back(atom[i]);
            nchi = 3;
        }
        // Chi4
        else if(name == "CZ" || name == "NZ"){
            four.push_back(atom[i]);
            nchi = 4;
        }
        // Chi5
        else if(name == "Nh1"){
            four.push_back(atom[i]);
            nchi = 5;
        }
        else{
            poop = false;
        }

        if(poop){
            int s = four.size();
            printf("%s-%d %s-%s-%s-%s chi%d %f\n",
                   atom[i].residue,
                   atom[i].aa_id,
                   four[s-4].atom_name,
                   four[s-3].atom_name,
                   four[s-2].atom_name,
                   four[s-1].atom_name,
                   nchi,
                   torsion_angle(four[s-4].pos,
                                 four[s-3].pos,
                                 four[s-2].pos,
                                 four[s-1].pos)
                );
            poop = false;
        }
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


void print_usage(){
    fprintf(stderr, "Usage: subcommand [arguments] < STDIN\n");
    fprintf(stderr, "./a.out bstat < in.pdb\n");
    fprintf(stderr, "./a.out rotate [residue number] [phi|psi] [angle] < in.pdb\n");
    fprintf(stderr, "./a.out mindist < in.pdb\n");
    fprintf(stderr, "./a.out chi < in.pdb\n");
}

int main(int argc, char* argv[])
{
    vector<struct Atom> atoms;
    string subcommand;

    if(argc > 1){
        subcommand = argv[1];
    } else {
        print_usage();
        return 0;
    }

    // Prints the bond lengths and backbone angles
    if(subcommand == "bstat"){
        atoms = load_pdb_file(); 
        vector<struct Peptide> backbone = get_backbone(atoms);
        print_backbone_statistics(backbone);
    } 
    // Rotates about the psi or phi bond of the given residue
    else if(subcommand == "rotate"){ 
        if(argc == 5){
            if(! print_rotated(atoi(argv[2]), argv[3], atof(argv[4])) ){
                return 1;
            }
        }
        else {
            fprintf(stderr, "'rotate' requires 3 arguments:\n");
            fprintf(stderr, "1) residue number\n2) bond['psi' or 'phi']\n3) rotation angle\n");
            return 1;
        }
    }
    // Calculates chi torsion angles for all side chains
    else if(subcommand == "chi"){ 
        atoms = load_pdb_file(); 
        print_chi(atoms);
    }
    // Calculates each atom's nearest neighbor
    else if(subcommand == "mindist"){ 
        atoms = load_pdb_file(); 
        print_mindist(atoms);
    }
    // If you want help ...
    else if(subcommand == "-h" || subcommand == "-help" || subcommand == "-help"){
        print_usage();
    }
    else {
        print_usage();
    }
    return 0;
}
