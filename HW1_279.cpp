#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <iomanip>
using namespace std;

class Au{
    private:
        // Keeping these variable private ensures that users won't be able to tamper with their values
        const double Epsilon_Au = 5.29; // (kcal/mol)
        const double Sigma_Au = 2.951;  // (Angstrom)
        const double tol = 1e-5; // Tolerance for convergance in line search SD

    public:
    
        // Function to calculate the Lennard Jones potential energy between 2 atoms
        double Calculate_LJ(double distance) const {
            double sig_dist = Sigma_Au / distance; // Calculating the Rij term
            double term6 = pow(sig_dist , 6);
            double term12 = pow(sig_dist , 12); 
            return (Epsilon_Au * (term12 - (2 * term6) )); // Formula for LJ Energy (Au)
        }   


        // Function to calculate the distance between 2 atoms
        double Calculate_Distance(const vector<double> & a1 , const vector<double> & a2) const {
            // Distance is calculated by finding the squared difference between each of the coordinates

            return sqrt(pow( (a2[0] - a1[0]) , 2 ) + 
                        pow( (a2[1] - a1[1]) , 2 ) +
                        pow( (a2[2] - a1[2]) , 2 ) );
        }


        // Function to read in the text files that contain relevant info regarding our atoms
        vector<vector<double>> Read_Atoms(const string & filename) const {
            ifstream file(filename);
            if (!file){
                throw runtime_error("Unable to open file!"); // If specified file isn't opened properly, then program will error out
            }

            int rows;
            file >> rows;
            file.ignore(numeric_limits<streamsize>::max() , '\n');
            vector<vector<double>> data; // Creating a vector to store atomic coordinates

            string line;
            while (getline(file , line)){
                istringstream iss(line);
                vector<double> row;
                double coord; // Initialzing vector that will contain coordinate data

                int atomic_num;
                iss >> atomic_num;
                if(atomic_num != 79){
                    throw runtime_error("Error: Atom other than Gold(Au) is present in the data!"); // If another element is detected, program will error out
                }

                while (iss >> coord){
                    row.push_back(coord); // vector is filled with coordinate data
                }
                data.push_back(row); // vector from above is appended to another vector
            }

            return data; // Returning the nested vectors that hold the relevant atomic data
    }


    // Function that calculates the total energy of a cluster of atoms. Utilizes the 'Calculate_Distance' & 'Calculate_LJ' function from above
    double Calculate_Cluster_Energy(const vector<vector<double>> & atoms){
        double cluster_energy = 0.0; // Initialzing the variable that'll hold the LJ Energy for the entire cluster

        for(size_t i = 0; i < atoms.size(); i++){ // Outer loop should run for as many iterations as there are atoms in 'atoms' vector
            for(size_t j = i + 1; j < atoms.size(); j++){ // Only atoms in which i < j shall be used in calculations
                double distance = Calculate_Distance(atoms[i] , atoms[j]); // Calculating distance for atoms in which i < j
                double LJ = Calculate_LJ(distance); // Calculating LJ Energy for i < j atoms
                cluster_energy += LJ; // Riemann Sum of the LJ Energies to get the whole cluster
            }
        }

        return cluster_energy; // Returning the Riemann Sum of the LJ Energy
    }   


    // Function that calculates the Analytical Force- which is the derivative of the Lennard Jones Energy formula
    double Calculate_Analytical_Force(double distance) const {
        double term1 = -12 * pow(Sigma_Au, 12) / pow(distance, 13); // Derivative of the 1st term from the LJ formula
        double term2 = 12 * pow(Sigma_Au, 6) / pow(distance, 7); // Derivative of the 2nd term from the LJ formula 
        return (-1 * Epsilon_Au * (term1 + term2)); // Analytical force formula for Au
    }

    // Function that calculates Analytical Force for an entire cluster of atoms
    vector<vector<double>> Calculate_Analytical_Force_Cluster(const vector<vector<double>> & atoms) const {
        vector<vector<double>> forces(atoms.size(), vector<double>(3, 0.0));

        // For-loop structure is very similar to one used to calculate LJ Energy
        for(size_t i = 0; i < atoms.size(); i++){ // Outer loop should run for as many iterations as there are atoms in 'atoms' vector
            for(size_t j = i + 1; j < atoms.size(); j++){
                double dist = Calculate_Distance(atoms[j] , atoms[i]); // Like with LJ Energy, distance is needed before further Energy calculation
                double force_magnitude = Calculate_Analytical_Force(dist);
                
                // Calculating the difference in positions x,y,z between i < j atoms
                double dx = atoms[j][0] - atoms[i][0];
                double dy = atoms[j][1] - atoms[i][1];
                double dz = atoms[j][2] - atoms[i][2];

                // Normalizing the distance to avoid the possibility of dividing by zero 
                // By taking inverse of the distance, we properly scale the directional force when
                // calculating the force component
                double norm_dist = 1.0 / dist;

                // Calculating force components in all the directions 
                double fx = force_magnitude * dx * norm_dist;
                double fy = force_magnitude * dy * norm_dist;
                double fz = force_magnitude * dz * norm_dist;

                
                // Forces on atoms i & j are opposite but equal (Newton's 3rd Law)
                forces[i][0] -= fx;
                forces[i][1] -= fy;
                forces[i][2] -= fz;

                forces[j][0] += fx;
                forces[j][1] += fy;
                forces[j][2] += fz;
            }
            

        }
        return forces; // Returning the nested vector that contains info regarding forces acting on each atom in each direction
    }


    // Function that calculates the forward difference
    double Calculate_Foward_Difference(const vector<vector<double>> & atoms, double h, int axis, int atom_index){
        vector<vector<double>> shifted_atoms = atoms; // Creating a copy of the atoms 

        shifted_atoms[atom_index][axis] += h; // Shifting the atom's position by +h. 'axis' denotes x,y,z (x=0, y=1, z=2)

        double shifted_energy = Calculate_Cluster_Energy(shifted_atoms); // Calculating LJ energy for forward shited atoms
        double original_energy = Calculate_Cluster_Energy(atoms); // Calculating LJ energy for original position atoms

        return -(shifted_energy - original_energy)/h; // Formula for Forward Difference
    }


    // Function that calculates the forward difference on each atom in a cluster of Au atoms
    vector<vector<double>> Calculate_Forward_Difference_Cluster(const vector<vector<double>> & atoms, double h){
        vector<vector<double>> forward_forces(atoms.size(), vector<double>(3, 0.0));

        for(size_t i = 0; i < atoms.size(); i++){ // Outer loop should run for as many iterations as there are atoms in 'atoms' vector
            for(int axis = 0; axis < 3; axis++){ // Loop is set to 3 iterations since there are 3 axis' possible (x,y,z)
                forward_forces[i][axis] = Calculate_Foward_Difference(atoms, h, axis, i);
            }
        }
        return forward_forces; // Returns nested vector that contains forward forces acting on each atom on each axis plane
    }


    // Function that calculates the central difference
    double Calculate_Central_Difference(const vector<vector<double>> & atoms, double h, int axis, int atom_index){

        // Making a copy of input vector so it can be shifted in later code
        vector<vector<double>> plus_h_atoms = atoms;
        vector<vector<double>> minus_h_atoms = atoms;

        plus_h_atoms[atom_index][axis] += h; // Shifting the atom's position by +h. 'axis' denotes x,y,z (x=0, y=1, z=2)
        minus_h_atoms[atom_index][axis] -= h; // Shifting the atom's position by -h. 'axis' denotes x,y,z (x=0, y=1, z=2)

        double plus_energy = Calculate_Cluster_Energy(plus_h_atoms); // Calculating LJ energy for forward shited atoms
        double minus_energy = Calculate_Cluster_Energy(minus_h_atoms); // Calculating LJ energy for backward shited atoms

        return -((plus_energy - minus_energy) / (2 * h)); // Formula for calculating central difference
    }


    // Function that calculates the central difference on each atom in a cluster of Au atoms
    vector<vector<double>> Calculate_Central_Difference_Cluster(const vector<vector<double>> & atoms, double h){
        vector<vector<double>> central_forces(atoms.size(), vector<double>(3, 0.0));

        for(size_t i = 0; i < atoms.size(); i++){ // Outer loop should run for as many iterations as there are atoms in 'atoms' vector
            for(int axis = 0; axis < 3; axis++){ // Loop is set to 3 iterations since there are 3 axis' possible (x,y,z)
                central_forces[i][axis] = Calculate_Central_Difference(atoms, h, axis, i); // Calculating Central Difference
            }
        }
        return central_forces; // Returns nested vector that contains central forces acting on each atom on each axis plane
    }


    // Function that formats the forces that were calculated and outputs them to console
    void Print_Forces(const vector<vector<double>> & forces){
        for(int i = 0; i < 3; i++){ // Loop is set to 3 iterations since there are 3 axis' possible (x,y,z)
            for(size_t j = 0; j < forces.size(); j++){
                cout << scientific << setprecision(4) << setw(10)
                    << "" << forces[j][i];
            }
            cout << endl;
        }
    }

    // Function that performs Steepest Descent with Line Search
    void Steepest_Descent_With_Line_Search(vector<vector<double>> &atoms, double step_size, double conv_threshold){
        double current_energy = Calculate_Cluster_Energy(atoms); // Calculating LJ Energy for the cluster
        vector<vector<double>> forces = Calculate_Central_Difference_Cluster(atoms, 0.001); // Using central difference for force

        int iterations = 0; // Variable to count the number of iterations our SD will take to converge
        while(true){ // Starting an infinite loop that will run until convergance is reached
            
            // Norm of the forces to check convergence
            double force_norm = 0.0;
            for (auto atom_force : forces){
                for (auto f : atom_force){
                    force_norm += f * f; // Squaring each component of force for each atom and adding it all up
                }
            }
            force_norm = sqrt(force_norm); // Getting total norm requires us to take square root of the total

            if (force_norm < conv_threshold) { // If the force_norm is less than the threshold, loop stops
                break;
            }

            // Apply golden section line search to find optimal step size
            double alpha = Golden_Section_Search(atoms, forces, step_size);
            cout << "new_point" << endl;

            for (const auto &atom : atoms) {
                cout << "79(" << fixed << setprecision(4) << atom[0] << ", " 
                << atom[1] << ", " << atom[2] << ")" << endl;
            }


            // Nested for-loops to update the positions based on forces & alpha value
            for (size_t i = 0; i < atoms.size(); i++){ // Outer loop should run for as many iterations as there are atoms in 'atoms' vector
                for (size_t j = 0; j < 3; j++){ // Loop is set to 3 iterations since there are 3 axis' possible (x,y,z)
                    atoms[i][j] += alpha * forces[i][j];
                }
            }

            // Recalculate the energy and forces
            current_energy = Calculate_Cluster_Energy(atoms);
            forces = Calculate_Central_Difference_Cluster(atoms, 0.001);

            cout << "current energy: " << scientific << setprecision(4) << current_energy << endl;
            cout << "Central Difference Force" << endl;
            Print_Forces(forces);

            iterations++; // Increasing iterations each time we reach the end of the loop
        }

        // Information to be printed at the end, once convergance has been reached
        cout << "Total iterations: " << iterations << endl;
        cout << "Final energy: " << scientific << setprecision(4) << current_energy << endl;
        cout << "Optimized structure:" << endl;

        // Outputting the final positions of our Gold atom cluster
        for (const auto &atom : atoms) {
            cout << "79(" << fixed << setprecision(4) << atom[0] << ", " 
            << atom[1] << ", " << atom[2] << ")" << endl;
        }

    }


    double Golden_Section_Search(const vector<vector<double>> &atoms, const vector<vector<double>> &forces, double step_size){
        double a = 0.0;
        double b = step_size;
        double gr = (1 + sqrt(5)) / 2; // golden ratio

        // Below 2 lines calculate 2 new points within the interval of interest
        double c = b - (b - a) / gr;
        double d = a + (b - a) / gr;

        vector<vector<double>> temp_atoms = atoms; // Making a copy of input atoms in order to run series of calculations later
        
        // Updating the position and calculating the LJ Energy of the updated cluster
        Update_Positions(temp_atoms, forces, c);
        double energy_c = Calculate_Cluster_Energy(temp_atoms);

        temp_atoms = atoms;

        // Updating the position and calculating the LJ Energy of the updated cluster
        Update_Positions(temp_atoms, forces, d);
        double energy_d = Calculate_Cluster_Energy(temp_atoms);

        // As long as difference between c & d is greater than the tolerance, this block will always run
        while (fabs(c - d) > tol) { // Choose fabs() as we are dealing with decimal places in these calculations

            // If  c is lower than d, the minimum likely lies between a and d.
            if (energy_c < energy_d) {
                b = d;
            } else {
                a = c; // Otherwise lower bound a is set to c
            }

            // Recalculating c and d after adjusting the interval
            c = b - (b - a) / gr;
            d = a + (b - a) / gr;

            // Making a copy of atoms, updating the positions, and finding its cluster energy
            temp_atoms = atoms;
            Update_Positions(temp_atoms, forces, c);
            energy_c = Calculate_Cluster_Energy(temp_atoms);

            temp_atoms = atoms;
            Update_Positions(temp_atoms, forces, d);
            energy_d = Calculate_Cluster_Energy(temp_atoms);
        }

        return (b + a) / 2; // Return the midpoint of the final interval
    }


    // Function that updates the coordinates of our atoms based on results from SD function
    void Update_Positions(vector<vector<double>> &atoms, const vector<vector<double>> &forces, double step_size){
        for(size_t i = 0; i < atoms.size(); i++){ // Outer loop should run for as many iterations as there are atoms in 'atoms' vector
            for(size_t j = 0; j < 3; j++){ // Loop is set to 3 iterations since there are 3 axis' possible (x,y,z)
                atoms[i][j] += step_size * forces[i][j]; // New position is updated to current position plus product of step size & force values
            }
        }
    }

        
};




int main(){
    // Using 1.txt from 'Energy/'
    Au atom1;
    vector<vector<double>> coords = atom1.Read_Atoms("sample_input/Energy/1.txt");
    for(auto row : coords){
        cout << "79(" << row[0] << ", " << row[1] << ", " << row[2] << ")" << endl; 
    }
    double total_energy = atom1.Calculate_Cluster_Energy(coords);
    cout << "E_LJ = " << total_energy << endl;
    cout << endl << endl;


    // Using 2.txt from 'Energy/'
    Au atom2;
    vector<vector<double>> coords2 = atom2.Read_Atoms("sample_input/Energy/2.txt");
    for(auto row : coords2){
        cout << "79(" << row[0] << ", " << row[1] << ", " << row[2] << ")" << endl;
    }
    double total = atom2.Calculate_Cluster_Energy(coords2);
    cout << "E_LJ = " << total << endl;
    cout << endl << endl;


    // Using 3.txt from 'Energy/'
    Au atom3;
    vector<vector<double>> coords3 = atom3.Read_Atoms("sample_input/Energy/3.txt");
    for(auto row : coords3){
        cout << "79(" << row[0] << ", " << row[1] << ", " << row[2] << ")" << endl;
    }
    double t = atom3.Calculate_Cluster_Energy(coords3);
    cout << "E_LJ = " << t << endl;
    cout << endl << endl;


    // Using my_testcase.txt from 'Energy/'
    Au atom4;
    vector<vector<double>> coords_my_testcase = atom3.Read_Atoms("sample_input/Energy/my_testcase.txt");
    for(auto row : coords_my_testcase){
        cout << "79(" << row[0] << ", " << row[1] << ", " << row[2] << ")" << endl;
    }
    double my_t_E = atom4.Calculate_Cluster_Energy(coords_my_testcase);
    cout << "E_LJ = " << my_t_E << endl;
    cout << endl << endl;


    // // Using 1.txt from '/Force'
    Au atom1A;
    vector<vector<double>> coords1A = atom1A.Read_Atoms("sample_input/Force/1.txt");

    double total_energy1A = atom1A.Calculate_Cluster_Energy(coords1A);
    cout << "E_LJ = " << total_energy1A << endl;
    vector<vector<double>> forces = atom1A.Calculate_Analytical_Force_Cluster(coords1A);

    cout << "F_LJ_analytical" << endl;
    atom1A.Print_Forces(forces);




    cout << "Stepsize for finite difference:0.1" << endl;

    cout << "F_LJ_forward_difference" << endl;
    vector<vector<double>> forward_force = atom1A.Calculate_Forward_Difference_Cluster(coords1A, 0.1);
    atom1A.Print_Forces(forward_force);

    vector<vector<double>> central_force = atom1A.Calculate_Central_Difference_Cluster(coords1A, 0.1);
    cout << "F_LJ_central_difference" << endl;
    atom1A.Print_Forces(central_force);




    cout << "Stepsize for finite difference:0.01" << endl;

    cout << "F_LJ_forward_difference" << endl;
    vector<vector<double>> forward_force2 = atom1A.Calculate_Forward_Difference_Cluster(coords1A, 0.01);
    atom1A.Print_Forces(forward_force2);

    vector<vector<double>> central_force2 = atom1A.Calculate_Central_Difference_Cluster(coords1A, 0.01);
    cout << "F_LJ_central_difference" << endl;
    atom1A.Print_Forces(central_force2);




    cout << "Stepsize for finite difference:0.001" << endl;

    cout << "F_LJ_forward_difference" << endl;
    vector<vector<double>> forward_force3 = atom1A.Calculate_Forward_Difference_Cluster(coords1A, 0.001);
    atom1A.Print_Forces(forward_force3);

    vector<vector<double>> central_force3 = atom1A.Calculate_Central_Difference_Cluster(coords1A, 0.001);
    cout << "F_LJ_central_difference" << endl;
    atom1A.Print_Forces(central_force3);




    cout << "Stepsize for finite difference:0.0001" << endl;

    cout << "F_LJ_forward_difference" << endl;
    vector<vector<double>> forward_force4 = atom1A.Calculate_Forward_Difference_Cluster(coords1A, 0.0001);
    atom1A.Print_Forces(forward_force4);

    vector<vector<double>> central_force4 = atom1A.Calculate_Central_Difference_Cluster(coords1A, 0.0001);
    cout << "F_LJ_central_difference" << endl;
    atom1A.Print_Forces(central_force4);

    cout << endl << endl;


    // Using 1.txt from '/SD_with_line_search'
    Au atom1B;
    vector<vector<double>> coords1B = atom1B.Read_Atoms("sample_input/SD_with_line_search/1.txt");
    double total_E = atom1B.Calculate_Cluster_Energy(coords1B);

    vector<vector<double>> central_force_SD = atom1B.Calculate_Central_Difference_Cluster(coords1B, 0.0001);

    double step_size = 0.3;
    double convergance_threshold = 0.01;

    cout << "start steepest descent with golden section line search" << endl;
    cout << "Initial energy:" << total_E << endl;
    cout << "Stepsize for central difference is:0.0001;Initial stepsize for line search is:0.3;Threshold for convergence in force is:0.01" << endl;
    cout << "Central Difference Force" << endl;
    atom1B.Print_Forces(central_force_SD);

    atom1B.Steepest_Descent_With_Line_Search(coords1B, step_size, convergance_threshold);


    // Using 2.txt from '/SD_with_line_search'
    // Au atom3;
    // vector<vector<double>> coords3 = atom3.Read_Atoms("sample_input/SD_with_line_search/2.txt");
    // double total_E2 = atom3.Calculate_Cluster_Energy(coords3);
    // vector<vector<double>> central_force_SD2 = atom3.Calculate_Central_Difference_Cluster(coords3, 0.0001);

    // cout << "start steepest descent with golden section line search" << endl;
    // cout << "Initial energy:" << total_E2 << endl;
    // cout << "Stepsize for central difference is:0.0001;Initial stepsize for line search is:0.3;Threshold for convergence in force is:0.01" << endl;
    // cout << "Central Difference Force" << endl;
    // atom2.Print_Forces(central_force_SD2);

    // atom2.Steepest_Descent_With_Line_Search(coords3, step_size, convergance_threshold);




    return 0;
}