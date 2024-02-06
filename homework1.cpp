#include <vector>
#include <array>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <stdio.h>

using namespace std;

#define X = 0;
#define Y = 1;
#define Z = 2;

class Lennard_Jones
{

public:
    Lennard_Jones(int atom, double epsilon, double sigma);
    ~Lennard_Jones();
    vector<vector<double>> read_xyz(string file_path);
    double calculate_distance(vector<double> coord_1, vector<double> coord_2);
    double calculate_energy(vector<double> atom1, vector<double> atom2);
    double total_energy(vector<vector<double>> coords);
    void output_file(string file_name);
    void output_force_file(string file_name);
    void print_output_force();
    void print_output();
    double calculate_analytic_force(int component, int atom_index, vector<vector<double>> coords);
    double calculate_forward_difference_force(int component, int atom_index, double step_size, vector<vector<double>> coords);
    double calculate_finite_difference_force(int component, int atom_index, double step_size, vector<vector<double>> coords);
    vector<vector<double>> return_analytic_force(vector<vector<double>> coords);
    vector<vector<double>> return_forward_difference_force(vector<vector<double>> coords, double stepsize);
    vector<vector<double>> return_central_difference_force(vector<vector<double>> coords, double stepsize);
    vector<vector<vector<double>>> optimize(vector<vector<double>> coords, double threshold);

private:
    double epsilon;
    double sigma;
    int atomic_num;
    vector<vector<double>> atom_coords;
    vector<int> elements;
    int num_of_atoms;
    int error;
};

Lennard_Jones::~Lennard_Jones()
{
}

Lennard_Jones::Lennard_Jones(int atom, double binding_energy, double eq_distance)
{
    atomic_num = atom;
    epsilon = binding_energy;
    sigma = eq_distance;
    error = 0;
}

vector<vector<double>> Lennard_Jones::read_xyz(string file_path)
{
    ifstream myfile(file_path);
    double current;
    int atom;

    if (myfile.is_open())
    {

        myfile >> num_of_atoms;
        while (myfile >> atom)
        {
            vector<double> coords;
            elements.push_back(atom);
            if (atom != atomic_num)
            {

                cout << "All elements must have atomic number " << atomic_num << endl;
                error = -1;
            }

            for (int i = 0; i < 3; i++)
            {
                if (myfile >> current)
                {
                    // cout << current << " "; uncomment if you want them printed to command line each time the function is ran
                    coords.push_back(current);
                }
            }

            if (coords.size() == 3)
            {
                atom_coords.push_back(coords);
            }

            // cout << endl;
        }

        return atom_coords;
    }
}

double Lennard_Jones::calculate_distance(vector<double> coord_1, vector<double> coord_2)
{

    return sqrt(pow(coord_1[0] - coord_2[0], 2) + pow(coord_1[1] - coord_2[1], 2) + pow(coord_1[2] - coord_2[2], 2));
}

double Lennard_Jones::calculate_energy(vector<double> atom1, vector<double> atom2)
{
    if (!(error < 0))
    {
        double energy = 0;
        double distance;
        distance = calculate_distance(atom1, atom2);
        energy = epsilon * (pow(sigma / distance, 12) - 2 * pow(sigma / distance, 6));
        return energy;
    }
}

double Lennard_Jones::total_energy(vector<vector<double>> coords)
{
    double energy = 0;
    if (!(error < 0))
    {

        for (int i = 0; i < coords.size(); i++)
        {
            for (int j = i + 1; j < coords.size(); j++)
            {
                energy += calculate_energy(coords[i], coords[j]);
            }
        }

        return energy;
    }
}
// void Lennard_Jones::output_file(string file_name)
// {

//     ofstream myfile(file_name);
//     if (myfile.is_open())
//     {
//         for (int i = 0; i < num_of_atoms; i++)
//         {
//             myfile << elements[i] << " (" << atom_coords[i][0] << ", " << atom_coords[i][1]
//                    << ", " << atom_coords[i][2] << ")" << endl;
//         }
//         myfile << "E_LJ = " << total_energy() << endl;
//     }
// }

// void Lennard_Jones::print_output()
// {
//     for (int i = 0; i < num_of_atoms; i++)
//     {
//         cout << elements[i] << " (" << atom_coords[i][0] << ", " << atom_coords[i][1]
//              << ", " << atom_coords[i][2] << ")" << endl;
//     }
//     cout << "E_LJ = " << total_energy() << endl;
// }

double Lennard_Jones::calculate_analytic_force(int component, int atom_index, vector<vector<double>> coords)
{
    if (!(error < 0))
    {
        int k = atom_index;
        double force = 0;
        double distance;

        for (int i = 0; i < coords.size(); i++)
        {
            if (i != k)
            {

                distance = calculate_distance(coords[k], coords[i]);

                force += (epsilon * ((12 * pow(sigma, 12) / pow(distance, 13)) - (12 * pow(sigma, 6) / pow(distance, 7)))) * ((coords[k][component] - coords[i][component]) / distance);
            }
        }

        return force;
    }
}

double Lennard_Jones::calculate_forward_difference_force(int component, int atom_index, double step_size, vector<vector<double>> coords)
{
    if (!(error < 0))
    {
        int k = atom_index;
        double h = step_size;
        double force = 0;
        double energy = 0;
        double energy_h = 0;

        vector<double> plus_h;

        for (int i = 0; i < 3; i++)
        {
            if (i == component)
            {
                plus_h.push_back(coords[k][i] + h);
                continue;
            }
            plus_h.push_back(coords[k][i]);
        }

        for (int i = 0; i < coords.size(); i++)
        {

            if (i != k)
            {
                energy_h += calculate_energy(plus_h, coords[i]);
                energy += calculate_energy(coords[k], coords[i]);
            }
        }
        force = -(energy_h - energy) / h;

        return force;
    }
}

double Lennard_Jones::calculate_finite_difference_force(int component, int atom_index, double step_size, vector<vector<double>> coords)
{
    if (!(error < 0))
    {
        int k = atom_index;
        double h = step_size;
        double force = 0;
        double energy_plus = 0;
        double energy_minus = 0;

        vector<double> plus_h;
        vector<double> minus_h;
        for (int i = 0; i < 3; i++)
        {
            if (i == component)
            {
                plus_h.push_back(coords[k][i] + h);
                minus_h.push_back(coords[k][i] - h);
                continue;
            }
            plus_h.push_back(coords[k][i]);
            minus_h.push_back(coords[k][i]);
        }

        for (int i = 0; i < coords.size(); i++)
        {

            if (i != k)
            {
                energy_plus += calculate_energy(plus_h, coords[i]);
                energy_minus += calculate_energy(minus_h, coords[i]);
            }
        }
        force = -(energy_plus - energy_minus) / (2 * h);

        return force;
    }
}

// void Lennard_Jones::output_force_file(string file_name)
// {

//     ofstream myfile(file_name);
//     if (myfile.is_open())
//     {
//         myfile << "E_LJ = " << total_energy() << endl;
//         myfile << "F_LJ analytical" << endl;
//         for (int i = 0; i < 3; i++)
//         {
//             for (int j = 0; j < num_of_atoms; j++)
//             {
//                 myfile << calculate_analytic_force(i, j) << " ";
//             }
//             myfile << endl;
//         }

//         myfile << "Stepsize for finite difference:0.1" << endl;
//         myfile << "F_LJ forward difference" << endl;
//         for (int i = 0; i < 3; i++)
//         {
//             for (int j = 0; j < num_of_atoms; j++)
//             {
//                 myfile << calculate_forward_difference_force(i, j, 0.1) << " ";
//             }
//             myfile << endl;
//         }

//         myfile << "F_LJ central difference" << endl;
//         for (int i = 0; i < 3; i++)
//         {
//             for (int j = 0; j < num_of_atoms; j++)
//             {
//                 myfile << calculate_finite_difference_force(i, j, 0.1) << " ";
//             }
//             myfile << endl;
//         }

//         for (double h = 0.1; h > 0.0001; h * 0.1)
//         {

//             myfile << "Stepsize for finite difference: " << h << endl;
//             myfile << "F_LJ forward difference" << endl;
//             for (int i = 0; i < 3; i++)
//             {
//                 for (int j = 0; j < num_of_atoms; j++)
//                 {
//                     myfile << calculate_forward_difference_force(i, j, h) << " ";
//                 }
//                 myfile << endl;
//             }

//             myfile << "F_LJ central difference" << endl;
//             for (int i = 0; i < 3; i++)
//             {
//                 for (int j = 0; j < num_of_atoms; j++)
//                 {
//                     myfile << calculate_finite_difference_force(i, j, h) << " ";
//                 }
//                 myfile << endl;
//             }
//         }
//     }
// }

// void Lennard_Jones::print_output_force()
// {

//     cout << "E_LJ = " << total_energy() << endl;
//     cout << "F_LJ analytical" << endl;
//     for (int i = 0; i < 3; i++)
//     {
//         for (int j = 0; j < num_of_atoms; j++)
//         {
//             cout << calculate_analytic_force(i, j) << " ";
//         }
//         cout << endl;
//     }

// cout << "Stepsize for finite difference:0.1" << endl;
// cout << "F_LJ forward difference" << endl;
// for (int i = 0; i < 3; i++)
// {
//     for (int j = 0; j < num_of_atoms; j++)
//     {
//         cout << calculate_forward_difference_force(i, j, 0.1) << " ";
//     }
//     cout << endl;
// }

// cout << "F_LJ central difference" << endl;
// for (int i = 0; i < 3; i++)
// {
//     for (int j = 0; j < num_of_atoms; j++)
//     {
//         cout << calculate_finite_difference_force(i, j, 0.1) << " ";
//     }
//     cout << endl;
// }
// double h = 0.1;
// while (h > 0.0001)
// {

//     cout << "Stepsize for finite difference: " << h << endl;
//     cout << "F_LJ forward difference" << endl;
//     for (int i = 0; i < 3; i++)
//     {
//         for (int j = 0; j < num_of_atoms; j++)
//         {
//             cout << calculate_forward_difference_force(i, j, h) << " ";
//         }
//         cout << endl;
//     }

//     cout << "F_LJ central difference" << endl;
//     for (int i = 0; i < 3; i++)
//     {
//         for (int j = 0; j < num_of_atoms; j++)
//         {
//             cout << calculate_finite_difference_force(i, j, h) << " ";
//         }
//         cout << endl;
//     }
//     h = h * 0.1;
// }}

vector<vector<double>> Lennard_Jones::return_analytic_force(vector<vector<double>> coords)
{
    if (!(error < 0))
    {
        vector<vector<double>> forces;

        // vector<double> force(3);
        for (int j = 0; j < coords.size(); j++)
        {
            vector<double> force;
            for (int i = 0; i < 3; i++)
            {
                force.push_back(calculate_analytic_force(i, j, coords));
            }
            forces.push_back(force);
        }

        return forces;
    }
}
vector<vector<double>> Lennard_Jones::return_central_difference_force(vector<vector<double>> coords, double stepsize)
{
    if (!(error < 0))
    {
        vector<vector<double>> forces;
        double h = stepsize;

        // vector<double> force(3);
        for (int j = 0; j < coords.size(); j++)
        {
            vector<double> force;
            for (int i = 0; i < 3; i++)
            {
                force.push_back(calculate_finite_difference_force(i, j, h, coords));

                // forces[j][i] = calculate_analytic_force(i, j, coords);
            }
            forces.push_back(force);
        }

        return forces;
    }
}

vector<vector<double>> Lennard_Jones::return_forward_difference_force(vector<vector<double>> coords, double stepsize)
{
    if (!(error < 0))
    {
        vector<vector<double>> forces;
        double h = stepsize;

        for (int j = 0; j < coords.size(); j++)
        {
            vector<double> force;
            for (int i = 0; i < 3; i++)
            {
                force.push_back(calculate_forward_difference_force(i, j, h, coords));
            }
            forces.push_back(force);
        }

        return forces;
    }
}

vector<vector<vector<double>>> Lennard_Jones::optimize(vector<vector<double>> coords, double threshold)
{
    if (!(error < 0))
    {
        vector<vector<double>> current_pos = coords;
        vector<vector<double>> next_pos = coords;
        bool is_converged = false;
        double difference = 1;
        int count = 0;
        while (difference > threshold)
        {

            double current_energy = total_energy(current_pos);

            vector<vector<double>> current_forces = return_central_difference_force(current_pos, 0.0001);

            vector<vector<double>> normalized_force = current_forces;

            for (int i = 0; i < current_pos.size(); i++)
            {

                double magnitude = calculate_distance({0.0, 0.0, 0.0}, current_forces[i]);
                if (magnitude == 0)
                {
                    normalized_force[i] = {0, 0, 0};
                }
                else
                {
                    normalized_force[i][0] /= magnitude;
                    normalized_force[i][1] /= magnitude;
                    normalized_force[i][2] /= magnitude;
                }
            }

            double h = 0.3;
            vector<vector<vector<double>>> visted;
            double a;
            double b;
            double c;
            double energy_a;
            double energy_b;
            double energy_c;

            for (int i = 0; i < coords.size(); i++)
            {

                for (int j = 0; j < 3; j++)
                {
                    a = current_pos[i][j];
                    energy_a = current_energy;
                    b = a - normalized_force[i][j] * h;
                    current_pos[i][j] = b;
                    energy_b = total_energy(current_pos);
                    c = b - normalized_force[i][j] * h;
                    current_pos[i][j] = c;
                    energy_c = total_energy(current_pos);

                    while ((energy_b > energy_a))
                    {
                        h /= 2;
                        b = a + normalized_force[i][j] * h;
                        current_pos[i][j] = b;
                        energy_b = total_energy(current_pos);
                    }

                    h = 0.3;
                    while (energy_c < energy_b)
                    {
                        h *= 1.1;
                        c = b + normalized_force[i][j] * h * 1.1;
                        current_pos[i][j] = c;
                        energy_c = total_energy(current_pos);
                    }

                    double tau = (sqrt(5.0) - 1.0) / 2.0;
                    double x1, x2, energy_x1, energy_x2;

                    double x = b;
                    double energy_x = energy_b;

                    while (std::abs(energy_c - energy_a) > threshold)
                    {

                        x1 = a + (1 - tau) * (c - a);
                        x2 = a + tau * (c - a);
                        current_pos[i][j] = x1;
                        energy_x1 = total_energy(current_pos);
                        current_pos[i][j] = x2;
                        energy_x2 = total_energy(current_pos);

                        if (energy_x1 < energy_x2)
                        {
                            c = x2;
                            energy_c = total_energy(current_pos);
                            x = x1;
                            energy_x = energy_x1;
                        }
                        else
                        {
                            a = x1;
                            energy_a = total_energy(current_pos);
                            x = x2;
                            energy_x = energy_x2;
                        }
                    }
                    count++;

                    current_pos[i][j] = next_pos[i][j];
                    next_pos[i][j] = x;
                }
            }
            double energy_final = total_energy(next_pos);
            difference = std::abs(energy_final - current_energy);
            current_pos = next_pos;
            visted.push_back(current_pos);
            return visted;
        }
    }
}

int main(int argc, char *argv[])
{
    string file_path1 = argv[1];
    int atomic_num = 79;
    Lennard_Jones gold_cluster_1(atomic_num, 5.29, 2.951);
    vector<vector<double>> coords_1 = gold_cluster_1.read_xyz(file_path1);

    ofstream myfile1("output_1.txtx");
    if (myfile1.is_open())
    {
        for (int i = 0; i < coords_1.size(); i++)
        {
            myfile1 << atomic_num << " (" << coords_1[i][0] << ", " << coords_1[i][1]
                    << ", " << coords_1[i][2] << ")" << endl;
        }
        myfile1 << "E_LJ = " << gold_cluster_1.total_energy(coords_1) << endl;
    }

    string file_path2 = argv[2];
    Lennard_Jones gold_cluster_2(atomic_num, 5.29, 2.951);
    vector<vector<double>> coords_2 = gold_cluster_2.read_xyz(file_path2);

    ofstream myfile2("output_2.txt");
    ofstream myfileg("output_g.csv");

    // double force;
    vector<vector<double>> analytic_force;
    vector<vector<vector<double>>> forward_diff;
    vector<vector<vector<double>>> central_diff;

    if (myfile2.is_open())
    {
        myfile2 << "E_LJ = " << gold_cluster_2.total_energy(coords_2) << endl;
        myfile2 << "F_LJ analytical" << endl;
        analytic_force = gold_cluster_2.return_analytic_force(coords_2);

        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < coords_2.size(); j++)
            {
                if (i == 2 && j == coords_2.size() - 1)
                {
                    myfileg << analytic_force[j][i];
                    continue;
                }
                myfile2 << analytic_force[j][i] << " ";
                myfileg << analytic_force[j][i] << ",";
            }
            myfile2 << endl;
        }
        myfileg << endl;
        analytic_force = gold_cluster_2.return_analytic_force(coords_2);

        // myfile << "Stepsize for finite difference:0.1" << endl;
        // myfile << "F_LJ forward difference" << endl;
        // for (int i = 0; i < 3; i++)
        // {
        //     for (int j = 0; j < num_of_atoms; j++)
        //     {

        //         myfile << gold_cluster_2.calculate_forward_difference_force(i, j, 0.1) << " ";
        //     }
        //     myfile << endl;
        // }

        // myfile << "F_LJ central difference" << endl;
        // for (int i = 0; i < 3; i++)
        // {
        //     for (int j = 0; j < num_of_atoms; j++)
        //     {
        //         myfile << gold_clust_2.calculate_finite_difference_force(i, j, 0.1) << " ";
        //     }
        //     myfile << endl;
        // }
        int count = 0;
        double h = 0.1;
        vector<vector<double>> forces;
        for (int k = 0; k < 4; k++)
        {

            myfile2 << "Stepsize for finite difference: " << h << endl;
            myfile2 << "F_LJ forward difference" << endl;
            forces = gold_cluster_2.return_forward_difference_force(coords_2, h);
            forward_diff.push_back(forces);
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < coords_2.size(); j++)
                {
                    if (i == 2 && j == coords_2.size() - 1)
                    {
                        myfileg << forces[j][i];
                        continue;
                    }
                    myfile2 << forces[j][i] << " ";
                    myfileg << forces[j][i] << ",";
                }
                myfile2 << endl;
            }
            myfileg << endl;

            myfile2 << "F_LJ central difference" << endl;
            forces = gold_cluster_2.return_central_difference_force(coords_2, h);
            central_diff.push_back(forces);
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < coords_2.size(); j++)
                {
                    if (i == 2 && j == coords_2.size() - 1)
                    {
                        myfileg << forces[j][i];
                        continue;
                    }
                    myfile2 << forces[j][i] << " ";
                    myfileg << forces[j][i] << ",";
                }
                myfile2 << endl;
            }
            myfileg << endl;
            h *= 0.1;
        }
    }

    string file_path3 = argv[3];
    Lennard_Jones gold_cluster_3(atomic_num, 5.29, 2.951);
    vector<vector<double>> coords_3 = gold_cluster_3.read_xyz(file_path3);

    ofstream myfile3("output_3.txt");

    if (myfile3.is_open())
    {
        int count = 0;
        myfile3 << "start steepest descent with golden section line search" << endl;
        myfile3 << "Initial energy: " << gold_cluster_3.total_energy(coords_3) << endl;
        myfile3 << "Stepsize for central difference is:0.0001";
        myfile3 << "Initial stepsize for line search is:0.3";
        myfile3 << "Threshold for convergence in force is:0.01" << endl;
        myfile3 << "Central Difference Force" << endl;
        vector<vector<double>> central_forces;
        vector<vector<double>> forces;

        central_forces = gold_cluster_3.return_central_difference_force(coords_3, 0.0001);
        central_diff[count] = central_forces;
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < coords_3.size(); j++)
            {

                myfile3 << central_forces[j][i] << " ";
            }
            myfile3 << endl;
        }

        myfile3 << " Start steepest descent with golden section line search using central difference force" << endl;
        myfile3 << "Start golden section search" << endl;
        vector<vector<vector<double>>> visited = gold_cluster_3.optimize(coords_3, 0.01);
        myfile3 << "new_point" << endl;
        for (int i = 0; i < visited.size(); i++)
        {
            myfile3 << "new_point" << endl;
            for (int k = 0; k < 3; k++)
            {
                for (int j = 0; j < visited[i].size(); j++)
                {
                    myfile3 << visited[i][j][k];
                }
                myfile3 << endl;
            }

            myfile3 << "Central Difference Force" << endl;

            central_forces = gold_cluster_3.return_central_difference_force(visited[i], 0.0001);

            for (int k = 0; k < 3; k++)
            {
                for (int j = 0; j < visited[i].size(); j++)
                {
                    myfile3 << central_forces[j][k] << " ";
                }
                myfile3 << endl;
            }
        }

        myfile3 << "Total iterations: " << visited.size() << endl;
        myfile3 << "Final energy: " << gold_cluster_3.total_energy(visited[visited.size() - 1]) << endl;
        myfile3 << "Otpized Structure:" << endl;
        for (int n = 0; n < visited[visited.size() - 1].size(); n++)
        {
            myfile3 << atomic_num << "(";

            for (int z = 0; z < 3; z++)
            {
                myfile3 << visited[visited.size() - 1][n][z] << " ";
            }
            myfile3 << ")" << endl;
        }
        return 0;
    }
}
