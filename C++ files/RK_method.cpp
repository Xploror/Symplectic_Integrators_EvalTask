// RK_method.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include "cross_p.h"
#include "integrators.h"

using namespace std;

double q = -16 * pow(10, -20);
double m = 9.109 * pow(10, -31);
double step = pow(10, -12);
double dopri_step = pow(10, -14);
vector<long double> E = { 0, 0, 0 };
vector<long double> B = { 0, 0, 1 };

#define f(v) mult(q/m,(add(E, cross_product(v,B))))

int main()
{
    long double c = 3*pow(10,8);
    vector<long double> k1, k2, k3, k4;
    vector<vector<long double>> v = { { 0, 0.9 * c, 0 } };
    vector<vector<long double>> p = { mult(m, v[0]) };
    vector<vector<long double>> x = { {0, 0, 0} };
  
    // Function RK to implement Runge Kutta Order 4  
    RK(v, x, p, E, B, step, m, q);

    //Implementing Leapfrog integrator
    //LeapFrog(v, x, p, E, B, step, m, q);

    //Implementing DormantPrince Method
    //DoPri(v, x, p, E, B, dopri_step, m, q);
    
    // Output position values at each timestep
    for (unsigned int i = 0; i < x.size(); i++)
    {
        cout << "Position (x, y, z) : ";
        for (unsigned int j = 0; j < x[i].size(); j++)
        {
            cout<<x[i][j] << "   ";
        }

        cout << "     Momentum (x, y, z) : ";
        for (unsigned int j = 0; j < p[i].size(); j++)
        {
            cout << p[i][j] << "   ";
        }
        cout << endl;
    }

    // Radius calculation
    double v_mag = sqrt(pow(v[0][0], 2) + pow(v[0][1], 2) + pow(v[0][2], 2));
    vector<long double> F = add(mult(q, E), mult(q, cross_product(v[0], B)));
    double f_mag = sqrt(pow(F[0], 2) + pow(F[1], 2) + pow(F[2], 2));
    double R = m * pow(v_mag, 2) / f_mag;
    cout << "\nRadius : " << R;
}
