// RK_method.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include "cross_p.h"

using namespace std;

double q = -16 * pow(10, -20);
double m = 9.109 * pow(10, -31);
double step = pow(10, -12);
vector<long double> E = { 0, 0, 0 };
vector<long double> B = { 0, 0, 1 };

#define f(v) mult(q/m,(add(E, cross_product(v,B))))

int main()
{
    long double c = 3*pow(10,8);
    vector<long double> k1, k2, k3, k4, x_k;
    vector<vector<long double>> v = { { 0, 0.9 * c, 0 } };
    vector<vector<long double>> p = { mult(m, v[0]) };
    vector<vector<long double>> x = { {0, 0, 0} };
    
    for (double i = 0; i < 5*pow(10,-10); i+=step)
    {
        
        k1 = mult(step, f(v[int(i/step)]));
        k2 = mult(step, f(add(v[int(i / step)], mult(0.5,k1))));
        k3 = mult(step, f(add(v[int(i / step)], mult(0.5, k2))));
        k4 = mult(step, f(add(v[int(i / step)], k3)));
        v.push_back(add(v[int(i/step)], mult(0.16667, add(k1, mult(2, k2), mult(2, k3), k4))));
        p.push_back(mult(m, v[1 + int(i / step)]));
        x_k = add(x[int(i / step)], mult(step, v[1 + int(i / step)]));
        x.push_back(x_k);
    }
    
    
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
