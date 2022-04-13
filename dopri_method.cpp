// RK_method.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <iomanip>
#include <vector>
#include "cross_p.h"

using namespace std;

double q = -16 * pow(10, -20);
double m = 9.109 * pow(10, -31);
double step = pow(10, -14);
vector<long double> E = { 0, 0, 0 };
vector<long double> B = { 0, 0, 1 };

#define f(v) mult(q/m,(add(E, cross_product(v,B))))

int main()
{
    long double c = 3 * pow(10, 8);
    long double err_x, err_p;
    double eps = pow(10, -7);
    vector<long double> k1, k2, k3, k4, k5, k6, k7, x_k, p_k;
    vector<vector<long double>> v = { { 0, 0.9 * c, 0 } };
    vector<vector<long double>> z = { { 0, 0.9 * c, 0} };
    vector<vector<long double>> zx = { { 0, 0, 0} };
    vector<vector<long double>> zp = { mult(m, z[0]) };
    vector<vector<long double>> p = { mult(m, v[0]) };
    vector<vector<long double>> x = { {0, 0, 0} };

    for (double i = 0; i < 5 * pow(10, -11); i += step)
    {

        k1 = mult(step, f(v[int(i / step)]));
        k2 = mult(step, f(add(v[int(i / step)], mult(0.2, k1))));
        k3 = mult(step, f(add(v[int(i / step)], add(mult(3/double(40), k1), mult(9/double(40), k2)))));
        k4 = mult(step, f(add(v[int(i / step)], add(mult(44/double(45), k1), mult(-56/double(15), k2), mult(32/double(9), k3)))));
        k5 = mult(step, f(add(v[int(i / step)], add(mult(19372 / double(6561), k1), mult(-25360 / double(2187), k2), mult(64448 / double(6561), k3), mult(-212 / double(729), k4)))));
        k6 = mult(step, f(add(v[int(i / step)], add(mult(9017 / double(3168), k1), mult(-355 / double(33), k2), mult(46732 / double(5247), k3), mult(49 / double(176), k4), mult(-5103 / double(18656), k5)))));
        k7 = mult(step, f(add(v[int(i / step)], add(mult(35 / double(384), k1), mult(500 / double(1113), k3), mult(125 / double(192), k4), mult(-2187 / double(6784), k5), mult(11 / double(84), k6)))));
        v.push_back(add(v[int(i / step)], mult(4.68*pow(10,-8), add(mult(1921409, k1), mult(9690880, k3), mult(13122270, k4), mult(-5802111, k5), mult(1902912, k6), mult(534240, k7)))));
        z.push_back(add(z[int(i / step)], add(mult(5179 / double(57600), k1), mult(7571 / double(16695), k3), mult(393 / double(640), k4), mult(-92097 / double(339200), k5), mult(187 / double(2100), k6), mult(1 / double(40), k7))));
        p.push_back(mult(m, v[1 + int(i / step)]));
        x_k = add(x[int(i / step)], mult(step, v[1 + int(i / step)]));
        x.push_back(x_k);
        zx.push_back(add(zx[int(i / step)], mult(step, z[1 + int(i / step)])));
        zp.push_back(mult(m, z[1 + int(i / step)]));

        err_x = sqrt(pow(zx[i][0] - x[i][0], 2) + pow(zx[i][1] - x[i][1], 2) + pow(zx[i][2] - x[i][2], 2));
        err_p = sqrt(pow(zp[i][0] - p[i][0], 2) + pow(zp[i][1] - p[i][1], 2) + pow(zp[i][2] - p[i][2], 2));
        if (err_x<eps*sqrt(pow(x[i][0],2) + pow(x[i][1],2) + pow(x[i][2],2)) && err_p<eps*sqrt(pow(p[i][0], 2) + pow(p[i][1], 2) + pow(p[i][2], 2)))
        {
            step = pow((eps * step / (2 * err_x)), 0.25) * step;
        }
        //cout << step;
    }


    // Output position values at each timestep
    for (unsigned int i = 0; i < x.size(); i++)
    {
        cout << "Position (x, y, z) : ";
        for (unsigned int j = 0; j < x[i].size(); j++)
        {
            cout << x[i][j] << "   ";
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
    cout << "\nRadius : " <<R;
}
