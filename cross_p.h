#pragma once
#include <vector>

using namespace std;

vector<long double> cross_product(vector<long double> a, vector<long double> b)
{
	vector<long double> e;
	e.push_back(a[1] * b[2] - b[1] * a[2]);
	e.push_back(-a[0] * b[2] + b[0] * a[2]);
	e.push_back(a[0] * b[1] - b[0] * a[1]);

	return e;
}

vector<long double> add(vector<long double> a, vector<long double> b)
{
	vector<long double> c;
	for (unsigned int i = 0; i < a.size(); i++)
	{
		c.push_back(a[i] + b[i]);
	}

	return c;
}

vector<long double> add(vector<long double> a, vector<long double> b, vector<long double> c)
{
	vector<long double> e;
	for (unsigned int i = 0; i < a.size(); i++)
	{
		e.push_back(a[i] + b[i] + c[i]);
	}

	return e;
}

vector<long double> add(vector<long double> a, vector<long double> b, vector<long double> c, vector<long double> d)
{
	vector<long double> e;
	for (unsigned int i = 0; i < a.size(); i++)
	{
		e.push_back(a[i] + b[i] + c[i] + d[i]);
	}

	return e;
}

vector<long double> add(vector<long double> a, vector<long double> b, vector<long double> c, vector<long double> d, vector<long double> f)
{
	vector<long double> e;
	for (unsigned int i = 0; i < a.size(); i++)
	{
		e.push_back(a[i] + b[i] + c[i] + d[i] + f[i]);
	}

	return e;
}

vector<long double> add(vector<long double> a, vector<long double> b, vector<long double> c, vector<long double> d, vector<long double> f, vector<long double> g)
{
	vector<long double> e;
	for (unsigned int i = 0; i < a.size(); i++)
	{
		e.push_back(a[i] + b[i] + c[i] + d[i] + f[i] + g[i]);
	}

	return e;
}

vector<long double> mult(double a, vector<long double> b)
{
	vector<long double> c;
	for (unsigned int i = 0; i < b.size(); i++)
	{
		c.push_back(a * b[i]);
	}

	return c;
}
