#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
using namespace std;


// tablea

class GaussQuadrature {
public:
    int N;
    vector<double> xi;
    vector<double> w;

    GaussQuadrature(int n) {
        N = n;
        switch (n) {
        case 1: {
            xi = { 0.0 };
            w = { 2.0 };
            break;
        }
        case 2: {
            double a = sqrt(1.0 / 3.0);
            xi = { -a, a };
            w = { 1.0, 1.0 };
            break;
        }
        case 3: {
            double a = sqrt(3.0 / 5.0);
            xi = { -a, 0.0, a };
            w = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };
            break;
        }
        case 4: {
            double a = sqrt((3.0 / 7.0) - (2.0 / 7.0) * sqrt(6.0 / 5.0));
            double b = sqrt((3.0 / 7.0) + (2.0 / 7.0) * sqrt(6.0 / 5.0));
            double w1 = (18.0 - sqrt(30.0)) / 36.0;
            double w2 = (18.0 + sqrt(30.0)) / 36.0;
            xi = { -b, -a, a, b };
            w = { w1, w2, w2, w1 };
            break;
        }
        case 5: {
            double a = (1.0 / 3.0) * sqrt(5.0 - 2.0 * sqrt(10.0 / 7.0));
            double b = (1.0 / 3.0) * sqrt(5.0 + 2.0 * sqrt(10.0 / 7.0));
            double w1 = (322.0 - 13.0 * sqrt(70.0)) / 900.0;
            double w2 = 128.0 / 225.0;
            double w3 = (322.0 + 13.0 * sqrt(70.0)) / 900.0;
            xi = { -b, -a, 0.0, a, b };
            w = { w1, w3, w2, w3, w1 };
            break;
        }
        default:
            cerr << "Nieobslugiwany schemat Gaussa: N=" << n << endl;
            break;
        }
    }
};


double f1(double x) {
    return 5 * x * x + 3 * x + 6;
}

double f2(double x, double y) {
    return 5 * x * x * y * y + 3 * x * y + 6;
}


//  1D 

double gauss1D(double (*f)(double), const GaussQuadrature& G) {
    double sum = 0.0;
    for (int i = 0; i < G.N; ++i)
        sum += G.w[i] * f(G.xi[i]);
    return sum;
}


//  2D

double gauss2D(double (*f)(double, double), const GaussQuadrature& G) {
    double sum = 0.0;
    for (int i = 0; i < G.N; ++i)
        for (int j = 0; j < G.N; ++j)
            sum += G.w[i] * G.w[j] * f(G.xi[i], G.xi[j]);
    return sum;
}



int main() {
    cout << fixed << setprecision(6);

   
    GaussQuadrature G2(2);
    GaussQuadrature G3(3);

 
    cout << "1D" << endl;
    cout << " 2-punktowy " << gauss1D(f1, G2) << endl;
    cout << " 3-punktowy " << gauss1D(f1, G3) << endl;


   
    cout << "\n2D" << endl;
    cout << " 2x2  " << gauss2D(f2, G2) << endl;
    cout << " 3x3  " << gauss2D(f2, G3) << endl;




    return 0;
}
