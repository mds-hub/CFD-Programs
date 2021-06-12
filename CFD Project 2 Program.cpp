#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>
using namespace std;
double gamma = 0.25; //convection-diffusion equation variables
double phizero = 2, phiL = 1;
double rho = 1;
double u{};
double L = 1, N{};
int i=1;
int z;


void analytical(){ //void function for calculating analytical solution to convection-diffusion equation
printf("\nANALYTICAL SOLUTION\n");
double dx = L/N;
i=1;
double x = dx/2;
double PeL = (rho*u*L)/gamma;
while(i < (N+1)){
cout << "x" << i << ": " << x << endl;
double PeX = (rho*u*x)/gamma;
cout << "PeX: " << PeX << endl;
double phi = (((exp(PeX)-1)/(exp(PeL)-1))*(phiL-phizero))+phizero;
    ofstream outData;
    outData.open("Anlytphi.csv", ios::app);
    outData << phi << endl;
cout << "phi" << i << ": " << phi << endl;
i++;
ofstream outDatadist;
    outDatadist.open("distances.csv", ios::app);
    outDatadist << x << endl;
        x = x+dx;
}
}

void central(){ //void function for numerically solving convection-diffusion equation using the central method to approximate the convection term
printf("Please Input the Number of Cells\n");
printf("N: ");
cin >> N;
printf("Please Input the velocity (u)\n");
printf("u: ");
cin >> u;
printf("Matrix Size: ");
cin >> z;
float dx = L/N;
double F = rho*u;
double y[z][z];
double q[z];
printf("\nTRIDIAGONAL COEFFICIENTS (Central Method)\n"); //the following code computes the tridiagonal coefficients for the convection-diffusion equation and solves the tridiagonal matrix using the Thomas Algorithm
while(i < (N+1)){
    if(i == 1){
        double b = ((F/2)+((3*gamma)/dx));
        double c = (F/2)-(gamma/dx);
        double d = ((F)+((2*gamma)/dx))*phizero;
        cout << endl << "i: " << i << endl << "B: " << b << endl << "C: " << c << endl << "D: " << d << endl;
        i++;
        y[0][0] = b;
        y[0][1] = c;
        q[0] = d;
        ofstream triData;
        triData.open("tridicoeff.csv", ios::app);
        triData << 0 << endl << b << endl << c << endl << d << endl;
    }
        else if(i >= 2 && i < N){   //Interior Cells

        double a = ((-F/2)-((gamma)/dx));
        double b = 2*(gamma/dx);
        double c = (F/2)-(gamma/dx);;
        double d = 0;
        cout << endl << "i: " << i << endl << "A: " << a << endl << "B: " << b << endl << "C: " << c << endl << "D: " << d << endl;

        y[i-1][i-2] = a;
        y[i-1][i-1] = b;
        y[i-1][i] = c;
        q[i-1] = d;

        i++;
        ofstream triData;
        triData.open("tridicoeff.csv", ios::app);
        triData << endl << a << endl << b << endl << c << endl << d << endl;
    }
        else{                            //Right Boundary Condition

        double a = ((-F/2)-((gamma)/dx));
        double b = ((-F/2)+(3*(gamma)/dx));
        double d = ((2*gamma/dx)-(F))*phiL;
        cout << endl << "i: " << i << endl << "A: " << a << endl << "B: " << b << endl << "D: " << d <<endl;

        y[i-1][i-2] = a;
        y[i-1][i-1] = b;
        q[i-1] = d;

        i++;
        ofstream triData;
        triData.open("tridicoeff.csv", ios::app);
        triData << endl << a << endl << b << endl << 0 << endl << d << endl;
    }
}



cout << endl;
double tac[z];
double tad[z];
double final[z];
    //forward sweep
    i=1;
    while(i < z+1){
    if(i == 1){
    double cprime = (y[i-1][i])/(y[i-1][i-1]);
    double dprime = (q[i-1])/(y[i-1][i-1]);
    tac[i-1] = cprime;
    tad[i-1] = dprime;

    i++;
    }
    else if(i >= 2 && i < z){
     double cprime = (y[i-1][i])/(y[i-1][i-1]-(y[i-1][i-2]*tac[i-2]));
     double dprime = (q[i-1]-(y[i-1][i-2]*tad[i-2]))/(y[i-1][i-1]-(y[i-1][i-2]*tac[i-2]));
     tac[i-1] = cprime;
     tad[i-1] = dprime;

     i++;
    }
    else{
        double cprime = tac[i-2];
        double dprime = (q[i-1]-(y[i-1][i-2]*tad[i-2]))/(y[i-1][i-1]-(y[i-1][i-2]*tac[i-2]));
         tad[i-1] = dprime;
         i++;
    }
    }

    i=z;
    double phi = tad[i-1];
    ofstream outData;
    outData.open("NumTemps.csv", ios::app);
    outData << phi << endl;
    printf("\nNUMERICAL SOLUTION (Central Method)\n");
    cout << "Phi" << i << ": " << phi << endl;



    while(i > 1){
    phi = tad[i-2]-tac[i-2]*phi;
    ofstream outData;
    outData.open("NumTemps.csv", ios::app);
    outData << phi << endl;
    cout << "Phi" << i-1 << ": " << phi << endl;
    i--;
    }



}

void upwind(){ //void function that similarly solves the convection-diffusion equation but using the 1st Order Upwind method to approximate the convection term
float dx = L/N;
double F = rho*u;
double y[z][z];
double q[z];
printf("\nTRIDIAGONAL COEFFICIENTS (1st Order Upwind Method)\n");
while(i < (N+1)){
    if(i == 1){
        double b = (F+((3*gamma)/dx));
        double c = -(gamma/dx);
        double d = ((F)+((2*gamma)/dx))*phizero;
        cout << endl << "i: " << i << endl << "B: " << b << endl << "C: " << c << endl << "D: " << d << endl;
        i++;
        y[0][0] = b;
        y[0][1] = c;
        q[0] = d;
        ofstream triDataupwind;
        triDataupwind.open("tridicoeffupwind.csv", ios::app);
        triDataupwind << 0 << endl << b << endl << c << endl << d << endl;
    }
        else if(i >= 2 && i < N){   //Interior Cells

        double a = -(F+(gamma/dx));
        double b = F+((2*gamma)/dx);
        double c = -(gamma/dx);;
        double d = 0;
        cout << endl << "i: " << i << endl << "A: " << a << endl << "B: " << b << endl << "C: " << c << endl << "D: " << d << endl;

        y[i-1][i-2] = a;
        y[i-1][i-1] = b;
        y[i-1][i] = c;
        q[i-1] = d;

        i++;
        ofstream triDataupwind;
        triDataupwind.open("tridicoeffupwind.csv", ios::app);
        triDataupwind << endl << a << endl << b << endl << c << endl << d << endl;
    }
        else{                            //Right Boundary Condition

        double a = (-F-(gamma/dx));
        double b = ((F)+(3*(gamma)/dx));
        double d = ((2*gamma/dx))*phiL;
        cout << endl << "i: " << i << endl << "A: " << a << endl << "B: " << b << endl << "D: " << d <<endl;

        y[i-1][i-2] = a;
        y[i-1][i-1] = b;
        q[i-1] = d;

        i++;
        ofstream triDataupwind;
        triDataupwind.open("tridicoeffupwind.csv", ios::app);
        triDataupwind << endl << a << endl << b << endl << 0 << endl << d << endl;
    }
}



cout << endl;
double tac[z];
double tad[z];
double final[z];
    //forward sweep
    i=1;
    while(i < z+1){
    if(i == 1){
    double cprime = (y[i-1][i])/(y[i-1][i-1]);
    double dprime = (q[i-1])/(y[i-1][i-1]);
    tac[i-1] = cprime;
    tad[i-1] = dprime;

    i++;
    }
    else if(i >= 2 && i < z){
     double cprime = (y[i-1][i])/(y[i-1][i-1]-(y[i-1][i-2]*tac[i-2]));
     double dprime = (q[i-1]-(y[i-1][i-2]*tad[i-2]))/(y[i-1][i-1]-(y[i-1][i-2]*tac[i-2]));
     tac[i-1] = cprime;
     tad[i-1] = dprime;

     i++;
    }
    else{
        double cprime = tac[i-2];
        double dprime = (q[i-1]-(y[i-1][i-2]*tad[i-2]))/(y[i-1][i-1]-(y[i-1][i-2]*tac[i-2]));
         tad[i-1] = dprime;
         i++;
    }
    }

    i=z;
    double phi = tad[i-1];
    ofstream outDataupwind;
    outDataupwind.open("NumTempsupwind.csv", ios::app);
    outDataupwind << phi << endl;
    printf("\nNUMERICAL SOLUTION (1st Order Upwind Method)\n");
    cout << "Phi" << i << ": " << phi << endl;



    while(i > 1){
    phi = tad[i-2]-tac[i-2]*phi;
    ofstream outDataupwind;
    outDataupwind.open("NumTempsupwind.csv", ios::app);
    outDataupwind << phi << endl;
    cout << "Phi" << i-1 << ": " << phi << endl;
    i--;
    }



}

int main(){ //outputs all computed values 

    central();
    upwind();
    analytical();
    system("pause");
    return 0;
}
