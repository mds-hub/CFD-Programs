#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>
using namespace std;
double k = 100, h = 50; //diffusion equation variables
double Tinf = 100, Tb = 200;
double T, perim, area, d = 0.08;
double L=1, N{};
int i=1;
int z;

void analytical(){ //void function that computes the analytical solution for the 1-D Diffusion Equation
printf("\nANALYTICAL SOLUTION\n");
double dx = L/N;
i=1;
double msq = (h*perim)/(k*area);
double m = sqrt(msq);
double x = dx/2;

cout << "x" << i << ": " << x << endl;
while(i < (N+1)){
double theta = cosh(m*(L-x))/cosh(m*L);
    ofstream outData;
    outData.open("AnlytTheta.csv", ios::app);
    outData << theta << endl;
        T = theta*(Tb-Tinf)+Tinf;
    outData.open("AnlytTheta.csv", ios::app);
    outData << theta << endl;
    cout << "T" << i << ": " << T << endl;
    outData.open("AnlytTemps.csv", ios::app);
    outData << T << endl;
i++;
        x = x+dx;
    outData.open("distances.csv", ios::app);
    outData << x << endl;
    cout << "x" << i << ": " << x << endl;
}
}

void setup(){ //void function that calculates the tridiagonal coefficients for each mesh size then computes the cell centered values for that mesh size by solving the tridiagonal matrix using the Thomas Algorithm
perim = (M_PI)*d;
area = (M_PI_4)*d*d;
printf("Please Input the Number of Cells\n");
printf("N: ");
cin >> N;
cout << "perim equals: " << perim << endl << "area equals: " << area << endl;
printf("Matrix Size: ");
cin >> z;
float dx = L/N;
double NumTheta{};
double y[z][z];
double q[z];
printf("\nTRIDIAGONAL COEFFICIENTS\n");
while(i < (N+1)){
    if(i == 1){
        double b = (-2*k/dx)-(k/dx)-(((h*perim)/area)*dx);
        double c = k/dx;
        double d = ((-(h*perim)/area)*Tinf*dx)-(2*k/dx)*Tb;
        cout << endl << "i: " << i << endl << "B: " << b << endl << "C: " << c << endl << "D: " << d << endl;
        i++;
        y[0][0] = b;
        y[0][1] = c;
        q[0] = d;
        ofstream triData; 
        triData.open("tridicoeff.csv", ios::app); //opens a csv file entitled "tridicoeff" in order to write data to the file and manipulate it in excel
        triData << 0 << endl << b << endl << c << endl << d << endl; //writes output data to the csv file
    }
        else if(i >= 2 && i < N){   //Interior Cells

        double a = k/dx;
        double b = ((-k/dx)-(k/dx))-(((h*perim)/area)*dx);
        double c = k/dx;
        double d = ((-(h*perim)/area)*Tinf*dx);
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

        double a = (k/dx);
        double b = (-k/dx)+(-(h*perim*dx)/area);
        double d = (-((h*perim)/area))*Tinf*dx;
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
    double T = tad[i-1];
    ofstream outData;
    outData.open("NumTemps.csv", ios::app);
    outData << T << endl;
    printf("\nNUMERICAL SOLUTION\n");
    cout << "T" << i << ": " << T << ", " << "Numerical Theta " << i << ": " << NumTheta << endl;
    NumTheta = (T-Tinf)/(Tb-Tinf);


    while(i > 1){
    T = tad[i-2]-tac[i-2]*T;
    NumTheta = (T-Tinf)/(Tb-Tinf);
    ofstream outData;
    outData.open("NumTemps.csv", ios::app);
    outData << endl << T << endl;
    cout << "T" << i-1 << ": " << T << ", " << "Numerical Theta " << i-1 << ": " << NumTheta << endl;
    NumTheta = (T-Tinf)/(Tb-Tinf);
    i--;
    }



}



int main(){ //outputs all calculated values

    setup();
    analytical();
    system("pause");
    return 0;
}
