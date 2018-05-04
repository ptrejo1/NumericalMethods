//
//
//  multi grid
//
//  Phoenix Trejo
//  AOS 180     304482946
//
//  Unfortunately the accompanying report
//  providing an overview and analysis of
//  the multigrid solver was deleted.
//

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <ctime>

using namespace std;

ofstream out;

/* Kahan summation algorithim
 ** input: vector
 ** return: sum
 */
double KahanSum(const vector<vector<double>> vor, int n_x, int n_y)
{
    double sum = 0.0;
    double c = 0.0;
    
    for (int a = 1; a <= n_y; a++) {
        for (int b = 1; b <= n_x; b++) {
            double y = vor[a][b]-c;
            double t = sum + y;
            c = (t-sum)-y;
            sum = t;
        }
    }
    return sum;
}

/*
 ** Kahan summation algorithim of absolutue value of value
 */
double KahanSumAbs(const vector<vector<double>> psi, const int a1, const int b1)
{
    double sum = 0.0;
    double c = 0.0;
    
    for (int a = 1; a <= a1; a++) {
        for (int b = 1; b <= b1; b++) {
            double y = abs(psi[a][b])-c;
            double t = sum + y;
            c = (t-sum)-y;
            sum = t;
        }
    }
    return sum;
}

/*
 ** Remove average from grid
 */
void RemoveAvg(const double ksum, vector<vector<double>>& vor_avg, int n_x, int n_y)
{
    double avg = ksum / (n_x*n_y);
    for (int a = 0; a <= n_y+1; a++) {
        for (int b = 0; b <= n_x+1; b++) {
            vor_avg[a][b] -= avg;
        }
    }
}

/*
 ** Compute the norm infinite of grid
 */
double NormInf(vector<vector<double>> vec, int n_x, int n_y)
{
    double max1, max2;
    vector<double> max(n_x);
    
    for (int a = 1; a <= n_y; a++) {
        for (int b = 1; b <= n_x; b++) {
            vec[a][b] = abs(vec[a][b]);
        }
    }
    
    for (int z = 1; z <= n_x; z++) {
        max1 = *max_element(vec[z].begin(), vec[z].end());
        max[z] = max1;
    }
    max2 = *max_element(max.begin(), max.end());
    return max2;
}

/*
 ** Negate values of grid
 */
void NegW(vector<vector<double>>& w, int n_x, int n_y){
    for (int j = 0; j <= n_y; j++) {
        for (int i = 0; i <= n_x; i++) {
            w[j][i] *= -1;
        }
    }
}

/*
 ** Gauss Seidel algorithim
 */
void GSFunction(vector<vector<double>>& res, vector<vector<double>> vor, vector<vector<double>>& psi,int n_x, int n_y, int iy, double delta_x, double delta_y, int sweeps){
    
    int m = 0;
    do {
        // update psi values using SOR method
        for (int j = iy; j <= n_y-1; j++) {
            for (int k = iy; k <= n_x-1; k++) {
                
                res[j][k] = (1.0/pow(delta_x, 2.0))*(psi[j][k+1]-2*psi[j][k]+psi[j][k-1])+(1.0/pow(delta_y, 2.0))*(psi[j+1][k]-2*psi[j][k]+psi[j-1][k])-vor[j][k];
                
                psi[j][k] += (res[j][k])/(2*( (1/pow(delta_x, 2.0)) + (1/pow(delta_y, 2.0)) ));
            }
        }
        m++;
    } while (m < sweeps);
}

/*
 ** Function for interpolation of points between grids
 */
void Interpolation(vector<vector<double>>& Correction, vector<vector<double>> Psi2, int n_x, int n_y){
    for (int j = 1; j <= n_y/2+1; j++) {
        for (int i = 1; i <= n_x/2+1; i++){
            Correction[2*j-1][2*i-1] = Psi2[j][i];
        }
    }
    for (int j = 1; j <= n_y/2+1; j++) {
        for (int i = 1; i <= n_x/2+1; i++){
            Correction[2*j-1][2*i] = 0.5*(Correction[2*j-1][2*i-1]+Correction[2*j-1][2*i+1]);
        }
    }
    
    for (int j = 1; j <= n_y/2; j++) {
        for (int i = 1; i <= n_x; i++){
            Correction[2*j][i] = 0.5*(Correction[2*j-1][i]+Correction[2*j+1][i]);
        }
    }
}

/*
 ** Function to calculate L2 norm error
 */
double L2Norm(vector<vector<double>> Solution, vector<vector<double>> Psi1, int n_x, int n_y){
    
    double error;
    double non;
    double accum = 0;
    for (int h = 1; h <= n_y; h++) {
        for (int n = 1; n < n_x; n++) {
            non = pow(abs(Solution[h][n]-Psi1[h][n]),2);
            accum += non;
        }
    }
    return error = sqrt(accum);
}

/*
 ** function for SOR algorithim for comparison against multigrid
 */
void SORFunction(vector<vector<double>>& res, vector<vector<double>> vor, vector<vector<double>>& psi,int n_x, int n_y, int iy, double delta_x, double delta_y, int sweeps, double alpha, double res_norm, double psi_norm, double vor_norm, double& epsilon){
    
    int m = 0;
    do {
        // update psi values using SOR method
        for (int j = iy; j <= n_y-1; j++) {
            for (int k = iy; k <= n_x-1; k++) {
                
                res[j][k] = (1.0/pow(delta_x, 2.0))*(psi[j][k+1]-2*psi[j][k]+psi[j][k-1])+(1.0/pow(delta_y, 2.0))*(psi[j+1][k]-2*psi[j][k]+psi[j-1][k])-vor[j][k];
                
                psi[j][k] += alpha*(res[j][k])/(2*( (1/pow(delta_x, 2.0)) + (1/pow(delta_y, 2.0)) ));
            }
        }
        
        m++;
        
    } while (m < sweeps);
    
    cout << epsilon << ", " << m << ", ";
}

/*
 ** Multigrid V_cycle algorithim
 */
void MGFunction(vector<vector<double>>& Res1, vector<vector<double>> vor, vector<vector<double>>& Psi1, int n_x, int n_y, double delta_x, double delta_y, vector<vector<double>>& Res2, vector<vector<double>>& Res_sub, vector<vector<double>>& Psi2, double res_norm, double psi_norm, double vor_norm, vector<vector<double>>& Correction){
    
    double epsilon;
    double tol = 1e-7;
    
    do {
        GSFunction(Res1, vor, Psi1, n_x, n_y, 2, delta_x, delta_y, 5);
        
        
        for (int g = 1; g <= n_y/2+1; g++) {
            for (int h = 1; h <= n_x/2+1; h++) {
                Res2[g][h] =  Res1[2*g-1][2*h-1];
            }
        }
        
        NegW(Res2, n_x/2+3, n_y/2+3);
        GSFunction(Res_sub, Res2, Psi2, n_x/2+1, n_y/2+1, 2, delta_x, delta_y, 20);
        
        Interpolation(Correction, Psi2, n_x, n_y);
        
        for (int j = 1; j <= n_y; j++) {
            for (int i = 1; i <= n_x; i++){
                Psi1[j][i] += Correction[j][i];
            }
        }
        
        SORFunction(Res1,vor,Psi1,n_x,n_y,2,delta_x,delta_y, 10, 1.0,res_norm,psi_norm,vor_norm,epsilon);
        
        // calculate residual norm and psi norm to be used in finding epsilon
        res_norm = NormInf(Res1, n_x, n_y);
        psi_norm = KahanSumAbs(Psi1, n_x-1, n_y-1);
        
        // calculate epsilon and read out the value
        epsilon = res_norm/((2/pow(delta_x, 2.0)+2/pow(delta_y, 2.0))*psi_norm+vor_norm);

        
    } while (epsilon > tol);
    
}

/*
 ** Function to print solution
 */
void PrintSoln(vector<vector<double>> Psi1, int n_x, int n_y){
    out.open("Soln.txt");
    for (int i = 1; i <= n_y; i++) {
        out << "[";
        for (int u = 1; u <= n_x; u++) {
            if (u != n_x) {
                out << Psi1[i][u] << ", ";
            } else {
                out << Psi1[i][u];
            }
        }
        out << "], " << endl;
    }
    out.close();
}

int main() {

    //declare variables and grids
    int L_x = 2000;
    int N_x = 201;
    int N_y = 201;
    double delta_x = L_x/(N_x-1);
    double delta_y = L_x/(N_y-1);
    const double pi = acos(-1);
    double k1 = 2*pi/2000;
    double sum;
    double res_norm;
    double psi_norm;
   
    vector<double> x1(N_x+2), y1(N_y+2);
    vector<vector<double>> XY;
    XY.push_back(x1);
    XY.push_back(y1);
    vector<double> vor1(N_x+2);
    vector<vector<double>> vor(N_x+2);
    vector<double> psi1(N_x+2);
    vector<vector<double>> Psi1(N_y+2);
    vector<vector<double>> Psi2(N_x/2+3);
    vector<vector<double>> Res1(N_x+2);
    vector<double> res2(N_x/2+3);
    vector<vector<double>> Res2(N_x/2+3);
    vector<vector<double>> Res_sub(N_x/2+3);
    vector<vector<double>> Correction(N_y+2);
    vector<double> soln(N_x+2);
    vector<vector<double>> Solution(N_y+2);
    
    // set smaller grids
    for (int g = 0; g <= N_y/2+3; g++) {
        Res2[g] = res2;
    }
    Psi2 = Res2;
    Res_sub = Res2;
    
    // set position grid
    for (int even = 0; even <= N_x+1; even++){
        XY[0][even] = even*delta_x;
        XY[1][even] = even*delta_y;
    }
    
    //set vorticity grid and others
    for (int i = 0; i <= N_y; i++) {
        for (int u = 0; u <= N_x; u++){
            vor1[u+1] = -2*pow(k1, 2.0)*sin(k1*XY[0][u])*sin(k1*XY[1][i]);
            soln[u+1] = sin(k1*XY[0][u])*sin(k1*XY[1][i]);
        }
        vor1[N_x+1] = 0;
        soln[N_x+1] = 0;
        vor[i+1] = vor1;
        Solution[i+1] = soln;
        Psi1[i] = psi1;
    }
    vor[0] = psi1;
    vor[N_y+1] = psi1;
    Solution[0] = psi1;
    Solution[N_y+1] = psi1;
    Psi1[N_y+1]= psi1;
    Res1 = Psi1;
    Correction = Psi1;
    
    // satisfy solvability condition
    sum = KahanSum(vor, N_x, N_y);
    RemoveAvg(sum, vor, N_x, N_y);
    double vor_norm = NormInf(vor, N_x, N_y);
    
    //  Multigrid with clock set to time run
    clock_t begin = clock();
    
    MGFunction(Res1, vor, Psi1, N_x, N_y, delta_x, delta_y, Res2, Res_sub, Psi2, res_norm, psi_norm, vor_norm, Correction);
    
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << endl << elapsed_secs;
    
    // l2 norm error and solution output
    double err = L2Norm(Solution, Psi1, N_x, N_y);
    PrintSoln(Psi1, N_x, N_y);
    
    return 0;
}
