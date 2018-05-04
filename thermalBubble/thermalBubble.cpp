//
//
//  aos180
//
//  Phoenix Trejo
//
//

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <numeric>

using namespace std;

/* Kahan summation algorithim
 ** input: vector
 ** return: sum
 */
double KahanSum(const vector<vector<double>> vor)
{
    double sum = 0.0;
    double c = 0.0;
    
    for (int a = 1; a <= 201; a++) {
        for (int b = 1; b <= 200; b++) {
            double y = vor[a][b]-c;
            double t = sum + y;
            c = (t-sum)-y;
            sum = t;
        }
    }
    return sum;
}

/* Kahan summation algorithim with absolute value
 ** input: psi, a1, a1
 ** return: sum
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

/* Remove average from grid points
 ** input: vor, ksum
 */
void RemoveAvg(const double ksum, vector<vector<double>>& vor_avg)
{
    double avg = ksum / 40401.0;
    for (int a = 0; a <= 202; a++) {
        for (int b = 0; b <= 202; b++) {
            vor_avg[a][b] -= avg;
        }
    }
}

/* Find infinite norm
 ** input: vec
 ** return: max2
 */
double NormInf(vector<vector<double>> vec)
{
    double max1, max2;
    vector<double> max(201);
    
    for (int a = 1; a <= 201; a++) {
        for (int b = 1; b <= 201; b++) {
            vec[a][b] = abs(vec[a][b]);
        }
    }
    
    for (int z = 1; z <= 201; z++) {
        max1 = *max_element(vec[z].begin(), vec[z].end());
        max[z] = max1;
    }
    max2 = *max_element(max.begin(), max.end());
    return max2;
}

/* Set psi back to zerp
 ** input: psi
 */
void ZeroPsi(vector<vector<double>>& psi){
    for (int i = 0; i <= 202; i++) {
        for (int j = 0; j <= 202; j++) {
            psi[i][j] = 0;
        }
    }
}

/* Function for solving psi with SOR
 */
void SORFunction(vector<vector<double>>& res, vector<vector<double>> vor, vector<vector<double>>& psi, double alpha, double delta_x, double delta_y, double res_norm, double psi_norm, double vor_norm, double epsilon, int m, double tol){
    
    do {
        // update psi values using SOR method
        for (int j = 2; j <= 200; j++) {
            for (int k = 1; k <= 201; k++) {
                
                res[j][k] = (1.0/pow(delta_x, 2.0))*(psi[j][k+1]-2*psi[j][k]+psi[j][k-1])+(1.0/pow(delta_y, 2.0))*(psi[j+1][k]-2*psi[j][k]+psi[j-1][k])-vor[j][k];
                
                psi[j][k] += alpha*(res[j][k])/(2*( (1/pow(delta_x, 2.0)) + (1/pow(delta_y, 2.0)) ));
            }
        }
        
        // calculate residual norm and psi norm to be used in finding epsilon
        res_norm = NormInf(res);
        psi_norm = KahanSumAbs(psi, 201, 200);
        
        // calculate epsilon and read out the value
        epsilon = res_norm/((2/pow(delta_x, 2.0)+2/pow(delta_y, 2.0))*psi_norm+vor_norm);
        
        m++;
    } while (epsilon > tol && m < 10000);
    cout << m << ", ";
}

/* Function for making w negative to solve poisson eqn
 */
void NegW(vector<vector<double>>& w){
    for (int j = 0; j <= 202; j++) {
        for (int i = 0; i <= 202; i++) {
            w[j][i] *= -1;
        }
    }
}

int main() {
    
    // declaration of neccessary variables
    int temp_0 = 300;
    int x_0 = 1000;
    int z_0 = 260;
    int r_0 = 250;
    int N_x = 203;
    int N_y = 203;
    int delta_t = 2;
    int T_tot = 620;
    int N_t = T_tot/delta_t;
    int delta_x = 10;
    int delta_z = 10;
    int m = 0;
    double max1;
    double K = 0.1;
    double delta_temp = 0.5;
    double gravity = 9.81;
    double sum;
    double w_norm;
    double alpha = 1.9;
    double epsilon;
    double psi_norm;
    double res_norm;
    double tol = 1e-7;
    vector<double> help(N_y,300);
    
    vector<double> x1(N_y), z1(N_y);
    vector<vector<double>> XZ;
    XZ.push_back(x1);
    XZ.push_back(z1);
    
    vector<double> psi1(N_y);
    vector<vector<double>> psi(N_y);
    vector<vector<double>> res(N_y);
    
    vector<double> temp1(N_y);
    vector<vector<double>> temp(N_y);
    vector<vector<double>> temp_current(N_y);
    vector<vector<double>> temp_new(N_y);
    
    vector<vector<double>> w(N_y);
    vector<vector<double>> w_current(N_y);
    vector<vector<double>> w_new(N_y);
    
    vector<vector<double>> jac_w(N_y);
    vector<vector<double>> jac_w2(N_y);
    vector<vector<double>> jac_w3(N_y);
    vector<vector<double>> jac_t(N_y);
    vector<vector<double>> jac_t2(N_y);
    vector<vector<double>> jac_t3(N_y);

    //output streams
    ofstream out;
    ofstream fout;
    fout.open("max_temp.txt");
    
    // setting position grids
    for (int even = 0; even <= N_x-1; even++) {
        XZ[0][even] = even*delta_x;
        XZ[1][even] = even*delta_z;
    }
    
    // set intitial temperature conditions and set all grids needed for problem
    for (int i = 0; i <= N_y-2; i++) {
        for (int u = 0; u <= N_x-2; u++) {
            if (sqrt(pow((XZ[0][u]-x_0), 2.0)+pow((XZ[1][i]-z_0), 2.0)) <= r_0) {
                temp1[u+1] = temp_0+delta_temp;
            } else {
                temp1[u+1] = temp_0;
            }
        }
        temp1[0] = 300;
        temp1[N_y-1] = 300;
        temp[i+1] = temp1;
        w[i] = psi1;
    }
    temp[0] = help;
    temp[N_y-1] = help;
    temp_current = temp;
    temp_new = temp;
    w[N_y-1] = psi1;
    w_current = w;
    w_new = w;
    psi = w;
    res = w;
    jac_t = w;
    jac_t2 = w;
    jac_t3 = w;
    jac_w = w;
    jac_w2 = w;
    jac_w3 = w;

    // set intitial w grid
    for (int j = 2; j <= N_y-3; j++) {
        for (int i = 1; i <= N_x-2; i++) {
            w[j][i] = (-1*gravity/temp_0)*1/delta_x*(temp[j][i]-temp[j][i-1]);
        }
    }
    // prepare grids for poisson solver and satisfy solvability condition
    NegW(w);
    sum = KahanSum(w);
    RemoveAvg(sum, w);
    w_norm = NormInf(w);
    
    // Poisson Solver
    SORFunction(res, w, psi, alpha, delta_x, delta_z, res_norm, psi_norm, w_norm, epsilon, m, tol);
    NegW(w);
    
    // First time advancement for vorticity and temp
    for (int j = 2; j <= N_y-3; j++) {
        for (int i = 1; i <= N_x-2; i++) {
            jac_t[j][i] = (1/(2*delta_x))*(psi[j][i+1]-psi[j][i-1])*(1/(2*delta_z))*(temp[j+1][i]-temp[j-1][i])-(1/(2*delta_z))*(psi[j+1][i]-psi[j-1][i])*(1/(2*delta_x))*(temp[j][i+1]-temp[j][i-1]);
            
            temp_current[j][i] = temp[j][i]-delta_t*jac_t[j][i];
            
            jac_w[j][i] = (1/(2*delta_x))*(psi[j][i+1]-psi[j][i-1])*(1/(2*delta_z))*(w[j+1][i]-w[j-1][i])-(1/(2*delta_z))*(psi[j+1][i]-psi[j-1][i])*(1/(2*delta_x))*(w[j][i+1]-w[j][i-1]);
            
            w_current[j][i] = w[j][i]-delta_t*(jac_w[j][i]+gravity/temp_0*(temp[j][i+1]-temp[j][i-1])/(2.0*delta_x) );
        }
    }
    
    // Poisson Solver
    NegW(w_current);
    SORFunction(res, w_current, psi, alpha, delta_x, delta_z, res_norm, psi_norm, w_norm, epsilon, m, tol);
    NegW(w_current);
    
    // MAIN LOOP
    for (int a = 2; a <= 30; a++) {
        
        // Advance vorticity and temperature in time
        for (int j = 2; j <= N_y-3; j++) {
            for (int i = 1; i <= N_x-2; i++) {
                
                jac_t[j][i] = (1/(4*pow(delta_x, 2.0)))*((psi[j][i+1]-psi[j][i-1])*(temp_current[j+1][i]-temp_current[j-1][i])-(psi[j+1][i]-psi[j-1][i]) * (temp_current[j][i+1]-temp_current[j][i-1]));
                
                jac_t2[j][i] = (1/(4*pow(delta_x, 2.0)))*(psi[j][i+1]*(temp_current[j+1][i+1]-temp_current[j-1][i+1])-psi[j][i-1]*(temp_current[j+1][i-1]-temp_current[j-1][i-1])-psi[j+1][i]*(temp_current[j+1][i+1]-temp_current[j+1][i-1])+psi[j-1][i]*(temp_current[j-1][i+1]-temp_current[j-1][i-1]));
                
                jac_t3[j][i] = (1/(4*pow(delta_x, 2.0)))*(temp_current[j+1][i]*(psi[j+1][i+1]-psi[j+1][i-1])-temp_current[j-1][i]*(psi[j-1][i+1]-psi[j-1][i-1])-temp_current[j][i+1]*(psi[j+1][i+1]-psi[j-1][i+1])+temp_current[j][i-1]*(psi[j+1][i-1]-psi[j-1][i-1]));
                
                temp_new[j][i] = temp[j][i]-2*delta_t*( (jac_t[j][i]+jac_t2[j][i]+jac_t3[j][i])/3 )+K*(2*delta_t/pow(delta_x, 2.0))*(temp_current[j][i+1]-2*temp_current[j][i]+temp_current[j][i-1])+K*(2*delta_t/pow(delta_z, 2.0))*(temp_current[j+1][i]-2*temp_current[j][i]+temp_current[j-1][i]);
                
                jac_w[j][i] = (1/(4*pow(delta_x, 2.0)))*((psi[j][i+1]-psi[j][i-1])*(w_current[j+1][i]-w_current[j-1][i])-(psi[j+1][i]-psi[j-1][i]) * (w_current[j][i+1]-w_current[j][i-1]));
                
                jac_w2[j][i] = (1/(4*pow(delta_x, 2.0)))*(psi[j][i+1]*(w_current[j+1][i+1]-w_current[j-1][i+1])-psi[j][i-1]*(w_current[j+1][i-1]-w_current[j-1][i-1])-psi[j+1][i]*(w_current[j+1][i+1]-w_current[j+1][i-1])+psi[j-1][i]*(w_current[j-1][i+1]-w_current[j-1][i-1]));
                
                jac_w3[j][i] = (1/(4*pow(delta_x, 2.0)))*(w_current[j+1][i]*(psi[j+1][i+1]-psi[j+1][i-1])-w_current[j-1][i]*(psi[j-1][i+1]-psi[j-1][i-1])-w_current[j][i+1]*(psi[j+1][i+1]-psi[j-1][i+1])+w_current[j][i-1]*(psi[j+1][i-1]-psi[j-1][i-1]));
                
                w_new[j][i] = w[j][i]-2*delta_t*( ((jac_w[j][i]+jac_w2[j][i]+jac_w3[j][i])/3) + gravity/temp_0*(temp_current[j][i+1]-temp_current[j][i-1])/(2.0*delta_x))+K*(2*delta_t/pow(delta_x, 2.0))*(w_current[j][i+1]-2*w_current[j][i]+w_current[j][i-1])+K*(2*delta_t/pow(delta_z, 2.0))*(w_current[j+1][i]-2*w_current[j][i]+w_current[j-1][i]);

            }
            // boundary conditions
            temp_new[j][N_x-1] = temp_new[j][2];
            temp_new[j][0] = temp_new[j][N_x-3];
            w_new[j][N_x-1] = w_new[j][2];
            w_new[j][0] = w_new[j][N_x-3];
        }
        
        // Poisson Solver
        NegW(w_new);
        SORFunction(res, w_new, psi, alpha, delta_x, delta_z, res_norm, psi_norm, w_norm, epsilon, m, tol);
        NegW(w_new);
        
        // Update neccessary grids
        temp = temp_current;
        temp_current = temp_new;
        w = w_current;
        w_current = w_new;
        
        // Find max temperature
        max1 = NormInf(temp_current);
        fout << max1 << ", ";
        
        // readout solution
        string path = "temp"+to_string(a)+".txt";
        out.open(path);
        for (int i = 1; i <= N_y-2; i++) {
            out << "[";
            for (int u = 1; u <= N_x-2; u++) {
                if (u != N_x-2) {
                    out << temp_current[i][u] << ", ";
                } else {
                    out << temp_current[i][u];
                }
            }
            out << "], " << endl;
        }
        out.close();
    }
    
    return 0;
}
