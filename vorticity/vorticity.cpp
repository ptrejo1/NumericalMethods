//
//
//  aos180hw6
//
//  Phoenix Trejo
//  304482946
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
    
    for (int a = 1; a <= 101; a++) {
        for (int b = 1; b <= 101; b++) {
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
    double nodes = 101.0 * 101.0;
    double avg = ksum / nodes;
    for (int a = 0; a <= 102; a++) {
        for (int b = 0; b <= 102; b++) {
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
    vector<double> max(101);
    
    for (int a = 1; a <= 101; a++) {
        for (int b = 1; b <= 101; b++) {
            vec[a][b] = abs(vec[a][b]);
        }
    }
    
    for (int z = 1; z <= 101; z++) {
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
    for (int i = 0; i <= 102; i++) {
        for (int j = 0; j <= 102; j++) {
            psi[i][j] = 0;
        }
    }
}

/* Function for solving psi with SOR
 */
void SORFunction(vector<vector<double>>& res, vector<vector<double>> vor, vector<vector<double>>& psi,int iy, int ey, double alpha, double delta_x, double delta_y, double res_norm, double psi_norm, double vor_norm, double epsilon, int m, double tol){
    
    do {
        // update psi values using SOR method
        for (int j = iy; j <= ey; j++) {
            for (int k = iy; k <= ey; k++) {
                
                res[j][k] = (1.0/pow(delta_x, 2.0))*(psi[j][k+1]-2*psi[j][k]+psi[j][k-1])+(1.0/pow(delta_y, 2.0))*(psi[j+1][k]-2*psi[j][k]+psi[j-1][k])-vor[j][k];
                
                psi[j][k] += alpha*(res[j][k])/(2*( (1/pow(delta_x, 2.0)) + (1/pow(delta_y, 2.0)) ));
            }
        }
        
        // calculate residual norm and psi norm to be used in finding epsilon
        res_norm = NormInf(res);
        psi_norm = KahanSumAbs(psi, 100, 100);
        
        // calculate epsilon and read out the value
        epsilon = res_norm/((2/pow(delta_x, 2.0)+2/pow(delta_y, 2.0))*psi_norm+vor_norm);
        
        m++;
    } while (epsilon > tol && m < 10000);
    cout << m << ", ";
}

/* Solve for enstrophy
 ** input: vor
 ** return: sum2
 */
double enstrophy(vector<vector<double>> vor){
    double sum1, sum2 = 0;
    for (int j = 1; j <= 101; j++) {
        for (int k = 1; k <= 101; k++) {
            vor[j][k] = pow(vor[j][k], 2.0);
        }
    }
    
    for (int t = 1; t <= 100; t++) {
        sum1 = accumulate(vor[t].begin()+1, vor[t].end()-1, 0.0);
        sum2 += sum1;
    }
    return sum2;
}

int main() {
    
    double a = 180.0;
    int L_x = 2000;
    
    // number of grid nodes
    int N_x = 103;
    int N_y = 103;
    
    // Centers of vortices
    int x_01 = 1300;
    int x_02 = 700;
    int y_0 = 1000;
    
    // time step information
    double K = 0.1;
    int delta_t = 1000;
    double delta_x = 5000;
    double delta_y = 5000;
    double alpha = 1.9;
    double gamma_term = 8e-5;
    
    // vortex values
    double sum;
    double vor_norm;
    double epsilon;
    double psi_norm;
    double res_norm;
    double tol = 1e-7;
    
    // node helpers
    int ix = 1;
    int iy = 2;
    int ex = 101;
    int ey = 100;
    
    vector<double> x1, y1;
    vector<vector<double>> XY;
    XY.push_back(x1);
    XY.push_back(y1);
    
    vector<double> psi1(N_x);
    vector<vector<double>> psi(N_x);
    vector<double> vor1(N_x, 0);
    
    vector<vector<double>> vor(N_x);
    vector<vector<double>> vor_current(N_x);
    vector<vector<double>> vor_new(N_x);
    vector<vector<double>> res(N_x);
    vector<vector<double>> u_g(N_x);
    vector<vector<double>> jac(N_x);
    vector<vector<double>> jac2(N_x);
    vector<vector<double>> jac3(N_x);
    
    ofstream out;
    out.open("/Users/phoenixtrejo/Desktop/num/vorticity.txt");
    
    for (int even = 0; even <= L_x+40; even+=20) {
        XY[0].push_back(even);
        XY[1].push_back(even);
    }
    
    // loop for setting vorticity field and other grids
    vor[N_x-1] = vor1;
    vor[0] = vor1;
    u_g[N_x-1] = psi1;
    psi[N_x-1] = psi1;
    res[N_x-1] = psi1;
    jac[N_x-1] = psi1;
    jac2[N_x-1] = psi1;
    jac3[N_x-1] = psi1;
    vor_current[N_x-1] = psi1;
    vor_new[N_x-1] = psi1;
    for (int i = 0; i <= N_y-2; i++) {
        for (int u = 0; u <= N_x-2; u++) {
            double vparam1 = (pow((XY[0][u]-x_01), 2.0)+pow((XY[1][i]-y_0), 2.0))/pow(a, 2.0);
            double vparam2 = (pow((XY[0][u]-x_02), 2.0)+pow((XY[1][i]-y_0), 2.0))/pow(a, 2.0);
            vor1[u+1] = gamma_term*exp(-1*vparam1)+gamma_term*exp(-1*vparam2);
        }
        vor1[N_x-1] = 0;
        u_g[i] = psi1;
        res[i] = psi1;
        psi[i] = psi1;
        jac[i] = psi1;
        jac2[i] = psi1;
        jac3[i] = psi1;
        vor[i+1] = vor1;
        vor_current[i] = psi1;
        vor_new[i] = psi1;
    }
    // remove vorticity average and solve for vorticity normal
    sum = KahanSum(vor);
    RemoveAvg(sum, vor);
    vor_norm = NormInf(vor);
    
    // Solve for psi
    int m = 0;
    SORFunction(res, vor, psi, iy, ey, alpha, delta_x, delta_y, res_norm, psi_norm, vor_norm, epsilon, m, tol);
    
    // advance vorticity with EF
    for (int j = ix; j <= ex; j++) {
        for (int k = ix; k <= ex; k++) {
            jac[j][k] = (1/(2*delta_x))*(psi[j][k+1]-psi[j][k-1])*(1/(2*delta_y))*(vor[j+1][k]-vor[j-1][k])-(1/(2*delta_y))*(psi[j+1][k]-psi[j-1][k])*(1/(2*delta_x))*(vor[j][k+1]-vor[j][k-1]);
            vor_current[j][k] = vor[j][k]-delta_t*jac[j][k];
        }
        vor_current[j][1] = 0.0;
        vor_current[j][101] = 0.0;
    }
    vor1.assign(N_x,0);
    vor_current[1] = vor1;
    vor_current[101] = vor1;
    
    // Solve for new psi
    ZeroPsi(psi);
    SORFunction(res, vor_current, psi, iy, ey, alpha, delta_x, delta_y, res_norm, psi_norm, vor_norm, epsilon, m, tol);
    
    for (int i = 2; i < 30; i++) {
        
        //loop for solving jacobians and advancing vorticity in time
        for (int j = ix; j <= ex; j++) {
            for (int k = ix; k <= ex; k++) {
                
                jac2[j][k] = (1/(4*pow(delta_x, 2.0)))*(psi[j][k+1]*(vor_current[j+1][k+1]-vor_current[j-1][k+1])-psi[j][k-1]*(vor_current[j+1][k-1]-vor_current[j-1][k-1])-psi[j+1][k]*(vor_current[j+1][k+1]-vor_current[j+1][k-1])+psi[j-1][k]*(vor_current[j-1][k+1]-vor_current[j-1][k-1]));
                
                jac3[j][k] = (1/(4*pow(delta_x, 2.0)))*(vor_current[j+1][k]*(psi[j+1][k+1]-psi[j+1][k-1])-vor_current[j-1][k]*(psi[j-1][k+1]-psi[j-1][k-1])-vor_current[j][k+1]*(psi[j+1][k+1]-psi[j-1][k+1])+vor_current[j][k-1]*(psi[j+1][k-1]-psi[j-1][k-1]));
                
                jac[j][k] = (1/(4*pow(delta_x, 2.0)))*((psi[j][k+1]-psi[j][k-1])*(vor_current[j+1][k]-vor_current[j-1][k])-(psi[j+1][k]-psi[j-1][k]) * (vor_current[j][k+1]-vor_current[j][k-1]));
                
                vor_new[j][k] = vor[j][k]-2*delta_t*( (jac[j][k]+jac2[j][k]+jac3[j][k])/3 )+K*(2*delta_t/pow(delta_x, 2.0))*(vor_current[j][k+1]-2*vor_current[j][k]+vor_current[j][k-1])+K*(2*delta_t/pow(delta_y, 2.0))*(vor_current[j+1][k]-2*vor_current[j][k]+vor_current[j-1][k]);
                
            }
            // boundary conditions
            vor_new[j][1] = 0;
            vor_new[j][101] = 0;
        }
        // boundary conditions
        vor_new[1] = vor1;
        vor_new[101] = vor1;
        
        //ZeroPsi(psi);
        SORFunction(res, vor_new, psi, iy, ey, alpha, delta_x, delta_y, res_norm, psi_norm, vor_norm, epsilon, m, tol);
        
        // update values
        vor = vor_current;
        vor_current = vor_new;
    }
    
    // output for vorticity
    for (int i = 1; i <= 101; i++) {
        for (int u = 1; u <= 101; u++) {
            if (u != 101) {
                out << vor_current[i][u] << ",";
            } else {
                out << vor_current[i][u];
            }
        }
        out << endl;
    }
    
    return 0;
}
