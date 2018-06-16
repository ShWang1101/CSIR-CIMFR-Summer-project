/*|-----------------------------------------------|
  | Problem solved by Priyank Jairaj, BITS Pilani |
  | Github: priyankjairaj100, StopStalk: la_flame |
  |-----------------------------------------------|*/

#include <bits/stdc++.h>
using namespace std;
const long long MOD = 1000000007;
#define ll long long
#define el "\n"
#define fr(i, n) for (ll i = 0; i < n; ++i)
#define fr1(i, n) for (ll i = 1; i <= n; ++i)

#define N 4

void getCofactor(ll A[N][N], ll temp[N][N], ll p, ll q, ll n)
{
    ll i = 0, j = 0;

    // Looping for each element of the matrix
    for (ll row = 0; row < n; row++)
    {
        for (ll col = 0; col < n; col++)
        {
            //  Copying into temporary matrix only those element
            //  which are not in given row and column
            if (row != p && col != q)
            {
                temp[i][j++] = A[row][col];

                // Row is filled, so increase row index and
                // reset col index
                if (j == n - 1)
                {
                    j = 0;
                    i++;
                }
            }
        }
    }
}

/* Recursive function for finding determinant of matrix.
   n is current dimension of A[][]. */
ll determinant(ll A[N][N], ll n)
{
    ll D = 0; // Initialize result

    //  Base case : if matrix contains single element
    if (n == 1)
        return A[0][0];

    ll temp[N][N]; // To store cofactors

    ll sign = 1; // To store sign multiplier

    // Iterate for each element of first row
    for (ll f = 0; f < n; f++)
    {
        // Getting Cofactor of A[0][f]
        getCofactor(A, temp, 0, f, n);
        D += sign * A[0][f] * determinant(temp, n - 1);

        // terms are to be added with alternate sign
        sign = -sign;
    }

    return D;
}

// Function to get adjoint of A[N][N] in adj[N][N].
void adjoint(ll A[N][N], ll adj[N][N])
{
    if (N == 1)
    {
        adj[0][0] = 1;
        return;
    }

    // temp is used to store cofactors of A[][]
    ll sign = 1, temp[N][N];

    for (ll i = 0; i < N; i++)
    {
        for (ll j = 0; j < N; j++)
        {
            // Get cofactor of A[i][j]
            getCofactor(A, temp, i, j, N);

            // sign of adj[j][i] positive if sum of row
            // and column indexes is even.
            sign = ((i + j) % 2 == 0) ? 1 : -1;

            // Interchanging rows and columns to get the
            // transpose of the cofactor matrix
            adj[j][i] = (sign) * (determinant(temp, N - 1));
        }
    }
}

// Function to calculate and store inverse, returns false if
// matrix is singular
bool inverse(ll A[N][N], float inverse[N][N])
{
    // Find determinant of A[][]
    ll det = determinant(A, N);
    if (det == 0)
    {
        cout << "Singular matrix, can't find its inverse";
        return false;
    }

    // Find adjoint
    ll adj[N][N];
    adjoint(A, adj);

    // Find Inverse using formula "inverse(A) = adj(A)/det(A)"
    for (ll i = 0; i < N; i++)
        for (ll j = 0; j < N; j++)
            inverse[i][j] = adj[i][j] / float(det);

    return true;
}

void display(float A[N][N])
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
            cout << A[i][j] <<setprecision(2)<< " ";
        cout << endl;
    }
}

int main()
{
    ios::sync_with_stdio(0);
    cin.tie(NULL);
    ll nn = 4;
    ll ndof = 2 * nn;
    ll ne = 6;
    ll node[4][2];
    node[0][0] = 0;
    node[0][1] = 0;
    node[1][0] = 20;
    node[1][1] = 0;
    node[2][0] = 20;
    node[2][1] = 15;
    node[3][0] = 0;
    node[3][1] = 15;
    fr(i, 4) fr(j, 2) node[i][j] *= 12;
    ll conn[6][2] = {{1, 2}, {1, 3}, {2, 4}, {2, 3}, {4, 3}, {1, 4}};
    double A = 1.0;
    ll E = 10000000;
    ll p = 1000;
    ll f[8] = {0, 0, 0, 0, 0, p, 0, 0};//f
    ll isol[4] = {3, 4, 5, 6};//isol
    ll k[ndof][ndof];
    fr(i, ndof) fr(j, ndof) k[i][j] = 0;
    float d[ndof];
    fr(i, ndof) d[i] = 0;
    for (ll e = 0; e < ne; ++e)
    {
        ll n1 = conn[e][0];
        ll n2 = conn[e][1];
        ll x1 = node[n1-1][0];
        ll y1 = node[n1-1][1];
        ll x2 = node[n2-1][0];
        ll y2 = node[n2-1][1];
        double L = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
        double C = (x2 - x1) / L;
        double S = (y2 - y1) / L;
        double C2 = C * C;
        double S2 = S * S;
        double CS = C * S;
        double ke1[nn][nn] = {{C2, CS, -C2, -CS},
                              {CS, S2, -CS, -S2},
                              {-C2, -CS, C2, CS},
                              {-CS, -S2, CS, S2}};
        double ke[nn][nn];
    
        fr(i, nn) fr(j, nn) ke[i][j] = (A * E / L) * ke1[i][j];//ke matrix computed.

        ll sctr[4] = {2 * n1 - 1, 2 * n1, 2 * n2 - 1, 2 * n2};//sctr vector calculated.

        for (ll i = 0; i < nn; ++i)
            for (ll j = 0; j < nn; ++j){
                k[sctr[i] - 1][sctr[j] - 1] = k[sctr[i] - 1][sctr[j] - 1] + ke[i][j];}
                //calculation for k 8x8 matrix.
    }
    ll kx[4][4]; fr(i,4)fr(j,4)kx[i][j]=k[isol[i]-1][isol[j]-1]; float kxinv[4][4];

    ll fx[4][1]; fr(i,4)fx[i][0]= f[isol[i]-1];
    fr(i,4)fr(j,4)kxinv[i][j]=0;
    //calculating the inverse of kx.

    inverse(kx,kxinv);
        //display(kxinv);

    float mult[4][1];//r1xc2 r1=4 c1=4 r2=4 c2=1
    fr(i,4)mult[i][0]=0;

    for (ll i = 0; i < 4; ++i)//r1
        for (ll j = 0; j <1 ; ++j)
            for (ll k = 0; k < 4; ++k)//c1=r2
                mult[i][j] += kxinv[i][k] * (float)fx[k][j];

    fr(i,4)d[isol[i]-1]=mult[i][0];
    cout<<el<<"Nodal Displacements"<<el;
    cout<<el<<"NID  X-DISP  Y-DISP"<<el;
    fr1(i,nn)cout<<i<<" "<<d[2*(i-1)]<<"    "<<d[2*i-1]<<el;
    //calculating k*d;
    float f1[8]; fr(i,8)f1[i]=0;

    for (ll i = 0; i < 8; ++i) //r1
            for (ll j = 0; j < 8; ++j) //c1=r2
                f1[i] += k[i][j] * d[j];

    cout<<el<<"External Forces"<<el;
    cout<<el<<"NID  X-FORCE Y-FORCE"<<el;
    fr1(i,4)
        cout<<el<<i<<"  "<<f1[2*(i-1)]<<"  "<<f1[2*i-1]<<el;

    cout<<el<<"Element Stress"<<el;
    cout<<el<<"EID  STRAIN  STRESS"<<el;
    fr1(e,6){
        ll n1 = conn[e-1][0];
        ll n2 = conn[e-1][1];
        ll x1 = node[n1-1][0];
        ll y1 = node[n1-1][1];
        ll x2 = node[n2-1][0];
        ll y2 = node[n2-1][1];
        float L = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
        float C = (x2 - x1) / L;
        float S = (y2 - y1) / L;
        float B[4] = {-C/L,-S/L,C/L,S/L};
        ll sctr[4] = {2*n1-1,2*n1,2*n2-1,2*n2};
        float dsctr[4];
        fr(i,4)dsctr[i]=d[sctr[i]-1];
        float strain=0;
        for (ll i = 0; i < 4; ++i) //r1
                    strain += B[i] * dsctr[i];
        float stress;
        stress=E*strain;
        cout<<el<<e<<"  "<<strain<<"    "<<stress<<el;
    }
    return (0);
}