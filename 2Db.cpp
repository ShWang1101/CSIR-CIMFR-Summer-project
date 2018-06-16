#include <bits/stdc++.h>
using namespace std;
#define ll long long
#define size 4

void getCofactor(double A[size][size], double temp[size][size], ll p, ll q, ll n, ll N)
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
double determinant(double A[size][size], ll n, ll N)
{
    double D = 0; // Initialize result

    //  Base case : if matrix contains single element
    if (n == 1)
        return A[0][0];

    double temp[size][size]; // To store cofactors

    ll sign = 1; // To store sign multiplier

    // Iterate for each element of first row
    for (ll f = 0; f < n; f++)
    {
        // Getting Cofactor of A[0][f]
        getCofactor(A, temp, 0, f, n, N);
        D += sign * A[0][f] * determinant(temp, n - 1, N);

        // terms are to be added with alternate sign
        sign = -sign;
    }

    return D;
}

// Function to get adjoint of A[N][N] in adj[N][N].
void adjoint(double A[size][size], double adj[size][size], ll N)
{
    if (N == 1)
    {
        adj[0][0] = 1;
        return;
    }

    // temp is used to store cofactors of A[][]
    ll sign = 1;
    double temp[size][size];

    for (ll i = 0; i < N; i++)
    {
        for (ll j = 0; j < N; j++)
        {
            // Get cofactor of A[i][j]
            getCofactor(A, temp, i, j, N, N);

            // sign of adj[j][i] positive if sum of row
            // and column indexes is even.
            sign = ((i + j) % 2 == 0) ? 1 : -1;

            // Interchanging rows and columns to get the
            // transpose of the cofactor matrix
            adj[j][i] = (sign) * (determinant(temp, N - 1, N));
        }
    }
}

// Function to calculate and store inverse, returns false if
// matrix is singular
bool inverse(double A[size][size], double inverse[size][size], ll N)
{
    // Find determinant of A[][]
    double det = determinant(A, N, N);
    if (det == 0)
    {
        cout << "Singular matrix, can't find its inverse";
        return false;
    }

    // Find adjoint
    double adj[size][size];
    adjoint(A, adj, N);

    // Find Inverse using formula "inverse(A) = adj(A)/det(A)"
    for (ll i = 0; i < N; i++)
        for (ll j = 0; j < N; j++)
            inverse[i][j] = adj[i][j] / float(det);

    return true;
}

void display(double A[size+1][size+1]){
    cout << "______________________________________"<<endl;
    for (ll i = 1; i <= size; ++i)
    {
        for(ll j = 1; j <= size; ++j){
            cout<<setprecision(3)<<A[i][j]<<"    ";
        }
        cout<<endl;
    }
    cout << "______________________________________"<<endl;
}

int main(){
    cout << endl
         << "______________________________________" << endl<<
            "Solving the 3x3 Heat transfer problem."
         << endl
         << "______________________________________" << endl;
    ll NN = 4, N = 3;
    cout<<"Enter the length of the side."<<endl;
    double L; cin>>L;
    double dx = L/N, dy = L/N;
    double K[size + 1][size + 1] = {{0, 0, 0, 0, 0},
                         {0, (dy / (3 * dx) + dx / (3 * dy)), (-dy / (3 * dx) + dx / (6 * dy)), (-dy / (6 * dx) - dx / (6 * dy)), (dy / (6 * dx) - dx / (3 * dy))},
                         {0, (-dy / (3 * dx) + dx / (6 * dy)), (dy / (3 * dx) + dx / (3 * dy)), (dy / (6 * dx) - dx / (3 * dy)), (-dy / (6 * dx) - dx / (6 * dy))},
                         {0, (-dy / (6 * dx) - dx / (6 * dy)), (dy / (6 * dx) - dx / (3 * dy)), (dy / (3 * dx) + dx / (3 * dy)), (-dy / (3 * dx) + dx / (6 * dy))},
                         {0, (dy / (6 * dx) - dx / (3 * dy)), (-dy / (6 * dx) - dx / (6 * dy)), (-dy / (3 * dx) + dx / (6 * dy)), (dy / (3 * dx) + dx / (3 * dy))}};

    //display(KM);
    double KM[size + 1][size + 1] = {{0, 0, 0, 0, 0},
                             {0, (K[1][1] + K[2][2] + K[3][3] + K[4][4]), (K[4][3] + K[1][2]), (K[2][3] + K[1][4]), (K[1][3])},
                             {0, (K[3][4] + K[2][1]), (K[1][1] + K[2][2] + K[3][3] + K[4][4]), (K[2][4]), (K[2][3] + K[1][4])},
                             {0, (K[3][2] + K[4][1]), (K[4][2]), (K[1][1] + K[2][2] + K[3][3] + K[4][4]), (K[4][3] + K[1][2])},
                             {0, (K[3][1]), (K[3][2] + K[4][1]), (K[3][4] + K[2][1]), (K[1][1] + K[2][2] + K[3][3] + K[4][4])}};
    
    //display(KM);
    double F[size+1] = {0,dx*dy/4,dx*dy/4,dx*dy/4,dx*dy/4};
    cout << endl
         << "Enter the known temperatures." << endl
         << "______________________________________"<< endl;
    double t[size + 1];
    for( ll i = 1; i <= size; ++i) cin >> t[i];
    double T[(size - 1) * (size - 1) + 1][size + 1] = {{0,0,0,0,0},
                                                       {0, (t[1]+t[4])/2, (t[1]), (-1), (t[4])},
                                                       {0, (t[1]), (t[1]), (-1), (-1)},
                                                       {0, (t[1]), (t[1]+t[2])/2, (t[2]), (-1)},
                                                       {0, (t[4]), (-1), (-1), (t[4])},
                                                       {0, (-1), (-1), (-1), (-1)},
                                                       {0, (-1), (t[2]), (t[2]), (-1)},
                                                       {0, (t[4]), (-1), (t[3]), (t[3]+t[4])/2},
                                                       {0, (-1), (-1), (t[3]), (t[3])},
                                                       {0, (-1), (t[2]), (t[2]+t[3])/2, (t[3])}};

    double FM[size + 1] = {0,
                           (F[1] + F[2] + F[3] + F[4] - (K[3][1] * T[1][1] + K[3][2] * T[1][2] + K[3][4] * T[1][4] + K[4][1] * T[2][1] + K[4][2] * T[2][2] + K[2][1] * T[4][1] + K[2][4] * T[4][4])),
                           (F[1] + F[2] + F[3] + F[4] - (K[3][1] * T[2][1] + K[3][2] * T[2][2] + K[4][1] * T[3][1] + K[4][2] * T[3][2] + K[4][3] * T[3][3] + K[1][2] * T[6][2] + K[1][3] * T[6][3])),
                           (F[1] + F[2] + F[3] + F[4] - (K[1][3] * T[8][3] + K[1][4] * T[8][4] + K[2][3] * T[7][3] + K[2][4] * T[7][4] + K[2][1] * T[7][1] + K[3][1] * T[4][1] + K[3][4] * T[4][4])),
                           (F[1] + F[2] + F[3] + F[4] - (K[4][2] * T[6][2] + K[4][3] * T[6][3] + K[1][2] * T[9][2] + K[1][3] * T[9][3] + K[1][4] * T[9][4] + K[2][3] * T[8][3] + K[2][4] * T[8][4]))};

    double kmc[size][size], tc[size], fc[size], kmcinv[size][size];
    for( ll i = 1; i<=size; ++i) for( ll j = 1; j<=size; ++j) kmc[i-1][j-1] = KM[i][j];
    for( ll i = 1; i<=size; ++i) fc[i-1] = FM[i];
    for (ll i = 1; i <= size; ++i) tc[i-1] = 0;
    if(inverse(kmc,kmcinv,size));
//  display(kmcinv)
    for (ll i = 0; i < size; ++i)
        for (ll j = 0; j < size; ++j)
            tc[i] += fc[j] * kmcinv[i][j];

    cout<<"______________________________________"<<endl<<"Unknown temperature for the fifth element are:"<<endl<<endl;
    for (ll i = 0; i < size; ++i)
        cout<<"Tempaerature of element #"<<i+1<<" is:    "<<tc[i]<<endl;
    return 0;
}