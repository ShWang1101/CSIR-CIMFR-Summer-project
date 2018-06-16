#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <bits/stdc++.h>
using namespace std;
#define ll long long
#define MAX 20
#define el endl
#define MASTER 0
#define FROM_MASTER 1
#define FROM_WORKER 2

double km[MAX][MAX],km1[MAX][MAX],fm[MAX],kminv[MAX][MAX], T[MAX][1], fm1[MAX][1];
ll N;
MPI_Status status;

void getCofactor(double A[MAX][MAX], double temp[MAX][MAX], ll p, ll q, ll n, ll N)
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
double determinant(double A[MAX][MAX], ll n, ll N)
{
    double D = 0; // Initialize result

    //  Base case : if matrix contains single element
    if (n == 1)
        return A[0][0];

    double temp[MAX][MAX]; // To store cofactors

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
void adjoint(double A[MAX][MAX], double adj[MAX][MAX], ll N)
{
    if (N == 1)
    {
        adj[0][0] = 1;
        return;
    }

    // temp is used to store cofactors of A[][]
    ll sign = 1;double temp[MAX][MAX];

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
bool inverse(double A[MAX][MAX], double inverse[MAX][MAX], ll N)
{
    // Find determinant of A[][]
    double det = determinant(A, N, N);
    if (det == 0)
    {
        cout << "Singular matrix, can't find its inverse";
        return false;
    }

    // Find adjoint
    double adj[MAX][MAX];
    adjoint(A, adj, N);

    // Find Inverse using formula "inverse(A) = adj(A)/det(A)"
    for (ll i = 0; i < N; i++)
        for (ll j = 0; j < N; j++)
            inverse[i][j] = adj[i][j] / float(det);

    return true;
}

void display(double A[MAX][MAX], ll N)
{
    for (ll i = 0; i < N; i++)
    {
        for (ll j = 0; j < N; j++)
            cout << A[i][j] << setprecision(2) << " ";
        cout << endl;
    }
}



int main(int argc, char **argv){
    ios::sync_with_stdio(0);
    cin.tie(NULL);

    int task_id, num_tasks;
    ll averow, extra, offset, m_type, rows, source;
    double startwtime, endwtime;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &task_id);
    ll num_workers = num_tasks - 1;

/*----------------------------------------------------------------------------------------*/

    cout<<"Task ID = "<<task_id<<el;

/*-------------------------------MASTER TASK BEGIN-------------------------------------*/

    if(task_id == MASTER)
    {
        cout<<"Number of workers ="<<num_workers<<el;
        ll NN;
        cout<<"Enter the number of nodes in the one dimensional rod."<<el;
        cin>>NN;
        N = NN-2; 
        if(NN<3){cout<<"Atleast two nodes must be present."<<el; return 0;}
        double L;
        cout<<"Enter the length of the rod"<<el;
        cin>>L;
        double dx = L/(NN-1);
        double k11 = 1/dx, k12 = -1/dx, k21 = -1/dx, k22=1/dx;
        double f1 = dx/2, f2 = dx/2;
        /***********************************/
        for(ll i = 1; i <= N ; ++i){
            for(ll j = 1; j <= N; ++j){
                if(i==j)
                    km[i][j] = k22 + k11;
                else if(abs(i-j)==1){
                    if(i>j)
                        km[i][j] = k21;
                    else
                        km[i][j] = k12;
                }
                else
                    km[i][j] = 0;
            }
        }
        /************************************/
        for(ll i = 1; i <= N; ++i)
            fm[i] = f2 + f1;
        cout<<"Enter the temperature at the first node."<<el;
        cin>>T[0][0];
        cout<<"Enter the temperature at the last node."<<el;
        cin>>T[N+1][0];
        fm[1] -= k21*T[0][0];
        fm[N] -= k12*T[N+1][0];
        for(ll i=0; i<N; ++i) fm1[i][0] = fm[i+1];
        for(ll i=0;i<N;++i)for(ll j=0;j<N;++j)km1[i][j]=km[i+1][j+1];
        //if(inverse(km1,kminv,N));

        double t = 0;
        for (ll i = 0; i < N; i++)
        {
            for (ll j = N; j < 2 * N; j++)
            {
                if (i == j - N)
                    km1[i][j] = 1;
                else
                    km1[i][j] = 0;
            }
        }
        for (ll i = 0; i < N; i++)
        {
            t = km1[i][i];
            for (ll j = i; j < 2 * N; j++)
            {
                //if(!t)continue;
                km1[i][j] = km1[i][j] / t;
            }
            for (ll j = 0; j < N; j++)
            {
                if (i != j)
                {
                    t = km1[j][i];
                    for (ll k = 0; k < 2 * N; k++)
                        km1[j][k] = km1[j][k] - t * km1[i][k];
                }
            }
        }
        for (ll i = 0; i < N; ++i)
            for (ll j = 0; j < N; ++j)
                    kminv[i][j] = km1[i][j];
                    
                for (ll i = 1; i <= N; ++i)
                    T[i][0] = 0;

        startwtime = MPI_Wtime();
        /*Send matrix data to worker tasks.*/
        averow = N / num_workers;
        extra = N % num_workers;
        offset = 0;

        m_type = FROM_MASTER;
        for (ll dest = 1; dest <= num_workers; dest++)
        {
            rows = (dest <= extra) ? averow + 1 : averow;
            printf("   sending %d rows to task %d\n", rows, dest);
            MPI_Send(&offset, 1, MPI_INT, dest, m_type, MPI_COMM_WORLD);
            MPI_Send(&N, 1, MPI_INT, dest, m_type, MPI_COMM_WORLD);
            MPI_Send(&rows, 1, MPI_INT, dest, m_type, MPI_COMM_WORLD);
            MPI_Send(&kminv[offset][0], rows * MAX , MPI_DOUBLE, dest, m_type,
                     MPI_COMM_WORLD);
            MPI_Send(&fm1, MAX , MPI_DOUBLE, dest, m_type, MPI_COMM_WORLD);

            offset = offset + rows;
        }

        /* wait for results from all worker tasks */
        m_type = FROM_WORKER;
        for (ll i=1; i<=num_workers; i++)
        {
            source = i;
            MPI_Recv(&offset, 1, MPI_INT, source, m_type, MPI_COMM_WORLD, &status);
            MPI_Recv(&rows, 1, MPI_INT, source, m_type, MPI_COMM_WORLD, &status);
            MPI_Recv(&T[offset+1][0], rows , MPI_DOUBLE, source, m_type, MPI_COMM_WORLD, &status);
        }
        endwtime=MPI_Wtime();
        printf("mxm time = %f  %f Mflops\n",endwtime-startwtime, 2.*N*N/(endwtime-startwtime)*1.e-6);;
        /* print results */
        for(ll i =0;i<N+2;++i)cout<<T[i][0]<<"  ";
    }

    /*-------------------------------MASTER TASK ENDS--------------------------------*/

    /*-------------------------------WORKER TASK STARTS------------------------------*/

    if(task_id > MASTER){
        m_type = FROM_MASTER;
        MPI_Recv(&offset, 1, MPI_INT, MASTER, m_type, MPI_COMM_WORLD, &status);
        MPI_Recv(&N, 1, MPI_INT, MASTER, m_type, MPI_COMM_WORLD, &status);
        MPI_Recv(&rows, 1, MPI_INT, MASTER, m_type, MPI_COMM_WORLD, &status);
        MPI_Recv(&kminv, rows * MAX, MPI_DOUBLE, MASTER, m_type, MPI_COMM_WORLD, &status);
        MPI_Recv(&fm1, MAX * 1, MPI_DOUBLE, MASTER, m_type, MPI_COMM_WORLD, &status);

        for (ll i=0; i<rows; i++)
        {
            for (ll j=0; j<N; j++)
                T[i+1][0] += kminv[i][j] * fm1[j][0];
        }     
        m_type = FROM_WORKER;
        MPI_Send(&offset, 1, MPI_INT, MASTER, m_type, MPI_COMM_WORLD);
        MPI_Send(&rows, 1, MPI_INT, MASTER, m_type, MPI_COMM_WORLD);
        MPI_Send(&T[1][0], rows, MPI_DOUBLE, MASTER, m_type, MPI_COMM_WORLD);
    }
    /*-------------------------------WORKER TASK ENDS--------------------------------*/
    MPI_Finalize();
    return 0;
}