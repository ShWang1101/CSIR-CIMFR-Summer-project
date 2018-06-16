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


int main()
{
    ll nn;
    ll ne;
    cout<<"Enter the number of nodes."<<el;
    cin>>nn;
    cout<<"Enter the number of elements."<<el;
    cin>>ne;
    ll node[nn][2];
    fr(i,nn){cout<<"Enter the co-ordinates the "<<i+1<<"the element in X Y format."<<el;cin>>node[nn][0]>>node[nn][1];}
    ll ndof = 2 * nn;
    ll conn[ne][2];
    fr(i,ne){cout<<"Enter the nodes between which there is an edge in X Y format."<<el;cin>>conn[i][0]>>conn[i][1];}
    cout<<"Enter the value of A E P"<<el;
    double A;
    ll E;
    ll p;
    cin>>A>>E>>p;
    cout<<"Enter the number of node on which force is exerted."<<el;
    ll x; cin>>x;
    ll f[ndof]; fr(i,ndof)f[i]=0;
    f[x-1]=p;
    cout<<"Enter the number of free nodes."<<el;
    ll numf; cin>>numf;
    ll isol[numf]; cout<<"Enter the free nodes separated by space."<<el;
    fr(i,numf)cin>>isol[i];
    ll k[ndof][ndof];
    fr(i, ndof) fr(j, ndof) k[i][j] = 0;
    ll d[ndof];
    fr(i, ndof) d[i] = 0;
    for (ll e = 0; e < ne; ++e)
    {
        ll n1 = conn[e][0];
        ll n2 = conn[e][1];
        ll x1 = node[n1][0];
        ll y1 = node[n1][1];
        ll x2 = node[n2][0];
        ll y2 = node[n2][1];
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

        ll sctr[nn] = {2 * n1 - 1, 2 * n1, 2 * n2 - 1, 2 * n2};//sctr vector calculated.

        cout<<"Iteration nmber:"<<e<<el;
        for (ll i = 0; i < nn; ++i)
            for (ll j = 0; j < nn; ++j){
                k[sctr[i] - 1][sctr[j] - 1] = k[sctr[i] - 1][sctr[j] - 1] + ke[i][j];}
                //calculation for k 8x8 matrix.


        fr(i, ndof)
        {
            fr(j, ndof) cout << k[i][j] << " ";
            cout << el;
        }
        cout << "************************************" << el;
        fr(i, nn)
            d[isol[i] - 1] = abs(k[isol[i] - 1][isol[i] - 1] - f[isol[i] - 1]);

        //cout << el << "External Forces" << el;
        //cout << "NID  X-FORCE Y-FORCE" << el;
        fr1(i, nn)
        {
            //   cout << i << "  " << d[2 * i - 2] << "  " << d[2 * i - 1] << el;
        }
        //cout<<"Element Stress"<<el;
    }
    return (0);
}