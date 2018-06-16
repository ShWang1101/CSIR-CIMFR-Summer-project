/*|-----------------------------------------------|
  | Problem solved by Priyank Jairaj, BITS Pilani |
  | Github: priyankjairaj100, StopStalk: la_flame |
  |-----------------------------------------------|*/

#include<bits/stdc++.h>
using namespace std;
const long long MOD=1000000007;
#define ll long long
#define el "\n"
ll modl(const ll input, const ll ceil) {
    return input >= ceil ? input % ceil : input;}
#define fr(i,n) for(ll i=0;i<n;++i)
#define fr1(i,n) for(ll i=1;i<=n;++i)
bool myfn(int i, int j) { return i < j; }

ll gcd(ll a, ll b)
{
    if (a == 0 || b == 0)
        return 0;
    if (a == b)
        return a;
    if (a > b)
        return gcd(a - b, b);
    return gcd(a, b - a);
}

void dfs1(ll x,ll c,vector <ll> g[],ll vis[],ll& maxc){
    vis[x]=1;
    c++;
    for(auto i:g[x]){
        if(!vis[i]){
            if(c>=maxc){maxc=c;}//t=i;}//for depth
            dfs1(i,c,g,vis,maxc);
            }
    }
}

void dfs(ll x,vector <ll> g[],ll& maxc){
    ll n;//comment this out
    ll vis[n+1]; fr(i,n+1)vis[i]=0;
    ll c=0;
    dfs1(x,c,g,vis,maxc);
}

int main()
{
    ios::sync_with_stdio(0);
    cin.tie(NULL);
    
    return(0);
}