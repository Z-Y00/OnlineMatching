#pragma once
#include "base.h"

struct Point2d{
    double x, y;
    Point2d() = default;
    Point2d(double xx, double yy);
    static double dis(Point2d a, Point2d b);
    void print(){
        printf("(%f, %f)\n", x, y);
    }
};

struct Edge2d{
    int st, ed;
    double dis;
    Edge2d() = default;
    Edge2d(int s, int t, double d);
};

template<int maxn>
class Graph{
   std::vector<int> g[maxn]; 
   std::vector<Edge2d>e;
public:
   Graph() = default;
   void add_edge(int s, int t, double d);
   std::vector<Edge2d> get_edges(int s);
   Edge2d get_edge(int s, int t);
};

template<int N>
class km{
    constexpr static double inf = 100000000.0;
    int n;
    double val[N][N];//n*n对点之间的距离
    double lx[N],ly[N];//lx:server的偏好值 ly:client的偏好值
    int linky[N];
    int pre[N];
    bool vis[N];//记录点是否被用到
    double slack[N];
    void bfs(int k);
public:
    std::vector<std::pair<Point2d, Point2d> > KM(std::vector<Point2d> s, std::vector<Point2d> r);
};


template<int maxn>
class Model{
    std::vector<Point2d> requests, servers;
    int servers_match[maxn];
    double y[maxn<<1];
    km<maxn+10> kkm;
    double t;
public:
    Model();
    double get_sum(std::vector<std::pair<Point2d, Point2d>>);
    std::vector<std::pair<Point2d, Point2d>> add_requests(std::vector<Point2d> requests_batch);
    void add_servers(std::vector<Point2d> servers_batch);
};


Point2d::Point2d(double xx, double yy):
    x(xx), y(yy){}

double Point2d::dis(Point2d a, Point2d b){
    return sqrt((a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y));
}

Edge2d::Edge2d(int s, int t, double d):
    st(s), ed(t), dis(d){}

template<int maxn>
void Graph<maxn>::add_edge(int s, int t, double d){
    // std::cout << "build " << s << " to " << t << " dis = " << d << std::endl;
    g[s].push_back(e.size());
    e.push_back(Edge2d(s, t, d));
}

template<int maxn>
std::vector<Edge2d> Graph<maxn>::get_edges(int s){
    std::vector<Edge2d> e_ans;
    for(auto e_idx : g[s])
        e_ans.push_back(e[e_idx]);
    return e_ans;
}

template<int maxn>
Edge2d Graph<maxn>::get_edge(int s, int t){
    for(auto e_idx : g[s]){
        if(e[e_idx].ed == t){
            return e[e_idx];
        }
    }
}

template<int maxn>
Model<maxn>::Model(){
    std::fill(servers_match, servers_match+maxn, -1);
    t = maxn*maxn + 1;
    memset(y, 0, sizeof(y));
}
template<int maxn>
std::vector<std::pair<Point2d, Point2d>> Model<maxn>::add_requests(std::vector<Point2d> requests_batch){

    // Save k servers until KM
    std::vector<int> server_cand;

    // Calculate per request
    for(auto new_r : requests_batch){
        new_r.print();
        requests.push_back(new_r);

        // Build weight residual graph
        // Notice in rg & y, request's index add servers.size()
        Graph< maxn<<1 > rg;
        for(int i = 0; i < servers.size(); ++i)
            for(int j = 0; j < requests.size(); ++j){
                if(servers_match[i] == j+servers.size()){
                    // Edge in M* only starts from a server
                    // servers[i].print();
                    double d = Point2d::dis(servers[i], requests[j]) - y[i] - y[j+servers.size()];
                    // std::cout << y[i] << " " << y[j+servers.size()] << std::endl;
                    rg.add_edge(i, j+servers.size(), d);
                }
                else{
                    // Edge NOT in M* only starts from a request 
                    // servers[i].print();
                    // std::cout << y[i] << " " << y[j+servers.size()] << std::endl;
                    double d = t * Point2d::dis(servers[i], requests[j]) - y[i] - y[j+servers.size()];
                    rg.add_edge(j+servers.size(), i, d);
                }
            }

        // Run Dijkstra
        // Get the free request's index
        int st = requests.size() + servers.size() - 1;

        // Min distance
        double d[maxn<<1], min_d;
        int server_m = -1;

        // Save father node
        int fa[maxn<<1];

        // Init
        std::fill(d, d+(maxn<<1), 1000000.0);
        d[st] = 0;
        fa[st] = -1;

        // Data structure for heap
        struct Disp{
            int n;
            double d;
            bool operator<(const Disp &a) const{
                return d > a.d;
            }
            Disp() = default;
            Disp(int nn, double dd): n(nn), d(dd) {}
        };
        std::priority_queue<Disp> pq;
        pq.push(Disp(st, 0));
        while(!pq.empty()){
            Disp now_dp = pq.top();pq.pop();
            // std::cout << "nowd = " << now_dp.d << std::endl;
            int now_n = now_dp.n;
            double now_d = now_dp.d;

            if(now_d > d[now_n]) continue;

            auto nb_e = rg.get_edges(now_dp.n);
            for(auto e_ : nb_e){
                // Loop for every edge
                int to = e_.ed;
                if(d[to] > d[now_n] + e_.dis){
                    // Relax
                    d[to] = d[now_n] + e_.dis;
                    pq.push(Disp(to, d[to]));
                    fa[to] = now_n;
                }
            }
        }
        min_d = 10000000.0;
        // Search P_min_d which ends at a free server
        for(int i = 0; i < servers.size(); ++i){
            if(servers_match[i] != -1) continue;
            if(min_d > d[i]){
                min_d = d[i];
                server_m = i;
            }
        }

        // Update M* using P_server_m
        std::cout << "P* : " << std::endl;
        int now_nd = server_m;
        bool now_nd_is_server = true;
        while(now_nd != -1){
            std::cout << now_nd << std::endl; 
            if(now_nd_is_server){
                servers_match[now_nd] = fa[now_nd];
                y[fa[now_nd]] -= (t-1.0) * Point2d::dis(requests[fa[now_nd]-servers.size()], servers[now_nd]);
                // std::cout << fa[now_nd] <<  " " << now_nd << std::endl;
                // std::cout << y[fa[now_nd]] << std::endl;
            }
            now_nd_is_server = !now_nd_is_server;
            now_nd = fa[now_nd];
        }

        // Update y for EACH node
        std::cout << "mind = " << min_d << std::endl;
        for(int i = 0; i < servers.size()+requests.size(); ++i){
            if(d[i] >= min_d){
                // Do NOTHING
                ;
            }
            else{
                // When dis_i < d
                if(i < servers.size()){
                    // For server
                    y[i] -= (min_d - d[i]);
                }
                else{
                    // For request
                    y[i] += (min_d - d[i]);
                }
            }
        }

        // Finally, push server_m into server cacndidates
        server_cand.push_back(server_m);
        // printf("serverm = %d\n", server_m);
    }
    std::cout << "M*:" << std::endl;
    for(int i = 0 ; i < servers.size(); ++i){
        if(servers_match[i] != -1){
            std::cout << "s(" << i << ") <-> r(" << servers_match[i] <<  ") " << std::endl;
        }
    }

    std::vector<std::pair<Point2d, Point2d>> m_ori, m_km;
    std::vector<Point2d> servers_batch;
    for(int i = 0; i < server_cand.size(); ++i){
        servers_batch.push_back(servers[server_cand[i]]);
        m_ori.push_back(std::make_pair(servers[server_cand[i]], requests_batch[i]));
    }
    
    m_km = kkm.KM(servers_batch, requests_batch);
    std::cout << "Sum of orignial match: " << get_sum(m_ori) << std::endl;
    std::cout << "Sum of km match: " << get_sum(m_km) << std::endl;
    return m_ori;
    // Use KM to determine final match 
    // return km.match(server_cand, requests_batch);

}


template<int maxn>
void Model<maxn>::add_servers(std::vector<Point2d> servers_batch){
    servers = servers_batch;
    t = servers.size()*servers.size() + 1;
}

template<int maxn>
double Model<maxn>::get_sum(std::vector<std::pair<Point2d, Point2d>> m){
    double sum = 0;
    for(auto _p : m){
        // std::cout << "now dis = " << Point2d::dis(_p.first, _p.second) << std::endl;
        // _p.first.print();
        // _p.second.print();
        sum += Point2d::dis(_p.first, _p.second);
    }
    return sum;
}

template<int N>
void km<N>::bfs(int k){
    int px, py = 0,yy = 0;
    double d;
    memset(pre, 0, sizeof(pre));
    memset(slack, inf, sizeof(slack));
    linky[py]=k;
    do{
        px = linky[py],d = inf, vis[py] = 1;
        for(int i = 1; i <= n; i++)
            if(!vis[i]){
                if(slack[i] > lx[px] + ly[i] - val[px][i])
                    slack[i] = lx[px] + ly[i] -val[px][i], pre[i]=py;
                if(slack[i]<d) d=slack[i],yy=i;
            }
        for(int i = 1; i <= n; i++)
            if(vis[i]) lx[linky[i]] -= d, ly[i] += d;
            else slack[i] -= d;
        py = yy;
    }while(linky[py]);
    while(py) linky[py] = linky[pre[py]] , py=pre[py];
}
template<int N>
std::vector<std::pair<Point2d, Point2d> > km<N>::KM(std::vector<Point2d> s, std::vector<Point2d> r)
{
    int nn=s.size();
    n=nn;
    memset(val,0,sizeof(val));
    for(int i=0;i<nn;i++)
    {
        for(int j=0;j<n;j++)
        {
            val[i+1][j+1]=-Point2d::dis(s[i],r[j]);
        }
    }
    std::vector<std::pair<Point2d, Point2d> > result;
    memset(lx, 0, sizeof(lx));
    memset(ly, 0, sizeof(ly));
    memset(linky, 0, sizeof(linky));
    for(int i = 1; i <= nn; i++)
        memset(vis, 0, sizeof(vis)), bfs(i);
    double ans = 0.0;
    for(int i = 1; i <= nn; ++i)
        ans += lx[i] + ly[i];
    for(int i=0;i<nn;i++)
    {
        result.push_back(std::make_pair(s[linky[i]],r[i]));
    }
    std::cout << "ans = " << ans << std::endl;
    return result;
}
