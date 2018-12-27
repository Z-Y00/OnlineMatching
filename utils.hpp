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

// Thanks to https://www.jianshu.com/p/36c3de94efd7
template<int MAXN>
class km{
    static constexpr int INF=0x3f3f3f3f;
    double love[MAXN][MAXN],slack[MAXN];
    int match[MAXN],ans[MAXN];
    int visb[MAXN],visg[MAXN];
    double eb[MAXN],eg[MAXN];
    int n;
    double eps=1e-10;
    bool DFS(int girl);
public:
    std::vector<std::pair<Point2d, Point2d> > KM(std::vector<Point2d> s, std::vector<Point2d> r);
};


template<int maxn>
class Model{
    std::vector<Point2d> requests, servers;
    int servers_match[maxn];
    double y[maxn<<1];
    km<(maxn<<1)+10> kkm;
    double t;
public:
    double rm_cost;
    double our_cost;
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
Model<maxn>::Model(): rm_cost(0), our_cost(0){
    std::fill(servers_match, servers_match+maxn, -1);
    t = maxn*maxn + 1;
    memset(y, 0, sizeof(y));
}
template<int maxn>
std::vector<std::pair<Point2d, Point2d>> Model<maxn>::add_requests(std::vector<Point2d> requests_batch){
    static int rounds = 1;
    std::cout << "-------------------Batch " << rounds++ << " --------------------"<< std::endl;
    // Save k servers until KM
    std::vector<int> server_cand;

    // Calculate per request
    for(auto new_r : requests_batch){
        // new_r.print();
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
        // std::cout << "P* : " << std::endl;
        int now_nd = server_m;
        bool now_nd_is_server = true;
        while(now_nd != -1){
            // std::cout << now_nd << std::endl; 
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
        // std::cout << "mind = " << min_d << std::endl;
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
    std::cout << "Total M*:" << std::endl;
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
    double rm_this_cost = get_sum(m_ori), our_this_cost = get_sum(m_km);
    rm_cost += rm_this_cost;
    our_cost += our_this_cost;
    std::cout << "DASHBOARD:" << std::endl;
    std::cout << "K = " << requests_batch.size() << std::endl;
    std::cout << "Cost of orignial match (this batch): " << rm_this_cost << std::endl;
    std::cout << "Cost of km match       (this batch): " << our_this_cost << std::endl;
    std::cout << "Cost of orignial match      (total): " << rm_cost << std::endl;
    std::cout << "Cost of orignial match      (total): " << our_cost << std::endl;
    std::cout << "Total Improvement: " << ((rm_cost/our_cost)-1.0) * 100.0 << "%" << std::endl;
    return m_km;
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


template<int MAXN>
bool km<MAXN>::DFS(int girl)
{
    visg[girl]=true;
    for(int boy=1;boy<=n;boy++)
    {
        if(visb[boy]) continue;
        double gap=abs(eb[boy]+eg[girl]-love[girl][boy]);
        if(gap<=eps)
        {
            visb[boy]=1;
            if(match[boy]==-1||DFS(match[boy]))
            {
                match[boy]=girl;
                return true;
            }
        }
        else slack[boy]=std::min(slack[boy],gap);
    }
    return false;
}
template<int MAXN>
std::vector<std::pair<Point2d, Point2d> > km<MAXN>::KM(std::vector<Point2d> s, std::vector<Point2d> r)
{
    memset(match,-1,sizeof(match));
    memset(eb,0,sizeof(eb));
    memset(love,0,sizeof(love));
    std::vector<std::pair<Point2d, Point2d> > result;
    for(int i = 0; i < MAXN; ++i)
        for(int j = 0; j < MAXN; ++j)
            love[i][j] = 0;
    n=s.size();
    for(int i=1;i<=n;i++)
    {
        for(int j=1;j<=n;j++)
        {
            love[i][j]=-Point2d::dis(s[i-1],r[j-1]);
        }
    }

    for(int i=1;i<=n;i++)
    {
        eg[i]=love[i][1];
        for(int j=2;j<=n;j++)
        {
            eg[i]=std::max(eg[i],love[i][j]);
        }
    }
    for(int i=1;i<=n;i++)
    {
        std::fill(slack+1,slack+1+n,10000000.0);
        while(true)
        {
            memset(visg,0,sizeof(visg));
            memset(visb,0,sizeof(visb));
            if(DFS(i)) break;
            double d=100000000.0;
            for(int i=1;i<=n;i++)
            {
                if(!visb[i]) d=std::min(d,slack[i]);
            }
            for(int i=1;i<=n;i++)
            {
                if(visg[i]) eg[i]-=d;
                if(visb[i]) eb[i]+=d;
                else slack[i]-=d;
            }
        }
    }
    for(int i=1;i<=n;i++)
    {
        result.push_back(std::make_pair(s[match[i]-1],r[i-1]));
    }
    return result;
}