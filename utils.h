#pragma once
#include "base.h"

struct Point2d{
    double x, y;
    Point2d() = default;
    Point2d(double xx, double yy);
    static double dis(Point2d a, Point2d b);
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
};

template<int maxn>
class KM{
public:
    std::vector<std::pair<Point2d, Point2d>> 
        match(std::vector<Point2d> s, std::vector<Point2d> r);
};


template<int maxn>
class Model{
    std::vector<Point2d> requests, servers;
    int servers_match[maxn];
    int y[maxn<<1];
    KM<maxn> km;
    int t;
public:
    Model();
    std::vector<std::pair<Point2d, Point2d>> add_requests(std::vector<Point2d> requests_batch);
    void add_servers(std::vector<Point2d> servers_batch);
};

