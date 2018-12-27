#include <iostream>
#include <vector>
#include <string>
#include "utils.hpp"
using namespace std;
#define  K 20
#define  SERVER_NUM 100

int picNum=0;

void show(){
    picNum++;
    string s;
    s="./show.sh "+std::to_string(picNum);
    system(s.c_str());
}

void display_add(const std::vector<std::pair<Point2d, Point2d>> & results){
    for (auto result:results){
        string s;
        double height=result.second.y-result.first.y;
        double width=result.second.x-result.first.x;
        s="echo "+std::to_string(result.first.x)+" "+std::to_string(result.first.y)+" "
          +std::to_string(height)+" "+std::to_string(width)+" "+" >>out.dat";
        system(s.c_str());
    }
}

int main() {
    int x,y;
    std::vector<Point2d>  servers_batch;
    for (int i = 0; i < SERVER_NUM; ++i) {
        x=rand() % 10;
        y=rand() % 10;
        Point2d tmp(x,y);
        servers_batch.push_back(tmp);
    }
    string s;
    s="rm out.dat; \nrm servers.dat; pwd;";
        system(s.c_str());
    //init display
    for (auto point:servers_batch) {
        string s;
        double height=point.x;
        double width=point.y;

        s="echo "
          +std::to_string(height)+" "+std::to_string(width)+" "+" >>servers.dat";
        system(s.c_str());
    }
    show();

    //input servers
    Model<SERVER_NUM> match;
    match.add_servers(servers_batch);
    std::vector<std::pair<Point2d, Point2d>>  result;

    for (int j = 0; j < SERVER_NUM/K ; ++j) {
    //request, get result
        vector<Point2d>  request;
        for (int i = 0; i < K; ++i) {
            x=rand() % 10;
            y=rand() % 10;
            Point2d tmp(x,y);
            request.push_back(tmp);
        }
        result=match.add_requests(request);
        display_add(result);
        result.clear();
        show();
    }


    return 0;
}