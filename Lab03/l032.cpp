#include <algorithm>
#include <iostream> 
#include <iterator>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <chrono>
#include <math.h>
#include <ctime>
#include <list>

using namespace std;

class Point {
    public:
        double x;
        double y;
        Point() {
            x = 0.0;
            y = 0.0;
        }
        Point(double valx, double valy) {
            x = valx;
            y = valy;
        }
        double getx() {
            return x;
        }
        double gety() {
            return y;
        }
        bool compare(Point a) {
            return x == a.getx() && y == a.gety();
        }
};

struct xvalsort {
    bool operator() (Point i, Point j) {
        return (i.getx() < j.getx());
    }
};

int LIMIT = 1000;
int draw[800][800][3];
Point recurPts[2];
double globalDist = sqrt(2);
ofstream ppm("output.ppm");
ofstream datafile("results.txt");

bool valid(int a) {
    return a >= 0 && a < 800;
}
double calcDistance(Point a, Point b) { 
    return sqrt(( (b.getx() - a.getx()) * (b.getx() - a.getx()) ) + ( (b.gety() - a.gety()) * (b.gety() - a.gety()) )); 
}
double recur(vector<Point> v) {
    double d1, d2, d3;
    if(v.size() == 2)
    {
        d1 = calcDistance(v.at(0), v.at(1));
        if(d1 < globalDist)
        {
            globalDist = d1;
            recurPts[0] = v.at(0);
            recurPts[1] = v.at(1);
        }
        return d1;
    }
    else if(v.size() == 3)
    {
        d1 = calcDistance(v.at(0), v.at(1));
        d2 = calcDistance(v.at(1), v.at(2));
        d3 = calcDistance(v.at(2), v.at(0));
        if(d1 < d2)
        {
            if(d1 < d3)
            {
                if(d1 < globalDist)
                {
                    globalDist = d1;
                    recurPts[0] = v.at(0);
                    recurPts[1] = v.at(1);
                }
                return d1;
            }
            else
            {
                if(d3 < globalDist)
                {
                    globalDist = d3;
                    recurPts[0] = v.at(2);
                    recurPts[1] = v.at(0);
                }
                return d3;
            }
        }
        else
        {
            if(d2 < d3)
            {
                if(d2 < globalDist)
                {
                    globalDist = d2;
                    recurPts[0] = v.at(1);
                    recurPts[1] = v.at(2);
                }
                return d2;
            }
            else
            { 
                if(d3 < globalDist)
                {
                    globalDist = d3;
                    recurPts[0] = v.at(2);
                    recurPts[1] = v.at(0);
                }
                return d3;
            }
        }
    }
    else
    {
        vector<Point> left(v.begin(), v.begin() + v.size() / 2);
        vector<Point> right(v.begin() + v.size() / 2, v.end());
        d1 = recur(left);
        d2 = recur(right);
        return min(d1, d2);
    }
}

void part1() {
    srand(time(NULL));
    ppm << "P3\n800 800\n255" << endl;
    double timer, d = sqrt(2);
    Point a1, a2;

    datafile << "Brute Force Algorithm" << endl;
    datafile << "------------------------------" << endl;

    list<Point> randoms;
    ofstream outfile("points.txt");
    double xval, yval;
    for(int i = 0; i < LIMIT; i ++) {
        xval = ((double) rand())/RAND_MAX;
        yval = ((double) rand())/RAND_MAX;
        randoms.push_back(Point(xval, yval));
        outfile << xval << "  " << yval << endl;
    }
    outfile.close();
    
    bool n;
    double temp;
    Point p1, p2;

    chrono::time_point<chrono::system_clock> start = chrono::system_clock::now();
    for(list<Point>::iterator it=randoms.begin(); it != randoms.end(); ++it) {
        p1 = *it;
        for(list<Point>::iterator it2=randoms.begin(); it2 != randoms.end(); ++it2) {
            p2 = *it2;
            n = p2.compare(p1); 
            if(n)
                break;
            temp = calcDistance(p1, p2);
            if(temp <= d) {
                a1.x = p1.getx();
                a2.x = p2.getx();
                a1.y = p1.gety();
                a2.y = p2.gety();
                d = temp;
            }
        }
    }
    chrono::time_point<chrono::system_clock> end = chrono::system_clock::now();

    datafile << "Point1: (" << a1.x << ", " << a1.y << ")" << endl;
    datafile << "Point2: (" << a2.x << ", " << a2.y << ")" << endl;
    datafile << "Distance: " << d << endl;

    chrono::duration<float> difference = end - start;
    timer = (double) difference.count();
    datafile << "Time: " << timer << "s" << endl;

    int a0, b0, b1, b1_new, ty, x, y, r;
    double amax;
    Point curr;
    for(list<Point>::iterator it=randoms.begin(); it != randoms.end(); ++it) {
        curr = *it;
        x = (int) (curr.getx()*800);
        y = (int) (curr.gety()*800);
        r = 2;
        amax = (r * 2.0*acos(0.0) / 4.0);
        b0 = r;
        b1 = b0 * b0;
        ty = (2 * b0) - 1;
        b1_new = b1;
        for (a0 = 0; a0 <= amax; a0++)
        {
            if ((b1 - b1_new) >= ty)
            {
                b1 -= ty;
                b0 -= 1;
                ty -= 2;
            }
            if(valid(a0+x) && valid(b0+y))
            {
                draw[a0+x][b0+y][1] = 255;
                draw[a0+x][b0+y][2] = 255;
                if(!curr.compare(a1) && !curr.compare(a2)) {
                    draw[a0+x][b0+y][0] = 255;
                }
            }
            if(valid(a0+x) && valid(-b0+y))
            {
                draw[a0+x][-b0+y][1] = 255;
                draw[a0+x][-b0+y][2] = 255;
                if(!curr.compare(a1) && !curr.compare(a2)) {
                    draw[a0+x][-b0+y][0] = 255;
                }
            }
            if(valid(-a0+x) && valid(b0+y)) 
            {
                draw[-a0+x][b0+y][1] = 255;
                draw[-a0+x][b0+y][2] = 255;
                if(!curr.compare(a1) && !curr.compare(a2)) {
                    draw[-a0+x][b0+y][0] = 255;
                }
            }
            if(valid(-a0+x) && valid(-b0+y))   
            {
                draw[-a0+x][-b0+y][1] = 255;
                draw[-a0+x][-b0+y][2] = 255;
                if(!curr.compare(a1) && !curr.compare(a2)) {
                    draw[-a0+x][-b0+y][0] = 255;
                }
            }
            if(valid(a0+y) && valid(b0+x))
            {
                draw[b0+x][a0+y][1] = 255;
                draw[b0+x][a0+y][2] = 255;
                if(!curr.compare(a1) && !curr.compare(a2)) {
                    draw[b0+x][a0+y][0] = 255;
                }
            }
            if(valid(-a0+y) && valid(b0+x))
            {
                draw[b0+x][-a0+y][1] = 255;
                draw[b0+x][-a0+y][2] = 255;
                if(!curr.compare(a1) && !curr.compare(a2)) {
                    draw[b0+x][-a0+y][0] = 255;
                }
            }
            if(valid(a0+y) && valid(-b0+x))
            {
                draw[-b0+x][a0+y][1] = 255;
                draw[-b0+x][a0+y][2] = 255;
                if(!curr.compare(a1) && !curr.compare(a2)) {
                    draw[-b0+x][a0+y][0] = 255;
                }
            }
            if(valid(-a0+y) && valid(-b0+x))
            {
                draw[-b0+x][-a0+y][1] = 255;
                draw[-b0+x][-a0+y][2] = 255;
                if(!curr.compare(a1) && !curr.compare(a2)) {
                    draw[-b0+x][-a0+y][0] = 255;
                }
            }
            b1_new -= (2 * a0) - 3;
        }
    }

    for(int i = 0; i < 800; i++)
    {
        for(int j = 0; j < 800; j++)
        {
            ppm << (255-draw[i][j][0]) << " " << (255-draw[i][j][1]) << " " << (255-draw[i][j][2]) << "\t";
        }
        ppm << "\n";
    }

    ppm.close();
}

void part2() {
    srand(time(NULL));
    double timer, d = sqrt(2), d1, d2;
    Point a1, a2, a3, a4;

    datafile << "\nPreliminary Recursive Algorithm" << endl;
    datafile << "------------------------------" << endl;

    vector<Point> points;
    ifstream infile("points.txt");
    double xval, yval;
    string xval_raw, yval_raw;
    for(int i = 0; i < LIMIT; i ++) {
        infile >> xval_raw;
        infile >> yval_raw;
        xval = stod(xval_raw);
        yval = stod(yval_raw);
        points.push_back(Point(xval, yval));
    }

    sort(points.begin(), points.end(), xvalsort());

    chrono::time_point<chrono::system_clock> start = chrono::system_clock::now();

    vector<Point> left(points.begin(), points.begin() + points.size() / 2);
	vector<Point> right(points.begin() + points.size() / 2, points.end());
    Point l1, l2, r1, r2;
    d1 = recur(left);
    l1 = Point(recurPts[0].getx(), recurPts[0].gety());
    l2 = Point(recurPts[1].getx(), recurPts[1].gety());
    d2 = recur(right);
    r1 = Point(recurPts[0].getx(), recurPts[0].gety());
    r2 = Point(recurPts[1].getx(), recurPts[1].gety());
    if(d1 < d2)
    {
        a1 = Point(l1.getx(), l1.gety());
        a2 = Point(l2.getx(), l2.gety());
        d = d1;
    }
    else
    {
        a1 = Point(r1.getx(), r1.gety());
        a2 = Point(r2.getx(), r2.gety());
        d = d2;
    }

    vector<Point> strip;
    Point mid = right.at(0);
    strip.push_back(Point(mid.getx(), mid.gety()));
    for(auto& it : points)
        if(calcDistance(it, mid))
            strip.push_back(Point(it.getx(), it.gety()));
    
    bool n;
    double temp, dt = sqrt(2);
    Point p1, p2;

    for(auto& it : points) {
        p1 = it;
        for(auto& it2 : points) {
            p2 = it2;
            n = p2.compare(p1); 
            if(n)
                break;
            temp = calcDistance(p1, p2);
            if(temp < dt) {
                a3.x = p1.getx();
                a4.x = p2.getx();
                a3.y = p1.gety();
                a4.y = p2.gety();
                dt = temp;
            }
        }
    }
    
    if(dt < d)
    {
        a1.x = a3.getx();
        a1.y = a3.gety();
        a2.x = a4.getx();
        a2.y = a4.gety();
        d = dt;
    }

    chrono::time_point<chrono::system_clock> end = chrono::system_clock::now();
    
    datafile << "Point1: (" << a1.x << ", " << a1.y << ")" << endl;
    datafile << "Point2: (" << a2.x << ", " << a2.y << ")" << endl;
    datafile << "Distance: " << d << endl;

    chrono::duration<float> difference = end - start;
    timer = (double) difference.count();
    datafile << "Time: " << timer << "s" << endl;
}

int main() {
    part1();
    part2();
    return 0;
}