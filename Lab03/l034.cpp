#include <unordered_map>
#include <algorithm>
#include <iostream> 
#include <iterator>
#include <fstream>
#include <cstdlib>
#include <utility>
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
        Point copy() {
            return Point(x, y);
        }
};

struct xvalsort {
    bool operator() (Point i, Point j) {
        return (i.getx() < j.getx());
    }
};

struct yvalsort {
    bool operator() (Point i, Point j) {
        return (i.gety() < j.gety());
    }
};

struct hash_pair { 
    template <class T1, class T2> 
    size_t operator()(const pair<T1, T2>& p) const
    { 
        auto hash1 = hash<T1>{}(p.first); 
        auto hash2 = hash<T2>{}(p.second); 
        return hash1 ^ hash2; 
    } 
}; 

int draw[800][800][3];
Point recurPts[2], part2final[2], part3final[2], part4final[2];
double globalDist = sqrt(2), d2final, d3final, d4final, time2final, time3final, time4final;
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
    for(int i = 0; i < 1000; i ++) {
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
    Point a1, a2;

    vector<Point> points;
    ifstream infile("points.txt");
    double xval, yval;
    string xval_raw, yval_raw;
    double oldxval, oldyval;
    string oldxval_raw, oldyval_raw;
    while(true) {
        infile >> xval_raw;
        infile >> yval_raw;
        xval = stod(xval_raw);
        yval = stod(yval_raw);
        if(xval_raw == oldxval_raw && yval_raw == oldyval_raw && xval == oldxval && yval == oldyval)
            break;
        points.push_back(Point(xval, yval));
        oldxval_raw = xval_raw;
        oldyval_raw = yval_raw;
        oldxval = xval;
        oldyval = yval;
    }

    sort(points.begin(), points.end(), xvalsort());

    chrono::time_point<chrono::system_clock> start = chrono::system_clock::now();

    vector<Point> left(points.begin(), points.begin() + points.size() / 2);
	vector<Point> right(points.begin() + points.size() / 2, points.end());
    d1 = recur(left);
    d2 = recur(right);
    a1 = Point(recurPts[0].getx(), recurPts[0].gety());
    a2 = Point(recurPts[1].getx(), recurPts[1].gety());
    d = min(d1, d2);

    vector<Point> strip;
    Point p, mid = right.at(0);
    for(auto& it : points) {
        p = Point(it.getx(), it.gety());
        if(abs(p.getx()-mid.getx()) <= d) {
            strip.push_back(Point(p.getx(), p.gety()));
        }
    }
    
    bool n;
    double temp;
    Point p1, p2;

    for(auto& it : strip) {
        p1 = it;
        for(auto& it2 : strip) {
            p2 = it2;
            temp = calcDistance(p1, p2);
            if(temp == 0)
                break;
            if(temp < d) {
                a1.x = p1.getx();
                a2.x = p2.getx();
                a1.y = p1.gety();
                a2.y = p2.gety();
                d = temp;
            }
        }
    }

    chrono::time_point<chrono::system_clock> end = chrono::system_clock::now();
    
    part2final[0] = Point(a1.getx(), a1.gety());
    part2final[1] = Point(a2.getx(), a2.gety());
    d2final = d;

    chrono::duration<float> difference = end - start;
    timer = (double) difference.count();
    time2final = timer;
}

void part3(vector<Point> points) {
    srand(time(NULL));
    double timer, d = sqrt(2), d1, d2;
    Point a1, a2;

    sort(points.begin(), points.end(), xvalsort());

    chrono::time_point<chrono::system_clock> start = chrono::system_clock::now();

    vector<Point> left(points.begin(), points.begin() + points.size() / 2);
	vector<Point> right(points.begin() + points.size() / 2, points.end());
    d1 = recur(left);
    d2 = recur(right);
    a1 = Point(recurPts[0].getx(), recurPts[0].gety());
    a2 = Point(recurPts[1].getx(), recurPts[1].gety());
    d = min(d1, d2);

    vector<Point> strip;
    Point p, mid = right.at(0);
    for(auto& it : points) {
        p = Point(it.getx(), it.gety());
        if(abs(p.getx()-mid.getx()) <= d) {
            strip.push_back(Point(p.getx(), p.gety()));
        }
    }
    
    sort(strip.begin(), strip.end(), yvalsort());

    bool n;
    double temp;
    Point p1, p2;

    for(auto& it : strip) {
        p1 = it;
        n = -1;
        for(auto& it2 : strip) {
            if(n >= 15)
                break;
            p2 = it2;
            temp = calcDistance(p1, p2);
            if(temp == 0) {
                n = 0;
                continue;
            }
            if(n >= 0 && temp < d) {
                a1.x = p1.getx();
                a2.x = p2.getx();
                a1.y = p1.gety();
                a2.y = p2.gety();
                d = temp;
            }
            n += 1;
            
        }
    }

    chrono::time_point<chrono::system_clock> end = chrono::system_clock::now();

    part3final[0] = Point(a1.getx(), a1.gety());
    part3final[1] = Point(a2.getx(), a2.gety());
    d3final = d;

    chrono::duration<float> difference = end - start;
    timer = (double) difference.count();
    time3final = timer;
}

void part4(vector<Point> points) {
    srand(time(NULL));
    double timer, d;
    Point a1, a2;

    chrono::time_point<chrono::system_clock> start = chrono::system_clock::now();

    for(int a = points.size()-1; a >= 0; a--) {
        int x = rand() % (a+1);
        Point temp = points[a].copy();  
        points[a] = points[x].copy();  
        points[x] = temp;
    }

    // datafile << "CREATING MAP" << endl;
    d = calcDistance(points[0], points[1]);
    a1 = Point(points[0].getx(), points[0].gety());
    a2 = Point(points[1].getx(), points[1].gety());
    unordered_map<pair<unsigned long long, unsigned long long>, Point, hash_pair> subSquares = {};
    // datafile << "INITIALIZED" << endl;
    int xval, yval, a, b, width = floor(2.0/d)+1;
    double px, py, dtemp;
    pair<unsigned long long, unsigned long long> tpair;
    Point hold;

    // datafile << "Starting Loop..." << endl;
    for(int x = 0; x < points.size(); x++) {
        px = points[x].getx();
        py = points[x].gety();
        // datafile << "\nPoint: (" << px << ", " << py << ")" << endl;
        xval = floor((2*px)/d);
        yval = floor((2*py)/d);
        bool changed = false;
        // datafile << "Going through subSquares: (" << xval << ", " << yval << ")" << endl;
        for(int i = -2; i < 3; i++) {
            a = xval + i;
            // datafile << "[" << width << "]A: " << a << endl;
            if(a < 0)
                continue;
            if(a >= width)
                break;
            for(int j = -2; j < 3; j++) {
                b = yval + j;
                // datafile << "[" << width << "]A: " << a << "\tB: " << b << endl;
                if(b < 0)
                    continue;
                if(b >= width)
                    break;
                tpair = make_pair(a, b);
                // datafile << "subSquare (" << a << ", " << b << ")" << endl;

                if(subSquares.find(tpair) != subSquares.end()) {
                    hold = subSquares[tpair];
                    dtemp = calcDistance(points[x], hold);
                    // datafile << "(" << px << ", " << py << "), (" << hold.x << ", " << hold.y << ")" << endl;
                    // datafile << "Distance: " << dtemp << endl;
                    // datafile << "Current: " << d << endl;
                    if(dtemp < d) {
                        // datafile << "CHANGE" << endl;
                        d = dtemp;
                        a1 = Point(px, py);
                        a2 = Point(hold.getx(), hold.gety());
                        changed = true;
                    }
                }
            }
        }
        if(changed) {
            subSquares = {};
            width = floor(2.0/d)+1;
            // datafile << "New Dictionary [" << width << ", " << width << "]" << endl;
            for(int z = 0; z < x+1; z++) {
                px = points[z].getx();
                py = points[z].gety();
                xval = floor((2*px)/d);
                yval = floor((2*py)/d);
                subSquares[make_pair(xval, yval)] = Point(px, py);
            }
        }
        else {
            subSquares[make_pair(xval, yval)] = Point(px, py);
        }
    }
    // datafile << "Finished Loop..." << endl;

    chrono::time_point<chrono::system_clock> end = chrono::system_clock::now();

    part4final[0] = Point(a1.getx(), a1.gety());
    part4final[1] = Point(a2.getx(), a2.gety());
    d4final = d;

    chrono::duration<float> difference = end - start;
    timer = (double) difference.count();
    time4final = timer;
}

int main() {
    vector<Point> points, points3, points4;
    ifstream infile("points.txt");
    double xval, yval;
    string xval_raw, yval_raw;
    double oldxval, oldyval;
    string oldxval_raw, oldyval_raw;
    while(true) {
        infile >> xval_raw;
        infile >> yval_raw;
        xval = stod(xval_raw);
        yval = stod(yval_raw);
        if(xval_raw == oldxval_raw && yval_raw == oldyval_raw && xval == oldxval && yval == oldyval)
            break;
        points.push_back(Point(xval, yval));
        points3.push_back(Point(xval, yval));
        points4.push_back(Point(xval, yval));
        oldxval_raw = xval_raw;
        oldyval_raw = yval_raw;
        oldxval = xval;
        oldyval = yval;
    }

    part3(points3);
    part4(points4);

    datafile << "\nComplete Recursive Algorithm" << endl;
    datafile << "------------------------------" << endl;
    datafile << "Point1: (" << part3final[0].x << ", " << part3final[0].y << ")" << endl;
    datafile << "Point2: (" << part3final[1].x << ", " << part3final[1].y << ")" << endl;
    datafile << "Distance: " << d3final << endl;
    datafile << "Time: " << time3final << "s" << endl;

    datafile << "\nRandomized Algorithm" << endl;
    datafile << "------------------------------" << endl;
    datafile << "Point1: (" << part4final[0].x << ", " << part4final[0].y << ")" << endl;
    datafile << "Point2: (" << part4final[1].x << ", " << part4final[1].y << ")" << endl;
    datafile << "Distance: " << d4final << endl;
    datafile << "Time: " << time4final << "s" << endl;


    ppm << "P3\n800 800\n255" << endl;

    int a0, b0, b1, b1_new, ty, x, y, r;
    double amax;
    Point curr;
    for(auto& it : points) {
        curr = it;
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
                if(!curr.compare(part3final[0]) && !curr.compare(part3final[1])) {
                    draw[a0+x][b0+y][0] = 255;
                }
            }
            if(valid(a0+x) && valid(-b0+y))
            {
                draw[a0+x][-b0+y][1] = 255;
                draw[a0+x][-b0+y][2] = 255;
                if(!curr.compare(part3final[0]) && !curr.compare(part3final[1])) {
                    draw[a0+x][-b0+y][0] = 255;
                }
            }
            if(valid(-a0+x) && valid(b0+y)) 
            {
                draw[-a0+x][b0+y][1] = 255;
                draw[-a0+x][b0+y][2] = 255;
                if(!curr.compare(part3final[0]) && !curr.compare(part3final[1])) {
                    draw[-a0+x][b0+y][0] = 255;
                }
            }
            if(valid(-a0+x) && valid(-b0+y))   
            {
                draw[-a0+x][-b0+y][1] = 255;
                draw[-a0+x][-b0+y][2] = 255;
                if(!curr.compare(part3final[0]) && !curr.compare(part3final[1])) {
                    draw[-a0+x][-b0+y][0] = 255;
                }
            }
            if(valid(a0+y) && valid(b0+x))
            {
                draw[b0+x][a0+y][1] = 255;
                draw[b0+x][a0+y][2] = 255;
                if(!curr.compare(part3final[0]) && !curr.compare(part3final[1])) {
                    draw[b0+x][a0+y][0] = 255;
                }
            }
            if(valid(-a0+y) && valid(b0+x))
            {
                draw[b0+x][-a0+y][1] = 255;
                draw[b0+x][-a0+y][2] = 255;
                if(!curr.compare(part3final[0]) && !curr.compare(part3final[1])) {
                    draw[b0+x][-a0+y][0] = 255;
                }
            }
            if(valid(a0+y) && valid(-b0+x))
            {
                draw[-b0+x][a0+y][1] = 255;
                draw[-b0+x][a0+y][2] = 255;
                if(!curr.compare(part3final[0]) && !curr.compare(part3final[1])) {
                    draw[-b0+x][a0+y][0] = 255;
                }
            }
            if(valid(-a0+y) && valid(-b0+x))
            {
                draw[-b0+x][-a0+y][1] = 255;
                draw[-b0+x][-a0+y][2] = 255;
                if(!curr.compare(part3final[0]) && !curr.compare(part3final[1])) {
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
    return 0;
}