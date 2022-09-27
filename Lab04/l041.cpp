#include <iostream> 
#include <fstream>
#include <cstdlib>
#include <vector>
#include <math.h>
#include <ctime>
#include <list>
#include <map>

using namespace std;

/************************
 * You can set the number of clusters and points down below.
 * This lab supports up to 20 different clusters and N points.
 * There are 20 hardcoded colors for the user to see the
 * algorithm at work as this was deemed the maximum possible 
 * colors for them to all appear unique and discernible.
************************/

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
class Color {
    public:
        int r;
        int g;
        int b;
        Color() {
            r = 0;
            g = 0;
            b = 0;
        }
        Color(int valr, int valg, int valb) {
            r = valr;
            g = valg;
            b = valb;
        }
        double getr() {
            return 255-r;
        }
        double getg() {
            return 255-g;
        }
        double getb() {
            return 255-b;
        }
};

int CLUSTERS = 5, POINTS = 100;
int draw[800][800][3];
ofstream ppm("clusters.ppm");

bool valid(int a) {
    return a >= 0 && a < 800;
}
double distance(Point a, Point b) { 
    return sqrt(( (b.getx() - a.getx()) * (b.getx() - a.getx()) ) + ( (b.gety() - a.gety()) * (b.gety() - a.gety()) )); 
}

void part1() {
    srand(time(NULL));
    ppm << "P3\n800 800\n255" << endl;

    vector<Color> colors;
    colors.push_back(Color(255, 0, 0));
    colors.push_back(Color(0, 255, 0));
    colors.push_back(Color(0, 0, 255));
    colors.push_back(Color(255, 255, 0));
    colors.push_back(Color(255, 0, 255));
    colors.push_back(Color(0, 255, 255));
    colors.push_back(Color(255, 102, 0));
    colors.push_back(Color(168, 0, 255));
    colors.push_back(Color(255, 0, 168));
    colors.push_back(Color(0, 127, 255));
    colors.push_back(Color(0, 255, 153));
    colors.push_back(Color(127, 0, 0));
    colors.push_back(Color(0, 127, 0));
    colors.push_back(Color(0, 0, 127));
    colors.push_back(Color(127, 127, 0));
    colors.push_back(Color(127, 0, 127));
    colors.push_back(Color(0, 127, 127));
    colors.push_back(Color(84, 0, 127));
    colors.push_back(Color(127, 0, 84));
    colors.push_back(Color(127, 102, 0));
    colors.push_back(Color(127, 127, 127));

    vector<Point> points;
    double xval, yval;
    for(int i = 0; i < POINTS; i ++) {
        xval = ((double) rand())/RAND_MAX;
        yval = ((double) rand())/RAND_MAX;
        points.push_back(Point(xval, yval));
    }
    
    vector<Point> centroids;
    for(int i = 0; i < CLUSTERS; i++) {
        xval = ((double) rand())/RAND_MAX;
        yval = ((double) rand())/RAND_MAX;
        centroids.push_back(Point(xval, yval));
    }

    map<int, vector<Point>> kmeans;
    
    bool condition = true;
    double d, min, totalx, totaly;
    int cnum, pnum, count;
    vector<int> oldlens(CLUSTERS, 0);
    vector<int> lens(CLUSTERS, 0);
    Point c, p;
    while(condition) {
        for(int i = 0; i < CLUSTERS; i++) {
            vector<Point> temp;
            kmeans[i] = temp;
        }
        for(vector<Point>::iterator it=points.begin(); it != points.end(); ++it) {
            p = *it;
            min = sqrt(2);
            count = 0;
            for(vector<Point>::iterator it2=centroids.begin(); it2 != centroids.end(); ++it2) {
                c = *it2;
                d = distance(p, c);
                if(d < min) {
                    cnum = count;
                    min = d;
                }
                count++;
            }
            kmeans[cnum].push_back(p);
        }
        for(int i = 0; i < CLUSTERS; i++) {
            lens[i] = kmeans[i].size();
        }
        vector<Point> newc;
        for(int i = 0; i < CLUSTERS; i++) {
            pnum = 0, totalx = 0, totaly = 0;
            for(vector<Point>::iterator it=kmeans[i].begin(); it != kmeans[i].end(); ++it) {
                p = *it;
                totalx += p.getx();
                totaly += p.gety();
                pnum++;
            }
            xval = totalx / pnum;
            yval = totaly / pnum;
            newc.push_back(Point(xval, yval));
        }
        condition = false;
        for(int i = 0; i < CLUSTERS; i++) {
            if(oldlens[i] != lens[i])
                condition = true;
        }
        centroids = newc;
        for(int i = 0; i < CLUSTERS; i++) {
            oldlens[i] = lens[i];
        }
    }

    int a0, b0, b1, b1_new, ty, x, y, r;
    double amax;
    Point curr;
    for(vector<Point>::iterator it=centroids.begin(); it != centroids.end(); ++it) {
        curr = *it;
        x = (int) (curr.getx()*800);
        y = (int) (curr.gety()*800);
        r = 3;
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
                draw[a0+x][b0+y][0] = 255;
                draw[a0+x][b0+y][1] = 255;
                draw[a0+x][b0+y][2] = 255;
            }
            if(valid(a0+x) && valid(-b0+y))
            {
                draw[a0+x][-b0+y][0] = 255;
                draw[a0+x][-b0+y][1] = 255;
                draw[a0+x][-b0+y][2] = 255;
            }
            if(valid(-a0+x) && valid(b0+y)) 
            {
                draw[-a0+x][b0+y][0] = 255;
                draw[-a0+x][b0+y][1] = 255;
                draw[-a0+x][b0+y][2] = 255;
            }
            if(valid(-a0+x) && valid(-b0+y))   
            {
                draw[-a0+x][-b0+y][0] = 255;
                draw[-a0+x][-b0+y][1] = 255;
                draw[-a0+x][-b0+y][2] = 255;
            }
            if(valid(a0+y) && valid(b0+x))
            {
                draw[b0+x][a0+y][0] = 255;
                draw[b0+x][a0+y][1] = 255;
                draw[b0+x][a0+y][2] = 255;
            }
            if(valid(-a0+y) && valid(b0+x))
            {
                draw[b0+x][-a0+y][1] = 255;
                draw[b0+x][-a0+y][2] = 255;
                draw[b0+x][-a0+y][0] = 255;
            }
            if(valid(a0+y) && valid(-b0+x))
            {
                draw[-b0+x][a0+y][0] = 255;
                draw[-b0+x][a0+y][1] = 255;
                draw[-b0+x][a0+y][2] = 255;
            }
            if(valid(-a0+y) && valid(-b0+x))
            {
                draw[-b0+x][-a0+y][0] = 255;
                draw[-b0+x][-a0+y][1] = 255;
                draw[-b0+x][-a0+y][2] = 255;
            }
            b1_new -= (2 * a0) - 3;
        }
    }
    for(int i = 0; i < CLUSTERS; i++) {
        for(vector<Point>::iterator it=kmeans[i].begin(); it != kmeans[i].end(); ++it) {
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
                    draw[a0+x][b0+y][0] = colors[i].getr();
                    draw[a0+x][b0+y][1] = colors[i].getg();
                    draw[a0+x][b0+y][2] = colors[i].getb();
                }
                if(valid(a0+x) && valid(-b0+y))
                {
                    draw[a0+x][-b0+y][0] = colors[i].getr();
                    draw[a0+x][-b0+y][1] = colors[i].getg();
                    draw[a0+x][-b0+y][2] = colors[i].getb();
                }
                if(valid(-a0+x) && valid(b0+y)) 
                {
                    draw[-a0+x][b0+y][0] = colors[i].getr();
                    draw[-a0+x][b0+y][1] = colors[i].getg();
                    draw[-a0+x][b0+y][2] = colors[i].getb();
                }
                if(valid(-a0+x) && valid(-b0+y))   
                {
                    draw[-a0+x][-b0+y][0] = colors[i].getr();
                    draw[-a0+x][-b0+y][1] = colors[i].getg();
                    draw[-a0+x][-b0+y][2] = colors[i].getb();
                }
                if(valid(a0+y) && valid(b0+x))
                {
                    draw[b0+x][a0+y][0] = colors[i].getr();
                    draw[b0+x][a0+y][1] = colors[i].getg();
                    draw[b0+x][a0+y][2] = colors[i].getb();
                }
                if(valid(-a0+y) && valid(b0+x))
                {
                    draw[b0+x][-a0+y][0] = colors[i].getr();
                    draw[b0+x][-a0+y][1] = colors[i].getg();
                    draw[b0+x][-a0+y][2] = colors[i].getb();
                }
                if(valid(a0+y) && valid(-b0+x))
                {
                    draw[-b0+x][a0+y][0] = colors[i].getr();
                    draw[-b0+x][a0+y][1] = colors[i].getg();
                    draw[-b0+x][a0+y][2] = colors[i].getb();
                }
                if(valid(-a0+y) && valid(-b0+x))
                {
                    draw[-b0+x][-a0+y][0] = colors[i].getr();
                    draw[-b0+x][-a0+y][1] = colors[i].getg();
                    draw[-b0+x][-a0+y][2] = colors[i].getb();
                }
                b1_new -= (2 * a0) - 3;
            }
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
}

int main() {
    part1();
    return 0;
}