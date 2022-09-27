#include <iostream> 
#include <fstream>
#include <cstdlib>
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

int draw[800][800][3];
ofstream ppm("output.ppm");

bool valid(int a) {
    return a >= 0 && a < 800;
}
double distance(Point a, Point b) { 
    return sqrt(( (b.getx() - a.getx()) * (b.getx() - a.getx()) ) + ( (b.gety() - a.gety()) * (b.gety() - a.gety()) )); 
}

void part1() {
    srand(time(NULL));
    ofstream outfile("points.txt");
    ppm << "P3\n800 800\n255" << endl;
    
    list<Point> randoms;

    double xval, yval;
    for(int i = 0; i < 50; i ++) {
        xval = ((double) rand())/RAND_MAX;
        yval = ((double) rand())/RAND_MAX;
        randoms.push_back(Point(xval, yval));
        outfile << xval << "  " << yval << endl;
    }
    
    bool n;
    double temp, d = sqrt(2);
    Point p1, p2, a1, a2;
    for(list<Point>::iterator it=randoms.begin(); it != randoms.end(); ++it) {
        p1 = *it;
        for(list<Point>::iterator it2=randoms.begin(); it2 != randoms.end(); ++it2) {
            p2 = *it2;
            n = p2.compare(p1); 
            if(n)
                break;
            temp = distance(p1, p2);
            // outfile << endl;
            // outfile << (int) (p1.getx()*800.0) << "  " << (int) (p1.gety()*800.0) << endl;
            // outfile << (int) (p2.getx()*800.0) << "  " << (int) (p2.gety()*800.0) << endl;
            // outfile << temp << "\t" << n << "\t" << d << endl;
            if(temp <= d) {
                a1.x = p1.getx();
                a2.x = p2.getx();
                a1.y = p1.gety();
                a2.y = p2.gety();
                d = temp;
                // outfile << "CHANGE" << endl;
            }
        }
    }

    // outfile << endl;
    // outfile << (int) (a1.getx()*800.0) << "  " << (int) (a1.gety()*800.0) << endl;
    // outfile << (int) (a2.getx()*800.0) << "  " << (int) (a2.gety()*800.0) << endl;
    // outfile << d << endl;

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
                draw[a0+x][b0+y][0] += 255;
                if(!curr.compare(a1) && !curr.compare(a2)) {
                    draw[a0+x][b0+y][1] += 255;
                    draw[a0+x][b0+y][2] += 255;
                }
            }
            if(valid(a0+x) && valid(-b0+y))
            {
                draw[a0+x][-b0+y][0] += 255;
                if(!curr.compare(a1) && !curr.compare(a2)) {
                    draw[a0+x][-b0+y][1] += 255;
                    draw[a0+x][-b0+y][2] += 255;
                }
            }
            if(valid(-a0+x) && valid(b0+y)) 
            {
                draw[-a0+x][b0+y][0] += 255;
                if(!curr.compare(a1) && !curr.compare(a2)) {
                    draw[-a0+x][b0+y][1] += 255;
                    draw[-a0+x][b0+y][2] += 255;
                }
            }
            if(valid(-a0+x) && valid(-b0+y))   
            {
                draw[-a0+x][-b0+y][0] += 255;
                if(!curr.compare(a1) && !curr.compare(a2)) {
                    draw[-a0+x][-b0+y][1] += 255;
                    draw[-a0+x][-b0+y][2] += 255;
                }
            }
            if(valid(a0+y) && valid(b0+x))
            {
                draw[b0+x][a0+y][0] += 255;
                if(!curr.compare(a1) && !curr.compare(a2)) {
                    draw[b0+x][a0+y][1] += 255;
                    draw[b0+x][a0+y][2] += 255;
                }
            }
            if(valid(-a0+y) && valid(b0+x))
            {
                draw[b0+x][-a0+y][0] += 255;
                if(!curr.compare(a1) && !curr.compare(a2)) {
                    draw[b0+x][-a0+y][1] += 255;
                    draw[b0+x][-a0+y][2] += 255;
                }
            }
            if(valid(a0+y) && valid(-b0+x))
            {
                draw[-b0+x][a0+y][0] += 255;
                if(!curr.compare(a1) && !curr.compare(a2)) {
                    draw[-b0+x][a0+y][1] += 255;
                    draw[-b0+x][a0+y][2] += 255;
                }
            }
            if(valid(-a0+y) && valid(-b0+x))
            {
                draw[-b0+x][-a0+y][0] += 255;
                if(!curr.compare(a1) && !curr.compare(a2)) {
                    draw[-b0+x][-a0+y][1] += 255;
                    draw[-b0+x][-a0+y][2] += 255;
                }
            }
            b1_new -= (2 * a0) - 3;
        }
    }

    for(int i = 0; i < 800; i++)
    {
        for(int j = 0; j < 800; j++)
        {
            ppm << draw[i][j][0] << " " << draw[i][j][1] << " " << draw[i][j][2] << "\t";
        }
        ppm << "\n";
    }

    outfile.close();
    ppm.close();
}

int main() {
    part1();
    return 0;
}