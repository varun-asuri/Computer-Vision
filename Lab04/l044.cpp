#include <algorithm>
#include <iostream> 
#include <fstream>
#include <cstdlib>
#include <utility>
#include <string>
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
    private:
        double x;
        double y;
    public:
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
class Pixel {
    public:
        int r;
        int g;
        int b;
        int c;
        Pixel() {
            r = 0;
            g = 0;
            b = 0;
            c = 0;
        }
        Pixel(int valr, int valg, int valb, int valc) {
            r = valr;
            g = valg;
            b = valb;
            c = valc;
        }
        int getr() {
            return 255-r;
        }
        int getg() {
            return 255-g;
        }
        int getb() {
            return 255-b;
        }
        int getc() {
            return c;
        }
};
struct Element { 
    Point p;
    struct Element* left;
    struct Element* right;
    Element(Point valp) { 
        p = valp; 
        left = right = NULL; 
    } 
    double getx() {
        return p.getx();
    }
    double gety() {
        return p.gety();
    }
};

vector<int> kmct;
vector<int> kmot;
vector<double> kmtx;
vector<double> kmty;
vector<Point> kmcd;
map<pair<double, double>, int> kmcr;
int CLUSTERS = 5, POINTS = 100;
int draw[800][800][3];
double xlimtop, xlimbot, ylimtop, ylimbot;
bool ccode = true;
// ofstream debug("debug.txt");

bool valid(int a) {
    return a >= 0 && a < 800;
}
double distance(Point a, Point b) { 
    return sqrt(( (b.getx() - a.getx()) * (b.getx() - a.getx()) ) + ( (b.gety() - a.gety()) * (b.gety() - a.gety()) )); 
}
double difference(Pixel a, Color b) {
    return sqrt(( (b.getr() - a.getr()) * (b.getr() - a.getr()) ) + ( (b.getg() - a.getg()) * (b.getg() - a.getg()) ) + ( (b.getb() - a.getb()) * (b.getb() - a.getb()) )); 
}
Element* insert(Element* pointer, Point p) {
    double tempx, tempy;
    if(pointer == NULL)
        return new Element(p);
    tempx = pointer->p.getx();
    tempy = pointer->p.gety();
    if(ccode) {
        if(p.getx() > tempx) {
            xlimbot = tempx;
            ccode = false;
            pointer->right = insert(pointer->right, p);
        } else {
            xlimtop = tempx;
            ccode = false;
            pointer->left = insert(pointer->left, p);
        }
    } else {
        if(p.gety() > tempy) {
            ylimbot = tempy;
            ccode = true;
            pointer->left = insert(pointer->left, p);
        } else {
            ylimtop = tempy;
            ccode = true;
            pointer->right = insert(pointer->right, p);
        }
    }
    return pointer;
}
void update(Element* pointer, bool iter) {
    if(pointer == NULL)
        return;

    vector<int> dist;
    int i, n, c;
    double tempx, tempy, d=sqrt(2);
    tempx = pointer->p.getx();
    tempy = pointer->p.gety();
    pair<double, double> xy = make_pair(tempx, tempy);

    for(i = 0; i < 5; i++) {
        double temp = distance(pointer->p, kmcd[i]);
        dist.push_back(temp);
        if(temp < d) {
            d = temp;
            n = i;
        }
    }
    // debug << "\n" << n << " - " << tempx << " " << tempy;

    if(iter) {
        c = kmcr[xy];
        // debug << " " << c;
        kmct[c] = kmct[c]-1;
    }
    kmcr[xy] = n;
    kmct[n] = kmct[n]+1;
    kmtx[n] += tempx;
    kmty[n] += tempy;

    update(pointer->left, iter);
    update(pointer->right, iter);
}

void part1() {
    srand(time(NULL));
    ofstream ppm("clusters.ppm");
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
void part2() {
    srand(time(NULL));
    ifstream read("input.ppm");
    ofstream ppm("output.ppm");

    string fill;
    int width, height;
    read >> fill;
    read >> fill;
    width = stoi(fill);
    read >> fill;
    height = stoi(fill);
    read >> fill;

    ppm << "P3\n" << width << " " << height << "\n255" << endl;

    int r, g, b;
    vector<Pixel> image;
    while(!read.eof()) {
        read >> fill;
        r = stoi(fill);
        read >> fill;
        g = stoi(fill);
        read >> fill;
        b = stoi(fill);
        image.push_back(Pixel(r, g, b, 0));
    }
    read.close();

    vector<Color> centroids;
    centroids.push_back(Color(0, 0, 0));
    centroids.push_back(Color(85, 85, 85));
    centroids.push_back(Color(170, 170, 170));
    centroids.push_back(Color(255, 255, 255));

    map<int, vector<Pixel>> kmeans;
    
    Pixel c1, p;
    Color c2;
    bool condition = true;
    int cnum, pnum, count;
    double d, min, totalr, totalg, totalb, rval, gval, bval;
    vector<int> oldlens(4, 0);
    vector<int> lens(4, 0);

    while(condition) {
        count = 0;
        for(int i = 0; i < 4; i++) {
            kmeans[i].clear();
            kmeans[i].resize(0);
        }
        kmeans.clear();
        lens.clear();
        lens.resize(4);
        for(auto it : image) {
            c1 = it;
            min = 442.0;
            for(int i = 0; i < 4; i++) {
                c2 = centroids[i];
                d = difference(c1, c2);
                if(d < min) {
                    cnum = i;
                    min = d;
                }
            }
            image[count].c = cnum;
            kmeans[cnum].push_back(c1);
            lens[cnum] += 1;
            count++;
        }
        vector<Color> newc;
        for(int i = 0; i < 4; i++) {
            totalr = 0, totalg = 0, totalb = 0;
            pnum = 0;
            for(auto it : kmeans[i]) {
                p = it;           
                totalr += p.r;
                totalg += p.g;
                totalb += p.b;
                pnum++;
            }
            rval = totalr / pnum;
            gval = totalg / pnum;
            bval = totalb / pnum;
            newc.push_back(Color(rval, gval, bval));
        }
        
        condition = false;
        for(int i = 0; i < 4; i++) {
            if(oldlens[i] != lens[i])
                condition = true;
        }
        centroids = newc;
        for(int i = 0; i < 4; i++) {
            oldlens[i] = lens[i];
        }
    }

    for(auto it : image) {
        p = it;
        ppm << centroids[p.c].r << " " << centroids[p.c].g << " " << centroids[p.c].b << " ";
    }
    ppm.close();
}
void part3() {
    srand(time(NULL));
    ofstream ppm("diagram.ppm");
    ppm << "P3\n800 800\n255" << endl;
    
    string response;
    cout << "Would you like to generate points?\t";
    cin >> response;
    transform(response.begin(), response.end(), response.begin(), ::tolower);

    if(response.compare("yes") == 0) {
        ofstream outfile("points.txt");
        double xval, yval;
        for(int i = 0; i < 10; i ++) {
            xval = ((double) rand())/RAND_MAX;
            yval = ((double) rand())/RAND_MAX;
            outfile << xval << "  " << yval << endl;
        }
        outfile.close();
    }

    vector<Point> points;
    ifstream infile("points.txt");
    double xval, yval, oldxval, oldyval;
    string xval_raw, yval_raw, oldxval_raw, oldyval_raw;
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
    infile.close();

    Point curr;
    struct Element* head = NULL;
    double vx, vy, tempx, tempy;
    for(auto& it : points) {
        ccode = true;
        curr = it;
        vx = curr.getx();
        vy = curr.gety();
        xlimtop = 1, xlimbot = 0, ylimtop = 1, ylimbot = 0;

        head = insert(head, curr);

        double dx, dy;
        int x0, x1, y0, y1;
        if(ccode) {
            x0 = (int) (vx*800);
            x1 = (int) (vx*800);
            y0 = (int) (ylimbot*800);
            y1 = (int) (ylimtop*800);

            dx = x1 - x0;
            dy = y1 - y0;
            
            if (abs(dx) >= abs(dy))
            {
                double e = 0;
                double r = dy/dx;
                while(x0 != x1)
                {
                    draw[x0][y0][1] = 255;
                    draw[x0][y0][2] = 255;
                    if(e > 0)
                    {
                        if(dy < 0)
                            y0 -= 1;
                        else
                            y0 += 1;
                        e -= 1;
                    }
                    e += abs(r);
                    if(dx < 0)
                        x0 -= 1;
                    else
                        x0 += 1;
                }
            }
            else
            {
                double e = 0;
                double r = dx/dy;
                while(y0 != y1)
                {
                    draw[x0][y0][1] = 255;
                    draw[x0][y0][2] = 255;
                    if(e > 0)
                    {
                        if(dx < 0)
                            x0 -= 1;
                        else
                            x0 += 1;
                        e -= 1;
                    }
                    e += abs(r);
                    if(dy < 0)
                        y0 -= 1;
                    else
                        y0 += 1;
                }
            }
        }
        else {
            x0 = (int) (xlimbot*800);
            x1 = (int) (xlimtop*800);
            y0 = (int) (vy*800);
            y1 = (int) (vy*800);

            dx = x1 - x0;
            dy = y1 - y0;
            
            if (abs(dx) >= abs(dy))
            {
                double e = 0;
                double r = dy/dx;
                while(x0 != x1)
                {
                    draw[x0][y0][0] += 255;
                    draw[x0][y0][1] += 255;
                    if(e > 0)
                    {
                        if(dy < 0)
                            y0 -= 1;
                        else
                            y0 += 1;
                        e -= 1;
                    }
                    e += abs(r);
                    if(dx < 0)
                        x0 -= 1;
                    else
                        x0 += 1;
                }
            }
            else
            {
                double e = 0;
                double r = dx/dy;
                while(y0 != y1)
                {
                    draw[x0][y0][0] += 255;
                    draw[x0][y0][1] += 255;
                    if(e > 0)
                    {
                        if(dx < 0)
                            x0 -= 1;
                        else
                            x0 += 1;
                        e -= 1;
                    }
                    e += abs(r);
                    if(dy < 0)
                        y0 -= 1;
                    else
                        y0 += 1;
                }
            }
        }

        int a0, b0, b1, b1_new, ty, x, y, r;
        double amax;

        x = (int) (vx*800);
        y = (int) (vy*800);
        r = 2;

        amax = (r * 2*acos(0.0) / 4);
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
                draw[b0+x][-a0+y][0] = 255;
                draw[b0+x][-a0+y][1] = 255;
                draw[b0+x][-a0+y][2] = 255;
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
void part4() {
    srand(time(NULL));
    ofstream ppm("clusters.ppm");
    ppm << "P3\n800 800\n255" << endl;

    vector<Color> colors;
    colors.push_back(Color(255, 0, 0));
    colors.push_back(Color(0, 255, 0));
    colors.push_back(Color(0, 0, 255));
    colors.push_back(Color(255, 0, 255));
    colors.push_back(Color(127, 127, 0));
    
    string response;
    cout << "Would you like to generate points?\t";
    cin >> response;
    transform(response.begin(), response.end(), response.begin(), ::tolower);

    if(response.compare("yes") == 0) {
        ofstream outfile("points.txt");
        double xval, yval;
        for(int i = 0; i < 50; i ++) {
            xval = ((double) rand())/RAND_MAX;
            yval = ((double) rand())/RAND_MAX;
            outfile << xval << "  " << yval << endl;
        }
        outfile.close();
    }

    vector<Point> points;
    ifstream infile("points.txt");
    double xval, yval, oldxval, oldyval;
    string xval_raw, yval_raw, oldxval_raw, oldyval_raw;
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
    infile.close();

    Point curr;
    double vx, vy;
    struct Element* head = NULL;
    for(auto& it : points) {
        curr = it;
        vx = curr.getx();
        vy = curr.gety();
        head = insert(head, curr);
    }

    kmcd.push_back(Point(.2, .2));
    kmcd.push_back(Point(.2, .8));
    kmcd.push_back(Point(.5, .5));
    kmcd.push_back(Point(.8, .2));
    kmcd.push_back(Point(.8, .8));
    for(int a = 0; a < 5; a++) {
        kmct.push_back(0);
        kmot.push_back(0);
        kmtx.push_back(0);
        kmty.push_back(0);
    }

    update(head, false);
    // debug << "\n\nPRE" << endl;

    while(kmct[0] != kmot[0] || kmct[1] != kmot[1] || kmct[2] != kmot[2] || kmct[3] != kmot[3] || kmct[4] != kmot[4]) {
        for(int m = 0; m < 5; m++) {
            // debug << kmct[m] << " - " << kmcd[m].getx() << " " << kmcd[m].gety() << endl;
            kmot[m] = kmct[m];
            kmcd[m] = Point(kmtx[m]/((double)kmct[m]), kmty[m]/((double)kmct[m]));
            kmtx[m] = 0;
            kmty[m] = 0;
        }
        update(head, true);
        // debug << "\n\nPOST" << endl;
    }
    for(int p = 0; p < 5; p++) {
        // debug << kmct[p] << " - " << kmcd[p].getx() << " " << kmcd[p].gety() << endl;
        kmcd[p] = Point(kmtx[p]/((double)kmct[p]), kmty[p]/((double)kmct[p]));
    }

    int a0, b0, b1, b1_new, ty, x, y, r, val;
    double amax;
    pair<double, double> cp;
    for(map<pair<double, double>, int>::iterator it = kmcr.begin(); it != kmcr.end(); ++it) {
        cp = it->first;
        val = it->second;
        // debug << cp.first << " " << cp.second << " " << val << endl;
        x = (int) (cp.first*800);
        y = (int) (cp.second*800);
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
                draw[a0+x][b0+y][0] = colors[val].getr();
                draw[a0+x][b0+y][1] = colors[val].getg();
                draw[a0+x][b0+y][2] = colors[val].getb();
            }
            if(valid(a0+x) && valid(-b0+y))
            {
                draw[a0+x][-b0+y][0] = colors[val].getr();
                draw[a0+x][-b0+y][1] = colors[val].getg();
                draw[a0+x][-b0+y][2] = colors[val].getb();
            }
            if(valid(-a0+x) && valid(b0+y)) 
            {
                draw[-a0+x][b0+y][0] = colors[val].getr();
                draw[-a0+x][b0+y][1] = colors[val].getg();
                draw[-a0+x][b0+y][2] = colors[val].getb();
            }
            if(valid(-a0+x) && valid(-b0+y))   
            {
                draw[-a0+x][-b0+y][0] = colors[val].getr();
                draw[-a0+x][-b0+y][1] = colors[val].getg();
                draw[-a0+x][-b0+y][2] = colors[val].getb();
            }
            if(valid(a0+y) && valid(b0+x))
            {
                draw[b0+x][a0+y][0] = colors[val].getr();
                draw[b0+x][a0+y][1] = colors[val].getg();
                draw[b0+x][a0+y][2] = colors[val].getb();
            }
            if(valid(-a0+y) && valid(b0+x))
            {
                draw[b0+x][-a0+y][0] = colors[val].getr();
                draw[b0+x][-a0+y][1] = colors[val].getg();
                draw[b0+x][-a0+y][2] = colors[val].getb();
            }
            if(valid(a0+y) && valid(-b0+x))
            {
                draw[-b0+x][a0+y][0] = colors[val].getr();
                draw[-b0+x][a0+y][1] = colors[val].getg();
                draw[-b0+x][a0+y][2] = colors[val].getb();
            }
            if(valid(-a0+y) && valid(-b0+x))
            {
                draw[-b0+x][-a0+y][0] = colors[val].getr();
                draw[-b0+x][-a0+y][1] = colors[val].getg();
                draw[-b0+x][-a0+y][2] = colors[val].getb();
            }
            b1_new -= (2 * a0) - 3;
        }
    }
    for(vector<Point>::iterator it = kmcd.begin(); it != kmcd.end(); ++it) {
        curr = *it;
        x = (int) (curr.getx()*800);
        y = (int) (curr.gety()*800);
        r = 5;
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
                draw[b0+x][-a0+y][0] = 255;
                draw[b0+x][-a0+y][1] = 255;
                draw[b0+x][-a0+y][2] = 255;
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

    // debug << "END" << endl;

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

int main() {
    //part1();
    //part2();
    //part3();
    part4();
    return 0;
}