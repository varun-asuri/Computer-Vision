#include <iostream> 
#include <fstream>
#include <cstdlib>
#include <string>
#include <cstdio>
#include <math.h>
#include <ctime>

using namespace std;

class Point {
    public:
        double x;
        double y;
        Point(double valx, double valy) {
            x = valx;
            y = valy;
        }
};

int draw[800][800][3];
ofstream ppm("output.ppm");

bool valid(int a) {
    if(a >= 0 && a < 800)
        return true;
    return false;
}

double distance(double x1, double y1, double x2, double y2)  { 
    return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2)); 
}

double areaTriangle(double x1, double y1, double x2, double y2, double x3, double y3) {
    return abs(x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2))/2.0;
}

void part1() { 
    srand(time(NULL));
    ofstream outfile("points.txt");
    ppm << "P3\n800 800\n255" << endl;

    int point1[2] = {rand() % 800, rand() % 800};
    int point2[2] = {rand() % 800, rand() % 800};
    int point3[2] = {rand() % 800, rand() % 800};
    int point4[2] = {rand() % 800, rand() % 800};
    double A = areaTriangle(point1[0], point1[1], point2[0], point2[1], point3[0], point3[1]);
    double A1 = areaTriangle(point1[0], point1[1], point2[0], point2[1], point4[0], point4[1]);
    double A2 = areaTriangle(point1[0], point1[1], point4[0], point4[1], point3[0], point3[1]);
    double A3 = areaTriangle(point4[0], point4[1], point2[0], point2[1], point3[0], point3[1]);
    while((A1+A2+A3) == A) {
        point4[0] = rand() % 800;
        point4[1] = rand() % 800;
        A1 = areaTriangle(point1[0], point1[1], point2[0], point2[1], point4[0], point4[1]);
        A2 = areaTriangle(point1[0], point1[1], point4[0], point4[1], point3[0], point3[1]);
        A3 = areaTriangle(point4[0], point4[1], point2[0], point2[1], point3[0], point3[1]);
    }
    outfile << "(" << ((double)point1[0])/800.0 << "," << ((double)point1[1])/800.0 << ") , (" << ((double)point2[0])/800.0 << "," << ((double)point2[1])/800.0 << ") , (" << ((double)point3[0])/800.0 << "," << ((double)point3[1])/800.0 << ") , (" << ((double)point4[0])/800.0 << "," << ((double)point4[1])/800.0 << ")" << endl; 

    double dx;
    double dy;
    int x0;
    int x1;
    int y0;
    int y1;

    x0 = point1[0];
    x1 = point2[0];
    y0 = point1[1];
    y1 = point2[1];

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
            draw[x0][y0][2] += 255;
            if(e > 0)
            {
                if(dy < 0)
                {
                    y0 -= 1;
                }
                else
                {
                    y0 += 1;
                }
                e -= 1;
            }
            e += abs(r);
            if(dx < 0)
            {
                x0 -= 1;
            }
            else
            {
                x0 += 1;
            }
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
            draw[x0][y0][2] += 255;
            if(e > 0)
            {
                if(dx < 0)
                {
                    x0 -= 1;
                }
                else
                {
                    x0 += 1;
                }
                e -= 1;
            }
            e += abs(r);
            if(dy < 0)
            {
                y0 -= 1;
            }
            else
            {
                y0 += 1;
            }
        }
    }
    
    x0 = point2[0];
    x1 = point3[0];
    y0 = point2[1];
    y1 = point3[1];

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
            draw[x0][y0][2] += 255;
            if(e > 0)
            {
                if(dy < 0)
                {
                    y0 -= 1;
                }
                else
                {
                    y0 += 1;
                }
                e -= 1;
            }
            e += abs(r);
            if(dx < 0)
            {
                x0 -= 1;
            }
            else
            {
                x0 += 1;
            }
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
            draw[x0][y0][2] += 255;
            if(e > 0)
            {
                if(dx < 0)
                {
                    x0 -= 1;
                }
                else
                {
                    x0 += 1;
                }
                e -= 1;
            }
            e += abs(r);
            if(dy < 0)
            {
                y0 -= 1;
            }
            else
            {
                y0 += 1;
            }
        }
    }

    x0 = point3[0];
    x1 = point1[0];
    y0 = point3[1];
    y1 = point1[1];

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
            draw[x0][y0][2] += 255;
            if(e > 0)
            {
                if(dy < 0)
                {
                    y0 -= 1;
                }
                else
                {
                    y0 += 1;
                }
                e -= 1;
            }
            e += abs(r);
            if(dx < 0)
            {
                x0 -= 1;
            }
            else
            {
                x0 += 1;
            }
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
            draw[x0][y0][2] += 255;
            if(e > 0)
            {
                if(dx < 0)
                {
                    x0 -= 1;
                }
                else
                {
                    x0 += 1;
                }
                e -= 1;
            }
            e += abs(r);
            if(dy < 0)
            {
                y0 -= 1;
            }
            else
            {
                y0 += 1;
            }
        }
    }

    draw[point4[0]][point4[1]][0] += 255;
    draw[point4[0]][point4[1]][1] += 255;
    draw[point4[0]][point4[1]][2] += 255;
    
    for(int i = 0; i < 800; i++)
    {
        for(int j = 0; j < 800; j++)
        {
            ppm << draw[i][j][0] << " " << draw[i][j][1] << " " << draw[i][j][2] << "\t";
        }
        ppm << "\n";
    }

    ppm.close();
}

void part2() {
    srand(time(NULL));
    ifstream infile("points.txt");
    ofstream outfile("output.txt");
    ppm << "P3\n800 800\n255" << endl;

    int comma;
    string junk;
    string point1;
    string point2;
    string point3;
    string point4;

    infile >> point1;
    infile >> junk;
    infile >> point2;
    infile >> junk;
    infile >> point3;
    infile >> junk;
    infile >> point4;

    comma = point1.find(",");
    double pointA[2] = { stof(point1.substr(1, comma-1)), stof(point1.substr(comma+1, point1.length()-comma-2)) };
    comma = point2.find(",");
    double pointB[2] = { stof(point2.substr(1, comma-1)), stof(point2.substr(comma+1, point2.length()-comma-2)) };
    comma = point3.find(",");
    double pointC[2] = { stof(point3.substr(1, comma-1)), stof(point3.substr(comma+1, point3.length()-comma-2)) };
    comma = point4.find(",");
    double pointD[2] = { stof(point4.substr(1, comma-1)), stof(point4.substr(comma+1, point4.length()-comma-2)) };

    outfile << "(" << pointA[0] << "," << pointA[1] << ") , (" << pointB[0] << "," << pointB[1] << ") , (" << pointC[0] << "," << pointC[1] << ") , (" << pointD[0] << "," << pointD[1] << ")\n";

    int a0, b0, amax, b1, b1_new, ty, x, y, r;

    x = (int) (pointA[0]*800);
    y = (int) (pointA[1]*800);
    r = 2;
    amax = (int) (r * 2*acos(0.0) / 4);
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
        }
        if(valid(a0+x) && valid(-b0+y))
        {
            draw[a0+x][-b0+y][0] += 255;
        }
        if(valid(-a0+x) && valid(b0+y)) 
        {
            draw[-a0+x][b0+y][0] += 255;
        }
        if(valid(-a0+x) && valid(-b0+y))   
        {
            draw[-a0+x][-b0+y][0] += 255;
        }
        if(valid(a0+y) && valid(b0+x))
        {
            draw[b0+x][a0+y][0] += 255;
        }
        if(valid(-a0+y) && valid(b0+x))
        {
            draw[b0+x][-a0+y][0] += 255;
        }
        if(valid(a0+y) && valid(-b0+x))
        {
            draw[-b0+x][a0+y][0] += 255;
        }
        if(valid(-a0+y) && valid(-b0+x))
        {
            draw[-b0+x][-a0+y][0] += 255;
        }
        b1_new -= (2 * a0) - 3;
    }

    x = (int) (pointB[0]*800);
    y = (int) (pointB[1]*800);
    r = 2;
    amax = (int) (r * 2*acos(0.0) / 4);
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
        }
        if(valid(a0+x) && valid(-b0+y))
        {
            draw[a0+x][-b0+y][0] += 255;
        }
        if(valid(-a0+x) && valid(b0+y)) 
        {
            draw[-a0+x][b0+y][0] += 255;
        }
        if(valid(-a0+x) && valid(-b0+y))   
        {
            draw[-a0+x][-b0+y][0] += 255;
        }
        if(valid(a0+y) && valid(b0+x))
        {
            draw[b0+x][a0+y][0] += 255;
        }
        if(valid(-a0+y) && valid(b0+x))
        {
            draw[b0+x][-a0+y][0] += 255;
        }
        if(valid(a0+y) && valid(-b0+x))
        {
            draw[-b0+x][a0+y][0] += 255;
        }
        if(valid(-a0+y) && valid(-b0+x))
        {
            draw[-b0+x][-a0+y][0] += 255;
        }
        b1_new -= (2 * a0) - 3;
    }

    x = (int) (pointC[0]*800);
    y = (int) (pointC[1]*800);
    r = 2;
    amax = (int) (r * 2*acos(0.0) / 4);
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
        }
        if(valid(a0+x) && valid(-b0+y))
        {
            draw[a0+x][-b0+y][0] += 255;
        }
        if(valid(-a0+x) && valid(b0+y)) 
        {
            draw[-a0+x][b0+y][0] += 255;
        }
        if(valid(-a0+x) && valid(-b0+y))   
        {
            draw[-a0+x][-b0+y][0] += 255;
        }
        if(valid(a0+y) && valid(b0+x))
        {
            draw[b0+x][a0+y][0] += 255;
        }
        if(valid(-a0+y) && valid(b0+x))
        {
            draw[b0+x][-a0+y][0] += 255;
        }
        if(valid(a0+y) && valid(-b0+x))
        {
            draw[-b0+x][a0+y][0] += 255;
        }
        if(valid(-a0+y) && valid(-b0+x))
        {
            draw[-b0+x][-a0+y][0] += 255;
        }
        b1_new -= (2 * a0) - 3;
    }

    x = (int) (pointD[0]*800);
    y = (int) (pointD[1]*800);
    r = 2;
    amax = (int) (r * 2*acos(0.0) / 4);
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
        }
        if(valid(a0+x) && valid(-b0+y))
        {
            draw[a0+x][-b0+y][0] += 255;
        }
        if(valid(-a0+x) && valid(b0+y)) 
        {
            draw[-a0+x][b0+y][0] += 255;
        }
        if(valid(-a0+x) && valid(-b0+y))   
        {
            draw[-a0+x][-b0+y][0] += 255;
        }
        if(valid(a0+y) && valid(b0+x))
        {
            draw[b0+x][a0+y][0] += 255;
        }
        if(valid(-a0+y) && valid(b0+x))
        {
            draw[b0+x][-a0+y][0] += 255;
        }
        if(valid(a0+y) && valid(-b0+x))
        {
            draw[-b0+x][a0+y][0] += 255;
        }
        if(valid(-a0+y) && valid(-b0+x))
        {
            draw[-b0+x][-a0+y][0] += 255;
        }
        b1_new -= (2 * a0) - 3;
    }
    
    double areas[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double square[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double smallest[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double measures[16] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double Ax;
    double Ay;
    double Bx;
    double By;
    double Cx;
    double Cy;
    double Dx;
    double Dy;
    double Ex;
    double Ey;
    double m;
    double m_raw;
    double d;
    double dx;
    double dy;
    double A_i;
    double B_i;
    double C_i;
    double D_i;
    double AC_s;
    double BD_s;
    int minimum = 1;

    Ax = pointA[0];
    Ay = pointA[1];
    Bx = pointB[0];
    By = pointB[1];
    Cx = pointC[0];
    Cy = pointC[1];
    Dx = pointD[0];
    Dy = pointD[1];

    m_raw = (Cy - Ay)/(Cx - Ax);
    d = distance(Ax, Ay, Cx, Cy);
    m = -1/m_raw;
    dx = d / sqrt(1 + (m*m));
    dy = (d*m) / sqrt(1 + (m*m));

    Ex = Bx + dx;
    Ey = By + dy;

    BD_s = (Ey - Dy)/(Ex - Dx);
    D_i = Ey-(BD_s*Ex);
    AC_s = -1/BD_s;
    A_i = Ay-(AC_s*Ax);
    C_i = Cy-(AC_s*Cx);
    B_i = By-(BD_s*Bx);

    square[0] = (B_i-A_i)/(AC_s-BD_s);
    square[1] = AC_s*square[0]+A_i;
    square[2] = (C_i-B_i)/(BD_s-AC_s);
    square[3] = BD_s*square[2]+B_i;
    square[4] = (D_i-C_i)/(AC_s-BD_s);
    square[5] = AC_s*square[4]+C_i;
    square[6] = (A_i-D_i)/(BD_s-AC_s);
    square[7] = BD_s*square[6]+D_i;

    outfile << "(" << square[0] << "," << square[1] << ") , (" << square[2] << "," << square[3] << ") , (" << square[4] << "," << square[5] << ") , (" << square[6] << "," << square[7] << ")\n";
    areas[0] = distance(square[0], square[1], square[2], square[3])*distance(square[0], square[1], square[2], square[3]);
    outfile << "Area: " << areas[0] << "\n";
    if(areas[0] < minimum) {
        minimum = areas[0];
        for(int i = 0; i < 8; i++) {
            smallest[i] = square[i];
        }
        measures[0] = AC_s;
        measures[1] = BD_s;
        measures[2] = A_i;
        measures[3] = B_i;
        measures[4] = C_i;
        measures[5] = D_i;
        measures[6] = Ax;
        measures[7] = Ay;
        measures[8] = Bx;
        measures[9] = By;
        measures[10] = Cx;
        measures[11] = Cy;
        measures[12] = Dx;
        measures[13] = Dy;
        measures[14] = Ex;
        measures[15] = Ey;
    }
    
    Ax = pointA[0];
    Ay = pointA[1];
    Bx = pointB[0];
    By = pointB[1];
    Cx = pointD[0];
    Cy = pointD[1];
    Dx = pointC[0];
    Dy = pointC[1];
    
    m_raw = (Cy - Ay)/(Cx - Ax);
    d = distance(Ax, Ay, Cx, Cy);
    m = -1/m_raw;
    dx = d / sqrt(1 + (m*m));
    dy = (d*m) / sqrt(1 + (m*m));

    Ex = Bx + dx;
    Ey = By + dy;

    BD_s = (Ey - Dy)/(Ex - Dx);
    D_i = Ey-(BD_s*Ex);
    AC_s = -1/BD_s;
    A_i = Ay-(AC_s*Ax);
    C_i = Cy-(AC_s*Cx);
    B_i = By-(BD_s*Bx);

    square[0] = (B_i-A_i)/(AC_s-BD_s);
    square[1] = AC_s*square[0]+A_i;
    square[2] = (C_i-B_i)/(BD_s-AC_s);
    square[3] = BD_s*square[2]+B_i;
    square[4] = (D_i-C_i)/(AC_s-BD_s);
    square[5] = AC_s*square[4]+C_i;
    square[6] = (A_i-D_i)/(BD_s-AC_s);
    square[7] = BD_s*square[6]+D_i;

    outfile << "(" << square[0] << "," << square[1] << ") , (" << square[2] << "," << square[3] << ") , (" << square[4] << "," << square[5] << ") , (" << square[6] << "," << square[7] << ")\n";
    areas[1] = distance(square[0], square[1], square[2], square[3])*distance(square[0], square[1], square[2], square[3]);
    outfile << "Area: " << areas[1] << "\n";
    if(areas[1] < minimum) {
        minimum = areas[1];
        for(int i = 0; i < 8; i++) {
            smallest[i] = square[i];
        }
        measures[0] = AC_s;
        measures[1] = BD_s;
        measures[2] = A_i;
        measures[3] = B_i;
        measures[4] = C_i;
        measures[5] = D_i;
        measures[6] = Ax;
        measures[7] = Ay;
        measures[8] = Bx;
        measures[9] = By;
        measures[10] = Cx;
        measures[11] = Cy;
        measures[12] = Dx;
        measures[13] = Dy;
        measures[14] = Ex;
        measures[15] = Ey;
    }
    
    Ax = pointA[0];
    Ay = pointA[1];
    Bx = pointC[0];
    By = pointC[1];
    Cx = pointB[0];
    Cy = pointB[1];
    Dx = pointD[0];
    Dy = pointD[1];

    m_raw = (Cy - Ay)/(Cx - Ax);
    d = distance(Ax, Ay, Cx, Cy);
    m = -1/m_raw;
    dx = d / sqrt(1 + (m*m));
    dy = (d*m) / sqrt(1 + (m*m));

    Ex = Bx + dx;
    Ey = By + dy;

    BD_s = (Ey - Dy)/(Ex - Dx);
    D_i = Ey-(BD_s*Ex);
    AC_s = -1/BD_s;
    A_i = Ay-(AC_s*Ax);
    C_i = Cy-(AC_s*Cx);
    B_i = By-(BD_s*Bx);

    square[0] = (B_i-A_i)/(AC_s-BD_s);
    square[1] = AC_s*square[0]+A_i;
    square[2] = (C_i-B_i)/(BD_s-AC_s);
    square[3] = BD_s*square[2]+B_i;
    square[4] = (D_i-C_i)/(AC_s-BD_s);
    square[5] = AC_s*square[4]+C_i;
    square[6] = (A_i-D_i)/(BD_s-AC_s);
    square[7] = BD_s*square[6]+D_i;

    outfile << "(" << square[0] << "," << square[1] << ") , (" << square[2] << "," << square[3] << ") , (" << square[4] << "," << square[5] << ") , (" << square[6] << "," << square[7] << ")\n";
    areas[2] = distance(square[0], square[1], square[2], square[3])*distance(square[0], square[1], square[2], square[3]);
    outfile << "Area: " << areas[2] << "\n";
    if(areas[2] < minimum) {
        minimum = areas[2];
        for(int i = 0; i < 8; i++) {
            smallest[i] = square[i];
        }
        measures[0] = AC_s;
        measures[1] = BD_s;
        measures[2] = A_i;
        measures[3] = B_i;
        measures[4] = C_i;
        measures[5] = D_i;
        measures[6] = Ax;
        measures[7] = Ay;
        measures[8] = Bx;
        measures[9] = By;
        measures[10] = Cx;
        measures[11] = Cy;
        measures[12] = Dx;
        measures[13] = Dy;
        measures[14] = Ex;
        measures[15] = Ey;
    }
    
    Ax = pointA[0];
    Ay = pointA[1];
    Bx = pointC[0];
    By = pointC[1];
    Cx = pointD[0];
    Cy = pointD[1];
    Dx = pointB[0];
    Dy = pointB[1];

    m_raw = (Cy - Ay)/(Cx - Ax);
    d = distance(Ax, Ay, Cx, Cy);
    m = -1/m_raw;
    dx = d / sqrt(1 + (m*m));
    dy = (d*m) / sqrt(1 + (m*m));

    Ex = Bx + dx;
    Ey = By + dy;
    
    BD_s = (Ey - Dy)/(Ex - Dx);
    D_i = Ey-(BD_s*Ex);
    AC_s = -1/BD_s;
    A_i = Ay-(AC_s*Ax);
    C_i = Cy-(AC_s*Cx);
    B_i = By-(BD_s*Bx);

    square[0] = (B_i-A_i)/(AC_s-BD_s);
    square[1] = AC_s*square[0]+A_i;
    square[2] = (C_i-B_i)/(BD_s-AC_s);
    square[3] = BD_s*square[2]+B_i;
    square[4] = (D_i-C_i)/(AC_s-BD_s);
    square[5] = AC_s*square[4]+C_i;
    square[6] = (A_i-D_i)/(BD_s-AC_s);
    square[7] = BD_s*square[6]+D_i;

    outfile << "(" << square[0] << "," << square[1] << ") , (" << square[2] << "," << square[3] << ") , (" << square[4] << "," << square[5] << ") , (" << square[6] << "," << square[7] << ")\n";
    areas[3] = distance(square[0], square[1], square[2], square[3])*distance(square[0], square[1], square[2], square[3]);
    outfile << "Area: " << areas[3] << "\n";
    if(areas[3] < minimum) {
        minimum = areas[3];
        for(int i = 0; i < 8; i++) {
            smallest[i] = square[i];
        }
        measures[0] = AC_s;
        measures[1] = BD_s;
        measures[2] = A_i;
        measures[3] = B_i;
        measures[4] = C_i;
        measures[5] = D_i;
        measures[6] = Ax;
        measures[7] = Ay;
        measures[8] = Bx;
        measures[9] = By;
        measures[10] = Cx;
        measures[11] = Cy;
        measures[12] = Dx;
        measures[13] = Dy;
        measures[14] = Ex;
        measures[15] = Ey;
    }
    
    Ax = pointA[0];
    Ay = pointA[1];
    Bx = pointD[0];
    By = pointD[1];
    Cx = pointB[0];
    Cy = pointB[1];
    Dx = pointC[0];
    Dy = pointC[1];

    m_raw = (Cy - Ay)/(Cx - Ax);
    d = distance(Ax, Ay, Cx, Cy);
    m = -1/m_raw;
    dx = d / sqrt(1 + (m*m));
    dy = (d*m) / sqrt(1 + (m*m));

    Ex = Bx + dx;
    Ey = By + dy;
    
    BD_s = (Ey - Dy)/(Ex - Dx);
    D_i = Ey-(BD_s*Ex);
    AC_s = -1/BD_s;
    A_i = Ay-(AC_s*Ax);
    C_i = Cy-(AC_s*Cx);
    B_i = By-(BD_s*Bx);

    square[0] = (B_i-A_i)/(AC_s-BD_s);
    square[1] = AC_s*square[0]+A_i;
    square[2] = (C_i-B_i)/(BD_s-AC_s);
    square[3] = BD_s*square[2]+B_i;
    square[4] = (D_i-C_i)/(AC_s-BD_s);
    square[5] = AC_s*square[4]+C_i;
    square[6] = (A_i-D_i)/(BD_s-AC_s);
    square[7] = BD_s*square[6]+D_i;

    outfile << "(" << square[0] << "," << square[1] << ") , (" << square[2] << "," << square[3] << ") , (" << square[4] << "," << square[5] << ") , (" << square[6] << "," << square[7] << ")\n";
    areas[4] = distance(square[0], square[1], square[2], square[3])*distance(square[0], square[1], square[2], square[3]);
    outfile << "Area: " << areas[4] << "\n";
    if(areas[4] < minimum) {
        minimum = areas[4];
        for(int i = 0; i < 8; i++) {
            smallest[i] = square[i];
        }
        measures[0] = AC_s;
        measures[1] = BD_s;
        measures[2] = A_i;
        measures[3] = B_i;
        measures[4] = C_i;
        measures[5] = D_i;
        measures[6] = Ax;
        measures[7] = Ay;
        measures[8] = Bx;
        measures[9] = By;
        measures[10] = Cx;
        measures[11] = Cy;
        measures[12] = Dx;
        measures[13] = Dy;
        measures[14] = Ex;
        measures[15] = Ey;
    }
    
    Ax = pointA[0];
    Ay = pointA[1];
    Bx = pointD[0];
    By = pointD[1];
    Cx = pointC[0];
    Cy = pointC[1];
    Dx = pointB[0];
    Dy = pointB[1];

    m_raw = (Cy - Ay)/(Cx - Ax);
    d = distance(Ax, Ay, Cx, Cy);
    m = -1/m_raw;
    dx = d / sqrt(1 + (m*m));
    dy = (d*m) / sqrt(1 + (m*m));

    Ex = Bx + dx;
    Ey = By + dy;

    BD_s = (Ey - Dy)/(Ex - Dx);
    D_i = Ey-(BD_s*Ex);
    AC_s = -1/BD_s;
    A_i = Ay-(AC_s*Ax);
    C_i = Cy-(AC_s*Cx);
    B_i = By-(BD_s*Bx);

    square[0] = (B_i-A_i)/(AC_s-BD_s);
    square[1] = AC_s*square[0]+A_i;
    square[2] = (C_i-B_i)/(BD_s-AC_s);
    square[3] = BD_s*square[2]+B_i;
    square[4] = (D_i-C_i)/(AC_s-BD_s);
    square[5] = AC_s*square[4]+C_i;
    square[6] = (A_i-D_i)/(BD_s-AC_s);
    square[7] = BD_s*square[6]+D_i;

    outfile << "(" << square[0] << "," << square[1] << ") , (" << square[2] << "," << square[3] << ") , (" << square[4] << "," << square[5] << ") , (" << square[6] << "," << square[7] << ")\n";
    areas[5] = distance(square[0], square[1], square[2], square[3])*distance(square[0], square[1], square[2], square[3]);
    outfile << "Area: " << areas[5] << "\n";
    if(areas[5] < minimum) {
        minimum = areas[5];
        for(int i = 0; i < 8; i++) {
            smallest[i] = square[i];
        }
        measures[0] = AC_s;
        measures[1] = BD_s;
        measures[2] = A_i;
        measures[3] = B_i;
        measures[4] = C_i;
        measures[5] = D_i;
        measures[6] = Ax;
        measures[7] = Ay;
        measures[8] = Bx;
        measures[9] = By;
        measures[10] = Cx;
        measures[11] = Cy;
        measures[12] = Dx;
        measures[13] = Dy;
        measures[14] = Ex;
        measures[15] = Ey;
    }

    // x = (int) (measures[14]*800);
    // y = (int) (measures[15]*800);
    // r = 2;
    // amax = (int) (r * 2*acos(0.0) / 4);
    // b0 = r;
    // b1 = b0 * b0;
    // ty = (2 * b0) - 1;
    // b1_new = b1;
    // for (a0 = 0; a0 <= amax; a0++)
    // {
    //     if ((b1 - b1_new) >= ty)
    //     {
    //         b1 -= ty;
    //         b0 -= 1;
    //         ty -= 2;
    //     }
    //     if(valid(a0+x) && valid(b0+y))
    //     {
    //         draw[a0+x][b0+y][1] += 255;
    //     }
    //     if(valid(a0+x) && valid(-b0+y))
    //     {
    //         draw[a0+x][-b0+y][1] += 255;
    //     }
    //     if(valid(-a0+x) && valid(b0+y)) 
    //     {
    //         draw[-a0+x][b0+y][1] += 255;
    //     }
    //     if(valid(-a0+x) && valid(-b0+y))   
    //     {
    //         draw[-a0+x][-b0+y][1] += 255;
    //     }
    //     if(valid(a0+y) && valid(b0+x))
    //     {
    //         draw[b0+x][a0+y][1] += 255;
    //     }
    //     if(valid(-a0+y) && valid(b0+x))
    //     {
    //         draw[b0+x][-a0+y][1] += 255;
    //     }
    //     if(valid(a0+y) && valid(-b0+x))
    //     {
    //         draw[-b0+x][a0+y][1] += 255;
    //     }
    //     if(valid(-a0+y) && valid(-b0+x))
    //     {
    //         draw[-b0+x][-a0+y][1] += 255;
    //     }
    //     b1_new -= (2 * a0) - 3;
    // }

    int x0;
    int x1;
    int y0;
    int y1;

    double x0temp;
    double y0temp;
    double x1temp;
    double y1temp;
    double h;
    
    // x0 = (int) (measures[6]*800);
    // y0 = (int) (measures[7]*800);
    // x1 = (int) (measures[10]*800);
    // y1 = (int) (measures[11]*800);

    // dx = x1 - x0;
    // dy = y1 - y0;
    
    // if (abs(dx) >= abs(dy))
    // {
    //     double e = 0;
    //     double r = dy/dx;
    //     while(x0 != x1)
    //     {
    //         draw[x0][y0][0] += 255;
    //         draw[x0][y0][1] += 255;
    //         draw[x0][y0][2] += 255;
    //         if(e > 0)
    //         {
    //             if(dy < 0)
    //             {
    //                 y0 -= 1;
    //             }
    //             else
    //             {
    //                 y0 += 1;
    //             }
    //             e -= 1;
    //         }
    //         e += abs(r);
    //         if(dx < 0)
    //         {
    //             x0 -= 1;
    //         }
    //         else
    //         {
    //             x0 += 1;
    //         }
    //     }
    // }
    // else
    // {
    //     double e = 0;
    //     double r = dx/dy;
    //     while(y0 != y1)
    //     {
    //         draw[x0][y0][0] += 255;
    //         draw[x0][y0][1] += 255;
    //         draw[x0][y0][2] += 255;
    //         if(e > 0)
    //         {
    //             if(dx < 0)
    //             {
    //                 x0 -= 1;
    //             }
    //             else
    //             {
    //                 x0 += 1;
    //             }
    //             e -= 1;
    //         }
    //         e += abs(r);
    //         if(dy < 0)
    //         {
    //             y0 -= 1;
    //         }
    //         else
    //         {
    //             y0 += 1;
    //         }
    //     }
    // }
    
    // x0 = (int) (measures[8]*800);
    // y0 = (int) (measures[9]*800);
    // x1 = (int) (measures[14]*800);
    // y1 = (int) (measures[15]*800);

    // dx = x1 - x0;
    // dy = y1 - y0;
    
    // if (abs(dx) >= abs(dy))
    // {
    //     double e = 0;
    //     double r = dy/dx;
    //     while(x0 != x1)
    //     {
    //         draw[x0][y0][0] += 255;
    //         draw[x0][y0][1] += 255;
    //         draw[x0][y0][2] += 255;
    //         if(e > 0)
    //         {
    //             if(dy < 0)
    //             {
    //                 y0 -= 1;
    //             }
    //             else
    //             {
    //                 y0 += 1;
    //             }
    //             e -= 1;
    //         }
    //         e += abs(r);
    //         if(dx < 0)
    //         {
    //             x0 -= 1;
    //         }
    //         else
    //         {
    //             x0 += 1;
    //         }
    //     }
    // }
    // else
    // {
    //     double e = 0;
    //     double r = dx/dy;
    //     while(y0 != y1)
    //     {
    //         draw[x0][y0][0] += 255;
    //         draw[x0][y0][1] += 255;
    //         draw[x0][y0][2] += 255;
    //         if(e > 0)
    //         {
    //             if(dx < 0)
    //             {
    //                 x0 -= 1;
    //             }
    //             else
    //             {
    //                 x0 += 1;
    //             }
    //             e -= 1;
    //         }
    //         e += abs(r);
    //         if(dy < 0)
    //         {
    //             y0 -= 1;
    //         }
    //         else
    //         {
    //             y0 += 1;
    //         }
    //     }
    // }
    
    x0temp = smallest[0];
    y0temp = smallest[1];
    x1temp = smallest[2];
    y1temp = smallest[3];
    
    dx = x1temp - x0temp;
    dy = y1temp - y0temp;
    h = sqrt((dx*dx)+(dy*dy));

    while(x0temp >= (dx/h) && x0temp < (1+(dx/h)) && y0temp > (dy/h) && y0temp < (1+(dx/h))) {
        x0temp = x0temp - (dx/h);
        y0temp = y0temp - (dy/h);
    }
    while(x1temp >= (-dx/h) && x1temp < (1-(dx/h)) && y1temp > (-dy/h) && y1temp < (1-(dx/h))) {
        x1temp = x1temp + (dx/h);
        y1temp = y1temp + (dy/h);
    }
    
    x0 = (int) (x0temp*800);
    y0 = (int) (y0temp*800);
    x1 = (int) (x1temp*800);
    y1 = (int) (y1temp*800);

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
            draw[x0][y0][2] += 255;
            if(e > 0)
            {
                if(dy < 0)
                {
                    y0 -= 1;
                }
                else
                {
                    y0 += 1;
                }
                e -= 1;
            }
            e += abs(r);
            if(dx < 0)
            {
                x0 -= 1;
            }
            else
            {
                x0 += 1;
            }
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
            draw[x0][y0][2] += 255;
            if(e > 0)
            {
                if(dx < 0)
                {
                    x0 -= 1;
                }
                else
                {
                    x0 += 1;
                }
                e -= 1;
            }
            e += abs(r);
            if(dy < 0)
            {
                y0 -= 1;
            }
            else
            {
                y0 += 1;
            }
        }
    }

    x0temp = smallest[2];
    y0temp = smallest[3];
    x1temp = smallest[4];
    y1temp = smallest[5];
    
    dx = x1temp - x0temp;
    dy = y1temp - y0temp;
    h = (dx*dx)+(dy*dy);

    while(x0temp >= (dx/h) && x0temp < (1+(dx/h)) && y0temp > (dy/h) && y0temp < (1+(dx/h))) {
        x0temp = x0temp - (dx/h);
        y0temp = y0temp - (dy/h);
    }
    while(x1temp >= (-dx/h) && x1temp < (1-(dx/h)) && y1temp > (-dy/h) && y1temp < (1-(dx/h))) {
        x1temp = x1temp + (dx/h);
        y1temp = y1temp + (dy/h);
    }
    
    x0 = (int) (x0temp*800);
    y0 = (int) (y0temp*800);
    x1 = (int) (x1temp*800);
    y1 = (int) (y1temp*800);

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
            draw[x0][y0][2] += 255;
            if(e > 0)
            {
                if(dy < 0)
                {
                    y0 -= 1;
                }
                else
                {
                    y0 += 1;
                }
                e -= 1;
            }
            e += abs(r);
            if(dx < 0)
            {
                x0 -= 1;
            }
            else
            {
                x0 += 1;
            }
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
            draw[x0][y0][2] += 255;
            if(e > 0)
            {
                if(dx < 0)
                {
                    x0 -= 1;
                }
                else
                {
                    x0 += 1;
                }
                e -= 1;
            }
            e += abs(r);
            if(dy < 0)
            {
                y0 -= 1;
            }
            else
            {
                y0 += 1;
            }
        }
    }
    
    x0temp = smallest[4];
    y0temp = smallest[5];
    x1temp = smallest[6];
    y1temp = smallest[7];
    
    dx = x1temp - x0temp;
    dy = y1temp - y0temp;
    h = (dx*dx)+(dy*dy);

    while(x0temp >= (dx/h) && x0temp < (1+(dx/h)) && y0temp > (dy/h) && y0temp < (1+(dx/h))) {
        x0temp = x0temp - (dx/h);
        y0temp = y0temp - (dy/h);
    }
    while(x1temp >= (-dx/h) && x1temp < (1-(dx/h)) && y1temp > (-dy/h) && y1temp < (1-(dx/h))) {
        x1temp = x1temp + (dx/h);
        y1temp = y1temp + (dy/h);
    }
    
    x0 = (int) (x0temp*800);
    y0 = (int) (y0temp*800);
    x1 = (int) (x1temp*800);
    y1 = (int) (y1temp*800);

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
            draw[x0][y0][2] += 255;
            if(e > 0)
            {
                if(dy < 0)
                {
                    y0 -= 1;
                }
                else
                {
                    y0 += 1;
                }
                e -= 1;
            }
            e += abs(r);
            if(dx < 0)
            {
                x0 -= 1;
            }
            else
            {
                x0 += 1;
            }
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
            draw[x0][y0][2] += 255;
            if(e > 0)
            {
                if(dx < 0)
                {
                    x0 -= 1;
                }
                else
                {
                    x0 += 1;
                }
                e -= 1;
            }
            e += abs(r);
            if(dy < 0)
            {
                y0 -= 1;
            }
            else
            {
                y0 += 1;
            }
        }
    }
    
    x0temp = smallest[6];
    y0temp = smallest[7];
    x1temp = smallest[0];
    y1temp = smallest[1];
    
    dx = x1temp - x0temp;
    dy = y1temp - y0temp;
    h = (dx*dx)+(dy*dy);

    while(x0temp >= (dx/h) && x0temp < (1+(dx/h)) && y0temp > (dy/h) && y0temp < (1+(dx/h))) {
        x0temp = x0temp - (dx/h);
        y0temp = y0temp - (dy/h);
    }
    while(x1temp >= (-dx/h) && x1temp < (1-(dx/h)) && y1temp > (-dy/h) && y1temp < (1-(dx/h))) {
        x1temp = x1temp + (dx/h);
        y1temp = y1temp + (dy/h);
    }
    
    x0 = (int) (x0temp*800);
    y0 = (int) (y0temp*800);
    x1 = (int) (x1temp*800);
    y1 = (int) (y1temp*800);

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
            draw[x0][y0][2] += 255;
            if(e > 0)
            {
                if(dy < 0)
                {
                    y0 -= 1;
                }
                else
                {
                    y0 += 1;
                }
                e -= 1;
            }
            e += abs(r);
            if(dx < 0)
            {
                x0 -= 1;
            }
            else
            {
                x0 += 1;
            }
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
            draw[x0][y0][2] += 255;
            if(e > 0)
            {
                if(dx < 0)
                {
                    x0 -= 1;
                }
                else
                {
                    x0 += 1;
                }
                e -= 1;
            }
            e += abs(r);
            if(dy < 0)
            {
                y0 -= 1;
            }
            else
            {
                y0 += 1;
            }
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
    //part1();
    part2();
    return 0;
}