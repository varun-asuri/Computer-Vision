#include <iostream> 
#include <fstream>
#include <cstdlib>
#include <math.h>
#include <ctime>

using namespace std;

int draw[800][800][3];
ofstream ppm("triangle.ppm");

/*************************
/White - Triangle        /
/Red - Nine-Point Circle /
/Green - Incircle        /
/Blue - Circumcircle     /
/Yellow - Euler’s Line   /
/Purple - Orthocircle    /
*************************/

bool valid(double a)
{
    if(a >= 0.0 && a < 800.0)
        return true;
    return false;
}

double distance(double x1, double y1, double x2, double y2) 
{ 
    return sqrt(( (x2 - x1) * (x2 - x1) ) + ( (y2 - y1) * (y2 - y1) )); 
}

void drawLine(int x0, int y0, int x1, int y1, int c)
{
    double dx;
    double dy;
    
    dx = x1 - x0;
    dy = y1 - y0;
    
    if (abs(dx) >= abs(dy))
    {
        double e = 0;
        double r = dy/dx;
        while(x0 != x1)
        {
            if(c == 0) {
                draw[x0][y0][0] += 255;
                draw[x0][y0][1] += 255;
                draw[x0][y0][2] += 255;
            } else if(c == 1) {
                draw[x0][y0][0] += 255;
            } else if(c == 2) {
                draw[x0][y0][1] += 255;
            } else if(c == 3) {
                draw[x0][y0][2] += 255;
            } else if(c == 4) {
                draw[x0][y0][0] += 255;
                draw[x0][y0][1] += 255;
            } else if(c == 5) {
                draw[x0][y0][0] += 255;
                draw[x0][y0][2] += 255;
            } else if(c == 6) {
                draw[x0][y0][1] += 255;
                draw[x0][y0][2] += 255;
            }
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
            if(c == 0) {
                draw[x0][y0][0] += 255;
                draw[x0][y0][1] += 255;
                draw[x0][y0][2] += 255;
            } else if(c == 1) {
                draw[x0][y0][0] += 255;
            } else if(c == 2) {
                draw[x0][y0][1] += 255;
            } else if(c == 3) {
                draw[x0][y0][2] += 255;
            } else if(c == 4) {
                draw[x0][y0][0] += 255;
                draw[x0][y0][1] += 255;
            } else if(c == 5) {
                draw[x0][y0][0] += 255;
                draw[x0][y0][2] += 255;
            } else if(c == 6) {
                draw[x0][y0][1] += 255;
                draw[x0][y0][2] += 255;
            }
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

void drawCircle(int x, int y, int r, int c)
{
    int a0, b0, b1, b1_new, ty;
    double amax;

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
            if(c == 0) {
                draw[a0+x][b0+y][0] += 255;
                draw[a0+x][b0+y][1] += 255;
                draw[a0+x][b0+y][2] += 255;
            } else if(c == 1) {
                draw[a0+x][b0+y][0] += 255;
            } else if(c == 2) {
                draw[a0+x][b0+y][1] += 255;
            } else if(c == 3) {
                draw[a0+x][b0+y][2] += 255;
            } else if(c == 4) {
                draw[a0+x][b0+y][0] += 255;
                draw[a0+x][b0+y][1] += 255;
            } else if(c == 5) {
                draw[a0+x][b0+y][0] += 255;
                draw[a0+x][b0+y][2] += 255;
            } else if(c == 6) {
                draw[a0+x][b0+y][1] += 255;
                draw[a0+x][b0+y][2] += 255;
            }
        }
        if(valid(a0+x) && valid(-b0+y))
        {
            if(c == 0) {
                draw[a0+x][-b0+y][0] += 255;
                draw[a0+x][-b0+y][1] += 255;
                draw[a0+x][-b0+y][2] += 255;
            } else if(c == 1) {
                draw[a0+x][-b0+y][0] += 255;
            } else if(c == 2) {
                draw[a0+x][-b0+y][1] += 255;
            } else if(c == 3) {
                draw[a0+x][-b0+y][2] += 255;
            } else if(c == 4) {
                draw[a0+x][-b0+y][0] += 255;
                draw[a0+x][-b0+y][1] += 255;
            } else if(c == 5) {
                draw[a0+x][-b0+y][0] += 255;
                draw[a0+x][-b0+y][2] += 255;
            } else if(c == 6) {
                draw[a0+x][-b0+y][1] += 255;
                draw[a0+x][-b0+y][2] += 255;
            }
        }
        if(valid(-a0+x) && valid(b0+y)) 
        {
            if(c == 0) {
                draw[-a0+x][b0+y][0] += 255;
                draw[-a0+x][b0+y][1] += 255;
                draw[-a0+x][b0+y][2] += 255;
            } else if(c == 1) {
                draw[-a0+x][b0+y][0] += 255;
            } else if(c == 2) {
                draw[-a0+x][b0+y][1] += 255;
            } else if(c == 3) {
                draw[-a0+x][b0+y][2] += 255;
            } else if(c == 4) {
                draw[-a0+x][b0+y][0] += 255;
                draw[-a0+x][b0+y][1] += 255;
            } else if(c == 5) {
                draw[-a0+x][b0+y][0] += 255;
                draw[-a0+x][b0+y][2] += 255;
            } else if(c == 6) {
                draw[-a0+x][b0+y][1] += 255;
                draw[-a0+x][b0+y][2] += 255;
            }
        }
        if(valid(-a0+x) && valid(-b0+y))   
        {
            if(c == 0) {
                draw[-a0+x][-b0+y][0] += 255;
                draw[-a0+x][-b0+y][1] += 255;
                draw[-a0+x][-b0+y][2] += 255;
            } else if(c == 1) {
                draw[-a0+x][-b0+y][0] += 255;
            } else if(c == 2) {
                draw[-a0+x][-b0+y][1] += 255;
            } else if(c == 3) {
                draw[-a0+x][-b0+y][2] += 255;
            } else if(c == 4) {
                draw[-a0+x][-b0+y][0] += 255;
                draw[-a0+x][-b0+y][1] += 255;
            } else if(c == 5) {
                draw[-a0+x][-b0+y][0] += 255;
                draw[-a0+x][-b0+y][2] += 255;
            } else if(c == 6) {
                draw[-a0+x][-b0+y][1] += 255;
                draw[-a0+x][-b0+y][2] += 255;
            }
        }
        if(valid(a0+y) && valid(b0+x))
        {
            if(c == 0) {
                draw[b0+x][a0+y][0] += 255;
                draw[b0+x][a0+y][1] += 255;
                draw[b0+x][a0+y][2] += 255;
            } else if(c == 1) {
                draw[b0+x][a0+y][0] += 255;
            } else if(c == 2) {
                draw[b0+x][a0+y][1] += 255;
            } else if(c == 3) {
                draw[b0+x][a0+y][2] += 255;
            } else if(c == 4) {
                draw[b0+x][a0+y][0] += 255;
                draw[b0+x][a0+y][1] += 255;
            } else if(c == 5) {
                draw[b0+x][a0+y][0] += 255;
                draw[b0+x][a0+y][2] += 255;
            } else if(c == 6) {
                draw[b0+x][a0+y][1] += 255;
                draw[b0+x][a0+y][2] += 255;
            }
        }
        if(valid(-a0+y) && valid(b0+x))
        {
            if(c == 0) {
                draw[b0+x][-a0+y][0] += 255;
                draw[b0+x][-a0+y][1] += 255;
                draw[b0+x][-a0+y][2] += 255;
            } else if(c == 1) {
                draw[b0+x][-a0+y][0] += 255;
            } else if(c == 2) {
                draw[b0+x][-a0+y][1] += 255;
            } else if(c == 3) {
                draw[b0+x][-a0+y][2] += 255;
            } else if(c == 4) {
                draw[b0+x][-a0+y][0] += 255;
                draw[b0+x][-a0+y][1] += 255;
            } else if(c == 5) {
                draw[b0+x][-a0+y][0] += 255;
                draw[b0+x][-a0+y][2] += 255;
            } else if(c == 6) {
                draw[b0+x][-a0+y][1] += 255;
                draw[b0+x][-a0+y][2] += 255;
            }
        }
        if(valid(a0+y) && valid(-b0+x))
        {
            if(c == 0) {
                draw[-b0+x][a0+y][0] += 255;
                draw[-b0+x][a0+y][1] += 255;
                draw[-b0+x][a0+y][2] += 255;
            } else if(c == 1) {
                draw[-b0+x][a0+y][0] += 255;
            } else if(c == 2) {
                draw[-b0+x][a0+y][1] += 255;
            } else if(c == 3) {
                draw[-b0+x][a0+y][2] += 255;
            } else if(c == 4) {
                draw[-b0+x][a0+y][0] += 255;
                draw[-b0+x][a0+y][1] += 255;
            } else if(c == 5) {
                draw[-b0+x][a0+y][0] += 255;
                draw[-b0+x][a0+y][2] += 255;
            } else if(c == 6) {
                draw[-b0+x][a0+y][1] += 255;
                draw[-b0+x][a0+y][2] += 255;
            }
        }
        if(valid(-a0+y) && valid(-b0+x))
        {
            if(c == 0) {
                draw[-b0+x][-a0+y][0] += 255;
                draw[-b0+x][-a0+y][1] += 255;
                draw[-b0+x][-a0+y][2] += 255;
            } else if(c == 1) {
                draw[-b0+x][-a0+y][0] += 255;
            } else if(c == 2) {
                draw[-b0+x][-a0+y][1] += 255;
            } else if(c == 3) {
                draw[-b0+x][-a0+y][2] += 255;
            } else if(c == 4) {
                draw[-b0+x][-a0+y][0] += 255;
                draw[-b0+x][-a0+y][1] += 255;
            } else if(c == 5) {
                draw[-b0+x][-a0+y][0] += 255;
                draw[-b0+x][-a0+y][2] += 255;
            } else if(c == 6) {
                draw[-b0+x][-a0+y][1] += 255;
                draw[-b0+x][-a0+y][2] += 255;
            }
        }
        b1_new -= (2 * a0) - 3;
    }
}

int main() { 
    srand(time(NULL));
    ppm << "P3\n800 800\n255" << endl;

    double point1[2] = {(rand() / (RAND_MAX + 1.0))*800.0, (rand() / (RAND_MAX + 1.0))*800.0};
    double point2[2] = {(rand() / (RAND_MAX + 1.0))*800.0, (rand() / (RAND_MAX + 1.0))*800.0};
    double point3[2] = {(rand() / (RAND_MAX + 1.0))*800.0, (rand() / (RAND_MAX + 1.0))*800.0};
    // ppm << point1[0] << " " << point1[1] << endl; 
    // ppm << point2[0] << " " << point2[1] << endl; 
    // ppm << point3[0] << " " << point3[1] << endl; 

    draw[(int)point1[0]][(int)point1[1]][0] += 255;
    draw[(int)point2[0]][(int)point2[1]][1] += 255;
    draw[(int)point3[0]][(int)point3[1]][2] += 255;
    draw[(int)point1[0]][(int)point1[1]][0] += 255;
    draw[(int)point2[0]][(int)point2[1]][1] += 255;
    draw[(int)point3[0]][(int)point3[1]][2] += 255;
    draw[(int)point1[0]][(int)point1[1]][0] += 255;
    draw[(int)point2[0]][(int)point2[1]][1] += 255;
    draw[(int)point3[0]][(int)point3[1]][2] += 255;
    
    double dist1 = distance(point1[0], point1[1], point2[0], point2[1]);
    double dist2 = distance(point2[0], point2[1], point3[0], point3[1]);
    double dist3 = distance(point3[0], point3[1], point1[0], point1[1]);
    
    double semi = (dist1 + dist2 + dist3) / 2.0;
    double r_incircle = sqrt(((semi - dist1) * (semi - dist2) * (semi - dist3)) / semi);
    double r_circumcircle = (dist1 * dist2 * dist3) / (4.0 * semi * r_incircle);
    double r_npcircle = r_circumcircle / 2.0;

    double slope1 = (point3[1] - point2[1]) / (point3[0] - point2[0]);
    double slope2 = (point1[1] - point3[1]) / (point1[0] - point3[0]);
    double slope3 = (point2[1] - point1[1]) / (point2[0] - point1[0]);

    double slope1_o = -1.0 * (1.0 / slope1);
    double slope2_o = -1.0 * (1.0 / slope2);
    double slope3_o = -1.0 * (1.0 / slope3);

    double slope1_c = -1.0 * (1.0 / slope1);
    double slope2_c = -1.0 * (1.0 / slope2);
    double slope3_c = -1.0 * (1.0 / slope3);

    double intercept1 = point2[1] - (slope1 * point2[0]);
    double intercept2 = point3[1] - (slope2 * point3[0]);
    double intercept3 = point1[1] - (slope3 * point1[0]);

    double o_intercept1 = point1[1] - (slope1_o * point1[0]);
    double o_intercept2 = point2[1] - (slope2_o * point2[0]);

    double c_intercept1 = ((point3[1] + point2[1]) / 2.0) - (slope1_c * ((point3[0] + point2[0]) / 2.0));
    double c_intercept2 = ((point1[1] + point3[1]) / 2.0) - (slope2_c * ((point1[0] + point3[0]) / 2.0));

    double circum_x = (c_intercept2 - c_intercept1) / (slope1_c - slope2_c);
    double circum_y = slope1_c * circum_x + c_intercept1;
    double ortho_x = (o_intercept2 - o_intercept1) / (slope1_o - slope2_o);
    double ortho_y = slope1_o * ortho_x + o_intercept1;
    double in_x = ((dist2 * point1[0]) + (dist3 * point2[0]) + (dist1 * point3[0])) / (dist1 + dist2 + dist3);
    double in_y = ((dist2 * point1[1]) + (dist3 * point2[1]) + (dist1 * point3[1])) / (dist1 + dist2 + dist3);

    double np_x = (circum_x + ortho_x) / 2.0;
    double np_y = (circum_y + ortho_y) / 2.0;

    double e_slope = (ortho_y - circum_y) / (ortho_x - circum_x);
    double e_intercept = np_y - e_slope * np_x;
    double e_x0 = e_intercept;
    double e_y0 = (-1 * e_intercept) / e_slope;
    double e_x800 = 800.0 * e_slope + e_intercept;
    double e_y800 = (800.0 - e_intercept) / e_slope;
    int e_x1, e_x2, e_y1, e_y2;

    if(valid(e_x0) && valid(e_y0))
    {
        e_x1 = 0.0;
        e_y1 = e_x0;
        e_x2 = e_y0;
        e_y2 = 0;
    }
    else if(valid(e_x0) && valid(e_x800))
    {
        e_x1 = 0.0;
        e_y1 = e_x0;
        e_x2 = 800.0;
        e_y2 = e_x800;
    }
    else if(valid(e_x0) && valid(e_y800))
    {
        e_x1 = 0.0;
        e_y1 = e_x0;
        e_x2 = e_y800;
        e_y2 = 800.0;
    }
    else if(valid(e_y0) && valid(e_x800))
    {
        e_x1 = e_y0;
        e_y1 = 0.0;
        e_x2 = 800.0;
        e_y2 = e_x800;
    }
    else if(valid(e_y0) && valid(e_y800))
    {
        e_x1 = e_y0;
        e_y1 = 0.0;
        e_x2 = e_y800;
        e_y2 = 800.0;
    }
    else if(valid(e_x800) && valid(e_y800))
    {
        e_x1 = 800.0;
        e_y1 = e_x800;
        e_x2 = e_y800;
        e_y2 = 800.0;
    }
    
    // ppm << np_x << " " << np_y << " " << r_npcircle << endl;

    double dx;
    double dy;
    int x0;
    int x1;
    int y0;
    int y1;

    x0 = (int)point1[0];
    x1 = (int)point2[0];
    y0 = (int)point1[1];
    y1 = (int)point2[1];

    drawLine(x0, y0, x1, y1, 0);
    
    // ppm << endl;

    x0 = (int)point2[0];
    x1 = (int)point3[0];
    y0 = (int)point2[1];
    y1 = (int)point3[1];

    drawLine(x0, y0, x1, y1, 0);

    // ppm << endl;

    x0 = (int)point3[0];
    x1 = (int)point1[0];
    y0 = (int)point3[1];
    y1 = (int)point1[1];

    drawLine(x0, y0, x1, y1, 0);

    x0 = (int)e_x1;
    x1 = (int)e_x2;
    y0 = (int)e_y1;
    y1 = (int)e_y2;

    drawLine(x0, y0, x1, y1, 4);

    int a0, b0, b1, b1_new, ty, x, y, r;
    double amax;
    
    x = (int)np_x;
    y = (int)np_y;
    r = (int)r_npcircle;
    
    drawCircle(x, y, r, 1);

    // ppm << endl;

    x = (int)circum_x;
    y = (int)circum_y;
    r = (int)r_circumcircle;
    
    drawCircle(x, y, r, 3);

    // ppm << endl;

    x = (int)in_x;
    y = (int)in_y;
    r = (int)r_incircle;
    
    drawCircle(x, y, r, 2);

    // ppm << endl;

    for(int i = 0; i < 800; i++)
    {
        for(int j = 0; j < 800; j++)
        {
            ppm << draw[i][j][0] << " " << draw[i][j][1] << " " << draw[i][j][2] << "\t";
        }
        ppm << "\n";
    }

    ppm.close();
    return 0;
}