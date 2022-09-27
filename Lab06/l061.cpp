// Varun Asuri

#include <algorithm>
#include <iostream> 
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <math.h>
#include <queue>
#include <cmath>
#include <ctime>
#include <list>
#include <map>
#include <set>

using namespace std;

int width, height;
const double PI = 3.141592653589793238462643383279502884;
vector<int> shade;
// ofstream debug("debug.txt");

int xGradient(int n)
{
    int total = 0;
    if(n-1-width >= 0)
        total += shade[n-1-width];
    if(n-1 >= 0)
        total += 2*shade[n-1];
    if(n+width-1 < width*height)
        total += shade[n+width-1];
    if(n-width+1 >= 0)
        total -= shade[n-width+1];
    if(n+1 < width*height)
        total -= 2*shade[n+1];
    if(n+width+1 < width*height)
        total -= shade[n+width+1];
    return total;
}
int yGradient(int n)
{
    int total = 0;
    if(n-1-width >= 0)
        total += shade[n-1-width];
    if(n-width >= 0)
        total += 2*shade[n-width];
    if(n-width+1 >= 0)
        total += shade[n-width+1];
    if(n+width-1 < width*height)
        total -= shade[n+width-1];
    if(n+width < width*height)
        total -= 2*shade[n+width];
    if(n+width+1 < width*height)
        total -= shade[n+width+1];
    return total;
}
// void drawLine(int x0, int y0, int x1, int y1)
// {
//     double dx;
//     double dy;
//   
//     dx = x1 - x0;
//     dy = y1 - y0;
//   
//     if (abs(dx) >= abs(dy))
//     {
//         double e = 0;
//         double r = dy/dx;
//         while(x0 != x1)
//         {
//             draw[x0][y0][0] += 255;
//             draw[x0][y0][1] += 255;
//             draw[x0][y0][2] += 255;
//             if(e > 0)
//             {
//                 if(dy < 0)
//                     y0 -= 1;
//                 else
//                     y0 += 1;
//                 e -= 1;
//             }
//             e += abs(r);
//             if(dx < 0)
//                 x0 -= 1;
//             else
//                 x0 += 1;
//         }
//     }
//     else
//     {
//         double e = 0;
//         double r = dx/dy;
//         while(y0 != y1)
//         {
//             draw[x0][y0][0] += 255;
//             draw[x0][y0][1] += 255;
//             draw[x0][y0][2] += 255;
//             if(e > 0)
//             {
//                 if(dx < 0)
//                     x0 -= 1;
//                 else
//                     x0 += 1;
//                 e -= 1;
//             }
//             e += abs(r);
//             if(dy < 0)
//                 y0 -= 1;
//             else
//                 y0 += 1;
//         }
//     }
// }

void part1() {
    srand(time(NULL));
    ifstream read("image.ppm");
    ofstream ppm1("imagef.ppm");
    ofstream ppm2("imageCC.ppm");

    string fill;
    read >> fill; read >> fill;
    width = stoi(fill);
    read >> fill;
    height = stoi(fill);
    read >> fill;

    ppm1 << "P3\n" << width << " " << height << "\n255" << endl;
    ppm2 << "P3\n" << width << " " << height << "\n255" << endl;
    
    int r, g, b, a;
    while(!read.eof()) {
        read >> fill;
        r = stoi(fill);
        read >> fill;
        g = stoi(fill);
        read >> fill;
        b = stoi(fill);
        a = (r + g + b)/3;
        shade.push_back(a);
    }

    vector<int> dthresh;
    vector<int> nmsgrid;
    vector<int> nmsvals;
    vector<int> canedge;
    vector<int> votesct;
    vector<double> angles;

    double sum, angle;
    int xg = 0, yg = 0;
    for(int w = 0; w < width*height; w++) {
        xg = pow(xGradient(w), 2);
        yg = pow(yGradient(w), 2);
        sum = sqrt(xg+yg);
        angle = atan( ((double)yg) / ((double)xg) ) * 180/PI;
        nmsgrid.push_back(sum);
        angles.push_back(angle);

        if(sum >= 300)
            sum = 255;
        else if(sum < 50)
            sum = 0;
        else
            sum = 127;

        dthresh.push_back(sum);
    }

    double replace;
    int round, n1, n2, output;
    map<int, pair<int, int>> angloc;
    angloc[0] = make_pair(1, -1);
    angloc[45] = make_pair(1-width, width-1);
    angloc[90] = make_pair(width, -width);
    angloc[135] = make_pair(-1-width, width+1);
    for(int q = 0; q < width*height; q++) {
        replace = fmod(angles[q]+180, 180);
        round = ((int)((replace+22.5)/45))*45;
        pair<int, int> dir = angloc[round];

        if(q+dir.first >= 0 && q+dir.first < width*height)
            n1 = nmsgrid[q+dir.first];
        else
            n1 = -1;

        if(q+dir.second >= 0 && q+dir.second < width*height)
            n2 = nmsgrid[q+dir.second];
        else
            n2 = -1;
        
        if(nmsgrid[q] > n1 && nmsgrid[q] > n2)
            output = 255;
        else
            output = 0;
        
        nmsvals.push_back(output);
    }

    int curr = 0;
    bool condition;
    queue<int> itfill;
    vector<int> check;
    for(int z = 0; z < width*height; z++) {
        if(dthresh[z] != 127)
            continue;
        itfill.push(z);
        condition = false;
        check.push_back(z);
        while(!itfill.empty()) {
            curr = itfill.front();
            itfill.pop();
            if(dthresh[curr] == 0)
                continue;
            if(dthresh[curr] == 255) {
                condition = true;
                continue;
            }
            if(curr-1 >= 0 && count(check.begin(), check.end(), curr-1) == 0) {
                itfill.push(curr-1);
                check.push_back(curr-1);
            } if(curr-width >= 0 && count(check.begin(), check.end(), curr-width) == 0) {
                itfill.push(curr-width);
                check.push_back(curr-width);
            } if(curr+1 < width*height && count(check.begin(), check.end(), curr+1) == 0) {
                itfill.push(curr+1);
                check.push_back(curr+1);
            } if(curr+width < width*height && count(check.begin(), check.end(), curr+width) == 0) {
                itfill.push(curr+width);
                check.push_back(curr+width);
            } if(curr-1-width >= 0 && count(check.begin(), check.end(), curr-1-width) == 0) {
                itfill.push(curr-1-width);
                check.push_back(curr-1-width);
            } if(curr+1-width >= 0 && count(check.begin(), check.end(), curr+1-width) == 0) {
                itfill.push(curr+1-width);
                check.push_back(curr+1-width);
            } if(curr-1+width >= 0 && count(check.begin(), check.end(), curr-1+width) == 0) {
                itfill.push(curr-1+width);
                check.push_back(curr-1+width);
            } if(curr+1+width >= 0 && count(check.begin(), check.end(), curr+1+width) == 0) {
                itfill.push(curr+1+width);
                check.push_back(curr+1+width);
            }
        }
        
        if(condition) {
            for(int b = 0; b < check.size(); b++)
                dthresh[check[b]] = 255;
        } else {
            for(int b = 0; b < check.size(); b++)
                dthresh[check[b]] = 0;
        }
        check.clear();
    }

    int temp = 0;
    for(int p = 0; p < width*height; p++) {
        if(dthresh[p] == 255 && nmsvals[p] == 255)
            temp = 255;
        else
            temp = 0;
        canedge.push_back(temp);
        votesct.push_back(0);
        ppm1 << temp << " " << temp << " " << temp << "\t";
        if((p+1)%width==0)
            ppm1 << "\n";
    }

    for(int w = 0; w < width*height; w++) {
        if(canedge[w] == 255) {
            xg = pow(xGradient(w), 2);
            yg = pow(yGradient(w), 2);
            angle = atan( ((double)yg) / ((double)xg) );
            
            double dx;
            double dy;
            int x0, y0;
        
            dx = cos(angle);
            dy = sin(angle);
            x0 = w % width;
            y0 = w / width;

            // debug << dx << " " << dy << endl;

            if(abs(dx) >= abs(dy))
            {
                double e = 0;
                double r = dy/dx;
                while(x0 != 0 && x0 != width && y0 != 0 && y0 != height)
                {
                    votesct[y0*width + x0] += 1;
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
                while(x0 != 0 && x0 != width && y0 != 0 && y0 != height)
                {
                    votesct[(y0*width)+x0] += 1;
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

            dx = -cos(angle);
            dy = -sin(angle);
            x0 = w % width;
            y0 = w / width;
        
            if (abs(dx) >= abs(dy))
            {
                double e = 0;
                double r = dy/dx;
                while(x0 != 0 && x0 != width && y0 != 0 && y0 != height)
                {
                    votesct[(y0*width)+x0] += 1;
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
                while(x0 != 0 && x0 != width && y0 != 0 && y0 != height)
                {
                    votesct[y0*width + x0] += 1;
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
    }

    for(int z = 0; z < width*height; z++) {
        if(votesct[z] > 150)
            ppm2 << "255 0 0\t";
        else
            ppm2 << votesct[z] << " " << votesct[z] << " " << votesct[z] << "\t";
        if((z+1)%width==0)
            ppm2 << "\n";
    }

    read.close();
    ppm1.close();
    ppm2.close();
    // debug.close();
}

int main() {
    part1();
    return 0;
}
