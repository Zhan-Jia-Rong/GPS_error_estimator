#include <tuple>
#include "stdio.h"
#include "math.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include "vector"
#include "iomanip"

#define M_PI 3.14159265358979323846

using namespace std;

std::vector<double> vinc(double latp, double latc, double longp, double longc)
{
    constexpr double req = 6378137.0;             //Radius at equator
    constexpr double flat = 1 / 298.257223563;    //flattening of earth
    constexpr double rpol = (1 - flat) * req;

    double sin_sigma, cos_sigma, sigma, sin_alpha, cos_sq_alpha, cos2sigma;
    double C, lam_pre;

    // convert to radians
    latp = M_PI * latp / 180.0;
    latc = M_PI * latc / 180.0;
    longp = M_PI * longp / 180.0;
    longc = M_PI * longc / 180.0;

    const double u1 = atan((1 - flat) * tan(latc));
    const double u2 = atan((1 - flat) * tan(latp));

    double lon = longp - longc;
    double lam = lon;
    double tol = pow(10.0,-12.0); // iteration tolerance
    double diff = 1.;

    while (abs(diff) > tol)
    {
        sin_sigma = sqrt(pow((cos(u2) * sin(lam)), 2.) + pow(cos(u1) * sin(u2) - sin(u1) * cos(u2) * cos(lam), 2.));
        cos_sigma = sin(u1) * sin(u2) + cos(u1) * cos(u2) * cos(lam);
        sigma = atan(sin_sigma / cos_sigma);
        sin_alpha = (cos(u1) * cos(u2) * sin(lam)) / sin_sigma;
        cos_sq_alpha = 1 - pow(sin_alpha, 2.);
        cos2sigma = cos_sigma - ((2 * sin(u1) * sin(u2)) / cos_sq_alpha);
        C = (flat / 16) * cos_sq_alpha * (4 + flat * (4 - 3 * cos_sq_alpha));
        lam_pre = lam;
        lam = lon + (1 - C) * flat * sin_alpha * (sigma + C * sin_sigma * (cos2sigma + C * cos_sigma * (2 * pow(cos2sigma, 2.) - 1)));
        diff = abs(lam_pre - lam);
    }
    const double usq = cos_sq_alpha * ((pow(req, 2.) - pow(rpol, 2.)) / pow(rpol, 2.));
    const double A = 1 + (usq / 16384) * (4096 + usq * (-768 + usq * (320 - 175 * usq)));
    const double B = (usq / 1024) * (256 + usq * (-128 + usq * (74 - 47 * usq)));
    const double delta_sig = B * sin_sigma * (cos2sigma + 0.25 * B * (cos_sigma * (-1 + 2 * pow(cos2sigma, 2.)) -
        (1 / 6) * B * cos2sigma * (-3 + 4 * pow(sin_sigma, 2.)) *
        (-3 + 4 * pow(cos2sigma, 2.))));
    const double dis = rpol * A * (sigma - delta_sig);
    const double azi1 = atan2((cos(u2) * sin(lam)), (cos(u1) * sin(u2) - sin(u1) * cos(u2) * cos(lam)));
    std::vector<double> v(2);
    v[0] = dis;
    v[1] = azi1;
    std::cout << "dis=" << v[0] << "m  deg=" << v[1] * 180 / M_PI << std::endl;
    return v;
}
void read_GPS()
{
    ifstream sr("gps_test.csv", ios::in);
    if (!sr)
    {
        cout << "開啟檔案失敗！" << endl;
        exit(1);
    }
    string line;
    vector<string> ours_data(2);
    char del=',';
    int b = 0;
    while (getline(sr, line)) //read csv 
    {
        int a = line.find(",");
        if (a == line.length() - 1 || a + 1 >= line.length())
        {
            cout << line <<" <-Error input" << endl;
            break;
        }
        if (b == 1)
        {
            double od = stold(ours_data[0]);
            double od1 = stold(line.substr(0, a));
            cout << ours_data[0] << " " << ours_data[1] << endl;
            cout << od << " " << od1 << endl;

            //vinc(od, stod(ours_data[1]), od1 , stod(line.substr(a + 1, line.length())));
            //cout << stod(ours_data[0]) << " " << ours_data[1] << "\n" << line.substr(0, a) << " " << line.substr(a + 1, line.length()) << "\n\n" << endl;
            b = 0;
        }
        ours_data[0] = line.substr(0,a);
        ours_data[1] = line.substr(a+1, line.length());
        b = 1;
        
        
    }
}
int main()
{
    double lat1, lat2, long1, long2;
    lat1 = 24.79573446;
    lat2 = 24.79576845;
    long1 = 121.0367003;
    long2 = 121.0366977;
    
    string str1="24.7950928944444";
    string str2="24.7954168571428";
    string str3="121.036531572222";
    string str4="121.03661195238";

    double ans1=79573446;
    double ans2=79576845;
    double ans3=367003;
    double ans4=366977;
    ans1=ans1/100000000;
    ans2=ans2/100000000;
    ans3=ans3/10000000;
    ans4=ans4/10000000;

    vinc(stod(str1),stod(str2),stod(str3),stod(str4));
    //vinc(lat1, lat2, long1, long2);
    //read_GPS();
    return 0;
}