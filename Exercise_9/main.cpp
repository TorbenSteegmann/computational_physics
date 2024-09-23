#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>

#define WITHOUT_NUMPY
#include "../Matplotlib/matplotlibcpp.h"

using namespace std;

// controll configuration usage
#define CONFIGURATION 13

namespace e1
{
// preliminaries
double constexpr pi = M_PI;
double constexpr f_0 = 4;
double constexpr f_1 = 0.25;
double constexpr B_0 = 2 * pi * f_0;
double constexpr h = 2 * pi * f_1;
double constexpr w_0 = B_0;
double constexpr y = 1;
std::vector<double> M_1 = {0, 0, -1};
double constexpr tau = 0.01;
double constexpr t_0 = 0.;


#if CONFIGURATION == 0
std::vector<double> M_0 = {0, 1, 0};
double constexpr t_1 = 1.;
double constexpr T_1 = 0.;
double constexpr T_2 = 0.;
double phi = 0.;
#endif

#if CONFIGURATION == 1
std::vector<double> M_0 = {0, 1, 0};
double constexpr t_1 = 4.;
double constexpr T_1 = 0.;
double constexpr T_2 = 0.;
double phi = 0.;
#endif
#if CONFIGURATION == 2
std::vector<double> M_0 = {0, 1, 0};
double constexpr t_1 = 4.;
double constexpr T_1 = 0.;
double constexpr T_2 = 1.;
double phi = 0.;
#endif
#if CONFIGURATION == 3
std::vector<double> M_0 = {0, 1, 0};
double constexpr t_1 = 4.;
double constexpr T_1 = 1.;
double constexpr T_2 = 0.;
double phi = 0.;
#endif
#if CONFIGURATION == 4
std::vector<double> M_0 = {0, 1, 0};
double constexpr t_1 = 4.;
double constexpr T_1 = 1.;
double constexpr T_2 = 1.;
double phi = 0.;
#endif
#if CONFIGURATION == 5
std::vector<double> M_0 = {0, 1, 0};
double constexpr t_1 = 4.;
double constexpr T_1 = 0.;
double constexpr T_2 = 0.;
double phi = 0.5 * pi;
#endif
#if CONFIGURATION == 6
std::vector<double> M_0 = {0, 1, 0};
double constexpr t_1 = 4.;
double constexpr T_1 = 0.;
double constexpr T_2 = 0.;
double phi = 0.25 * pi;
#endif
#if CONFIGURATION == 7 // example from exercise
std::vector<double> M_0 = {1, 0, 1};
double constexpr t_1 = 4.;
double constexpr T_1 = 0.;
double constexpr T_2 = 1.;
double phi = 0.;
#endif
#if CONFIGURATION == 8
std::vector<double> M_0 = {1, 0, 1};
double constexpr t_1 = 4.;
double constexpr T_1 = 0.;
double constexpr T_2 = 1.;
double phi = 0.5 * pi;
#endif
#if CONFIGURATION == 9
std::vector<double> M_0 = {1, 0, 1};
double constexpr t_1 = 4.;
double constexpr T_1 = 0.;
double constexpr T_2 = 1.;
double phi = 0.25 * pi;
#endif
#if CONFIGURATION == 10
std::vector<double> M_0 = {1, 0, 1};
double constexpr t_1 = 4.;
double constexpr T_1 = 0.;
double constexpr T_2 = 1.;
double phi = 0.5 * pi;
#endif
#if CONFIGURATION == 11
std::vector<double> M_0 = {1, 0, 0};
double constexpr t_1 = 4.;
double constexpr T_1 = 0.;
double constexpr T_2 = 1.;
double phi = 0.25 * pi;
#endif
#if CONFIGURATION == 12
std::vector<double> M_0 = {0, 0, 1};
double constexpr t_1 = 4.;
double constexpr T_1 = 0.;
double constexpr T_2 = 1.;
double phi = 0.25 * pi;
#endif
#if CONFIGURATION == 13
std::vector<double> M_0 = {1, 0, 1};
double constexpr t_1 = 4.;
double constexpr T_1 = 0.;
double constexpr T_2 = 1.;
double phi = 0.5 * pi;
#endif
#if CONFIGURATION == 14
std::vector<double> M_0 = {1, 0, 1};
double constexpr t_1 = 4.;
double constexpr T_1 = 0.;
double constexpr T_2 = 1.;
double phi = 0.5 * pi;
#endif

std::vector<double> B(double t) { return std::vector{h * std::cos(w_0 * (t - 0.5 * tau) + phi), -h * std::sin(w_0 * (t - 0.5 * tau) + phi), B_0}; };

void NMR()
{
    double t = t_0;
    std::vector<double> M_t = M_0;
    std::vector<double> time_x((t_1 / tau) + 1);
    std::vector<double> time_y((t_1 / tau) + 1);
    std::vector<double> time_z((t_1 / tau) + 1);

    time_x[0] = M_0[0];
    time_y[0] = M_0[1];
    time_z[0] = M_0[2];
    // iterate over all time-steps
    for (int i = 1; i < (t_1 / tau) + 1; ++i)
    {
        t += tau;

        // calculation of e^(tau C / 2) * M
        M_t[0] *= std::exp(-tau * 0.5 * T_2);
        M_t[1] *= std::exp(-tau * 0.5 * T_2);
        M_t[2] *= std::exp(-tau * 0.5 * T_1);

        // calc current B
        std::vector<double> B_t = B(t);
        std::vector<std::vector<double>> e_B(3, std::vector<double>(3));
        // calc current omega vals
        double Omega_squared = (B_t[0] * B_t[0]) + (B_t[1] * B_t[1]) + (B_t[2] * B_t[2]);
        double Omega = std::sqrt(Omega_squared);

        // calculation of e^(tau B)
        e_B[0][0] = ((B_t[0] * B_t[0]) + (((B_t[1] * B_t[1]) + (B_t[2] * B_t[2])) * std::cos(Omega * tau * y))) / (Omega_squared);
        e_B[0][1] = ((B_t[0] * B_t[1] * (1 - std::cos(Omega * tau * y))) + (Omega * B_t[2] * std::sin(Omega * tau * y))) / (Omega_squared);
        e_B[0][2] = ((B_t[0] * B_t[2] * (1 - std::cos(Omega * tau * y))) - (Omega * B_t[1] * std::sin(Omega * tau * y))) / (Omega_squared);

        e_B[1][0] = ((B_t[0] * B_t[1] * (1 - std::cos(Omega * tau * y))) - (Omega * B_t[2] * std::sin(Omega * tau * y))) / (Omega_squared);
        e_B[1][1] = ((B_t[1] * B_t[1]) + (((B_t[0] * B_t[0]) + (B_t[2] * B_t[2])) * std::cos(Omega * tau * y))) / (Omega_squared);
        e_B[1][2] = ((B_t[1] * B_t[2] * (1 - std::cos(Omega * tau * y))) + (Omega * B_t[0] * std::sin(Omega * tau * y))) / (Omega_squared);

        e_B[2][0] = ((B_t[0] * B_t[2] * (1 - std::cos(Omega * tau * y))) + (Omega * B_t[1] * std::sin(Omega * tau * y))) / (Omega_squared);
        e_B[2][1] = ((B_t[1] * B_t[2] * (1 - std::cos(Omega * tau * y))) - (Omega * B_t[0] * std::sin(Omega * tau * y))) / (Omega_squared);
        e_B[2][2] = ((B_t[2] * B_t[2]) + (((B_t[1] * B_t[1]) + (B_t[0] * B_t[0])) * std::cos(Omega * tau * y))) / (Omega_squared);

        // calculation of e^(tau B) * M
        auto M_old = M_t;
        M_t[0] = e_B[0][0] * M_old[0] + e_B[0][1] * M_old[1] + e_B[0][2] * M_old[2];
        M_t[1] = e_B[1][0] * M_old[0] + e_B[1][1] * M_old[1] + e_B[1][2] * M_old[2];
        M_t[2] = e_B[2][0] * M_old[0] + e_B[2][1] * M_old[1] + e_B[2][2] * M_old[2];

        // calculation of e^(tau C / 2) * M
        M_t[0] *= std::exp(-tau * 0.5 * T_2);
        M_t[1] *= std::exp(-tau * 0.5 * T_2);
        M_t[2] *= std::exp(-tau * 0.5 * T_1);

        // store for plotting
        time_x[i] = M_t[0];
        time_y[i] = M_t[1];
        time_z[i] = M_t[2];
    }

    // plotting
    namespace plt = matplotlibcpp;

    plt::plot(time_x, {{"label", "$M^x(t)$"}});
    plt::plot(time_y, {{"label", "$M^y(t)$"}});
    plt::plot(time_z, {{"label", "$M^z(t)$"}});

    std::string title = "Nuclear Magnetic Resonance, Configuration: ";
    title = title.append(std::to_string(CONFIGURATION));
    plt::title(title);


#if CONFIGURATION != 0
    std::vector<double> ticks = {0, 50, 100, 150, 200, 250, 300, 350, 400};
    std::vector<std::string> labels = {"0", "0.5", "1", "1.5", "2", "2.5", "3", "3.5", "4"};
    plt::xticks(ticks, labels);
#endif
#if CONFIGURATION == 0
    std::vector<double> ticks = {0, 20, 40, 60, 80, 100};
    std::vector<std::string> labels = {"0", "0.2", "0.4", "0.6", "0.8", "1"};
    plt::xticks(ticks, labels);
#endif

    plt::xlabel("t");
    plt::ylabel("Magnetization");
    plt::grid(true);
    plt::legend();
    plt::show();
}
} // namespace e1

int main()
{
    e1::NMR();
    return 0;
}
