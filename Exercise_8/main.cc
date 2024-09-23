#define _USE_MATH_DEFINES
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

#define WITHOUT_NUMPY
#include "../Matplotlib/matplotlibcpp.h"

using namespace std::complex_literals;
typedef std::complex<double> complex;

#define BARRIER 1 // Controll whether or not a Barrier is used

namespace e1
{
// Configuration
double const pi = M_PI;
double const sigma = 3.;
double const sigma_squared = 9.;
double const x_0 = 20.;   // starting position
double const q = 1;       // direction
double const Delta = 0.1; // spacial resolution
double const tau = 0.001; // temporal resolution
double const M = 1.;
double const h = 1.;
double const barrier_start = 50.;
double const barrier_end = 50.5;
int const timesteps = 50001;
int const L = 1001;

// V(x) barrier function calculation
double V(float x)
{
#if !BARRIER
    return 0;
#endif
    if (x >= barrier_start && x <= barrier_end)
    {
        return 2.;
    }
    return 0.;
}

// Calculate Starting Value
std::vector<complex> Phi_0()
{
    std::vector<complex> positions(L);
    double x = 0.;
    for (int i = 0; i < L; ++i)
    {
        double x_Delta = x - x_0;
        double const prefix = std::pow(2 * pi * sigma_squared, -0.25) * std::exp(-(x_Delta * x_Delta) / (4 * sigma_squared));

        double a = prefix * std::cos(q * x_Delta);
        double b = prefix * std::sin(q * x_Delta);

        positions[i] = a + 1i * b;

        x += Delta;
    }

    return positions;
}

// Main TDSE function
void TDSE()
{
    complex const c = std::cos(tau / (4 * Delta * Delta));
    complex const is = 1i * std::sin(tau / (4 * Delta * Delta));
    std::vector<complex> positions(L);
    std::vector<std::vector<complex>> time(timesteps, positions);

    time[0] = Phi_0();                  // Set Phi(x, t=0)
    for (int i = 1; i < timesteps; ++i) // Loop over all timesteps
    {
        for (int j = 0; j < L; j += 2) // Loop over all grid points for 1. part of product formula
        {
            positions[j] = time[i - 1][j] * c + is * time[i - 1][j + 1];
            positions[j + 1] = time[i - 1][j] * is + c * time[i - 1][j + 1];
        }
        positions[L - 1] = time[i - 1][0]; // K1

        positions[0] *= 1;
        for (int j = 1; j < L; j += 2) // Loop over all grid points for 2. part of product formula
        {
            complex pj = positions[j] * c + is * positions[j + 1];
            complex pj1 = positions[j] * is + c * positions[j + 1];
            positions[j] = pj;
            positions[j + 1] = pj1;
        } // K2

        double x = 0.;
        for (int j = 0; j < L; ++j) // Loop over all grid points for 3. part of product formula
        {
            complex v_pot = -(1i * tau * (1. + (Delta * Delta) * V(x)) / (Delta * Delta));
            complex e_v = std::exp(v_pot);

            positions[j] *= e_v;
            x += Delta;
        } // V

        positions[0] *= 1;
        for (int j = 1; j < L; j += 2) // Loop over all grid points for 4. part of product formula
        {
            complex pj = positions[j] * c + is * positions[j + 1];
            complex pj1 = positions[j] * is + c * positions[j + 1];
            positions[j] = pj;
            positions[j + 1] = pj1;
        } // K2

        for (int j = 0; j < L; j += 2) // Loop over all grid points for 5. part of product formula
        {
            complex pj = positions[j] * c + is * positions[j + 1];
            complex pj1 = positions[j] * is + c * positions[j + 1];
            positions[j] = pj;
            positions[j + 1] = pj1;
        }
        positions[L - 1] *= 1; // K1

        time[i] = positions; // We only store this for simpler plotting, we could and should not store the entire matrix for real simulations
    }

    // plotting starts here
    namespace plt = matplotlibcpp;

    std::vector<double> probabilities_0(L);
    std::vector<double> probabilities_1(L);
    std::vector<double> probabilities_2(L);
    std::vector<double> probabilities_3(L);
    std::vector<double> probabilities_4(L);
    double post_barrier_0 = 0.;
    double post_barrier_1 = 0.;
    double post_barrier_2 = 0.;
    double post_barrier_3 = 0.;
    double post_barrier_4 = 0.;

    for (int i = 0; i < L; ++i)
    {
        probabilities_0[i] = std::sqrt(std::pow(time[0][i].real(), 2) + std::pow(time[0][i].imag(), 2));
        probabilities_0[i] *= std::sqrt(std::pow(time[0][i].real(), 2) + std::pow(time[0][i].imag(), 2));
        probabilities_0[i] *= Delta;

        probabilities_1[i] = std::sqrt(std::pow(time[(5 / tau)][i].real(), 2) + std::pow(time[(5 / tau)][i].imag(), 2));
        probabilities_1[i] *= std::sqrt(std::pow(time[(5 / tau)][i].real(), 2) + std::pow(time[(5 / tau)][i].imag(), 2));
        probabilities_1[i] *= Delta;

        probabilities_2[i] = std::sqrt(std::pow(time[(40 / tau)][i].real(), 2) + std::pow(time[(40 / tau)][i].imag(), 2));
        probabilities_2[i] *= std::sqrt(std::pow(time[(40 / tau)][i].real(), 2) + std::pow(time[(40 / tau)][i].imag(), 2));
        probabilities_2[i] *= Delta;

        probabilities_3[i] = std::sqrt(std::pow(time[(45 / tau)][i].real(), 2) + std::pow(time[(45 / tau)][i].imag(), 2));
        probabilities_3[i] *= std::sqrt(std::pow(time[(45 / tau)][i].real(), 2) + std::pow(time[(45 / tau)][i].imag(), 2));
        probabilities_3[i] *= Delta;

        probabilities_4[i] = std::sqrt(std::pow(time[(50 / tau)][i].real(), 2) + std::pow(time[(50 / tau)][i].imag(), 2));
        probabilities_4[i] *= std::sqrt(std::pow(time[(50 / tau)][i].real(), 2) + std::pow(time[(50 / tau)][i].imag(), 2));
        probabilities_4[i] *= Delta;

        if (i * Delta < barrier_end)
        {
            continue;
        }
        post_barrier_0 += probabilities_0[i];
        post_barrier_1 += probabilities_1[i];
        post_barrier_2 += probabilities_2[i];
        post_barrier_3 += probabilities_3[i];
        post_barrier_4 += probabilities_4[i];
    }

#if BARRIER
    // normalize probabilities by maximum of all probabilities post barrier
    double post_max = std::max(post_barrier_0, std::max(post_barrier_1, std::max(post_barrier_2, std::max(post_barrier_3, post_barrier_4))));
    for (int i = 0; i < L; ++i)
    {
        if (i * Delta < barrier_end)
            continue;

        probabilities_0[i] /= post_max;
        probabilities_1[i] /= post_max;
        probabilities_2[i] /= post_max;
        probabilities_3[i] /= post_max;
        probabilities_4[i] /= post_max;
    }
#endif

    plt::plot(probabilities_0, {{"label", "t = 0"}});
    plt::plot(probabilities_1, {{"label", "t = 5"}});
    plt::plot(probabilities_2, {{"label", "t = 40"}});
    plt::plot(probabilities_3, {{"label", "t = 45"}});
    plt::plot(probabilities_4, {{"label", "t = 50"}});

    plt::title("TDSE Probability Distribution without Potential Barrier");
#if BARRIER

    plt::title("TDSE Probability Distribution with Potential Barrier");
    plt::axvspan(barrier_start / Delta, barrier_end / Delta);
#endif
    plt::xlabel("x");
    plt::ylabel("P(x,t)");
    std::vector<double> ticks = {0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000};
    std::vector<std::string> labels = {"0", "10", "20", "30", "40", "50", "60", "70", "80", "90", "100"};
    plt::xticks(ticks, labels); // Set custom ticks and labels
    plt::grid(true);
    plt::legend();
    plt::show();
}

} // namespace e1

int main()
{
    e1::TDSE();

    return 0;
}
