#define _USE_MATH_DEFINES
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

#define WITHOUT_NUMPY
#include "../Matplotlib/matplotlibcpp.h"

namespace exam
{
using namespace std::complex_literals;
namespace plt = matplotlibcpp;
typedef std::complex<double> complex;

double const Delta = 0.025;
int const L = 1201;
double const tau = 0.00025;
int const m = 40000;
double const pi = M_PI;

// controlls configurations
double const Omega = 1.;
double const sigma = 1.;
double const x_0 = 0.;

double const sigma_squared = sigma * sigma;

std::vector<complex> Phi_0() // Initial Wave Function
{
    std::vector<complex> phi_0(L);
    for (int i = 0; i < L; ++i)
    {
        complex x = (i - L / 2) * Delta;
        complex term_1 = std::pow(pi * sigma_squared, -0.25);
        complex term_2 = std::exp(-(((x - x_0) * (x - x_0)) / (2. * sigma_squared)));
        phi_0[i] = term_1 * term_2;
    };

    return phi_0;
}

std::vector<complex> V() // (Omega^2 / 2) * x^2
{
    std::vector<complex> v(L);
    for (int i = 0; i < L; ++i)
    {
        auto x = (i - L / 2) * Delta;
        v[i] = ((Omega * Omega) / 2) * (x * x);
    }
    return v;
}

void TDSE() // driving function, implementation of the product forumla approach
{
    std::vector<std::vector<complex>> phi_t; // matrix to store phi at different time steps for plotting
    std::vector<complex> phi = Phi_0();      // initialize phi as described in exercise
    phi_t.emplace_back(phi);                 // add phi(t=0) to plot later
    std::vector<double> probabilities(L);    // for plotting
    std::vector<double> xv1(m + 1);          // <x>
    double xv2 = 0.;                         // <x^2>
    std::vector<double> variance_n(m + 1);   // variance vector numerical
    std::vector<double> expectation_a(m + 1); // variance vector analytical
    std::vector<double> variance_a(m + 1);    // variance vector analytical
    for (int i = 0; i < L; ++i)              // variance at t=0
    {
        double x = (i - L / 2) * Delta;
        double probability;
        probability = std::sqrt(std::pow(phi[i].real(), 2) + std::pow(phi[i].imag(), 2));
        probability *= std::sqrt(std::pow(phi[i].real(), 2) + std::pow(phi[i].imag(), 2));
        probability *= Delta;
        xv1[0] += probability * x;
        xv2 += (probability * (x * x));
    }
    variance_n[0] = (xv2 - (xv1[0] * xv1[0]));
    xv2 = 0;
    expectation_a[0] = (x_0)*std::cos(0);
    variance_a[0] = (0.5 * sigma_squared) * (std::cos(Omega * 0) * std::cos(Omega * 0))
                    + ((std::sin(Omega * 0) * std::sin(Omega * 0)) / (2 * sigma_squared * Omega * Omega));

    // pre-computations
    complex const c = std::cos(tau / (4. * Delta * Delta));
    complex const is = 1i * std::sin(tau / (4. * Delta * Delta));
    auto v = V();
    // Product formula implementation
    for (int j = 1; j < m + 1; ++j)
    {
        double time = j * tau;
        complex p1;
        complex p2;
        for (int i = 0; i + 1 < L; i += 2) // K1/2
        {
            p1 = c * phi[i] + is * phi[i + 1];
            p2 = is * phi[i] + c * phi[i + 1];
            phi[i] = p1;
            phi[i + 1] = p2;
        }
        for (int i = 1; i < L; i += 2) // K2/2
        {
            p1 = c * phi[i] + is * phi[i + 1];
            p2 = is * phi[i] + c * phi[i + 1];
            phi[i] = p1;
            phi[i + 1] = p2;
        }
        for (int i = 0; i < L; ++i) // V
        {
            phi[i] *= std::exp(-(1i * tau * (1. + (Delta * Delta) * v[i] / (Delta * Delta))));
        }
        for (int i = 1; i < L; i += 2) // K2/2
        {
            p1 = c * phi[i] + is * phi[i + 1];
            p2 = is * phi[i] + c * phi[i + 1];
            phi[i] = p1;
            phi[i + 1] = p2;
        }
        for (int i = 0; i + 1 < L; i += 2) // K1/2
        {
            p1 = c * phi[i] + is * phi[i + 1];
            p2 = is * phi[i] + c * phi[i + 1];
            phi[i] = p1;
            phi[i + 1] = p2;
        }

        // Variance
        for (int i = 0; i < L; ++i)
        {
            double x = (i - L / 2) * Delta;
            double probability;
            probability = std::sqrt(std::pow(phi[i].real(), 2) + std::pow(phi[i].imag(), 2));
            probability *= std::sqrt(std::pow(phi[i].real(), 2) + std::pow(phi[i].imag(), 2));
            probability *= Delta;
            xv1[j] += probability * x;
            xv2 += probability * (x * x);
        }
        variance_n[j] = (xv2 - (xv1[j] * xv1[j]));
        xv2 = 0;
        expectation_a[j] = (x_0)*std::cos(Omega * time);
        variance_a[j] = (0.5 * sigma_squared) * (std::cos(Omega * time) * std::cos(Omega * time))
                        + ((std::sin(Omega * time) * std::sin(Omega * time)) / (2 * sigma_squared * Omega * Omega));

        if (time != 2 && time != 4 && time != 6 && time != 8 && time != 10)
            continue;

        phi_t.emplace_back(phi);
    }

    // Probability calculations
    for (int i = 0; i < phi_t.size(); ++i)
    {
        for (int j = 0; j < L; ++j)
        {
            probabilities[j] = std::sqrt(std::pow(phi_t[i][j].real(), 2) + std::pow(phi_t[i][j].imag(), 2));
            probabilities[j] *= std::sqrt(std::pow(phi_t[i][j].real(), 2) + std::pow(phi_t[i][j].imag(), 2));
            probabilities[j] *= Delta;
        }
        plt::plot(probabilities, {{"label", "t = " + std::to_string(2 * i)}});
    }

    // Plotting
    plt::title("Quantum Harmonic Oscillator Position \nk=" + std::to_string((int)Omega) + ", s=" + std::to_string((int)sigma)
               + ", $x_0$=" + std::to_string((int)x_0));
    plt::xlabel("x");
    plt::ylabel("P(x,t)");
    std::vector<double> xticks = {0, 200, 400, 600, 800, 1000, 1200};
    std::vector<std::string> labels = {"-15", "-10", "-5", "0", "5", "10", "15"};
    plt::xticks(xticks, labels); // Set custom ticks and labels
    plt::grid(true);
    plt::legend();
    plt::show();

    plt::title("Quantum Harmonic Oscillator Variance \nk=" + std::to_string((int)Omega) + ", s=" + std::to_string((int)sigma)
               + ", $x_0$=" + std::to_string((int)x_0));
    plt::xlabel("t");
    plt::ylabel("variance");
    plt::plot(variance_n, {{"label", "$<x^2> - <x>^2$ numerical"}});
    plt::plot(xv1, {{"label", "<x> numerical"}});
    plt::plot(expectation_a, {{"label", "<x> analytical"}});
    plt::plot(variance_a, {{"label", "$<x^2> - <x>^2$ analytical"}});
    xticks = {0, 8000, 16000, 24000, 32000, 40000};
    labels = {"0", "2", "4", "6", "8", "10"};
    plt::xticks(xticks, labels);
    plt::grid(true);
    plt::legend();
    plt::show();
};
} // namespace exam


int main()
{
    exam::TDSE();
    return 0;
}
