#include <cmath>
#include <iostream>
#include <vector>

#define WITHOUT_NUMPY
#include "../Matplotlib/matplotlibcpp.h"

double const EulerConstant = std::exp(1.0);

using namespace std;

// helper struct
struct results
{
    std::vector<float> U;
    std::vector<float> M;
    std::vector<float> C;
};
namespace task_1
{
// main simulation for the 1D-Lattice
results simulation(double N, double N_Samples)
{
    results r;
    // make random spin config
    std::vector<float> spins(N);
    for (auto& element : spins)
    {
        element = -1 + 2 * (std::rand() % 2);
    }

    // calc E of configuration
    auto E = 0;
    for (int j = 0; j < N - 1; ++j)
        E += spins[j] * spins[j + 1];
    E *= -1;

    // iterate over temperature, double inacc doesnt matter here
    double T = 4.0;
    for (; T > 0; T -= 0.2)
    {
        // recalc beta
        float Beta = 1 / T;
        float E_sum = 0;
        float E_sum_squared = 0;
        // calc for one temperature
        for (int k = 0; k < N_Samples; ++k)
        {
            for (int i = 0; i < N; ++i)
            {
                int j = (int)((static_cast<float>(rand()) / (static_cast<float>(RAND_MAX))) * (N - 1));
                int E_prime = 0;
                int weight_prior = 0;
                int weight_after = 0;
                // special weighting for first and last element
                if (j == 0)
                {
                    weight_prior = spins[j] * spins[j + 1];
                    weight_after = (-1 * spins[j]) * spins[j + 1];
                }
                else if (j == N - 1)
                {
                    weight_prior = spins[j - 1] * spins[j];
                    weight_after = spins[j - 1] * (-1 * spins[j]);
                }
                else
                {
                    weight_prior = spins[j - 1] * spins[j] + spins[j] * spins[j + 1];
                    weight_after = spins[j - 1] * (-1 * spins[j]) + (-1 * spins[j]) * spins[j + 1];
                }
                E_prime = -weight_after + weight_prior;
                auto q = std::pow(EulerConstant, -1 * Beta * E_prime);
                if (q > (static_cast<float>(rand()) / (static_cast<float>(RAND_MAX))))
                {
                    spins[j] *= -1;
                    E = E + E_prime;
                    continue;
                }
            }
            // add E to sum
            E_sum += E;
            E_sum_squared += (E * E);
        }
        // average over all sample iterations
        double U = E_sum / N_Samples;
        double C = (Beta * Beta) * ((E_sum_squared / N_Samples) - (U * U));
        // average over all spins
        U /= N;
        C /= N;
        // add to result to plot later
        r.U.emplace_back(U);
        r.C.emplace_back(C);
    }
    return r;
}


void task_1()
{
    // run simulations for different lattice and sample sizes
    results r1 = simulation(10, 1000);
    results r2 = simulation(10, 10000);
    results r3 = simulation(100, 1000);
    results r4 = simulation(100, 10000);
    results r5 = simulation(1000, 1000);
    results r6 = simulation(1000, 10000);

    // calc theoretical values
    std::vector<float> T_plot;
    std::vector<float> U_theory_plot;
    std::vector<float> C_theory_plot;
    double T = 4.0;
    for (; T > 0; T -= 0.2)
    {
        double beta = 1.0 / T;
        T_plot.emplace_back(T);
        double U_Theory = -((1000. - 1.) / (1000.)) * std::tanh(beta);
        double C_Theory = ((1000. - 1.) / (1000.)) * ((beta / std::cosh(beta)) * (beta / std::cosh(beta)));
        U_theory_plot.emplace_back(U_Theory);
        C_theory_plot.emplace_back(C_Theory);
    }
    // plot it out
    namespace plt = matplotlibcpp;
    plt::plot(T_plot, r1.U, {{"label", "N = " + std::to_string((int)10) + ", $N_{Samples} = $" + std::to_string((int)1000)}});
    plt::plot(T_plot, r2.U, {{"label", "N = " + std::to_string((int)10) + ", $N_{Samples} = $" + std::to_string((int)10000)}});
    plt::plot(T_plot, r3.U, {{"label", "N = " + std::to_string((int)100) + ", $N_{Samples} = $" + std::to_string((int)1000)}});
    plt::plot(T_plot, r4.U, {{"label", "N = " + std::to_string((int)100) + ", $N_{Samples} = $" + std::to_string((int)10000)}});
    plt::plot(T_plot, r5.U, {{"label", "N = " + std::to_string((int)1000) + ", $N_{Samples} = $" + std::to_string((int)1000)}});
    plt::plot(T_plot, r6.U, {{"label", "N = " + std::to_string((int)1000) + ", $N_{Samples} = $" + std::to_string((int)10000)}});
    plt::plot(T_plot, U_theory_plot, {{"label", "U_Theory for N = 1000"}});
    // plot settings
    plt::title("Average Energy per Spin at Temperature T");
    plt::xlabel("Temperature T");
    plt::ylabel("Average Heat per spin $U$ / $N$");
    plt::grid(true);
    plt::legend();
    plt::show();

    plt::plot(T_plot, r1.C, {{"label", "N = " + std::to_string((int)10) + ", $N_{Samples} = $" + std::to_string((int)1000)}});
    plt::plot(T_plot, r2.C, {{"label", "N = " + std::to_string((int)10) + ", $N_{Samples} = $" + std::to_string((int)10000)}});
    plt::plot(T_plot, r3.C, {{"label", "N = " + std::to_string((int)100) + ", $N_{Samples} = $" + std::to_string((int)1000)}});
    plt::plot(T_plot, r4.C, {{"label", "N = " + std::to_string((int)100) + ", $N_{Samples} = $" + std::to_string((int)10000)}});
    plt::plot(T_plot, r5.C, {{"label", "N = " + std::to_string((int)1000) + ", $N_{Samples} = $" + std::to_string((int)1000)}});
    plt::plot(T_plot, r6.C, {{"label", "N = " + std::to_string((int)1000) + ", $N_{Samples} = $" + std::to_string((int)10000)}});
    plt::plot(T_plot, C_theory_plot, {{"label", "C_Theory for N = 1000"}});
    // plot settings
    plt::title("Average Heat per Spin at Temperature T");
    plt::xlabel("Temperature T");
    plt::ylabel("Average Heat per spin $C$ / $N^2$");
    plt::grid(true);
    plt::legend();
    plt::show();
}
} // namespace task_1

namespace task_2
{
// main simulation loop for 2D-Lattice
results simulation(double N, double N_Samples)
{
    results r;
    double N_Relaxation = N_Samples * 0.1;
    double M = 0;
    // make random spin config
    std::vector<std::vector<int>> spins(N, std::vector<int>(N));
    for (auto& row : spins)
    {
        for (auto& element : row)
        {
            element = -1 + 2 * (std::rand() % 2);
            M += element;
        }
    }
    // calc E of configuration
    int E_left = 0;
    for (int i = 0; i < N - 1; ++i)
        for (int j = 0; j < N; ++j)
            E_left += spins[i][j] * spins[i + 1][j];
    int E_right = 0;
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N - 1; ++j)
            E_right += spins[i][j] * spins[i][j + 1];
    int E = -E_left - E_right;

    // Iterate over Temperature
    double T = 4.0;
    for (; T > 0; T -= 0.2)
    {
        // reset values
        float Beta = 1 / T;
        float E_sum = 0;
        float E_sum_squared = 0;
        double M_sum = 0;
        // calc for one temperature
        for (int k = 0; k < N_Samples + N_Relaxation; ++k)
        {
            for (int n = 0; n < N * N; ++n)
            {
                int i = (int)((static_cast<float>(rand()) / (static_cast<float>(RAND_MAX))) * (N - 1));
                int j = (int)((static_cast<float>(rand()) / (static_cast<float>(RAND_MAX))) * (N - 1));

                int weight_prior = 0;
                int weight_after = 0;
                if (i + 1 < N)
                {
                    weight_prior -= spins[i][j] * spins[i + 1][j];
                    weight_after -= -spins[i][j] * spins[i + 1][j];
                }
                if (i - 1 >= 0)
                {
                    weight_prior -= spins[i - 1][j] * spins[i][j];
                    weight_after -= spins[i - 1][j] * -spins[i][j];
                }
                if (j + 1 < N)
                {
                    weight_prior -= spins[i][j] * spins[i][j + 1];
                    weight_after -= -spins[i][j] * spins[i][j + 1];
                }
                if (j - 1 >= 0)
                {
                    weight_prior -= spins[i][j - 1] * spins[i][j];
                    weight_after -= spins[i][j - 1] * -spins[i][j];
                }
                int E_prime = -weight_prior + weight_after;

                auto q = std::pow(EulerConstant, -1 * Beta * E_prime);
                if (q > (static_cast<float>(rand()) / (static_cast<float>(RAND_MAX))))
                {
                    spins[i][j] *= -1;
                    E += E_prime;
                    M += 2 * spins[i][j];
                }
            }
            if (k < N_Relaxation)
                continue;
            E_sum += E;
            E_sum_squared += (E * E);
            M_sum += M;
        }
        // average over all sample size
        auto U = E_sum / N_Samples;
        auto C = (Beta * Beta) * ((E_sum_squared / N_Samples) - (U * U));
        M_sum /= N_Samples;
        // average over spins
        M_sum /= N * N;
        U /= (N * N);
        C /= (N * N);
        r.U.emplace_back(U);
        r.C.emplace_back(C);
        r.M.emplace_back(M_sum);
    }
    return r;
}

void task_2()
{
    // run simulation for different lattice and sample sizes
    results r1 = simulation(10, 1000);
    results r2 = simulation(10, 10000);
    results r3 = simulation(50, 1000);
    results r4 = simulation(50, 10000);
    results r5 = simulation(100, 1000);
    results r6 = simulation(100, 10000);

    // calc theory values
    std::vector<float> T_plot;
    std::vector<float> M_theory_plot;
    double T = 4.f;
    for (; T > 0; T -= 0.2)
    {
        auto beta = 1 / T;
        T_plot.emplace_back(T);
        double M_Theory = 0;
        if (T < ((2) / (std::log1p(std::sqrt(2)))))
        {
            M_Theory = std::pow(1.0 - std::pow(std::sinh(2.0 * beta), -4.0), (1.0 / 8.0));
        }
        M_theory_plot.emplace_back(M_Theory);
    }

    // plot it out
    namespace plt = matplotlibcpp;
    plt::plot(T_plot, r1.U, {{"label", "N = " + std::to_string((int)10) + ", $N_{Samples} = $" + std::to_string((int)1000)}});
    plt::plot(T_plot, r2.U, {{"label", "N = " + std::to_string((int)10) + ", $N_{Samples} = $" + std::to_string((int)10000)}});
    plt::plot(T_plot, r3.U, {{"label", "N = " + std::to_string((int)50) + ", $N_{Samples} = $" + std::to_string((int)1000)}});
    plt::plot(T_plot, r4.U, {{"label", "N = " + std::to_string((int)50) + ", $N_{Samples} = $" + std::to_string((int)10000)}});
    plt::plot(T_plot, r5.U, {{"label", "N = " + std::to_string((int)100) + ", $N_{Samples} = $" + std::to_string((int)1000)}});
    plt::plot(T_plot, r6.U, {{"label", "N = " + std::to_string((int)100) + ", $N_{Samples} = $" + std::to_string((int)10000)}});
    // plot settings
    plt::title("Average Energy per Spin at Temperature T");
    plt::xlabel("Temperature T");
    plt::ylabel("Average Energy per Spin $U$ / $N^2$");
    plt::grid(true);
    plt::legend();
    plt::show();

    plt::plot(T_plot, r1.C, {{"label", "N = " + std::to_string((int)10) + ", $N_{Samples} = $" + std::to_string((int)1000)}});
    plt::plot(T_plot, r2.C, {{"label", "N = " + std::to_string((int)10) + ", $N_{Samples} = $" + std::to_string((int)10000)}});
    plt::plot(T_plot, r3.C, {{"label", "N = " + std::to_string((int)50) + ", $N_{Samples} = $" + std::to_string((int)1000)}});
    plt::plot(T_plot, r4.C, {{"label", "N = " + std::to_string((int)50) + ", $N_{Samples} = $" + std::to_string((int)10000)}});
    plt::plot(T_plot, r5.C, {{"label", "N = " + std::to_string((int)100) + ", $N_{Samples} = $" + std::to_string((int)1000)}});
    plt::plot(T_plot, r6.C, {{"label", "N = " + std::to_string((int)100) + ", $N_{Samples} = $" + std::to_string((int)10000)}});
    // plot settings
    plt::title("Average Heat per Spin at Temperature T");
    plt::xlabel("Temperature T");
    plt::ylabel("Average Heat per spin $C$ / $N^2$");
    plt::grid(true);
    plt::legend();
    plt::show();

    plt::plot(T_plot, r1.M, {{"label", "N = " + std::to_string((int)10) + ", $N_{Samples} = $" + std::to_string((int)1000)}});
    plt::plot(T_plot, r2.M, {{"label", "N = " + std::to_string((int)10) + ", $N_{Samples} = $" + std::to_string((int)10000)}});
    plt::plot(T_plot, r3.M, {{"label", "N = " + std::to_string((int)50) + ", $N_{Samples} = $" + std::to_string((int)1000)}});
    plt::plot(T_plot, r4.M, {{"label", "N = " + std::to_string((int)50) + ", $N_{Samples} = $" + std::to_string((int)10000)}});
    plt::plot(T_plot, r5.M, {{"label", "N = " + std::to_string((int)100) + ", $N_{Samples} = $" + std::to_string((int)1000)}});
    plt::plot(T_plot, r6.M, {{"label", "N = " + std::to_string((int)100) + ", $N_{Samples} = $" + std::to_string((int)10000)}});
    plt::plot(T_plot, M_theory_plot, {{"label", "Analytical result of M for the infinite system"}});
    // plot settings
    plt::title("Average Magnetization per Spin at Temperature T");
    plt::xlabel("Temperature T");
    plt::ylabel("Average Magnetization per spin $M$ / $N^2$");
    plt::grid(true);
    plt::legend();
    plt::show();
}
} // namespace task_2

int main()
{
    std::srand(5378);
    task_1::task_1();
    task_2::task_2();
    return 0;
}
