#include <iostream>
#include <vector>

// plotting library available here: https://github.com/lava/matplotlib-cpp
#define WITHOUT_NUMPY // setting for matplotlib
#include "../Matplotlib/matplotlibcpp.h"

namespace task_2 {
void task_2()
{
    namespace plt = matplotlibcpp; // convenience

    std::srand(3227); // seed
    int N_part = 10000;
    int N = 1000;
    std::vector<int> lattice(2 * N + 1); // 0 initialized lattice
    lattice[N] = N_part;

    std::vector<float> variance_per_N;
    for (int j = 0; j < N + 1; ++j)
    {
        auto lattice_copy = lattice;
        for (int i = 0; i < lattice.size(); ++i)
        {
            int it = lattice[i];
            for (int k = 0; k < it; ++k)
            {
                int direction = -1 + 2 * (std::rand() % 2); // choose direction
                --lattice_copy[i];                          // move particle in lattice
                ++lattice_copy[i + direction];
            }
        }
        lattice = lattice_copy;
        // square up distances
        int total_distance_squared = 0;
        int total_distance = 0;
        for (int i = 0; i < lattice.size(); ++i)
        {
            total_distance_squared += (lattice[i] * ((i - N) * (i - N)));
            total_distance += (lattice[i] * (i - N));
        }
        float mean_squared = ((float)total_distance_squared) / ((float)N_part);
        float squared_mean = ((float)total_distance) / ((float)N_part);
        squared_mean *= squared_mean;
        float variance = mean_squared - squared_mean;
        variance_per_N.emplace_back(variance);
    }

    std::vector<int> x_axis(2 * N + 1);
    std::vector<float> y_axis = variance_per_N;

    x_axis[0] = -N;
    for (int i = 1; i < x_axis.size(); ++i)
    {
        x_axis[i] = x_axis[i - 1] + 1;
    }

    // 2.: Plot Polynomials
    plt::plot(x_axis, lattice);

    // plot settings
    plt::title("Random Walk");
    plt::xlabel("Final position x");
    plt::ylabel("Number of particels");
    plt::grid(true);
    plt::show();
}
} // namespace task_2

int main()
{
    task_2::task_2();

    return 0;
}
