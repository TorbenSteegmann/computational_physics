#define WITHOUT_NUMPY
#include "../Matplotlib/matplotlibcpp.h"

#include <iostream>
#include <vector>

namespace exercise_1 {
// Implement the Euler algorithm
// Harmonic Oscillator (Spring, no friction)
// m = 1, k = 1
// h = 0.1, 0.01, 0.001
// 1 <= j <= 10000 / h
// r(t + h) = r(t) + v(t) * h
// v(t + h) = v(t) + a(t) * h
// a(t) = 1 * F
// F = k * (x_0 - x)
// Discribe why this is not usefull (it explodes)
float r_update(float r, float v, float h)
{
    r += v * h;
    return r;
};
float v_update(float v, float a, float h)
{
    v += a * h;
    return v;
};
float a_update(float r_0, float r, float const K)
{
    return K * (r_0 - r); // M^-1 = 1
}
void spring_simulation(float h)
{
    float r_0 = 0.f;
    float r = r_0;
    float v = 1.f;
    float const K = 1.f;
    float a = a_update(r_0, r, K);
    int iterations = (10000 / h) + 1;
    int test = 10 / h;
    //int iterations = 100;
    std::vector<float> numerical_solution(iterations);
    std::vector<float> analytical_solution(iterations);
    std::vector<float> time(iterations);

    for (int j = 0; j < iterations; ++j) {
        if (j == test)
            std::cout << r << std::endl;
        numerical_solution[j] = r;
        analytical_solution[j] = std::sin(j * h);
        time[j] = j * h;
        r = r_update(r, v, h);
        v = v_update(v, a, h);
        a = a_update(r_0, r, K);
    }

    namespace plt = matplotlibcpp;

    plt::plot(time, numerical_solution, {{"label", "Numerical Solution"}});
    plt::plot(time, analytical_solution, {{"label", "Analytical Solution"}});

    plt::title("Oscillator Simulation using the Euler Algorithm");
    plt::xlabel("Time t");
    plt::ylabel("Position x of Oscillator");
    plt::grid(true);
    plt::legend();
    plt::show();
};
} // namespace exercise_1

namespace exercise_2 {
// Implement the Euler-Cromer algorithm
// Harmonic Oscillator (Spring, no friction)
// m = 1, k = 1
// h = 0.01
// 1 <= j <= 1000

// (a) v(t + h) = v(t) - F(r(t)) * h
// (a) r(t + h) = r(t) + v(t + h) * h

// (b) r(t + h) = r(t) + v(t) * h
// (b) v(t + h) = v(t) - F(r(t+h)) * h

// a(t) = 1 * F
// F = k * (x_0 - x)
float r_update(float r, float v, float h)
{
    r += v * h;
    return r;
};
float v_update(float v, float r, float h)
{
    v -= r * h;
    return v;
};
float a_update(float r_0, float r, float const K)
{
    return K * (r_0 - r); // M^-1 = 1
}
void spring_simulation(float h)
{
    float r_0 = 0.f;
    float r = r_0;
    float v = 1.f;
    int iterations = 10001;
    std::vector<float> numerical_solution(iterations);
    std::vector<float> analytical_solution(iterations);
    std::vector<float> energy(iterations);
    std::vector<float> time(iterations);

    for (int j = 0; j < iterations; ++j) {
        energy[j] = ((v * v) / 2) + (0.5 * r * r);
        numerical_solution[j] = r;
        analytical_solution[j] = std::sin(j * h);
        time[j] = j * h;
        r = r_update(r, v, h);
        v = v_update(v, r, h);
    }

    namespace plt = matplotlibcpp;

    plt::plot(time, energy, {{"label", "Numerical Solution"}});
    plt::plot(time, analytical_solution, {{"label", "Analytical Solution"}});

    plt::title("Oscillator Simulation using the Euler Algorithm");
    plt::xlabel("Time t");
    plt::ylabel("Position x of Oscillator");
    plt::grid(true);
    plt::legend();
    plt::show();
};
} // namespace exercise_2

namespace exercise_3
{
// Implement the velocity Verlet algorithm
// Harmonic Oscillator (Spring, no friction)
// m = 1, k = 1
// h = 0.1, 0.01
// 1 <= j <= 1000

// r(t + h)       = r(t) + v(t) * h + 0.5 * a(t) * h^2
// v(t + h)       = v(t)  + 0.5 * (a(t) + a(t + h)) * h

// a(t) = 1 * F
// F = k * (x_0 - x)
float r_update(float r, float v, float a, float h)
{
    r += v * h + 0.5 * a * h * h;
    return r;
};
float v_update(float v, float a_old, float a_new, float h)
{
    v += 0.5 * (a_old + a_new) * h;
    return v;
};
float a_update(float r_0, float r, float const K)
{
    return K * (r_0 - r); // M^-1 = 1
}
void spring_simulation(float h)
{
    float r_0 = 0.f;
    float r = r_0;
    float v = 1.f;
    float a = 0.f;
    int iterations = 1001;
    float const K = 1.f;
    std::vector<float> numerical_position_solution(iterations);
    std::vector<float> numerical_velocity_solution(iterations);
    std::vector<float> analytical_solution(iterations);
    std::vector<float> energy(iterations);
    std::vector<float> energy_analytical(iterations);
    std::vector<float> time(iterations);

    for (int j = 0; j < iterations; ++j)
    {
        energy[j] = ((v * v) / 2) + (0.5 * r * r);
        energy_analytical[j] = 0.5f;
        numerical_position_solution[j] = r;
        numerical_velocity_solution[j] = v;
        analytical_solution[j] = std::sin(j * h);
        time[j] = j * h;
        r = r_update(r, v, a, h);
        float a_new = a_update(r_0, r, K);
        v = v_update(v, a, a_new, h);
        a = a_new;
    }

    namespace plt = matplotlibcpp;

    plt::plot(time, numerical_position_solution, {{"label", "Position: Numerical Solution"}});
    plt::plot(time, numerical_velocity_solution, {{"label", "Velocity: Numerical Solution"}});
    plt::plot(time, analytical_solution, {{"label", "Position: Analytical Solution"}});
    plt::plot(time, energy, {{"label", "Energy: Numerical Solution"}});
    plt::plot(time, energy_analytical, {{"label", "Energy: Analytical Solution"}});

    plt::title("Oscillator Simulation using the Euler Algorithm");
    plt::xlabel("Time t");
    plt::ylabel("Position x of Oscillator");
    plt::grid(true);
    plt::legend();
    plt::show();
};
} // namespace exercise_3

namespace exercise_4
{
// Implement the velocity Verlet algorithm
// Many Coupled Oscillators (Springs, no friction)
// N = 4, 16, 128
// m = 1, k = 1
// h = 0.1, 0.01
// 1 <= j <= 1000 ?

// r(t + h)       = r(t) + v(t) * h + 0.5 * a(t) * h^2
// v(t + h)       = v(t)  + 0.5 * (a(t) + a(t + h)) * h

// a(t) = 1 * F
// F = k * (x_0 - x)

// Initial Configurations:
// (a) all v = 0
// (a) all x = 0, expect x_(N/2) = 1

// (b) all v = 0
// (b) all x = sin((pi * j * k) / (N + 1)), k = 1,...N, j = 1, N/2
float r_update(float r, float v, float a, float h)
{
    r += v * h + 0.5f * a * h * h;
    return r;
};
float v_update(float v, float a_old, float a_new, float h)
{
    v += 0.5f * (a_old + a_new) * h;
    return v;
};
float a_update(float r_0, float r, float r_2)
{
    return 0 - (2 * r - r_0 - r_2); // M^-1 = 1, K = 1
}
float a_update(float r_0, float r_1)
{
    return 0 - (r_0 - r_1); // M^-1 = 1, K = 1
}
void spring_simulation(float h, int N)
{
    std::vector<float> r(N);
    r[N / 2] = 1;
    std::vector<float> v(N);
    std::vector<float> a(N);
    int iterations = 1001;
    std::vector<std::vector<float>> numerical_position_solution(N, std::vector<float>(iterations));
    std::vector<float> time(iterations);

    for (int j = 0; j < iterations; ++j)
    {
        time[j] = j * h;
        for (int k = 0; k < N; ++k)
        {
            numerical_position_solution[k][j] = r[k];
            r[k] = r_update(r[k], v[k], a[k], h); // change

            float a_new = 0;
            if (k == 0)
                a_new = a_update(r[k], r[k + 1]);
            else if (k == N - 1)
                a_new = a_update(r[k], r[k - 1]);
            else
                a_new = a_update(r[k - 1], r[k], r[k + 1]);


            v[k] = v_update(v[k], a[k], a_new, h);
            a[k] = a_new;
        }
    }

    namespace plt = matplotlibcpp;
    for (int k = 0; k < N; ++k)
    {
        plt::plot(time, numerical_position_solution[k], {{"label", "Position of Element " + std::to_string(k) + ": Numerical Solution"}});
    }

    plt::title("Oscillator Simulation using the Euler Algorithm");
    plt::xlabel("Time t");
    plt::ylabel("Position x of Oscillator");
    plt::grid(true);
    // plt::legend(); // legend turned off for visual clarity
    plt::show();
};
} // namespace exercise_4

int main()
{
    std::cout << "Hello World!" << std::endl;

    // exercise_1::spring_simulation(0.01f);
    // exercise_2::spring_simulation(0.001f);
    // exercise_3::spring_simulation(0.001f);
    exercise_4::spring_simulation(0.01f, 128);

    return 0;
}
