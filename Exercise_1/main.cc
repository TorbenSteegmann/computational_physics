#include <iostream>
#include <vector>

// plotting library available here: https://github.com/lava/matplotlib-cpp
#define WITHOUT_NUMPY // setting for matplotlib
#include "../Matplotlib/matplotlibcpp.h"

namespace task_1
{
struct Index
{
    int row = 0;
    int column = 0;
};

struct Matrix
{
    // Constructor for a row_length x column_length matrix with random number in its cells
    Matrix(int row_length, int column_length) : row_length(row_length), column_length(column_length)
    {
        // 3: create vectors, initialized at -5 because that is the smallest number
        largest_row = std::vector<float>(column_length, -5);
        largest_column = std::vector<float>(row_length, -5);
        for (int i = 0; i < row_length; ++i)
        {
            matrix.emplace_back();
            for (int j = 0; j < column_length; ++j)
            {
                // 1.:create random element between -5.0 and 5.0
                float element = -5.f + static_cast<float>(rand()) / (static_cast<float>(RAND_MAX / 11));
                // place element into matrix
                matrix[i].emplace_back(element);
                std::cout << matrix[i][j] << ',';
                // 3.: check if element is largest found in column j so far
                if (element > largest_row[j])
                {
                    // 3.: if so place replace current element
                    largest_row[j] = element;
                }
                // 3.: check if element smaller or equal to the currently largest element in row i
                if (element <= largest_column[i])
                {
                    // if so there is no need for the last check as we know it cant be the largest in the matrix
                    continue;
                }
                // 3.: else replace currently largest row element
                largest_column[i] = element;
                // 2.: check if element is largest in matrix
                if (element > matrix[largest_element_index.row][largest_element_index.column])
                {
                    // 2.: if so replace the index of the largest element
                    largest_element_index.row = i;
                    largest_element_index.column = j;
                }
            }
            std::cout << ";\n";
        }

        // only do the next part if it is mathematically possible (for this exercise it always does)
        if (row_length != column_length)
            return;

        for (int i = 0; i < row_length; ++i)
        {
            // 3.: multiply row with column vector
            largest_scalar += largest_row[i] * largest_column[i];
        }
    };

    Matrix(Matrix A, Matrix B) // constructor for 4.
    {
        // make sure it is mathematically possible
        if (A.column_length != B.row_length)
        {
            std::cout << "Invalid Matrix sizes.\n";
            return;
        }

        for (int k = 0; k < A.row_length; ++k)
        {
            matrix.emplace_back();
            for (int i = 0; i < A.column_length; ++i)
            {
                matrix[k].emplace_back();
                for (int j = 0; j < B.row_length; ++j)
                {
                    // standard matrix multiplication
                    matrix[k][i] += A.matrix[k][j] * B.matrix[j][i];
                }
                std::cout << matrix[k][i] << ',';
            }
            std::cout << ";\n";
        }
    };

    // output function for 2. (value never changes in exercise -> no need to recalculate)
    void largest_element()
    {
        std::cout << "The largest element has a value of " << matrix[largest_element_index.row][largest_element_index.column] << " and is in row "
                  << largest_element_index.row << " and column " << largest_element_index.column << ".\n";
    };

    // output for 3. (value never changes in exercise -> no need to recalculate)
    void largest_scalar_result() { std::cout << "The result of the vector multiplication is " << largest_scalar << ".\n"; };

    std::vector<std::vector<float>> matrix;
    int largest_scalar = 0;
    int row_length;
    int column_length;
    Index largest_element_index;
    std::vector<float> largest_row;
    std::vector<float> largest_column;
};

void task_1()
{
    std::srand(5378);

    // 1.: create matrix A
    std::cout << "Matrix A: \n";
    Matrix A(6, 6);
    // 2.: output largest element
    A.largest_element();
    // 3.: output largest scalar
    A.largest_scalar_result();

    // 4.: Create Matrix B
    std::cout << "\nMatrix B: \n";
    Matrix B(6, 6);

    // 4.: Create Matrix C as A * B
    std::cout << "\nMatrix C: \n";
    Matrix C(A, B);

    // 4.: Create Matrix D as B * A
    std::cout << "\nMatrix D: \n";
    Matrix D(B, A);
}
} // namespace task_1

namespace task_2
{
// 1.: cheby function
auto cheby(std::vector<float> x, int N)
{
    std::vector<std::vector<float>> result(N + 1); // N + 1 rows

    std::vector<float> chebn(N + 1);
    // solves the values of 1 datapoint element for the N+1 chebychev Polynomials recursively
    auto chebx = [&chebn](float element, int n) mutable -> void
    {
        // level of indirection for recursion
        auto chebx_impl = [&chebn](std::vector<float>& p_chebn, float element, int n, auto& chebx_ref) mutable -> void
        {
            if (n == 0)
            {
                chebn[n] = 1;
                return;
            }
            if (n == 1)
            {
                chebx_ref(chebn, element, n - 1, chebx_ref);
                chebn[n] = element;
                return;
            }

            chebx_ref(chebn, element, n - 1, chebx_ref);
            chebn[n] = 2 * element * chebn[n - 1] - chebn[n - 2];
            return;
        };
        chebx_impl(chebn, element, n, chebx_impl);
        return;
    };

    std::cout << "\nChebyshev Polynomial Matrix:\n";

    for (int i = 0; i < x.size(); ++i)
    {
        // chebx gets called for every element
        chebx(x[i], N);
        for (int j = 0; j <= N; ++j)
        {
            // place value of polynomial j in row j
            result[j].emplace_back(chebn[j]);
        }
    }

    for (auto const& row : result)
    {
        for (auto const& element : row)
        {
            std::cout << element << '\t';
        }
        std::cout << '\n';
    }

    return result;
}

void task_2()
{
    namespace plt = matplotlibcpp; // convenience

    // create a vector x with 201 data points
    std::vector<float> x(201);
    float e = -1.f;
    for (int i = 0; i < 201; ++i) {
        x[i] = e;
        e += 0.01f;
    }
    // 1.: create matrix
    auto cheb = task_2::cheby(x, 5);

    // 2.: Plot Polynomials
    for (int i = 0; i < cheb.size(); ++i)
    {
        std::vector<float> y = cheb[i];

        plt::plot(x, y, {{"label", "N = " + std::to_string(i)}});
    }

    // plot settings
    plt::title("Chebyshev Polynomials");
    plt::xlabel("x");
    plt::ylabel("$T_N$(x)");
    plt::grid(true);
    plt::legend();
    plt::show();
}
} // namespace task_2

int main()
{
    task_1::task_1();

    task_2::task_2();

    return 0;
}
