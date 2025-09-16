#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>
#include <iomanip>
#include <map>
#include <vector>
#include <algorithm>

using namespace std;

const int PRECISION = 10;
const int K = 7; // Number of points to fit polynomial of degree 6

// Convert a string in given base to decimal
long double convertToDecimal(const string& value, int base) {
    long double result = 0;
    for (char c : value) {
        int digit = isdigit(c) ? c - '0' : tolower(c) - 'a' + 10;
        result = result * base + digit;
    }
    return result;
}

// Compute log base b of value
long double logBase(long double value, int base) {
    return log(value) / log(base);
}

// Solve Ax = b using Gaussian elimination
vector<long double> solvePolynomial(const vector<long double>& x, const vector<long double>& y) {
    int n = x.size();
    vector<vector<long double>> A(n, vector<long double>(n));
    vector<long double> b = y;

    // Build Vandermonde matrix
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            A[i][j] = pow(x[i], j);

    // Gaussian elimination
    for (int i = 0; i < n; ++i) {
        // Pivot
        for (int k = i + 1; k < n; ++k) {
            if (abs(A[k][i]) > abs(A[i][i])) {
                swap(A[i], A[k]);
                swap(b[i], b[k]);
            }
        }

        // Eliminate
        for (int k = i + 1; k < n; ++k) {
            long double factor = A[k][i] / A[i][i];
            for (int j = i; j < n; ++j)
                A[k][j] -= factor * A[i][j];
            b[k] -= factor * b[i];
        }
    }

    // Back-substitution
    vector<long double> coeffs(n);
    for (int i = n - 1; i >= 0; --i) {
        coeffs[i] = b[i];
        for (int j = i + 1; j < n; ++j)
            coeffs[i] -= A[i][j] * coeffs[j];
        coeffs[i] /= A[i][i];
    }

    return coeffs;
}

int main() {
    ifstream file("input.json");
    if (!file.is_open()) {
        cerr << "Error: Could not open input.json" << endl;
        return 1;
    }

    string line, key;
    int base = 0;
    string value;
    vector<long double> y_values;

    cout << fixed << setprecision(PRECISION);

    while (getline(file, line)) {
        line.erase(remove(line.begin(), line.end(), ' '), line.end());

        if (line.find("\"base\"") != string::npos) {
            size_t start = line.find("\"", line.find(":") + 1) + 1;
            size_t end = line.find("\"", start);
            base = stoi(line.substr(start, end - start));
        }

        if (line.find("\"value\"") != string::npos) {
            size_t start = line.find("\"", line.find(":") + 1) + 1;
            size_t end = line.find("\"", start);
            value = line.substr(start, end - start);

            long double decimalValue = convertToDecimal(value, base);
            long double logResult = logBase(decimalValue, base);
            y_values.push_back(logResult);
            cout << "log base " << base << " of " << value << " = " << logResult << endl;
        }
    }

    // Use first K points to fit polynomial
    if (y_values.size() < K) {
        cerr << "Not enough points to fit polynomial of degree " << K - 1 << endl;
        return 1;
    }

    vector<long double> x_values;
    for (int i = 1; i <= K; ++i)
        x_values.push_back(i);

    vector<long double> y_subset(y_values.begin(), y_values.begin() + K);
    vector<long double> coeffs = solvePolynomial(x_values, y_subset);

    cout << "\nPolynomial Coefficients (lowest degree to highest):" << endl;
    for (int i = 0; i < coeffs.size(); ++i)
        cout << "a[" << i << "] = " << coeffs[i] << endl;

    return 0;
}
