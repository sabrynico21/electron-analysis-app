// =============================================================================
// _mb_core.cpp — compiled MB coordinate-descent core.
//
// Near-literal translation of `huge_core.cpp` (huge 2.0.1):
//   * mb_partial_residual  (huge_core.cpp:640-649)
//   * mb_refine_active     (huge_core.cpp:651-677)
//   * mb                   (huge_core.cpp:679-810)
//
// Only the solver is compiled. Everything else (CLR, correlation matrix, lambda
// path, StARS, selection, refit) stays in Python, unchanged.
//
// Fidelity notes — the C++ keeps exactly the operations of the original:
//   * column-major accessor `cm`, same as huge (the input array is requested in
//     Fortran order, so the same expression indexes the same element);
//   * scalar accumulation loops for the partial residuals (no BLAS);
//   * coordinate sweep order j = 0..d-1 skipping the response m;
//   * sequential strong rule 2*lambda[i] - lambda[i-1] (2*lambda[0] for i = 0);
//   * active-set refinement with relative coefficient change `thol`;
//   * pruning pass, then the full KKT certification pass that re-runs the
//     strong loop on violations;
//   * same soft-threshold (`threshold_l1`) and the same `thol` / MAX_ITER.
//
// Deliberate differences from the original:
//   * the OpenMP `#pragma omp for` over the response index m is not used (the
//     caller parallelises independent subsamples instead), which does not change
//     any per-m result;
//   * `collect_sorted` is not reproduced: the Python side consumes a dense
//     (n_lambda, d, d) array, and the values are identical regardless of the
//     emission order;
//   * `hit_max_iter` is returned to the caller instead of being dropped.
//
// Build flags (see setup.py): -O2 -ffp-contract=off, no -ffast-math and no
// -march=native, so the floating-point behaviour matches R's build (-O2).
// =============================================================================

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace py = pybind11;

namespace {

// huge_core.cpp:16 — soft threshold.
inline double threshold_l1(double x, double thr) {
    if (x > thr) return x - thr;
    if (x < -thr) return x + thr;
    return 0.0;
}

// huge_core.cpp:23 — column-major accessor for a raw const pointer.
inline double cm(const double* data, int nrows, int r, int c) {
    return data[static_cast<std::size_t>(c) * static_cast<std::size_t>(nrows) +
                static_cast<std::size_t>(r)];
}

// huge_core.cpp:640 — residual gradient of coordinate j for the m-th nodewise
// lasso: r = S[m,j] - sum_k S[j, a_k] * w[a_k] over the active set.
inline double mb_partial_residual(const double* S_data, int d, int m, int j,
                                  const int* idx_a, int size_a, const double* w) {
    double r = cm(S_data, d, m, j);
    for (int k = 0; k < size_a; k++) r -= cm(S_data, d, j, idx_a[k]) * w[idx_a[k]];
    return r;
}

// huge_core.cpp:651 — cycle coordinate descent over the active set until the
// relative coefficient change drops below thol.
void mb_refine_active(const double* S_data, int d, int m, double ilambda,
                      double thol, int max_iter, const int* idx_a, int size_a,
                      double* w0, double* w1) {
    double gap_int = 1;
    int iter_int = 0;
    while (gap_int > thol && iter_int < max_iter) {
        double change_sum = 0, coef_sum = 0;
        for (int j = 0; j < size_a; j++) {
            int w_idx = idx_a[j];
            double r = cm(S_data, d, m, w_idx) + w0[w_idx];
            for (int k = 0; k < size_a; k++)
                r -= cm(S_data, d, w_idx, idx_a[k]) * w0[idx_a[k]];

            w1[w_idx] = threshold_l1(r, ilambda);
            coef_sum += std::fabs(w1[w_idx]);
            change_sum += std::fabs(w1[w_idx] - w0[w_idx]);
            w0[w_idx] = w1[w_idx];
        }
        gap_int = (coef_sum > 0) ? change_sum / coef_sum : 0;
        iter_int++;
    }
}

// huge_core.cpp:679 — MB graph estimation over the whole lambda path,
// warm-started from one lambda to the next.
py::tuple solve_mb_path(
    py::array_t<double, py::array::f_style | py::array::forcecast> S_in,
    py::array_t<double, py::array::c_style | py::array::forcecast> lambda_in,
    double thol, int max_iter) {
    if (S_in.ndim() != 2) throw std::invalid_argument("S must be a 2-D array.");
    const py::ssize_t d_ss = S_in.shape(0);
    if (S_in.shape(1) != d_ss) throw std::invalid_argument("S must be square.");
    if (d_ss < 2)
        throw std::invalid_argument("At least two features are required.");
    if (lambda_in.ndim() != 1 || lambda_in.shape(0) < 1)
        throw std::invalid_argument("lambda_path must be a non-empty 1-D array.");
    if (max_iter < 1) throw std::invalid_argument("max_iter must be >= 1.");

    const int d = static_cast<int>(d_ss);
    const int nlambda = static_cast<int>(lambda_in.shape(0));
    const double* S_data = S_in.data();
    const double* lambda = lambda_in.data();

    py::array_t<double> out(
        {static_cast<py::ssize_t>(nlambda), d_ss, d_ss});
    double* out_data = out.mutable_data();
    std::fill(out_data, out_data + static_cast<std::size_t>(nlambda) *
                                     static_cast<std::size_t>(d) * d,
              0.0);

    // Per-response scratch, reused across m exactly as in the original.
    std::vector<double> w0(d, 0.0), w1(d, 0.0);
    std::vector<int> idx_a(d), idx_i(d);
    std::vector<unsigned char> strong(d, 0);
    std::vector<double> grad_abs(d, 0.0);

    bool hit_max_iter = false;

    // The numeric core touches no Python object: release the GIL.
    {
        py::gil_scoped_release release;

        for (int m = 0; m < d; m++) {
            bool hit_any_max_iter = false;
            idx_i[m] = 0;
            for (int j = 0; j < m; j++) idx_i[j] = 1;
            for (int j = m + 1; j < d; j++) idx_i[j] = 1;

            int size_a = 0;
            std::fill(w0.begin(), w0.end(), 0.0);
            std::fill(w1.begin(), w1.end(), 0.0);

            // Gradient at w = 0 is the m-th row of S.
            for (int j = 0; j < d; j++)
                grad_abs[j] = std::fabs(cm(S_data, d, m, j));

            for (int i = 0; i < nlambda; i++) {
                const double ilambda = lambda[i];
                const double strong_thr = (i > 0)
                    ? 2.0 * lambda[i] - lambda[i - 1]
                    : 2.0 * lambda[i];
                for (int j = 0; j < d; j++)
                    strong[j] =
                        (idx_i[j] == 1 && grad_abs[j] > strong_thr) ? 1 : 0;

                int iter_ext = 0;
                bool certified = false;
                while (!certified && iter_ext < max_iter) {
                    int gap_ext = 1;
                    while (gap_ext != 0 && iter_ext < max_iter) {
                        const int size_a_prev = size_a;
                        for (int j = 0; j < d; j++) {
                            if (idx_i[j] == 1 && strong[j]) {
                                double r = mb_partial_residual(
                                    S_data, d, m, j, idx_a.data(), size_a,
                                    w0.data());
                                grad_abs[j] = std::fabs(r);
                                w1[j] = threshold_l1(r, ilambda);
                                if (w1[j] != 0) {
                                    idx_a[size_a++] = j;
                                    idx_i[j] = 0;
                                }
                                w0[j] = w1[j];
                            }
                        }
                        gap_ext = size_a - size_a_prev;

                        mb_refine_active(S_data, d, m, ilambda, thol, max_iter,
                                         idx_a.data(), size_a, w0.data(),
                                         w1.data());
                        int junk_a = 0;
                        for (int j = 0; j < size_a; j++) {
                            int w_idx = idx_a[j];
                            if (w1[w_idx] == 0) {
                                junk_a++;
                                idx_i[w_idx] = 1;
                                strong[w_idx] = 1;  // recently active: keep sweeping
                            } else {
                                idx_a[j - junk_a] = w_idx;
                            }
                        }
                        size_a -= junk_a;
                        iter_ext++;
                    }

                    // Certification: one full pass over every inactive coordinate.
                    int violations = 0;
                    for (int j = 0; j < d; j++) {
                        if (idx_i[j] == 1) {
                            double r = mb_partial_residual(
                                S_data, d, m, j, idx_a.data(), size_a, w0.data());
                            grad_abs[j] = std::fabs(r);
                            w1[j] = threshold_l1(r, ilambda);
                            if (w1[j] != 0) {
                                idx_a[size_a++] = j;
                                idx_i[j] = 0;
                                strong[j] = 1;
                                violations++;
                            }
                            w0[j] = w1[j];
                        }
                    }
                    certified = (violations == 0);
                }
                if (!certified) hit_any_max_iter = true;

                double* row =
                    out_data + (static_cast<std::size_t>(i) * d + m) * d;
                for (int k = 0; k < size_a; k++) row[idx_a[k]] = w1[idx_a[k]];
            }

            if (hit_any_max_iter) hit_max_iter = true;
        }
    }

    return py::make_tuple(out, hit_max_iter);
}

}  // namespace

PYBIND11_MODULE(_mb_core, m) {
    m.doc() =
        "Compiled MB coordinate-descent core (translation of huge_core.cpp::mb).";
    m.def("solve_mb_path", &solve_mb_path, py::arg("S"), py::arg("lambda_path"),
          py::arg("thol") = 1e-4, py::arg("max_iter") = 10000,
          "Warm-started MB nodewise Lasso over the whole lambda path.\n\n"
          "S must be the d x d correlation matrix (Fortran order is used\n"
          "internally, as in huge). Returns (coefficients, hit_max_iter) where\n"
          "coefficients has shape (n_lambda, d, d) with [i, node, predictor].");
}
