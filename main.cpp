#include <iostream>
#include <cstddef>
#include <utility>
struct Mtx {
    size_t rows;
    size_t cols;
    double* data;
    Mtx() : rows(0), cols(0), data(nullptr) {
    }
    Mtx(size_t r, size_t c, double v = 0.0) : rows(r), cols(c), data(nullptr) {
        if (rows * cols > 0) {
            data = new double[rows * cols];
            for (size_t i = 0; i < rows * cols; ++i) {
                data[i] = v;
            }
        }
    }
    ~Mtx() {
        delete[] data;
    }
    Mtx(const Mtx& o) : rows(o.rows), cols(o.cols), data(nullptr) {
        if (rows * cols > 0) {
            data = new double[rows * cols];
            for (size_t i = 0; i < rows * cols; ++i) {
                data[i] = o.data[i];
            }
        }
    }
    Mtx(Mtx&& o) noexcept : rows(o.rows), cols(o.cols), data(o.data) {
        o.rows = 0;
        o.cols = 0;
        o.data = nullptr;
    }
    Mtx& operator=(const Mtx& o) {
        if (this == &o) {
            return *this;
        }
        delete[] data;
        rows = o.rows;
        cols = o.cols;
        data = nullptr;
        if (rows * cols > 0) {
            data = new double[rows * cols];
            for (size_t i = 0; i < rows * cols; ++i) {
                data[i] = o.data[i];
            }
        }
        return *this;
    }
    Mtx& operator=(Mtx&& o) noexcept {
        if (this != &o) {
            delete[] data;
            rows = o.rows;
            cols = o.cols;
            data = o.data;
            o.rows = 0;
            o.cols = 0;
            o.data = nullptr;
        }
        return *this;
    }
    double* operator[](size_t r) {
        return data + r * cols;
    }
    const double* operator[](size_t r) const {
        return data + r * cols;
    }
};
struct Mtxs {
    size_t k;
    Mtx* mtxs;
    Mtxs() : k(0), mtxs(nullptr) {
    }
    explicit Mtxs(size_t n) : k(n), mtxs(nullptr) {
        if (k > 0) {
            mtxs = new Mtx[k];
        }
    }
    ~Mtxs() {
        delete[] mtxs;
    }
    Mtxs(const Mtxs& o) : k(o.k), mtxs(nullptr) {
        if (k > 0) {
            mtxs = new Mtx[k];
            for (size_t i = 0; i < k; ++i) {
                mtxs[i] = o.mtxs[i];
            }
        }
    }
    Mtxs(Mtxs&& o) noexcept : k(o.k), mtxs(o.mtxs) {
        o.k = 0;
        o.mtxs = nullptr;
    }
    Mtxs& operator=(const Mtxs& o) {
        delete[] mtxs;
        k = o.k;
        mtxs = nullptr;
        if (k > 0) {
            mtxs = new Mtx[k];
            for (size_t i = 0; i < k; ++i) {
                mtxs[i] = o.mtxs[i];
            }
        }
        return *this;
    }
    Mtxs& operator=(Mtxs&& o) noexcept {
        if (this != &o) {
            delete[] mtxs;
            k = o.k;
            mtxs = o.mtxs;
            o.k = 0;
            o.mtxs = nullptr;
        }
        return *this;
    }
    Mtx conv(Mtx m, Mtx(*t)(Mtx, Mtx)) {
        for (size_t i = 0; i < k; ++i) {
            m = t(m, mtxs[i]);
        }
        return m;
    }
};
Mtx mul(Mtx a, Mtx b) {
    Mtx r(a.rows, b.cols);
    for (size_t i = 0; i < a.rows; ++i) {
        for (size_t j = 0; j < b.cols; ++j) {
            for (size_t x = 0; x < a.cols; ++x) {
                r[i][j] += a[i][x] * b[x][j];
            }
        }
    }
    return r;
}
void print(const Mtxs& s) {
    std::cout << "k = " << s.k << ":";
    for (size_t i = 0; i < s.k; ++i) {
        std::cout << " " << s.mtxs[i][0][0];
    }
    std::cout << "\n";
}
int main() {
    Mtxs a(3);
    for (size_t i = 0; i < a.k; ++i) {
        a.mtxs[i] = Mtx(2, 2, double(i + 1));
    }
    print(a);
    Mtxs b = a;
    Mtxs c = std::move(b);
    print(c);
    std::cout << "conv: " << a.conv(Mtx(2, 2, 1), mul)[0][0] << "\n";
    return 0;
}
