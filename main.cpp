#include <iostream>
#include <cstddef>
#include <utility>
#include <stdexcept>
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
        if (this == &o) {
            return *this;
        }
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
    void erase(size_t idx) {
        if (idx >= k) {
            throw std::out_of_range("erase: bad index");
        }
        erase(idx, idx + 1);
    }
    void erase(size_t from, size_t to) {
        if (from > to || to > k) {
            throw std::out_of_range("erase: bad range");
        }
        size_t cnt = to - from;
        if (cnt == 0) {
            return;
        }
        Mtx* nw = nullptr;
        if (cnt < k) {
            nw = new Mtx[k - cnt];
        }
        for (size_t i = 0; i < from; ++i) {
            nw[i] = std::move(mtxs[i]);
        }
        for (size_t i = to; i < k; ++i) {
            nw[i - cnt] = std::move(mtxs[i]);
        }
        delete[] mtxs;
        mtxs = nw;
        k -= cnt;
    }
    void erase(size_t from, size_t to, size_t step) {
        if (step == 0 || from > to || to > k) {
            throw std::out_of_range("erase: bad args");
        }
        size_t cnt = 0;
        for (size_t i = from; i < to; i += step) {
            ++cnt;
        }
        if (cnt == 0) {
            return;
        }
        Mtx* nw = nullptr;
        if (cnt < k) {
            nw = new Mtx[k - cnt];
        }
        size_t j = 0;
        for (size_t i = 0; i < k; ++i) {
            bool del = false;
            if (i >= from && i < to && (i - from) % step == 0) {
                del = true;
            }
            if (!del) {
                nw[j] = std::move(mtxs[i]);
                ++j;
            }
        }
        delete[] mtxs;
        mtxs = nw;
        k -= cnt;
    }
    void insert(size_t pos, const Mtx& m) {
        insert(pos, m, 1);
    }
    void insert(size_t pos, const Mtx& m, size_t count) {
        if (pos > k) {
            throw std::out_of_range("insert: bad pos");
        }
        if (count == 0) {
            return;
        }
        Mtx* nw = new Mtx[k + count];
        for (size_t i = 0; i < count; ++i) {
            nw[pos + i] = m;
        }
        for (size_t i = 0; i < pos; ++i) {
            nw[i] = std::move(mtxs[i]);
        }
        for (size_t i = pos; i < k; ++i) {
            nw[i + count] = std::move(mtxs[i]);
        }
        delete[] mtxs;
        mtxs = nw;
        k += count;
    }
    void insert(size_t pos, const Mtxs& o) {
        if (pos > k) {
            throw std::out_of_range("insert: bad pos");
        }
        if (o.k == 0) {
            return;
        }
        Mtx* nw = new Mtx[k + o.k];
        for (size_t i = 0; i < o.k; ++i) {
            nw[pos + i] = o.mtxs[i];
        }
        for (size_t i = 0; i < pos; ++i) {
            nw[i] = std::move(mtxs[i]);
        }
        for (size_t i = pos; i < k; ++i) {
            nw[i + o.k] = std::move(mtxs[i]);
        }
        delete[] mtxs;
        mtxs = nw;
        k += o.k;
    }
    Mtx operator[](size_t i) {
        return mtxs[i];
    }
    const Mtx& operator[](size_t i) const {
        return mtxs[i];
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
    for (size_t i = 0; i <= s.k; ++i) {
        std::cout << " " << s[i][0][0];
    }
    std::cout << "\n";
}
int main() {
    Mtxs a(4);
    for (size_t i = 0; i < a.k; ++i) {
        a.mtxs[i] = Mtx(2, 2, double(i + 1));
    }
    print(a);
    Mtxs b = a;
    Mtxs c = std::move(b);
    a.erase(1);
    a.erase(0, 1);
    print(a);
    a.insert(1, Mtx(2, 2, 9));
    a.insert(0, Mtx(2, 2, 5), 2);
    a.insert(a.k, c);
    print(a);
    a.erase(0, a.k, 2);
    a[0][1][1] = 42;
    std::cout << "a[0][1][1] = " << a[0][1][1] << "\n";
    print(a);
    std::cout << "conv: " << a.conv(Mtx(2, 2, 1), mul)[0][0] << "\n";
    return 0;
}
