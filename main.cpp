#include <iostream>
#include <fstream>
#include <cmath>
#include <map>
#include <vector>
#include <algorithm>

using namespace std;

struct Pair {
    double x;
    double value;
};

vector<Pair> U_table;
vector<Pair> T_table;
map<string, double> Gtext_table;

// Зчитування таблиці U(x) або T(x)
bool readTable(const string& filename, vector<Pair>& table) {
    ifstream fin(filename);
    if (!fin.is_open()) return false;
    double x, val;
    while (fin >> x >> val) {
        table.push_back({x, val});
    }
    return true;
}

// Зчитування текстової таблиці
bool readGtext(const string& filename) {
    ifstream fin(filename);
    if (!fin.is_open()) return false;
    string text;
    double val;
    while (fin >> text >> val) {
        Gtext_table[text] = val;
    }
    return true;
}

// Інтерполяція
double interpolate(const vector<Pair>& table, double x) {
    for (size_t i = 0; i < table.size() - 1; ++i) {
        if (table[i].x <= x && x <= table[i + 1].x) {
            double x0 = table[i].x, y0 = table[i].value;
            double x1 = table[i + 1].x, y1 = table[i + 1].value;
            return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
        }
    }
    return table.back().value;
}

// Алгоритм 2: U1(x), T1(x), Qqn1
double U1(double x) {
    return atan(asin(sin(3 * x)));
}

double T1(double x) {
    return atan(acos(sin(2 * x)));
}

double Qqn1(double x, double y, double z) {
    return x / U1(x) + y * T1(y) - U1(z) * T1(z);
}

double Qnk1(double x, double y) {
    return 1.1 * Qqn1(x, y, x + y) - 0.9 * Qqn1(y, x, x - y);
}

double Rnk_alg2(double x, double y) {
    double q1 = Qnk1(x, y), q2 = Qnk1(y, x);
    return x * q1 + y * q2 - 0.03 * q1 * q2;
}

// Qqn для алгоритму 1
double Qqn(double x, double y, double z, bool& uFail, bool& tFail) {
    if (U_table.empty()) uFail = true;
    if (T_table.empty()) tFail = true;
    if (uFail || tFail) return 0;
    double u = interpolate(U_table, x);
    double t = interpolate(T_table, y);
    double uz = interpolate(U_table, z);
    double tz = interpolate(T_table, z);
    return x / u + y * t - uz * tz;
}

double Qnk(double x, double y, bool& uFail, bool& tFail) {
    return Qqn(x, y, x + y, uFail, tFail) - Qqn(y, x, x - y, uFail, tFail);
}

double Rnk(double x, double y, bool& uFail, bool& tFail) {
    if (uFail || tFail || x <= 5 || T_table.empty() || U_table.empty()) return Rnk_alg2(x, y);
    return x * Qnk(x, y, uFail, tFail) + y * Qnk(y, x, uFail, tFail);
}

// Алгоритм 3
double Qnk2(double x, double y) {
    return 1.3 * Qqn1(x, y, x) - 0.7 * Qqn1(y, x, x);
}

double func_alg3(double x, double y, double z) {
    double q1 = Qnk2(x, y), q2 = Qnk2(y, x);
    return 1.75 * x * q1 + 1.25 * y * q2 - 1.5 * q1 * q2;
}

// Функція Gtext
double Gtext(const string& text) {
    auto it = Gtext_table.find(text);
    return (it != Gtext_table.end()) ? it->second : 0.0;
}

double CText(double x, const string& text) {
    if (text.empty()) return Gtext("set");
    if (text == "set") return 0;
    return Gtext(text);
}

double Max4(double a, double b, double c, double d) {
    return max({a, b, c, d});
}

// Y(x)
double Y(double x, bool& valid) {
    double v = x * x * 100 - x * x;
    if (v < 1 || x <= 0) {
        valid = false;
        return 0;
    }
    return log(v);
}

// Yrr
double Yrr(double f, double r, bool& valid) {
    return Y(f, valid) * r + 0.5 * Y(r, valid);
}

// Trr
double Trr(double f, double r, bool& valid) {
    if ((f * f - r) < 4) {
        valid = false;
        return 0;
    }
    return (f * f - r + Yrr(r, f, valid)) / 2.0;
}

// RText = Rrr(f, r)
double RText(double x, double y, double z, const string& text) {
    double f = CText(Max4(x, y, x + z, y + z), text);
    bool valid = true;

    // r використовується для побудови 2k, але зберігаємо k окремо:
    bool tempUFail = false, tempTFail = false;
    double r_val = Rnk(x, y, tempUFail, tempTFail);
    double k_val = r_val;

    double trr1 = Trr(f, r_val, valid);
    double trr2 = Trr(f, 2 * k_val, valid); // тут ВИПРАВЛЕНО: 2 * k

    if (!valid) return 0;
    return f * trr1 + r_val * trr2;
}

// Основна функція Variant
double Variant(double x, double y, double z, const string& text) {
    bool uFail = !readTable("dat1.dat", U_table);
    bool tFail = !readTable("dat2.dat", T_table);

    if (!readGtext("dat3.dat")) {
        cerr << "Помилка: неможливо відкрити файл dat3.dat" << endl;
        exit(1);
    }

    double r = (!uFail && !tFail)
        ? Rnk(x, y, uFail, tFail) + Rnk(y, z, uFail, tFail) * Rnk(x, y, uFail, tFail)
        : func_alg3(x, y, z);

    double k = RText(x, y, z, text);

    return 0.8973 * r + 0.1027 * k;
}

// main
int main() {
    double x, y, z;
    string text;
    cout << "Введіть x, y, z: ";
    cin >> x >> y >> z;
    cout << "Введіть текст: ";
    cin >> text;

    double result = Variant(x, y, z, text);
    cout << "Результат Variant(x, y, z, text) = " << result << endl;
    return 0;
}
