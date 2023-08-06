#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <iostream>
#include <windows.h>
#include <cmath>
#include <fstream>
#include <string>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}


using namespace std;
const string path = "C:/Games/Out"; //Путь к выводу файла
const double f = 0.0001;
const double rb = 0.0001;
const double dt = 0.01; //разница времен ???
const int nx = 385; //сетка по х
const int ny = 420; //сетка по y
const double g = 9.8; //скорость свободного падения
const int dx = 2; // разница горизонтальная координата
const int tn = 100000;
const int kx = 1000;
const int ro = 1000; //плотность воды
const double pi = 3.1415926535897932384626433832795; //число Пи
QString nolab;
QString n;
QString no;
typedef double tt1[nx][ny]; //двумерный массив 385 на 420

double ww;
double hs, k0, rs, mt, bb, c, l, k, T0, w00, w0, h01, hh, w1, r; //k -волновое число
int nx1, i, j, t, j0_0; //j0_0 - переменная j0, t - Время
string ts;
int Xg[ny];
tt1 Lm0, Lm, Ks, Cs, h, h0, u, v, fx, fy; //сет
POINT gr1[nx], gr2[nx], gr3[nx];

double Zaph()
{
    cout << "\nИсходная матрица: \n";
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            h[i][j] = (rand() % 2101) - 100;
            cout << h[i][j] << " ";
        }

        cout << "\n";
    }
    cout << "\n";
    return 0;
}

double ZapCs()
{
    cout << "\nИсходная матрица: \n";
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            Cs[i][j] = (rand() % 2101) - 100;
            cout << Cs[i][j] << " ";
        }

        cout << "\n";
    }
    cout << "\n";
    return 0;
}

double ZapKs()
{
    cout << "\nИсходная матрица: \n";
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            Ks[i][j] = (rand() % 2101) - 100;
            cout << Ks[i][j] << " ";
        }

        cout << "\n";
    }
    cout << "\n";
    return 0;
}

double Zapfx()
{
    cout << "\nИсходная матрица: \n";
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            fx[i][j] = (rand() % 2101) - 100;
            cout << fx[i][j] << " ";
        }

        cout << "\n";
    }
    cout << "\n";
    return 0;
}

double Zapfy()
{
    cout << "\nИсходная матрица: \n";
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            fy[i][j] = (rand() % 2101) - 100;
            cout << fy[i][j] << " ";
        }

        cout << "\n";
    }
    cout << "\n";
    return 0;
}

double Zaph0()
{
    cout << "\nИсходная матрица: \n";
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            h0[i][j] = (rand() % 2101) - 100;
            cout << h0[i][j] << " ";
        }

        cout << "\n";
    }
    cout << "\n";
    return 0;
}

double Zapv()
{
    cout << "\nИсходная матрица: \n";
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            h0[i][j] = (rand() % 2101) - 100;
            cout << h0[i][j] << " ";
        }

        cout << "\n";
    }
    cout << "\n";
    return 0;
}

double Zapu()
{
    cout << "\nИсходная матрица: \n";
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            u[i][j] = (rand() % 2101) - 100;
            cout << u[i][j] << " ";
        }

        cout << "\n";
    }
    cout << "\n";
    return 0;
}

double ZapLm()
{
    cout << "\nИсходная матрица: \n";
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            Lm[i][j] = (rand() % 2101) - 100;
            cout << Lm[i][j] << " ";
        }
        cout << "\n";
    }
    cout << "\n";
    return 0;
}
double ZapLm0()
{
    cout << "\nИсходная матрица: \n";
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            Lm0[i][j] = (rand() % 2101) - 100;
            cout << Lm0[i][j] << " ";
        }
        cout << "\n";
    }
    cout << "\n";
    return 0;
}

double th(double a1)
{
    return (exp(a1) - exp(-a1)) / (exp(a1) + exp(-a1));
}

// функция для нахождения KK
double KW(double tt, double h0, double& KK)
{
    double KK1, a, b;
    KK = 1;
    b = KK * h0;
    KK1 = 4 * pi * pi / tt / tt / g / th(b); //TT - возможно период
    a = fabs((KK - KK1) / KK);
    if (a > 0.01)
    {
        KK = KK1;
    }
    nolab = QString::number(a,'g',15);
    return a;
}

double plan(tt1 zz) {
    int  i, j;
    double d, hmax, hmin;
    hmax = 0;
    hmin = 0;
    for (i = 2; i <= (nx - 1); i++) {
        for (j = 2; j <= (ny - 1); j++) {
            if (zz[i][j] > hmax) {
                hmax = zz[i][j];
            }
            if (zz[i][j] < hmin) {
                hmin = zz[i][j];
            }
        }
    }
    d = round(10 * (hmax - hmin));
    //cout << "Минимальная высота" << hmin << endl;
        //cout << "Максимальная высота" << hmax << endl;
    //cout << "Разница в высоте" << d << endl;
    n = QString::number(d,'g',15);
        return d;
}
double MAX(int t) {
    int i, j;
    double tt;

    tt = t % 10000;

    if ((abs(tt) < 0.001) || (t == 0)) {
        for (i = 2; i <= nx - 1; i++) {
            for (j = 2; j <= ny - 1; j++) {
                Lm0[i][j] = Lm[i][j];
            }
        }
        for (i = 2; i <= nx - 1; i++) {
            for (j = 2; j <= ny - 1; j++) {
                Lm[i][j] = 0;
            }
        }
    }

    for (i = 2; i <= nx - 1; i++) {
        for (j = 2; j <= ny - 1; j++) {
            if (Lm[i][j] < abs(h0[i][j])) {
                Lm[i][j] = abs(h0[i][j]);
            }
        }
    }
    c = 2 * pi / T0 / k;
    string c1 = to_string(c);
    r = rb / hh;
    for (j = 1; j <= ny; j++) {
        if (j < 280) {
            Xg[j] = 105 + (j - 1);
        }
        else {
            Xg[j] = 385;
        }
    }
    for (i = 1; i < nx; i++) {
        for (j = 1; j < ny; j++) {
            h0[i][j] = 0;
            Lm[i][j] = 0;
        }
    }
    for (i = 1; i <= nx; i++) {
        for (j = 1; j <= ny; j++) {
            if (i <= Xg[j]) {
                h[i][j] = hh * (1 - pow((j - 210), 2) / 210 / 210) * (1 - i * i / Xg[j] / Xg[j]) + 1;
            }
            else {
                h[i][j] = 0;
            }
            if (j == 1 || j == ny || i == nx) {
                h[i][j] = 0;
            }
            if ((nx - i) + j < 245) {
                h[i][j] = 0;
            }
            u[i][j] = 0; //u{0}[i][j] = 0
            v[i][j] = 0; //v{ 0 }[i][j]
        }
    }
    for (i = 1; i < nx; i++) {
        for (j = 1; j < ny; j++) {
            hs = h[i][j];
            if (hs > 0) {
                KW(T0, hs, k); // ввод переменных в функцию кв
                Ks[i][j] = k;
                Cs[i][j] = 2 * pi / T0 / Ks[i][j];
            }
            else {
                Cs[i][j] = 0;
            }
        }
        for (i = 2; i < (nx - 1); i++) {
            for (j = 2; j < (ny - 1); j++) {
                fx[i][j] = -g * (h0[i][j] - h0[i - 1][j]) / dx;
                fy[i][j] = -g * (h0[i][j] - h0[i][j - 1]) / dx;
            }
        }
        for (i = 2; i < (nx - 1); i++) {
            for (j = 2; j < (ny - 1); j++) {
                u[i][j] = u[i][j] * exp(-r * dt) + fx[i][j] / r * (1 - exp(-r * dt));
                v[i][j] = v[i][j] * exp(-r * dt) + fy[i][j] / r * (1 - exp(-r * dt));
            }
        }

        for (j = 1; j < ny; j++) {
            u[1][j] = -g / c * h0[1][j];
            u[nx][j] = 0;
        }

        for (i = 1; i < nx; i++) {
            v[i][1] = 0;
            v[i][ny] = 0;
        }

        for (j = j0_0; j < (ny - 1); j++) {
            u[20][j] = 0;
            u[21][j] = 0;
        }
        cout << "h" << h[i][j] << endl;
        cout << "Ks" << Ks[i][j] << endl;
        for (i = 1; i < (nx - 1); i++) {
            for (j = 1; j < (ny - 1); j++) {
                cout << "bb1231 = " << bb << endl;
                cout << "h = " << h[i][j] << endl;
                cout << "Ks" << Ks[i][j] << endl;
                if (h[i][j] > 0) {
                    w0 = 0;
                    bb = Ks[i][j] * h[i][j];
                    if (i == 3) {
                        w0 = w00 * sin(2 * pi / T0 * (t - 1) * dt);
                    }
                    k0 = Ks[i][j];
                    ww = w0 - (u[i + 1][j] - u[i][j]) / dx * th(bb) / k0 - (v[i][j + 1] - v[i][j]) / dx * th(bb) / k0;
                    h0[i][j] = h0[i][j] + ww * dt;
                }
            }
        }
    }
    cout << "Xg" << Xg[j] << endl;
    cout << "Lm" << Lm[i][j] << endl;
    cout << "Lm0" << Lm0[i][j] << endl;
    cout << "h" << h[i][j] << endl;
    cout << "u" << u[i][j] << endl;
    cout << "v" << v[i][j] << endl;
    cout << "hs" << hs << endl;
    cout << "Ks" << Ks[i][j] << endl;
    cout << "Cs" << Cs[i][j] << endl;
    cout << "bb" << bb << endl;
    cout << "k0" << k0 << endl;
    cout << "ww" << ww << endl;
    cout << "h0" << h0[i][j] << endl;
    no = QString::number(h0[i][j],'g',15);
    return h0[i][j];

}
void MainWindow::on_pushButton_clicked()
{
    double T0,w00,rs,j0_0,hh;
    T0 = (ui ->lineEdit->text()).toDouble();
    w00 = (ui ->lineEdit_2->text()).toDouble();
    hh = (ui ->lineEdit_3->text()).toDouble();
    rs = (ui ->lineEdit_4->text()).toDouble();
    j0_0= (ui ->lineEdit_5->text()).toDouble();
    cout << Zaph();
    cout << ZapCs();
    cout << ZapKs();
    cout << Zapfx();
    cout << Zapfy();
    cout << Zaph0();
    cout << Zapv();
    cout << Zapu();
    cout << ZapLm();
    cout << ZapLm0();
    cout << th;
    cout << KW(T0, w00, hh);
    cout << plan(h0);
    cout << MAX(T0);
    ui->label_10->setText(no);
    ui->label_11->setText(nolab);
    ui->label_12->setText(n);

}

