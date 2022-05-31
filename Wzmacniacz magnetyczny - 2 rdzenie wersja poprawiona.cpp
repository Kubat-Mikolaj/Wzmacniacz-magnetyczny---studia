#include <stdio.h>
#include <conio.h>
#include <iostream>                                                                            
#include <iomanip>                                                        
#include <fstream>                                                    
#include <cmath>
#define N 3
#define M 2

using namespace std;

// Dane:

double alfa11 = 98, alfa21 = 98, alfa1S = 140, alfa2S = 140;
double alfa[2], alfaS[2];
double a, FI[2];
double r11 = 3.5, r21 = 3.5, r1s = 4.6, r2s = 4.6;
double Rs = 3, ro = 30;
double C = 0.00002;
double w = 314.159265;
double us = 36;
double UM = 240;
double u1;

double tk = 1., T = 0.02;
int ilosc = 1;

// Dane krzywa magnesowania:

double a1 = 1, a2 = 50;
double a0 = 45;
double PS1 = 0.5, PS2 = 1.1;
double FI1 = 0.5, FI2 = 10;
double A0, A1, A2, A3;

/////////////////////////////////////////////////////////////////////
//Funkcje:

// Mnożenie macierzy a*b=c
void mm(double a[][M], double b[][M], double c[][M])
{
    int i, j, k;
    for (i = 0; i < M; i++)
        for (j = 0; j < M; j++)
        {
            c[i][j] = 0.;
            for (k = 0; k < M; k++)
                c[i][j] = c[i][j] + a[i][k] * b[k][j];
        }
}
///////////////////////////////////////////////////////

// Mnożenie macierzy na kolumne a*b=c
void mk(double a[][M], double b[], double c[])
{
    int i, j, k;
    for (i = 0; i < M; i++)
    {
        c[i] = 0.;
        for (j = 0; j < M; j++)
            c[i] = c[i] + a[i][j] * b[j];
    }
}
///////////////////////////////////////////////////////


void aproks(double y[], double fi[])
{
    double ps[2], psi[2];
    int i;
    psi[0] = y[0];
    psi[1] = y[1];
    ps[0] = fabs(psi[0]);
    ps[1] = fabs(psi[1]);
    for (i = 0; i < 2; i++)
    {
        if (ps[i] < PS1)
        {
            alfa[i] = a1;
            alfaS[i] = a1;
        }
        else
            if (ps[i] < PS2)
            {
                alfa[i] = A1 + 2. * A2 * ps[i] + 3. * A3 * ps[i] * ps[i];
                alfaS[i] = A0 / ps[i] + A1 + A2 * ps[i] + A3 * ps[i] * ps[i];
            }
            else
            {
                alfa[i] = a2;
                alfaS[i] = a2 - a0 / ps[i];
            }
        fi[i] = alfaS[i] * psi[i];
        FI[i] = fi[i];
    }
}

void dydt(double t, double y[], double dy[])
{
    //   psi1=y[0], psi2=y[1], uc=y[2]
    double r1, i1, is, rs, alfa1, alfaS, a1, a2, b, D[2][2], P[2][2], LAMBDA[2][2], delta, DPS[2];
    double fi[2], dp[2];
    u1 = UM * sin(w * t + 1.5);
    alfa1 = alfa11 * alfa21 / (alfa11 + alfa21);  // 7.5
    alfaS = alfa1S * alfa2S / (alfa1S + alfa2S);
    a1 = alfa[0] + alfa1 + alfaS;  // 7.8
    a2 = alfa[1] + alfa1 + alfaS;
    b = alfa1 - alfaS;
    delta = a1 * a2 - b * b; //  7.12
    P[0][0] = a2 / delta;
    P[0][1] = b / delta;
    P[1][0] = P[0][1];
    P[1][1] = a1 / delta;
    LAMBDA[0][0] = alfa1;  // 7.9
    LAMBDA[0][1] = -alfaS;
    LAMBDA[1][0] = -alfa1;
    LAMBDA[1][1] = -alfaS;
    mm(P, LAMBDA, D);  // mnożenie macierzy (P * LAMBDA = D)
    r1 = r11 + r21 + ro;  // 7.10
    rs = r1s + r2s + Rs;
    aproks(y, fi);
    i1 = (fi[0] - fi[1]) / 2.;  //7.17
    is = -(fi[0] + fi[1]) / 2.;
    DPS[0] = u1 - y[2] - r1 * i1;  // 7.10
    DPS[1] = us - rs * is;
    mk(D, DPS, dp);  //  mnorzenie macierzy na kolumne (D * DPS = dp)
    dy[0] = dp[0];
    dy[1] = dp[1];
    dy[2] = i1 / C;  //  7.15
}

// Metoda Rynge-Kutty a).
void runge_kutta(double t, double y[], double dy[]) {
    double h = 0.02 / 720.;
    double k1[N], k2[N], k3[N], n[N] = { 0.0 };
    fstream file;
    fstream file1;
    fstream file2;
    file.open("psi1.txt", ios::out | ios::trunc);
    file1.open("psi2.txt", ios::out | ios::trunc);
    file2.open("fi.txt", ios::out | ios::trunc);
    do {
        dydt(t, y, dy);
        for (int i = 0; i < N; i++)//pкtla dla k1
            k1[i] = h * dy[i];
        t += h / 2;
        for (int i = 0; i < N; i++)//petla dla k2
            y[i] = n[i] + k1[i] / 2.;
        dydt(t, y, dy);
        for (int i = 0; i < N; i++)//petla dla k2
            k2[i] = h * dy[i];
        t += h / 2;
        for (int i = 0; i < N; i++)//petla dla k3
            y[i] = n[i] + 2. * k2[i] - k1[i];
        dydt(t, y, dy);
        for (int i = 0; i < N; i++)//petla dla k3
            k3[i] = h * dy[i];
        for (int i = 0; i < N; i++)//zmiana y
            y[i] = n[i] + (k1[i] + 4. * k2[i] + k3[i]) / 6.;
        if (ilosc % 10 == 0)
        {
            file << setw(13) << t << setw(13) << FI[0] << endl;
            file1 << setw(13) << t << setw(13) << FI[1] << endl;
            file2 << setw(13) << t << setw(13) << -(FI[0] + FI[1]) / 2. << endl;
            cout << setw(13) << t << setw(13) << y[0] << setw(13) << y[1] << endl;
        }
        ilosc++;
        for (int i = 0; i < N; i++)n[i] = y[i];
    } while (t < tk);
    file.close();
}

int main()
{
    double y[N] = { 0.0 }, dy[N];
    double t = 0.0;
    double c0, c2, c3, d0, d1, d2, d3;
    // Metoda eliminacji gausa
    c0 = (FI1 - FI2) / (PS1 - PS2);
    c2 = (PS1 + PS2);
    c3 = PS1 * PS1 + PS1 * PS2 + PS2 * PS2;
    d0 = (c0 - a1) / (c2 - 2. * PS1);
    d1 = (c3 - 3. * PS1 * PS1) / (c2 - 2. * PS1);
    d2 = (c0 - a2) / (c2 - 2. * PS2);
    d3 = (c3 - 3. * PS2 * PS2) / (c2 - 2. * PS2);
    A3 = (d0 - d2) / (d1 - d3);
    A2 = d0 - A3 * d1;
    A1 = c0 - A2 * c2 - A3 * c3;
    A0 = FI1 - A1 * PS1 - A2 * PS1 * PS1 - A3 * PS1 * PS1 * PS1;
    runge_kutta(t, y, dy);
    return 0;
}