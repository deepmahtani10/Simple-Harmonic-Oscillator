#include <stdio.h>
#include <cmath>
#include <fstream>
#include <iostream>

// No AI used, all 100% my brain.

double dt = 0.10;
double disp = 0.00;
double mass = 0.00;
double k_const = 15; // given as 15N/m
double v = 1.00;
double t = 0.00;
double v_increase = 0.00;
double disp_increase = 0.00;
double Energy_tot = 0.00;
double a_new = 0;
int i = 0;

void RK4(double &disp, double &v, double &t, double &dt)
{
    double v_init;
    double k1X, k2X, k3X, k4X = 0.00;
    double k1V, k2V, k3V, k4V = 0.00;
    double const omega = sqrt(k_const / mass);
    double a_RK4 = (k_const / mass) * disp;
    int j = 0;
    auto dXdt = [&](double disp, double v)
    {
        return (v);
    };
    auto dVdt = [&](double disp)
    {
        return (-pow(omega, 2) * disp);
    };

    std::ofstream file("SimpleHarmonicOut_RK4.txt", std::ios::out);
    if (!file.is_open())
    {
        std::cerr << "Error: could not open SimpleHarmonicOut.txt for writing.\n";
        return;
    }
    for (j = 0; j < 151; j++)
    {
        if (j == 0)
        {
            // v_init = v;
            file << "m = " << mass << " k = " << k_const << " x init =" << disp << " step = " << dt << " v_init =" << v << " time =" << t << "\n";
            printf("m = %lf, k = %lf, x = %lf, v = %lf, t = %lf, step = %lf, E = %lf\n", mass, k_const, disp, v, t, dt, Energy_tot);
            t += dt;
            continue;
        }
        // K1 according to RK4
        k1X = dt * dXdt(disp, v);
        k1V = dt * dVdt(disp);

        // k2 according to RK4
        k2X = dt * dXdt(disp + 0.5 * k1X, v + 0.5 * k1V);
        k2V = dt * dVdt(disp + (0.5 * k1X));

        // k3 according to RK4
        k3X = dt * dXdt((disp + 0.5 * k2X), (v + (0.5 * k2V)));
        k3V = dt * dVdt(disp + (0.5 * k2X));

        // k4
        k4X = dt * dXdt((disp + k3X), (v + k3V));
        k4V = dt * dVdt((disp + k3X));

        // Updating the values of x and v and a
        a_RK4 = dVdt(disp);
        v += (k1V + (2 * k2V) + (2 * k3V) + k4V) / 6.00;
        disp += (k1X + (2 * k2X) + (2 * k3X) + k4X) / 6.00;
        Energy_tot = (0.5 * mass * pow(v, 2)) + (0.5 * k_const * pow(disp, 2));

        // Need to write to file I'll fill it later
        file << "m = " << mass << " k = " << k_const << " x init =" << disp << " step = " << dt << " v_init =" << v << " time =" << t << "\n";
        printf("m = %lf, k = %lf, x = %lf, v = %lf, t = %lf, step = %lf, E = %lf\n", mass, k_const, disp, v, t, dt, Energy_tot);
        t += dt;
    }
}

void Euler()
{
    double a_euler = -(k_const / mass) * disp;
    std::ofstream file("SimpleHarmonicOut_e.txt", std::ios::out);
    if (!file.is_open())
    {
        std::cerr << "Error: could not open SimpleHarmonicOut.txt for writing.\n";
        return;
    }
    for (i = 0; i < 151; i++) // we give it a 150 iterations which gives us 15seconds of the simulation
    {
        if (i == 0)
        {
            Energy_tot = ((mass * v * v) + (k_const * disp * disp)) / 2.0;
            file << "m = " << mass << " k = " << k_const << " x init =" << disp << " step = " << dt << " v_init =" << v << " time =" << t << "\n";
            printf("m = %lf, k = %lf, x = %lf, v = %lf, t = %lf, step = %lf, E = %lf\n", mass, k_const, disp, v, t, dt, Energy_tot);
            t += dt;
        }
        else
        {
            v += (a_euler * dt);
            disp += v * dt;
            t += dt;
            Energy_tot = ((mass * v * v) + (k_const * disp * disp)) / 2.0;
            file << "m = " << mass << " k = " << k_const << " x init =" << disp << " step = " << dt << " v_init =" << v << " time =" << t << "\n";
            printf("m = %lf, k = %lf, x = %lf, v = %lf, t = %lf, step = %lf, E = %lf\n", mass, k_const, disp, v, t, dt, Energy_tot);
        }
    }
}
void Velver()
{
    double a = (k_const / mass) * disp;
    // Velocity Verlet
    std::ofstream file("SimpleHarmonicOut_v.txt", std::ios::out);
    if (!file.is_open())
    {
        std::cerr << "Error: could not open SimpleHarmonicOut.txt for writing.\n";
        return;
    }
    for (i = 0; i < 151; i++) // we give it a 150 iterations which gives us 15seconds of the simulation
    {
        if (i == 0)
        {
            Energy_tot = ((mass * v * v) + (k_const * disp * disp)) / 2.0;
            file << "m = " << mass << " k = " << k_const << " x init =" << disp << " step = " << dt << " v_init =" << v << " time =" << t << "\n";
            printf("m = %lf, k = %lf, x = %lf, v = %lf, t = %lf, step = %lf, E = %lf\n", mass, k_const, disp, v, t, dt, Energy_tot);
            t += dt;
        }
        else
        {
            disp_increase = disp + (dt * v) + (0.5 * a * dt * dt);
            a_new = -(k_const / mass) * disp_increase;
            v_increase = v + (0.5 * (a + a_new) * dt);
            Energy_tot = ((mass * v_increase * v_increase) + (k_const * disp_increase * disp_increase)) / 2.0;
            v = v_increase;
            disp = disp_increase;
            t += dt;
            file << "x = " << disp << " velocity = " << v << " time = " << t << "\n";
            printf("m = %lf, k = %lf, x = %lf, v = %lf, t = %lf, step = %lf E = %lf\n", mass, k_const, disp, v, t, dt, Energy_tot);
        }
    }
    std::cout << "\nFile written.";
}

int main()
{
    double original_v;
    std::cout << "\nPleasse type in the mass of the oscillator:";
    std::cin >> mass;
    std::cout << "\nPlease type in the speed of the oscillator, default is 1:";
    std::cin >> v;
    if (std::cin.fail())
    {
        v = 1.00;
        std::cout << "\nInvalid input, using v=1.0.";
    }
    else
    {
        original_v = v;
    }
    disp = 0.00;
    t = 0.00;
    v = original_v; // Resetting data to original working state intended to passed to functions to avoid global variable data corrpution on each function call
    Euler();
    disp = 0.00;
    t = 0.00;
    v = original_v;
    Velver();
    disp = 0.00;
    t = 0.00;
    v = original_v;
    RK4(disp, v, t, dt);
    return 0;
}

// Euler diverges, Verlet conserves E.