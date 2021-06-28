#include "PathSimulator2F.h"
#include "RandomNormalGenerator.h"
#include <algorithm>
#include <iostream>


PathSimulator2F::PathSimulator2F(const Pair &initial_values, const Vector& time_points, const ModelHeston& model) :
    _initial_values(initial_factors), _time_points(time_points), _model(model)
{}

PathSimulator2F::PathSimulator2F(const PathSimulator2F& path_simulator) :
    _initial_factors(path_simulator._initial_values), _time_points(path_simulator._time_points),
    _model(path_simulator._model)
{}


PathSimulator2F& PathSimulator2f::operator=(const PathSimulator2F& path_simulator) {
    if (!(this == &path_simulator)) {
        
        _initial_values = path_simulator._initial_values;
        _time_points = path_simulator._time_points;
    }
    return *this;					
}

PathSimulator2F::~PathSimulator2F() {

}


Vector PathSimulator2F:path() const {
    
}


SchemaTG::SchemaTG(Pair initial_values, const Vector& time_points, const ModelHeston& model)
    :PathSimulator2F(initial_values, time_points, model) 
{
    
}


/*
 * Descente de gradient
 */
double SchemaTG::2DRootSearch(double psi, double r_0, double eps, int maxIter) {
    int numIter = 0;
    double r_1, r_2, g_psi, g_psi_d, r_1_square;
    r_1 = r_0;
    do {
        numIter++;
        r_1_square = r_1 * r_1;
        double phi_r = (1 / sqrt(2 * M_PI)) * exp(-r_1_square / 2);
        double Phi_r = 0.5 * erfc(-r_1 * M_SQRT1_2);
        // Computation of g_psi
        g_psi = r_1 * phi_r + Phi_r * (1 + r_1_square) - (1 + psi) * pow(phi_r + r_1 * Phi_r, 2);
        // Computation of g_psi_d (derivative of g_psi)
        g_psi_d = 2 * phi_r + 2 * r_1 * Phi_r - 2 * (1 + psi) * Phi_r * (phi_r + r_1 * Phi_r);
        r_2 = r_1 - g_psi / g_psi_d;
        r_1 = r_2;

    } while (abs(g_psi) > eps && numIter < maxIter);
    
    return r_1;
}


Pair SchemaTG::nextStep(int current_index, Pair current_values) const
{
    double time_gap = _time_points[current_index + 1] - _time_points[current_index];

    Pair nextStep;
    double rho = _model.get_correlation();
    double vol_of_vol = _model.get_vol_of_vol();

    double z2 = RandomNormalGenerator::normalRandom();
    // computation of  nextStep Second = positive part of (mu+sigma*N(0,1))
    double k = _model.get_mean_reversion_speed();
    double theta = _model.get_mean_reversion_level();
    double m = theta + (current_values.second - theta) * std::exp(-k * time_gap);

    double s_2 = (current_values.second * pow(vol_of_vol, 2) * std::exp(-k * time_gap))
        / k * (1 - std::exp(-k * time_gap))
        + (theta * pow(vol_of_vol, 2) / 2 * k) * pow(1 - std::exp(-k * time_gap), 2);
    double psi = s_2 / pow(m, 2);

 
    // root search
    double r_0 = 10.0; // be careful to find the best r_0, eps, maxIter !!!
    double eps = 10e-15;
    int maxIter = 100;//(int)pow(10, 2);
    double r_psi = SchemaTG::rootSearch(psi, r_0, eps, maxIter);
 


    // calculating sigma and mu 
    double phi_r = (1 / sqrt(2 * M_PI)) * exp(-pow(r_psi, 2) / 2);
    double Phi_r = 0.5 * erfc(-r_psi * sqrt(0.5));
    double mu = m / (Phi_r + phi_r / r_psi);
    double sigma = mu / r_psi;
    nextStep.second = std::max(mu + sigma * z2, (double)0);

    double variance = time_gap * (0.5 * current_values.second + 0.5 * nextStep.second);
    double z1 = sqrt(variance) * RandomNormalGenerator::normalRandom();// normal random variable (m =0 ; var = time_gap*(0.5*current_values.second+0.5.nextStep.second)
    // nextStep.first refers to ln(X) 
    nextStep.first = current_values.first + rho / vol_of_vol *
        (nextStep.second - current_values.second -
            k * theta * time_gap)
        + (k * rho / vol_of_vol - 0.5) * variance + sqrt(1 - pow(rho, 2)) * z1;

    return nextStep;
}


SchemaTG* SchemaTG::clone() const
{
    return  new SchemaTG(*this); 
}


