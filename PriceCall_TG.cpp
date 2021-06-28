/*
 * Test unitaire pour tester le schema TG sur un call europeen
 */ 



int main(){
  unsigned num_sims = 100000;   // Number of simulated asset paths
  unsigned num_intervals = 1000;  // Number of intervals for the asset path to be sampled 

  double S_0 = 100.0;    // Initial spot price
  double K = 100.0;      // Strike price
  double r = 0.0319;     // Risk-free rate
  double v_0 = 0.010201; // Initial volatility 
  double T = 1.00;       // One year until expiry

  double rho = -0.7;     // Correlation of asset and volatility
  double kappa = 6.21;   // Mean-reversion rate
  double theta = 0.019;  // Long run average volatility
  double xi = 0.61;      // "Vol of vol"


}