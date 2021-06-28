
#include "ModelHeston.h"
#include <cmath>


ModelHeston::ModelHeston(double drift, double correlation, 
                double mean_reversion_speed, 
                double mean_reversion_level,
                double vol_of_vol): 
    drift(drift), _correlation(correlation), _mean_reversion_speed(mean_reversion_speed),
    _mean_reversion_level(mean_reversion_level), _vol_of_vol(vol_of_vol)
{}



/*
 * returns a pair (drift asset, drift diffusion)
 */ 
Pair ModelHeston::drift_pair(const double &time,const Pair & spot_variance) const{
    Pair drift_pair;
    drift_pair.first = _drift*spot_variance.first;
    drift_pair.second = _mean_reversion_speed*(_mean_reversion_level-spot_variance.second);
    return drift_pair;
}


/*
 * returns a pair (volatility asset, vol_of_vol)
 */
Pair ModelHeston::volatility_term(const double &time,const Pair & spot_variance) const{
	Pair volatility_pair;
    volatility_pair.first = sqrt(factors.second)*factors.first;
    volatility_pair.second = _vol_of_vol*sqrt(factors.second);
    return volatility_pair;
}


double ModelHeston::correlation() const
{
    return _correlation;
}


double ModelHeston::vol_of_vol() const {
    return _vol_of_vol;
}

double ModelHeston::mean_reversion_speed() const {
    return _mean_reversion_speed;
}


double ModelHeston::mean_reversion_level() const 
{

    return _mean_reversion_level;
}


double ModelHeston::drift() const 

{
    return _drift;
}


