#ifndef HESTONMODEL_H
#define HESTONMODEL_H
#include <vector>

using Pair = std::pair<double, double>;

using namespace std;

class HestonModel
{
	public:
	
		HestonModel(double drift, double correlation, double mean_reversion_speed, double mean_reversion_level,
                            double vol_of_vol);

		Pair drift_pair(const double &time,const Pair & spot_variance) const;
		Pair diffusion_pair(const double &time,const Pair & spot_variance) const;
        double correlation() const;
        double vol_of_vol() const;
        double mean_reversion_speed() const;
        double mean_reversion_level() const;
        double drift() const;

    ~HestonModel() = default;

	private:

        
    		double _correlation;
    		double _drift;
    		double _mean_reversion_speed;
    		double _mean_reversion_level;
    		double _vol_of_vol;
};
#endif