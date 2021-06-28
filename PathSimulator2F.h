#ifndef PathSimulator2F_H
#define PathSimulator2F_H



#include "ModelHeston.h"
#include <cmath>
#include <vector>


using Vector = std::vector<double>;
using Vector_Pair = std::vector<pair<double, double> >;


using namespace std;
class PathSimulator2F
{
	public:
		PathSimulator2F( Pair initial_values, const Vector& time_points, const ModelHeston& model);
		PathSimulator2F(const PathSimulator2F& path_simulator);
		PathSimulator2F & operator=(const PathSimulator2F& path_simulator);
		

		
		
	
        //===============Destructor=======================
        virtual ~PathSimulator2F();

		//============return a vector of pair==============
		Vector path() const;

		//============== clone ===========================
		virtual PathSimulator2F*  clone() const = 0; 

        //==================Getters====================
        Pair getInitialValues() const;
        Vector getTimePoints() const;
        ModelHeston* GetModel() const;

	protected:
		
        //==========Members
        Pair _initial_values;		
		Vector _time_points;		
		const ModelHeston* _model;

        //===========get the pair at t+ delta
        virtual Pair nextStep(int current_index, Pair current_factors)  const = 0;


};



class SchemaTG : public PathSimulator2F
{
	public:
	
		SchemaTG(Pair initial_values,
			 const Vector& time_points,
			 const ModelHeston& model);

		
        static double 2DRootSearch(double psi, double r_0, double eps, int maxIter);

		Pair nextStep(int current_index, Pair current_values) const override;


		SchemaTG* clone() const override;

};
#endif
