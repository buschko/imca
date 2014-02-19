/*
 * lp.h
 *
 * Interface for LP solvers
 *
 *  Created on: Jul 17, 2013
 *      Author: harold
 */

#ifndef LP_H_
#define LP_H_

#include "real.h"

// choose between SoPlex and LP_solve
#ifdef __SOPLEX__

#include "soplex.h"
//using namespace soplex;
using soplex::SoPlex;
using soplex::DSVector;
using soplex::LPRow;
using soplex::LPCol;
using soplex::Real;

#else

#include "lp_lib.h"

#endif

#include <vector>
#include <map>

#ifdef __SOPLEX__
	typedef SoPlex LPModel;
#else
	typedef lprec* LPModel;
#endif //__SOPLEX__

class LPObjective {
public:
	LPObjective(unsigned long nrCols);

	virtual ~LPObjective() {};

	//bound = upper bound (lower = 0.0)
	void setCol(unsigned long index, Real value, Real bound);

	unsigned long getCols() { return m_cols; }

	void addToModel(LPModel model);
private:
	std::vector<unsigned long> m_indices;
	std::vector< std::pair<Real,Real> > m_values;

	unsigned long m_cols;
};

class LPConstraint {
public:
	enum Type {
		LEQ,
		GEQ
	};

	LPConstraint();

	virtual ~LPConstraint();

	void setValue(Real value);
	void setType(Type type);

	//index must be larger than of any previous call to setCol
	void setCol(unsigned long index, Real value);

	void addToModel(LPModel model);
private:
	std::vector<unsigned long> m_indices;
	std::vector<Real> m_values;

	Real m_value;
	Type m_type;
};

// Wrapper for LP solvers, designed specifically for LRA, so reuseability low
class LP {
public:
	// Last index of cols assumed to be the objective
	// Rows only for constraints, objective will be generated
	LP(unsigned long rows, unsigned long cols, bool max);

	virtual ~LP();

	void setObj(LPObjective objective);
	void addRow(LPConstraint constraint);

	void buildModel();

	bool solve();

	Real getObjective();
	Real getPrimal(unsigned long col);
private:
	LPModel m_model;
	unsigned long m_rows;
	unsigned long m_cols;

	bool m_maximize;

	Real m_result;
	std::vector< Real > m_primals;

	LPObjective m_objective;
	std::vector< LPConstraint > m_constraints;
};


#endif /* LP_H_ */
