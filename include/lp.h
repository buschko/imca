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

#define _SOPLEX_ 1
#define _LPSOLVE_ 2
#define _GLPK_ 3

#include "real.h"

// choose between SoPlex and LP_solve
#ifdef __LPSOLVER__

#if __LPSOLVER__!=_SOPLEX_ && __LPSOLVER__!=_LPSOLVE_ && __LPSOLVER__!=_GLPK_
#error "Only SOPLEX, LPSOLVE and GLPK solvers supported"
#endif


#if __LPSOLVER__==_SOPLEX_
#include "soplex.h"
//using namespace soplex;
using soplex::SoPlex;
using soplex::DSVector;
using soplex::LPRow;
using soplex::LPCol;
using soplex::Real;

#elif __LPSOLVER__==_LPSOLVE_

extern "C" {
#include "lp_lib.h"
}

#elif __LPSOLVER__==_GLPK_

extern "C" {
#include "glpk.h"
}

#endif

#include <vector>
#include <map>

#if __LPSOLVER__==_SOPLEX_
	typedef SoPlex LPModel;
#elif __LPSOLVER__==_LPSOLVE_
	typedef lprec* LPModel;
#elif __LPSOLVER__==_GLPK_
	typedef glp_prob* LPModel;
#endif //__LPSOLVER__==__SOPLEX__

class LPObjective {
	friend class LP;
public:
	struct Col {
		Col(unsigned long index_, Real value_, Real upperbound_, Real lowerbound_) :
			index(index_), value(value_), upperbound(upperbound_), lowerbound(lowerbound_) {}
		unsigned long index;
		Real value;
		Real upperbound;
		Real lowerbound;
	};
public:
	LPObjective();
	virtual ~LPObjective();

	//bound = upper bound (lower = 0.0)
	void setCol(unsigned long index, Real value, Real upperbound, Real lowerbound=0.0);

	unsigned long getMaxCol() const { return m_maxCol; }

	void addToModel(LPModel& model, unsigned long maxCol = 0);
private:
	std::vector< Col > m_cols;
	unsigned long m_maxCol;
};

class LPConstraint {
	friend class LP;
public:
	typedef std::pair<unsigned long, Real> Col;
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

	unsigned long getMaxCol() const { return m_maxCol; }

	void addToModel(LPModel& model);
private:
	std::vector< Col > m_cols;
	unsigned long m_maxCol;

	Real m_value;
	Type m_type;
};

// Wrapper for LP solvers, designed specifically for LRA, so reuseability low
class LP {
public:
	LP(bool max, Real delta=1e-6);
	virtual ~LP();

	LPObjective& getObj() { return m_objective; }
	//void setObj(LPObjective objective);
	void addRow(LPConstraint constraint);

	void printModel();

	bool solve();

	Real getObjective();
	Real getPrimal(unsigned long col);
private:
	void buildModel();

	LPModel m_model;

	bool m_maximize;

	Real m_result;
	std::vector< Real > m_primals;

	LPObjective m_objective;
	std::vector< LPConstraint > m_constraints;
};

#endif //__LPSOLVER__

#endif /* LP_H_ */
