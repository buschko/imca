/*
 * lp.cpp
 *
 *  Created on: Feb 19, 2014
 *      Author: hbruintjes
 */

#include "lp.h"

#include <algorithm>
#include <assert.h>

// Definme the actual Real infinity here
const Real infinity = 1e100;

struct CmpObjCol {
	bool operator() (const LPObjective::Col& a, const LPObjective::Col& b) const {
		return a.index < b.index;
	}
};

struct CmpConstrCol {
	bool operator() (const LPConstraint::Col& a, const LPConstraint::Col& b) const {
		return a.first < b.first;
	}
};

LPObjective::LPObjective(unsigned long nrCols) :
	m_nrCols(nrCols)
{
}

//bound = upper bound (lower = 0.0)
void LPObjective::setCol(unsigned long index, Real value, Real bound) {
	m_cols.push_back({index, value, bound});
}

void LPObjective::addToModel(LPModel& model) {
	//sort the vectors
	std::sort(m_cols.begin(), m_cols.end(), CmpObjCol());

#ifdef __SOPLEX__
	DSVector dummycol(0);
	std::vector<Col>::const_iterator cols_it = m_cols.begin();
	for(unsigned long i = 0; i < m_nrCols; i++) {
		if ((*cols_it).index == i) {
			model.addCol(LPCol((*cols_it).value, dummycol, (*cols_it).bound, 0));
			++cols_it;
		} else {
			model.addCol(LPCol(0.0, dummycol, 0, 0));
		}
	}
#else
	Real* row = new Real[m_cols.size()];
	int* colno = new int[m_cols.size()];

	std::vector<Col>::const_iterator cols_it = m_cols.begin();
	unsigned int index = 0;
	while(cols_it != m_cols.end()) {
		row[index] = (*cols_it).value;
		colno[index] = (*cols_it).index + 1;
		++cols_it;
		index++;
	}
	set_obj_fnex(model, m_cols.size(), row, colno);

	cols_it = m_cols.begin();
	while(cols_it != m_cols.end()) {
		set_upbo(model, index, (*cols_it++).bound);
	}

	delete [] row;
	delete [] colno;
#endif
}


LPConstraint::LPConstraint() :
	m_value(0.0),
	m_type(GEQ)
{

}

LPConstraint::~LPConstraint() {

}

void LPConstraint::setValue(Real value) {
	m_value = value;
}
void LPConstraint::setType(LPConstraint::Type type) {
	m_type = type;
}

void LPConstraint::setCol(unsigned long index, Real value) {
	m_cols.push_back(Col(index, value));
}

void LPConstraint::addToModel(LPModel& model) {
	//sort the vectors
	//std::sort(m_cols.begin(), m_cols.end(), CmpConstrCol());

#ifdef __SOPLEX__
	DSVector row(m_cols.size());
	std::vector<Col>::const_iterator cols_it = m_cols.begin();
	while(cols_it != m_cols.end()) {
		row.add((*cols_it).first, (*cols_it).second);
		++cols_it;
	}

	if (m_type == GEQ) {
		model.addRow(LPRow(row, LPRow::GREATER_EQUAL, m_value));
	} else {
		model.addRow(LPRow(row, LPRow::LESS_EQUAL, m_value));
	}
#else
	Real* row = new Real[m_cols.size()];
	int* colno = new int[m_cols.size()];
	int constr_type = (m_type == GEQ)?2:1;

	std::vector<Col>::const_iterator cols_it = m_cols.begin();
	unsigned int index = 0;
	while(cols_it != m_cols.end()) {
		row[index] = (*cols_it).second;
		colno[index] = (*cols_it).first + 1;
		++cols_it;
		index++;
	}

	add_constraintex(model, (int)m_cols.size(), row, colno, constr_type, m_value);

	delete [] row;
	delete [] colno;
#endif
}


//TODO: max is for min or max lra, not LP solving, so it is inverted. Fix this in API
LP::LP(unsigned long rows, unsigned long cols, bool max) :
	m_rows(rows), m_cols(cols), m_maximize(max), m_result(0.0), m_objective(cols) {
#ifdef __SOPLEX__
	// Default constructor works
	//m_model = SoPlex();
	if (max) {
		m_model.changeSense(soplex::SPxLP::MINIMIZE);
	} else {
		m_model.changeSense(soplex::SPxLP::MAXIMIZE);
	}
#else
	m_model = make_lp(rows, cols);
	set_add_rowmode(m_model, TRUE);
	set_sense(m_model, max?FALSE:TRUE);
#endif
}

LP::~LP() {
#ifdef __SOPLEX__
	// Default destructor works
	//delete m_model;
#else
	delete_lp(m_model);
#endif
}

void LP::setObj(LPObjective objective) {
	m_objective = objective;
}

void LP::addRow(LPConstraint constraint) {
	m_constraints.push_back(constraint);
}

void LP::printModel() {
	buildModel();

	printf("p =  %.1f x%lu",m_objective.m_cols[0].value, m_objective.m_cols[0].index);
	for(unsigned int i = 1 ;i < m_objective.m_cols.size(); i++) {
		printf(" + %.1f x%lu",m_objective.m_cols[i].value, m_objective.m_cols[i].index);
	}
	printf("\n    such that\n");

	std::vector< LPConstraint >::iterator constraint_it = m_constraints.begin();
	while(constraint_it != m_constraints.end()) {
		LPConstraint& constr = *constraint_it;
		printf("    %.1f %s=  %.1f x%lu", constr.m_value, constr.m_type==LPConstraint::GEQ?"<":">",constr.m_cols[0].second, constr.m_cols[0].first);
		for(unsigned int i = 1 ;i < constr.m_cols.size(); i++) {
			printf(" + %.1f x%lu",constr.m_cols[i].second, constr.m_cols[i].first);
		}
		printf("\n");

		++constraint_it;
	}
	printf("\n");
}

void LP::buildModel() {
	// Some checking
	//assert (m_rows == m_constraints.size());
	assert (m_cols == m_objective.getCols());

	if (m_maximize) {
		printf("Minimize ");
	} else {
		printf("Maximize ");
	}

	// Build the objective
	m_objective.addToModel(m_model);


	// Build the constraints
#ifdef __SOPLEX__
#else
	// Adjust row number estimate
	resize_lp(m_model, m_constraints.size(), get_Ncolumns(m_model));
#endif
	std::vector< LPConstraint >::iterator constraint_it = m_constraints.begin();
	while(constraint_it != m_constraints.end()) {
		(*constraint_it).addToModel(m_model);
		++constraint_it;
	}
}

bool LP::solve() {
	//TODO: optimize the row cols number mess
	m_primals.clear();
	m_result = 0.0;

	buildModel();

#ifdef __SOPLEX__
	m_model.setDelta(1e-6);

	soplex::SPxSolver::Status solve_res =soplex::SPxSolver::UNKNOWN;

	try {
		solve_res = m_model.solve();
	} catch(soplex::SPxStatusException& e) {
		printf("exception %s\n", e.what().c_str());
		return false;
	}

	if (solve_res == soplex::SPxSolver::OPTIMAL) {
		m_result = m_model.objValue();
		soplex::DVector primals(m_cols);
		m_model.getPrimal(primals);
		for(unsigned long i = 0; i < m_cols; i++) {
			m_primals.push_back(primals[i]);
		}
	}

	return (solve_res == soplex::SPxSolver::OPTIMAL);
#else
	int solve_res = ::solve(m_model);

	if (solve_res == OPTIMAL) {
		m_result = get_objective(m_model);
		Real* primals = new Real[1 + m_rows + m_cols];
		for(unsigned long i = 1 + m_rows; i < 1 + m_rows + m_cols; i++) {
			m_primals.push_back(primals[i]);
		}
		delete [] primals;
	}

	return (solve_res == OPTIMAL);
#endif
}

Real LP::getObjective() {
	return m_result;
}

Real LP::getPrimal(unsigned long col) {
	return m_primals[col];
}
