/*
 * lp.cpp
 *
 *  Created on: Feb 19, 2014
 *      Author: hbruintjes
 */

#ifndef __LPSOLVER__
#error "LPSOLVER must be specified for LP classes"
#endif

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

LPObjective::LPObjective() :
	m_maxCol(0)
{
}

LPObjective::~LPObjective() {
}

//bound = upper bound (lower = 0.0)
void LPObjective::setCol(unsigned long index, Real value, Real bound) {
	if ((index+1) > m_maxCol)
		m_maxCol = index+1;
	if (value == 0.0)
		return;
	m_cols.push_back({index, value, bound});
}

void LPObjective::addToModel(LPModel& model) {
	//sort the vectors
	std::sort(m_cols.begin(), m_cols.end(), CmpObjCol());

#if __LPSOLVER__==_SOPLEX_
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
#elif __LPSOLVER__==_LPSOLVE_
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
		set_upbo(model, (*cols_it).index + 1, (*cols_it).bound);
		++cols_it;
	}

	delete [] row;
	delete [] colno;
#endif
}


LPConstraint::LPConstraint() :
	m_value(0.0),
	m_type(GEQ),
	m_maxCol(0)
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
	if ((index+1) > m_maxCol)
		m_maxCol = index + 1;
	if (value == 0.0)
		return;
	m_cols.push_back(Col(index, value));
}

void LPConstraint::addToModel(LPModel& model) {
	//sort the vectors
	std::sort(m_cols.begin(), m_cols.end(), CmpConstrCol());

#if __LPSOLVER__==_SOPLEX_
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
#elif __LPSOLVER__==_LPSOLVE_
	Real* row = new Real[m_cols.size()];
	int* colno = new int[m_cols.size()];
	int constr_type = (m_type == GEQ)?GE:LE;

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
LP::LP(bool max, Real delta) :
	m_maximize(max), m_result(0.0), m_objective() {
#if __LPSOLVER__==_SOPLEX_
	// Default constructor works
	//m_model = SoPlex();
	if (max) {
		m_model.changeSense(soplex::SPxLP::MINIMIZE);
	} else {
		m_model.changeSense(soplex::SPxLP::MAXIMIZE);
	}
	m_model.setDelta(delta);
#elif __LPSOLVER__==_LPSOLVE_
	m_model = make_lp(0, 0);
	set_add_rowmode(m_model, TRUE);
	set_sense(m_model, max?FALSE:TRUE);
#endif
}

LP::~LP() {
#if __LPSOLVER__==_SOPLEX_
	// Default destructor works
	//delete m_model;
#elif __LPSOLVER__==_LPSOLVE_
	delete_lp(m_model);
#endif
}

void LP::addRow(LPConstraint constraint) {
	m_constraints.push_back(constraint);
}

void LP::printModel() {
	buildModel();

	if (m_maximize) {
		printf("Minimize ");
	} else {
		printf("Maximize ");
	}

	printf("p = %.1f x%lu",m_objective.m_cols[0].value, m_objective.m_cols[0].index);
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
	unsigned int cols = m_objective.getMaxCol();
	std::vector< LPConstraint >::iterator constraint_it = m_constraints.begin();
	while(constraint_it != m_constraints.end()) {
		const unsigned long& max = (*constraint_it).getMaxCol();
		if (cols < max)
			cols = max;
		++constraint_it;
	}
	unsigned int rows = m_constraints.size();

	// Clear the model
#if __LPSOLVER__==_SOPLEX_
	m_model.clear();
#elif __LPSOLVER__==_LPSOLVE_
	delete_lp(m_model);
	m_model = make_lp(rows, cols);
	set_sense(m_model, m_maximize?FALSE:TRUE);

	//TODO: below setting rowmode triggers a bug in solving method
	//set_add_rowmode(m_model, TRUE);
#endif

	// Build the objective
	m_objective.addToModel(m_model);

	// Build the constraints
	constraint_it = m_constraints.begin();
	while(constraint_it != m_constraints.end()) {
		(*constraint_it).addToModel(m_model);
		++constraint_it;
	}
}

bool LP::solve() {
	m_primals.clear();
	m_result = 0.0;

	buildModel();
	//printModel();

#if __LPSOLVER__==_SOPLEX_
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
#elif __LPSOLVER__==_LPSOLVE_
	// TODO: This hack prevents some segfault inside lp_solve
	// if row mode is enabled
	//LPModel cp = copy_lp(m_model);
	//delete_lp(m_model);
	//m_model = cp;

	set_verbose(m_model, IMPORTANT);

	int solve_res = ::solve(m_model);

	// requery row (and col) counts, they might change as a result of solving
	int lprows = get_Nrows(m_model);
	int lpcols = get_Ncolumns(m_model);

	if (solve_res == OPTIMAL) {
		m_result = get_objective(m_model);
		Real* primals = new Real[1 + lprows + lpcols];
		get_primal_solution(m_model, primals);
		for(unsigned long i = 1 + lprows; i < 1 + lprows + lpcols; i++) {
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

