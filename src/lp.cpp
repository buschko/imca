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
#include <stdio.h>

// Define the actual Real infinity here
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

void LPObjective::setCol(unsigned long index, Real value, Real upperbound, Real lowerbound) {
	if ((index+1) > m_maxCol)
		m_maxCol = index+1;
	m_cols.push_back(Col(index, value, upperbound, lowerbound));
}

void LPObjective::addToModel(LPModel& model, unsigned long maxCol) {
	//sort the vectors
	std::sort(m_cols.begin(), m_cols.end(), CmpObjCol());

#if __LPSOLVER__==_SOPLEX_
	DSVector dummycol(0);
	std::vector<Col>::const_iterator cols_it = m_cols.begin();
	if (maxCol == 0) maxCol = m_maxCol;

	for(unsigned long i = 0; i < maxCol; i++) {
		if ((*cols_it).index == i) {
			model.addCol(LPCol((*cols_it).value, dummycol, (*cols_it).upperbound, (*cols_it).lowerbound));
			++cols_it;
		} else {
			model.addCol(LPCol(0.0, dummycol, 0.0, 0.0));
		}
	}
#elif __LPSOLVER__==_LPSOLVE_
	Real* row = new Real[m_cols.size()];
	int* colno = new int[m_cols.size()];

	std::vector<Col>::const_iterator cols_it = m_cols.begin();
	unsigned int index = 0;
	while(cols_it != m_cols.end()) {
		row[index] = (*cols_it).value;
		colno[index] = (int)(*cols_it).index + 1;
		++cols_it;
		index++;
	}

	set_obj_fnex(model, (int)m_cols.size(), row, colno);

	cols_it = m_cols.begin();
	while(cols_it != m_cols.end()) {
		set_upbo(model, (int)(*cols_it).index + 1, (*cols_it).upperbound);
		set_lowbo(model, (int)(*cols_it).index + 1, (*cols_it).lowerbound);
		++cols_it;
	}

	delete [] row;
	delete [] colno;
#elif __LPSOLVER__==_GLPK_
	std::vector<Col>::const_iterator cols_it = m_cols.begin();
	while(cols_it != m_cols.end()) {
		glp_set_obj_coef(model, (int)(*cols_it).index + 1, (*cols_it).value);

		int type = GLP_FR;
		if ((*cols_it).lowerbound == -infinity) {
			if ((*cols_it).upperbound == infinity) {
				type = GLP_FR;
			} else {
				type = GLP_UP;
			}
		} else {
			if ((*cols_it).upperbound == infinity) {
				type = GLP_LO;
			} else if ((*cols_it).upperbound == (*cols_it).lowerbound) {
				type = GLP_FX;
			} else {
				type = GLP_DB;
			}
		}
		glp_set_col_bnds(model, (int)(*cols_it).index + 1, type, (*cols_it).lowerbound, (*cols_it).upperbound);

		++cols_it;
	}
#endif
}


LPConstraint::LPConstraint() :
	m_maxCol(0),
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
		colno[index] = (int)(*cols_it).first + 1;
		++cols_it;
		index++;
	}

	add_constraintex(model, (int)m_cols.size(), row, colno, constr_type, m_value);

	delete [] row;
	delete [] colno;
#elif __LPSOLVER__==_GLPK_
	Real* row = new Real[m_cols.size()+1];
	int* colno = new int[m_cols.size()+1];

	std::vector<Col>::const_iterator cols_it = m_cols.begin();
	unsigned int index = 0;
	while(cols_it != m_cols.end()) {
		row[index+1] = (*cols_it).second;
		colno[index+1] = (int)(*cols_it).first + 1;
		++cols_it;
		index++;
	}

	//TODO: inefficient row adding, because index otherwise unknown
	int rowindex = glp_get_num_rows(model);
	glp_add_rows(model, 1);
	glp_set_mat_row(model, rowindex, (int)m_cols.size(), colno, row);

	if (m_type == GEQ) {
		glp_set_row_bnds(model, rowindex, GLP_LO, m_value, m_value);
	} else {
		glp_set_row_bnds(model, rowindex, GLP_UP, m_value, m_value);
	}

	delete [] row;
	delete [] colno;
#endif
}


//TODO: max is for min or max lra, not LP solving, so it is inverted. Fix this in API
LP::LP(bool max, Real delta) :
	m_maximize(max), m_delta(delta), m_result(0.0), m_objective() {
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
	//set_add_rowmode(m_model, TRUE);
	set_sense(m_model, max?FALSE:TRUE);
	set_epsd(m_model, delta);
#elif __LPSOLVER__==_GLPK_
	m_model = glp_create_prob();
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
	unsigned long cols = m_objective.getMaxCol();
	std::vector< LPConstraint >::iterator constraint_it = m_constraints.begin();
	while(constraint_it != m_constraints.end()) {
		const unsigned long& max = (*constraint_it).getMaxCol();
		if (cols < max)
			cols = max;
		++constraint_it;
	}

	// Clear the model
#if __LPSOLVER__==_SOPLEX_
	m_model.clear();
	if (m_maximize) {
		m_model.changeSense(soplex::SPxLP::MINIMIZE);
	} else {
		m_model.changeSense(soplex::SPxLP::MAXIMIZE);
	}
	m_model.setDelta(delta);
#elif __LPSOLVER__==_LPSOLVE_
	delete_lp(m_model);
	m_model = make_lp((int)m_constraints.size(), (int)cols);
	set_sense(m_model, m_maximize?FALSE:TRUE);

	//TODO: below setting rowmode triggers a bug in solving method
	//set_add_rowmode(m_model, TRUE);

	//TODO: set delta
#elif __LPSOLVER__==_GLPK_
	glp_delete_prob(m_model);
	m_model = glp_create_prob();
	glp_add_rows(m_model, (int)m_constraints.size());
	glp_add_cols(m_model, (int)cols);
	glp_set_obj_dir(m_model, m_maximize?GLP_MIN:GLP_MAX);
#endif

	// Build the objective
	m_objective.addToModel(m_model, cols);

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
	soplex::SPxSolver::Status solve_res =soplex::SPxSolver::ERROR;

	try {
		solve_res = m_model.solve();
	} catch(soplex::SPxStatusException& e) {
		printf("exception %s\n", e.what().c_str());
		return false;
	}

	if (solve_res == soplex::SPxSolver::OPTIMAL) {
		m_result = m_model.objValue();
		int cols = m_model.nCols();
		soplex::DVector primals(cols);
		m_model.getPrimal(primals);
		for(int i = 0; i < cols; i++) {
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
		for(int i = 1 + lprows; i < 1 + lprows + lpcols; i++) {
			m_primals.push_back(primals[i]);
		}
		delete [] primals;
	}

	return (solve_res == OPTIMAL);
#elif __LPSOLVER__==_GLPK_
	glp_smcp parms;
	glp_init_smcp(&parms);
	//parms.meth = GLP_DUAL;
	parms.msg_lev = GLP_MSG_ERR;

	//TODO: parameters for tolerances need to take delta into account
	int solve_res = glp_simplex(m_model, &parms);

	// requery col counts
	int lpcols = glp_get_num_cols(m_model);

	if (solve_res == 0) {
		m_result = glp_get_obj_val(m_model);

		for(int i = 0; i < lpcols; i++) {
			Real prim = glp_get_col_prim(m_model, i+1);
			m_primals.push_back(prim);
		}
	}

	return (solve_res == 0);
#endif
}

Real LP::getObjective() {
	return m_result;
}

Real LP::getPrimal(unsigned long col) {
	return m_primals[col];
}

