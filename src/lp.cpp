/*
 * lp.cpp
 *
 *  Created on: Feb 19, 2014
 *      Author: hbruintjes
 */

#include "lp.h"

#include <assert.h>

LPObjective::LPObjective(unsigned long nrCols) :
	m_cols(nrCols)
{
}

//bound = upper bound (lower = 0.0)
void LPObjective::setCol(unsigned long index, Real value, Real bound) {

}

void LPObjective::addToModel(LPModel model) {
#ifdef __SOPLEX__
	DSVector dummycol(0);
	std::vector<unsigned long>::const_iterator indic_it = m_indices.begin();
	std::vector< std::pair<Real,Real> >::iterator value_it = m_values.begin();
	for(unsigned long i = 0; i < m_cols; i++) {
		if (*indic_it == i) {
			model.addCol(LPCol((*value_it).first, dummycol, (*value_it).second, 0));
			++indic_it;
			++value_it;
		} else {
			model.addCol(LPCol(0.0, dummycol, 0, 0));
		}
	}
#else
	Real* row = new Real[m_values.size()];
	int* colno = new int[m_indices.size()];

	std::vector<unsigned long>::const_iterator indic_it = m_indices.begin();
	std::vector< std::pair<Real,Real> >::const_iterator value_it = m_values.begin();
	unsigned int index = 0;
	while(indic_it != m_indices.end()) {
		row[index] = (*value_it).first;
		colno[index] = *indic_it + 1;
		++indic_it;
		++value_it;
		index++;
	}
	set_obj_fnex(model, m_values.size(), row, colno);

	value_it = m_values.begin();
	while(indic_it != m_indices.end()) {
		set_upbo(model, index, (*value_it++).second);
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

//index must be larger than of any previous call to setCol
void LPConstraint::setCol(unsigned long index, Real value) {
	assert (*(m_indices.end()-1) < index);
	m_indices.push_back(index);
	m_values.push_back(value);
}

void LPConstraint::addToModel(LPModel model) {
#ifdef __SOPLEX__
	DSVector row(m_values.size());
	std::vector<unsigned long>::const_iterator indic_it = m_indices.begin();
	std::vector<Real>::const_iterator value_it = m_values.begin();
	while(indic_it != m_indices.end()) {
		row.add(*indic_it, *value_it);
		++indic_it;
		++value_it;
	}

	if (m_type == GEQ) {
		model.addRow(LPRow(row, LPRow::GREATER_EQUAL, m_value));
	} else {
		model.addRow(LPRow(row, LPRow::LESS_EQUAL, m_value));
	}
#else
	Real* row = new Real[m_values.size()];
	int* colno = new int[m_indices.size()];
	int constr_type = (m_type == GEQ)?2:1;

	std::vector<unsigned long>::const_iterator indic_it = m_indices.begin();
	std::vector<Real>::const_iterator value_it = m_values.begin();
	unsigned int index = 0;
	while(indic_it != m_indices.end()) {
		row[index] = *value_it;
		colno[index] = *indic_it + 1;
		++indic_it;
		++value_it;
		index++;
	}

	add_constraintex(model, (int)m_values.size(), row, colno, constr_type, m_value);

	delete [] row;
	delete [] colno;
#endif
}



LP::LP(unsigned long rows, unsigned long cols, bool max) :
	m_rows(rows), m_cols(cols), m_maximize(max), m_result(0.0), m_objective(0) {
#ifdef __SOPLEX__
	// Default constructor works
	//m_model = SoPlex();
#else
	m_model = make_lp(rows, cols);
	set_add_rowmode(m_model, TRUE);
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

void LP::buildModel() {
	// Some checking
	assert (m_rows == m_constraints.size());
	assert (m_cols == m_objective.getCols());

	// Build the objective
	m_objective.addToModel(m_model);

	// Build the constraints
	std::vector< LPConstraint >::iterator constraint_it = m_constraints.begin();
	while(constraint_it != m_constraints.end()) {
		(*constraint_it).addToModel(m_model);
		++constraint_it;
	}
}

bool LP::solve() {
#ifdef __SOPLEX__
	m_model.setDelta(1e-6);

	soplex::SPxSolver::Status solve_res = m_model.solve();
	m_primals.clear();

	if (solve_res == soplex::SPxSolver::OPTIMAL) {
		m_result = m_model.objValue();
	} else {
		m_result = 0.0;
		soplex::DVector primals(m_cols);
		m_model.getPrimal(primals);
		for(int i = 0; i < m_cols; i++) {
			m_primals.push_back(primals[i]);
		}
	}

	return (solve_res == soplex::SPxSolver::OPTIMAL);
#else
	int solve_res = ::solve(m_model);
	m_primals.clear();
	if (solve_res == OPTIMAL) {
		m_result = get_objective(m_model);
		Real* primals = new Real[1 + m_rows + m_cols];
		for(int i = 1 + m_rows; i < 1 + m_rows + m_cols; i++) {
			m_primals.push_back(primals[i]);
		}
		delete [] primals;
	} else {
		m_result = 0.0;
	}
	return (OPTIMAL == solve_res);
#endif
}

Real LP::getObjective() {
	return m_result;
}

Real LP::getPrimal(unsigned long col) {
	return m_primals[col];
}
