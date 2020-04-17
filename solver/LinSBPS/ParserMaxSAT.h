/*!
 * \author Ruben Martins - ruben@sat.inesc-id.pt
 *
 * @section LICENSE
 *
 * MiniSat,  Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
 *           Copyright (c) 2007-2010, Niklas Sorensson
 * Open-WBO, Copyright (c) 2013-2015, Ruben Martins, Vasco Manquinho, Ines Lynce
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 */

#ifndef ParserMaxSAT_h
#define ParserMaxSAT_h

#include <stdio.h>

#include "MaxSATFormula.h"
#include "core/SolverTypes.h"
#include "utils/ParseUtils.h"

using NSPACE::mkLit;
using NSPACE::StreamBuffer;

namespace openwbo {

//=================================================================================================
// DIMACS Parser:

template <class B> static uint64_t parseWeight(B &in) {
  uint64_t val = 0;
  while ((*in >= 9 && *in <= 13) || *in == 32)
    ++in;
  if (*in < '0' || *in > '9')
    fprintf(stderr, "PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
  while (*in >= '0' && *in <= '9')
    val = val * 10 + (*in - '0'), ++in;
  return val;
}

template <class B, class MaxSATFormula>
static uint64_t readClause(B &in, MaxSATFormula *maxsat_formula,
                           vec<Lit> &lits) {
  int parsed_lit, var;
  int64_t weight = 1;
  lits.clear();
  if (maxsat_formula->getProblemType() == _WEIGHTED_)
    weight = parseWeight(in);
  assert(weight > 0);

  for (;;) {
    parsed_lit = parseInt(in);
    if (parsed_lit == 0)
      break;
    var = abs(parsed_lit) - 1;
    while (var >= maxsat_formula->nVars())
      maxsat_formula->newVar();
    lits.push((parsed_lit > 0) ? mkLit(var) : ~mkLit(var));
  }
  return weight;
}

template <class B, class MaxSATFormula>
static void parseMaxSAT(B &in, MaxSATFormula *maxsat_formula) {
  vec<Lit> lits;
  uint64_t hard_weight = UINT64_MAX;
  for (;;) {
    skipWhitespace(in);
    if (*in == EOF)
      break;
    else if (*in == 'p') {
      if (eagerMatch(in, "p cnf")) {
        parseInt(in); // Variables
        parseInt(in); // Clauses
      } else if (eagerMatch(in, "wcnf")) {
        maxsat_formula->setProblemType(_WEIGHTED_);
        parseInt(in); // Variables
        parseInt(in); // Clauses
        if (*in != '\r' && *in != '\n') {
          hard_weight = parseWeight(in);
          maxsat_formula->setHardWeight(hard_weight);
        }
      } else
        printf("c PARSE ERROR! Unexpected char: %c\n", *in),
            printf("s UNKNOWN\n"), exit(_ERROR_);
    } else if (*in == 'c' || *in == 'p')
      skipLine(in);
    else {
      uint64_t weight = readClause(in, maxsat_formula, lits);

      if (weight != hard_weight ||
          maxsat_formula->getProblemType() == _UNWEIGHTED_) {
        assert(weight > 0);
	
        // Updates the maximum weight of soft clauses.
        maxsat_formula->setMaximumWeight(weight);
        // Updates the sum of the weights of soft clauses.
        maxsat_formula->updateSumWeights(weight);
        maxsat_formula->addSoftClause(weight, lits);
      } else
        maxsat_formula->addHardClause(lits);
    }
  }
}
/*
template <class B, class MaxSATFormula>
static void parseMaxSAT(B &in, MaxSATFormula *maxsat_formula) {
  vec<Lit> lits;
  uint64_t hard_weight = UINT64_MAX;
  for (;;) {
    skipWhitespace(in);
    if (*in == EOF)
      break;
    else if (*in == 'p') {
      if (eagerMatch(in, "p cnf")) {
        parseInt(in); // Variables
        parseInt(in); // Clauses
      } else if (eagerMatch(in, "wcnf")) {
        maxsat_formula->setProblemType(_WEIGHTED_);
        parseInt(in); // Variables
        parseInt(in); // Clauses
        if (*in != '\r' && *in != '\n') {
          hard_weight = parseWeight(in);
          maxsat_formula->setHardWeight(hard_weight);
        }
      } else
        printf("c PARSE ERROR! Unexpected char: %c\n", *in),
            printf("s UNKNOWN\n"), exit(_ERROR_);
    } else if (*in == 'c' || *in == 'p')
      skipLine(in);
    else {
      uint64_t weight = readClause(in, maxsat_formula, lits);

      if (weight != hard_weight ||
          maxsat_formula->getProblemType() == _UNWEIGHTED_) {
        assert(weight > 0);
	
        // Updates the maximum weight of soft clauses.
        maxsat_formula->setMaximumWeight(weight);
        // Updates the sum of the weights of soft clauses.
        maxsat_formula->updateSumWeights(weight);
        maxsat_formula->addSoftClause(weight, lits);
      } else
        maxsat_formula->addHardClause(lits);
    }
  }
}*/

// Inserts problem into solver.
//
template <class MaxSATFormula>
static void parseMaxSATFormula(gzFile input_stream,
                               MaxSATFormula *maxsat_formula) {
  StreamBuffer in(input_stream);
  parseMaxSAT(in, maxsat_formula);
  if (maxsat_formula->getMaximumWeight() == 1)
    maxsat_formula->setProblemType(_UNWEIGHTED_);
  else
    maxsat_formula->setProblemType(_WEIGHTED_);

  // maxsat_formula->setInitialVars(maxsat_formula->nVars());
}

//add a heuristic algorithm to the mix

//check threshold applicability
//get the maximum number of digits
//get the starting point
//process the ORIGINAL formula with s

//if solved...

//save solution
//restore ORIGINAL formula
//process formula with s-1 if it is not already at 1 
//put the previous solution as a starting point
//solve

//process the formula with respect to the starting point

template <class MaxSATFormula> // = 500000
static bool isCandidateForIncreasePrecisionStrategy(MaxSATFormula *maxsat_formula, uint64_t threshold) {
	uint64_t sum(0);
	for (int i = 0; i < maxsat_formula->nSoft(); i++) {
		sum += maxsat_formula->getSoftClause(i).weight;
	} //I suppose sum_soft_weight would do the trick, check
	return sum >= threshold;
}

template <class MaxSATFormula>
static int getStartingPrecision(MaxSATFormula *maxsat_formula, double contribution_threshold) {
	vec<int> digit_counter(20);
	vec<uint64_t> digit_sum(20);


	uint64_t sum(0);
	for (int i = 0; i < maxsat_formula->nSoft(); i++) {
		sum += maxsat_formula->getSoftClause(i).weight;
	} //I suppose sum_soft_weight would do the trick, check



	int max_number_of_digits = -1;

	for (int i = 0; i < maxsat_formula->nSoft(); i++) {
		int number_of_digits = 0;
		uint64_t n = maxsat_formula->getSoftClause(i).weight;
		while (n > 0) {
			n /= 10;
			number_of_digits++;
		}

		max_number_of_digits = std::max(max_number_of_digits, number_of_digits);

		digit_counter[number_of_digits]++;
		digit_sum[number_of_digits] += maxsat_formula->getSoftClause(i).weight;
	}

	int starting_precision = max_number_of_digits;
	for (int d = 1; d < 20; d++) {
		if (digit_counter[d] > 0) {
			double contribution_factor = double(digit_sum[d]) / sum;
			if (contribution_factor >= contribution_threshold) {
				starting_precision = d;
				break;
			}
		}
	}

	return starting_precision;
}

template <class MaxSATFormula>
static void saveOriginalSoftClauses(MaxSATFormula *maxsat_formula) {
	assert(maxsat_formula->soft_clauses_original.size() == 0); //should be called only once for now

	for (int i = 0; i < maxsat_formula->nSoft(); i++) {
		Soft &sc = maxsat_formula->getSoftClause(i);
		maxsat_formula->soft_clauses_original.push();
		new (&maxsat_formula->soft_clauses_original[maxsat_formula->soft_clauses_original.size() - 1])
			Soft(sc.clause, sc.weight, sc.assumption_var, sc.relaxation_vars);
	}
}

static int gcd(uint64_t a, uint64_t b) {
	return b == 0 ? a : gcd(b, a % b);
}

template <class MaxSATFormula>
static void setFormulaToPrecision(MaxSATFormula *maxsat_formula, int precision) {
	assert(maxsat_formula->soft_clauses_original.size() != 0);
	assert(precision >= 1);

	uint64_t division_value = pow(10, precision - 1);
	printf("c div value: %d\n", division_value);

	vec<Soft> new_softies;
	uint64_t max_soft_weight_new = -1;
	uint64_t sum_soft_weights_new = 0;

	for (int i = 0; i < maxsat_formula->soft_clauses_original.size(); i++) {
		Soft &ori_soft = maxsat_formula->soft_clauses_original[i];
		uint64_t w = ori_soft.weight / division_value;
		
		if (w > 0) {
			sum_soft_weights_new += w;

			if (max_soft_weight_new < w) {
				max_soft_weight_new = w;
			}

			new_softies.push();
			new (&new_softies[new_softies.size() - 1])
				Soft(ori_soft.clause, w, ori_soft.assumption_var, ori_soft.relaxation_vars);
		}
	}

	//see if the common divisor can help
	int common_factor = new_softies[0].weight;
	for (int i = 1; i < new_softies.size(); i++) {
		common_factor = gcd(common_factor, new_softies[i].weight);
	}

	if (common_factor > 1) {
		max_soft_weight_new /= common_factor;
		sum_soft_weights_new /= common_factor;

		for (int i = 0; i < new_softies.size(); i++) {
			new_softies[i].weight /= common_factor;
		}
	}

	maxsat_formula->sum_soft_weight = sum_soft_weights_new;
	maxsat_formula->max_soft_weight = max_soft_weight_new;

	if (maxsat_formula->getMaximumWeight() == 1) {
		maxsat_formula->setProblemType(_UNWEIGHTED_);
	}
	else {
		maxsat_formula->setProblemType(_WEIGHTED_);
	}

	maxsat_formula->soft_clauses.growTo(new_softies.size());
	maxsat_formula->soft_clauses.shrink_(maxsat_formula->soft_clauses.size() - new_softies.size());
	maxsat_formula->n_soft = new_softies.size();

	for (int i = 0; i < maxsat_formula->n_soft; i++) {
		Soft &sc = maxsat_formula->getSoftClause(i);
		Soft &new_sc = new_softies[i];

		new_sc.clause.copyTo(sc.clause);
		sc.weight = new_sc.weight;
		sc.assumption_var = new_sc.assumption_var;
		new_sc.relaxation_vars.copyTo(sc.relaxation_vars);
	}

	printf("c proc %d done\n", division_value);
}

template <class MaxSATFormula> // = 500000
static void preprocessMaxSATformula(MaxSATFormula *maxsat_formula, int precision, uint64_t threshold) {
	/*
	precisionData returnData;
	returnData.maximum_precision = -1;
	returnData.starting_precision = -1;

	uint64_t sum(0);
	for (int i = 0; i < maxsat_formula->nSoft(); i++) {
		sum += maxsat_formula->getSoftClause(i).weight;
	} //I suppose sum_soft_weight would do the trick

	if (sum >= threshold) {
		printf("c im in\n");
		vec<int> digit_counter(20);
		vec<uint64_t> digit_sum(20);

		for (int i = 0; i < maxsat_formula->nSoft(); i++) {
			int number_of_digits = 0;
			uint64_t n = maxsat_formula->getSoftClause(i).weight;
			while (n > 0) {
				n /= 10;
				number_of_digits++;
			}

			returnData.maximum_precision = max(returnData.maximum_precision, number_of_digits);

			digit_counter[number_of_digits]++;
			digit_sum[number_of_digits] += maxsat_formula->getSoftClause(i).weight;
		}

		int digit_threshold = -1;
		for (int d = 1; d < 20; d++) {
			if (digit_counter[d] > 0) {
				double contribution_factor = double(digit_sum[d]) / sum;
				if (contribution_factor >= 0.05) {
					digit_threshold = d;
					break;
				}
			}
		}
		
		assert(digit_threshold != -1);
		uint64_t division_value = 1;
		for (int i = 1; i < digit_threshold; i++) {
			division_value *= 10;
		}

		printf("c div thresh and div value: %d, %d\n", digit_threshold, division_value);
		
		vec<Soft> new_softies;
		uint64_t max_soft_weight_new = -1;
		uint64_t sum_soft_weights_new = 0;

		//save original soft clauses before processing
		for (int i = 0; i < maxsat_formula->nSoft(); i++) {
			Soft &sc = maxsat_formula->getSoftClause(i);
			maxsat_formula->soft_clauses_original.push();
			new (&maxsat_formula->soft_clauses_original[maxsat_formula->soft_clauses_original.size() - 1])
				Soft(sc.clause, sc.weight, sc.assumption_var, sc.relaxation_vars);		
		}
		
		for (int i = 0; i < maxsat_formula->nSoft(); i++) {
			Soft &ori_soft = maxsat_formula->getSoftClause(i);
			uint64_t w = ori_soft.weight;
			w /= division_value;

			if (w > 0) {
				sum_soft_weights_new += w;
				
				if (max_soft_weight_new < w) {
					max_soft_weight_new = w;
				}

				new_softies.push();
				new (&new_softies[new_softies.size() - 1])
					Soft(ori_soft.clause, w, ori_soft.assumption_var, ori_soft.relaxation_vars);	
				//printf("c oerfirned\n");
			}
		}
		
		maxsat_formula->sum_soft_weight = sum_soft_weights_new;
		maxsat_formula->max_soft_weight = max_soft_weight_new;

		maxsat_formula->soft_clauses.shrink_(maxsat_formula->soft_clauses.size() - new_softies.size());
		maxsat_formula->n_soft = new_softies.size();
		for (int i = 0; i < maxsat_formula->n_soft; i++) {
			Soft &sc = maxsat_formula->getSoftClause(i);
			Soft &new_sc = new_softies[i];

			new_sc.clause.copyTo(sc.clause);
			sc.weight = new_sc.weight;
			sc.assumption_var = new_sc.assumption_var;
			new_sc.relaxation_vars.copyTo(sc.relaxation_vars);
		}
	}

	printf("c proc done %\n", );

	return returnData;*/
}

//=================================================================================================
} // namespace openwbo

#endif
