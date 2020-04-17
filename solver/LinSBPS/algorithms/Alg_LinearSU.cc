/*!
 * \author Ruben Martins - ruben@sat.inesc-id.pt
 *
 * @section LICENSE
 *
 * Open-WBO, Copyright (c) 2013-2017, Ruben Martins, Vasco Manquinho, Ines Lynce
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

#include "Alg_LinearSU.h"

using namespace openwbo;

/************************************************************************************************
 //
 // Linear Search Algorithm with Boolean Multilevel Optimization (BMO)
 //
 ************************************************************************************************/

/*_________________________________________________________________________________________________
  |
  |  bmoSearch : [void] ->  [void]
  |
  |  Description:
  |
  |    Linear search algorithm with lexicographical optimization.
  |
  |  For further details see:
  |    * Joao Marques-Silva, Josep Argelich, Ana Graça, Ines Lynce: Boolean
  |      lexicographic optimization: algorithms & applications. Ann. Math.
  |      Artif. Intell. 62(3-4): 317-343 (2011)
  |
  |  Post-conditions:
  |    * 'lbCost' is updated.
  |    * 'ubCost' is updated.
  |    * 'nbSatisfiable' is updated.
  |    * 'nbCores' is updated.
  |
  |________________________________________________________________________________________________@*/

StatusCode LinearSU::bmoSearch() {

  assert(orderWeights.size() > 0);
  lbool res = l_True;

  initRelaxation();

  uint64_t currentWeight = orderWeights[0];
  uint64_t minWeight = orderWeights[orderWeights.size() - 1];
  int posWeight = 0;

  vec<vec<Lit>> functions;
  vec<int> weights;

  solver = rebuildBMO(functions, weights, currentWeight);

  bool cardinality_used = false;

  uint64_t localCost = 0;
  ubCost = 0;

  for (;;) {

    vec<Lit> dummy;
    // Do not use preprocessing for linear search algorithm.
    // NOTE: When preprocessing is enabled the SAT solver simplifies the
    // relaxation variables which leads to incorrect results.
    res = searchSATSolver(solver, dummy);

    if (res == l_True) {
      nbSatisfiable++;

      uint64_t newCost = computeCostModel(solver->model, currentWeight);
      if (currentWeight == minWeight) {
        // If current weight is the same as the minimum weight, then we are in
        // the last lexicographical function.
        saveModel(solver->model);
        printf("o %" PRId64 "\n", newCost + lbCost + off_set); //pls
        printf("c t = %" PRId64 "\n", (time(0) - start_time));

        ubCost = newCost + lbCost;
      } else {
        if (verbosity > 0){
          printf("c BMO-UB : %-12" PRIu64 "\t (Function %d/%d)\n", newCost,
                 posWeight + 1, (int)orderWeights.size());
          printf("c t = %" PRId64 "\n", (time(0) - start_time));
        }
      }

      if (newCost == 0 && currentWeight == minWeight) {
        // Optimum value has been found.
        printAnswer(_OPTIMUM_);
        return _OPTIMUM_;
      } else {

        if (newCost == 0) {

          functions.push();
          new (&functions[functions.size() - 1]) vec<Lit>();
          objFunction.copyTo(functions[functions.size() - 1]);

          localCost = newCost;
          weights.push(localCost / currentWeight);

          posWeight++;
          currentWeight = orderWeights[posWeight];
          localCost = 0;

          delete solver;
          solver = rebuildBMO(functions, weights, currentWeight);

          if (verbosity > 0)
            printf("c LB : %-12" PRIu64 "\n", lbCost);
        } else {

          // Optimization of the current lexicographical function.
          if (localCost == 0)
            encoder.encodeCardinality(solver, objFunction,
                                      newCost / currentWeight - 1);
          else
            encoder.updateCardinality(solver, newCost / currentWeight - 1);


		  if (solver->_obj_vars.size() == 0) {
			  int *temp = new int[objFunction.size()];
			  for (int i = 0; i < objFunction.size(); i++) {
				  temp[i] = var(objFunction[i]);
			  }
			  std::sort(&temp[0], &temp[0] + objFunction.size());
			  for (int i = 0; i < objFunction.size(); i++) {
				  solver->_obj_vars.push(temp[i]);
			  }

			  delete[] temp;
		  }




          localCost = newCost;
        }
      }
    } else {
      nbCores++;

      if (currentWeight == minWeight) {
        // There are no more functions to be optimized.

        if (model.size() == 0) {
          assert(nbSatisfiable == 0);
          // If no model was found then the MaxSAT formula is unsatisfiable
          printAnswer(_UNSATISFIABLE_);
          return _UNSATISFIABLE_;
        } else {
          printAnswer(_OPTIMUM_);
          return _OPTIMUM_;
        }
      } else {

        // The current lexicographical function has been optimize. Go to the
        // next lexicographical function.
        functions.push();
        new (&functions[functions.size() - 1]) vec<Lit>();
        objFunction.copyTo(functions[functions.size() - 1]);

        weights.push(localCost / currentWeight);
        lbCost += localCost;

        posWeight++;
        currentWeight = orderWeights[posWeight];
        localCost = 0;

        delete solver;
        solver = rebuildBMO(functions, weights, currentWeight);

        if (verbosity > 0)
          printf("c LB : %-12" PRIu64 "\n", lbCost);
      }
    }
  }
}

#include <iostream>

/*_________________________________________________________________________________________________
  |
  |  normalSearch : [void] ->  [void]
  |
  |  Description:
  |
  |    Linear search algorithm.
  |
  |  For further details see:
  |    *  Daniel Le Berre, Anne Parrain: The Sat4j library, release 2.2. JSAT
  |       7(2-3): 59-6 (2010)
  |    *  Miyuki Koshimura, Tong Zhang, Hiroshi Fujita, Ryuzo Hasegawa: QMaxSAT:
  |       A Partial Max-SAT Solver. JSAT 8(1/2): 95-100 (2012)
  |
  |  Post-conditions:
  |    * 'ubCost' is updated.
  |    * 'nbSatisfiable' is updated.
  |    * 'nbCores' is updated.
  |
  |________________________________________________________________________________________________@*/

#include <algorithm>
#include <fstream>

StatusCode LinearSU::normalSearch() {

	//uint64_t bestUB_true = UINT64_MAX;

  lbool res = l_True;

  initRelaxation();
  
  if (solver != NULL) {
	  delete solver;
	  solver = NULL;
  }
  solver = rebuildSolver();

  bool processing_happened = (_lns_params._eproc && maxsat_formula->soft_clauses_original.size() > 0);

  bool cardinality_used = false;


  if (processing_happened && initial_model.size() > 0 ) {
	  model.clear();
	  for (int i = 0; i < initial_model.size(); i++) {
		  model.push(initial_model[i]);
	  }	  
  }

  vec<char> old_decision;

  if (_use_only_original_vars) {
	  for (int i = maxsat_formula->n_initial_vars; i < solver->nVars(); i++) {
		  old_decision.push(solver->decision[i]);
		  solver->decision[i] = 0;
	  }
  }

  while (res == l_True) {

    vec<Lit> dummy;
	
    // Do not use preprocessing for linear search algorithm.
    // NOTE: When preprocessing is enabled the SAT solver simplifies the
    // relaxation variables which leads to incorrect results.
    res = searchSATSolver(solver, dummy);	

	_use_only_original_vars = false;
	if (old_decision.size() > 0) {
		for (int i = 0; i < old_decision.size(); i++) {
			solver->decision[i + maxsat_formula->n_initial_vars] = old_decision[i];
		}
		old_decision.clear();
	}

    if (res == l_True) {
      nbSatisfiable++;
	  uint64_t newCost = computeCostModel(solver->model); 

//save model only if its true cost is better 
	  uint64_t trueCost;
	  bool wasBetter = false;
	  
	  if (processing_happened) {
		  trueCost = computeCostModel_beforePreproc(solver->model);		 

		  //printf("c c = %" PRId64 "\n", trueCost);

		  if (trueCost < bestUB_true) { //it could be that the new solution is worse than the previous one if processing is used
			  saveModel(solver->model);
			  bestUB_true = trueCost;
			  wasBetter = true;
			}
	  }
	  else {
		  saveModel(solver->model); //standard way, without proc is must be better
	  }


      if (maxsat_formula->getFormat() == _FORMAT_PB_) {
        // optimization problem
        if (maxsat_formula->getObjFunction() != NULL) {

			if (processing_happened) {
				printf("error: pb unexpected\n");
				std::ofstream output_emergency("problems.txt", std::ios::app);
				output_emergency << "Problem!\n";
                return _ERROR_;
			}

			printf("o %" PRId64 "\n", newCost + off_set);
			printf("c t = %" PRId64 "\n", (time(0) - start_time));			
        }
	  }
	  else {
		  if (processing_happened) {
			  if (wasBetter) { //do not display if it was worse!
				  printf("o %" PRId64 "\n", trueCost + off_set);
				  printf("c t = %" PRId64 "\n", (time(0) - start_time));
			  }
			  printf("c m %" PRId64 "\n", newCost);
		  }
		  else { //standard display
			  printf("o %" PRId64 "\n", newCost + off_set);
			  printf("c t = %" PRId64 "\n", (time(0) - start_time));
		  }
	  }

      if (newCost == 0) {
        // If there is a model with value 0 then it is an optimal model
        ubCost = newCost;

        if (maxsat_formula->getFormat() == _FORMAT_PB_ &&
            maxsat_formula->getObjFunction() == NULL) {
			
			if (processing_happened) {
				printf("c ??? pb and processing?\n");
				//exit(1);
			}

          printAnswer(_SATISFIABLE_);
          return _SATISFIABLE_;
        } else {
	          if (!processing_happened) {
				printAnswer(_OPTIMUM_);
				return _OPTIMUM_;
			}
			  else {
				  break;
			  }			
        }

      } else {
        if (maxsat_formula->getProblemType() == _WEIGHTED_) {
			if (!encoder.hasPBEncoding()) {
				printf("c first encoding\n");
				cardinality_used = false;
				encoder.encodePB(solver, objFunction, coeffs, newCost - 1);
				printf("c done with encoding\n");
			}
          else
            encoder.updatePB(solver, newCost - 1);
        } else {
          // Unweighted.
			if (!encoder.hasCardEncoding()) {
				encoder.encodeCardinality(solver, objFunction, newCost - 1);
				cardinality_used = true;
			}
          else
            encoder.updateCardinality(solver, newCost - 1);
        }

		if (solver->_obj_vars.size() == 0) {
			int *temp = new int[objFunction.size()];
			for (int i = 0; i < objFunction.size(); i++) {
				temp[i] = var(objFunction[i]);
			}
			std::sort(&temp[0], &temp[0] + objFunction.size());
			for (int i = 0; i < objFunction.size(); i++) {
				solver->_obj_vars.push(temp[i]);
			}
			
			delete[] temp;
		}

        ubCost = newCost;
      }
    } else {
      nbCores++;
      if (model.size() == 0) {
        assert(nbSatisfiable == 0);
        // If no model was found then the MaxSAT formula is unsatisfiable
        printAnswer(_UNSATISFIABLE_);
        return _UNSATISFIABLE_;
      } else {
		  if (processing_happened) {
			  printf("c done with this round\n");
			  break;
			}
		  else {
			  printAnswer(_OPTIMUM_);
			  return _OPTIMUM_;
		  }
      }
    }
  }

  if (processing_happened) {
	  if (cardinality_used) {
		  encoder.deleteCardinality();
	  }
	  else {
		  encoder.deletePB();
	  }
  } 

}

// Public search method
StatusCode LinearSU::search() {

	bool processing_happened = (_lns_params._eproc && maxsat_formula->soft_clauses_original.size() > 0);

  if (maxsat_formula->getProblemType() == _WEIGHTED_)
    is_bmo = isBMO();

  printConfiguration(is_bmo, maxsat_formula->getProblemType());

  if (maxsat_formula->getProblemType() == _WEIGHTED_) {
    if (bmoMode && is_bmo && !processing_happened) //make bmo work with preprocessing...but its unlikely bmo would take place after preprocessing
      bmoSearch();
    else
      normalSearch();
  } else
    normalSearch();
}

/************************************************************************************************
 //
 // Rebuild MaxSAT solver
 //
 ************************************************************************************************/

/*_________________________________________________________________________________________________
  |
  |  rebuildSolver : (minWeight : int)  ->  [Solver *]
  |
  |  Description:
  |
  |    Rebuilds a SAT solver with the current MaxSAT formula.
  |    If a weight is specified, then it only considers soft clauses with weight
  |    smaller than the specified weight.
  |    NOTE: a weight is specified in the 'bmo' approach.
  |
  |________________________________________________________________________________________________@*/
Solver *LinearSU::rebuildSolver(uint64_t min_weight) {

  vec<bool> seen;
  seen.growTo(maxsat_formula->nVars(), false);

  Solver *S = newSATSolver();

  for (int i = 0; i < maxsat_formula->nVars(); i++)
    newSATVariable(S);

  for (int i = 0; i < maxsat_formula->nHard(); i++)
    S->addClause(maxsat_formula->getHardClause(i).clause);

  for (int i = 0; i < maxsat_formula->nPB(); i++) {
    Encoder *enc = new Encoder(_INCREMENTAL_NONE_, _CARD_MTOTALIZER_,
                               _AMO_LADDER_, _PB_GTE_);

    // Make sure the PB is on the form <=
    // if (maxsat_formula->getPBConstraint(i)->_sign)
    //  maxsat_formula->getPBConstraint(i)->changeSign();
    assert(maxsat_formula->getPBConstraint(i)->_sign);

    enc->encodePB(S, maxsat_formula->getPBConstraint(i)->_lits,
                  maxsat_formula->getPBConstraint(i)->_coeffs,
                  maxsat_formula->getPBConstraint(i)->_rhs);

    delete enc;
  }

  for (int i = 0; i < maxsat_formula->nCard(); i++) {
    Encoder *enc = new Encoder(_INCREMENTAL_NONE_, _CARD_MTOTALIZER_,
                               _AMO_LADDER_, _PB_GTE_);

    if (maxsat_formula->getCardinalityConstraint(i)->_rhs == 1) {
      enc->encodeAMO(S, maxsat_formula->getCardinalityConstraint(i)->_lits);
    } else {

      enc->encodeCardinality(S,
                             maxsat_formula->getCardinalityConstraint(i)->_lits,
                             maxsat_formula->getCardinalityConstraint(i)->_rhs);
    }

    delete enc;
  }

  vec<Lit> clause;
  for (int i = 0; i < maxsat_formula->nSoft(); i++) {
    if (maxsat_formula->getSoftClause(i).weight < min_weight)
      continue;

    clause.clear();
    maxsat_formula->getSoftClause(i).clause.copyTo(clause);

    for (int j = 0; j < maxsat_formula->getSoftClause(i).relaxation_vars.size();
         j++) {
      clause.push(maxsat_formula->getSoftClause(i).relaxation_vars[j]);
    }

    S->addClause(clause);
  }

  return S;
}

/*_________________________________________________________________________________________________
  |
  |  rebuildBMO : (functions : int)  ->  [Solver *]
  |
  |  Description:
  |
  |    Rebuilds a SAT solver with the current MaxSAT formula.
  |    Only considers soft clauses with the weight of the current
  |    lexicographical optimization weight ('currentWeight')
  |    For each function already computed in the BMO algorithm it encodes the
  |    respective cardinality constraint.
  |
  |________________________________________________________________________________________________@*/
Solver *LinearSU::rebuildBMO(vec<vec<Lit>> &functions, vec<int> &rhs,
                             uint64_t currentWeight) {

  assert(functions.size() == rhs.size());

  Solver *S = rebuildSolver(currentWeight);

  objFunction.clear();
  coeffs.clear();
  for (int i = 0; i < maxsat_formula->nSoft(); i++) {
    if (maxsat_formula->getSoftClause(i).weight == currentWeight) {
      objFunction.push(maxsat_formula->getSoftClause(i).relaxation_vars[0]);
      coeffs.push(maxsat_formula->getSoftClause(i).weight);
    }
  }

  for (int i = 0; i < functions.size(); i++)
    encoder.encodeCardinality(S, functions[i], rhs[i]);

  return S;
}

/************************************************************************************************
 //
 // Other protected methods
 //
 ************************************************************************************************/

/*_________________________________________________________________________________________________
  |
  |  initRelaxation : (objective : vec<Lit>&) (weights : vec<int>&)  ->  [void]
  |
  |  Description:
  |
  |    Initializes the relaxation variables by adding a fresh variable to the
  |    'relaxationVars' of each soft clause.
  |
  |  Post-conditions:
  |    * 'objFunction' contains all relaxation variables that were added to soft
  |       clauses.
  |    * 'coeffs' contains the weights of all soft clauses.
  |
  |________________________________________________________________________________________________@*/
void LinearSU::initRelaxation() {

	//in case we did processing ensure they are empty
	objFunction.shrink(objFunction.size());
	coeffs.shrink(coeffs.size());

	for (int i = 0; i < maxsat_formula->nSoft(); i++) {
		Lit l = maxsat_formula->newLiteral();
		maxsat_formula->getSoftClause(i).relaxation_vars.push(l); //this seems harmless for the proc version, could be unnecessary but does not hurt? todo check
		objFunction.push(l);
		coeffs.push(maxsat_formula->getSoftClause(i).weight);
	}
	
}

// Print LinearSU configuration.
void LinearSU::print_LinearSU_configuration() {
  printf("c |  Algorithm: %23s                                             "
         "                      |\n",
         "LinearSU");

  if (maxsat_formula->getProblemType() == _WEIGHTED_) {
    if (bmoMode)
      printf("c |  BMO strategy: %20s                      "
             "                                             |\n",
             "On");
    else
      printf("c |  BMO strategy: %20s                      "
             "                                             |\n",
             "Off");

    if (bmoMode) {
      if (is_bmo)
        printf("c |  BMO search: %22s                      "
               "                                             |\n",
               "Yes");
      else
        printf("c |  BMO search: %22s                      "
               "                                             |\n",
               "No");
    }
  }
}
