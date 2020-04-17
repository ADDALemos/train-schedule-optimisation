/*!
 * \author Ruben Martins - ruben@sat.inesc-id.pt
 *
 * @section LICENSE
 *
 * MiniSat,  Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
 *           Copyright (c) 2007-2010, Niklas Sorensson
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

#include <cstddef>

#include "LNSparameters.h"

#include "utils/Options.h"
#include "utils/ParseUtils.h"
#include "utils/System.h"
#include <errno.h>
#include <signal.h>
#include <zlib.h>

#include <fstream>
#include <iostream>
#include <map>
#include <stdlib.h>
#include <string>
#include <vector>
#include <fstream>

#ifdef SIMP
#include "simp/SimpSolver.h"
#else
#include "core/Solver.h"
#endif

#include "MaxSAT.h"
#include "MaxTypes.h"
#include "ParserMaxSAT.h"
#include "ParserPB.h"

// Algorithms
#include "algorithms/Alg_LinearSU.h"
#include "algorithms/Alg_MSU3.h"
#include "algorithms/Alg_OLL.h"
#include "algorithms/Alg_PartMSU3.h"
#include "algorithms/Alg_WBO.h"

#define VER1_(x) #x
#define VER_(x) VER1_(x)
#define SATVER VER_(SOLVERNAME)
#define VER VER_(VERSION)

using NSPACE::cpuTime;
using NSPACE::OutOfMemoryException;
using NSPACE::IntOption;
using NSPACE::BoolOption;
using NSPACE::StringOption;
using NSPACE::IntRange;
using NSPACE::parseOptions;
using namespace openwbo;

//=================================================================================================

static MaxSAT *mxsolver;

bool g_should_print_at_the_end = true;

static void SIGINT_exit(int signum) {
  if (g_should_print_at_the_end) mxsolver->printAnswer(_UNKNOWN_);
  exit(_UNKNOWN_);
}

//=================================================================================================
// Main:

int main(int argc, char **argv) {


  printf(
	"c\nc LinSBPS\t an experimental Open-WBO-based MaxSAT solver by Emir Demirovic and Prof Peter J. Stuckey.\n");

  printf(
      "c\nc Open-WBO:\t a Modular MaxSAT Solver -- based on %s (%s version)\n",
      SATVER, VER);
  printf("c Version:\t 2017 -- Release: 2.0\n");
  printf("c Authors:\t Ruben Martins, Vasco Manquinho, Ines Lynce\n");
  printf("c Contributors:\t Miguel Neves, Saurabh Joshi, Mikolas Janota\n");
  printf("c Contact:\t open-wbo@sat.inesc-id.pt -- "
         "http://sat.inesc-id.pt/open-wbo/\nc\n");
  try {
    NSPACE::setUsageHelp("c USAGE: %s [options] <input-file>\n\n");

#if defined(__linux__)
    fpu_control_t oldcw, newcw;
    _FPU_GETCW(oldcw);
    newcw = (oldcw & ~_FPU_EXTENDED) | _FPU_DOUBLE;
    _FPU_SETCW(newcw);
    printf(
        "c WARNING: for repeatability, setting FPU to use double precision\n");
#endif

    BoolOption printmodel("Open-WBO", "print-model", "Print model.\n", true);

    IntOption verbosity("Open-WBO", "verbosity",
                        "Verbosity level (0=minimal, 1=more).\n", 0,
                        IntRange(0, 1));

    IntOption algorithm("Open-WBO", "algorithm",
                        "Search algorithm "
                        "(0=wbo,1=linear-su,2=msu3,3=part-msu3,4=oll,5=best)."
                        "\n",
                        5, IntRange(0, 5));

    IntOption partition_strategy("PartMSU3", "partition-strategy",
                                 "Partition strategy (0=sequential, "
                                 "1=sequential-sorted, 2=binary)"
                                 "(only for unsat-based partition algorithms).",
                                 2, IntRange(0, 2));

    IntOption graph_type("PartMSU3", "graph-type",
                         "Graph type (0=vig, 1=cvig, 2=res) (only for unsat-"
                         "based partition algorithms).",
                         2, IntRange(0, 2));

    BoolOption bmo("Open-WBO", "bmo", "BMO search.\n", true);

    IntOption cardinality("Encodings", "cardinality",
                          "Cardinality encoding (0=cardinality networks, "
                          "1=totalizer, 2=modulo totalizer).\n",
                          1, IntRange(0, 2));

    IntOption amo("Encodings", "amo", "AMO encoding (0=Ladder).\n", 0,
                  IntRange(0, 0));

    IntOption pb("Encodings", "pb", "PB encoding (0=SWC,1=GTE).\n", 1,
                 IntRange(0, 0));

    IntOption formula("Open-WBO", "formula",
                      "Type of formula (0=WCNF, 1=OPB).\n", 0, IntRange(0, 1));

    IntOption weight(
        "WBO", "weight-strategy",
        "Weight strategy (0=none, 1=weight-based, 2=diversity-based).\n", 2,
        IntRange(0, 2));

    BoolOption symmetry("WBO", "symmetry", "Symmetry breaking.\n", true);

    IntOption symmetry_lim(
        "WBO", "symmetry-limit",
        "Limit on the number of symmetry breaking clauses.\n", 500000,
        IntRange(0, INT32_MAX));


	BoolOption phase_saving_solution_based("Open-WBO", "phase-saving-solution-based", "Phase saving solution based.\n", false);
	BoolOption backjump_strat("Open-WBO", "backjump-strat", "Backjump strategy to sbps until first conflict.\n", false);
	BoolOption backjump_alternate("Open-WBO", "backjump-alternate", "Set objective variables to one in phase saving.\n", false);

	BoolOption alternate_sbps_ps("Open-WBO", "alternate-sbps-ps", "Alternate between sbps and ps.\n", false);
	BoolOption obj_phasing_zero("Open-WBO", "obj-ps-zero", "Set objective variables to zero in phase saving.\n", false);
	BoolOption obj_phasing_one("Open-WBO", "obj-ps-one", "Set objective variables to one in phase saving.\n", false);

	IntOption should_print_end("Open-WBO", "end-print", "indicates whether the solver should print the solution at the end\n", 1, IntRange(0, 1));

	IntOption neg_sbps("Open-WBO", "negsbps",
		"Search algorithm "
		"(0=wbo,1=linear-su,2=msu3,3=part-msu3,4=oll,5=best)."
		"\n",
		0, IntRange(0, 100));

	IntOption sbps_chance("Open-WBO", "sbps-chance",
		"Search algorithm "
		"(0=wbo,1=linear-su,2=msu3,3=part-msu3,4=oll,5=best)."
		"\n",
		0, IntRange(0, 100));

	IntOption eproc_thresh("Open-WBO", "eprocthreshold",
		"the sum of soft weights that trigger eproc.\n", 500000,
		IntRange(1, 50000000));

		IntOption eproc_contribution_factor(
			"Open-WBO", "eprocfactor",
			"fraction of contribution necessary to be qualified for a starting precision in eproc.\n", 5,
			IntRange(0, 100));

	BoolOption emir_preprocessing("Open-WBO", "eproc", "Reduce precision and increase with search.\n", false);

	
    parseOptions(argc, argv, true);

	g_should_print_at_the_end = should_print_end;

    double initial_time = cpuTime();
    MaxSAT *S = NULL;

	LNSparameters lns_params;

	

	lns_params._eproc = emir_preprocessing;
	lns_params._sbps = phase_saving_solution_based;
	lns_params._backjump_shift_strat = backjump_strat;
	lns_params._obj_phasing_zero = obj_phasing_zero;
	lns_params._obj_phasing_one = obj_phasing_one;
	lns_params._negative_sbps_chance = neg_sbps;
	lns_params._alternate_sbps_ps = alternate_sbps_ps;
	lns_params._alternate_backjumps = backjump_alternate;

	lns_params._sbps_chance = sbps_chance;

	if (lns_params._negative_sbps_chance > 0 && lns_params._sbps == false) {
		printf("sbps chance non zero but sbps is off\n");
		exit(1);
	}

    switch ((int)algorithm) {
    case _ALGORITHM_WBO_:
      S = new WBO(verbosity, weight, symmetry, symmetry_lim);
      break;

    case _ALGORITHM_LINEAR_SU_:
      S = new LinearSU(verbosity, bmo, cardinality, pb, lns_params);
      break;

    case _ALGORITHM_PART_MSU3_:
      S = new PartMSU3(verbosity, partition_strategy, graph_type, cardinality);
      break;

    case _ALGORITHM_MSU3_:
      S = new MSU3(verbosity);
      break;

    case _ALGORITHM_OLL_:
      S = new OLL(verbosity, cardinality);
      break;

    case _ALGORITHM_BEST_:
      break;

    default:
      printf("c Error: Invalid MaxSAT algorithm.\n");
      printf("s UNKNOWN\n");
      exit(_ERROR_);
    }

    signal(SIGXCPU, SIGINT_exit);
    signal(SIGTERM, SIGINT_exit);

    if (argc == 1) {
      printf("c Error: no filename.\n");
      printf("s UNKNOWN\n");
      exit(_ERROR_);
    }

    gzFile in = (argc == 1) ? gzdopen(0, "rb") : gzopen(argv[1], "rb");
    if (in == NULL)
      printf("c ERROR! Could not open file: %s\n",
             argc == 1 ? "<stdin>" : argv[1]),
          printf("s UNKNOWN\n"), exit(_ERROR_);

    MaxSATFormula *maxsat_formula = new MaxSATFormula();
	MaxSATFormula *maxsat_formula_eproc;

	int starting_precision = -1;

	time_t my_start_time = time(0);

    if ((int)formula == _FORMAT_MAXSAT_) {
      parseMaxSATFormula(in, maxsat_formula);

	  if (lns_params._eproc && maxsat_formula->getProblemType() == _WEIGHTED_) {

		  if (isCandidateForIncreasePrecisionStrategy(maxsat_formula, eproc_thresh)) {
			  double frac = double(eproc_contribution_factor) / 100;

			  starting_precision = getStartingPrecision(maxsat_formula, frac);
			  saveOriginalSoftClauses(maxsat_formula);
			  setFormulaToPrecision(maxsat_formula, starting_precision);
		  }
	  }

      maxsat_formula->setFormat(_FORMAT_MAXSAT_);
    } else {
      ParserPB *parser_pb = new ParserPB();
      parser_pb->parsePBFormula(argv[1], maxsat_formula);
      maxsat_formula->setFormat(_FORMAT_PB_);
    }
    gzclose(in);

    printf("c |                                                                "
           "                                       |\n");
    printf("c ========================================[ Problem Statistics "
           "]===========================================\n");
    printf("c |                                                                "
           "                                       |\n");

    if (maxsat_formula->getFormat() == _FORMAT_MAXSAT_)
      printf(
          "c |  Problem Format:  %17s                                         "
          "                          |\n",
          "MaxSAT");
    else
      printf(
          "c |  Problem Format:  %17s                                         "
          "                          |\n",
          "PB");

    if (maxsat_formula->getProblemType() == _UNWEIGHTED_)
      printf("c |  Problem Type:  %19s                                         "
             "                          |\n",
             "Unweighted");
    else
      printf("c |  Problem Type:  %19s                                         "
             "                          |\n",
             "Weighted");

    printf("c |  Number of variables:  %12d                                    "
           "                               |\n",
           maxsat_formula->nVars());
    printf("c |  Number of hard clauses:    %7d                                "
           "                                   |\n",
           maxsat_formula->nHard());
    printf("c |  Number of soft clauses:    %7d                                "
           "                                   |\n",
           maxsat_formula->nSoft());
    printf("c |  Number of cardinality:     %7d                                "
           "                                   |\n",
           maxsat_formula->nCard());
    printf("c |  Number of PB :             %7d                                "
           "                                   |\n",
           maxsat_formula->nPB());
    double parsed_time = cpuTime();

    printf("c |  Parse time:           %12.2f s                                "
           "                                 |\n",
           parsed_time - initial_time);
    printf("c |                                                                "
           "                                       |\n");

    if (algorithm == _ALGORITHM_BEST_) {
      assert(S == NULL);

      if (maxsat_formula->getProblemType() == _UNWEIGHTED_) {
        // Unweighted
        S = new PartMSU3(_VERBOSITY_MINIMAL_, _PART_BINARY_, RES_GRAPH,
                         cardinality);
        S->loadFormula(maxsat_formula);

        if (((PartMSU3 *)S)->chooseAlgorithm() == _ALGORITHM_MSU3_) {
          // FIXME: possible memory leak
          S = new MSU3(_VERBOSITY_MINIMAL_);
        }

      } else {
        // Weighted
        S = new OLL(_VERBOSITY_MINIMAL_, cardinality);
      }
    }

    if (S->getMaxSATFormula() == NULL)
      S->loadFormula(maxsat_formula);


    S->setPrintModel(printmodel);
    S->setInitialTime(initial_time);
	mxsolver = S;
    
	
	//insert heuristic here?

	printf("c var %d vs ini vars %d\n", maxsat_formula->nVars(), maxsat_formula->n_initial_vars);

	int n_ini_vars = maxsat_formula->n_initial_vars;


	while (1 == 1) {
		mxsolver->search(); 

		vec<lbool> previous_model;
		for (int i = 0; i < mxsolver->model.size(); i++) {
			previous_model.push(mxsolver->model[i]);
		}
		uint64_t oldUB = mxsolver->bestUB_true;

		starting_precision--;
		if (starting_precision > 0) {
			setFormulaToPrecision(maxsat_formula, starting_precision);

			for (int i = 0; i < maxsat_formula->soft_clauses.size(); i++) {
				maxsat_formula->soft_clauses[i].relaxation_vars.clear();
			}

			S = new LinearSU(verbosity, bmo, cardinality, pb, lns_params);
			S->setPrintModel(printmodel);
			S->setInitialTime(initial_time);
			S->start_time = my_start_time;
			mxsolver = S;
			mxsolver->loadFormula(maxsat_formula);

			//for (int i = 0; i < previous_model.size(); i++) {
			for(int i = 0; i < n_ini_vars; i++){
				mxsolver->model.push(previous_model[i]);
			}
			mxsolver->bestUB_true = oldUB;

			mxsolver->_use_only_original_vars = true;

		}
		else { //we done
			if (g_should_print_at_the_end) mxsolver->printAnswer(_UNKNOWN_);
			exit(_UNKNOWN_);
		}
	}
	
	

	//solve the problem
	//if the precision is not at the maximum
	//unload the current formula, preprocess the original formula with higher precision
	//load the formula
	//set the previous solution as a hot start
	//solve again

	//remember to forbid saying that a solution is optimal when this precision style is used
	//essentially, when you find the optimum solution, don't print it
	//let the print function at the end display the solution

  } catch (OutOfMemoryException &) {
    //sleep(1);
    printf("c Error: Out of memory.\n");

	if (mxsolver && mxsolver->model.size() > 0) {
		if (g_should_print_at_the_end) mxsolver->printAnswer(_UNKNOWN_);
		exit(_UNKNOWN_);
	}
	else {
		printf("s UNKNOWN\n");
		exit(_ERROR_);
	}

    
  }
}
