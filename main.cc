/*!
 * \author Ruben Martins - ruben@sat.inesc-id.pt
 *
 * @section LICENSE
 *
  * MiniSat,  Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
 *           Copyright (c) 2007-2010, Niklas Sorensson
 * Open-WBO, Copyright (c) 2013-2017, Ruben Martins, Vasco Manquinho, Ines Lynce
 * Open-WBO-Inc Copyright (c) 2018  Saurabh Joshi, Prateek Kumar, Ruben Martins, Sukrut Rao
 * TT-Open-WBO-Inc Copyright (c) 2019 Alexander Nadel
 * Timetabler Copyright (c) 2019 Alexandre Lemos, Pedro T Monteiro, Ines Lynce
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
#include <cstring>

#ifdef SIMP
#include "simp/SimpSolver.h"
#else

#include "core/Solver.h"

#endif

#include "MaxSAT.h"
#include "MaxSATFormula.h"
#include "MaxTypes.h"
#include "ParserMaxSAT.h"
#include "ParserPB.h"
#include "Torc.h"

// Algorithms
#include "algorithms/Alg_LinearSU.h"
#include "algorithms/Alg_LinearSU_Clustering.h"
#include "algorithms/Alg_LinearSU_Mod.h"
#include "algorithms/Alg_MSU3.h"
#include "algorithms/Alg_OLL.h"
#include "algorithms/Alg_OLL_Mod.h"
#include "algorithms/Alg_PartMSU3.h"
#include "algorithms/Alg_WBO.h"
#include "algorithms/Alg_OBV.h"
#include "algorithms/Alg_BLS.h"

//RapidJSON reader
#include "rapidjson/reader.h"
#include "rapidjson/document.h"
#include "rapidjson/istreamwrapper.h"
#include "rapidjson/stringbuffer.h"
#include "rapidjson/writer.h"



//Problem Domain
#include "problem/Instance.h"
#include "problem/Train.h"
#include "problem/Resource.h"
#include "problem/Route.h"
#include "problem/route_path.h"
#include "problem/route_section.h"



#define VER1_(x) #x
#define VER_(x) VER1_(x)
#define SATVER VER_(SOLVERNAME)
#define VER VER_(VERSION)


using NSPACE::BoolOption;
using NSPACE::IntOption;
using NSPACE::IntRange;
using NSPACE::OutOfMemoryException;
using NSPACE::StringOption;
using NSPACE::cpuTime;
using NSPACE::parseOptions;
using namespace openwbo;

//=================================================================================================

static MaxSAT *mxsolver;
#include "Test.h"

//Print Solver stats
void printSolverStats(MaxSATFormula*maxsat_formula,double initial_time);

int getVariableID(std::string varName,MaxSATFormula*maxsat_formula);

static void SIGINT_exit(int signum) {
    mxsolver->printAnswer(_UNKNOWN_);
    exit(_UNKNOWN_);
}


int size=-1;
Instance readJSONFile(char *);

void outputJSONFile(Instance instance);


void printPairs(
        vector<list<pair<route_section, route_section>, allocator<pair<route_section, route_section>>>, allocator<list<pair<route_section, route_section>, allocator<pair<route_section, route_section>>>>> &Pairs);

using namespace rapidjson;
using namespace std;

int main(int argc, char **argv) {
    double initial_time = cpuTime();

    try {
#if defined(__linux__)
        fpu_control_t oldcw, newcw;
        _FPU_GETCW(oldcw);
        newcw = (oldcw & ~_FPU_EXTENDED) | _FPU_DOUBLE;
        _FPU_SETCW(newcw);
        printf(
            "c WARNING: for repeatability, setting FPU to use double precision\n");
#endif


        BoolOption printmodel("Open-WBO", "print-model", "Print model.\n", true);

        IntOption num_tests("Test", "num_tests", "Number of tests\n", 0,
                            IntRange(0, 10000000));

        IntOption test_rhs("Test", "test_rhs",
                           "RHS for a custom encoding test\n", 0,
                           IntRange(0, 10000000));

        IntOption test_rhs2("Test", "test_rhs2",
                            "RHS for a custom encoding test for the second tree\n", 0,
                            IntRange(0, 10000000));

        IntOption test_nsoft("Test", "test_nsoft",
                             "Nsoft for a custom encoding test\n", 0,
                             IntRange(0, 10000000));

        IntOption test_join("Test", "test_join",
                            "Join for a custom encoding test\n", 0, IntRange(0, 1));

        IntOption verbosity("Open-WBO", "verbosity",
                            "Verbosity level (0=minimal, 1=more).\n", 0,
                            IntRange(0, 1));

        IntOption algorithm("Open-WBO", "algorithm",
                            "Search algorithm "
                                    "(_ALGORITHM_WBO_ = 0,_ALGORITHM_LINEAR_SU_,_ALGORITHM_MSU3_,"
                                    "_ALGORITHM_PART_MSU3_,_ALGORITHM_OLL_,_ALGORITHM_BEST_,_ALGORITHM_LSU_CLUSTER_,"
                                    "_ALGORITHM_LSU_MRSBEAVER_,_ALGORITHM_LSU_MCS_\n",
                            6, IntRange(0, 8));

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

        IntOption pb("Encodings", "pb", "PB encoding (0=SWC,1=GTE,2=GTECluster).\n",
                     1, IntRange(0, 2));

        IntOption formula("Open-WBO", "formula",
                          "Type of formula (0=WCNF, 1=OPB).\n", 0, IntRange(0, 1));
        IntOption cpu_lim("Open-WBO", "cpu-lim",
                          "Limit on CPU time allowed in seconds.\n", 0,
                          IntRange(0, INT_MAX));

        IntOption weight(
                "WBO", "weight-strategy",
                "Weight strategy (0=none, 1=weight-based, 2=diversity-based).\n", 2,
                IntRange(0, 2));

        BoolOption symmetry("WBO", "symmetry", "Symmetry breaking.\n", true);

        IntOption symmetry_lim(
                "WBO", "symmetry-limit",
                "Limit on the number of symmetry breaking clauses.\n", 500000,
                IntRange(0, INT32_MAX));

        IntOption cluster_algorithm("Clustering", "ca",
                                    "Clustering algorithm "
                                            "(0=none, 1=DivisiveMaxSeparate)",
                                    1, IntRange(0, 1));
        IntOption num_clusters("Clustering", "c", "Number of agglomerated clusters",
                               100000, IntRange(1, INT_MAX));

        IntOption rounding_strategy(
                "Clustering", "rs",
                "Statistic used to select"
                        " common weights in a cluster (0=Mean, 1=Median, 2=Min)",
                0, IntRange(0, 2));


        IntOption num_conflicts(
                "Incomplete", "conflicts", "Limit on the number of conflicts.\n", 10000,
                IntRange(0, INT32_MAX));

        IntOption num_iterations(
                "Incomplete", "iterations", "Limit on the number of iterations.\n", 100000,
                IntRange(0, INT32_MAX));

        BoolOption local("Incomplete", "local", "Local limit on the number of conflicts.\n", false);

        BoolOption polConservative("TorcOpenWbo", "conservative", "Apply conservative polarity heuristic?\n", true);
        BoolOption conservativeUseAllVars("TorcOpenWbo", "conservative_use_all_vars",
                                          "Re-use the polarity of all the variables within the conservative approach (or, otherwise, only the initial once)?\n",
                                          true);
        BoolOption polOptimistic("TorcOpenWbo", "optimistic", "Set target variables' polarity to the optimum?\n", true);
        IntOption targetVarsBumpVal("TorcOpenWbo", "target_vars_bump_val",
                                    "Bump factor of the activity of the targets at the beginning\n", 113);
        BoolOption targetVarsBumpRelWeights("TorcOpenWbo", "target_vars_bump_rel_weights",
                                            "Bump the variable scores, where the bump value is relative to the weights?\n",
                                            true);

        IntOption targetVarsBumpMaxRandVal("TorcOpenWbo", "target_vars_bump_max_rand_val",
                                           "Maximal random bump factor\n", 552);
        BoolOption optC1("Timetabler", "opt-allocation",
                         "Optimality for Allocation?\n",
                         false);

        parseOptions(argc, argv, true);

        if ((int) num_tests) {
            if ((int) test_join) {
                for (int i = 0; i < (int) num_tests; i++) {
                    test_encoding_join();
                }
            } else {
                for (int i = 0; i < (int) num_tests; i++) {
                    test_encoding();
                }
            }

            return 0;
        }

        Torc::Instance()->SetPolConservative(polConservative);
        Torc::Instance()->SetConservativeAllVars(conservativeUseAllVars);
        Torc::Instance()->SetPolOptimistic(polOptimistic);
        Torc::Instance()->SetTargetVarsBumpVal(targetVarsBumpVal);
        Torc::Instance()->SetBumpRelWeights(targetVarsBumpRelWeights);
        Torc::Instance()->SetTargetBumpMaxRandVal(targetVarsBumpMaxRandVal);

        MaxSAT *S = NULL;

        Statistics rounding_statistic =
                static_cast<Statistics>((int) rounding_strategy);

        switch ((int) algorithm) {
            case _ALGORITHM_WBO_:
                S = new WBO(verbosity, weight, symmetry, symmetry_lim);
                break;

            case _ALGORITHM_LINEAR_SU_:
                if ((int) (cluster_algorithm) == 1) {
                    S = new LinearSUMod(verbosity, bmo, cardinality, pb,
                                        ClusterAlg::_DIVISIVE_, rounding_statistic,
                                        (int) (num_clusters));
                } else {
                    S = new LinearSU(verbosity, bmo, cardinality, pb);
                }
                break;

            case _ALGORITHM_PART_MSU3_:
                S = new PartMSU3(verbosity, partition_strategy, graph_type, cardinality);
                break;

            case _ALGORITHM_MSU3_:
                S = new MSU3(verbosity);
                break;

            case _ALGORITHM_LSU_CLUSTER_:
                S = new LinearSUClustering(verbosity, bmo, cardinality, pb,
                                           ClusterAlg::_DIVISIVE_, rounding_statistic,
                                           (int) (num_clusters));
                break;

            case _ALGORITHM_LSU_MRSBEAVER_:
                S = new OBV(verbosity, cardinality, num_conflicts, num_iterations, local);
                break;

            case _ALGORITHM_LSU_MCS_:
                S = new BLS(verbosity, cardinality, num_conflicts, num_iterations, local);
                break;

            case _ALGORITHM_OLL_:
                if ((int) (cluster_algorithm) == 1) {
                    S = new OLLMod(verbosity, cardinality, ClusterAlg::_DIVISIVE_,
                                   rounding_statistic, (int) (num_clusters));
                } else {
                    S = new OLL(verbosity, cardinality);
                }
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

        MaxSATFormula *maxsat_formula = new MaxSATFormula();
        maxsat_formula->setFormat(_FORMAT_PB_);

        Instance instance= readJSONFile(argv[1]);

        for (int i = 0; i < instance.train.size() ; ++i) {
            for (int j = 0; j < instance.route[instance.train[i].route].totalSeq; ++j) {
                getVariableID("t_"+std::to_string(instance.train[i].id)+"_"+std::to_string(j),maxsat_formula);
            }
        }


        /*compact graph
       for (std::list<Route>::iterator it=instance.route.begin(); it != instance.route.end(); ++it)
            for (std::list<route_path>::iterator it1=it->route_path.begin(); it1 != it1=it->route_path.end(); ++it)
                for (std::list<route_section>::iterator it1=it->route_path.begin(); it1 != it1=it->route_path.end(); ++it)

                }

                }

        }*/


        if (cpu_lim != 0) ;



        S->loadFormula(maxsat_formula);
        printSolverStats(maxsat_formula,initial_time);


        if ((int) (cluster_algorithm) == 1) {
            switch ((int) algorithm) {
                case _ALGORITHM_LINEAR_SU_:
                    static_cast<LinearSUMod *>(S)->initializeCluster();
                    break;
                case _ALGORITHM_OLL_:
                    static_cast<OLLMod *>(S)->initializeCluster();
                    break;
                case _ALGORITHM_LSU_CLUSTER_:
                    static_cast<LinearSUClustering *>(S)->initializeCluster();
                    break;
            }
        }

        S->search();


        outputJSONFile(instance);



    } catch (OutOfMemoryException &) {
        sleep(1);
        printf("c Error: Out of memory.\n");
        printf("s UNKNOWN\n");
        exit(_ERROR_);
    }
}




void printSolverStats(MaxSATFormula*maxsat_formula,double initial_time){
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

}

// Get the variable identifier corresponding to a given name. If the
// variable does not exist, a new identifier is created.
int getVariableID(std::string varName,MaxSATFormula*maxsat_formula) {
    char *cstr = new char[varName.length() + 1];
    strcpy(cstr, varName.c_str());
    int id = maxsat_formula->varID(cstr);
    if (id == var_Undef)
        id = maxsat_formula->newVarName(cstr);
    delete[] cstr;
    return id;
}


void outputJSONFile(Instance instance) {
    StringBuffer s;
    Writer<StringBuffer> writer(s);
    writer.StartObject();               // Between StartObject()/EndObject(),
    writer.Key("problem_instance_label");                // output a key,
    writer.String(instance.label.c_str());             // follow by a value.
    writer.Key("problem_instance_hash");                // output a key,
    writer.Int(instance.hash);             // follow by a value.
    writer.Key("hash");                // output a key,
    writer.Int(42);             // follow by a value.
    writer.Key("train_runs");
    writer.StartArray();
    for(int i=0;i<1;i++){
        writer.Key("service_intention_id");
        writer.Int(2);
        writer.Key("train_run_sections");
        writer.StartArray();
        for(int j=0;j<1;j++){
            writer.Key("entry_time");
            writer.String("");
            writer.Key("exit_time");
            writer.String("");
            writer.Key("route");
            writer.String("");

            writer.Key("route_section_id");
            writer.String("");

            writer.Key("sequence_number");
            writer.String("");

            writer.Key("route_path");
            writer.String("");

            writer.Key("section_requirement");
            writer.String("");




        }
        writer.EndArray();


    }



    writer.EndArray();

    writer.EndObject();







    //Solution to file
    ofstream myfile;
    myfile.open (instance.label+".json");
    myfile << s.GetString();
    myfile.close();




}

Instance readJSONFile(char* local) {

    ifstream ifs(local);
    IStreamWrapper isw(ifs);
    Document d;
    d.ParseStream(isw);

    Instance Instance;

    Instance.hash=d["hash"].GetInt();
    Instance.label=d["label"].GetString();
    std::vector<Train> tt;

    for (int i = 0; i < d["service_intentions"].GetArray().Size(); ++i) {
        Train train;
        train.id=d["service_intentions"].GetArray()[i]["id"].GetInt();
        train.route=d["service_intentions"].GetArray()[i]["route"].GetInt();
        if(train.id!=train.route)
            cerr << "ERR" << endl;
        for (int j = 0; j <d["service_intentions"].GetArray()[i]["section_requirements"].GetArray().Size() ; ++j) {
            list<Requirement> re;
            int id=-1,delay=-1;
            string entry_ea="",exit_earliest="",type="",min_stopping_time="",marker="",exit_latest="",entry_latest="";
            if(d["service_intentions"].GetArray()[i]["section_requirements"].GetArray()[j].HasMember("entry_latest"))
                entry_latest=d["service_intentions"].GetArray()[i]["section_requirements"].GetArray()[j]["entry_latest"].GetString();
            if(d["service_intentions"].GetArray()[i]["section_requirements"].GetArray()[j].HasMember("exit_latest"))
                exit_latest=d["service_intentions"].GetArray()[i]["section_requirements"].GetArray()[j]["exit_latest"].GetString();

            if(d["service_intentions"].GetArray()[i]["section_requirements"].GetArray()[j].HasMember("entry_earliest"))
                entry_ea=d["service_intentions"].GetArray()[i]["section_requirements"].GetArray()[j]["entry_earliest"].GetString();
            if(d["service_intentions"].GetArray()[i]["section_requirements"].GetArray()[j].HasMember("exit_earliest"))
                exit_earliest=d["service_intentions"].GetArray()[i]["section_requirements"].GetArray()[j]["exit_earliest"].GetString();
            if(d["service_intentions"].GetArray()[i]["section_requirements"].GetArray()[j].HasMember("section_marker"))
                marker=d["service_intentions"].GetArray()[i]["section_requirements"].GetArray()[j]["section_marker"].GetString();
            if(d["service_intentions"].GetArray()[i]["section_requirements"].GetArray()[j].HasMember("type"))
                type=d["service_intentions"].GetArray()[i]["section_requirements"].GetArray()[j]["type"].GetString();
            if(d["service_intentions"].GetArray()[i]["section_requirements"].GetArray()[j].HasMember("min_stopping_time"))
                min_stopping_time=d["service_intentions"].GetArray()[i]["section_requirements"].GetArray()[j]["min_stopping_time"].GetString();
            if(d["service_intentions"].GetArray()[i]["section_requirements"].GetArray()[j].HasMember("sequence_number"))
                id=d["service_intentions"].GetArray()[i]["section_requirements"].GetArray()[j]["sequence_number"].GetInt();
            if(d["service_intentions"].GetArray()[i]["section_requirements"].GetArray()[j].HasMember("entry_delay_weight"))
                delay=d["service_intentions"].GetArray()[i]["section_requirements"].GetArray()[j]["entry_delay_weight"].GetInt();

            list<connection> clist;
            if(!d["service_intentions"].GetArray()[i]["section_requirements"].GetArray()[j]["connections"].IsNull())
            for (int k = 0; k < d["service_intentions"].GetArray()[i]["section_requirements"].GetArray()[j]["connections"].GetArray().Size(); ++k) {
                connection c = connection(d["service_intentions"].GetArray()[i]["section_requirements"].GetArray()[j]["connections"].GetArray()[k]["onto_service_intention"].GetInt(),
                                                d["service_intentions"].GetArray()[i]["section_requirements"].GetArray()[j]["connections"].GetArray()[k]["onto_section_marker"].GetString(),
                                                d["service_intentions"].GetArray()[i]["section_requirements"].GetArray()[j]["connections"].GetArray()[k]["min_connection_time"].GetString());
                clist.push_back(c);
            }

            if(id!=-1) {
                Requirement r = Requirement(id,
                                            marker,
                                            type,
                                            min_stopping_time,
                                            entry_ea,
                                            delay,
                                            exit_earliest,entry_latest,exit_latest);
                r.connections = clist;
                //std::cout << r << std::endl;
                //r.toString();
                re.push_back(r);
            }
        }
        tt.push_back(train);
    }
    Instance.train=tt;
    std::map<int,Route> rr;
    std::vector<std::map<std::string,std::vector<std::pair<route_section, route_section>>>> Pairs;

    for (int m = 0; m < d["routes"].GetArray().Size(); ++m) {
        int nSeq=0;
       Route r;
       r.id=d["routes"].GetArray()[m]["id"].GetInt();
       std::list<route_path> rpl;
       route_path rp;std::map<std::string,std::vector<std::pair<route_section, route_section>>> pairs;

        for (int i = 0; i < d["routes"].GetArray()[m]["route_paths"].GetArray().Size(); ++i) {
            if(d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["id"].IsInt())
                rp.id=d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["id"].GetInt()+"";
            else
                rp.id=d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["id"].GetString();
            std::list<route_section> rsl;
            route_section rs, rs1;
            for (int j = 0; j < d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"].GetArray().Size(); j++) {
                nSeq++;
                size++;
                rs.sequence_number = d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["sequence_number"].GetInt();
                std::list<std::string> temp;
                if (d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j].HasMember(
                        "route_alternative_marker_at_entry")) {
                    for (int k = 0; k <
                                    d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["route_alternative_marker_at_entry"].GetArray().Size(); ++k) {
                        temp.push_front(
                                d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["route_alternative_marker_at_entry"].GetArray()[k].GetString());
                    }
                }
                rs.route_alternative_marker_at_entry = temp;
                temp.clear();
                if (d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j].HasMember(
                        "route_alternative_marker_at_exit")) {
                    for (int k = 0; k <
                                    d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["route_alternative_marker_at_exit"].GetArray().Size(); ++k) {
                        temp.push_front(
                                d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["route_alternative_marker_at_exit"].GetArray()[k].GetString());
                    }

                }
                rs.route_alternative_marker_at_exit = temp;
                temp.clear();
                if (d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j].HasMember(
                        "section_marker")) {
                    for (int k = 0; k <
                                    d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["section_marker"].GetArray().Size(); ++k) {
                        temp.push_front(
                                d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["section_marker"].GetArray()[k].GetString());
                    }
                }
                rs.section_marke = temp;
                if (d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j].HasMember(
                        "resource_occupations")) {
                    std::list<Resource> tempR;
                    for (int k = 0; k <
                                    d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["resource_occupations"].GetArray().Size(); ++k) {
                        Resource r;
                        if (d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["resource_occupations"].GetArray()[k]["occupation_direction"].IsString())
                            r = Resource(
                                    d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["resource_occupations"].GetArray()[k]["resource"].GetString(),
                                    d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["resource_occupations"].GetArray()[k]["occupation_direction"].GetString());
                        else
                            r = Resource(
                                    d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["resource_occupations"].GetArray()[k]["resource"].GetString());

                        tempR.push_front(r);
                    }
                    rs.resource_occupations = tempR;
                } else {
                    std::list<Resource> tempR;
                    rs.resource_occupations = tempR;
                }
                if (!d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["penalty"].IsNull())
                    rs.penalty = d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["penalty"].GetDouble();
                else
                    rs.penalty = 0;
                rs.starting_point = d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["starting_point"].GetString();
                rs.minimum_running_time = d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["minimum_running_time"].GetString();
                rs.minimum_running_time = rs.minimum_running_time.substr(2, 2);
                rs.ending_point = d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["ending_point"].GetString();
                if (j != 0) {
                    std::pair<route_section, route_section> p(rs1, rs);
                    std::string name(rp.id + std::to_string(rs.sequence_number));
                    auto it = pairs.find(name);
                    if (it != pairs.end())
                        it->second.push_back(p);
                    else {
                        std::vector<std::pair<route_section, route_section>> v;
                        v.push_back(p);
                        pairs.insert(
                                std::pair<std::string, std::vector<std::pair<route_section, route_section>>>(
                                        rp.id + std::to_string(rs.sequence_number), v));
                    }
                }



                if(i!=0){
                    if(j==0 || j==d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"].GetArray().Size()-1){
                       // std::cout<<rs.sequence_number<<" c "<<*rs.route_alternative_marker_at_exit.begin()<<" "<<rp.route_section.size()<<std::endl;
                        for (auto l= rpl.begin(); l != rpl.end(); ++l) {
                            for (auto k = l->route_section.begin(); k != l->route_section.end(); ++k) {
                                if (k->route_alternative_marker_at_entry.size() != 0 &&
                                    rs.route_alternative_marker_at_exit.size() != 0) {
                                    std::string a = *rs.route_alternative_marker_at_exit.begin();
                                    if (a.compare(*k->route_alternative_marker_at_entry.begin()) == 0) {
                                        std::pair<route_section, route_section> p(rs, *k);
                                        auto it = pairs.find(rp.id + std::to_string(k->sequence_number));
                                        if (it != pairs.end())
                                            it->second.push_back(p);
                                        else {
                                            std::vector<std::pair<route_section, route_section>> v;
                                            v.push_back(p);
                                            pairs.insert(
                                                    std::pair<std::string, std::vector<std::pair<route_section, route_section>>>(
                                                            rp.id + std::to_string(k->sequence_number),
                                                            std::vector<std::pair<route_section, route_section> >(v)));
                                        }                                       // std::cout << k->sequence_number << " " << rs.sequence_number << " A "
                                        //       << *k->route_alternative_marker_at_entry.begin() << " " << a << std::endl;

                                    }
                                }
                                if (k->route_alternative_marker_at_exit.size() != 0 &&
                                    rs.route_alternative_marker_at_entry.size() != 0) {
                                    std::string a = *rs.route_alternative_marker_at_entry.begin();
                                    if (a.compare(*k->route_alternative_marker_at_exit.begin()) == 0) {
//                                            std::cout<<k->sequence_number<<" "<<rs.sequence_number<<" "<<*k->route_alternative_marker_at_entry.begin()<<" "<<*rs.route_alternative_marker_at_exit.begin()<<std::endl;
                                        std::pair<route_section, route_section> p(*k,rs);
                                        auto it = pairs.find(rp.id + std::to_string(rs.sequence_number));
                                        if (it != pairs.end())
                                            it->second.push_back(p);
                                        else {
                                            std::vector<std::pair<route_section, route_section>> v;
                                            v.push_back(p);
                                            pairs.insert(
                                                    std::pair<std::string, std::vector<std::pair<route_section, route_section>>>(
                                                            rp.id + std::to_string(rs.sequence_number),
                                                            std::vector<std::pair<route_section, route_section> >(v)));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                rs1=rs;
                rsl.push_front(rs1);
            }

            rp.route_section=rsl;
            rpl.push_front(rp);

        }
        Pairs.push_back(pairs);
       r.route_path=rpl;
        r.totalSeq=nSeq;
        rr.insert(std::pair<int,Route>(r.id,r));
    }
    Instance.route=rr;

    Instance.pair=Pairs;
    //printPairs(Pairs);

    std::list<Resource> reso;
    for (int l = 0; l < d["resources"].GetArray().Size(); ++l) {
        Resource resource = Resource(d["resources"].GetArray()[l]["id"].GetString(),d["resources"].GetArray()[l]["release_time"].GetString(),d["resources"].GetArray()[l]["following_allowed"].GetBool());
      //  std::cout<<resource<<std::endl;
        reso.push_front(resource);

    }
    Instance.resource=reso;
    Instance.maxBandabweichung=d["parameters"].GetObject()["maxBandabweichung"].GetString();

    return Instance;
}

void printPairs(
        vector<list<pair<route_section, route_section>, allocator<pair<route_section, route_section>>>, allocator<list<pair<route_section, route_section>, allocator<pair<route_section, route_section>>>>> &Pairs) {
    for (int n = 0; n < Pairs.size(); ++n) {
        for (auto i = Pairs.at(n).begin(); i != Pairs.at(n).end() ; ++i) {
            cout << i->first.sequence_number << " " << i->second.sequence_number << endl;
        }

    }
}