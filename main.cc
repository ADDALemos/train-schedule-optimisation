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
#include <limits.h>
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
#include "rapidjson/prettywriter.h"


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
#include "../../../../Library/Developer/CommandLineTools/usr/include/c++/v1/climits"

//Print Solver stats
void printSolverStats(MaxSATFormula*maxsat_formula,double initial_time);

int getVariableID(std::string varName,MaxSATFormula*maxsat_formula);

static void SIGINT_exit(int signum) {
    mxsolver->printAnswer(_UNKNOWN_);
    exit(_UNKNOWN_);
}


int size=-1;
Instance readJSONFile(char *);
Instance readOutputJSONFile(char*);
void outputJSONFile(Instance instance);



using namespace rapidjson;
using namespace std;

int main(int argc, char **argv) {
    readJSONFile(argv[1]);
//    readOutputJSONFile(argv[1]);
    std::exit(1);
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

        /*for (int i = 0; i < instance.train.size() ; ++i) {
            for (int j = 0; j < instance.route[instance.train[i].route].totalSeq; ++j) {
                getVariableID("t^"+instance.train[i].id+"^"+std::to_string(j),maxsat_formula);
            }
        }*/
        std::map<std::string, std::map<int,std::vector<route_section*>>>::iterator it = instance.end.begin();;

        while (it != instance.end.end()) {
            std::map<int,std::vector<route_section*>>::iterator it1 = it->second.begin();
            while (it1 != it->second.end()) {
                if(it1->second[0]->route_alternative_marker_at_entry.size()==0){
                    vec<Lit> lit;
                    lit.push(~mkLit(getVariableID("t^"+it->first+"^"+std::to_string(it1->first),maxsat_formula)));
                    //printf("~%s ",("t^"+it->first+"^"+std::to_string(it1->first)).c_str());
                    for (int i = 1; i < it1->second.size(); ++i) {
                        lit.push(mkLit(getVariableID("t^"+it->first+"^"+std::to_string(it1->second[i]->sequence_number),maxsat_formula)));
                        //printf("%s ",("t^"+it->first+"^"+std::to_string(it1->second[i]->sequence_number)).c_str());

                    }
                    //printf("\n");
                    maxsat_formula->addHardClause(lit);
                    lit.clear();
                }
                it1++;

            }
            it++;

        }
        printf("splits\n");
        std::string delimiter = "^";
        std::map<std::string,std::vector<route_section*>> ::iterator it2 = instance.entryMap.begin();;

        while (it2 != instance.entryMap.end()) {
            for(int y=0; y<it2->second.size();y++) {
                vec <Lit> lit;
                std::string rid = it2->first.substr(it2->first.find(delimiter) + 1, it2->first.size());
                if(instance.exitMap[it2->first].size()>0) {
                    lit.push(~mkLit(getVariableID("t^" + rid + "^" + std::to_string(it2->second[y]->sequence_number),
                                                  maxsat_formula)));
                    // printf("~%s ", ("t^" + rid + "^" + std::to_string(it2->second[y]->sequence_number)).c_str());
                    for (int i = 0; i < instance.exitMap[it2->first].size(); ++i) {
                        lit.push(mkLit(getVariableID(
                                "t^" + rid + "^" + std::to_string(instance.exitMap[it2->first][i]->sequence_number),
                                maxsat_formula)));
                        //printf("%s ", ("t^" + rid + "^" + std::to_string(instance.exitMap[it2->first][i]->sequence_number)).c_str());

                    }
                    //printf("\n");
                    maxsat_formula->addHardClause(lit);
                    lit.clear();
                }
            }
            it2++;



        }

        printf("musts\n");
        for (int j = 0; j < instance.train.size(); ++j) {

            for(Requirement *r: instance.train[j].t){

                vec<Lit> lit;
                    for(int k=0; k<instance.markerMap[instance.train[j].id+"^"+r->section_marker].size();k++){
                        lit.push(mkLit(getVariableID(
                                "t^" + instance.train[j].id + "^" + std::to_string(instance.markerMap[instance.train[j].id+"^"+r->section_marker][k]->sequence_number),maxsat_formula)));
                    }

                    maxsat_formula->addHardClause(lit);
                    lit.clear();

            }

        }
        printf("Opt\n");
        std::map<std::string, double >::iterator itpen = instance.route_pen.begin();;
        PBObjFunction *of = new PBObjFunction();
        while (itpen != instance.route_pen.end()) {
            //vec<Lit> litpen;
            std::string rid = itpen->first.substr(0, itpen->first.find(delimiter));
            std::string section = itpen->first.substr(itpen->first.find(delimiter) + 1, itpen->first.size());
            //litpen.push(mkLit(getVariableID("t^" + rid + "^" + section,maxsat_formula)));

            //printf("%f %s \n",itpen->second,("t^" + rid + "^" + section).c_str());
            of->addProduct(mkLit(getVariableID(
                    "t^" + rid + "^" + section,maxsat_formula)),ceil(itpen->second));
            //maxsat_formula->addSoftClause(100,litpen);
            //litpen.clear();
            itpen++;
        }
        maxsat_formula->addObjFunction(of);

        printf("Time\n");

        for (int j = 0; j < instance.train.size(); ++j) {
            printf("%d\n",j);
            for(Requirement *r: instance.train[j].t){
                printf("ee: %d el: %d xe: %d xl: %d\n",r->sec_entry_earliest,r->sec_entry_latest,
                       r->sec_exit_earliest,r->sec_exit_latest);

            }

        }
        //for all path
        //30000+ d+ d <= 31800
        //30000+ d >= 3600

        /*compact graph
       for (std::list<Route>::iterator it=instance.route.begin(); it != instance.route.end(); ++it)
            for (std::list<route_path>::iterator it1=it->route_path.begin(); it1 != it1=it->route_path.end(); ++it)
                for (std::list<route_section>::iterator it1=it->route_path.begin(); it1 != it1=it->route_path.end(); ++it)

                }

                }

        }*/


        if (cpu_lim != 0)
            ;



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

        if(S->search()) {
            for (int i = 0; i < S->model.size(); i++) {
                indexMap::const_iterator iter = maxsat_formula->getIndexToName().find(i);
                if (iter != maxsat_formula->getIndexToName().end()) {
                    if (S->model[i] != l_False) {
                        std::string id =iter->second.substr(iter->second.find(delimiter) + 1, iter->second.size());
                        std::string sid = id.substr(id.find(delimiter) + 1, id.size());
                        std::string rid = id.substr(0,id.find(delimiter));
                        train_run_sections * trs = new train_run_sections();
                        trs->entry_time="";
                        trs->exit_time="";
                        trs->route=rid;
                        trs->route_section_id=rid+"#"+sid;
                        trs->route_path=instance.sectionMap[rid][std::stoi(sid)]->route_pathName;
                        for (int j = 0; j < instance.train.size(); ++j) {
                            if(instance.train[j].id.compare(rid)!=0)
                                continue;
                            for (Requirement *r: instance.train[j].t) {
                                for(int k=0; k<instance.markerMap[instance.train[j].id+"^"+r->section_marker].size();k++) {
                                      //  printf("%s %s %s %s \n",rid.c_str(),sid.c_str(),std::to_string(instance.markerMap[instance.train[j].id+"^"+r->section_marker][k]->sequence_number).c_str(),r->section_marker.c_str());
                                    if (std::to_string(instance.markerMap[instance.train[j].id+"^"+r->section_marker][k]->sequence_number).compare(sid) == 0) {
                                        trs->section_requirement=r->section_marker;
                                        break;
                                    }
                                }
                            }
                        }
                        if(instance.results.find(rid)!=instance.results.end())
                            instance.results[rid].insert(std::pair<int,train_run_sections*>(std::stoi(sid),trs));
                        else{
                            std::map<int,train_run_sections*> trsv;
                            trsv.insert(std::pair<int,train_run_sections*>(std::stoi(sid),trs));
                            instance.results.insert(std::pair<std::string,std::map<int,train_run_sections*>>(rid,trsv));
                        }


                    }

                }
            }
            outputJSONFile(instance);
        }




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

    PrettyWriter<StringBuffer> writer(s);
    writer.StartObject();               // Between StartObject()/EndObject(),
    writer.Key("problem_instance_label");                // output a key,
    writer.String(instance.label.c_str());             // follow by a value.
    writer.Key("problem_instance_hash");                // output a key,
    writer.Int(instance.hash);             // follow by a value.
    writer.Key("hash");                // output a key,
    writer.Int(42);             // follow by a value.
    writer.Key("train_runs");
    writer.StartArray();
    std::map<std::string,std::map<int,train_run_sections*>> ::iterator it = instance.results.begin();
    while (it != instance.results.end()) {
        writer.StartObject();
        writer.Key("service_intention_id");
        writer.String(it->first.c_str());
        writer.Key("train_run_sections");
        writer.StartArray();
        std::map<int,train_run_sections*>::iterator it1 = it->second.begin();
        int j=1;
        while (it1 != it->second.end()) {
            writer.StartObject();
            writer.Key("entry_time");
            writer.String(it1->second->entry_time.c_str());
            writer.Key("exit_time");
            writer.String(it1->second->exit_time.c_str());
            writer.Key("route");
            writer.String(it1->second->route.c_str());

            writer.Key("route_section_id");
            writer.String(it1->second->route_section_id.c_str());

            writer.Key("sequence_number");
            writer.Int(j);

            writer.Key("route_path");
            writer.String(it1->second->route_path.c_str());

            writer.Key("section_requirement");
            if(it1->second->section_requirement.size()==0)
                writer.Null();
            else
                writer.String(it1->second->section_requirement.c_str());


            writer.EndObject();
            it1++;
            j++;


        }
        it++;
        writer.EndArray();
        writer.EndObject();



    }



    writer.EndArray();

    writer.EndObject();







    //Solution to file

    ofstream myfile;
    myfile.open ("data/"+instance.label+".out.json");
    myfile << s.GetString();
    myfile.close();




}

Instance readOutputJSONFile(char* local) {
    ifstream ifs(local);
    IStreamWrapper isw(ifs);
    Document d;
    d.ParseStream(isw);

    Instance instance;

    instance.hash=d["problem_instance_hash"].GetInt();
    instance.solution_hash=d["hash"].GetInt();
    instance.label=d["problem_instance_label"].GetString();
    std::map<std::string,std::map<int,train_run_sections*>> results;
    int distance=0;
    for (int i = 0; i < d["train_runs"].GetArray().Size(); ++i) {
        std::string service_intention_id;
        if(d["train_runs"].GetArray()[i]["service_intention_id"].IsInt())
            service_intention_id = std::to_string(d["train_runs"].GetArray()[i]["service_intention_id"].GetInt());
        else
            service_intention_id = d["train_runs"].GetArray()[i]["service_intention_id"].GetString();
        std::map<int,train_run_sections*> res;
        int min=INT_MAX;
        int max=0;
        for (int j = 0; j < d["train_runs"].GetArray()[i]["train_run_sections"].GetArray().Size(); ++j) {
            int h1=0,m1=0,s1=0;
            train_run_sections* trs =new train_run_sections();
            trs->entry_time=d["train_runs"].GetArray()[i]["train_run_sections"].GetArray()[j]["entry_time"].GetString();
            sscanf(trs->entry_time.c_str(), "%d:%d:%d", &h1, &m1,&s1);
            if(((h1 * 60*60) + (m1*60)+s1)<min)
                min=((h1 * 60*60) + (m1*60)+s1);
            trs->exit_time=d["train_runs"].GetArray()[i]["train_run_sections"].GetArray()[j]["exit_time"].GetString();
            sscanf(trs->exit_time.c_str(), "%d:%d:%d", &h1, &m1,&s1);
            if(((h1 * 60*60) + (m1*60)+s1)>max)
                max=((h1 * 60*60) + (m1*60)+s1);
            if(d["train_runs"].GetArray()[i]["train_run_sections"].GetArray()[j]["route"].IsInt())
                trs->route=d["train_runs"].GetArray()[i]["train_run_sections"].GetArray()[j]["route"].GetInt();
            else
                trs->route=d["train_runs"].GetArray()[i]["train_run_sections"].GetArray()[j]["route"].GetString();

            if(d["train_runs"].GetArray()[i]["train_run_sections"].GetArray()[j]["route_section_id"].IsString())
                trs->route_section_id=d["train_runs"].GetArray()[i]["train_run_sections"].GetArray()[j]["route_section_id"].GetString();
            else
                trs->route_section_id=std::to_string(d["train_runs"].GetArray()[i]["train_run_sections"].GetArray()[j]["route_section_id"].GetInt());
            trs->route_path=d["train_runs"].GetArray()[i]["train_run_sections"].GetArray()[j]["route_path"].GetString();
            if(d["train_runs"].GetArray()[i]["train_run_sections"].GetArray()[j].HasMember("section_requirement")){
                if(!d["train_runs"].GetArray()[i]["train_run_sections"].GetArray()[j]["section_requirement"].IsNull()){
                    trs->section_requirement=d["train_runs"].GetArray()[i]["train_run_sections"].GetArray()[j]["section_requirement"].GetString();
                }
            }
            if(trs->route_section_id.find("#")!= std::string::npos)
                res.insert(std::pair<int,train_run_sections*>(std::stoi(trs->route_section_id.substr(
                        trs->route_section_id.find("#")+1,trs->route_section_id.size())),trs));
            else
                res.insert(std::pair<int,train_run_sections*>(std::stoi(trs->route_section_id),trs));
        }
        if((max-min)>distance)
            distance=(max-min);
        results.insert(std::pair<std::string,std::map<int,train_run_sections*>>(service_intention_id,res));
    }
    printf("%d\n",distance);
    return instance;



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
    std::map<std::string, double > route_pen;
    for (int i = 0; i < d["service_intentions"].GetArray().Size(); ++i) {
        Train train;
        if(d["service_intentions"].GetArray()[i]["id"].IsInt())
            train.id=std::to_string(d["service_intentions"].GetArray()[i]["id"].GetInt());
        else
            train.id=d["service_intentions"].GetArray()[i]["id"].GetString();

        if(d["service_intentions"].GetArray()[i]["route"].IsInt())
            train.route=std::to_string(d["service_intentions"].GetArray()[i]["route"].GetInt());
        else
            train.route=d["service_intentions"].GetArray()[i]["route"].GetString();
        std::vector<Requirement*> re;

        for (int j = 0; j <d["service_intentions"].GetArray()[i]["section_requirements"].GetArray().Size() ; ++j) {
            string id="",delay="";
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
            if(d["service_intentions"].GetArray()[i]["section_requirements"].GetArray()[j].HasMember("sequence_number")){
               if(d["service_intentions"].GetArray()[i]["section_requirements"].GetArray()[j]["sequence_number"].IsInt())
                   id=std::to_string(d["service_intentions"].GetArray()[i]["section_requirements"].GetArray()[j]["sequence_number"].GetInt());
                else
                   id = d["service_intentions"].GetArray()[i]["section_requirements"].GetArray()[j]["sequence_number"].GetString();


            }

            if(d["service_intentions"].GetArray()[i]["section_requirements"].GetArray()[j].HasMember("entry_delay_weight")) {
                if(d["service_intentions"].GetArray()[i]["section_requirements"].GetArray()[j]["entry_delay_weight"].IsInt())
                    delay=std::to_string(d["service_intentions"].GetArray()[i]["section_requirements"].GetArray()[j]["entry_delay_weight"].GetInt());
                else
                    delay=std::to_string(d["service_intentions"].GetArray()[i]["section_requirements"].GetArray()[j]["entry_delay_weight"].GetFloat());

            }

            list<connection> clist;

            if(d["service_intentions"].GetArray()[i]["section_requirements"].GetArray()[j].HasMember("connections")){
                if(!d["service_intentions"].GetArray()[i]["section_requirements"].GetArray()[j]["connections"].IsNull()){
                    for (int k = 0; k < d["service_intentions"].GetArray()[i]["section_requirements"].GetArray()[j]["connections"].GetArray().Size(); ++k) {
                        if(!d["service_intentions"].GetArray()[i]["section_requirements"].GetArray()[j]["connections"].GetArray()[k].IsNull()) {
                            connection c = connection(
                                    d["service_intentions"].GetArray()[i]["section_requirements"].GetArray()[j]["connections"].GetArray()[k]["onto_service_intention"].GetInt(),
                                    d["service_intentions"].GetArray()[i]["section_requirements"].GetArray()[j]["connections"].GetArray()[k]["onto_section_marker"].GetString(),
                                    d["service_intentions"].GetArray()[i]["section_requirements"].GetArray()[j]["connections"].GetArray()[k]["min_connection_time"].GetString());
                            clist.push_back(c);
                        }

                    }
                }
            }



            if(id.compare("")!=0) {
                Requirement *r = new Requirement(id,
                                            marker,
                                            type,
                                            min_stopping_time,
                                            entry_ea,
                                            delay,
                                            exit_earliest,entry_latest,exit_latest);
                r->connections = clist;
                //printf("Marker!: %s \n",marker.c_str());
                //std::cout <<"now "<< *r << std::endl;
                //r->toString();
                if(re.size()>0){
                    //std::cout <<"old "<< *re[re.size()-1] << std::endl;
                    if(re[re.size()-1]->exit_latest.compare("")==0){
                        if(r->entry_earliest.compare("")!=0){
                            re[re.size()-1]->sec_exit_latest=r->sec_entry_earliest;//+re[re.size()-1]->min_stopping_time;
                            //printf("earl  %d\n",re[re.size()-1]->sec_exit_latest);
                        } else if(r->exit_latest.compare("")!=0){
                            re[re.size()-1]->sec_exit_latest=r->sec_exit_latest;//+re[re.size()-1]->min_stopping_time;
                            //printf("exit %d\n",re[re.size()-1]->sec_exit_latest);
                        } else {
                            re[re.size()-1]->sec_exit_latest=r->sec_exit_earliest;//+re[re.size()-1]->min_stopping_time;
                            //printf("exit %d\n",re[re.size()-1]->sec_exit_earliest);
                        }
                    }
                    if(r->entry_earliest.compare("")==0){
                        if(re[re.size()-1]->exit_latest.compare("")!=0){
                            r->sec_entry_earliest=re[re.size()-1]->sec_exit_latest;//+re[re.size()-1]->min_stopping_time;
                            //printf("exit ow %d\n",r->sec_entry_earliest);
                        } else  if(re[re.size()-1]->sec_entry_earliest!=-1){
                            r->sec_entry_earliest=re[re.size()-1]->sec_entry_earliest;//+re[re.size()-1]->min_stopping_time;
                            //printf("earl ow %d\n",r->sec_entry_earliest);
                        } else {
                            printf("shit\n");
                        }
                    }


                }
                re.push_back(r);

            }

            //if(std::stoi(id)==d["service_intentions"].GetArray()[i]["section_requirements"].GetArray().Size())
              //  printf("exa: %s l: %s\n",exit_earliest.c_str(),exit_latest.c_str());


        }
        train.t=re;


        tt.push_back(train);
    }


    Instance.train=tt;
    std::map<std::string,Route> rr;
    std::map<std::string, std::map<int,std::vector<route_section*>>> end1;
    std::map<std::string,std::vector<route_section*>> entryMap;
    std::map<std::string,std::vector<route_section*>> exitMap;
    std::map<std::string,std::vector<route_section*>> markerMap;

    std::map<std::string,std::map<int,route_section*>> secMap;


    for (int m = 0; m < d["routes"].GetArray().Size(); ++m) {
        int nSeq=0;
       Route r;
        if(d["routes"].GetArray()[m]["id"].IsInt())
            r.id=std::to_string(d["routes"].GetArray()[m]["id"].GetInt());

        else
            r.id=d["routes"].GetArray()[m]["id"].GetString();
        std::list<route_path> rpl;
        route_path rp;

        for (int i = 0; i < d["routes"].GetArray()[m]["route_paths"].GetArray().Size(); ++i) {
            if(d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["id"].IsInt())
                rp.id=std::to_string(d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["id"].GetInt());
            else
                rp.id=d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["id"].GetString();
            std::list<route_section*> rsl;
            route_section *rs= new route_section();
            route_section *rs1= new route_section();


            for (int j = 0; j < d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"].GetArray().Size(); j++) {
                nSeq++;
                size++;
                rs->sequence_number = d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["sequence_number"].GetInt();
                std::list<std::string> temp;
                if (d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j].HasMember(
                        "route_alternative_marker_at_entry")) {
                    for (int k = 0; k <
                                    d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["route_alternative_marker_at_entry"].GetArray().Size(); ++k) {
                        std::string e =d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["route_alternative_marker_at_entry"].GetArray()[k].GetString();
                        temp.push_front(e);
                    std::string c = e +"^"+ r.id;
//                            printf("Entry: %s s %d\n",c.c_str(),rs->sequence_number);
                        if(entryMap.find(c)!=entryMap.end()){
                            entryMap[c].push_back(rs);
                        } else {
                            std::vector<route_section*> rsV;
                            rsV.push_back(rs);
                            entryMap.insert(std::pair<std::string,std::vector<route_section*>>(c,rsV));
                        }
                    }
                }

                rs->route_alternative_marker_at_entry = temp;
                temp.clear();
                if (d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j].HasMember(
                        "route_alternative_marker_at_exit")) {
                    for (int k = 0; k <
                                    d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["route_alternative_marker_at_exit"].GetArray().Size(); ++k) {
                        std::string e =d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["route_alternative_marker_at_exit"].GetArray()[k].GetString();
                        temp.push_front(e);
                        std::string c = e +"^"+ r.id;
//                            printf("Exit: %s s %d\n",c.c_str(),rs->sequence_number);
                        if(exitMap.find(c)!=exitMap.end()){
                            exitMap[c].push_back(rs);
                        } else {
                            std::vector<route_section*> rsV;
                            rsV.push_back(rs);
                            exitMap.insert(std::pair<std::string,std::vector<route_section*>>(c,rsV));
                        }
                    }

                }
                rs->route_alternative_marker_at_exit = temp;
                temp.clear();
                if (d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j].HasMember(
                        "section_marker")) {
                    for (int k = 0; k <
                                    d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["section_marker"].GetArray().Size(); ++k) {
                        std::string e = d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["section_marker"].GetArray()[k].GetString();
                        temp.push_front(e);
                        std::string c = r.id +"^"+e;
//                        printf("Marker: %s s %d\n",c.c_str(),rs->sequence_number);
                        if(markerMap.find(c)!=exitMap.end()){
                            markerMap[c].push_back(rs);
                        } else {
                            std::vector<route_section*> rsV;
                            rsV.push_back(rs);
                            markerMap.insert(std::pair<std::string,std::vector<route_section*>>(c,rsV));
                        }


                    }
                }
                rs->section_marke = temp;
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
                    rs->resource_occupations = tempR;
                } else {
                    std::list<Resource> tempR;
                    rs->resource_occupations = tempR;
                }
                if(d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j].HasMember("penalty")) {
                    if (!d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["penalty"].IsNull())
                        rs->penalty = d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["penalty"].GetDouble();
                    else
                        rs->penalty = 0;
                } else
                    rs->penalty = 0;
                if(rs->penalty != 0)
                    route_pen.insert(std::pair<std::string, double>(r.id+"^"+std::to_string(rs->sequence_number),rs->penalty));
                rs->route_pathName=rp.id;
                rs->starting_point = d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["starting_point"].GetString();
                rs->minimum_running_time = d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["minimum_running_time"].GetString();
                rs->minimum_running_time = rs->minimum_running_time.substr(2, 2);
                rs->ending_point = d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["ending_point"].GetString();
                if (j > 0) {
                    //printf("train: %s origin %d dest %d\n",r.id.c_str(),rs1->sequence_number,rs->sequence_number);
                    auto it = end1.find(r.id);
                    if (it == end1.end()) {
                        std::map<int, std::vector < route_section * >> mapEnd;
                        std::vector<route_section*> t;
                        t.push_back(rs);
                        t.push_back(rs1);
                        mapEnd.insert(std::pair<int, std::vector < route_section * >>(rs->sequence_number,t));
                        end1.insert(std::pair < std::string, std::map < int,
                                    std::vector < route_section * >> > (r.id, mapEnd));
                    } else {
                        auto it1 = it->second.find(rs->sequence_number);
                        if (it1 == it->second.end()) {
                            std::vector<route_section*> t;
                            t.push_back(rs);
                            t.push_back(rs1);
                            it->second.insert(std::pair<int, std::vector < route_section * >>(rs->sequence_number,t));
                        } else {
                            it1->second.push_back(rs1);
                        }
                    }
                }

                rs1=rs;

                if(secMap.find(r.id)!=secMap.end()) {
                    if (secMap[r.id].find(rs->sequence_number) == secMap[r.id].end()) {
                        secMap[r.id].insert(std::pair<int,route_section*>(rs->sequence_number,rs));
                    } else{
                        printf("OPS: This should not happen: line 959\n");
                        std::exit(1);
                    }

                } else {
                    std::map<int,route_section*> mapT;
                    mapT.insert(std::pair<int,route_section*>(rs->sequence_number,rs));
                    secMap.insert(std::pair<std::string,std::map<int,route_section*>>(r.id,mapT));
                }
                rsl.push_front(rs1);
                rs= new route_section();


            }

            rp.route_section=rsl;
            rpl.push_front(rp);

        }
        r.route_path=rpl;
        r.totalSeq=nSeq;
        rr.insert(std::pair<std::string,Route>(r.id,r));
    }
    Instance.route=rr;
    Instance.exitMap=exitMap;
    Instance.entryMap=entryMap;
    Instance.markerMap=markerMap;
    Instance.route_pen=route_pen;
    Instance.sectionMap=secMap;
    Instance.end=end1;

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



