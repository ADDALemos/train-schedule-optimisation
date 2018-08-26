
#include "include/rapidjson/reader.h"
#include <iostream>
#include <fstream>
#include "include/rapidjson/document.h"
#include "include/rapidjson/istreamwrapper.h"
#include "include/rapidjson/stringbuffer.h"
#include "include/rapidjson/writer.h"

#include <fstream>
#include "Instance.h"
#include "Train.h"
#include "Resource.h"
#include "Route.h"
#include "route_path.h"
#include "route_section.h"

#ifndef IL_STD
#define IL_STD
#endif
#include <cstring>
#include <ilcplex/ilocplex.h>
#include <vector>

ILOSTLBEGIN




Instance readJSONFile();

void outputJSONFile(Instance instance);

void ILP(const Instance &instance);

using namespace rapidjson;
using namespace std;
const char *name="sample_scenario";
typedef IloArray<IloNumVarArray> NumVarMatrix;

int main() {


    Instance instance= readJSONFile();

    /*compact graph
   for (std::list<Route>::iterator it=instance.route.begin(); it != instance.route.end(); ++it)
        for (std::list<route_path>::iterator it1=it->route_path.begin(); it1 != it1=it->route_path.end(); ++it)
            for (std::list<route_section>::iterator it1=it->route_path.begin(); it1 != it1=it->route_path.end(); ++it)

            }

            }

    }*/
    //ILP(instance);


    outputJSONFile(instance);


    return 0;
}

void ILP(const Instance &instance) {
    IloEnv env;
    IloModel model(env);
    NumVarMatrix accept(env,instance.train.size());
    for(int i=0; i< instance.train.size(); i++) {
        accept[i] = IloNumVarArray(env, instance.route.size());
        for(int j=0; j< instance.route.size(); j++) {
            accept[i][j] = IloNumVar(env, 0.0, 1.0);

        }
    }

    IloCplex cplex(model);
}

void outputJSONFile(Instance instance) {
    StringBuffer s;
    Writer<StringBuffer> writer(s);
    writer.StartObject();               // Between StartObject()/EndObject(),
    writer.Key("problem_instance_label");                // output a key,
    writer.String(name);             // follow by a value.
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
    myfile.open ("../example.txt");
    myfile << s.GetString();
    myfile.close();




}

Instance readJSONFile() {
    std::string local="/Volumes/MAC/train-schedule-optimisation-challenge-starter-kit-master/sample_files/";//"../problem_instances/";
    local.append(name);
    local.append(".json");
    ifstream ifs(local);
    IStreamWrapper isw(ifs);
    Document d;
    d.ParseStream(isw);

    Instance Instance;

    Instance.hash=d["hash"].GetInt();
    Instance.label=d["label"].GetString();
    std::list<Train> tt;

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
              //  std::cout << r << std::endl;
                //r.toString();
                re.push_back(r);
            }
        }
        tt.push_front(train);
    }
    Instance.train=tt;
    std::list<Route> rr;
    std::vector<std::list<std::pair<route_section, route_section>>> Pairs;

    for (int m = 0; m < d["routes"].GetArray().Size(); ++m) {
       Route r;
       r.id=d["routes"].GetArray()[m]["id"].GetInt();
       std::list<route_path> rpl;
       route_path rp;std::list<std::pair<route_section, route_section>> pairs;

        for (int i = 0; i < d["routes"].GetArray()[m]["route_paths"].GetArray().Size(); ++i) {
            if(d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["id"].IsInt())
                rp.id=d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["id"].GetInt()+"";
            else
                rp.id=d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["id"].GetString();
            std::list<route_section> rsl;
            route_section rs, rs1;
            for (int j = 0; j < d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"].GetArray().Size(); j++) {
                rs.sequence_number=d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["sequence_number"].GetInt();
                std::list<std::string> temp;
                if(d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j].HasMember("route_alternative_marker_at_entry")) {
                    for (int k = 0; k < d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["route_alternative_marker_at_entry"].GetArray().Size(); ++k) {
                        temp.push_front(d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["route_alternative_marker_at_entry"].GetArray()[k].GetString());
                    }
                }
                rs.route_alternative_marker_at_entry=temp;
                temp.clear();
                if(d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j].HasMember("route_alternative_marker_at_exit")) {
                    for (int k = 0; k <
                                    d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["route_alternative_marker_at_exit"].GetArray().Size(); ++k) {
                        temp.push_front(
                                d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["route_alternative_marker_at_exit"].GetArray()[k].GetString());
                    }

                }
                rs.route_alternative_marker_at_exit = temp;
                temp.clear();
                if(d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j].HasMember("section_marker")){
                    for (int k = 0; k < d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["section_marker"].GetArray().Size(); ++k) {
                        temp.push_front(d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["section_marker"].GetArray()[k].GetString());
                    }
                }
                rs.section_marke=temp;
                if(d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j].HasMember("resource_occupations"))
                {
                    std::list<Resource> tempR;
                    for (int k = 0; k < d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["resource_occupations"].GetArray().Size(); ++k) {
                        Resource r;
                        if(d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["resource_occupations"].GetArray()[k]["occupation_direction"].IsString())
                            r= Resource(                        d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["resource_occupations"].GetArray()[k]["resource"].GetString()
                                    ,d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["resource_occupations"].GetArray()[k]["occupation_direction"].GetString());
                        else
                            r=Resource(                        d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["resource_occupations"].GetArray()[k]["resource"].GetString());

                        tempR.push_front(r);
                    }
                    rs.resource_occupations=tempR;
                } else{
                    std::list<Resource> tempR;
                    rs.resource_occupations = tempR;
                }
                if(!d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["penalty"].IsNull())
                    rs.penalty=d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["penalty"].GetDouble();
                else
                    rs.penalty=0;
                rs.starting_point=d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["starting_point"].GetString();
                rs.minimum_running_time=d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["minimum_running_time"].GetString();
                rs.minimum_running_time=rs.minimum_running_time.substr(2,2);
                rs.ending_point=d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["ending_point"].GetString();
                if(j!=0) {
                    std::pair<route_section, route_section> p(rs1, rs);
                    pairs.push_front(p);
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
                                        pairs.push_front(p);
                                       // std::cout << k->sequence_number << " " << rs.sequence_number << " A "
                                         //       << *k->route_alternative_marker_at_entry.begin() << " " << a << std::endl;

                                    }
                                }
                                if (k->route_alternative_marker_at_exit.size() != 0 &&
                                    rs.route_alternative_marker_at_entry.size() != 0) {
                                    std::string a = *rs.route_alternative_marker_at_entry.begin();
                                    if (a.compare(*k->route_alternative_marker_at_exit.begin()) == 0) {
//                                            std::cout<<k->sequence_number<<" "<<rs.sequence_number<<" "<<*k->route_alternative_marker_at_entry.begin()<<" "<<*rs.route_alternative_marker_at_exit.begin()<<std::endl;
                                        std::pair<route_section, route_section> p(*k,rs);
                                        pairs.push_front(p);
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

       rr.push_front(r);
    }
    for (int n = 0; n < 1; ++n) {
        for (auto i = Pairs.at(n).begin(); i != Pairs.at(n).end() ; ++i) {
            std::cout<<i->first.sequence_number<<" "<<i->second.sequence_number<<std::endl;
        }

    }

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