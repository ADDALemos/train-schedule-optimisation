
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
ILOSTLBEGIN




Instance readJSONFile();

void outputJSONFile(Instance instance);

void ILP(const Instance &instance);

using namespace rapidjson;
using namespace std;
const char *name="01_dummy";
typedef IloArray<IloNumVarArray> NumVarMatrix;

int main() {


    Instance instance= readJSONFile();

    //compact graph
    for (int i = 0; i < instance.route.size(); ++i) {

    }
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
    std::string local="../problem_instances/";
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
                std::cout << r << std::endl;
                //r.toString();
                re.push_back(r);
            }
        }
        tt.push_front(train);
    }
    Instance.train=tt;
    std::list<Route> rr;
   for (int m = 0; m < d["routes"].GetArray().Size(); ++m) {
       Route r;
       r.id=d["routes"].GetArray()[m]["id"].GetInt();
       std::list<route_path> rpl;
       route_path rp;
        for (int i = 0; i < d["routes"].GetArray()[m]["route_paths"].GetArray().Size(); ++i) {
            rp.id=d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["id"].GetString();
            std::list<route_section> rsl;
            route_section rs;
            for (int j = 0; j < d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"].GetArray().Size(); ++j) {
                rs.sequence_number=d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["sequence_number"].GetInt();
                if(d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j].HasMember("route_alternative_marker_at_entry")) {
                    std::list<std::string> temp;
                    for (int k = 0; k < d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["route_alternative_marker_at_entry"].GetArray().Size(); ++k) {
                        temp.push_front(d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["route_alternative_marker_at_entry"].GetArray()[k].GetString());
                    }
                    rs.route_alternative_marker_at_entry=temp;
                }
                if(d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j].HasMember("route_alternative_marker_at_exit")) {
                    std::list<std::string> temp;
                    for (int k = 0; k <
                                    d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["route_alternative_marker_at_exit"].GetArray().Size(); ++k) {
                        temp.push_front(
                                d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["route_alternative_marker_at_exit"].GetArray()[k].GetString());
                    }

                    rs.route_alternative_marker_at_exit = temp;
                }
                if(d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j].HasMember("section_marker")){
                    std::list<std::string> temp;
                    for (int k = 0; k < d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["section_marker"].GetArray().Size(); ++k) {
                        temp.push_front(d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["section_marker"].GetArray()[k].GetString());
                    }
                    rs.section_marke=temp;
                }
                if(d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j].HasMember("resource_occupations"))
                {
                    std::list<Resource> temp;
                    for (int k = 0; k < d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["resource_occupations"].GetArray().Size(); ++k) {
                        Resource r;
                        if(d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["resource_occupations"].GetArray()[k]["occupation_direction"].IsString())
                            r= Resource(                        d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["resource_occupations"].GetArray()[k]["resource"].GetString()
                                    ,d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["resource_occupations"].GetArray()[k]["occupation_direction"].GetString());
                        else
                            r=Resource(                        d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["resource_occupations"].GetArray()[k]["resource"].GetString());

                        temp.push_front(r);
                    }
                    rs.resource_occupations=temp;
                }
                if(!d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["penalty"].IsNull())
                    rs.penalty=d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["penalty"].GetDouble();
                rs.starting_point=d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["starting_point"].GetString();
                rs.minimum_running_time=d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["minimum_running_time"].GetString();
                rs.minimum_running_time=rs.minimum_running_time.substr(2,2);
                rs.ending_point=d["routes"].GetArray()[m]["route_paths"].GetArray()[i]["route_sections"][j]["ending_point"].GetString();


            }
            rsl.push_front(rs);
            rp.route_section=rsl;

        }
       rpl.push_front(rp);
       r.route_path=rpl;

       rr.push_front(r);
    }

    
    std::list<Resource> reso;
    for (int l = 0; l < d["resources"].GetArray().Size(); ++l) {
        Resource resource = Resource(d["resources"].GetArray()[l]["id"].GetString(),d["resources"].GetArray()[l]["release_time"].GetString(),d["resources"].GetArray()[l]["following_allowed"].GetBool());
        std::cout<<resource<<std::endl;
        reso.push_front(resource);

    }
    Instance.resource=reso;
    Instance.maxBandabweichung=d["parameters"].GetObject()["maxBandabweichung"].GetString();

    return Instance;
}