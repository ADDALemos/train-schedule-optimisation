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

Instance readJSONFile();

void outputJSONFile(Instance instance);

using namespace rapidjson;
using namespace std;
const char *name="02_a_little_less_dummy";

int main() {
    Instance instance= readJSONFile();
    outputJSONFile(instance);


    return 0;
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


    for (int i = 0; i < d["service_intentions"].GetArray().Size(); ++i) {
        Train train;
        train.id=d["service_intentions"].GetArray()[i]["id"].GetInt();
        train.route=d["service_intentions"].GetArray()[i]["route"].GetInt();
        if(train.id!=train.route)
            cerr << "ERR" << endl;
        for (int j = 0; j <d["service_intentions"].GetArray()[i]["section_requirements"].GetArray().Size() ; ++j) {
            list<Requirement> re;
            int id=-1,delay=-1;
            string entry_ea="",exit_earliest="",type="",min_stopping_time="",marker="";
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
                                            exit_earliest);
                r.connections = clist;
               // std::cout << r << std::endl;
                //r.toString();
                re.push_back(r);
            }
        }
    }

    return Instance;
}