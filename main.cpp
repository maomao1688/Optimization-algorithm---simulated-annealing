#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <tuple>
#include <algorithm>
#include <time.h>
#include <stdlib.h>

using namespace std;

#define INPUT "simple1.txt"
#define OUTPUT "output.txt"
#define UPTRY 500 //尝试布局次数上限
#define RANDOMTERM 1 //随机系数
#define RUNTIMES 10000     //优化次数

// 定义Pin类
class Pin {
public:
    std::string name;
    double x;
    double y;
};

// 定义FlipFlop类
class FlipFlop {
public:
    int bits;
    std::string name;
    double width;
    double height;
    int pinCount;
    std::vector<Pin> pins;
};

// 定义Gate类
class Gate {
public:
    std::string name;
    double width;
    double height;
    int pinCount;
    std::vector<Pin> pins;
};

// 定义Instance类
class Instance {
public:
    std::string name;
    std::string libCellName;
    double x;
    double y;
    bool ifFF;          //是否为FF，true为FF，false为Gate
};

// 定义Net类
class Net {
public:
    std::string name;
    int numPins;
    std::vector<std::string> pins;
};

// 定义PlacementRow类
class PlacementRow {
public:
    double startX;
    double startY;
    double siteWidth;
    double siteHeight;
    int totalNumOfSites;
};


// 定义Design类来封装整个设计
class Design {
public:
    double alpha;
    double beta;
    double gamma;
    double dieSize[4];
    std::vector<std::pair<std::string, std::pair<double, double>>> inputPins;
    std::vector<std::pair<std::string, std::pair<double, double>>> outputPins;
    std::vector<FlipFlop> flipFlopslib;
    std::map<string, Gate> gateslib;
    std::vector<Instance> instances;
    std::vector<Net> nets;
    std::vector<PlacementRow> placementRows; //每行网格的绝对大小和个数,目前简单认为均匀
    double displacementDelayCoefficient;
    std::map<std::string, double> qPinDelays;
    std::vector<std::tuple<std::string, std::string, double>> timingSlacks;
    std::map<std::string, double> gatePowers;

    void parseInputFile(const std::string& filename);
    void writeOutputFile(const std::string& filename, const std::vector<Instance>& outputInstances, const std::vector<pair<std::string, std::string>>& pinMappings) const;
};

class location
{
public:
    int type;           //instance的类型，具体储存为在FFlib中的序号
    int x;              //格点矩阵的x
    int y;              //格点矩阵的y
};

void error(string str) {
    std::cout << str << endl;
    system("PAUSE");
    exit(0);
}

// 定义解析输入文件的方法
void Design::parseInputFile(const std::string& filename) {
    std::ifstream inputFile(filename);
    std::string line;
    if (!inputFile.is_open()) {
        error("Can't find the input file!");
    }
    while (std::getline(inputFile, line)) {
        std::istringstream iss(line);
        std::string keyword;
        iss >> keyword;

        if (keyword == "Alpha") {
            iss >> alpha;
        }
        else if (keyword == "Beta") {
            iss >> beta;
        }
        else if (keyword == "Gamma") {
            iss >> gamma;
        }
        else if (keyword == "DieSize") {
            iss >> dieSize[0] >> dieSize[1] >> dieSize[2] >> dieSize[3];
        }
        else if (keyword == "NumInput") {
            int numInput;
            iss >> numInput;
            for (int i = 0; i < numInput; ++i) {
                std::getline(inputFile, line);
                std::istringstream pinIss(line);
                std::string pinName;
                double x, y;
                pinIss >> keyword >> pinName >> x >> y;
                inputPins.emplace_back(pinName, std::make_pair(x, y));
            }
        }
        else if (keyword == "NumOutput") {
            int numOutput;
            iss >> numOutput;
            for (int i = 0; i < numOutput; ++i) {
                std::getline(inputFile, line);
                std::istringstream pinIss(line);
                std::string pinName;
                double x, y;
                pinIss >> keyword >> pinName >> x >> y;
                outputPins.emplace_back(pinName, std::make_pair(x, y));
            }
        }
        else if (keyword == "FlipFlop") {
            FlipFlop ff;
            iss >> ff.bits >> ff.name >> ff.width >> ff.height >> ff.pinCount;
            for (int i = 0; i < ff.pinCount; ++i) {
                std::getline(inputFile, line);
                std::istringstream pinIss(line);
                Pin pin;
                pinIss >> keyword >> pin.name >> pin.x >> pin.y;
                ff.pins.push_back(pin);
            }
            flipFlopslib.push_back(ff);
        }
        else if (keyword == "Gate") {
            Gate gate;
            iss >> gate.name >> gate.width >> gate.height >> gate.pinCount;
            for (int i = 0; i < gate.pinCount; ++i) {
                std::getline(inputFile, line);
                std::istringstream pinIss(line);
                Pin pin;
                pinIss >> keyword >> pin.name >> pin.x >> pin.y;
                gate.pins.push_back(pin);
            }
            gateslib.insert(make_pair(gate.name, gate));
        }
        else if (keyword == "NumInstances") {
            int numInstances;
            iss >> numInstances;
            for (int i = 0; i < numInstances; ++i) {
                std::getline(inputFile, line);
                std::istringstream instIss(line);
                Instance inst;
                instIss >> keyword >> inst.name >> inst.libCellName >> inst.x >> inst.y;
                instances.push_back(inst);
            }
        }
        else if (keyword == "NumNets") {
            int numNets;
            iss >> numNets;
            for (int i = 0; i < numNets; ++i) {
                std::getline(inputFile, line);
                std::istringstream netIss(line);
                Net net;
                netIss >> keyword >> net.name >> net.numPins;
                for (int j = 0; j < net.numPins; ++j) {
                    std::getline(inputFile, line);
                    std::istringstream pinIss(line);
                    std::string pinName;
                    pinIss >> keyword >> pinName;
                    net.pins.push_back(pinName);
                }
                nets.push_back(net);
            }
        }
        else if (keyword == "PlacementRows") {
            PlacementRow row;
            iss >> row.startX >> row.startY >> row.siteWidth >> row.siteHeight >> row.totalNumOfSites;
            placementRows.push_back(row);
        }
        else if (keyword == "DisplacementDelay") {
            iss >> displacementDelayCoefficient;
        }
        else if (keyword == "QpinDelay") {
            std::string libCellName;
            double delay;
            iss >> libCellName >> delay;
            qPinDelays[libCellName] = delay;
        }
        else if (keyword == "TimingSlack") {
            std::string instanceCellName, pinName;
            double slack;
            iss >> instanceCellName >> pinName >> slack;
            timingSlacks.emplace_back(instanceCellName, pinName, slack);
        }
        else if (keyword == "GatePower") {
            std::string libCellName;
            double powerConsumption;
            iss >> libCellName >> powerConsumption;
            gatePowers[libCellName] = powerConsumption;
        }
    }
}

// 定义写入输出文件的方法
void Design::writeOutputFile(const std::string& filename, const std::vector<Instance>& outputInstances, const std::vector<pair<std::string, std::string>>& pinMappings) const {
    std::ofstream outputFile(filename);
    outputFile << "CellInst " << outputInstances.size() << "\n";
    for (const auto& inst : outputInstances) {
        outputFile << "Inst " << inst.name << " " << inst.libCellName << " " << inst.x << " " << inst.y << "\n";
    }
    outputFile << endl;
    for (int i = 0; i < pinMappings.size(); i++) {
        outputFile << pinMappings[i].first << "\tmap\t" << pinMappings[i].second << "\n";
    }
}

void printDesign(const Design& design) {
    std::cout << "Design Parameters:" << std::endl;
    std::cout << "Alpha: " << design.alpha << std::endl;
    std::cout << "Beta: " << design.beta << std::endl;
    std::cout << "Gamma: " << design.gamma << std::endl;

    std::cout << "Die Size: [" << design.dieSize[0] << ", " << design.dieSize[1] << ", " << design.dieSize[2] << ", " << design.dieSize[3] << "]" << std::endl;

    std::cout << "Input Pins:" << std::endl;
    for (const auto& inputPin : design.inputPins) {
        std::cout << "  Name: " << inputPin.first << ", Delay: " << inputPin.second.second << std::endl;
    }

    std::cout << "Output Pins:" << std::endl;
    for (const auto& outputPin : design.outputPins) {
        std::cout << "  Name: " << outputPin.first << ", Delay: " << outputPin.second.second << std::endl;
    }

    std::cout << "Flip Flops:" << std::endl;
    for (const auto& flipFlop : design.flipFlopslib) {
        std::cout << "  Name: " << flipFlop.name << ", Bits: " << flipFlop.bits << ", Width: " << flipFlop.width << ", Height: " << flipFlop.height << ", Pin Count: " << flipFlop.pinCount << std::endl;
        for (const auto& pin : flipFlop.pins) {
            std::cout << "    Pin Name: " << pin.name << ", X: " << pin.x << ", Y: " << pin.y << std::endl;
        }
    }

    std::cout << "Gates:" << std::endl;
    for (const auto& gate : design.gateslib) {
        std::cout << "  Name: " << gate.second.name << ", Width: " << gate.second.width << ", Height: " << gate.second.height << ", Pin Count: " << gate.second.pinCount << std::endl;
        for (const auto& pin : gate.second.pins) {
            std::cout << "    Pin Name: " << pin.name << ", X: " << pin.x << ", Y: " << pin.y << std::endl;
        }
    }

    std::cout << "Instances:" << std::endl;
    for (const auto& instance : design.instances) {
        std::cout << "  Name: " << instance.name << ", LibCellName: " << instance.libCellName << ", X: " << instance.x << ", Y: " << instance.y << std::endl;
    }

    std::cout << "Nets:" << std::endl;
    for (const auto& net : design.nets) {
        std::cout << "  Name: " << net.name << ", Num Pins: " << net.numPins << std::endl;
        for (const auto& pin : net.pins) {
            std::cout << "    Pin: " << pin << std::endl;
        }
    }

    std::cout << "Placement Rows:" << std::endl;
    for (const auto& row : design.placementRows) {
        std::cout << "  StartX: " << row.startX << ", StartY: " << row.startY << ", SiteWidth: " << row.siteWidth << ", SiteHeight: " << row.siteHeight << ", TotalNumOfSites: " << row.totalNumOfSites << std::endl;
    }

    std::cout << "Displacement Delay Coefficient: " << design.displacementDelayCoefficient << std::endl;

    std::cout << "Q Pin Delays:" << std::endl;
    for (const auto& qPinDelay : design.qPinDelays) {
        std::cout << "  Name: " << qPinDelay.first << ", Delay: " << qPinDelay.second << std::endl;
    }

    std::cout << "Timing Slacks:" << std::endl;
    for (const auto& slack : design.timingSlacks) {
        std::cout << "  Start Point: " << std::get<0>(slack) << ", End Point: " << std::get<1>(slack) << ", Slack: " << std::get<2>(slack) << std::endl;
    }

    std::cout << "Gate Powers:" << std::endl;
    for (const auto& gatePower : design.gatePowers) {
        std::cout << "  Gate: " << gatePower.first << ", Power: " << gatePower.second << std::endl;
    }
}

int findinlib(Design* design, Instance instance) {        //查找某个instance所属的FF在libvector中的序号，若不是FF，则返回-1
    for (int i = 0; i < design->flipFlopslib.size(); i++) {
        if (instance.libCellName.compare(design->flipFlopslib[i].name) == 0) return i;
    }
    return -1;
}//查找某个instance所属的FF在libvector中的序号

//根据名字查找instance
Instance* findinstance(string name, Design* design) {
    for (int i = 0; i < design->instances.size(); i++) {
        if (name.compare(design->instances[i].name) == 0) {
            return &design->instances[i];
        }
    }
}

location findlattice(Instance* instance, Design* design) {        //根据位置坐标求格点矩阵xy的函数，目前适用每行高相同且宽相同的情况
    location loc;
    loc.x = (instance->x) / (design->placementRows[0].siteWidth);
    loc.y = (instance->y) / (design->placementRows[0].siteHeight);
    loc.type = findinlib(design, *instance);
    return loc;
}

struct Pos {
    int x;
    int y;
};

class Mysort {
public:
    bool operator()(Pos p1, Pos p2) const {
        if (p1.x < p2.x) {
            return true;
        }
        else if (p1.x > p2.x) {
            return false;
        }
        else {
            if (p1.y < p2.y) {
                return true;
            }
            else {
                return false;
            }
        }
    }
};

//检查是否重叠
bool coverTest(vector<vector<int>>& lattice, int sx, int sy, int w, int h, int width, int height, int mw, int mh);

//矩阵拷贝覆盖
void copyLattice(vector<vector<int>>& lattice, vector<vector<int>>& lattice_save, int mw, int mh);

//生成latticeHeight*latticeWidth全矩阵格点坐标的随机序列
void randomPos_maker(vector<pair<int, int>>& randomPos, int mw, int mh);

//半矩形周长模型计算
int miniRectangle(const vector<pair<int, int>>& points);

int main() {
    srand((unsigned int)time(NULL));
    std::string inputFileName = INPUT;
    std::string outputFileName = OUTPUT;
    Design design;
    design.parseInputFile(inputFileName);

    printDesign(design);

    //对数据结构进行初始处理
    //对FFlib根据寄存器位数排序，由高到低,以便于优化时确定不同FF的数目
    sort(design.flipFlopslib.begin(), design.flipFlopslib.end(), [&](const FlipFlop& a, const FlipFlop& b)->bool {
        return a.bits > b.bits;
        });
    for (int i = 0; i < design.flipFlopslib.size(); i++) {
        cout << design.flipFlopslib[i].name << "\t" << design.flipFlopslib[i].bits << endl;
    }
    //考虑到placerows的起始xy可能不为0，故在所有处理之前将所有instance的xy减去起始xy
    for (int i = 0; i < design.instances.size(); i++) {
        design.instances[i].x -= design.placementRows[0].startX;
        design.instances[i].y -= design.placementRows[0].startY;
    }


    //初始情况
    //一、每种寄存器的采用数目
    vector<int> FFnum;              //每种FF所用数量,其顺序与lib vector的顺序相同
    for (int i = 0; i < design.flipFlopslib.size(); i++) {
        FFnum.push_back(0);
    }//初始化FFnum,全部置0
    for (int i = 0; i < design.instances.size(); i++) {
        design.instances[i].ifFF = false;   //先将所有instance的ifFF初始化为false
        int iffind = findinlib(&design, design.instances[i]);
        if (iffind >= 0) {      //若为FF
            design.instances[i].ifFF = true;
            FFnum[iffind]++;
        }
    }   //取得初始输入电路的FFnum情况

    for (int i = 0; i < FFnum.size(); i++) {
        cout << FFnum[i] << endl;
    }

    int bitsnum = 0;             //定义所有FF的总位数
    for (int i = 0; i < FFnum.size(); i++) {
        bitsnum += FFnum[i] * design.flipFlopslib[i].bits;
    }//求解总位数

    cout << "bitsnum=" << bitsnum << endl;

    //二、每个instance的位置
    int FFinstance_num = 0;             //定义FFinstance_num代表一共有多少个FF实例
    for (int i = 0; i < FFnum.size(); i++) {
        FFinstance_num += FFnum[i];
    }

    vector<location> FFlocation;  //代表每个FF实例的位置，以格点矩阵的x、y表示

    //接下来求输入初始电路的FFlocation
    for (int i = 0; i < design.instances.size(); i++) {
        if (design.instances[i].ifFF) {         //若为FF
            FFlocation.push_back(findlattice(&design.instances[i], &design));
        }
    }
    cout << "初始位置：" << endl;
    for (int i = 0; i < FFlocation.size(); i++) {
        cout << design.flipFlopslib[FFlocation[i].type].name << " " << FFlocation[i].x << "\t" << FFlocation[i].y << endl;
    }

    vector<pair<Pin, Pin>> FFpins;                  //存储初始电路的所有FFinstance的D、Q节点名称以及绝对位置,first为D,Second为Q
    Instance tmpinstance;
    string tmp_instancename, tmp_pinname;
    Pin tmppinD;
    Pin tmppinQ;
    for (int i = 0; i < design.nets.size(); i++) {
        for (int j = 0; j < design.nets[i].pins.size(); j++) {
            if (design.nets[i].pins[j].compare(design.nets[i].pins[j].size() - 1, 1, "D") == 0 || design.nets[i].pins[j].compare(design.nets[i].pins[j].size() - 2, 1, "D") == 0) {    //若某个pin的name的最后一位是D或倒数第二位是D
                tmp_instancename = design.nets[i].pins[j].substr(0, design.nets[i].pins[j].find("/"));  //instance名字
                tmpinstance = *findinstance(tmp_instancename, &design);
                tmp_pinname = design.nets[i].pins[j].substr(design.nets[i].pins[j].find("/") + 1);                       //去掉前缀的instance名字
                for (int k = 0; k < design.flipFlopslib[findinlib(&design, tmpinstance)].pins.size(); k++) {    //在器件库的模板里寻找这个引脚
                    if (design.flipFlopslib[findinlib(&design, tmpinstance)].pins[k].name.compare(tmp_pinname) == 0) {      //前面的name是器件库里面的，是不带instance名字的
                        tmppinD.x = design.flipFlopslib[findinlib(&design, tmpinstance)].pins[k].x + tmpinstance.x;
                        tmppinD.y = design.flipFlopslib[findinlib(&design, tmpinstance)].pins[k].y + tmpinstance.y;
                        //根据pinname找到他所属的instanace以获取器件的绝对坐标，进而找到FFlib以获取引脚相对于器件的相对坐标,二者相加，即为引脚的绝对坐标
                    }
                }
                tmppinD.name = design.nets[i].pins[j];
                //接下来将tmp_pinname中的D改成Q，再处理一遍
                tmp_pinname[0] = 'Q';
                for (int k = 0; k < design.flipFlopslib[findinlib(&design, tmpinstance)].pins.size(); k++) {    //在器件库的模板里寻找这个引脚
                    if (design.flipFlopslib[findinlib(&design, tmpinstance)].pins[k].name.compare(tmp_pinname) == 0) {      //前面的name是器件库里面的，是不带instance名字的
                        tmppinQ.x = design.flipFlopslib[findinlib(&design, tmpinstance)].pins[k].x + tmpinstance.x;
                        tmppinQ.y = design.flipFlopslib[findinlib(&design, tmpinstance)].pins[k].y + tmpinstance.y;
                        //根据pinname找到他所属的instanace以获取器件的绝对坐标，进而找到FFlib以获取引脚相对于器件的相对坐标,二者相加，即为引脚的绝对坐标
                    }
                }
                tmppinQ.name = tmp_instancename + "/" + tmp_pinname;
                FFpins.push_back(pair<Pin, Pin>(tmppinD, tmppinQ));
            }
        }

    }
    
    //求初始电路的maprelation
    vector<vector<int>> mapRelation;        //存储映射关系，存储数字代表FFpins向量的序号
    mapRelation.resize(FFinstance_num);
    int k = 0;
    for (int i = 0; i < FFnum.size(); i++) {
        for (int j = 0; j < FFnum[i]; j++) {
            mapRelation[k].resize(design.flipFlopslib[i].bits);     //二维vector中第一维存储所有FFinstance，第二维的容量为此触发器的位数
            k++;
        }
    }



	//此处计算初始解
	double cost_matric = 0;
	int tmplabel = 0;
    double TNS=0;
    tmplabel = 0;
    for (int i = 0; i < FFnum.size(); i++) {
        for (int j = 0; j < FFnum[i]; j++) {
            FFlocation[tmplabel].type = i;
            cost_matric += design.beta * design.gatePowers.find(design.flipFlopslib[FFlocation[tmplabel].type].name)->second;
            cost_matric += design.gamma * design.flipFlopslib[FFlocation[tmplabel].type].width * design.flipFlopslib[FFlocation[tmplabel].type].height;
        }
    }
    for (auto it = design.timingSlacks.begin(); it != design.timingSlacks.end(); it++) {
        cost_matric += get<2>(*it);
    }
    double cost_initial = cost_matric;
    double costfinal = cost_initial;
    cout << "初始电路的目标函数为："<< cost_initial << endl;


    //优化过程
    //数据结构初始处理
    int T0 = 1000;              //初始温度    
    int T = 1000;               //温度(不断下降)
   
    //格点矩阵
    int width = design.placementRows[0].siteWidth;
    int height = design.placementRows[0].siteHeight;
    int latticeHeight = design.placementRows.size() + 1;
    int latticeWidth = design.placementRows[0].totalNumOfSites + 1;
    vector<vector<int>> lattice(latticeHeight, vector<int>(latticeWidth, 0));       //格点矩阵lattice，初始化为0
    //初始化格点矩阵
    for (vector<Instance>::iterator it = design.instances.begin(); it != design.instances.end(); it++) {
        if (it->ifFF == true) {
            continue;
        }
        else { //是组合逻辑的instance
            int up = it->y + design.gateslib.find(it->libCellName)->second.height;
            int right = it->x + design.gateslib.find(it->libCellName)->second.width;
            for (int i = it->x; i <= right; i += width) {
                for (int j = it->y; j <= up; j += height) {
                    lattice[j / height][i / width] = 1;
                }
            }
        }
    }


    //开始优化
    for (int times = 0; times <= RUNTIMES; times++)
    {
        T=T-T0/RUNTIMES;
    restart://是否需要考虑重置的情况
        //1.确定不同FF数目
        //更新FFnum这个vector
        int remainbits = bitsnum;           //剩余的位数;
        int randchangepos = 0;
        int randchangeneg = 0;
        int maxnum = 0;
        for (int i = 0; i < FFnum.size(); i++) {
            if (i == FFnum.size() - 1) {
                FFnum[i] = remainbits;
            }
            else
            {

                if (remainbits > 0)
                {
                    maxnum = remainbits / design.flipFlopslib[i].bits;
                    if ((maxnum - FFnum[i]) > 0)
                    {
                        randchangepos = rand() % (maxnum - FFnum[i] + 1);
                    }
                    else {
                        if ((maxnum - FFnum[i] == 0)) {
                            randchangepos = 0;
                        }
                        else {
                            randchangepos = 0;
                        }
                    }
                    if (FFnum[i] > 0)
                    {
                        if (maxnum - FFnum[i] >= 0) randchangeneg = rand() % (FFnum[i] + 1);
                        else randchangeneg = FFnum[i] - maxnum + rand() % (maxnum + 1);

                    }
                    else {
                        randchangeneg = 0;
                    }
                    FFnum[i] = FFnum[i] + randchangepos * (T / T0) - randchangeneg * (T / T0);
                }
                else if (remainbits == 0) FFnum[i] = 0;
                else error("remianbits<0");
            }
            remainbits -= (FFnum[i] * design.flipFlopslib[i].bits);
            cout << "remianbits=" << remainbits << "\t" << "FFnum[i]=" << FFnum[i] << endl;

        }

        //更新FFinstance_num
        FFinstance_num = 0;
        for (int i = 0; i < FFnum.size(); i++) {
            FFinstance_num += FFnum[i];
        }
        cout << FFinstance_num << endl;


        //2.确定映射关系
        mapRelation.resize(FFinstance_num);
        k = 0;
        for (int i = 0; i < FFnum.size(); i++) {
            for (int j = 0; j < FFnum[i]; j++) {
                mapRelation[k].resize(design.flipFlopslib[i].bits);     //二维vector中第一维存储所有FFinstance，第二维的容量为此触发器的位数
                k++;
            }
        }

        k = 0;
        for (int i = 0; i < mapRelation.size(); i++) {
            for (int j = 0; j < mapRelation[i].size(); j++) {
                mapRelation[i][j] = k;
                k++;
            }
        }       //初始的映射为按顺序映射，存入的数字为FFpins的序号

        for (int i = 0; i < mapRelation.size(); i++) {
            for (int j = 0; j < mapRelation[i].size(); j++) {
                cout <<"初始的mapRelation:" << mapRelation[i][j] << "\t";
            }
            cout << endl;
        }

        //优化过程：随机抽取打乱
        int randi = 0;
        int randj = 0;
        int tmp = 0;
        for (int i = mapRelation.size() - 1; i >= 0; i--) {
            for (int j = mapRelation[i].size() - 1; j >= 0; j--) {
                randi = rand() % (i + 1);
                randj = rand() % (mapRelation[randi].size());
                tmp = mapRelation[i][j];
                mapRelation[i][j] = mapRelation[randi][randj];
                mapRelation[randi][randj] = tmp;
            }
        }
        for (int i = 0; i < mapRelation.size(); i++) {
            for (int j = 0; j < mapRelation[i].size(); j++) {
                cout << mapRelation[i][j] << "\t";
            }
            cout << endl;
        }
        cout << "接下来确定位置" << endl;

        //3.确定每个FFinstance的位置
        //备份格点矩阵
        vector<vector<int>> lattice_save(latticeHeight, vector<int>(latticeWidth, 0));
        copyLattice(lattice, lattice_save, latticeWidth, latticeHeight);
        //开始布局
        FFlocation.resize(FFinstance_num);
        int tryTime = 0; //尝试布局次数
        int index = 0;
        vector<pair<int, int>> randomPos;
    rePlacement:
        tmplabel = 0;
        copyLattice(lattice_save, lattice, latticeWidth, latticeHeight);
        for (int i = 0; i < FFnum.size(); i++) {
            for (int j = 0; j < FFnum[i]; j++) {
                FFlocation[tmplabel].type = i;
                index = 0;
                randomPos_maker(randomPos, latticeWidth, latticeHeight);
                FFlocation[tmplabel].x = randomPos[index].first;
                FFlocation[tmplabel].y = randomPos[index].second;
                while (coverTest(lattice, FFlocation[tmplabel].x, FFlocation[tmplabel].y, design.flipFlopslib[FFlocation[tmplabel].type].width, design.flipFlopslib[FFlocation[tmplabel].type].height, width, height, latticeWidth, latticeHeight) == false) {
                    index++;
                    FFlocation[tmplabel].x = randomPos[index].first;
                    FFlocation[tmplabel].y = randomPos[index].second;
                    if (index >= latticeWidth * latticeHeight - 1) {
                        tryTime++;
                        cout << "trytime is:" << tryTime << endl;
                        if (tryTime <= UPTRY) {
                            goto rePlacement;
                        }
                        else {
                            cout << "restart" << endl;
                            goto restart;

                        }
                    }
                }
                cout << "布局序号：" << tmplabel << endl;
                tmplabel++;
            }
        }
        cout << "开始计算" << endl;
        //此处计算目标函数
        TNS = 0;
        tmplabel = 0;
        cost_matric = 0;
        for (int i = 0; i < FFnum.size(); i++) {
            for (int j = 0; j < FFnum[i]; j++) {
                FFlocation[tmplabel].type = i;
                cost_matric += design.beta * design.gatePowers.find(design.flipFlopslib[FFlocation[tmplabel].type].name)->second;
                cost_matric += design.gamma * design.flipFlopslib[FFlocation[tmplabel].type].width * design.flipFlopslib[FFlocation[tmplabel].type].height;
                TNS = 0;
                //计算延迟模型
                int DQnum = design.flipFlopslib[FFlocation[tmplabel].type].pinCount / 2;//当前计数器包含的D-Q对数
                for (int k = 0; k <= DQnum - 1; k++) { //遍历每对D-Q
                    string Dname = FFpins[mapRelation[tmplabel][k]].first.name;//D的绝对名称
                    //现在查找连接D的另一节点或更多节点的绝对名称
                    vector<string> link_Dname;//连接D的另一节点或更多节点的绝对名称
                    link_Dname.clear(); //清零
                    for (vector<Net>::iterator nt = design.nets.begin(); nt != design.nets.end(); nt++) { //遍历所有互连线
                        if (find(nt->pins.begin(), nt->pins.end(), Dname) != nt->pins.end()) { //发现D是某一互连线的其中一个节点
                            link_Dname.insert(link_Dname.end(), nt->pins.begin(), nt->pins.end());
                            link_Dname.erase(find(link_Dname.begin(), link_Dname.end(), Dname));
                            break;
                        }
                    }
                    
                    //获取Dname的初始绝对坐标
                    int D_ix, D_iy, drx, dry;
                    int tmp = 0;
                    for (int i = 0; i < FFpins.size(); i++) {
                        if (FFpins[i].first.name.compare(Dname) == 0) {
                            D_ix = FFpins[i].first.x;
                            D_iy = FFpins[i].first.y;
                            tmp = i;
                        }
                    }
                   
                    //获取Dname的当前绝对坐标
                    int D_cx, D_cy;
                    int tmpi = 0;
                    int tmpj = 0;
                    for (int i = 0; i < mapRelation.size(); i++) {
                        for (int j = 0; j < mapRelation[i].size(); j++) {
                            if (mapRelation[i][j] == tmp) {
                                tmpi = i;
                                tmpj = j;
                            }
                        }
                    }
                    drx = design.flipFlopslib[FFlocation[tmpi].type].pins[tmpj].x;
                    dry = design.flipFlopslib[FFlocation[tmpi].type].pins[tmpj].y;
                    D_cx = FFlocation[tmpi].x + drx;
                    D_cy = FFlocation[tmpi].y + dry;


                    
                    //获取link_Dname的身份及其初始绝对坐标
                    vector<pair<int, int>> linkD_pos;  //组合逻辑的link_D的初始绝对坐标，与link_Dname一一对应
                    vector<pair<int, int>> linkD_pos1; //触发器的link_D的初始绝对坐标
                    vector<pair<int, int>> linkD_pos2; //触发器的link_D的当前绝对坐标
                    linkD_pos.clear();
                    linkD_pos1.clear();
                    linkD_pos2.clear();
                    bool isIO = false;
                    for (vector<string>::iterator st = link_Dname.begin(); st != link_Dname.end(); st++) {
                        //先考虑连接到输入、输出端
                        isIO = false;
                        for (vector<pair<string, pair<double, double>>>::iterator pt = design.inputPins.begin(); pt != design.inputPins.end(); pt++) {
                            if (pt->first == *st) { //匹配到输入端
                                linkD_pos.push_back(pt->second);
                                isIO = true;
                                break;
                            }
                        }
                        for (vector<pair<string, pair<double, double>>>::iterator pt = design.outputPins.begin(); pt != design.outputPins.end(); pt++) {
                            if (pt->first == *st) { //匹配到输出端
                                linkD_pos.push_back(pt->second);
                                isIO = true;
                                break;
                            }
                        }
                        if (isIO == false) {
                            //连接到instance的端口
                            for (vector<Instance>::iterator it = design.instances.begin(); it != design.instances.end(); it++) { //遍历所有instance
                                if (it->ifFF == false) { //组合逻辑
                                    Gate tmpGate = design.gateslib.find(it->libCellName)->second;
                                    vector<Pin>::iterator vP;
                                    for (vP = tmpGate.pins.begin(); vP != tmpGate.pins.end(); vP++) {
                                        if (vP->name == *st) {
                                            break;
                                        }
                                    }
                                    if (vP != tmpGate.pins.end()) { //定位到该组合逻辑模块中
                                        linkD_pos.push_back(make_pair(vP->x, vP->y)); //存储该节点的初始绝对坐标
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    
                    double timingStack_old; //当前D初始时所在触发器的timingStack
                    for (vector<tuple<string, string, double>>::iterator it = design.timingSlacks.begin(); it != design.timingSlacks.end(); it++) {
                        if (get<0>(*it) == Dname.substr(0, Dname.find('/')) && get<1>(*it) == Dname.substr(Dname.find('/') + 1, Dname.length())) {
                            timingStack_old = get<2>(*it);
                            break;
                        }
                    }
                    double HPWL_old, HPWL_new;
                    vector<pair<int, int>> transPos;
                    transPos.clear();
                    transPos.insert(transPos.end(), linkD_pos.begin(), linkD_pos.end());    //D连接组合逻辑节点
                    transPos.insert(transPos.end(), linkD_pos1.begin(), linkD_pos1.end());  //D连接初始触发器节点，一般为空
                    transPos.push_back(make_pair(D_ix, D_iy)); //D初始绝对坐标
                    HPWL_old = miniRectangle(transPos);
                    transPos.clear();
                    transPos.insert(transPos.end(), linkD_pos.begin(), linkD_pos.end());    //D连接组合逻辑节点
                    transPos.insert(transPos.end(), linkD_pos2.begin(), linkD_pos2.end());  //D连接当前触发器节点，一般为空
                    transPos.push_back(make_pair(D_cx, D_cy)); //D当前绝对坐标
                    HPWL_new = miniRectangle(transPos);
                    double timingStack_new; //当前D现在的timingStack
                    timingStack_new = timingStack_old + design.displacementDelayCoefficient * (HPWL_old - HPWL_new);
                    TNS += timingStack_new;
                    
                }
                cost_matric += design.alpha * (design.qPinDelays.find(design.flipFlopslib[FFlocation[tmplabel].type].name)->second + TNS);
                tmplabel++;
                
            }
        }

        srand((unsigned int)time(NULL));
        
        cout << "第" << times << "次，" << "此次优化的目标函数值为：" << cost_matric << endl;
        cout << endl;
        if (cost_matric < costfinal) costfinal = cost_matric;
        cout << "目前的最小目标函数值为" << costfinal << endl;
        cout << "初始目标函数值为：" << cost_initial << endl;
        cout << "共运行" << times << "次" << endl;
        
    }

    






    // 示例输出实例和引脚映射 (请根据具体需求修改)
    std::vector<Instance> outputInstances;

    string newInstance_name;
    string newInstance_type;
    Instance newInstance;
    for (int i = 0; i < FFlocation.size(); i++) {
        newInstance_name = "C" + to_string(design.instances.size() + i + 1);
        newInstance_type = design.flipFlopslib[FFlocation[i].type].name;
        newInstance.name = newInstance_name;
        newInstance.libCellName = newInstance_type;
        newInstance.x = design.placementRows[0].startX + design.placementRows[0].siteWidth * FFlocation[i].x;
        newInstance.y = design.placementRows[0].startY + design.placementRows[0].siteHeight * FFlocation[i].y;
        outputInstances.push_back(newInstance);
    }

    vector<pair<std::string, std::string>> pinMappings;

    string newpinname, oldpinname;
    for (int i = 0; i < mapRelation.size(); i++) {
        for (int j = 0; j < mapRelation[i].size(); j++)
        {
            if (mapRelation[i].size() == 1) {
                newpinname = "C" + to_string(design.instances.size() + i + 1) + "/D";
            }
            else {
                newpinname = "C" + to_string(design.instances.size() + i + 1) + "/D" + to_string(j);
            }

            pinMappings.push_back(pair<string, string>(FFpins[mapRelation[i][j]].first.name, newpinname));
            oldpinname = FFpins[mapRelation[i][j]].first.name;
            oldpinname.replace(oldpinname.find('D', 0), 1, "Q");
            newpinname.replace(newpinname.find('D', 0), 1, "Q");
            pinMappings.push_back(pair<string, string>(oldpinname, newpinname));
        }
    }
    cout << "运行时间为:" << clock() / CLOCKS_PER_SEC << "s" << endl;

    design.writeOutputFile(outputFileName, outputInstances, pinMappings);

    return 0;
}

int miniRectangle(const vector<pair<int, int>>& points) {
    if (points.empty()) {
        return 0;
    }

    int min_x = points[0].first;
    int max_x = points[0].first;
    int min_y = points[0].second;
    int max_y = points[0].second;

    for (const auto& point : points) {
        if (point.first < min_x) min_x = point.first;
        if (point.first > max_x) max_x = point.first;
        if (point.second < min_y) min_y = point.second;
        if (point.second > max_y) max_y = point.second;
    }

    int width = max_x - min_x;
    int height = max_y - min_y;

    return width + height;
}

//检查是否重叠
//sx,sy为instance锚点所在格点，w,h为instance大小，width、height为die的每个网格的绝对大小，mw，mh为die的矩阵行列数
bool coverTest(vector<vector<int>>& lattice, int sx, int sy, int w, int h, int width, int height, int mw, int mh) {//sx,sy为instance锚点绝对坐标，w,h为instance大小，width、height为die的每个网格的绝对大小，mw，mh为die的矩阵行列数
    int x, y; //当前待判定坐标
    //首先判断关键点位
    x = sx;
    y = sy;
    if (x >= mw || y >= mh) return false;
    if (lattice[y][x] == 1) {
        return false;
    }
    x = sx + w / width;
    y = sy + h / height;
    if (x >= mw || y >= mh) return false;
    if (lattice[y][x] == 1) {
        return false;
    }
    x = sx / width;
    y = sy + h / height;
    if (x >= mw || y >= mh) return false;
    if (lattice[y][x] == 1) {
        return false;
    }
    x = sx + w / width;
    y = sy / height;
    if (x >= mw || y >= mh) return false;
    if (lattice[y][x] == 1) {
        return false;
    }
    //完全扫描
    for (int i = sx ; i <= sx + w / width ; i++) {
        for (int j = sy; j <= sy + h / height ; j++) {
            if (i >= mw || j >= mh) return false;
            if (lattice[j][i] == 1) {
                return false;
            }
        }
    }
    return true;
}

//矩阵拷贝覆盖
void copyLattice(vector<vector<int>>& lattice, vector<vector<int>>& lattice_save, int mw, int mh) {
    for (int i = 0; i < mw; i++) {
        for (int j = 0; j < mh; j++) {
            lattice_save[j][i] = lattice[j][i];
        }
    }
}

//生成latticeHeight*latticeWidth全矩阵格点坐标的随机序列
void randomPos_maker(vector<pair<int, int>>& randomPos, int mw, int mh) {
    int index = 0, m, n;
    randomPos.clear();
    for (int i = 0; i < mw; i++) {
        for (int j = 0; j < mh; j++) {
            randomPos.push_back(make_pair(i, j));
        }
    }
    index = 1;
    pair<int, int> tmp;
    while (index++ <= mw * mh * RANDOMTERM) {
        //srand((unsigned int)time(NULL));
        m = rand() % mw;
        //srand((unsigned int)time(NULL));
        n = rand() % mh;
        tmp = randomPos[m];
        randomPos[m] = randomPos[n];
        randomPos[n] = tmp;
    }
}