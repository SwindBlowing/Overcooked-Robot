#include <enum.h>
#include <framework.h>
#include <decision.h>
#include <string>
#include <string.h>
#include <iostream>
#include <sstream>
#include <string>
#include <cassert>
#include <vector>
#include <string>
#include <queue>
#include <cmath>
#include <sstream>
#define inf 0x3f3f3f3f

const double eps = 1e-7;
const double sht = sqrt(2) / 2;
const double pdx[4 + 5] = {1, 1, -1, -1};
const double pdy[4 + 5] = {-1, 1, 1, -1};
const double dx[8 + 5] = {0, sht, 1, sht, 0, -sht, -1, -sht};
const double dy[8 + 5] = {-1, -sht, 0, sht, 1, sht, 0, -sht};
int nowframe = 0;
std::string directions[8 + 5] = {"U", "RU", "R", "RD", "D", "LD", "L", "LU"};

/* 按照读入顺序定义 */
int width, height;
char Map[20 + 5][20 + 5];
int IngredientCount;
struct Ingredient Ingredient[20 + 5];
int recipeCount;
struct Recipe Recipe[20 + 5];
int totalTime, randomizeSeed, totalOrderCount;
struct Order totalOrder[20 + 5];
int orderCount;
struct Order Order[20 + 5];
int k;
struct Player Players[2 + 5];
int entityCount;
struct Entity Entity[20 + 5];
int remainFrame, Fund;
int plateNum;
int maxWaitingTime = 0;
int emptyNum[2 + 5];
int mustAvoid = -1;

std::pair<double, double> panhome, pothome, chopplace, submitplace, dirtyplateshome, washingplace, deadpoint;// to be init
std::vector <std::pair<int, int>> concretePoints;

int findWork(int x);

bool checkValidPoint(int x, int y)
{
    return x >= 0 && x < width && y >= 0 && y < height;
}

double calDist(std::pair<double, double> a, std::pair<double, double> b)
{
    double ans = 0;
    ans += (a.first - b.first) * (a.first - b.first);
    ans += (a.second - b.second) * (a.second - b.second);
    return sqrt(ans);
}

bool reachTarget(std::pair<double, double> a, std::pair<double, double> b)
{
    return calDist(a, b) < 0.1;
}

bool nearCenter(std::pair<double, double> p1, std::pair<double, double> p2, double d)
{
    return fabs(p1.first - p2.first) < 0.35 && fabs(p1.second - p2.second) < d;
}

double calVecMul(std::pair<double, double> p1, std::pair<double, double> p2)
{
    return p1.first * p2.first + p1.second * p2.second;
}

struct Floyd
{
    double dist[20 + 5][20 + 5][20 + 5][20 + 5];
    
    void initFloyd()
    {
        for (int i1 = 0; i1 < width; i1++)
            for (int j1 = 0; j1 < height; j1++)
                for (int i2 = 0; i2 < width; i2++)
                    for (int j2 = 0; j2 < height; j2++) 
                        if (i1 == i2 && j1 == j2) dist[i1][j1][i2][j2] = 0;
                        else dist[i1][j1][i2][j2] = inf;
        for (int i1 = 0; i1 < width; i1++)
            for (int j1 = 0; j1 < height; j1++) {
                for (int k = 0; k < 8; k += 2) {
                    int i2 = i1 + dx[k], j2 = j1 + dy[k];
                    if (checkValidPoint(i2, j2) && Map[i2][j2] == '.' && Map[i1][j1] == '.') {
                        double addDist = calDist({(double)i1 + 0.5, (double)j1 + 0.5}, {deadpoint.first + 0.5, deadpoint.second + 0.5});
                        dist[i1][j1][i2][j2] = 1 + 1 / (fabs(addDist - 1) + 0.001);
                    }
                }
                /*for (int k = 0; k < 4; k++) {
                    int i2 = i1 + pdx[k], j2 = j1 + pdy[k];
                    if (checkValidPoint(i2, j2) && Map[i2][j2] == '.' && Map[i1][j1] == '.') {
                        double addDist = calDist({(double)i1 + 0.5, (double)j1 + 0.5}, {deadpoint.first + 0.5, deadpoint.second + 0.5});
                        dist[i1][j1][i2][j2] = sht * 2 + 1 / (fabs(addDist - 1) + 0.001);
                    }
                }*/
            }
        for (int k1 = 0; k1 < width; k1++)
            for (int k2 = 0; k2 < height; k2++)
                for (int i1 = 0; i1 < width; i1++)
                    for (int j1 = 0; j1 < height; j1++)
                        for (int i2 = 0; i2 < width; i2++)
                            for (int j2 = 0; j2 < height; j2++)
                                dist[i1][j1][i2][j2] = std::min(dist[i1][j1][i2][j2], dist[i1][j1][k1][k2] + dist[k1][k2][i2][j2]);

    }
}myFloyd;


struct Action
{
    int flag; //0 move, 1 interact, 2 putOrPick, 3 wait
    std::string EntityName;
    Action *pre = NULL, *next = NULL;
};

void InsertBefore(Action *nex, Action *now) {
    Action *pre = nex->pre;
    now->pre = pre;
    now->next = nex;
    pre->next = now;
    nex->pre = now;
}

double calFloyd(std::pair<int, int> a, std::pair<int, int> b)
{
    return myFloyd.dist[a.first][a.second][b.first][b.second];                                                                //TBD
}

struct Decision_High
{
    struct Order currentOrder;
    int currentOrderMaking;
    int recipeNum;
    struct Recipe recipeList[20 + 5];
    Action* ActionListHead;
    Action* ActionListTail;
    struct Entity fixedPlate;

    void ActionInit()
    {
        ActionListHead = new(Action);
        ActionListTail = new(Action);
        ActionListHead->next = ActionListTail;
        ActionListTail->pre = ActionListHead;
    }

    bool checkend(bool *flag)
    {
        if (ActionListHead->next == ActionListTail) {
            recipeNum--;
            *flag = 1;
            if (recipeNum == 0) {
                currentOrderMaking++;
                if (currentOrderMaking == currentOrder.recipe.size())
                    return 1;
                else setRecipe();
            }
            else setActionList();
        }
        return 0;
    }

    void setActionList()
    {
        ActionInit();
        if (recipeList[recipeNum].kind == "-ingre->") {
            Action *nowAction1 = new(Action);
            nowAction1->flag = 0; 
            nowAction1->EntityName = recipeList[recipeNum].nameBefore;
            InsertBefore(ActionListTail, nowAction1);

            Action *nowAction2 = new(Action);
            nowAction2->flag = 2; 
            nowAction2->EntityName = recipeList[recipeNum].nameBefore;
            InsertBefore(ActionListTail, nowAction2);
        }
        else if (recipeList[recipeNum].kind == "-chop->" || recipeList[recipeNum].kind == "-pan->" || recipeList[recipeNum].kind == "-pot->") {
            std::string st;
            if (recipeList[recipeNum].kind == "-chop->") st = "chop";
            if (recipeList[recipeNum].kind == "-pan->") st = "pan";
            if (recipeList[recipeNum].kind == "-pot->") st = "pot";

            Action *nowAction1 = new(Action);
            nowAction1->flag = 0; 
            nowAction1->EntityName = st;
            InsertBefore(ActionListTail, nowAction1);

            Action *nowAction2 = new(Action);
            nowAction2->flag = 2; 
            nowAction2->EntityName = st;
            InsertBefore(ActionListTail, nowAction2);

            if (st == "chop") {
                Action *nowAction3 = new(Action);
                nowAction3->flag = 1; 
                nowAction3->EntityName = st;
                InsertBefore(ActionListTail, nowAction3);
            }

            Action *nowAction4 = new(Action);
            nowAction4->flag = 3; 
            nowAction4->EntityName = st;
            InsertBefore(ActionListTail, nowAction4);

            Action *nowAction5 = new(Action);
            nowAction5->flag = 2; 
            nowAction5->EntityName = st;
            InsertBefore(ActionListTail, nowAction5);
        }
        else if (recipeList[recipeNum].kind == "-plate->") {
            Action *nowAction1 = new(Action);
            nowAction1->flag = 0; 
            nowAction1->EntityName = "plate";
            InsertBefore(ActionListTail, nowAction1);

            Action *nowAction2 = new(Action);
            nowAction2->flag = 2; 
            nowAction2->EntityName = "plate";
            InsertBefore(ActionListTail, nowAction2);
        }
        else if (recipeList[recipeNum].kind == "-submit->") {
            Action *nowAction1 = new(Action);
            nowAction1->flag = 2; 
            nowAction1->EntityName = "plate";
            InsertBefore(ActionListTail, nowAction1);

            Action *nowAction2 = new(Action);
            nowAction2->flag = 0; 
            nowAction2->EntityName = "submit";
            InsertBefore(ActionListTail, nowAction2);

            Action *nowAction3 = new(Action);
            nowAction3->flag = 2; 
            nowAction3->EntityName = "submit";
            InsertBefore(ActionListTail, nowAction3);
        }
        else assert(0);
    }

    void setRecipe()
    {
        //std::cerr << "Hello1!" << std::endl;
        recipeNum = 0;
        std::string str = currentOrder.recipe[currentOrderMaking];
        std::cerr << "My name is " << str << std::endl;
        if (str == "submit") {
            recipeList[++recipeNum].nameBefore = str;
            recipeList[recipeNum].kind = "-submit->";
            setActionList();
            return ;
        }
        recipeList[++recipeNum].nameBefore = str;
        recipeList[recipeNum].kind = "-plate->";
        //std::cerr << "Hello2!" << std::endl;
        bool breakflag = 0;
        while (1) {
            std::cerr << "My name is " << str << std::endl;
            for (int i = 0; i < IngredientCount; i++)
                if (Ingredient[i].name == str) {
                    recipeList[++recipeNum].nameBefore = str;
                    recipeList[recipeNum].kind = "-ingre->";
                    breakflag = 1;
                    break;
                }
            if (breakflag)
                break;
            for (int i = 0; i < recipeCount; i++)
                if (Recipe[i].nameAfter == str) {
                    str = Recipe[i].nameBefore;
                    recipeList[++recipeNum] = Recipe[i];
                }
        }
        setActionList();
    }

    void setWork(struct Order currOrder)
    {
        currentOrder = currOrder;
        currentOrderMaking = 0;
        currentOrder.recipe.push_back("submit");
        setRecipe();
        // set a special "final order" here
    }

};


bool fixedcp[3 + 5];// 0 chop, 1 pan, 2 pot
std::vector <std::pair<double, double>> fixedPlates[2];
std::vector <std::pair<double, double>> fixedDirtyPlates[2];
bool isInteracting[2 + 5];
bool isWashingDishes[2 + 5];
bool lastFlag[2 + 5];

struct Decision_Low
{
    int node;
    struct FullPlayer p;
    Decision_High fd;
    std::pair<double, double> actionTarget;
    std::pair<double, double> plateTarget;
    // action Entity name: 
    // "ingres", chop, pan, pot, plate, dirtyplate
    bool recipechanged;
    bool orderchanged;
    // add a self-actionlist here
    Action *selfListHead = NULL, *selfListTail = NULL;
    bool canDetect = 0;
    bool selfCanDetect = 0;
    ContainerKind containerKind;
    std::vector<std::string> myentity;
    int startWait = 3;
    Action *storedAction;

    void init(int node)
    {
        this->node = node;
        selfListHead = new(Action);
        selfListTail = new(Action);
        selfListHead->next = selfListTail;
        selfListTail->pre = selfListHead;
        canDetect = 0;
        selfCanDetect = 0;
        p.isfree = 1;
        recipechanged = 0;
        orderchanged = 0;
        containerKind = ContainerKind::None;
        myentity.clear();
        for (int i = 0; i < 3; i++)
            if (p.fixedcp[i]) fixedcp[i] = 0, p.fixedcp[i] = 0;
        fixedPlates[node].clear();
        isInteracting[node] = 0;
        isWashingDishes[node] = 0;
    }

    bool findplace(struct FullPlayer nowp, std::string st, std::pair<double, double> *placeTarget)
    {
        if (node == 0) {
            std::cerr << "findplace: " << st << std::endl;
        }
        // find a entity named st in map and fix it. if failed, just ingore the frame.
        for (int i = 0; i < IngredientCount; i++)
            if (Ingredient[i].name == st) {
                *placeTarget = std::make_pair(Ingredient[i].x, Ingredient[i].y);
                if (node == 0) {
                    std::cerr << "ingre: " << placeTarget->first << "," << placeTarget->second << std::endl;
                }
                return 1;
            }
        if (st == "chop") {
            if (fixedcp[0] && !p.fixedcp[0]) return 0;
            *placeTarget = chopplace;
            return 1;
        }
        else if (st == "pan") {
            if (fixedcp[1] && !p.fixedcp[1]) return 0;
            for (int i = 0; i < entityCount; i++)
                if (Entity[i].containerKind == ContainerKind::Pan) {
                    *placeTarget = std::make_pair(Entity[i].x, Entity[i].y);
                    return 1;
                }
            if (Players[node].containerKind == ContainerKind::Pan) 
                return 1;
            return 0;
        }
        else if (st == "pot") {
            if (fixedcp[2] && !p.fixedcp[2]) return 0;
            for (int i = 0; i < entityCount; i++)
                if (Entity[i].containerKind == ContainerKind::Pot) {
                    *placeTarget = std::make_pair(Entity[i].x, Entity[i].y);
                    return 1;
                }
            if (node == 0)
                std::cerr << "containerKind: " << (int)Players[node].containerKind << std::endl;
            if (Players[node].containerKind == ContainerKind::Pot) 
                return 1;
            return 0;
        }
        else if (st == "dirtyplatehome") {
            *placeTarget = dirtyplateshome;
            return 1;
        }
        else if (st == "plate") {
            std::vector <std::pair<double, double>> invalidPos; invalidPos.clear();
            int invalidNum[20 + 5][20 + 5] = {0};
            /*for (int i = 0; i < 2; i++)
                for (int j = 0; j < fixedPlates[i].size(); j++)
                    invalidPos.push_back(fixedPlates[i][j]);*/
            for (int i = 0; i < 2; i++)
                for (int j = 0; j < fixedPlates[i].size(); j++)
                    invalidNum[(int)fixedPlates[i][j].first][(int)fixedPlates[i][j].second]++;
            for (int i = 0; i < 2; i++)
                if (Players[i].containerKind == ContainerKind::Plate) 
                    invalidPos.push_back(std::make_pair(Players[i].x, Players[i].y));
            for (int i = 0; i < entityCount; i++)
                if (Entity[i].containerKind == ContainerKind::Plate) {
                    bool flag = true;
                    for (int j = 0; j < invalidPos.size(); j++)
                        if (fabs(Entity[i].x - invalidPos[j].first) < eps && fabs(Entity[i].y - invalidPos[j].second) < eps) {
                            flag = false;
                            break;
                        }
                    if (flag) {
                        if (node == 0) {
                            std::cerr << "invalid Plates " << (int)Entity[i].x << "," << (int)Entity[i].y
                                << " " << invalidNum[(int)Entity[i].x][(int)Entity[i].y];
                        }
                        if (invalidNum[(int)Entity[i].x][(int)Entity[i].y]) {
                            invalidNum[(int)Entity[i].x][(int)Entity[i].y]--;
                            flag = false;
                        }
                    }
                    if (flag) {
                        *placeTarget = std::make_pair(Entity[i].x, Entity[i].y);
                        return 1;
                    }
                }
            return 0;
        }
        else if (st == "pothome") {
            *placeTarget = pothome;
            return 1;
        }
        else if (st == "panhome") {
            *placeTarget = panhome;
            return 1;
        }
        else if (st == "submit") {
            *placeTarget = submitplace;
            if (node == 0) {
                    std::cerr << "submit: " << placeTarget->first << "," << placeTarget->second << std::endl;
                }
            return 1;
        }
        else if (st == "wash") {
            *placeTarget = washingplace;
            return 1;
        }
        else if (st == "dead") {
            *placeTarget = deadpoint;
            return 1;
        }
        else assert(0);
    }

    std::vector <std::pair<double, double>> findCandidates(std::pair<double, double> placeTarget)
    {
        std::vector <std::pair<double, double>> candidates;
        if (Map[(int)placeTarget.first][(int)placeTarget.second] == '.') 
            candidates.push_back({placeTarget.first, placeTarget.second});
        else {
            for (int i = 0; i < 8; i += 2) {
                int nowx = placeTarget.first + dx[i], nowy = placeTarget.second + dy[i];
                if (checkValidPoint(nowx, nowy) && Map[nowx][nowy] == '.')
                    candidates.push_back({nowx, nowy});
            }
        }
        return candidates;
    }

    std::string playerMove(struct FullPlayer nowp, std::pair<double, double> placeTarget, bool emergency)
    {
        // move to a point in map.
        // if placetarget can not be reached, then the destination is it's around cells.
        double choiceVal[8 + 5] = {0};

        if (emergency) {
            // avoid player crash
            double ndx = Players[node ^ 1].x - nowp.frameState.x, ndy = Players[node ^ 1].y - nowp.frameState.y;
            double dist = calDist({p.frameState.x, p.frameState.y}, {Players[node ^ 1].x, Players[node ^ 1].y});
            if (dist < 2.5) {
                double addNum[8 + 5] = {0};
                for (int i = 0; i < 8; i++)
                    addNum[i] = calVecMul({ndx, ndy}, {dx[i], dy[i]});
                for (int i = 0; i < 8; i++)
                    if (addNum[i] > 0)
                        choiceVal[i] += addNum[i] * 1 / (fabs(dist - 1) + eps);
            }

            // invalid directions
            for (int i = 0; i < 8; i++) {
                int nowx = nowp.frameState.x + dx[i], nowy = nowp.frameState.y + dy[i];
                if (!checkValidPoint(nowx, nowy) || Map[nowx][nowy] != '.') choiceVal[i] = 1000000;
            }

            int minn = 0;
            for (int i = 0; i < 8; i++)
                if (choiceVal[i] < choiceVal[minn]) minn = i;
            
            return directions[minn];
        }

        if (node == 0)
            std::cerr << "placeTarget: " << placeTarget.first << "," << placeTarget.second << std::endl;
        std::vector <std::pair<double, double>> candidates = findCandidates(placeTarget);
        if (node == 0) {
            for (int i = 0; i < candidates.size(); i++)
                std::cerr << "candidates" << i << ": " << candidates[i].first << "," << candidates[i].second << std::endl;
        }
        assert(candidates.size());

        for (int i = 0; i < candidates.size(); i++)
            if (calDist({p.frameState.x, p.frameState.y}, {candidates[i].first + 0.5, candidates[i].second + 0.5}) < 1) {
                double ndx = candidates[i].first + 0.5 - nowp.frameState.x, ndy = candidates[i].second + 0.5 - nowp.frameState.y;
                
                for (int i = 0; i < 8; i++)
                    choiceVal[i] = calVecMul({ndx, ndy}, {dx[i], dy[i]});

                int maxn = 0;
                for (int i = 0; i < 8; i++)
                    if (choiceVal[i] > choiceVal[maxn]) maxn = i;

                
                
                return directions[maxn];
            }
    
        for (int i = 0; i < 8; i++) choiceVal[i] = inf; // the smaller, the better
        for (int i = 0; i < 8; i++) {
            int nowx = nowp.frameState.x + dx[i], nowy = nowp.frameState.y + dy[i];
            /*if (node == 1) {
                std::cerr << "place " << nowp.frameState.x << "," << nowp.frameState.y << std::endl;
                std::cerr << nowx << "," << nowy << std::endl;
            }*/
            for (int j = 0; j < candidates.size(); j++) {
                double anotherDist = calDist({nowp.frameState.x, nowp.frameState.y}, {(double)nowx + 0.5, (double)nowy + 0.5});
                double fdist = calFloyd(candidates[j], {nowx, nowy});
                if (fdist + anotherDist < choiceVal[i])
                    choiceVal[i] = fdist + anotherDist;
                if (node == 0) 
                    std::cerr << anotherDist << "," << fdist << std::endl;
                }
        }

        double ndx, ndy, dist;

        // avoid player crash
        if (mustAvoid == node) {
            ndx = Players[node ^ 1].x - nowp.frameState.x, ndy = Players[node ^ 1].y - nowp.frameState.y;
            dist = calDist({p.frameState.x, p.frameState.y}, {Players[node ^ 1].x, Players[node ^ 1].y});
            if (dist < 2) {
                double addNum[8 + 5] = {0};
                for (int i = 0; i < 8; i++)
                    addNum[i] = calVecMul({ndx, ndy}, {dx[i], dy[i]});
                for (int i = 0; i < 8; i++)
                    if (addNum[i] > 0)
                        choiceVal[i] += addNum[i] * 5 / (fabs(dist - 1) + eps);
            }
        }

        // avoid concretePoints crash
        for (int i = 0; i < 8; i++) {
            double nowx = nowp.frameState.x + 0.2 * dx[i], nowy = nowp.frameState.y + 0.2 * dy[i];
            for (int j = 0; j < concretePoints.size(); j++) {
                double dist = calDist({nowx, nowy}, {(double)concretePoints[j].first + 0.5, (double)concretePoints[j].second + 0.5});
                if (dist < 2)
                    choiceVal[i] += (double)1 / std::pow(dist, 2);
            }
        }

        // avoid deadPoint
        bool willDead = 0;

        for (int j = 15; j; j--) {
            double nowx = nowp.frameState.x + 0.1 * j * p.frameState.X_Velocity, nowy = nowp.frameState.y + 0.1 * j * p.frameState.Y_Velocity;
            if (sameCell({nowx, nowy}, {deadpoint.first + 0.5, deadpoint.second + 0.5})) {
                willDead = 1;
                break;
            }
        }
        if (willDead) {
            dist = calDist({p.frameState.x, p.frameState.y}, {deadpoint.first + 0.5, deadpoint.second + 0.5});
            double v = sqrt(p.frameState.X_Velocity * p.frameState.X_Velocity + p.frameState.Y_Velocity * p.frameState.Y_Velocity);
            if (dist < 3 && v > dist * 2) return "";
        }

        for (int i = 0; i < 8; i++) {
            for (int j = 15; j; j--) {
                double nowx = nowp.frameState.x + 0.1 * j * dx[i], nowy = nowp.frameState.y + 0.1 * j * dy[i];
                if (sameCell({nowx, nowy}, {deadpoint.first + 0.5, deadpoint.second + 0.5}))
                    choiceVal[i] = j * 100000000;
            }
        }
        /*ndy = deadpoint.first + 0.5 - nowp.frameState.x, ndx = deadpoint.second + 0.5 - nowp.frameState.y;
        dist = calDist({p.frameState.x, p.frameState.y}, {deadpoint.first + 0.5, deadpoint.second + 0.5});
        if (dist < 3) {
            /*if (dist < 3) {
                double v = sqrt(p.frameState.X_Velocity * p.frameState.X_Velocity + p.frameState.Y_Velocity * p.frameState.Y_Velocity);
                if (v > 2) return "";
            }
            else *//*if (dist < 2) {
                double addNum[8 + 5] = {0};
                if (ndy > 0) addNum[4] += 3 * fabs(ndy) / fabs(dist - 1);
                else addNum[0] += 3 * fabs(ndy) / fabs(dist - 1);
                if (ndx > 0) addNum[2] += 3 * fabs(ndx) / fabs(dist - 1);
                else addNum[6] += 3 * fabs(ndx) / fabs(dist - 1);
                for (int i = 1; i < 8; i += 2)
                    addNum[i] = addNum[(i - 1 + 8) % 8] * sht + addNum[(i + 1) % 8] * sht;
                for (int i = 0; i < 8; i++)
                    choiceVal[i] += addNum[i];
            }
        }*/

        int minn = 0;
        for (int i = 0; i < 8; i++)
            if (choiceVal[i] < choiceVal[minn]) minn = i;

        /*if (node == 1) {
            for (int i = 0; i < 8; i++)
                std::cerr << "dist" << i << ": " << choiceVal[i] << std::endl;
        }*/
        
        return directions[minn];
    }

    std::string relativePos(struct FullPlayer nowp, std::pair<double, double> placeTarget)
    {
        placeTarget.first += 0.5; placeTarget.second += 0.5;
        if (node == 0) {
            std::cerr << "relativePos" << std::endl;
            std::cerr << nowp.frameState.x << "," << nowp.frameState.y << std::endl;
            std::cerr << placeTarget.first << "," << placeTarget.second << std::endl;
        }
        double dx = placeTarget.first - nowp.frameState.x; // right
        double dy = placeTarget.second - nowp.frameState.y; // down
        if (node == 0) {
            std::cerr << dx << "," << dy << std::endl;
        }
        if (fabs(dy) >= fabs(dx) && dy >= 0) return (std::string)"D";
        if (fabs(dy) >= fabs(dx) && dy < 0) return (std::string)"U";
        if (fabs(dy) < fabs(dx) && dx >= 0) return (std::string)"R";
        if (fabs(dy) < fabs(dx) && dx < 0) return (std::string)"L";
        assert(0);
    }

    bool refixActionTarget(Action *nowAction)
    {
        if (node == 0)
            std::cerr << nowAction->flag << "," << nowAction->EntityName << std::endl;
        if (nowAction->EntityName == "plate") {
            actionTarget = plateTarget;
            return true;
        }
        bool flag = findplace(p, nowAction->EntityName, &actionTarget);
        if (!flag) {
            if (node == 0) {
                std::cerr << "find failed!" << std::endl;
            }
            return 0;
        }
        if (nowAction->EntityName == "chop") fixedcp[0] = 1, p.fixedcp[0] = 1;
        else if (nowAction->EntityName == "pan") fixedcp[1] = 1, p.fixedcp[1] = 1;
        else if (nowAction->EntityName == "pot") fixedcp[2] = 1, p.fixedcp[2] = 1;
        return true;
    }

    bool refixPlate()
    {
        if (!findplace(p, "plate", &plateTarget)) {
            // ignore this now, the procedure is washing dirtyplates.
            // may need to reattribute plates.
            return false;
        }
        fixedPlates[node].push_back(plateTarget);
        return true;
    }

    bool actionFinished(Action *nowAction)
    {
        if (nowAction->flag == 0) {
            std::vector <std::pair<double, double>> candidates = findCandidates(actionTarget);
            int flag = 0;
            for (int i = 0; i < candidates.size(); i++) {
                double nowx = candidates[i].first + 0.5, nowy = candidates[i].second + 0.5;
                flag |= nearCenter(std::make_pair(p.frameState.x, p.frameState.y), {nowx, nowy}, 0.4);
            }
            return flag & (calDist({p.frameState.x, p.frameState.y}, {(double)actionTarget.first + 0.5, (double)actionTarget.second + 0.5}) < 1.3);
        }
        else if (nowAction->flag == 1) {
            for (int i = 0; i < entityCount; i++) {
                if (sameCell(actionTarget, std::make_pair(Entity[i].x, Entity[i].y))) {
                    if (Entity[i].totalFrame) return 1;
                    else return 0;
                }
            }
            return 0;
        }
        else if (nowAction->flag == 2) {
            bool flag = 0;
            if ((int)containerKind != (int)Players[node].containerKind) flag = 1;
            if (myentity.size() != Players[node].entity.size()) flag = 1;
            containerKind = Players[node].containerKind;
            myentity = Players[node].entity;
            return flag;
        }
        else if (nowAction->flag == 3) {
            if (startWait) {
                startWait--;
                return 0;
            }
            if (sameCell(actionTarget, washingplace)) {
                for (int i = 0; i < entityCount; i++) {
                    if (node == 0) {
                        std::cerr << "entities: " << (int)Entity[i].containerKind << " " << Entity[i].x
                            << " " << Entity[i].y << std::endl;
                    }
                    if (Entity[i].containerKind == ContainerKind::DirtyPlates 
                        && sameCell(washingplace, {Entity[i].x, Entity[i].y}))
                        return 0;
                }
                    
                isInteracting[node] = 0;
                return 1;
            }
            else {
                for (int i = 0; i < entityCount; i++) {
                    if (sameCell(actionTarget, std::make_pair(Entity[i].x, Entity[i].y))) {
                        if (Entity[i].totalFrame && Entity[i].currentFrame < Entity[i].totalFrame) 
                            return 0;
                        else {
                            isInteracting[node] = 0;
                            return 1;
                        }
                    }
                }
            }
            std::cerr << "assertion at line " << __LINE__ << " failed!" << std::endl;
            assert(0);
        }
        else assert(0);
        return 0;
    }

    bool adjacent(std::pair<double, double> a, std::pair<double, double> b)
    {
        if (calDist(a, b) > 1.35) return 0;
        for (int i = 0; i < 8; i += 2) {
            double nowx = b.first + dx[i], nowy = b.second + dy[i];
            if (checkValidPoint(nowx, nowy) && Map[(int)nowx][(int)nowy] == '.' && nearCenter(a, {nowx, nowy}, 0.4)) return true;
        }
        return false;
    }

    bool actionStep(Action *nowAction, std::string *frame_action)
    {
        storedAction = nowAction;
        frame_action->clear();
        if (node == 0) {
            std::cerr << "In action step, with task: " << nowAction->flag << "," << nowAction->EntityName << std::endl;
        }
        if (nowAction->flag == 0) {
            frame_action->append("Move");
            frame_action->append(" ");
            frame_action->append(playerMove(p, actionTarget, 0));
        }
        else if (nowAction->flag == 1) {
            if (!adjacent({p.frameState.x, p.frameState.y}, {actionTarget.first + 0.5, actionTarget.second + 0.5})) {
                Action *nowAction1 = new(Action);
                nowAction1->flag = 0;
                nowAction1->EntityName = nowAction->EntityName;
                InsertBefore(selfListHead->next, nowAction1);

                return false;
            }
            frame_action->append("Interact");
            frame_action->append(" ");
            frame_action->append(relativePos(p, actionTarget));
        }
        else if (nowAction->flag == 2) {
            if (!adjacent({p.frameState.x, p.frameState.y}, {actionTarget.first + 0.5, actionTarget.second + 0.5})) {
                if (node == 0)
                    std::cerr << "error: not adjacent!" << std::endl;
                Action *nowAction1 = new(Action);
                nowAction1->flag = 0;
                nowAction1->EntityName = nowAction->EntityName;
                InsertBefore(selfListHead->next, nowAction1);

                return false;
            }
            frame_action->append("PutOrPick");
            frame_action->append(" ");
            frame_action->append(relativePos(p, actionTarget));
        }
        else {
            if (!isInteracting[node]) {
                isInteracting[node] = 1;
                startWait = 3;
            }
            if (!adjacent({p.frameState.x, p.frameState.y}, {actionTarget.first + 0.5, actionTarget.second + 0.5})) {
                if (nowAction->EntityName != "wash" && nowAction->EntityName != "chop") ;
                else {
                    if (node == 0)
                        std::cerr << "error: not adjacent!" << std::endl;
                    Action *nowAction1 = new(Action);
                    nowAction1->flag = 0;
                    nowAction1->EntityName = nowAction->EntityName;
                    InsertBefore(selfListHead->next, nowAction1);

                    Action *nowAction2 = new(Action);
                    nowAction2->flag = 1;
                    nowAction2->EntityName = nowAction->EntityName;
                    InsertBefore(selfListHead->next->next, nowAction2);

                    return false;
                }
                
            }
            frame_action->append("Move ");
        }
        return true;
    }

    bool debuger(std::string *frame_action, std::string addit)
    {
        *frame_action = "Move ";
        *frame_action += addit;
        return 1;
    }

    void initSelfList()
    {
        selfListHead->next = selfListTail;
        selfListTail->pre = selfListHead;
    }

    void checkPlateInhand()
    {
        if (p.frameState.containerKind == ContainerKind::Plate) fixedPlates[node].clear();
    }

    bool makeDecision(std::string *frame_action)
    {
        p.frameFlush(Players[node]);
        storedAction = NULL;

        checkPlateInhand();

        //if (node == 1) return false;

        if (p.frameState.live) {
            init(node);
            for (int i = 0; i < 3; i++)
                if (p.fixedcp[i]) {
                    p.fixedcp[i] = 0;
                    fixedcp[i] = 0;
                }
            return false;
        }

        //if (node == 1) return debuger(frame_action, "U");
        
        
        //do the self list first
        if (node == 0)
            std::cerr << "check selfList" << std::endl;
        if (selfListHead->next != selfListTail) {
            Action *nowAction = selfListHead->next;
            if (selfCanDetect && actionFinished(nowAction)) {
                selfListHead->next = selfListHead->next->next;
                selfListHead->next->pre = selfListHead;

                if (selfListHead->next == selfListTail) 
                    return 0;
            }
            nowAction = selfListHead->next;
            if (!refixActionTarget(nowAction)) return 0;
            selfCanDetect = 1;
            return actionStep(nowAction, frame_action);
        }
        if (node == 0)
            std::cerr << "check empty!" << std::endl;

        isWashingDishes[node] = 0;

        if (p.isfree == 1) {
            // free the plate here.
            fixedPlates[node].clear();

            if (!refixPlate()) {
                assert(selfListHead->next == selfListTail);
                std::cerr << node << " plate not found!" << std::endl;
                if (!isWashingDishes[node ^ 1]) {
                    int dirtyPlateNum = 0;
                    for (int i = 0; i < entityCount; i++)
                        if (Entity[i].containerKind == ContainerKind::DirtyPlates) {
                            if (Entity[i].sum == 0) dirtyPlateNum++;
                            else dirtyPlateNum += Entity[i].sum;
                        }

                    //if (dirtyPlateNum < ceil((double)plateNum / 2)) return 0;
                    if (dirtyPlateNum == 0) return 0;


                    isWashingDishes[node] = 1;

                    Action *nowAction1 = new(Action);
                    nowAction1->flag = 0;
                    nowAction1->EntityName = "dirtyplatehome";
                    InsertBefore(selfListTail, nowAction1);

                    Action *nowAction2 = new(Action);
                    nowAction2->flag = 2;
                    nowAction2->EntityName = "dirtyplatehome";
                    InsertBefore(selfListTail, nowAction2);

                    Action *nowAction3 = new(Action);
                    nowAction3->flag = 0;
                    nowAction3->EntityName = "wash";
                    InsertBefore(selfListTail, nowAction3);

                    Action *nowAction4 = new(Action);
                    nowAction4->flag = 2;
                    nowAction4->EntityName = "wash";
                    InsertBefore(selfListTail, nowAction4);

                    Action *nowAction5 = new(Action);
                    nowAction5->flag = 1;
                    nowAction5->EntityName = "wash";
                    InsertBefore(selfListTail, nowAction5);

                    Action *nowAction6 = new(Action);
                    nowAction6->flag = 3;
                    nowAction6->EntityName = "wash";
                    InsertBefore(selfListTail, nowAction6);

                }
                /*else {
                    Action *nowAction1 = new(Action);
                    nowAction1->flag = 0;
                    nowAction1->EntityName = "wash";
                    InsertBefore(selfListTail, nowAction1);
                }*/
                return 0;
            }

            p.isfree = 0;
            fd.setWork(Order[findWork(node)]); // To Be Done Here
            orderchanged = 1;
            if (node == 0) {
                std::cerr << "reset order: ";
                for (int i = 0; i < fd.currentOrder.recipe.size(); i++)
                    std::cerr << fd.currentOrder.recipe[i] << " ";
                std::cerr << std::endl;
            }
        }

        if (orderchanged) {
            orderchanged = 0;
            recipechanged = 1;
        }

        if (recipechanged) {
            // if sometime an empty container is in hand, put it back.
            assert(selfListHead->next == selfListTail);
            if (selfListHead->next == selfListTail) {
                if (Players[node].containerKind == ContainerKind::Pan && Players[node].entity.empty()) {
                    initSelfList();
                    std::cerr << "I'm here!" << std::endl;
                    Action *nowAction1 = new(Action);
                    nowAction1->flag = 0;
                    nowAction1->EntityName = "panhome";
                    InsertBefore(selfListTail, nowAction1);

                    Action *nowAction2 = new(Action);
                    nowAction2->flag = 2;
                    nowAction2->EntityName = "panhome";
                    InsertBefore(selfListTail, nowAction2);

                    Action *nowAction3 = new(Action);
                    nowAction3->flag = 0;
                    nowAction3->EntityName = "plate";
                    InsertBefore(selfListTail, nowAction3);

                    return 0;
                }
                else if (Players[node].containerKind == ContainerKind::Pot && Players[node].entity.empty()) {
                    initSelfList();
                    std::cerr << "I'm here!" << std::endl;
                    Action *nowAction1 = new(Action);
                    nowAction1->flag = 0;
                    nowAction1->EntityName = "pothome";
                    InsertBefore(selfListTail, nowAction1);

                    Action *nowAction2 = new(Action);
                    nowAction2->flag = 2;
                    nowAction2->EntityName = "pothome";
                    InsertBefore(selfListTail, nowAction2);

                    Action *nowAction3 = new(Action);
                    nowAction3->flag = 0;
                    nowAction3->EntityName = "plate";
                    InsertBefore(selfListTail, nowAction3);

                    return 0;
                }
            }
            // free all fixed entities of p
            for (int i = 0; i < 3; i++)
                if (p.fixedcp[i]) {
                    p.fixedcp[i] = 0;
                    fixedcp[i] = 0;
                }
            recipechanged = 0;
        }

        Action *nowAction = fd.ActionListHead->next;
        if (!refixActionTarget(nowAction)) return 0;
        if (canDetect && actionFinished(nowAction)) {
            canDetect = 0;
            fd.ActionListHead->next = fd.ActionListHead->next->next;
            fd.ActionListHead->next->pre = fd.ActionListHead;
        }

        if (fd.checkend(&recipechanged)) {
            p.isfree = 1;
            return 0;
        }

        if (recipechanged) return 0;


        if (fd.ActionListHead->next == fd.ActionListTail) {
            std::cerr << "assertion at line " << __LINE__ << "failed!" << std::endl;
            assert(0);
        }
        nowAction = fd.ActionListHead->next;
        if (node == 0) {
            std::cerr << "refixActionTarget" << std::endl;
        }
        if (!refixActionTarget(nowAction)) return 0;
        if (node == 0)
            std::cerr << "actionStep" << std::endl;
        canDetect = 1;
        return actionStep(nowAction, frame_action);
    }

}DecisionSystem[2 + 5];

int preAct[2 + 5];
int clearTime[2 + 5];

bool makeDecision(int x, std::string *frameAction)
{
    std::cerr << "now emptyNums: " << emptyNum[0] << "," << emptyNum[1] << std::endl;
    emptyNum[0] = emptyNum[1] = 0;
    mustAvoid = -1;
    for (int i = 0; i < width; i++)
        for (int j = 0; j < height; j++)
            if (Map[i][j] == '.') {
                double dist0 = calDist({Players[0].x, Players[0].y}, {(double)i + 0.5, (double)j + 0.5});
                double dist1 = calDist({Players[1].x, Players[1].y}, {(double)i + 0.5, (double)j + 0.5});
                if (dist0 < dist1) emptyNum[0]++;
                else emptyNum[1]++;
            }
    mustAvoid =  (nowframe / 100) % 2;
    /*if (isInteracting[1]) mustAvoid = 0;
    if (!lastFlag[0]) mustAvoid = 0;
    if (emptyNum[0] > emptyNum[1] * 20) mustAvoid = 0;*/
    
    bool flag = DecisionSystem[x].makeDecision(frameAction);
    lastFlag[x] = flag;
    
    /*if (flag) {
        assert(DecisionSystem[x].storedAction != NULL);
        if (DecisionSystem[x].storedAction->flag != preAct[x]) clearTime[x] = 0;
        preAct[x] = DecisionSystem[x].storedAction->flag;
        if (DecisionSystem[x].storedAction->flag == 1 || DecisionSystem[x].storedAction->flag == 2 || DecisionSystem[x].storedAction->flag == 3) {
            clearTime[x]++;
            if (DecisionSystem[x].storedAction->flag != 3 && clearTime[x] >= 50) {
                clearTime[x] = 0;
                preAct[x] = 0;
                DecisionSystem[x].p.isfree = 1;
            }
            else if (DecisionSystem[x].storedAction->flag == 3 && clearTime[x] >= maxWaitingTime + 10) {
                clearTime[x] = 0;
                preAct[x] = 0;
                DecisionSystem[x].p.isfree = 1;
            }
        }
        if (*frameAction == "Move") clearTime[x]++;
        else {
            std::string preSt = frameAction->substr(1 + frameAction->find(' '));
            if (preSt == "PutOrPick" || preSt == "Interact") clearTime[x]++;
            else clearTime[x] = 0;
        }
        if (clearTime[x] >= 500) {
            DecisionSystem[x].init(x);

            Action *nowAction1 = new(Action);
            nowAction1->flag = 0;
            nowAction1->EntityName = "dead";
            InsertBefore(DecisionSystem[x].selfListTail, nowAction1);
        }
    }*/
    if (!flag) {
        // avoid crash
        if (calDist({Players[x ^ 1].x, Players[x ^ 1].y}, {DecisionSystem[x].p.frameState.x, DecisionSystem[x].p.frameState.y}) < 3) {
            frameAction->append("Move");
            frameAction->append(" ");
            frameAction->append(DecisionSystem[x].playerMove(DecisionSystem[x].p, DecisionSystem[x].actionTarget, 1));
            return 1;
        }
    }
    return flag;
}

int findWork(int x)
{
    int validWork[20 + 5];
    for (int i = 0; i < orderCount; i++)
        validWork[i] = 1;
    if (!DecisionSystem[x ^ 1].p.isfree) {
        std::vector <std::string> nowOrder; nowOrder.clear();
        for (int i = 0; i < DecisionSystem[x ^ 1].fd.currentOrder.recipe.size() - 1; i++)
            nowOrder.push_back(DecisionSystem[x ^ 1].fd.currentOrder.recipe[i]);
        for (int i = 0; i < orderCount; i++) {
            bool flag = 1;
            if (nowOrder.size() != Order[i].recipe.size()) flag = 0;
            else {
                for (int j = 0; j < nowOrder.size(); j++)
                    if (nowOrder[j] != Order[i].recipe[j]) flag = 0;
            }
            if (flag) {
                validWork[i] = 0;
                break;
            }
        }
    }
    int now = 0;
    while (!validWork[now]) now++;
    return now;
}

void init_all()
{
    for (int i = 0; i < 2; i++) 
        DecisionSystem[i].init(i);
    int countNum = 0;
    for (int i = 0; i < width; i++)
        for (int j = 0; j < height; j++) {
            if (Map[i][j] == 'c') chopplace = {i, j};
            else if (Map[i][j] == '$') submitplace = {i, j};
            else if (Map[i][j] == 'p') dirtyplateshome = {i, j};
            else if (Map[i][j] == 'k') washingplace = {i, j};
            else if (Map[i][j] == '_') deadpoint = {i, j};
            if (Map[i][j] != '.') 
                concretePoints.push_back({i, j});
        }
            
    for (int i = 0; i < entityCount; i++)
        if (Entity[i].entity[0] == "Pot") pothome = {Entity[i].x, Entity[i].y};
        else if (Entity[i].entity[0] == "Pan") panhome = {Entity[i].x, Entity[i].y};
        else if (Entity[i].entity[0] == "Plate") plateNum++;
    
    myFloyd.initFloyd();
    std::cerr << "attention!" << myFloyd.dist[2][7][1][4] << std::endl;
}

void init_read()
{
    std::string s;
    std::stringstream ss;
    int frame;

    /* 读取初始地图信息 */
    std::getline(std::cin, s, '\0');
    ss << s;

    /* 若按照该读入，访问坐标(x, y)等价于访问Map[y][x],你可按照自己的习惯进行修改 */
    // modified
    ss >> width >> height;
    std::cerr << "Map size: " << width << "x" << height << std::endl;
    for (int i = 0; i < height; i++)
        for (int j = 0; j < width; j++)
            ss >> Map[j][i];

    /* 读入原料箱：位置、名字、以及采购单价 */
    ss >> IngredientCount;
    for (int i = 0; i < IngredientCount; i++)
    {
        ss >> s;
        assert(s == "IngredientBox");
        ss >> Ingredient[i].x >> Ingredient[i].y >> Ingredient[i].name >> Ingredient[i].price;
    }

    /* 读入配方：加工时间、加工前的字符串表示、加工容器、加工后的字符串表示 */
    ss >> recipeCount;
    for (int i = 0; i < recipeCount; i++)
    {
        ss >> Recipe[i].time >> Recipe[i].nameBefore >> Recipe[i].kind >> Recipe[i].nameAfter;
        maxWaitingTime = std::max(maxWaitingTime, Recipe[i].time);
    }

    /* 读入总帧数、当前采用的随机种子、一共可能出现的订单数量 */
    ss >> totalTime >> randomizeSeed >> totalOrderCount;

    /* 读入订单的有效帧数、价格、权重、订单组成 */
    for (int i = 0; i < totalOrderCount; i++)
    {
        ss >> totalOrder[i].validFrame >> totalOrder[i].price >> totalOrder[i].frequency;
        getline(ss, s);
        std::stringstream tmp(s);
        while (tmp >> s)
            totalOrder[i].recipe.push_back(s);
    }

    /* 读入玩家信息：初始坐标 */
    ss >> k;
    assert(k == 2);
    for (int i = 0; i < k; i++)
    {
        ss >> Players[i].x >> Players[i].y;
        Players[i].containerKind = ContainerKind::None;
        Players[i].entity.clear();
    }

    /* 读入实体信息：坐标、实体组成 */
    ss >> entityCount;
    for (int i = 0; i < entityCount; i++)
    {
        ss >> Entity[i].x >> Entity[i].y >> s;
        Entity[i].entity.push_back(s);
    }
}

bool frame_read(int nowFrame)
{
    nowframe = nowFrame;
    std::string s;
    std::stringstream ss;
    int frame;
    std::getline(std::cin, s, '\0');
    ss.str(s);
    /*
      如果输入流中还有数据，说明游戏已经在请求下一帧了
      这时候我们应该跳过当前帧，以便能够及时响应游戏。
    */
    if (std::cin.rdbuf()->in_avail() > 0)
    {
        std::cerr << "Warning: skipping frame " << nowFrame
             << " to catch up with the game" << std::endl;
        return true;
    }
    ss >> s;
    assert(s == "Frame");
    int currentFrame;
    ss >> currentFrame;
    assert(currentFrame == nowFrame);
    ss >> remainFrame >> Fund;
    /* 读入当前的订单剩余帧数、价格、以及配方 */
    ss >> orderCount;
    for (int i = 0; i < orderCount; i++)
    {
        ss >> Order[i].validFrame >> Order[i].price;
        Order[i].recipe.clear();
        getline(ss, s);
        std::stringstream tmp(s);
        while (tmp >> s)
        {
            Order[i].recipe.push_back(s);
        }
    }
    ss >> k;
    assert(k == 2);
    /* 读入玩家坐标、x方向速度、y方向速度、剩余复活时间 */
    for (int i = 0; i < k; i++)
    {
        ss >> Players[i].x >> Players[i].y >> Players[i].X_Velocity >> Players[i].Y_Velocity >> Players[i].live;
        getline(ss, s);
        std::stringstream tmp(s);
        Players[i].containerKind = ContainerKind::None;
        Players[i].entity.clear();
        while (tmp >> s)
        {
            /*
                若若该玩家手里有东西，则接下来一个分号，分号后一个空格，空格后为一个实体。
                以下是可能的输入（省略前面的输入）：
                 ;  : fish
                 ; @  : fish
                 ; @ Plate : fish
                 ; Plate
                 ; DirtyPlates 1
                ...
            */

            /* 若你不需要处理这些，可直接忽略 */
            if (s == ";" || s == ":" || s == "@" || s == "*")
                continue;

            /* 
                Todo: 其他容器
            */
            if (s == "Plate")
                Players[i].containerKind = ContainerKind::Plate;
            else if (s == "Pan")
                Players[i].containerKind = ContainerKind::Pan;
            else if (s == "Pot")
                Players[i].containerKind = ContainerKind::Pot;
            else if (s == "DirtyPlates")
                Players[i].containerKind = ContainerKind::DirtyPlates;
            else
                Players[i].entity.push_back(s);
        }
    }

    ss >> entityCount;
    /* 读入实体坐标 */
    for (int i = 0; i < entityCount; i++)
    {
        ss >> Entity[i].x >> Entity[i].y;
        getline(ss, s);
        std::stringstream tmp(s);
        Entity[i].containerKind = ContainerKind::None;
        Entity[i].entity.clear();
        Entity[i].currentFrame = Entity[i].totalFrame = 0;
        Entity[i].sum = 1;
        while (tmp >> s)
        {
            /*
                读入一个实体，例子：
                DirtyPlates 2
                fish
                DirtyPlates 1 ; 15 / 180

            */

            /* 若你不需要处理这些，可直接忽略 */
            if (s == ":" || s == "@" || s == "*")
                continue;
            if (s == ";")
            {
                tmp >> Entity[i].currentFrame >> s >> Entity[i].totalFrame;
                assert(s == "/");
                break;
            }
            
            /* 
                Todo: 其他容器
            */
            if (s == "Plate")
                Entity[i].containerKind = ContainerKind::Plate;
            else if (s == "Pan")
                Entity[i].containerKind = ContainerKind::Pan;
            else if (s == "Pot")
                Entity[i].containerKind = ContainerKind::Pot;
            else if (s == "DirtyPlates")
            {
                Entity[i].containerKind = ContainerKind::DirtyPlates;
                tmp >> Entity[i].sum;
            }
            else
                Entity[i].entity.push_back(s);
        }
    }
    return false;
}