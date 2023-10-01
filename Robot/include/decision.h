#ifndef DECISION
#define DECISION

#include <string>
#include <vector>
#include <queue>
#include <framework.h>

bool sameCell(std::pair<double, double> p1, std::pair<double, double> p2)
{
    std::pair<int, int> q1 = {(int)p1.first, (int)p1.second};
    std::pair<int, int> q2 = {(int)p2.first, (int)p2.second};
    return q1 == q2;
}

struct FullPlayer
{
    Player frameState;
    bool isfree = 1;
    bool shouldMoving = 0;
    std::pair<double, double> movingPlace;
    bool fixedcp[3 + 5];


    void frameFlush(Player nowFrame)
    {
        frameState = nowFrame;
    }

    void stopWork()
    {
        isfree = 1;
    }

    bool setMoving(std::pair<double, double> p)
    {
        shouldMoving = 1;
        movingPlace = p;
        return 1;
    }

    bool arrived()
    {
        if (sameCell(movingPlace, std::make_pair(frameState.x, frameState.y)))
            shouldMoving = 0;
        return shouldMoving ^ 1;
    }
};

#endif