#include "header.hpp"

void printData ( Particle * pFirst )
{
  Particle * pTmp = pFirst;
  // 粒子データの出力
  for(int i = 0; i < N_PTCL; i++){
    std::cout << (pTmp + i)->getMass() << " "
              << (pTmp + i)->getPos()  << " "
              << (pTmp + i)->getVel()  << std::endl;
  }

}
