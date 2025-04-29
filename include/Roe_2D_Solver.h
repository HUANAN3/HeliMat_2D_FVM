#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <iomanip>
#include <omp.h>

#include "DataExporter.h"


#ifdef GLOBAL_ACCURACY
  #define REAL float
  #define ABS std::fabsf
  #define SQRT std::sqrtf
  #define TAN std::tanf
  #define COS std::cosf
#else
  #define REAL double
  #define ABS std::abs
  #define SQRT std::sqrt
  #define TAN std::tan
  #define COS std::cos
#endif

#define PI 3.1415926536

class Roe2DSolver {
public:
  struct Grid {
    std::vector<REAL> X;
    std::vector<REAL> Y;
    std::vector<REAL> Xnx;
    std::vector<REAL> Xny;
    std::vector<REAL> Ynx;
    std::vector<REAL> Yny;
    size_t numGridX;
    size_t numGridY;
    REAL XGridStart;
    REAL XGridEnd;
    REAL YGridStart;
    REAL YGridEnd;
    size_t numGrid;
  } grid;
  struct FlowPara {
    std::vector<REAL> p;
    std::vector<REAL> rho;
    std::vector<REAL> u;
    std::vector<REAL> v;
  } flowPara;
  struct Flux {
    std::vector<REAL> rhoV;
    std::vector<REAL> rhouVp;
    std::vector<REAL> rhovVp;
    std::vector<REAL> rhoHV;
  } F, G;
  struct W {
    std::vector<REAL> rho;
    std::vector<REAL> rhou;
    std::vector<REAL> rhov;
    std::vector<REAL> rhoE;
  } wn;
  Roe2DSolver(REAL Xs, REAL Xe, size_t Xnum, REAL Ys, REAL Ye, size_t Ynum, REAL Tt,REAL dt, std::string FN);
  ~Roe2DSolver();
  auto SetGrid() -> void;
  auto SetInitialValue() -> void;
  auto CalculateFluxField() -> void;
  auto RungeKuttaTimeAdvance() -> void;
  auto Solver() ->void;
  auto SetBC() -> void;
  auto qLMUSCL(REAL qW, REAL q, REAL qE) -> REAL;
  auto qRMUSCL(REAL q, REAL qE, REAL qEE) -> REAL;
  auto upDateWn() -> void;
  auto upDateFlux() -> void; 
  auto VanAlbada(REAL delplus, REAL delminus) -> REAL;
  auto linspace(REAL start, REAL end, size_t num) -> std::vector<REAL>;
  auto idx(int row, int col) -> size_t;
private:
  REAL TotalTime;
  const REAL gamma = 1.4;
  const REAL dt;
  const size_t nG = 2;
  size_t timeStep;
  std::string filename;
};
