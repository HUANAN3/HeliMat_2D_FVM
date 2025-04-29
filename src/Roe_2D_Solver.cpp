#include "Roe_2D_Solver.h"

Roe2DSolver::Roe2DSolver(REAL Xs, REAL Xe, size_t Xnum, 
                         REAL Ys, REAL Ye, size_t Ynum, 
                         REAL Tt,REAL dt, std::string FN) 
                        : TotalTime(Tt), dt(dt), filename(FN) {
                          grid.XGridStart = Xs;
                          grid.XGridEnd   = Xe;
                          grid.numGridX   = Xnum;
                          grid.YGridStart = Ys;
                          grid.YGridEnd   = Ye;
                          grid.numGridY   = Ynum;
                          std::cout << "**********************************" << std::endl;
                          std::cout << "                                  " << std::endl;
                          std::cout << "      HeliMat 2D FVM Solver       " << std::endl;
                          std::cout << "                                  " << std::endl;
                          std::cout << "     (c) J. Shen 2024-04-09       " << std::endl;
                          std::cout << "                                  " << std::endl;
                          std::cout << "**********************************" << std::endl;
                        }

Roe2DSolver::~Roe2DSolver() {}

auto Roe2DSolver::idx(int row, int col) -> size_t {  // IMPROTANT
  //return col + row * grid.X.size();
  return static_cast<size_t>(row) + nG + (static_cast<size_t>(col) + nG) * (grid.numGridX + 2*nG);
}

auto Roe2DSolver::SetGrid() -> void {
  grid.X = linspace(grid.XGridStart, grid.XGridEnd, grid.numGridX);
  grid.Y = linspace(grid.YGridStart, grid.YGridEnd, grid.numGridY);
  grid.numGrid = (grid.numGridX + 2*nG) * (grid.numGridY + 2*nG);
  grid.Xnx.resize(grid.numGrid, 0.0);
  grid.Xny.resize(grid.numGrid, 0.0);
  grid.Ynx.resize(grid.numGrid, 0.0);
  grid.Yny.resize(grid.numGrid, 0.0);
  //grid.numGrid = grid.X.size() * grid.Y.size();
  for (int i = -2; i < static_cast<int>(grid.numGridX) + 2; i++) {
    for (int j = -2; j < static_cast<int>(grid.numGridY) + 2; j++) {//// 
      grid.Xnx[idx(i,j)] = 1.0;
      grid.Xny[idx(i,j)] = 0.0;
      grid.Ynx[idx(i,j)] = 0.0;
      grid.Yny[idx(i,j)] = 1.0;
    }
  }
  std::cout << "Grid Generation Complete!" << std::endl;
  REAL max_speed = 1.5;
  auto CFL = dt*max_speed/(grid.XGridStart-grid.XGridEnd)*grid.numGridX;
  if (ABS(CFL) >= 1.0) throw std::invalid_argument("CFL条件不满足！");
}

auto Roe2DSolver::SetInitialValue() -> void {
  flowPara.p.resize(grid.numGrid, 0.0);       // 预先分配内存
  flowPara.rho.resize(grid.numGrid, 0.0);
  flowPara.u.resize(grid.numGrid, 0.0);
  flowPara.v.resize(grid.numGrid, 0.0);
  F.rhoV.resize(grid.numGrid, 0.0);
  F.rhouVp.resize(grid.numGrid, 0.0);
  F.rhovVp.resize(grid.numGrid, 0.0);
  F.rhoHV.resize(grid.numGrid, 0.0);
  G.rhoV.resize(grid.numGrid, 0.0);
  G.rhouVp.resize(grid.numGrid, 0.0);
  G.rhovVp.resize(grid.numGrid, 0.0);   
  G.rhoHV.resize(grid.numGrid, 0.0);
  wn.rho.resize(grid.numGrid, 0.0);
  wn.rhoE.resize(grid.numGrid, 0.0);
  wn.rhou.resize(grid.numGrid, 0.0);
  wn.rhov.resize(grid.numGrid, 0.0);
  for (size_t i = 0; i < grid.numGridX; i++) {
    for(size_t j = 0; j < grid.numGridY; j++) {
    // Initial Value Set Here
    //  if (grid.X[i] <= 1.0/6.0 + grid.Y[j]*TAN(PI/6.0)) {
    //   // if (grid.Y[j] > 0.1) {
    //     flowPara.u[idx(i, j)] = 7.145;        // 激波后
    //     flowPara.v[idx(i, j)] = -4.125;
    //     flowPara.rho[idx(i, j)] = 8.0;       //7.875
    //     flowPara.p[idx(i, j)] = 116.5;       //87.33
    //     } else {
    //     flowPara.u[idx(i, j)] = 0.0;          // 激波前
    //     flowPara.v[idx(i, j)] = 0.0;
    //     flowPara.rho[idx(i, j)] = 1.4;
    //     flowPara.p[idx(i, j)] = 1.0;
    //   }

    // if (grid.X[i] >= 0.75) {
    //   if (grid.Y[j] >= 0.75) {
    //     flowPara.u[idx(i, j)] = 0.0;
    //     flowPara.v[idx(i, j)] = 0.0;
    //     flowPara.rho[idx(i, j)] = 1.0;
    //     flowPara.p[idx(i, j)] = 1.0;
    //   } else {
    //     flowPara.u[idx(i, j)] = 0.0;
    //     flowPara.v[idx(i, j)] = 1.0;
    //     flowPara.rho[idx(i, j)] = 0.5;
    //     flowPara.p[idx(i, j)] = 0.25;
    //   }
    // } else {
    //   if (grid.Y[j] >= 0.75) {
    //     flowPara.u[idx(i, j)] = 1.0;
    //     flowPara.v[idx(i, j)] = 0.0;
    //     flowPara.rho[idx(i, j)] = 0.5;
    //     flowPara.p[idx(i, j)] = 0.25;
    //   } else {
    //     flowPara.u[idx(i, j)] = 0.0;
    //     flowPara.v[idx(i, j)] = 1.0;
    //     flowPara.rho[idx(i, j)] = 0.125;
    //     flowPara.p[idx(i, j)] = 0.25;
    //   }
    // }
    }
  }
  upDateWn();
  upDateFlux();
}

auto Roe2DSolver::SetBC() -> void {
  const REAL u0 = 0.0;      // 激波前参数
  const REAL v0 = 0.0;
  const REAL rho0 = 1.4;
  const REAL p0 = 1.0;
  const REAL u1 = 7.145;     // 激波后参数
  const REAL v1 = -4.125;
  const REAL rho1 = 8.0;
  const REAL p1 = 116.5;
  // 左边界条件
  for (size_t i = 0; i < nG; i++) {
    for (size_t j = 0; j < grid.numGridY; j++) {
        int im = - 1 - i;
        // flowPara.rho[idx(im,j)] = rho1;
        // flowPara.p[idx(im,j)] = p1;
        // flowPara.u[idx(im,j)] = u1;
        // flowPara.v[idx(im,j)] = v1;
        flowPara.rho[idx(im,j)] = flowPara.rho[idx(0,j)];
        flowPara.p[idx(im,j)] = flowPara.p[idx(0,j)];
        flowPara.u[idx(im,j)] = flowPara.u[idx(0,j)];
        flowPara.v[idx(im,j)] = flowPara.v[idx(0,j)];
    }
  }  
  // 右边界条件 
  for (size_t i = grid.numGridX - nG; i < grid.numGridX; i++) {
    for (size_t j = 0; j < grid.numGridY; j++) {
        int im = i + nG;
        flowPara.rho[idx(im,j)] = flowPara.rho[idx(grid.numGridX - 1,j)];
        flowPara.p[idx(im,j)] = flowPara.p[idx(grid.numGridX - 1,j)];
        flowPara.u[idx(im,j)] = flowPara.u[idx(grid.numGridX - 1,j)];
        flowPara.v[idx(im,j)] = flowPara.v[idx(grid.numGridX - 1,j)];
    }
  } 
  // 上边界条件
  for (size_t i = 0; i < grid.numGridX; i++) {
    for (size_t j = grid.numGridY - nG; j < grid.numGridY; j++) {
      int jm = j + nG;
    //   if (grid.X[i] <= 1.0/6.0 + grid.Y[j]*TAN(PI/6.0) + 10.0/COS(PI/6.0)*timeStep*dt){  // IMPORTANT
    //     flowPara.rho[idx(i,jm)] = rho1;
    //     flowPara.p[idx(i,jm)] = p1;
    //     flowPara.u[idx(i,jm)] = u1;
    //     flowPara.v[idx(i,jm)] = v1;
    //   } else {
    //     flowPara.rho[idx(i,jm)] = flowPara.rho[idx(i,grid.numGridY - 1)];
    //     flowPara.p[idx(i,jm)] = flowPara.p[idx(i,grid.numGridY - 1)];
    //     flowPara.u[idx(i,jm)] = flowPara.u[idx(i,grid.numGridY - 1)];
    //     flowPara.v[idx(i,jm)] = flowPara.v[idx(i,grid.numGridY - 1)];
    //   }
    // 
    flowPara.rho[idx(i, jm)] = flowPara.rho[idx(i, grid.numGridY - 1)];
    flowPara.p[idx(i, jm)] = flowPara.p[idx(i, grid.numGridY - 1)];
    flowPara.u[idx(i, jm)] = flowPara.u[idx(i, grid.numGridY - 1)];
    flowPara.v[idx(i, jm)] = flowPara.v[idx(i, grid.numGridY - 1)];
    }
  } 
  // 下边界条件
  for (size_t i = 0; i < grid.numGridX; i++) {
    for (size_t j = 0; j < nG; j++) {
  //     if (grid.X[i] <= 1.0/6.0) {
  //       int jm = - 1 - j;
  //       flowPara.rho[idx(i,jm)] = rho1;
  //       flowPara.p[idx(i,jm)] = p1;
  //       flowPara.u[idx(i,jm)] = u1;
  //       flowPara.v[idx(i,jm)] = v1;
  //     } else {
  //       int jm = - 1 - j;
  //       flowPara.rho[idx(i,jm)] = flowPara.rho[idx(i,0)];  // 固壁条件
  //       flowPara.p[idx(i,jm)] = flowPara.p[idx(i,0)];
  //       flowPara.u[idx(i,jm)] = flowPara.u[idx(i,0)];
  //       flowPara.v[idx(i,jm)] = -flowPara.v[idx(i,0)];     // IMPORTANT
  //     }
  //   }
  //  
  int jm = - 1 - j;
  flowPara.rho[idx(i, jm)] = flowPara.rho[idx(i, 0)];
  flowPara.p[idx(i, jm)] = flowPara.p[idx(i, 0)];
  flowPara.u[idx(i, jm)] = flowPara.u[idx(i, 0)];
  flowPara.v[idx(i, jm)] = flowPara.v[idx(i, 0)];
    }
  }
}


auto Roe2DSolver::CalculateFluxField() -> void {
  // BC
  SetBC();
  // InitialValue
  auto rhoLX = flowPara.rho;
  auto rhoLY = flowPara.rho;
  auto pLX = flowPara.p;
  auto pLY = flowPara.p;
  auto uLX = flowPara.u;
  auto uLY = flowPara.u;
  auto vLX = flowPara.v;
  auto vLY = flowPara.v;

  #pragma omp parallel for collapse(2)
  for (int i = -1; i < static_cast<int>(grid.numGridX); i++) {     // size_t -> int
    for(int j = -1; j < static_cast<int>(grid.numGridY); j++) { 
      rhoLX[idx(i,j)] = qLMUSCL(flowPara.rho[idx(i-1,j)], flowPara.rho[idx(i,j)], flowPara.rho[idx(i+1,j)]);
      rhoLY[idx(i,j)] = qLMUSCL(flowPara.rho[idx(i,j-1)], flowPara.rho[idx(i,j)], flowPara.rho[idx(i,j+1)]);
      pLX[idx(i,j)] = qLMUSCL(flowPara.p[idx(i-1,j)], flowPara.p[idx(i,j)], flowPara.p[idx(i+1,j)]);
      pLY[idx(i,j)] = qLMUSCL(flowPara.p[idx(i,j-1)], flowPara.p[idx(i,j)], flowPara.p[idx(i,j+1)]);
      uLX[idx(i,j)] = qLMUSCL(flowPara.u[idx(i-1,j)], flowPara.u[idx(i,j)], flowPara.u[idx(i+1,j)]);
      uLY[idx(i,j)] = qLMUSCL(flowPara.u[idx(i,j-1)], flowPara.u[idx(i,j)], flowPara.u[idx(i,j+1)]);
      vLX[idx(i,j)] = qLMUSCL(flowPara.v[idx(i-1,j)], flowPara.v[idx(i,j)], flowPara.v[idx(i+1,j)]);
      vLY[idx(i,j)] = qLMUSCL(flowPara.v[idx(i,j-1)], flowPara.v[idx(i,j)], flowPara.v[idx(i,j+1)]);
    }
  }

  auto rhoRX = flowPara.rho;
  auto rhoRY = flowPara.rho;
  auto pRX = flowPara.p;
  auto pRY = flowPara.p;
  auto uRX = flowPara.u;
  auto uRY = flowPara.u;
  auto vRX = flowPara.v;
  auto vRY = flowPara.v;

  #pragma omp parallel for collapse(2)
  for (int i = -1; i < static_cast<int>(grid.numGridX); i++) {        // size_t -> int
    for(int j = -1; j < static_cast<int>(grid.numGridY); j++) { 
      rhoRX[idx(i,j)] = qRMUSCL(flowPara.rho[idx(i,j)], flowPara.rho[idx(i+1,j)], flowPara.rho[idx(i+2,j)]);
      rhoRY[idx(i,j)] = qRMUSCL(flowPara.rho[idx(i,j)], flowPara.rho[idx(i,j+1)], flowPara.rho[idx(i,j+2)]);
      pRX[idx(i,j)] = qRMUSCL(flowPara.p[idx(i,j)], flowPara.p[idx(i+1,j)], flowPara.p[idx(i+2,j)]);
      pRY[idx(i,j)] = qRMUSCL(flowPara.p[idx(i,j)], flowPara.p[idx(i,j+1)], flowPara.p[idx(i,j+2)]);
      uRX[idx(i,j)] = qRMUSCL(flowPara.u[idx(i,j)], flowPara.u[idx(i+1,j)], flowPara.u[idx(i+2,j)]);
      uRY[idx(i,j)] = qRMUSCL(flowPara.u[idx(i,j)], flowPara.u[idx(i,j+1)], flowPara.u[idx(i,j+2)]);
      vRX[idx(i,j)] = qRMUSCL(flowPara.v[idx(i,j)], flowPara.v[idx(i+1,j)], flowPara.v[idx(i+2,j)]);
      vRY[idx(i,j)] = qRMUSCL(flowPara.v[idx(i,j)], flowPara.v[idx(i,j+1)], flowPara.v[idx(i,j+2)]);
    }
  }
  auto rhoAvX = flowPara.rho;
  auto rhoAvY = flowPara.rho;
  auto uAvX = flowPara.u;
  auto uAvY = flowPara.u;
  auto vAvX = flowPara.v;
  auto vAvY = flowPara.v;
  auto HLX = flowPara.p;
  auto HLY = flowPara.p;
  auto HRX = flowPara.p;
  auto HRY = flowPara.p;
  auto HAvX = flowPara.p;
  auto HAvY = flowPara.p;
  auto cAvX = flowPara.p;
  auto cAvY = flowPara.p;
  auto VAvX = flowPara.u;
  auto VAvY = flowPara.v;
  auto qAv2X = flowPara.u;
  auto qAv2Y = flowPara.v;
  auto FluxLX = F;
  auto FluxRX = F;
  auto FluxLY = G;
  auto FluxRY = G;
  auto VLX = flowPara.u;
  auto VLY = flowPara.u;
  auto VRX = flowPara.v;
  auto VRY = flowPara.v;

  auto Ac = [](REAL ac, REAL cav) -> REAL {
    auto del = 0.1*cav;
    return ac > del ? ac : ((ac*ac + del*del) / (2*del));
  };
  
  #pragma omp parallel for simd                      // 合并两层循环
  for (size_t idx = 0; idx < grid.numGrid; idx++) {
    rhoAvX[idx] = SQRT(rhoLX[idx] * rhoRX[idx]);
    rhoAvY[idx] = SQRT(rhoLY[idx] * rhoRY[idx]);
    uAvX[idx] = (uLX[idx] * SQRT(rhoLX[idx]) + uRX[idx] * SQRT(rhoRX[idx])) / (SQRT(rhoLX[idx]) + SQRT(rhoRX[idx]));
    uAvY[idx] = (uLY[idx] * SQRT(rhoLY[idx]) + uRY[idx] * SQRT(rhoRY[idx])) / (SQRT(rhoLY[idx]) + SQRT(rhoRY[idx]));
    vAvX[idx] = (vLX[idx] * SQRT(rhoLX[idx]) + vRX[idx] * SQRT(rhoRX[idx])) / (SQRT(rhoLX[idx]) + SQRT(rhoRX[idx]));
    vAvY[idx] = (vLY[idx] * SQRT(rhoLY[idx]) + vRY[idx] * SQRT(rhoRY[idx])) / (SQRT(rhoLY[idx]) + SQRT(rhoRY[idx]));
    VAvX[idx] = uAvX[idx] * grid.Xnx[idx] + vAvX[idx] * grid.Xny[idx];
    VAvY[idx] = uAvY[idx] * grid.Ynx[idx] + vAvY[idx] * grid.Yny[idx];
    qAv2X[idx] = uAvX[idx] * uAvX[idx] + vAvX[idx] * vAvX[idx];
    qAv2Y[idx] = uAvY[idx] * uAvY[idx] + vAvY[idx] * vAvY[idx];
    HLX[idx] = gamma/(gamma-1) * pLX[idx]/rhoLX[idx] + 0.5*uLX[idx]*uLX[idx] + 0.5*vLX[idx]*vLX[idx];
    HLY[idx] = gamma/(gamma-1) * pLY[idx]/rhoLY[idx] + 0.5*uLY[idx]*uLY[idx] + 0.5*vLY[idx]*vLY[idx];
    HRX[idx] = gamma/(gamma-1) * pRX[idx]/rhoRX[idx] + 0.5*uRX[idx]*uRX[idx] + 0.5*vRX[idx]*vRX[idx];
    HRY[idx] = gamma/(gamma-1) * pRY[idx]/rhoRY[idx] + 0.5*uRY[idx]*uRY[idx] + 0.5*vRY[idx]*vRY[idx];
    HAvX[idx] = (HLX[idx] * SQRT(rhoLX[idx]) + HRX[idx] * SQRT(rhoRX[idx])) / (SQRT(rhoLX[idx]) + SQRT(rhoRX[idx]));
    HAvY[idx] = (HLY[idx] * SQRT(rhoLY[idx]) + HRY[idx] * SQRT(rhoRY[idx])) / (SQRT(rhoLY[idx]) + SQRT(rhoRY[idx]));
    cAvX[idx] = SQRT((gamma-1) * (HAvX[idx] - 0.5 * qAv2X[idx]));
    cAvY[idx] = SQRT((gamma-1) * (HAvY[idx] - 0.5 * qAv2Y[idx]));
    VLX[idx] = uLX[idx] * grid.Xnx[idx] + vLX[idx] * grid.Xny[idx];
    VLY[idx] = uLY[idx] * grid.Ynx[idx] + vLY[idx] * grid.Yny[idx];
    VRX[idx] = uRX[idx] * grid.Xnx[idx] + vRX[idx] * grid.Xny[idx];
    VRY[idx] = uRY[idx] * grid.Ynx[idx] + vRY[idx] * grid.Yny[idx];
    FluxLX.rhoV[idx] = rhoLX[idx] * VLX[idx];
    FluxLX.rhouVp[idx] = rhoLX[idx] * uLX[idx] * VLX[idx] + grid.Xnx[idx] * pLX[idx];
    FluxLX.rhovVp[idx] = rhoLX[idx] * vLX[idx] * VLX[idx] + grid.Xny[idx] * pLX[idx];
    FluxLX.rhoHV[idx] = rhoLX[idx] * HLX[idx] * VLX[idx];

    FluxRX.rhoV[idx] = rhoRX[idx] * VRX[idx];
    FluxRX.rhouVp[idx] = rhoRX[idx] * uRX[idx] * VRX[idx] + grid.Xnx[idx] * pRX[idx];
    FluxRX.rhovVp[idx] = rhoRX[idx] * vRX[idx] * VRX[idx] + grid.Xny[idx] * pRX[idx];
    FluxRX.rhoHV[idx] = rhoRX[idx] * HRX[idx] * VRX[idx];

    FluxLY.rhoV[idx] = rhoLY[idx] * VLY[idx];
    FluxLY.rhouVp[idx] = rhoLY[idx] * uLY[idx] * VLY[idx] + grid.Ynx[idx] * pLY[idx];
    FluxLY.rhovVp[idx] = rhoLY[idx] * vLY[idx] * VLY[idx] + grid.Yny[idx] * pLY[idx];
    FluxLY.rhoHV[idx] = rhoLY[idx] * HLY[idx] * VLY[idx];

    FluxRY.rhoV[idx] = rhoRY[idx] * VRY[idx];
    FluxRY.rhouVp[idx] = rhoRY[idx] * uRY[idx] * VRY[idx] + grid.Ynx[idx] * pRY[idx];
    FluxRY.rhovVp[idx] = rhoRY[idx] * vRY[idx] * VRY[idx] + grid.Yny[idx] * pRY[idx];
    FluxRY.rhoHV[idx] = rhoRY[idx] * HRY[idx] * VRY[idx];
  }

  #pragma omp parallel for simd   
  for (size_t idx = 0; idx < grid.numGrid; idx++) {
    auto dpX = pRX[idx] - pLX[idx];
    auto dpY = pRY[idx] - pLY[idx];
    auto dVX = VRX[idx] - VLX[idx];
    auto dVY = VRY[idx] - VLY[idx];
    auto drhoX = rhoRX[idx] - rhoLX[idx];
    auto drhoY = rhoRY[idx] - rhoLY[idx];
    auto cAv2X = cAvX[idx] * cAvX[idx];
    auto cAv2Y = cAvY[idx] * cAvY[idx];
    auto duX = uRX[idx] - uLX[idx];
    auto duY = uRY[idx] - uLY[idx];
    auto dvX = vRX[idx] - vLX[idx];
    auto dvY = vRY[idx] - vLY[idx];

    auto f1termX = Ac(ABS(VAvX[idx]-cAvX[idx]), cAvX[idx]) * (dpX - rhoAvX[idx]*cAvX[idx]*dVX) / (2*cAv2X); 
    auto f1termY = Ac(ABS(VAvY[idx]-cAvY[idx]), cAvY[idx]) * (dpY - rhoAvY[idx]*cAvY[idx]*dVY) / (2*cAv2Y); 
    auto f234term1X = Ac(ABS(VAvX[idx]), cAvX[idx]) * (drhoX - dpX / cAv2X);
    auto f234term1Y = Ac(ABS(VAvY[idx]), cAvY[idx]) * (drhoY - dpY / cAv2Y);
    auto f234term2X = Ac(ABS(VAvX[idx]), cAvX[idx]) * rhoAvX[idx];
    auto f234term2Y = Ac(ABS(VAvY[idx]), cAvY[idx]) * rhoAvY[idx];
    auto f5termX = Ac(ABS(VAvX[idx] + cAvX[idx]), cAvX[idx]) * (dpX + rhoAvX[idx]*cAvX[idx]*dVX) / (2*cAv2X);
    auto f5termY = Ac(ABS(VAvY[idx] + cAvY[idx]), cAvY[idx]) * (dpY + rhoAvY[idx]*cAvY[idx]*dVY) / (2*cAv2Y);

    F.rhoV[idx] = 0.5*(FluxLX.rhoV[idx] + FluxRX.rhoV[idx] - f1termX*1 - f234term1X*1 - f5termX*1);
    G.rhoV[idx] = 0.5*(FluxLY.rhoV[idx] + FluxRY.rhoV[idx] - f1termY*1 - f234term1Y*1 - f5termY*1);

    F.rhouVp[idx] = 0.5*(FluxLX.rhouVp[idx] + FluxRX.rhouVp[idx] - f1termX*(uAvX[idx] - cAvX[idx]*grid.Xnx[idx]) - f234term1X*uAvX[idx] 
    - f234term2X*(duX - dVX*grid.Xnx[idx]) - f5termX*(uAvX[idx] + cAvX[idx]*grid.Xnx[idx]));
    G.rhouVp[idx] = 0.5*(FluxLY.rhouVp[idx] + FluxRY.rhouVp[idx] - f1termY*(uAvY[idx] - cAvY[idx]*grid.Ynx[idx]) - f234term1Y*uAvY[idx] 
    - f234term2Y*(duY - dVY*grid.Ynx[idx]) - f5termY*(uAvY[idx] + cAvY[idx]*grid.Ynx[idx]));

    F.rhovVp[idx] = 0.5*(FluxLX.rhovVp[idx] + FluxRX.rhovVp[idx] - f1termX*(vAvX[idx] - cAvX[idx]*grid.Xny[idx]) - f234term1X*vAvX[idx] 
    - f234term2X*(dvX - dVX*grid.Xny[idx]) - f5termX*(vAvX[idx] + cAvX[idx]*grid.Xny[idx]));
    G.rhovVp[idx] = 0.5*(FluxLY.rhovVp[idx] + FluxRY.rhovVp[idx] - f1termY*(vAvY[idx] - cAvY[idx]*grid.Yny[idx]) - f234term1Y*vAvY[idx] 
    - f234term2Y*(dvY - dVY*grid.Yny[idx]) - f5termY*(vAvY[idx] + cAvY[idx]*grid.Yny[idx]));

    F.rhoHV[idx] = 0.5*(FluxLX.rhoHV[idx] + FluxRX.rhoHV[idx] - f1termX*(HAvX[idx] - cAvX[idx]*VAvX[idx]) - f234term1X*0.5*qAv2X[idx] 
    - f234term2X*(uAvX[idx]*duX + vAvX[idx]*dvX - VAvX[idx]*dVX) - f5termX*(HAvX[idx] + cAvX[idx]*VAvX[idx])); 
    G.rhoHV[idx] = 0.5*(FluxLY.rhoHV[idx] + FluxRY.rhoHV[idx] - f1termY*(HAvY[idx] - cAvY[idx]*VAvY[idx]) - f234term1Y*0.5*qAv2Y[idx] 
    - f234term2Y*(uAvY[idx]*duY + vAvY[idx]*dvY - VAvY[idx]*dVY) - f5termY*(HAvY[idx] + cAvY[idx]*VAvY[idx])); 
  }
}

auto Roe2DSolver::RungeKuttaTimeAdvance() -> void {
  auto dx = grid.X[1] - grid.X[0];
  auto dy = grid.Y[1] - grid.Y[0];
  upDateWn();
  auto Un = wn;
  CalculateFluxField();
  auto U1 = wn;

  #pragma omp parallel for collapse(2)               // ！
  for (size_t i = 0; i < grid.numGridX; i++) {
    for (size_t j = 0; j < grid.numGridY; j++) {
      U1.rho[idx(i,j)] = Un.rho[idx(i,j)] - dt * (F.rhoV[idx(i,j)] - F.rhoV[idx(i-1,j)]) / dx
        - dt * (G.rhoV[idx(i,j)] - G.rhoV[idx(i,j-1)]) / dy;  
      U1.rhou[idx(i,j)] = Un.rhou[idx(i,j)] - dt * (F.rhouVp[idx(i,j)] - F.rhouVp[idx(i-1,j)]) / dx
        - dt * (G.rhouVp[idx(i,j)] - G.rhouVp[idx(i,j-1)]) / dy;
      U1.rhov[idx(i,j)] = Un.rhov[idx(i,j)] - dt * (F.rhovVp[idx(i,j)] - F.rhovVp[idx(i-1,j)]) / dx
        - dt * (G.rhovVp[idx(i,j)] - G.rhovVp[idx(i,j-1)]) / dy;
      U1.rhoE[idx(i,j)] = Un.rhoE[idx(i,j)] - dt * (F.rhoHV[idx(i,j)] - F.rhoHV[idx(i-1,j)]) / dx
        - dt * (G.rhoHV[idx(i,j)] - G.rhoHV[idx(i,j-1)]) / dy;

      flowPara.rho[idx(i,j)] = U1.rho[idx(i,j)];
      flowPara.u[idx(i,j)] = U1.rhou[idx(i,j)] / U1.rho[idx(i,j)];
      flowPara.v[idx(i,j)] = U1.rhov[idx(i,j)] / U1.rho[idx(i,j)];
      flowPara.p[idx(i,j)] = (gamma - 1) * (U1.rhoE[idx(i,j)] - 0.5 * U1.rhou[idx(i,j)] * U1.rhou[idx(i,j)] / U1.rho[idx(i,j)] 
        - 0.5 * U1.rhov[idx(i,j)] * U1.rhov[idx(i,j)] / U1.rho[idx(i,j)]);
    }
  }
  CalculateFluxField();
  upDateWn();
  auto U2 = wn;

  #pragma omp parallel for collapse(2) 
  for (size_t i = 0; i < grid.numGridX; i++) {
    for (size_t j = 0; j < grid.numGridY; j++) {
      U2.rho[idx(i,j)] = 0.75*Un.rho[idx(i,j)] + 0.25*U1.rho[idx(i,j)] - 0.25*dt * (F.rhoV[idx(i,j)] - F.rhoV[idx(i-1,j)]) / dx
        - dt * (G.rhoV[idx(i,j)] - G.rhoV[idx(i,j-1)]) / dy;  
      U2.rhou[idx(i,j)] = 0.75*Un.rhou[idx(i,j)] + 0.25*U1.rhou[idx(i,j)] - 0.25*dt * (F.rhouVp[idx(i,j)] - F.rhouVp[idx(i-1,j)]) / dx
        - dt * (G.rhouVp[idx(i,j)] - G.rhouVp[idx(i,j-1)]) / dy;
      U2.rhov[idx(i,j)] = 0.75*Un.rhov[idx(i,j)] + 0.25*U1.rhov[idx(i,j)] - 0.25*dt * (F.rhovVp[idx(i,j)] - F.rhovVp[idx(i-1,j)]) / dx
        - dt * (G.rhovVp[idx(i,j)] - G.rhovVp[idx(i,j-1)]) / dy;
      U2.rhoE[idx(i,j)] = 0.75*Un.rhoE[idx(i,j)] + 0.25*U1.rhoE[idx(i,j)] - 0.25*dt * (F.rhoHV[idx(i,j)] - F.rhoHV[idx(i-1,j)]) / dx
        - dt * (G.rhoHV[idx(i,j)] - G.rhoHV[idx(i,j-1)]) / dy;

      flowPara.rho[idx(i,j)] = U2.rho[idx(i,j)];
      flowPara.u[idx(i,j)] = U2.rhou[idx(i,j)] / U2.rho[idx(i,j)];
      flowPara.v[idx(i,j)] = U2.rhov[idx(i,j)] / U2.rho[idx(i,j)];
      flowPara.p[idx(i,j)] = (gamma - 1) * (U2.rhoE[idx(i,j)] - 0.5 * U2.rhou[idx(i,j)] * U2.rhou[idx(i,j)] / U2.rho[idx(i,j)] 
        - 0.5 * U2.rhov[idx(i,j)] * U2.rhov[idx(i,j)] / U2.rho[idx(i,j)]);
    }
  }
  CalculateFluxField();
  upDateWn();
  auto U3 = wn;

  #pragma omp parallel for collapse(2)  
  for (size_t i = 0; i < grid.numGridX; i++) {
    for (size_t j = 0; j < grid.numGridY; j++) {
      U3.rho[idx(i,j)] = 1.0/3.0*Un.rho[idx(i,j)] + 2.0/3.0*U2.rho[idx(i,j)] - 2.0/3.0*dt * (F.rhoV[idx(i,j)] - F.rhoV[idx(i-1,j)]) / dx
        - dt * (G.rhoV[idx(i,j)] - G.rhoV[idx(i,j-1)]) / dy;  
      U3.rhou[idx(i,j)] = 1.0/3.0*Un.rhou[idx(i,j)] + 2.0/3.0*U2.rhou[idx(i,j)] - 2.0/3.0*dt * (F.rhouVp[idx(i,j)] - F.rhouVp[idx(i-1,j)]) / dx
        - dt * (G.rhouVp[idx(i,j)] - G.rhouVp[idx(i,j-1)]) / dy;
      U3.rhov[idx(i,j)] = 1.0/3.0*Un.rhov[idx(i,j)] + 2.0/3.0*U2.rhov[idx(i,j)] - 2.0/3.0*dt * (F.rhovVp[idx(i,j)] - F.rhovVp[idx(i-1,j)]) / dx
        - dt * (G.rhovVp[idx(i,j)] - G.rhovVp[idx(i,j-1)]) / dy;
      U3.rhoE[idx(i,j)] = 1.0/3.0*Un.rhoE[idx(i,j)] + 2.0/3.0*U2.rhoE[idx(i,j)] - 2.0/3.0*dt * (F.rhoHV[idx(i,j)] - F.rhoHV[idx(i-1,j)]) / dx
        - dt * (G.rhoHV[idx(i,j)] - G.rhoHV[idx(i,j-1)]) / dy;

      flowPara.rho[idx(i,j)] = U3.rho[idx(i,j)];
      flowPara.u[idx(i,j)] = U3.rhou[idx(i,j)] / U3.rho[idx(i,j)];
      flowPara.v[idx(i,j)] = U3.rhov[idx(i,j)] / U3.rho[idx(i,j)];
      flowPara.p[idx(i,j)] = (gamma - 1) * (U3.rhoE[idx(i,j)] - 0.5 * U3.rhou[idx(i,j)] * U3.rhou[idx(i,j)] / U3.rho[idx(i,j)] 
        - 0.5 * U3.rhov[idx(i,j)] * U3.rhov[idx(i,j)] / U3.rho[idx(i,j)]);
    }
  }
}

auto Roe2DSolver::Solver() -> void {
  std::cout << "Calculation Start!" << std::endl;
  DataExporter dep(filename);
  int outPutNume = 0;
  for (this->timeStep = 0; timeStep <= TotalTime / dt; timeStep++) {
    RungeKuttaTimeAdvance();
    auto fp = flowPara;
    if (timeStep % 100 == 0) {
      dep.writeTimeStep2D(outPutNume,fp.p,fp.u,fp.v,fp.rho);
      outPutNume++;
    }
    if (timeStep % 50 == 0) {
      std::cout << std::fixed << std::setprecision(2) <<
        static_cast<REAL>(timeStep) / TotalTime * dt * 100 << "%" << std::endl;
    }
  }
  std::cout << "Calculation Complete!" << std::endl;
  std::cout << "Thank you for your help!" << std::endl;
}

auto Roe2DSolver::qLMUSCL(REAL qW, REAL q, REAL qE) -> REAL {
  REAL k = 1.0/3.0;
    auto delForw = qE - q;
    auto delBack= q - qW;
    auto psi = VanAlbada(delForw, delBack);
    return q + 0.25*psi*((1-k)*delBack + (1+k)*delForw);
}

auto Roe2DSolver::qRMUSCL(REAL q, REAL qE, REAL qEE) -> REAL {
  REAL k = 1.0/3.0;
    auto delForw = qEE - qE;
    auto delBack= qE - q;
    auto psi = VanAlbada(delForw, delBack);
    return qE - 0.25*psi*((1-k)*delForw + (1+k)*delBack);
}

auto Roe2DSolver::upDateWn() -> void {
  for (size_t idx = 0; idx < grid.numGrid; idx++) {
    wn.rho[idx] = flowPara.rho[idx];
    wn.rhou[idx] = flowPara.rho[idx] * flowPara.u[idx];
    wn.rhov[idx] = flowPara.rho[idx] * flowPara.v[idx];
    wn.rhoE[idx] = flowPara.rho[idx] * (flowPara.p[idx] / (gamma - 1) / flowPara.rho[idx] + 0.5 * flowPara.u[idx] * flowPara.u[idx] + 0.5 * flowPara.v[idx] * flowPara.v[idx]);
  }
}

auto Roe2DSolver::upDateFlux() -> void {
  for (size_t i = 0; i < grid.numGridX; i++) {     
    for(size_t j = 0; j < grid.numGridY; j++) { 
      auto Vx = flowPara.u[idx(i,j)] * grid.Xnx[idx(i,j)] + flowPara.v[idx(i,j)] * grid.Xny[idx(i,j)];
      auto Vy = flowPara.u[idx(i,j)] * grid.Ynx[idx(i,j)] + flowPara.v[idx(i,j)] * grid.Yny[idx(i,j)];

      F.rhoV[idx(i,j)] = Vx * flowPara.rho[idx(i,j)];
      F.rhouVp[idx(i,j)] = Vx * flowPara.u[idx(i,j)] * flowPara.rho[idx(i,j)] + flowPara.p[idx(i,j)] * grid.Xnx[idx(i,j)];
      F.rhovVp[idx(i,j)] = Vx * flowPara.v[idx(i,j)] * flowPara.rho[idx(i,j)] + flowPara.p[idx(i,j)] * grid.Xny[idx(i,j)];
      F.rhoHV[idx(i,j)] = flowPara.rho[idx(i,j)] * Vx * (gamma*flowPara.p[idx(i,j)]/flowPara.rho[idx(i,j)]/(gamma-1) 
        + 0.5*flowPara.u[idx(i,j)]*flowPara.u[idx(i,j)] + 0.5*flowPara.v[idx(i,j)]*flowPara.v[idx(i,j)]);

      G.rhoV[idx(i,j)] = Vy * flowPara.rho[idx(i,j)];
      G.rhouVp[idx(i,j)] = Vy * flowPara.u[idx(i,j)] * flowPara.rho[idx(i,j)] + flowPara.p[idx(i,j)] * grid.Xnx[idx(i,j)];
      G.rhovVp[idx(i,j)] = Vy * flowPara.v[idx(i,j)] * flowPara.rho[idx(i,j)] + flowPara.p[idx(i,j)] * grid.Xny[idx(i,j)];
      G.rhoHV[idx(i,j)] = flowPara.rho[idx(i,j)] * Vy * (gamma*flowPara.p[idx(i,j)]/flowPara.rho[idx(i,j)]/(gamma-1) 
        + 0.5*flowPara.u[idx(i,j)]*flowPara.u[idx(i,j)] + 0.5*flowPara.v[idx(i,j)]*flowPara.v[idx(i,j)]);
    }
  } 
}

auto Roe2DSolver::VanAlbada(REAL delplus, REAL delminus) -> REAL {
  return (2*delplus*delminus+1.0e-6)/(delplus*delplus+delminus*delminus+1.0e-6);
};

auto Roe2DSolver::linspace(REAL start, REAL end, size_t num) -> std::vector<REAL> {
  std::vector<REAL> result;
  if (num <= 0) return result;
  if (num == 1) {
    result.push_back(start);
    return result;
  }
  double step = (end - start) / num;
  for (size_t idx = 0; idx <= num; ++idx) {
    result.push_back(start + step * idx);
  }
  return result;
}