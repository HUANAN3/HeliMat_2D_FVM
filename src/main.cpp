#include "Roe_2D_Solver.h" 
auto main() -> int {
  try {
  omp_set_num_threads(24);  
  std::string filename = "D:/2024_Graduate/00_Code/C++/CFD_Project/result/Part2/Roe2D";
  // std::string filename = "R:/BangdaFreshOrange/Roe2D";
  //Roe2DSolver RSolver(0.0,4.0,4000,   0.0,1.0,1000,  0.3,0.00001,filename);
  Roe2DSolver RSolver(0.0,1.0,256,   0.0,1.0,256, 0.5,0.001,filename);
  RSolver.SetGrid();
  RSolver.SetInitialValue();
  RSolver.Solver();
  } catch (const std::exception& e) {
    std::cerr << "错误信息：" << e.what() << std::endl;
    return 1;
  }
  return 0;
}