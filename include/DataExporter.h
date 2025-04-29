#pragma once

#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <format>

class DataExporter {
 public:
  explicit DataExporter(const std::string& outputDir, bool clearDir = true);
  auto writeTimeStep(size_t timestep, const std::vector<double>& p, const std::vector<double>& u,
                     const std::vector<double>& rho) -> void;
  auto writeTimeStep2D(size_t timestep, const std::vector<double>& p, const std::vector<double>& u, const std::vector<double>& v,
                     const std::vector<double>& rho) -> void;

 private:
  std::filesystem::path outputDir;
  bool clearDir;
};