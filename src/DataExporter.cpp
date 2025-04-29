#include "DataExporter.h"

DataExporter::DataExporter(const std::string& outputDir, bool clearDir) 
  : outputDir(std::filesystem::path(outputDir)), clearDir(clearDir) {
  if (!std::filesystem::exists(outputDir)) {
    std::filesystem::create_directories(outputDir);
  } else if (clearDir) {
    for (const auto& entry : std::filesystem::directory_iterator(outputDir)) {
      std::filesystem::remove_all(entry.path()); // 删除目录内所有文件和子目录
    }
  }
}

auto DataExporter::writeTimeStep(size_t timestep, const std::vector<double>& p, const std::vector<double>& u,
  const std::vector<double>& rho) -> void {
  if (p.size() != u.size() || u.size() != rho.size()) {
    throw std::invalid_argument("向量长度不一致！");
  }
  auto filename = std::format("{}/frame_{:03d}.csv", outputDir.string(), timestep);
  std::ofstream file(filename, std::ios::trunc);
  if (!file) {
    throw std::runtime_error("无法写入文件: " + filename);
  }
  for (size_t i = 0; i < u.size(); ++i) {
    file << std::format("{},{},{}\n", p[i], u[i], rho[i]);
}
}

auto DataExporter::writeTimeStep2D(size_t timestep, const std::vector<double>& p, const std::vector<double>& u, const std::vector<double>& v,
    const std::vector<double>& rho) -> void {
  if (p.size() != u.size() || u.size() != rho.size()) {
    throw std::invalid_argument("向量长度不一致！");
  }
  auto filename = std::format("{}/frame_{:03d}.csv", outputDir.string(), timestep);
  std::ofstream file(filename, std::ios::trunc);
  if (!file) {
    throw std::runtime_error("无法写入文件: " + filename);
  }
  for (size_t i = 0; i < p.size(); ++i) {
    file << std::format("{},{},{},{}\n",  p[i], u[i], v[i], rho[i]);
  } 
}