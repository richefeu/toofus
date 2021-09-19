#include "../anova.hpp"

int main(int argc, char const *argv[]) {
  // for comparaison:
  // https://www.statskingdom.com/180Anova1way.html

  std::vector<std::vector<double>> samples;

  std::vector<double> group1;
  group1.push_back(-7.75405);
  group1.push_back(-7.70286);
  group1.push_back(-7.68725);
  group1.push_back(-7.61047);
  group1.push_back(-7.60942);
  samples.push_back(group1);

  std::vector<double> group2;
  group2.push_back(-7.35701);
  group2.push_back(-7.29485);
  group2.push_back(-7.28961);
  group2.push_back(-7.26047);
  samples.push_back(group2);

  anovaResult result;
  anova::test(samples, result);
  anova::print(result);

  return 0;
}