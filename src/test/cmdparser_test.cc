#include<iostream>

#include"../utility/cmdparser.h"

using namespace std;

void configure_parser(cli::Parser& parser) {
  parser.set_optional<std::string>("o", "output", "data", "Strings are naturally included.");
  parser.set_optional<int>("n", "number", 8, "Integers in all forms, e.g., unsigned int, long long, ..., are possible. Hexadecimal and Ocatl numbers parsed as well");
  parser.set_optional<cli::NumericalBase<int, 10>>("t", "temp", 0, "integer parsing restricted only to numerical base 10");
  parser.set_optional<double>("b", "beta", 11.0, "Also floating point values are possible.");
  parser.set_optional<bool>("a", "all", false, "Boolean arguments are simply switched when encountered, i.e. false to true if provided.");
  parser.set_required<std::vector<std::string>>("v", "values", "By using a vector it is possible to receive a multitude of inputs.");
  parser.set_required<std::vector<std::string>>("x", "xs", "By using a vector it is possible to receive a multitude of inputs.");
}

int main(int argc, char *argv[]) {
  //for (int i = 0; i < argc; ++i) {
  //  cout << i << ": "<< argv[i] << endl;
  //}
  cli::Parser parser(argc, argv);
  configure_parser(parser);
  parser.run_and_exit_if_error();

  auto number = parser.get<int>("n");
  //auto will be std::string
  auto output = parser.get<std::string>("o");

  //auto will be bool
  auto all = parser.get<bool>("a");

  //auto will be std::vector<short>
  auto values = parser.get<std::vector<std::string>>("v");

  auto xs = parser.get<std::vector<std::string>>("x");

  cout << number << endl;
  cout << output << endl;
  cout << all << endl;
  cout << "print out v, size is  " << values.size() << endl;
  for (auto v : values) {
    cout << "\t" << v << endl;
  }

  cout << "print out x, size is  " << xs.size() << endl;
  for (auto x : xs) {
    cout << "\t" << x << endl;
  }

  return 0;
}
