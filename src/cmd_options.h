#pragma once

#include "common.h"

#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <memory>
#include <cassert>
#include <regex>

class CMDOptions {
 private:
  std::string usageMessage;
  std::map<std::string, std::string> options;
  size_t verbose = 0;

  std::map<std::string, std::string> allowedOptions;
  std::vector<std::string> allowedOptionsOrder;
  std::map<std::string, std::string> defaultValues;
  std::map<std::string, std::vector<std::string> > allowedValues;

  CMDOptions(const CMDOptions&);
  CMDOptions& operator = (const CMDOptions&);
  CMDOptions() {}

 public:
  static std::unique_ptr<CMDOptions> Create() {
    return std::unique_ptr<CMDOptions>(new CMDOptions());
  }

  size_t ParseSafe(int argc, char** argv, const std::string& optionName, const size_t defaultValue) {
    for (int i = 1; i < argc; i++) {
      std::string s(argv[i]);
      size_t equalIndex = s.find('=');
      std::string name = s.substr(0, equalIndex);
      if (name == optionName) {
        std::string value = (equalIndex == std::string::npos ? "" : s.substr(equalIndex + 1));
        return to_sizet(value);
      }
    }
    return defaultValue;
  }

  void Parse(int argc, char** argv) {
    for (int i = 1; i < argc; i++) {
      std::string s(argv[i]);

      if (verbose) {
        if (s == "/?"  || s == "-?" || s == "--help" || s == "-help") {
          Usage(argv[0]);
          throw 0;
        }
      }

      SetOption(s);
    }
  }

  void SetUsageMessage(const std::string& msg) {
    usageMessage = msg;
  }

  void SetVerbose(size_t verbose_) {
    verbose = verbose_;
  }

  void SetOption(const std::string& s) {
    size_t equalIndex = s.find('=');
    std::string name = s.substr(0, equalIndex);

    if (!allowedOptions.count(name)) {
      if (equalIndex == std::string::npos && allowedOptions.count("")) {
        options[""] = name;
        return;
      }

      if (!verbose)
        return;
      UnrecognizedOption(name);
    }

    std::string value = (equalIndex == std::string::npos ? "" : s.substr(equalIndex + 1));

    if (!options.count(name) || (defaultValues.count(name) && options[name] == defaultValues[name])) {
      options[name] = value;
    }

    if (!allowedValues[name].empty()) {
      // check regexp
      bool found = false;
      for (std::string& allowedValue : allowedValues[name]) {
        if (allowedValue == value) {
          found = true;
          break;
        }
        // check regexp
        if (std::regex_match(value, std::regex(allowedValue))) {
          found = true;
          break;
        }
      }
      if (!found) {
        InvalidOption(name);
      }
    }
  }

  void AddAllowedOption(const std::string& optionName, const std::string& defaultValue, const std::string& description) {
    AddAllowedOption(optionName, description);
    options[optionName] = defaultValue;
    defaultValues[optionName] = defaultValue;
  }

  void AddAllowedOption(const std::string& optionName, const std::string& description) {
    assert(!allowedOptions.count(optionName));
    allowedOptions[optionName] = description;
    allowedOptionsOrder.push_back(optionName);
  }

  void AddAllowedValue(const std::string& optionName, const std::string& value) {
    assert(allowedOptions.count(optionName));
    allowedValues[optionName].push_back(value);
  }

  std::string getOption(const std::string& optionName) const {
    if (!options.count(optionName)) {
      if (allowedOptions.count(optionName)) {
        UnspecifiedOption(optionName);
      }

      UnrecognizedOption(optionName);
    }

    assert(options.count(optionName));
    return (*options.find(optionName)).second;
  }

  void setOption(const std::string& optionName, const std::string& value) {
    options[optionName] = value;
  }

  void set(const std::string& optionName, const std::string& value) {
    options[optionName] = value;
  }

  void setStr(const std::string& optionName, const std::string& value) {
    options[optionName] = value;
  }

  std::string getStr(const std::string& optionName) const {
    return getOption(optionName);
  }

  std::string get(const std::string& optionName) const {
    return getOption(optionName);
  }

  int getInt(const std::string& optionName) const {
    return to_int(getOption(optionName));
  }

  double getDouble(const std::string& optionName) const {
    return to_double(getOption(optionName));
  }

  void setInt(const std::string& optionName, int value) {
    if (!options.count(optionName)) {
      UnrecognizedOption(optionName);
    }

    assert(options.count(optionName));
    setOption(optionName, to_string(value));
  }

  bool getBool(const std::string& optionName) const {
    return getOption(optionName) != "false" && getOption(optionName) != "0";
  }

  void setBool(const std::string& optionName, bool value) {
    if (!options.count(optionName)) {
      UnrecognizedOption(optionName);
    }

    assert(options.count(optionName));
    setOption(optionName, value ? "true" : "false");
  }

  bool hasOption(const std::string& optionName) const {
    if (!allowedOptions.count(optionName)) {
      UnrecognizedOption(optionName);
    }

    return options.count(optionName) > 0;
  }

  void UnspecifiedOption(const std::string& optionName) const {
    if (verbose)
      std::cout << "required option \"" << optionName << "\" is not specified\n";
    throw 1;
  }

  void UnrecognizedOption(const std::string& optionName) const {
    if (verbose)
      std::cout << "unrecognized option \"" << optionName << "\"\n";
    throw 1;
  }

  void InvalidOption(const std::string& optionName) const {
    if (verbose)
      std::cout << "value \"" << getOption(optionName) << "\" is invalid for option \"" << optionName << "\"\n";
    throw 1;
  }

  void Usage(const std::string& program) const {
    if (usageMessage != "") {
      std::cout << usageMessage << "\n";
    } else {
      std::cout << "Usage: " << program << " [options]" << "\n";
    }

    std::cout << "Allowed options:";

    for (auto opt : allowedOptionsOrder) {
      std::string name = allowedOptions.find(opt)->first;

      if (name.length() == 0) {
        continue;
      }

      std::cout << "\n";
      std::cout << "  " << name;

      if (allowedValues.count(name)) {
        auto av = allowedValues.find(name)->second;

        if (!av.empty()) {
          std::cout << "=";
          bool first = true;

          for (std::string s : av)
            if (first)
            { std::cout << "[" << s; first = false; }
            else {
              std::cout << "|" << s;
            }

          std::cout << "]";
        }
      }

      std::cout << "\n";
      std::cout << "  " << allowedOptions.find(opt)->second << "\n";
    }
  }
};
