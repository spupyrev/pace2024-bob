#pragma once

#include <iostream>
#include <string>
#include <chrono>
#include <ctime>
#include <cstdarg>

// assertion (wild server error)
const int ERROR_EXIT_CODE = 99;
// application error
const int CHECK_EXIT_CODE = 40;

#define stringize(s) #s
#define XSTR(s) stringize(s)

#ifdef DEBUG

#define CHECK1(condition) \
  if (0 == (condition)) { \
    CHECK_LOG_EMPTY(XSTR(condition), __FILE__, __LINE__); \
    throw CHECK_EXIT_CODE; \
  }
#define CHECK2(condition, ...) \
  if (0 == (condition)) { \
    CHECK_LOG(XSTR(condition), __FILE__, __LINE__, __VA_ARGS__); \
    throw CHECK_EXIT_CODE; \
  }
#else

#define CHECK1(condition) {}
#define CHECK2(condition, ...) {}

#endif

#define GET_MACRO(_1,_2,_3,_4,_5,_6,_7,_8,NAME,...) NAME
  // TODO
#define CHECK(...) GET_MACRO(__VA_ARGS__, CHECK2, CHECK2, CHECK2, CHECK2, CHECK2, CHECK2, CHECK2, CHECK1)(__VA_ARGS__)

#define ERROR(message) \
  { \
    std::cerr << "\033[91m" << message << "\033[0m" << " [" << __FILE__ << ":" << __LINE__ << "]\n"; \
    throw ERROR_EXIT_CODE; \
  }

enum class TextColor {
  none,
  red,
  blue,
  green
};

inline void CHECK_LOG(const char* desc, const char* file, int line, const char* message, va_list args) {
#ifdef DEBUG
  char* buffer = new char[16384];
  setlocale(LC_NUMERIC, "");
  std::vsprintf(buffer, message, args);
  std::cerr << "\033[91m" << "assertion '" << desc << "' failed: " << std::string(buffer) << "\033[0m" << " [" << file << ":" << line << "]\n";
  delete[] buffer;
#endif
}

inline void CHECK_LOG(const char* desc, const char* file, int line, ...) {
#ifdef DEBUG
  va_list args;
  va_start(args, line);
  const char* message = va_arg(args, const char*);
  CHECK_LOG(desc, file, line, message, args);
  va_end(args);
#endif
}

inline void CHECK_LOG_EMPTY(const char* desc, const char* file, int line) {
#ifdef DEBUG
  std::cerr << "\033[91m" << "assertion '" << desc << "' failed" << "\033[0m" << " [" << file << ":" << line << "]\n";
#endif
}

inline void LOG(TextColor color, const char* message, va_list args) {
#ifdef DEBUG
  auto end = std::chrono::system_clock::now();
  auto time = std::chrono::system_clock::to_time_t(end);
  auto stime = std::string(std::ctime(&time));
  stime.erase(stime.find_last_not_of(" \n\r\t") + 1);
  stime = "\033[90m" + stime + ":\033[0m";
  char* buffer = new char[16384];
  setlocale(LC_NUMERIC, "");
  std::vsprintf(buffer, message, args);
  std::cerr << stime << " ";

  switch (color) {
    case TextColor::none: std::cerr << std::string(buffer); break;

    case TextColor::red : std::cerr << "\033[91m" << std::string(buffer) << "\033[0m"; break;

    case TextColor::blue: std::cerr << "\e[38;5;12m" << std::string(buffer) << "\e[0m"; break;

    case TextColor::green: std::cerr << "\e[38;5;34m" << std::string(buffer) << "\e[0m"; break;
  };

  std::cerr << "\n";
  delete[] buffer;
#endif
}

inline void LOG(const std::string& message) {
#ifdef DEBUG
  va_list args;
  LOG(TextColor::none, message.c_str(), args);
#endif
}

inline void LOG(const char* message, ...) {
#ifdef DEBUG
  va_list args;
  va_start(args, message);
  LOG(TextColor::none, message, args);
  va_end(args);
#endif
}

inline void LOG(TextColor color, const char* message, ...) {
#ifdef DEBUG
  va_list args;
  va_start(args, message);
  LOG(color, message, args);
  va_end(args);
#endif
}

inline void LOG_IF(TextColor color, bool condition, const char* message, ...) {
#ifdef DEBUG
  if (!condition)
    return;
  
  va_list args;
  va_start(args, message);
  LOG(color, message, args);
  va_end(args);
#endif
}

inline void LOG_IF(bool condition, const char* message, ...) {
#ifdef DEBUG
  if (!condition)
    return;

  va_list args;
  va_start(args, message);
  LOG(TextColor::none, message, args);
  va_end(args);
#endif
}
