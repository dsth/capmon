#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H

#include <stdexcept>

class CapServerError : public std::runtime_error {
  public:
    CapServerError(const std::string& msg = "") : std::runtime_error(msg) {}
};

class JsonError : public std::runtime_error {
  public:
    JsonError(const std::string& msg = "") : std::runtime_error(msg) {}
};

class MySqlError : public std::runtime_error {
  public:
    MySqlError(const std::string& msg = "") : std::runtime_error(msg) {}
};

class Sqlite3Error : public std::runtime_error {
  public:
    Sqlite3Error(const std::string& msg = "") : std::runtime_error(msg) {}
};

class FileTypeError : public std::runtime_error {
  public:
    FileTypeError(const std::string& msg = "") : std::runtime_error(msg) {}
};

class ExternalExecutable : public std::runtime_error {
  public:
    ExternalExecutable(const std::string& msg = "") : std::runtime_error(msg) {}
};

class MySqlConnError : public MySqlError {
  public:
    MySqlConnError(const std::string& msg = "") : MySqlError(msg) {}
};

class MyError {
    const char* const data;
  public:
    MyError(const char* const msg = 0) : data (msg) {}
};

#endif

