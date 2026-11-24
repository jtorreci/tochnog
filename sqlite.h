/*
    sqlite.h - Optional SQLite storage for the tabular export
    (control_print_tabular). Pure C++ (RAII), guarded by SQLITE_USE.

    Schema (normalised by the (node,t) composite key):
      primary(node, dof, t, value)              -- primary variables (long)
      derived(node, t, vmises, tresca, sig1..3) -- derived magnitudes (C++)
      user(node, t, ...)                        -- user variables (Python)
      meta(key, value)                          -- mesh, units, convention

    The long-format primary table lets any dof set be stored without ALTER
    TABLE; pandas/JOIN by (node,t) is used for analysis.

    If SQLITE_USE is 0 (or sqlite3.h is unavailable), the class methods are
    stubs that set `available=false`; the caller (print_tabular) warns and
    falls back to CSV.
*/

#ifndef TN_SQLITE_H
#define TN_SQLITE_H

#include <string>

#if SQLITE_USE
#include <sqlite3.h>
#endif

// RAII wrapper around the SQLite C API.
class SqliteDB {
public:
  bool available;
  std::string filename;

  SqliteDB();
  ~SqliteDB();
  bool open(const char* name);
  bool exec(const char* sql);
  void close();

private:
#if SQLITE_USE
  sqlite3* db_;
#else
  int db_; // unused placeholder
#endif
  void create_schema();
};

// Entry points used by print_tabular.
SqliteDB* sqlite_db_open( const char* filename );
void      sqlite_db_close( SqliteDB* db );

#endif // TN_SQLITE_H
