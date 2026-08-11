/*
    sqlite.cc - Optional SQLite storage for the tabular export
    (control_print_tabular). Pure C++ (RAII), guarded by SQLITE_USE.

    Schema (normalised by (node,t) composite key):
      primary(node, t, ux.., sigxx.., exx..)     -- primary variables
      derived(node, t, vmises, tresca, sig1..3)  -- derived magnitudes (C++)
      user(node, t, ...)                         -- user variables (Python)
      meta(key, value)                           -- mesh, units, convention

    If SQLITE_USE is 0 (or sqlite3.h is unavailable), the class methods are
    stubs that set `available=false`; the caller (print_tabular) warns and
    falls back to CSV.
*/

#include "tochnog.h"
#include <string>

#if SQLITE_USE
#include <sqlite3.h>
#endif

#if SQLITE_USE
// RAII wrapper around the SQLite C API.
class SqliteDB {
public:
  bool available;
  std::string filename;

  SqliteDB() : available(false), db_(nullptr) {}
  ~SqliteDB() { close(); }

  bool open(const char* name) {
    filename = name;
    if ( sqlite3_open(name, &db_) != SQLITE_OK ) {
      pri( "Error: cannot open SQLite database", name );
      close();
      return false;
    }
    available = true;
    create_schema();
    return true;
  }

  bool exec(const char* sql) {
    if ( !db_ ) return false;
    char* errmsg = nullptr;
    if ( sqlite3_exec(db_, sql, nullptr, nullptr, &errmsg) != SQLITE_OK ) {
      if ( errmsg ) {
        pri( "SQLite error: ", errmsg );
        sqlite3_free(errmsg);
      }
      return false;
    }
    return true;
  }

  void close() {
    if ( db_ ) { sqlite3_close(db_); db_ = nullptr; }
    available = false;
  }

private:
  sqlite3* db_;

  void create_schema() {
    exec( "CREATE TABLE IF NOT EXISTS meta (key TEXT PRIMARY KEY, value TEXT);" );
    exec( "CREATE TABLE IF NOT EXISTS primary_data ("
          " node INTEGER, t REAL,"
          " ux REAL, uy REAL, uz REAL,"
          " sigxx REAL, sigyy REAL, sigzz REAL, sigxy REAL, sigyz REAL, sigzx REAL,"
          " exx REAL, eyy REAL, ezz REAL, exy REAL, eyz REAL, ezx REAL,"
          " PRIMARY KEY (node, t));" );
    exec( "CREATE TABLE IF NOT EXISTS derived ("
          " node INTEGER, t REAL,"
          " vmises REAL, tresca REAL, sig1 REAL, sig2 REAL, sig3 REAL,"
          " PRIMARY KEY (node, t));" );
    exec( "CREATE TABLE IF NOT EXISTS user_data ("
          " node INTEGER, t REAL,"
          " PRIMARY KEY (node, t));" );
  }
};

#else
// Stub: no SQLite support in this build.
class SqliteDB {
public:
  bool available;
  std::string filename;
  SqliteDB() : available(false) {}
  ~SqliteDB() {}
  bool open(const char* name) { (void)name; available = false; return false; }
  bool exec(const char* sql) { (void)sql; return false; }
  void close() { available = false; }
};
#endif

// Entry point used by print_tabular: create/open the database and return
// the RAII object by pointer (the caller owns it).
SqliteDB* sqlite_db_open( const char* filename )
{
  SqliteDB* db = new SqliteDB();
  if ( !db->open(filename) ) {
    delete db;
    return nullptr;
  }
  return db;
}

void sqlite_db_close( SqliteDB* db )
{
  if ( db ) delete db;
}
