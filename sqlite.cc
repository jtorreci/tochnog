/*
    sqlite.cc - Implementation of SqliteDB (see sqlite.h).

    Schema (normalised by (node,t) composite key):
      primary(node, t, ux.., sigxx.., exx..)     -- primary variables
      derived(node, t, vmises, tresca, sig1..3)  -- derived magnitudes (C++)
      user(node, t, ...)                         -- user variables (Python)
      meta(key, value)                           -- mesh, units, convention
*/

#include "tochnog.h"
#include "sqlite.h"

#if SQLITE_USE
SqliteDB::SqliteDB() : available(false), db_(nullptr) {}
SqliteDB::~SqliteDB() { close(); }

bool SqliteDB::open(const char* name) {
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

bool SqliteDB::exec(const char* sql) {
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

void SqliteDB::close() {
  if ( db_ ) { sqlite3_close(db_); db_ = nullptr; }
  available = false;
}

void SqliteDB::create_schema() {
  exec( "CREATE TABLE IF NOT EXISTS meta (key TEXT PRIMARY KEY, value TEXT);" );
  exec( "CREATE TABLE IF NOT EXISTS primary_data ("
        " node INTEGER, dof TEXT, t REAL, value REAL,"
        " PRIMARY KEY (node, dof, t));" );
  exec( "CREATE TABLE IF NOT EXISTS derived ("
        " node INTEGER, t REAL,"
        " vmises REAL, tresca REAL, sig1 REAL, sig2 REAL, sig3 REAL,"
        " PRIMARY KEY (node, t));" );
  exec( "CREATE TABLE IF NOT EXISTS user_data ("
        " node INTEGER, t REAL,"
        " PRIMARY KEY (node, t));" );
}

#else
// Stub: no SQLite support in this build.
SqliteDB::SqliteDB() : available(false), db_(0) {}
SqliteDB::~SqliteDB() {}
bool SqliteDB::open(const char* name) { (void)name; available = false; return false; }
bool SqliteDB::exec(const char* sql) { (void)sql; return false; }
void SqliteDB::close() { available = false; }
void SqliteDB::create_schema() {}
#endif

// Entry points.
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
