import sqlite3
import shelve
import pathlib

from enum import IntEnum

import idseq_dag.util.log as log

DICT_DELIMITER = chr(1) # Delimiter for an array of values
class IdSeqDictValue(IntEnum):
    VALUE_TYPE_SCALAR = 1
    VALUE_TYPE_ARRAY = 2
SQLITE_TABLE_NAME = "idseq_dict"

class IdSeqDict(object):
    '''
        a key -> value permanent dictionary store with sqlite3 as the backend.
        key has to be a string
        value could be either array of strings or a string
    '''
    def __init__(self, db_path, value_type, read_only=True):
        self.db_path = db_path
        self.value_type = value_type
        self.read_only = read_only
        self.uri_db_path = pathlib.Path(self.db_path).as_uri()
        if self.read_only:
            self.uri_db_path += "?mode=ro&immutable=1&nolock=1&cache=private"
        with log.log_context(f"db_open", {"db_path": self.db_path, "read_only": self.read_only, "uri_db_path": self.uri_db_path}):
            self.db_conn = sqlite3.connect(self.uri_db_path, check_same_thread=False, uri=True) # make it thread safe
        self._open = True
        with log.log_context(f"db_assert_table", {"db_path": self.db_path, "read_only": self.read_only}):
            cursor = self.db_conn.cursor()
            if not self.read_only:
                res = cursor.execute(f"SELECT count(*) FROM sqlite_master WHERE type='table' AND name='{SQLITE_TABLE_NAME}';")
                v = res.fetchone()
                if v[0] == 0:
                    raise Exception(f"table {SQLITE_TABLE_NAME} doesn't exist in db {self.db_path}")
            else:
                cursor.execute(f"CREATE TABLE IF NOT EXISTS {SQLITE_TABLE_NAME} (dict_key VARCHAR(255) PRIMARY KEY, dict_value text)")
                self.db_conn.commit()
            cursor.close()

    def __enter__(self):
        ''' context manager enter '''
        return self

    def __exit__(self, type, value, traceback):
        ''' context manager exit '''
        self.close()
        return False

    def __del__(self):
        ''' destructor '''
        self.close()

    def close(self):
        ''' close db if it is open. If the db is not read-only, it also performs a commit before closing it. '''
        if self._open:
            if not self.read_only:
                with log.log_context("db_commit", {"db_path": self.db_path, "read_only": self.read_only}):
                    self.db_conn.commit()
            with log.log_context("db_close", {"db_path": self.db_path, "read_only": self.read_only}):
                self.db_conn.close()
            self._open = False

    def update(self, key, value):
        ''' Update a particular key value pair '''
        cursor = self.db_conn.cursor()
        val = value
        if self.value_type == IdSeqDictValue.VALUE_TYPE_ARRAY:
            val = DICT_DELIMITER.join([str(v) for v in value])
        cursor.execute(f"INSERT OR REPLACE INTO {SQLITE_TABLE_NAME} VALUES ('{key}', '{val}')")
        cursor.close()

    def batch_inserts(self, tuples):
        ''' Insert multiple records at a time '''
        if self.read_only:
            raise Exception(f"db {self.db_path} is open in read-only mode")
        if len(tuples) == 0:
            return
        def tuple_to_sql_str(user_tuple):
            val = user_tuple[1]
            if self.value_type == IdSeqDictValue.VALUE_TYPE_ARRAY:
                val = DICT_DELIMITER.join([str(v) for v in val])
            return (user_tuple[0], val)

        # Use a parameterized SQL string, so that special characters such as quotes are automatically escaped.
        value_arr = []
        for pair in tuples:
            values = tuple_to_sql_str(pair)
            value_arr += values

        parameter_str = ",".join(["(?,?)"] * len(tuples))

        cursor = self.db_conn.cursor()
        cursor.execute(f"INSERT OR REPLACE INTO {SQLITE_TABLE_NAME} VALUES {parameter_str}", value_arr)
        self.db_conn.commit()
        cursor.close()

    def get(self, key, default_value=None):
        ''' Emulate get as a python dictionary  '''
        cursor = self.db_conn.cursor()
        res = cursor.execute(f"SELECT dict_value FROM idseq_dict where dict_key = '{key}'")
        v = res.fetchone()
        if v is None:
            return default_value
        value = v[0]
        if self.value_type == IdSeqDictValue.VALUE_TYPE_ARRAY:
            return value.split(DICT_DELIMITER)
        cursor.close()
        return value

def open_file_db_by_extension(db_path, value_type=IdSeqDictValue.VALUE_TYPE_SCALAR, read_only=True):
    ''' if extension is .sqlite3 open it as an IdSeqDict, otherwise, open as shelve in read mode '''
    if db_path[-8:] == '.sqlite3': # sqlite3 format
        return IdSeqDict(db_path, value_type, read_only)
    # shelve format
    return shelve.open(db_path.replace('.db', ''), 'r')
