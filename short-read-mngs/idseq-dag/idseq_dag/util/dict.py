import sqlite3
import shelve
from enum import IntEnum

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
    def __init__(self, db_path, value_type):
        self.db_path = db_path
        self.value_type = value_type
        self.db_conn = sqlite3.connect(db_path, check_same_thread=False) # make it thread safe
        cursor = self.db_conn.cursor()
        cursor.execute(f"CREATE TABLE IF NOT EXISTS {SQLITE_TABLE_NAME} (dict_key VARCHAR(255) PRIMARY KEY, dict_value text)")
        self.db_conn.commit()

    def __del__(self):
        ''' destructor '''
        self.db_conn.commit()
        self.db_conn.close()

    def update(self, key, value):
        ''' Update a particular key value pair '''
        cursor = self.db_conn.cursor()
        val = value
        if self.value_type == IdSeqDictValue.VALUE_TYPE_ARRAY:
            val = DICT_DELIMITER.join([str(v) for v in value])
        cursor.execute(f"INSERT OR REPLACE INTO {SQLITE_TABLE_NAME} VALUES ('{key}', '{val}')")

    def batch_inserts(self, tuples):
        ''' Insert multiple records at a time '''
        if len(tuples) == 0:
            return
        def tuple_to_sql_str(user_tuple):
            val = user_tuple[1]
            if self.value_type == IdSeqDictValue.VALUE_TYPE_ARRAY:
                val = DICT_DELIMITER.join([str(v) for v in val])
            return f"('{user_tuple[0]}', '{val}')"
        value_str = ",".join(map(tuple_to_sql_str, tuples))
        cursor = self.db_conn.cursor()
        cursor.execute(f"INSERT OR REPLACE INTO {SQLITE_TABLE_NAME} VALUES {value_str}")
        self.db_conn.commit()

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
        return value

def open_file_db_by_extension(db_path, value_type=IdSeqDictValue.VALUE_TYPE_SCALAR):
    ''' if extension is .sqlite3 open it as an IdSeqDict, otherwise, open as shelve in read mode '''
    if db_path[-8:] == '.sqlite3': # sqlite3 format
        return IdSeqDict(db_path, value_type)
    # shelve format
    return shelve.open(db_path.replace('.db', ''), 'r')
