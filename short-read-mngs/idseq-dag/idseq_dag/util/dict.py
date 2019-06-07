import sqlite3
import shelve
import pathlib
import multiprocessing

from enum import IntEnum

import idseq_dag.util.log as log

DICT_DELIMITER = chr(1) # Delimiter for an array of values
class IdSeqDictValue(IntEnum):
    VALUE_TYPE_SCALAR = 1
    VALUE_TYPE_ARRAY = 2
SQLITE_TABLE_NAME = "idseq_dict"


class _IdSeqDictBase(object):
    '''
        Private superclass of IdSeqDict and IdSeqDictForUpdate.
        See IdSeqDict for more information.
    '''

    def __init__(self, db_path, value_type):
        self.db_path = db_path
        self.value_type = value_type
        # Private from the point of view of users.  Accessible to derived classes.
        self._db_conn = None
        self._lock = multiprocessing.RLock()
        self._lock.acquire()

    def _assert_lock_held(self):
        assert self._lock.acquire(block=False), "Sharing an IdSeqDict object among multiple threads or processes is not supported."

    def _uri_base(self):
        return pathlib.Path(self.db_path).as_uri()

    def _connect(self):
        self._assert_lock_held()
        uri_db_path = self._uri_base()
        with log.log_context(f"db_open", {"db_path": self.db_path, "uri_db_path": uri_db_path}):
            return sqlite3.connect(uri_db_path, uri=True)

    def _ensure_table_exists(self, _conn):
        assert False and self, "Abstract method.  Behavior is different for read-only vs read-write mode."

    def __enter__(self):
        ''' Open a connection, ensure table exists, and enter connection context (i.e. begin transaction). '''
        # Connection objects can be used as context managers that automatically commit or rollback transactions.
        # In the event of an exception, the transaction is rolled back; otherwise, the transaction is committed.
        # https://docs.python.org/2/library/sqlite3.html#using-the-connection-as-a-context-manager
        self._assert_lock_held()
        assert not self._db_conn
        conn = self._connect()
        try:
            # The connection is first entered and exited in ensure_table_exists(), then entered again
            # to start a new transaction for the context block being managed.
            self._ensure_table_exists(conn)
            conn.__enter__()
        except:
            # Try not to leak open connections.  If an exception is raised above,
            # the self.__exit__() method will never run, so we have to close
            # the connection here.
            conn.close()
            raise
        # Leaving this last ensures we don't __exit__ a conn that hasn't been __enter__'ed.
        self._db_conn = conn
        return self

    def __exit__(self, etype, evalue, etraceback):
        ''' If an exception has occurred, the connection will auto-rollback;  otherwise it will auto-commit. '''
        self._assert_lock_held()  # It would be extremely unusual for this to fail.
        conn = self._db_conn
        self._db_conn = None
        try:
            return conn.__exit__(etype, evalue, etraceback)
        finally:
            conn.close()


class IdSeqDict(_IdSeqDictBase):
    '''
        A key -> value permanent dictionary store with sqlite3 as the backend.

        The key has to be a string.

        The value could be either an array of strings or a string.

        Use only through a WITH context, like so:

           with IdSeqDict("somefile.sql3") as sf:
               print(sf.get("somekey"))

        This opens a read-only "connection" to the database file, and looks up
        the value for "somekey".  The connection is closed when the block exits.

        For write access, use IdSeqDictForUpdate.  That begins a transaction
        which is either committed or rolled back at block exit depending on
        whether exceptions have occurred in the block, as per sqlite3 doc.

        To open the same file from multiple threads or processes, each thread or
        process must execute its own WITH block with an IdSeqDict() constructor
        exactly as above.

        Accessing the same IdSeqDict instance/connection from multiple threads/
        processes is a checked runtime error.

        ILLEGAL ACCESSS

        Accessing "sf" outside of the WITH block will raise an assert.

        Accessing "sf" from a child thread or process will raise an assert.

        "Accessing" means calling get(), update(), batch_inserts(), or, trying
        to enter/exit as a context via "with".
    '''

    def __init__(self, db_path, value_type):
        super().__init__(db_path, value_type)

    def _get_uri(self):
        ''' Private: URI for readonly access '''
        return self._uri_base() + "?mode=ro&immutable=1&nolock=1&cache=private"

    def _ensure_table_exists(self, conn):
        ''' Private: Fail if the table does not exist.  Called when self._db_conn is still none. '''
        self._assert_lock_held()
        with log.log_context(f"db_assert_table", {"db_path": self.db_path}):
            with conn:
                res = conn.execute(f"SELECT count(*) FROM sqlite_master WHERE type='table' AND name='{SQLITE_TABLE_NAME}';")
                table_exists = res.fetchone()[0] != 0
            assert table_exists, f"table {SQLITE_TABLE_NAME} doesn't exist in db {self.db_path}"

    @staticmethod
    def _table_exists(conn, table_name):
        ''' Private: Return True just when table_name exists. '''
        with conn:
            res = conn.execute(f"SELECT count(*) FROM sqlite_master WHERE type='table' AND name='{table_name}';")
            return res.fetchone()[0] != 0

    def get(self, key, default_value=None):
        ''' Emulate get as a python dictionary  '''
        self._assert_lock_held()
        result = default_value
        v = self._db_conn.execute(f"SELECT dict_value FROM idseq_dict where dict_key = '{key}'").fetchone()
        if v:
            result = v[0]
            if self.value_type == IdSeqDictValue.VALUE_TYPE_ARRAY:
                result = result.split(DICT_DELIMITER)
        return result


class IdSeqDictForUpdate(IdSeqDict):
    '''
        A version of IdSeqDict that permits updates.  See IdSeqDict for more information.
    '''

    def __init__(self, db_path, value_type):
        super().__init__(db_path, value_type)

    def _get_uri(self):
        ''' URI for read and write access '''
        return self._uri_base()

    def _ensure_table_exists(self, conn):
        ''' Create writable table if one doesn't exist. '''
        self._assert_lock_held()
        with log.log_context(f"db_assert_table", {"db_path": self.db_path}):
            with conn:
                conn.execute(f"CREATE TABLE IF NOT EXISTS {SQLITE_TABLE_NAME} (dict_key VARCHAR(255) PRIMARY KEY, dict_value text)")

    def update(self, key, value):
        ''' Update a particular key value pair.  Will commit or roll back at exit of WITH block. '''
        self._assert_lock_held()
        val = value
        if self.value_type == IdSeqDictValue.VALUE_TYPE_ARRAY:
            val = DICT_DELIMITER.join(str(v) for v in value)
        self._db_conn.execute(f"INSERT OR REPLACE INTO {SQLITE_TABLE_NAME} VALUES ('{key}', '{val}')")

    def batch_inserts(self, tuples):
        ''' Insert multiple records at a time.  Will commit or roll back at exit of WITH block. '''
        self._assert_lock_held()

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

        # This adds to the pending transaction.
        self._db_conn.execute(f"INSERT OR REPLACE INTO {SQLITE_TABLE_NAME} VALUES {parameter_str}", value_arr)


def open_file_db_by_extension(db_path, value_type=IdSeqDictValue.VALUE_TYPE_SCALAR):
    ''' if extension is .sqlite3 open it as an IdSeqDict, otherwise, open as shelve in read mode '''
    if db_path[-8:] == '.sqlite3': # sqlite3 format
        return IdSeqDict(db_path, value_type)
    # shelve format;  shelve objects can be used as context managers
    return shelve.open(db_path.replace('.db', ''), 'r')
