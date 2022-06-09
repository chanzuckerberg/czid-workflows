import shelve

def open_file_db_by_extension(db_path):
    return shelve.open(db_path, 'r')
