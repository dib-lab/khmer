import DBConstants
import os
import sqlite3
import itertools

def create_db(filepath, fields, rcrditer):
    """
    Creates a screed database in the given filepath. Fields is a tuple
    specifying the names and relative order of attributes in a
    record. rcrditer is an iterator returning records over a
    sequence dataset. Records yielded are in dictionary form
    """
    if not filepath.endswith(DBConstants.fileExtension):
        filepath += DBConstants.fileExtension

    if os.path.exists(filepath): # Remove existing files
        os.unlink(filepath)

    con = sqlite3.connect(filepath)
    cur = con.cursor()

    # Sqlite PRAGMA settings for speed
    cur.execute("PRAGMA synchronous='OFF'")
    cur.execute("PRAGMA locking_mode=EXCLUSIVE")

    # Create the admin table
    cur.execute('CREATE TABLE %s (%s INTEGER PRIMARY KEY, '\
                '%s TEXT, %s TEXT)' % (DBConstants._SCREEDADMIN,
                                       DBConstants._PRIMARY_KEY,
                                       DBConstants._FIELDNAME,
                                       DBConstants._ROLENAME))
    query = 'INSERT INTO %s (%s, %s) VALUES (?, ?)' % \
            (DBConstants._SCREEDADMIN, DBConstants._FIELDNAME,
             DBConstants._ROLENAME)

    # Put the primary key in as an attribute
    cur.execute(query, (DBConstants._PRIMARY_KEY,
                        DBConstants._PRIMARY_KEY_ROLE))
    for attribute, role in fields:
        cur.execute(query, (attribute, role))

    # Setup the dictionary table creation field substring
    fieldsub = ','.join(['%s TEXT' % field for field, role in fields])

    # Create the dictionary table
    cur.execute('CREATE TABLE %s (%s INTEGER PRIMARY KEY, %s)' %
                (DBConstants._DICT_TABLE, DBConstants._PRIMARY_KEY,
                 fieldsub))

    # Setup the 'qmarks' sqlite substring
    qmarks = ','.join(['?' for i in range(len(fields))])

    # Setup the sql substring for inserting fields into database
    fieldsub = ','.join([fieldname for fieldname, role in fields])

    query = 'INSERT INTO %s (%s) VALUES (%s)' %\
            (DBConstants._DICT_TABLE, fieldsub, qmarks)
    # Pull data from the iterator and store in database
    # Commiting in batches seems faster than a single call to executemany
    data = (tuple(record[fieldname] for fieldname, role in fields) \
            for record in rcrditer)
    while True:
        batch = list(itertools.islice(data, 10000))
        if not batch: break
        cur.executemany(query, batch)
    con.commit()

    # Attribute to index
    queryby = fields[0][0] # Defaults to the first field
    for fieldname, role in fields:
        if role == DBConstants._INDEXED_TEXT_KEY:
            queryby = fieldname
            break

    # Make the index on the 'queryby' attribute
    cur.execute('CREATE UNIQUE INDEX %sidx ON %s(%s)' %
                (queryby, DBConstants._DICT_TABLE, queryby))

    con.commit()
    con.close()
