import os
import numpy as np
import pandas as pd

from collections import OrderedDict


def read_gdb(dbf_file, table_name, fields=None):
    """Reads the contents of a dbf table """
    import ogr

    # Initialize file
    driver = ogr.GetDriverByName("OpenFileGDB")
    gdb = driver.Open(dbf_file)

    # parsing layers by index
    tables = {gdb.GetLayerByIndex(i).GetName(): i for i in range(gdb.GetLayerCount())}
    table = gdb.GetLayer(tables[table_name])
    table_def = table.GetLayerDefn()
    if fields is None:
        fields = [table_def.GetFieldDefn(i).GetName() for i in range(table_def.GetFieldCount())]
    data = np.array([[row.GetField(f) for f in fields] for row in table])
    return pd.DataFrame(data=data, columns=fields)


def read_dbf(dbf_file):
    from dbfread import DBF, FieldParser

    class MyFieldParser(FieldParser):
        def parse(self, field, data):
            try:
                return FieldParser.parse(self, field, data)
            except ValueError:
                return None

    try:
        dbf = DBF(dbf_file)
        table = pd.DataFrame(iter(dbf))
    except ValueError:
        dbf = DBF(dbf_file, parserclass=MyFieldParser)
        table = pd.DataFrame(iter(dbf))
    table.rename(columns={column: column.lower() for column in table.columns}, inplace=True)
    return table


# NHD regions and the states that overlap
nhd_states = OrderedDict((('01', {"ME", "NH", "VT", "MA", "CT", "RI", "NY"}),
                          ('02', {"VT", "NY", "PA", "NJ", "MD", "DE", "WV", "DC", "VA"}),
                          ('03N', {"VA", "NC", "SC", "GA"}),
                          ('03S', {"FL", "GA"}),
                          ('03W', {"FL", "GA", "TN", "AL", "MS"}),
                          ('04', {"WI", "MN", "MI", "IL", "IN", "OH", "PA", "NY"}),
                          ('05', {"IL", "IN", "OH", "PA", "WV", "VA", "KY", "TN"}),
                          ('06', {"VA", "KY", "TN", "NC", "GA", "AL", "MS"}),
                          ('07', {"MN", "WI", "SD", "IA", "IL", "MO", "IN"}),
                          ('08', {"MO", "KY", "TN", "AR", "MS", "LA"}),
                          ('09', {"ND", "MN", "SD"}),
                          ('10U', {"MT", "ND", "WY", "SD", "MN", "NE", "IA"}),
                          ('10L', {"CO", "WY", "MN", "NE", "IA", "KS", "MO"}),
                          ('11', {"CO", "KS", "MO", "NM", "TX", "OK", "AR", "LA"}),
                          ('12', {"NM", "TX", "LA"}),
                          ('13', {"CO", "NM", "TX"}),
                          ('14', {"WY", "UT", "CO", "AZ", "NM"}),
                          ('15', {"NV", "UT", "AZ", "NM", "CA"}),
                          ('16', {"CA", "OR", "ID", "WY", "NV", "UT"}),
                          ('17', {"WA", "ID", "MT", "OR", "WY", "UT", "NV"}),
                          ('18', {"OR", "NV", "CA"})))

# All states
states = sorted(set().union(*nhd_states.values()))
