#!/usr/bin/env python

from MySQLdb import connect

CREATE_BIOM_TABLE = """CREATE TABLE %s (
    OBS_ID VARCHAR(32) NOT NULL,
    SAMPLE_ID VARCHAR(32) NOT NULL,
    VAL %s NOT NULL,
    PRIMARY KEY(OBS_ID, SAMPLE_ID))"""
    
COUNT_TABLES = """SELECT COUNT(TABLE_NAME) 
    FROM information_schema.tables 
    WHERE table_schema='%s' AND TABLE_NAME LIKE '%s'"""
CREATE_BIOM_TABLE_FROM = """CREATE TABLE %s
    SELECT * FROM %s"""
CREATE_BIOM_TABLE_FROM_SUBSET = """CREATE TABLE %s
    SELECT * FROM %s
    WHERE %s IN %s"""
GET_ITEM = """SELECT VAL 
    FROM %s
    WHERE OBS_ID='%s' AND SAMPLE_ID='%s'"""
SET_ITEM = """INSERT INTO %s(OBS_ID, SAMPLE_ID, VAL) 
    VALUES('%s','%s',%f)"""
UPDATE_ITEM = """UPDATE %s SET VAL=%d 
    WHERE OBS_ID='%s' AND SAMPLE_ID='%s'"""
SAMPLE_IN_OBS = "SELECT OBS_ID, SAMPLE_ID, VAL FROM %s WHERE OBS_ID='%s'"
OBS_IN_SAMPLE = "SELECT OBS_ID, SAMPLE_ID, VAL FROM %s WHERE SAMPLE_ID='%s'"
DELETE_ITEM = "DELETE FROM %s WHERE OBS_ID='%s' AND SAMPLE_ID='%s'"

# mysql doesnt implement MINUS... lame.
TABLE_VALUE_EQUALITY = """SELECT count(*) 
FROM %s a LEFT JOIN %s b on a.OBS_ID=b.OBS_ID AND a.SAMPLE_ID=b.SAMPLE_ID
WHERE a.VAL != b.VAL"""

TABLE_LABEL_EQUALITY = """SELECT count(*) 
FROM %s a LEFT JOIN %s b ON a.%s=b.%s 
WHERE b.%s IS NULL"""

DROP_TABLE = "DROP TABLE %s"                    
# http://weblogs.sqlteam.com/jeffs/archive/2004/11/10/2737.aspx
#TABLE_EQUALITY = """SELECT MIN(%s) as %s, OBS_ID, SAMPLE_ID, VAL
#FROM
#    (SELECT 'Table A' as %s, A.OBS_ID, A.SAMPLE_ID, A.VAL
#     FROM A
#     UNION ALL
#     SELECT 'Table B' as %S, B.OBS_ID, B.SAMPLE_ID, B.VAL
#     FROM B
#    ) tmp
#GROUP BY OBS_ID, SAMPLE_ID, VAL
#HAVING COUNT(*) = 1
#ORDER BY OBS_ID"""


class BiomDB(object):
    def __init__(self, hostname='localhost', user='mcdonald', password='stupid',
                 db='biom', table_basename='tmp_biomtable'):
        self.con = connect(hostname, user, password, db)
        self.db = db
        self._tmp_tables = []
        self._table_basename = table_basename
        self._num_tmp_tables = self._tmp_table_count()
        
    def _tmp_table_count(self):
        """Gets the count of tables in the db"""
        cursor = self.con.cursor()
        cursor.execute(COUNT_TABLES % (self.db, self._table_basename + '%'))
        return cursor.fetchone()[0]
        
    def _get_new_tableid(self):
        """Get a new unique table id"""
        self._num_tmp_tables += 1
        new_table_id = '_'.join([self._table_basename, str(self._num_tmp_tables)])
        return new_table_id
            
    def createTempTable(self, from_table=None, from_samples=None, 
            from_obs=None, dtype='INT'):
        """Create a temporary table
        
        if from_table is None, create empty structure else copy from_table
            - just pull samples or obs optionally
        """
        new_table_id = self._get_new_tableid()
        cursor = self.con.cursor()
        print from_table
        if from_table is None:
            cursor.execute(CREATE_BIOM_TABLE % (new_table_id, dtype))
        else:
            if from_samples is not None and from_obs is not None:
                raise ValueError, "Can only handle 1 axis at a time"
            if from_samples is not None:
                wherein = "('%s')" % "','".join(from_samples)
                cursor.execute(CREATE_BIOM_TABLE_FROM_SUBSET % (new_table_id, 
                                from_table, 'SAMPLE_ID', wherein))
            elif from_obs is not None:
                wherein = "('%s')" % "','".join(from_obs)
                cursor.execute(CREATE_BIOM_TABLE_FROM_SUBSET % (new_table_id, from_table, 'OBS_ID', wherein))
            else:
                print CREATE_BIOM_TABLE_FROM % (new_table_id, from_table)
                cursor.execute(CREATE_BIOM_TABLE_FROM % (new_table_id, from_table))
        self._tmp_tables.append(new_table_id)
        self.con.commit()
        return new_table_id

    def getItem(self, table_name, obs_id, sample_id):
        """Get a specific value from a table"""
        cursor = self.con.cursor()
        cursor.execute(GET_ITEM % (table_name, obs_id, sample_id))
        res = cursor.fetchone()
        if not res:
            return 0
        else:
            return res[0]

    def setItem(self, table_name, obs_id, sample_id, value):
        """Set a specific value in a table"""
        if value == 0:
            self.deleteItem(table_name, obs_id, sample_id)
            return
            
        cursor = self.con.cursor()
        cursor.execute(GET_ITEM % (table_name, obs_id, sample_id))
        res = cursor.fetchone()
        if not res:
            cursor.execute(SET_ITEM % (table_name, obs_id, sample_id, value))
        else:
            cursor.execute(UPDATE_ITEM % (table_name, value, obs_id, sample_id))
        self.con.commit()
        
    def getSampleByObs(self, table_name, obs_id):
        """Get all samples associated with an observation"""
        cursor = self.con.cursor()
        cursor.execute(SAMPLE_IN_OBS % (table_name, obs_id))
        return cursor.fetchall()    

    def getObsBySample(self, table_name, samp_id):
        """Get all observations associated with a sample"""
        cursor = self.con.cursor()
        cursor.execute(OBS_IN_SAMPLE % (table_name, samp_id))
        return cursor.fetchall()    

    def deleteItem(self, table_name, obs_id, sample_id):
        """Delete a specific item in a table"""
        cursor = self.con.cursor()
        cursor.execute(GET_ITEM % (table_name, obs_id, sample_id))
        res = cursor.fetchone()
        if not res:
            return
        else:
            cursor.execute(DELETE_ITEM % (table_name, obs_id, sample_id))
            self.con.commit()
            
    def tableEquality(self, table_a, table_b):
        """Check for table equality"""
        cursor = self.con.cursor()
        
        # UGH...
        
        cursor.execute(TABLE_VALUE_EQUALITY % (table_a, table_b))
        res = cursor.fetchone()[0]
        
        if res != 0:
            return False
        
        # mysql doesn't do full outer joins
        # table A -> table b
        cursor.execute(TABLE_LABEL_EQUALITY % (table_a, table_b, 'OBS_ID', 
                                                'OBS_ID','OBS_ID'))
        res = cursor.fetchone()[0]
        
        if res != 0:
            return False
        
        # table B -> table A
        cursor.execute(TABLE_LABEL_EQUALITY % (table_a, table_b, 'OBS_ID', 
                                                    'OBS_ID','OBS_ID'))
        res = cursor.fetchone()[0]

        if res != 0:
            return False
    
        # table A -> table B 
        cursor.execute(TABLE_LABEL_EQUALITY % (table_a, table_b, 'SAMPLE_ID', 
                                                    'SAMPLE_ID','SAMPLE_ID'))
        res = cursor.fetchone()[0]

        if res != 0:
            return False

        # table B -> table A 
        cursor.execute(TABLE_LABEL_EQUALITY % (table_b, table_a, 'SAMPLE_ID', 
                                                        'SAMPLE_ID','SAMPLE_ID'))
        res = cursor.fetchone()[0]

        if res != 0:
            return False

        return True
    
    def dropTable(self,table=None):
        """Drop a table or all tables if None"""
        cursor = self.con.cursor()
        
        if table is None:
            for t in self._tmp_tables:
                cursor.execute(DROP_TABLE % t)
            self._tmp_tables = []
        else:
            cursor.execute(DROP_TABLE % table)
            self._tmp_tables.pop(self._tmp_tables.index(table))
            
        self.con.commit()
            
            
