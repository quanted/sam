
from __future__ import absolute_import
import pymongo as pymongo
import json
import pandas as pd
import os
import numpy as np

import requests

#from flask_qed.pram_flask.tasks import sam_status
#from pram_flask.tasks import sam_status

IN_DOCKER = "False"


class SamPostprocessor(object):

    def __init__(self,task_id):
        self.task_id = task_id
        self.sam_data = None
        self.huc8_summary = None

    def connect_to_mongoDB(self):
        if IN_DOCKER == "False":
            # Dev env mongoDB
            mongo = pymongo.MongoClient(host='mongodb://localhost:27017/0')
            print("MONGODB: mongodb://localhost:27017/0")
        else:
            # Production env mongoDB
            mongo = pymongo.MongoClient(host='mongodb://mongodb:27017/0')
            print("MONGODB: mongodb://mongodb:27017/0")
        mongo_db = mongo['pram_tasks']
        mongo.pram_tasks.Collection.create_index([("date", pymongo.DESCENDING)], expireAfterSeconds=86400)
        # ALL entries into mongo.flask_hms must have datetime.utcnow() timestamp, which is used to delete the record after 86400
        # seconds, 24 hours.
        return mongo_db

    def get_sam_data(self):
        mongo_db = self.connect_to_mongoDB()
        posts = mongo_db.posts
        data = json.loads(posts.find_one({'_id': self.task_id})["data"])
        self.sam_data = data
        return

    def calc_huc_summary(self):
        sam_properties = [x["properties"] for x in self.sam_data['features']]
        data = pd.DataFrame(sam_properties)
        data['COMID'] = data['COMID'].astype(str)
        path_to_csv = os.path.join(os.path.dirname(__file__),'HUC12_comids.csv')
        huc_comid = pd.read_csv(path_to_csv)
        huc_comid[["HUC12","COMID"]] = huc_comid[["HUC12","COMID"]].astype(str)
        huc_comid['HUC12'] = huc_comid['HUC12'].apply(replace_leading_0)
        data = data.merge(huc_comid[["COMID", "HUC12"]], on="COMID")
        data["HUC8"] = data["HUC12"].str.slice(0,8)
        try:
            huc8_summary = data.groupby('HUC8').agg(['mean','max'])
            huc8_summary.columns = ["_".join(x) for x in huc8_summary.columns.ravel()]
            self.huc8_summary = huc8_summary
        except:
            self.huc8_summary = pd.DataFrame(columns=['HUC8','acute_human_mean', 'acute_human_max'])



    def append_sam_data(self):
        mongo_db = self.connect_to_mongoDB()
        posts = mongo_db.posts
        posts.update_one({'_id': self.task_id},{'$set': {'huc8_summary': self.huc8_summary.to_json(orient='index')}})
        return



def replace_leading_0(huc_str):
    if len(huc_str) == 12:
        return huc_str
    elif len(huc_str) == 11:
        return '0'+ huc_str
    else:
        raise NameError('Number that is neither 12 digits nor 11 digits!')



