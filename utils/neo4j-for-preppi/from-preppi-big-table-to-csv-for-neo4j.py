#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 18:44:56 2018

@author: av2729
"""

import sys, csv, time
from collections import defaultdict
import subprocess
#import numpy as np
# import pandas as pd

from neo4jrestclient import client

def parsePreppiBigTable(filename,n):
    prot_dict = defaultdict(lambda: defaultdict(float))

    with open(filename, "rt") as csvfile:
        datareader = csv.reader(csvfile,delimiter='\t')
        next(datareader) # Skipping header
        count = 0
        p1_count = 0
        p2_count = 0
        for row in datareader:
            if count < n:
                p = row[0].split("|")
#                print("- Protein 1 : " , p[0] )
#                print("- Protein 2 : " , p[1] )
#                print("-- SCORE: " , row[3] )
                count += 1 
                prot_dict[p[0]][p[1]] = row[3]
                
#                prot1_dict[p[0]] = row[3]
#                prot2_dict[p[1]] = row[3] 
#                if p1 in p[0]:
#                    p1_count += 1 
#                if p2 in p[1]:
#                    p2_count += 1

                if ((count % 1e6) == 0 ):
                    print(">>> *** Number of lines analyzed : " , count/1e6 , " Millions ***")                 
            else:
                break
#    print(prot1_dict.keys())

#    print( " - P1 as " , p1 , " found " , p1_count , " times" )
#    print( " - P2 as " , p2 , " found " , p2_count , " times" )
    print( "\n -- Read " , count , " lines" )
    print( " --- Size Hash P : " , sys.getsizeof(prot_dict) )
    # print( " --- Size Hash P2 : " , sys.getsizeof(prot2_dict) )
                
    # return( prot1_dict , prot2_dict )
    return( prot_dict )
                    
    
def parsePreppiTable_600(filename,n):
    prot_dict_total_score = defaultdict(lambda: defaultdict(list))

    with open(filename, "rt") as csvfile:
        datareader = csv.reader(csvfile,delimiter='\t')
        next(datareader) # Skipping header
        count = 0
        for row in datareader:
            if count < n:
                p1 = row[0]
                p2 = row[1]
                prot_dict_total_score[p1][p2].append( row[10] ) # PrePPI Total Score
                prot_dict_total_score[p1][p2].append( row[13] ) # PrePPI Experimental Score
                prot_dict_total_score[p1][p2].append( row[14] ) # PrePPI Final Score
                count += 1 

                if ((count % 1e3) == 0 ):
                    print(">>> *** Number of lines analyzed : " , count/1e3 , " Thousands ***")                 
            else:
                break
            
    print( "\n -- Read " , count , " lines" )
    print( " --- Size Hash P : " , sys.getsizeof(prot_dict_total_score) )
    
    return( prot_dict_total_score )    
           
##
# Main Function
# -------------
def main(argv):
    
    # preppi_big_table_filename = "/Users/av2729/Downloads/preppi-data/PrePPI_scores.txt"
    preppi_big_table_filename = "/ifs/scratch/c2b2/ac_lab/av2729/preppi-db-building/PrePPI_scores.txt"
    # preppi_table_600_filename = "/Users/av2729/Downloads/preppi-data/preppi_final600.txt"
    preppi_table_600_filename = "/ifs/scratch/c2b2/ac_lab/av2729/preppi-db-building/preppi_final600.txt"

    # preppi_big_table_Neo4j_nodes_filename = "/Users/av2729/Downloads/preppi-data/preppi-for-neo4j-nodes.csv"
    # preppi_big_table_Neo4j_relationships_filename = "/Users/av2729/Downloads/preppi-data/preppi-for-neo4j-relationships.csv"
    preppi_big_table_Neo4j_nodes_filename = "/ifs/scratch/c2b2/ac_lab/av2729/preppi-db-building/preppi-for-neo4j-nodes.csv"
    preppi_big_table_Neo4j_relationships_filename = "/ifs/scratch/c2b2/ac_lab/av2729/preppi-db-building/preppi-for-neo4j-relationships.csv"
    # preppi_table_600_Neo4j_relationships_filename = "/Users/av2729/Downloads/preppi-data/preppi-for-neo4j-table-600-relationships.csv"
    # csv_importer_script_filename = "/Users/av2729/Workspace/vaxtools/python/neo4j-for-preppi/import-database-from-csv-to-neo4j.sh"
    csv_importer_script_filename = "/ifs/scratch/c2b2/ac_lab/av2729/preppi-db-building/vaxtools/python/neo4j-for-preppi/import-database-from-csv-to-neo4j.sh"
    
    n_lines = 1e9
    # n_lines = 1e6
    
    # a = "P05067"
    # b = "O75509"
    
    print('---')
    print('>>> Starting reading PrePPI Big Table located in: ', preppi_big_table_filename )
    time_start = time.clock()
    # p1,p2 = parsePreppiTable( filename = preppi_big_table_filename , n = n_lines )
    x = parsePreppiBigTable( filename = preppi_big_table_filename , n = n_lines )
    time_elapsed = (time.clock() - time_start)
    print('>>> >> Reading Big Table Elapsed Time: ', time_elapsed )
    print('---')
    
    print('---')
    print('>>> Starting reading PrePPI (score > 600) Table located in: ', preppi_table_600_filename )
    time_start = time.clock()
    x_with_exp_scores = parsePreppiTable_600( filename = preppi_table_600_filename , n = n_lines )
    time_elapsed = (time.clock() - time_start)
    print('>>> >> Reading (score > 600) Table Elapsed Time: ', time_elapsed )
    print('---')       
    
    # from itertools import islice
    # list( islice( x_with_exp_scores.items() , 1 ) )

    
    
    
    # db = client.GraphDatabase( "http://localhost:7474/", username="preppi", password="preppi" )
    # q = "MATCH (n) DETACH DELETE n"
    # results = db.query(q, returns=(client.Node, str, client.Node))
    # shutil.rmtree( "/Users/av2729/Workspace/neo4j/databases/preppi.graphdb" , ignore_errors = False, onerror = None )    
    
    
    # protein = db.labels.create("Protein")
    # for p1_keys , p1_values in x.items():
    #     p1 = db.nodes.create(uniprot_id=p1_keys)
    #     protein.add(p1)
    #     for p2_keys , p2_values in p1_values.items():
    #         # print(p1_keys , " -> " , p2_keys , " : " , p2_values ) 
    #         p2 = db.nodes.create(uniprot_id=p2_keys)
    #         p1.relationships.create( "interact with", p2 , score = p2_values )
    #         protein.add(p2)
    
    print('---')
    print('>>> Starting writing Neo4J Nodes in: ', preppi_big_table_Neo4j_nodes_filename )  
    time_start = time.clock()
    with open(preppi_big_table_Neo4j_nodes_filename, 'wt') as csvfile:
        csvwriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        # csvwriter.writerow( ['NODE_ID','name',':LABEL'] )
        csvwriter.writerow( ['NODE_ID:ID','name'] )
        counter = 0
        unique_prot_ids = set()
        for p1_keys , p1_values in x.items():
            unique_prot_ids.add(p1_keys)
        for p in sorted(unique_prot_ids):
            counter += 1
            csvwriter.writerow( [p] + [p] )
            if ((counter % 1e3) == 0 ):
                print(">>> *** Number of lines wrote to file : " , counter/1e3 , " Thousands ***")
            
    # print('>>> Starting writing PrePPI Big Table Neo4J Relationships in: ', preppi_big_table_Neo4j_relationships_filename )              
    # with open(preppi_big_table_Neo4j_relationships_filename, 'wt') as csvfile:
    #     csvwriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    #     # csvwriter.writerow( [':START_ID','score:float',':END_ID',':TYPE'] )
    #     csvwriter.writerow( [':START_ID','score:float',':END_ID'] )
    #     counter = 0
    #     for p1_keys , p1_values in x.items():
    #         if ((counter % 1e3) == 0 ):
    #             print(">>> *** Number of protein relationships wrote to file : " , counter/1e3 , " Thousands ***")
    #         counter += 1
    #         for p2_keys , p2_values in p1_values.items():
    #             csvwriter.writerow( [p1_keys] + [1 if "NULL" in p2_values else p2_values] + [p2_keys] )
    #             # csvwriter.writerow( [p1_keys] + [p2_values] + [p2_keys] + ["INTERACT_WITH"] )            
                
    # print('>>> Starting writing PrePPI Table 600 Neo4J Relationships in: ', preppi_table_600_Neo4j_relationships_filename )              
    # with open(preppi_table_600_Neo4j_relationships_filename, 'wt') as csvfile:
    #     csvwriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    #     # csvwriter.writerow( [':START_ID','total_score:float','experimental_score:float','final_score:float',':END_ID',':TYPE'] )
    #     csvwriter.writerow( [':START_ID','total_score:float','experimental_score:float','final_score:float',':END_ID'] )
    #     counter = 0
    #     for p1_keys , p1_values in x_with_exp_scores.items():
    #         if ((counter % 1e3) == 0 ):
    #             print(">>> *** Number of protein relationships wrote to file : " , counter/1e3 , " Thousands ***")
    #         counter += 1
    #         for p2_keys , p2_values in p1_values.items():
    #             csvwriter.writerow( [p1_keys] + [1.0 if not aScore else float(aScore) for aScore in p2_values] + [p2_keys] )
        
    print('>>> Starting writing PrePPI Tables Neo4J Relationships in: ', preppi_big_table_Neo4j_relationships_filename )              
    with open(preppi_big_table_Neo4j_relationships_filename, 'wt') as csvfile:
        csvwriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        # csvwriter.writerow( [':START_ID','total_score:float','experimental_score:float','final_score:float',':END_ID',':TYPE'] )
        csvwriter.writerow( [':START_ID','total_score_B:float','total_score:float','experimental_score:float','final_score:float',':END_ID'] )
        counter = 0
        for p1_keys , p1_values in x.items():
            if ((counter % 1e3) == 0 ):
                print(">>> *** Number of protein relationships wrote to file : " , counter/1e3 , " Thousands ***")
            counter += 1
            for p2_keys , p2_values in p1_values.items():
                tmp = x_with_exp_scores[p1_keys][p2_keys]
                if tmp:
                    csvwriter.writerow( [p1_keys] + [1.0 if "NULL" in p2_values else float(p2_values)] + [1.0 if not aScore else float(aScore) for aScore in tmp] + [p2_keys] )
                else:
                    csvwriter.writerow( [p1_keys] + [1.0 if "NULL" in p2_values else float(p2_values)] + [0.0] + [0.0] + [0.0]+ [p2_keys] )


    # for p1_keys , p1_values in x.items():
    #     p1 = db.nodes.create(uniprot_id=p1_keys)
    #     protein.add(p1)
    #     for p2_keys , p2_values in p1_values.items():
    #         # print(p1_keys , " -> " , p2_keys , " : " , p2_values ) 
    #         p2 = db.nodes.create(uniprot_id=p2_keys)
    #         p1.relationships.create( "interact with", p2 , score = p2_values )
    #         protein.add(p2)    
        
    time_elapsed = (time.clock() - time_start)
    print('>>> >> Writing Neo4J CSV Tables Elapsed Time: ', time_elapsed )
    print('---')    
    
    
    # print('---')
    # print('>>> Starting importing Neo4J Tables ')
    # time_start = time.clock()
    # subprocess.call([csv_importer_script_filename])
    # time_elapsed = (time.clock() - time_start)
    # print('>>> >> Importing CSV into NEO4J DBMS Elapsed Time: ', time_elapsed )
    # print('---')    

if __name__ == "__main__":
   main(sys.argv[1:])
