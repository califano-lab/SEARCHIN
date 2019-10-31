#!/bin/bash

## To import the database
# rm -R "/Applications/Neo4j Community Edition 3.2.5.app/Contents/Resources/app/data/databases/preppi.graphdb"
PREPPI_NEO4J_DB="/Users/av2729/Workspace/neo4j/databases/preppi.graphdb"
echo -e "\033[1;97m\033[41m *** Deleting previous database *** \033[0m"
echo -e "\033[1;97m\033[41m --- ${PREPPI_NEO4J_DB} --- \033[0m"
rm -R ${PREPPI_NEO4J_DB}
NEO4J_HOME="/Users/av2729/Workspace/neo4j/"
/Applications/Neo4j\ Community\ Edition\ 3.2.5.app/Contents/Resources/app/bin/neo4j-admin import \
                --additional-config=/Users/av2729/Documents/Neo4j/.neo4j.conf \
                --database=preppi.graphdb --id-type string \
                 --nodes:Protein /Users/av2729/Downloads/preppi-data/preppi-for-neo4j-nodes.csv  \
                 --relationships:interacts_with /Users/av2729/Downloads/preppi-data/preppi-for-neo4j-relationships.csv \
                 --ignore-missing-nodes=true