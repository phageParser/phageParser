#!/bin/bash

read -p "This will create a new empty test database called phagedb_test. This will not work if the database phagedb_test already exists. Proceed? " -n 1 -r
echo    # (optional) move to a new line
if [[ $REPLY =~ ^[Yy]$ ]]
then
    pg_dump phagedb -s > phage_test_db.sql
    createdb phagedb_test
    psql -d phagedb_test -f phage_test_db.sql
    rm phage_test_db.sql
fi
