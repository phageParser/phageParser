#!/usr/bin/env sh
# Run the following, having activated a virtualenv
pip install --upgrade pip setuptools
pip install jupyter
pip install -r requirements.txt

# Setup Database and Schema
sudo -u postgres createdb -w testdb
./manage.py makemigrations restapi
./manage.py migrate
