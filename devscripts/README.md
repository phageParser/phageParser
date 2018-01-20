# Development Environment Setup Scripts

These scripts are for setting up a development environment using apt (e.g., Debian or Ubuntu Linux) and pip.

1. ```devscripts/bootstrap.sh``` uses apt to get required software including python, postgres, etc...
2. ```devscripts/server.sh``` starts the postgres service
3. ```source venv/bin/activate``` to start virtualenv
4. ```devscripts/setup.sh``` uses pip to get python dependencies, creates and migrates database structure

