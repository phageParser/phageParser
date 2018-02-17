#!/usr/bin/env sh
sudo apt install python3 python3-dev python3-venv postgresql postgresql-server-dev-9.5
python3 -m venv venv

# Note, if you're behind an SSL content inspecting firewall, try this in venv/pip.conf
#```
#[global]
#cert=/usr/local/share/ca-certificates/some_cert.pem
#```
# cp scripts/pip.conf venv/pip.conf

echo 'Next Step: source venv/bin/activate'
