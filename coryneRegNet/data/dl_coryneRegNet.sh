#!/bin/bash
curl -O http://coryneregnet.mmci.uni-saarland.de/v6/coryneregnet_database_content.zip
unzip coryneregnet_database_content.zip
rename -n "s/ /_/g" *
