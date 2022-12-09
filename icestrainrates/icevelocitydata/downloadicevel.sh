#!bin/bash

export SERVER=https://dataverse.geus.dk
export DOI=10.22008/promice/data/sentinel1icevelocity/greenlandicesheet

curl ${SERVER}/api/datasets/:persistentId?persistentId=doi:${DOI} > dv.json
cat dv.json | tr ',' '\n' | grep -E '"persistentId"' | cut -d'"' -f4 > urls.txt
while read -r PID; do
    curl -O -J $SERVER/api/access/datafile/:persistentId?persistentId=${PID}
done < urls.txt
rm dv.json urls.txt # cleanup
