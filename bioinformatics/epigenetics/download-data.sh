#/bin/bash

rootPath='/home/stanfish/Git/meissner-nature-2015/data'
uvx geofetch -i GSE61475 --just-metadata -m $rootPath
uvx geofetch -i sample-table.txt --processed --geo-folder $rootPath \
	&& mv sample-table $rootPath/GSE61475 \
	&& echo "complete"
