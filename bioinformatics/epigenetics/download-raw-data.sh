#/bin/bash

#this goes into ~/.local/bin
uv tool install geofetch
uv tool install looper
geofetch -i sample-table.txt --add-convert-modifier --discard-soft -m `pwd` \
	&& export SRARAW=`pwd` \
	&& mkdir -p fastqs \
	&& export SRAFQ=`pwd`/fastqs \
	&& looper run --config ./sample-table/looper_config.yaml -p local --output-dir . \
	&& echo "complete"
