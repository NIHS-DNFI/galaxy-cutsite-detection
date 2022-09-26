#!/usr/bin/env bash

echo "Checking /export/galaxy-central/database/ ..."

if [ ! -d /export/galaxy-central/database/shed_tools ]; then
	echo "Coping /local_tools/shed_tools to /export/galaxy-central/database"
	cp -r /local_tools/shed_tools /export/galaxy-central/database/
	chown -R galaxy.galaxy /export/galaxy-central/database/shed_tools
fi

echo "Starting /usr/bin/startup ..."
/usr/bin/startup
