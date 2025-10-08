#!/bin/bash
# get_scripts.sh - Create links and copy updated scripts
#
# SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum hereon GmbH
# SPDX-FileContributor: Kai Wirtz <kai.wirtz@hereon.de>
# SPDX-License-Identifier: CC0-1.0

POSTPROC_DIR="../../postprocessing"

# List of files to create symbolic links for (easily extensible)
LINK_FILES=("replace_y.sh" "read_nc_simple.m" "changeNsee.sh")

# Check if postprocessing directory exists
[ ! -d "$POSTPROC_DIR" ] && { echo "Error: $POSTPROC_DIR does not exist"; exit 1; }

#echo "Creating links to scripts in $POSTPROC_DIR..."

# Create symbolic links
for file in "${LINK_FILES[@]}"; do
    if [ -f "$POSTPROC_DIR/$file" ]; then
        [ -L "$file" ] && rm "$file"
        ln -s "$POSTPROC_DIR/$file" "$file"
        echo "Created link: $file"
    else
        echo "Warning: $POSTPROC_DIR/$file not found"
    fi
done

# Handle cmp_0Dres.m - backup and copy
if [ -f "$POSTPROC_DIR/cmp_0Dres.m" ]; then
    if [ -f "cmp_0Dres.m" ]; then
        mkdir -p old_scripts
        DATE_SUFFIX=$(date +%Y%m%d_%H%M)
        BACKUP_NAME="cmp_0Dres_${DATE_SUFFIX}.m"
        mv cmp_0Dres.m "old_scripts/$BACKUP_NAME"
        echo "Backed up as: old_scripts/$BACKUP_NAME"
    fi
    cp "$POSTPROC_DIR/cmp_0Dres.m" .
    echo "Copied new cmp_0Dres.m"
else
    echo "Error: $POSTPROC_DIR/cmp_0Dres.m not found"
    exit 1
fi
#echo "Done!"
