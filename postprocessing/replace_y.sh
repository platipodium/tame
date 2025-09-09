#!/bin/bash
#
# Script to update numerical values in YAML files while preserving comments
# Usage: ./update_yaml.sh <yaml_file> <token> <new_value>
#
# SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum hereon GmbH
# SPDX-FileContributor: Kai Wirtz <kai.wirtz@hereon.de>
# SPDX-License-Identifier: CC0-1.0
#

set -e  # Exit on any error

# Function to display usage information
show_usage() {
    echo "Usage: $0 <yaml_file> <token> <new_value>"
    echo ""
    echo "Arguments:"
    echo "  yaml_file  - Path to the YAML file to modify"
    echo "  token      - Name of the variable to update"
    echo "  new_value  - New numerical value to set"
    echo ""
    echo "Example:"
    echo "  $0 config.yaml port 8080"
    echo "  $0 settings.yml timeout 30"
    exit 1
}

# Check if correct number of arguments provided
if [ $# -ne 3 ]; then
    echo "Error: Invalid number of arguments"
    show_usage
fi

YAML_FILE="$1"
TOKEN="$2"
NEW_VALUE="$3"

# Validate inputs
if [ ! -f "$YAML_FILE" ]; then
    echo "Error: File '$YAML_FILE' does not exist"
    exit 1
fi

# Check if new value is numerical
if ! [[ "$NEW_VALUE" =~ ^-?[0-9]+\.?[0-9]*$ ]]; then
    echo "Error: '$NEW_VALUE' is not a valid numerical value"
    exit 1
fi

# Escape special characters in token for sed
ESCAPED_TOKEN=$(echo "$TOKEN" | sed 's/[[\.*^$()+?{|]/\\&/g')

# Create backup of original file
cp "$YAML_FILE" "${YAML_FILE}.backup"

# Use sed to replace the value while preserving comments
# This pattern matches:
# - Optional whitespace at start of line
# - The token followed by a colon
# - Optional whitespace
# - Any existing value (numbers, including decimals and negatives)
# - Captures any trailing content (comments, etc.)
sed -i.tmp "s/^\([[:space:]]*${ESCAPED_TOKEN}[[:space:]]*:[[:space:]]*\)[0-9.-]*\(.*\)$/\1${NEW_VALUE}\2/" "$YAML_FILE"

# Check if the replacement was successful by comparing files
if cmp -s "$YAML_FILE" "${YAML_FILE}.backup"; then
    echo "Warning: No changes made. Token '$TOKEN' may not exist or may not have a numerical value."
    echo "Please check that the token exists and has the format 'token: value'"
    rm "${YAML_FILE}.tmp" 2>/dev/null || true
    exit 1
else
    echo "Successfully updated '$TOKEN' to '$NEW_VALUE' in '$YAML_FILE'"
 #   echo "Backup saved as '${YAML_FILE}.backup'"
fi

# Clean up temporary file
rm "${YAML_FILE}.tmp" 2>/dev/null || true
